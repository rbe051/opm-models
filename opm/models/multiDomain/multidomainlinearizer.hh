#ifndef OPM_MULTI_DOMAIN_LINEARIZER_HH
#define OPM_MULTI_DOMAIN_LINEARIZER_HH

#include <opm/models/utils/propertysystem.hh>

#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/multiDomain/couplingelementcontext.hh>
#include <opm/models/multiDomain/multidomainproperties.hh>
#include <opm/material/densead/Math.hpp>

#include <dune/grid/common/intersection.hh>
#include <dune/common/indices.hh>

BEGIN_PROPERTIES
NEW_PROP_TAG(CouplerTypeTag);
END_PROPERTIES

namespace Opm
{

template <class TypeTag>
class MultiDomainLinearizer
{
    template <std::size_t i>
    using Coupler = typename GET_PROP_TYPE(TypeTag, CouplerTypeTag)::template Coupler<i>;
    template <std::size_t i>
    using Simulator = typename GET_PROP_TYPE(TypeTag, SubTypeTag)::template Simulator<i>;
    template <std::size_t i>
    using Stencil = typename GET_PROP_TYPE(TypeTag, SubTypeTag)::template Stencil<i>;
    template <std::size_t i>
    using GridView = typename GET_PROP_TYPE(TypeTag, SubTypeTag)::template GridView<i>;
    typedef typename GET_PROP_TYPE(TypeTag, SubTypeTag)::JacobianMatrix JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, SubTypeTag)::GlobalEqVector GlobalEqVector;

    template <std::size_t i>
    using Element = typename GridView<i>::template Codim<0>::Entity;
    using Simulators = typename GET_PROP_TYPE(TypeTag, SubTypeTag)::template TupleOfSharedPtr<Simulator>;
    using Couplers = typename GET_PROP_TYPE(TypeTag, CouplerTypeTag)::template TupleOfSharedPtr<Coupler>;

public:
    MultiDomainLinearizer()
    {
    }
    /*!
     * \brief Initialize the linearizer.
     *
     * At this point we can assume that all objects in the simulator
     * have been allocated. We cannot assume that they are fully
     * initialized, though.
     *
     * \copydetails Doxygen::simulatorParam
     */
    void init(Simulators &simulators, Couplers &couplers)
    {
        simulators_ = &simulators;
        couplers_ = &couplers;
    }

    void resetSystem_()
    {
        residual_ = 0.0;
        // zero all matrix entries
        jacobian_ = 0.0;
    }

    // Construct the BCRS matrix for the Jacobian of the residual function
    /*!
     * \brief Sets the jacobian build mode
     */
    void setJacobianBuildMode()
    {
        using namespace Dune::Hybrid;
        forEach(jacobian_, [](auto &jacRow) {
            forEach(jacRow, [](auto &jacBlock) {
                using BlockType = std::decay_t<decltype(jacBlock)>;
                if (jacBlock.buildMode() == BlockType::BuildMode::unknown)
                    jacBlock.setBuildMode(BlockType::BuildMode::random);
                else if (jacBlock.buildMode() != BlockType::BuildMode::random)
                    DUNE_THROW(Dune::NotImplemented, "Only BCRS matrices with random build mode are supported at the moment");
            });
        });
    }

    template <std::size_t i, std::size_t j, class Set>
    void reserve_(const std::vector<Set> &sparsityPattern)
    {
        Dune::index_constant<i> I;
        Dune::index_constant<j> J;
        auto rows_ = model_<i>().numTotalDof();
        auto cols_ = model_<j>().numTotalDof();

        // make sure sparsityPattern is consistent with number of rows
        assert(rows_ == sparsityPattern.size());

        // allocate space for the rows of the matrix
        for (size_t dofIdx = 0; dofIdx < rows_; ++dofIdx)
            jacobian_[I][J].setrowsize(dofIdx, sparsityPattern[dofIdx].size());
        jacobian_[I][J].endrowsizes();

        // fill the rows with indices. each degree of freedom talks to
        // all of its neighbors. (it also talks to itself since
        // degrees of freedom are sometimes quite egocentric.)
        for (size_t dofIdx = 0; dofIdx < rows_; ++dofIdx)
        {
            auto nIt = sparsityPattern[dofIdx].begin();
            auto nEndIt = sparsityPattern[dofIdx].end();
            for (; nIt != nEndIt; ++nIt)
                jacobian_[I][J].addindex(dofIdx, *nIt);
        }
        jacobian_[I][J].endindices();
    }

    void createMatrix()
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(jacobian_)), [&](const auto domainI) {
            forEach(integralRange(Dune::Hybrid::size(jacobian_[domainI])), [&](const auto domainJ) {
                setJacobianPattern_(domainI, domainJ);
            });
        });
    }

    template <std::size_t i, std::size_t j>
    void setJacobianPattern_(Dune::index_constant<i> I, Dune::index_constant<j> J)
    {
        typedef std::set<unsigned> NeighborSet;
        std::vector<NeighborSet> sparsityPattern(model_<i>().numTotalDof());

        auto rows_ = model_<i>().numTotalDof();
        auto cols_ = model_<j>().numTotalDof();

        jacobian_[I][J].setSize(rows_, cols_);
        residual_[I].resize(rows_);
        if (i == j)
        {
            addDiagonalPattern_<i>(sparsityPattern);
        }
        else
        {
            addOffDiagonalPattern_<i, j>(sparsityPattern);
        }

        // allocate raw matrix

        // create matrix structure based on sparsity pattern
        reserve_<i, j>(sparsityPattern);
    }

    template <std::size_t i, std::size_t j, class Set>
    void addOffDiagonalPattern_(std::vector<Set> &sparsityPattern)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(*couplers_)), [&](const auto K) {
            const auto I = std::get<K>(*couplers_)->template subDomainIndex<0>();
            const auto J = std::get<K>(*couplers_)->template subDomainIndex<1>();
            if (I == i && J == j)
                addCouplerPattern_<K>(false, sparsityPattern);
            if (I == j && J == i)
                addCouplerPattern_<K>(true, sparsityPattern);
        });
        return;
    }

    template <std::size_t k, class Set>
    void addCouplerPattern_(bool swap, Set &sparsityPattern)
    {
        const auto &model1 = coupler<k>().template model0();
        const auto &model2 = coupler<k>().template model1();
        const auto &mortarView = coupler<k>().mortarView();
        const auto &map = coupler<k>().map_;
        // for the main model, find out the global indices of the neighboring degrees of
        // freedom of each primary degree of freedom
        typedef std::set<unsigned> NeighborSet;
        int numDofs;
        if (!swap)
            numDofs = model1.numTotalDof();
        else
            numDofs = model2.numTotalDof();

        const auto &idxSet = mortarView.indexSet();
        for (const auto &elem : elements(mortarView))
        {
            const auto idx = idxSet.index(elem);
            //const auto face = model1.gridView().grid().entity(map_->[i][idx]);
            const auto element1 = map->template toElement<0>(elem);
            unsigned elIdx1 = model1.gridView().indexSet().index(element1);

            const auto element2 = map->template toElement<1>(elem);
            unsigned elIdx2 = model2.gridView().indexSet().index(element2);
            if (!swap)
                sparsityPattern[elIdx1].insert(elIdx2);
            else
                sparsityPattern[elIdx2].insert(elIdx1);
        }
    }

    template <std::size_t i, class Set>
    void addDiagonalPattern_(std::vector<Set> &sparsityPattern) const
    {
        // for the main model, find out the global indices of the neighboring degrees of
        // freedom of each primary degree of freedom
        Stencil<i> stencil(gridView_<i>(), model_<i>().dofMapper());

        auto elemIt = gridView_<i>().template begin<0>();
        const auto elemEndIt = gridView_<i>().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            const Element<i> &elem = *elemIt;
            stencil.update(elem);

            for (unsigned primaryDofIdx = 0; primaryDofIdx < stencil.numPrimaryDof(); ++primaryDofIdx)
            {
                unsigned myIdx = stencil.globalSpaceIndex(primaryDofIdx);

                for (unsigned dofIdx = 0; dofIdx < stencil.numDof(); ++dofIdx)
                {
                    unsigned neighborIdx = stencil.globalSpaceIndex(dofIdx);
                    sparsityPattern[myIdx].insert(neighborIdx);
                }
            }
        }

        // add the additional neighbors and degrees of freedom caused by the auxiliary
        // equations
        size_t numAuxMod = model_<i>().numAuxiliaryModules();
        for (unsigned auxModIdx = 0; auxModIdx < numAuxMod; ++auxModIdx)
            model_<i>().auxiliaryModule(auxModIdx)->addNeighbors(sparsityPattern);
    }

    void linearizeSubDomains()
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(jacobian_)), [&](const auto I) {
            std::get<I>(*simulators_)->model().linearizer().linearizeDomain();
            std::get<I>(*simulators_)->model().linearizer().linearizeAuxiliaryEquations();
            std::get<I>(*simulators_)->model().linearizer().finalize();
        });
        forEach(integralRange(Dune::Hybrid::size(*couplers_)), [&](const auto I) {
            std::get<I>(*couplers_)->linearize();
        });

        forEach(integralRange(Dune::Hybrid::size(jacobian_)), [&](const auto I) {
            jacobian_[I][I] += model_<I>().linearizer().jacobian().istlMatrix();
            residual_[I] += model_<I>().linearizer().residual();
        });

        forEach(integralRange(Dune::Hybrid::size(*couplers_)), [&](const auto K) {
            const auto I = std::get<K>(*couplers_)->template subDomainIndex<0>();
            const auto J = std::get<K>(*couplers_)->template subDomainIndex<1>();
            addToJacobian(I, J, std::get<K>(*couplers_)->jacobian());
            addToResidual(I, J, std::get<K>(*couplers_)->residual());
        });
    }

    template <std::size_t i, std::size_t j, class Jac>
    void addToJacobian(Dune::index_constant<i> I, Dune::index_constant<j> J, const Jac &localJac)
    {
        Dune::index_constant<0> _0;
        Dune::index_constant<1> _1;
        jacobian_[I][I] += localJac[_0][_0];
        jacobian_[I][J] += localJac[_0][_1];
        jacobian_[J][I] += localJac[_1][_0];
        jacobian_[J][J] += localJac[_1][_1];
    }

    template <std::size_t i, std::size_t j, class Res>
    void addToResidual(Dune::index_constant<i> I, Dune::index_constant<j> J, const Res &localRes)
    {
        Dune::index_constant<0> _0;
        Dune::index_constant<1> _1;
        residual_[I] += localRes[_0];
        residual_[J] += localRes[_1];
    }

    void linearize()
    {
        if (firstIteration_)
        {
            initFirstIteration_();
        }
        resetSystem_(); // reset the global linear system of equations.
        linearizeSubDomains();
    }

    void initFirstIteration_()
    {
        firstIteration_ = false;
        setJacobianBuildMode();
        createMatrix();
    }

    GlobalEqVector &residual()
    {
        return residual_;
    }

    template <std::size_t i>
    Simulator<i> &simulator() const
    {
        return *(std::get<i>(*simulators_));
    }

    template <std::size_t i>
    Coupler<i> &coupler() const
    {
        return *(std::get<i>(*couplers_));
    }
    template <std::size_t i>
    auto &model_() const
    {
        return simulator<i>().model();
    }

    template <std::size_t i>
    auto &gridView_() const
    {
        return simulator<i>().gridView();
    }

    JacobianMatrix &jacobian()
    {
        return jacobian_;
    }

    void printMatrixBlock(size_t i, size_t j)
    {
        std::stringstream streamMat;
        Dune::printmatrix(streamMat, jacobian_[i][j], "", "", 10, 2);
        std::cout << streamMat.str();
    }

    void printResidualBlock(size_t i)
    {
        std::cout << residual_[i] << std::endl;
    }

private:
    Simulators *simulators_;
    Couplers *couplers_;
    GlobalEqVector residual_;
    JacobianMatrix jacobian_;

    bool firstIteration_{true};
};

} // namespace Opm

#endif
