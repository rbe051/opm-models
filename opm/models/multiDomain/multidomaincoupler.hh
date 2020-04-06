#ifndef OPM_MULTI_DOMAIN_COUPLER_HH
#define OPM_MULTI_DOMAIN_COUPLER_HH

#include <opm/models/utils/propertysystem.hh>

#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/multiDomain/couplingelementcontext.hh>
#include <opm/models/multiDomain/ecfvcouplingstencil.hh>
#include <opm/models/multiDomain/multidomainproperties.hh>
#include <opm/models/multiDomain/multidomainmapper.hh>
#include <opm/material/densead/Math.hpp>

#include <dune/common/indices.hh>

namespace Opm
{
template <class TypeTag>
class DarcyCoupler;

template <class TypeTag>
class VoidClass;

} // namespace Opm

BEGIN_PROPERTIES
NEW_TYPE_TAG(DarcyCoupler, INHERITS_FROM(MultiDomain));
NEW_PROP_TAG(Coupler);
NEW_PROP_TAG(GridFile);
NEW_PROP_TAG(MappingFile);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(CouplingMapper);
NEW_PROP_TAG(SubTypeTag);
NEW_PROP_TAG(CouplingElementContext);
NEW_PROP_TAG(DomainI);
NEW_PROP_TAG(DomainJ);

SET_TYPE_PROP(DarcyCoupler, Coupler, Opm::DarcyCoupler<TypeTag>);
SET_TYPE_PROP(DarcyCoupler, CouplingMapper, Opm::FaceElementMapper<TypeTag>);
SET_TYPE_PROP(DarcyCoupler, CouplingElementContext, Opm::CouplingElementContext<TypeTag>);
SET_TYPE_PROP(DarcyCoupler, Simulator, Opm::DarcyCoupler<TypeTag>);

SET_INT_PROP(DarcyCoupler, TimeDiscHistorySize, 0);
SET_STRING_PROP(DarcyCoupler, MappingFile, "");

SET_PROP(DarcyCoupler, Stencil)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, CouplingMapper) Mapper;
    typedef typename GET_PROP_TYPE(TypeTag, SubTypeTag) SubTypeTag;

public:
    typedef Opm::EcfvMixedDimStencil<Scalar, GridView, Mapper, SubTypeTag> type;
};

END_PROPERTIES

namespace Opm
{
template <class TypeTag>
class VoidClass
{
};

template <class TypeTag>
class DarcyCoupler
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, DomainI) DomainI;
    typedef typename GET_PROP_TYPE(TypeTag, DomainJ) DomainJ;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) MortarGrid;
    typedef typename GET_PROP_TYPE(TypeTag, Vanguard) Vanguard;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) MortarView;
    typedef typename GET_PROP_TYPE(TypeTag, SubTypeTag) SubTypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, CouplingMapper) CouplingMapper;
    typedef typename GET_PROP_TYPE(TypeTag, CouplingElementContext) CouplingElementContext;

    template <std::size_t i>
    using Simulator = typename SubTypeTag::template Simulator<i>;
    template <std::size_t i>
    using Model = typename SubTypeTag::template Model<i>;
    template <std::size_t i>
    using Stencil = typename SubTypeTag::template Stencil<i>;
    template <std::size_t i>
    using Grid = typename SubTypeTag::template GridView<i>::Grid;
    template <std::size_t i>
    using ElementContext = typename GET_PROP_TYPE(TypeTag, SubTypeTag)::template ElementContext<i>;

    typedef typename GET_PROP_TYPE(typename SubTypeTag::template SubDomain<0>::TypeTag, Evaluation) Evaluation;

    typedef typename SubTypeTag::GlobalEqVector GlobalEqVector;
    typedef typename SubTypeTag::JacobianMatrix JacobianMatrix;

    enum
    {
        numPhases = GET_PROP_VALUE(typename SubTypeTag::template SubDomain<0>::TypeTag, NumPhases)
    };
    enum
    {
        dimWorld = MortarGrid::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Evaluation, dimWorld> EvalDimVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    DarcyCoupler(Simulator<0> &simulator0, Simulator<1> &simulator1)
        : simulator0_{simulator0}, simulator1_{simulator1}
    {
        if (GET_PROP_VALUE(typename SubTypeTag::template SubDomain<0>::TypeTag, NumPhases) !=
            GET_PROP_VALUE(typename SubTypeTag::template SubDomain<1>::TypeTag, NumPhases))
            throw std::runtime_error("Not implemented: Can only couple two models with the same number of faces ");
        if (!(std::is_same<typename GET_PROP_TYPE(typename SubTypeTag::template SubDomain<0>::TypeTag, Evaluation),
                           typename GET_PROP_TYPE(typename SubTypeTag::template SubDomain<1>::TypeTag, Evaluation)>::value))
            throw std::runtime_error("Not implemented: Can only couple two models with the same Evaluation");

        vanguard_.reset(new Vanguard(*this));
        finalizeInit_();
    }

    void volumeFlux(const CouplingElementContext &elemCtx)
    {
        const auto &stencil = elemCtx.stencil(/*timeIdx=*/0);
        const auto &face = stencil.template interiorFace<0>(0);
        auto focusDofIdx = elemCtx.focusDofIndex();
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            flux_[phaseIdx] = 0.0;
            Evaluation pL;
            Evaluation pR;
            if (focusDofIdx == 0)
            {
                pL = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).fluidState().pressure(phaseIdx);
                pR = Opm::getValue(elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0).fluidState().pressure(phaseIdx));
            }
            else if (focusDofIdx == 1)
            {
                pL = Opm::getValue(elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).fluidState().pressure(phaseIdx));
                pR = elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0).fluidState().pressure(phaseIdx);
            }
            else
                DUNE_THROW(Dune::NotImplemented, "Can only couple two degrees of freedom");

            auto deltay = pL - pR;

            const auto &interiorPos = stencil.template subControlVolume<0>(0).globalPos();
            const auto &exteriorPos = stencil.template subControlVolume<1>(0).globalPos();

            Scalar distSquared = 0.0;
            for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            {
                Scalar tmp = exteriorPos[dimIdx] - interiorPos[dimIdx];
                distSquared += tmp * tmp;
            }
            // divide the gradient by the squared distance between the centers of the
            // sub-control volumes: the gradient is the normalized directional vector between
            // the two centers times the ratio of the difference of the values and their
            // distance, i.e., d/abs(d) * delta y / abs(d) = d*delta y / abs(d)^2.
            EvalDimVector quantityGrad;
            for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            {
                Scalar tmp = exteriorPos[dimIdx] - interiorPos[dimIdx];
                quantityGrad[dimIdx] = (tmp / distSquared) * deltay; //deltay * (tmp / distSquared);
            }
            // determine the upstream and downstream DOFs
            const auto &faceNormal = face.normal();
            short upstreamDofIdx;
            short downstreamDofIdx;
            Evaluation tmp = 0.0;
            for (unsigned dimIdx = 0; dimIdx < faceNormal.size(); ++dimIdx)
                tmp += quantityGrad[dimIdx] * faceNormal[dimIdx];

            if (tmp > 0)
            {
                upstreamDofIdx = 0;
                downstreamDofIdx = 1;
            }
            else
            {
                upstreamDofIdx = 1;
                downstreamDofIdx = 0;
            }

            // we only carry the derivatives along if the upstream DOF is the one which
            // we currently focus on
            Evaluation mobility;
            if (upstreamDofIdx == 0)
            {
                if (static_cast<int>(focusDofIdx) == 0)
                    mobility = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).mobility(phaseIdx);
                else
                    mobility = Opm::getValue(elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).mobility(phaseIdx));
            }
            else
            {
                if (static_cast<int>(focusDofIdx) == 1)
                    mobility = elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0).mobility(phaseIdx);
                else
                    mobility = Opm::getValue(elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0).mobility(phaseIdx));
            }

            DimMatrix K;
            const auto K0 = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).intrinsicPermeability();
            const auto K1 = elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0).intrinsicPermeability();

            if (Grid<0>::dimension == Grid<1>::dimension && Grid<0>::dimension != 1)
            {
                // Entry-wise harmonic mean. this is almost certainly wrong if
                // you have off-main diagonal entries in your permeabilities!
                for (unsigned i = 0; i < dimWorld; ++i)
                    for (unsigned j = 0; j < dimWorld; ++j)
                        K[i][j] = Opm::harmonicMean(K0[i][j], K1[i][j]);
            }
            else
            {
                // Harmonic average with permeabilities scaled by element size
                K = K0;
                K *= K1[0][0];
                K *= 1;
                Scalar ndotDist = 0.0;
                for (unsigned dimIdx = 0; dimIdx < faceNormal.size(); ++dimIdx)
                    ndotDist += (exteriorPos[dimIdx] - interiorPos[dimIdx]) * faceNormal[dimIdx];
                K /= K1[0][0] + K0[0][0] * 1e-4 / 2 * ndotDist / distSquared;
            }
            EvalDimVector filterVelocity;
            K.mv(quantityGrad, filterVelocity);

            for (unsigned i = 0; i < faceNormal.size(); ++i)
                flux_[phaseIdx] += mobility * (filterVelocity[i] * faceNormal[i]);

            Scalar alpha = face.area() * elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).extrusionFactor();
            Opm::Valgrind::CheckDefined(alpha);
            assert(alpha > 0.0);
            assert(Opm::isfinite(alpha));
            flux_[phaseIdx] *= alpha;
        }
    }

    void resetSystem_()
    {
        residual_ = 0.0;
        // zero all matrix entries
        jacobian_ = 0.0;
    }

    void advectiveFluxCoupling()
    {
        const auto &idxSet = mortarView().indexSet();
        const auto &gridViewL = simulator0_.gridView();
        const auto &gridViewR = simulator1_.gridView();
        const auto simulators = std::forward_as_tuple(simulator0_, simulator1_);
        CouplingElementContext elemCtx(mortarView(), simulators, *map_);

        for (const auto &e : elements(mortarView()))
        {
            elemCtx.updateAll(e);
            // elemCtx.updateStencil(e);
            // elemCtx.updateAllIntensiveQuantities();
            // elemCtx.updateAllExtensiveQuantities();
            // compute the local residual and its Jacobian
            unsigned numPrimaryDof = elemCtx.numPrimaryDof(/*timeIdx=*/0);
            for (unsigned focusDofIdx = 0; focusDofIdx < numPrimaryDof; ++focusDofIdx)
            {
                elemCtx.setFocusDofIndex(focusDofIdx);
                elemCtx.updateAllExtensiveQuantities();

                // calculate the local residual
                volumeFlux(elemCtx);

                // Global here refers to global within each domain
                unsigned globI = elemCtx.template globalSpaceIndex<0>(/*spaceIdx=*/0, /*timeIdx=*/0);
                unsigned globJ = elemCtx.template globalSpaceIndex<1>(/*spaceIdx=*/0, /*timeIdx=*/0);
                for (unsigned eqIdx = 0; eqIdx < numPhases; ++eqIdx)
                {
                    assert(Opm::isfinite(flux_[eqIdx]));
                    if (focusDofIdx == 0)
                    {
                        residual_[_0][globI][eqIdx] += Opm::getValue(flux_[eqIdx]);
                        for (unsigned pvIdx = 0; pvIdx < numPhases; ++pvIdx)
                        {
                            jacobian_[_0][_0][globI][globI][eqIdx][pvIdx] += flux_[eqIdx].derivative(pvIdx);
                            jacobian_[_1][_0][globJ][globI][eqIdx][pvIdx] -= flux_[eqIdx].derivative(pvIdx);
                        }
                    }
                    else
                    {
                        residual_[_1][globJ][eqIdx] -= Opm::getValue(flux_[eqIdx]);

                        for (unsigned pvIdx = 0; pvIdx < numPhases; ++pvIdx)
                        {
                            jacobian_[_0][_1][globI][globJ][eqIdx][pvIdx] += flux_[eqIdx].derivative(pvIdx);
                            jacobian_[_1][_1][globJ][globJ][eqIdx][pvIdx] -= flux_[eqIdx].derivative(pvIdx);
                        }
                    }
                }
            }
        }
    }

    // Construct the BCRS matrix for the Jacobian of the residual function
    /*!
     * \brief Sets the jacobian build mode
     */
    void
    setJacobianBuildMode()
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
        auto rows_ = model<i>().numTotalDof();
        auto cols_ = model<j>().numTotalDof();

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
        std::vector<NeighborSet> sparsityPattern(model<i>().numTotalDof());

        auto rows_ = model<i>().numTotalDof();
        auto cols_ = model<j>().numTotalDof();

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
        if (0 == i && 1 == j)
            addCouplerPattern_(false, sparsityPattern);
        if (0 == j && 1 == i)
            addCouplerPattern_(true, sparsityPattern);

        return;
    }
    template <class Set>
    void addCouplerPattern_(bool swap, Set &sparsityPattern)
    {
        // for the main model, find out the global indices of the neighboring degrees of
        // freedom of each primary degree of freedom
        typedef std::set<unsigned> NeighborSet;
        int numDofs;
        if (!swap)
            numDofs = model0().numTotalDof();
        else
            numDofs = model1().numTotalDof();

        const auto &idxSet = mortarView().indexSet();
        for (const auto &elem : elements(mortarView()))
        {
            const auto idx = idxSet.index(elem);
            //const auto face = model1().gridView().grid().entity(map_->[i][idx]);
            const auto element1 = map_->template toElement<0>(elem);
            unsigned elIdx1 = model0().gridView().indexSet().index(element1);

            const auto element2 = map_->template toElement<1>(elem);
            unsigned elIdx2 = model1().gridView().indexSet().index(element2);
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
        Stencil<i> stencil(model<i>().gridView(), model<i>().dofMapper());

        auto elemIt = model<i>().gridView().template begin<0>();
        const auto elemEndIt = model<i>().gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            const auto &elem = *elemIt;
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
        size_t numAuxMod = model<i>().numAuxiliaryModules();
        for (unsigned auxModIdx = 0; auxModIdx < numAuxMod; ++auxModIdx)
            model<i>().auxiliaryModule(auxModIdx)->addNeighbors(sparsityPattern);
    }

    void linearize()
    {
        if (firstIteration_)
        {
            initFirstIteration_();
        }
        resetSystem_(); // reset the global linear system of equations.

        // Linearize coupling
        advectiveFluxCoupling();
    }

    void initFirstIteration_()
    {
        firstIteration_ = false;
        setJacobianBuildMode();
        createMatrix();
    }

    JacobianMatrix &jacobian()
    {
        return jacobian_;
    }

    GlobalEqVector &residual()
    {
        return residual_;
    }

    template <std::size_t i>
    void printResidualBlock(Dune::index_constant<i> I)
    {
        std::cout << residual_[I] << std::endl;
    }

    template <std::size_t k, typename std::enable_if_t<(k == 0), int> = 0>
    const auto subDomainIndex() const
    {
        DomainI K;
        return K;
    }
    template <std::size_t k, typename std::enable_if_t<(k == 1), int> = 0>
    const auto subDomainIndex() const
    {
        DomainJ K;
        return K;
    }
    // std::size_t subDomainIndex(std::size_t k)
    // {
    //     if (k == 0)
    //         return domainI_;
    //     else if (k == 1)
    //         return domainJ_;
    //     else
    //         assert(false);
    // }
    template <std::size_t k, typename std::enable_if_t<(k == 0), int> = 0>
    const auto &model() const
    {
        return model0();
    }
    template <std::size_t k, typename std::enable_if_t<(k == 1), int> = 0>
    const auto &model() const
    {
        return model1();
    }

    const auto &model0() const
    {
        return simulator0_.model();
    }
    const auto &model1() const
    {
        return simulator1_.model();
    }
    void finalizeInit_()
    {
        const std::string mappingFileName = EWOMS_GET_PARAM(TypeTag, std::string, MappingFile);

        if (mappingFileName.size() > 0)
            map_.reset(new CouplingMapper(mappingFileName, mortarView(), simulator0_.gridView(), simulator1_.gridView()));
        else
            map_.reset(new CouplingMapper(mortarView(), simulator0_.gridView(), simulator1_.gridView()));
    }

    /*!
     * \brief Return the grid view for which the simulation is done
     */
    const MortarView &mortarView() const
    {
        return vanguard_->gridView();
    }

    static void registerParameters()
    {
        Vanguard::registerParameters();
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, GridGlobalRefinements,
                             "The number of global refinements of the grid "
                             "executed after it was loaded");
        // EWOMS_REGISTER_PARAM(TypeTag, std::string, GridFile,
        //                      "The file name of the file to load");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, MappingFile,
                             "The file name of the mapping file to load");
    }
    //private:
    Dune::index_constant<0> _0;
    Dune::index_constant<1> _1;

    Simulator<0> &simulator0_;
    Simulator<1> &simulator1_;
    std::unique_ptr<CouplingMapper> map_;
    JacobianMatrix jacobian_;
    GlobalEqVector residual_;
    bool verbose_{true};
    std::unique_ptr<Vanguard> vanguard_;
    //std::unique_ptr<Dune::OneDGrid> onedgrid_;
    //std::unique_ptr<MortarView> mortarView_;

    Evaluation flux_[numPhases];

    bool firstIteration_{true};
}; // namespace Opm

} // namespace Opm

#endif
