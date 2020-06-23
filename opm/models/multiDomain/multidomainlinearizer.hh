// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::MultiDomainLinearizer
 */

#ifndef OPM_MULTI_DOMAIN_LINEARIZER_HH
#define OPM_MULTI_DOMAIN_LINEARIZER_HH

#include <opm/models/utils/propertysystem.hh>

#include <opm/material/densead/Math.hpp>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/multiDomain/couplingelementcontext.hh>
#include <opm/models/multiDomain/multidomainproperties.hh>

#include <dune/common/indices.hh>
#include <dune/grid/common/intersection.hh>

BEGIN_PROPERTIES
NEW_PROP_TAG(CouplerTypeTag);
END_PROPERTIES

namespace Opm {
/*!
 * \ingroup MultiDomainModel
 *
 * \brief Global linearizer that linearize the multidomain system
 * 
 * A multidomain model consist of several subdomains which are initiated as
 * standard OPM models. The global system is linearized by calling the
 * local linearizers for each subdomain, and the resulting matrix is assembled
 * into the global system based on the domain id.
 * 
*/
template <class TypeTag>
class MultiDomainLinearizer {
    using Simulator = typename GET_PROP_TYPE(TypeTag, Simulator);
    template <std::size_t i>
    using Coupler = typename GET_PROP_TYPE(TypeTag, CouplerTypeTag)::template Coupler<i>;
    template <std::size_t i>
    using SubSimulator = typename GET_PROP_TYPE(TypeTag, SubTypeTag)::template Simulator<i>;
    template <std::size_t i>
    using Stencil = typename GET_PROP_TYPE(TypeTag, SubTypeTag)::template Stencil<i>;
    template <std::size_t i>
    using GridView = typename GET_PROP_TYPE(TypeTag, SubTypeTag)::template GridView<i>;
    typedef typename GET_PROP_TYPE(TypeTag, SubTypeTag)::JacobianMatrix JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, SubTypeTag)::GlobalEqVector GlobalEqVector;

    template <std::size_t i>
    using Element = typename GridView<i>::template Codim<0>::Entity;
    using Simulators = typename GET_PROP_TYPE(TypeTag, SubTypeTag)::template TupleOfSharedPtr<SubSimulator>;
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
     */
    void init(Simulator& simulator, Couplers& couplers)
    {
        simulator_ = &simulator;
        couplers_ = &couplers;
    }
    /*!
     * \brief Register all run-time parameters for the multidomain linearizer.
     */
    static void registerParameters()
    {
    }

    /*!
     * \brief Causes the Jacobian matrix to be recreated from scratch before the next
     *        iteration.
     *
     * This method is usally called if the sparsity pattern has changed for some
     * reason. (e.g. by modifications of the grid or changes of the auxiliary equations.)
     */
    void eraseMatrix()
    {
        throw std::logic_error("Not implemented:eraseMatrix() method not implemented by the multidomainLinearizer");
    }

    /*!
     * \brief Linearize the full system of non-linear equations.
     *
     * This means the spatial subdomains plus all auxiliary equations.
     */
    void linearize()
    {
        linearizeSubModels();
    }
    /*!
     * \brief Linearize the full system of non-linear equations.
     *
     * That means that the global Jacobian of the residual is assembled and the residual
     * is evaluated for the current solution for all of the sub-domains.
     *
     * The current state of affairs (esp. the previous and the current solutions) is
     * represented by each of the submodel objects.
     */
    void linearizeSubModels()
    {
        // we defer the initialization of the Jacobian matrix until here because the
        // auxiliary modules usually assume the problem, model and grid to be fully
        // initialized...
        if (firstIteration_) {
            initFirstIteration_();
        }
        resetSystem_(); // reset the global linear system of equations.

        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(jacobian_)), [&](const auto I) {
            simulator_->template model<I>().linearizer().linearizeDomain();
            simulator_->template model<I>().linearizer().linearizeAuxiliaryEquations();
            simulator_->template model<I>().linearizer().finalize();
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

    /*!
     * \brief Return reference to global Jacobian matrix backend.
     */
    JacobianMatrix& jacobian()
    {
        return jacobian_;
    }

    /*!
     * \brief Return reference to global residual vector.
     */
    GlobalEqVector& residual()
    {
        return residual_;
    }

private:
    template <std::size_t i>
    Coupler<i>& coupler_() const
    {
        return *(std::get<i>(*couplers_));
    }
    template <std::size_t i>
    auto& model_() const
    {
        return simulator_->template model<i>();
    }

    template <std::size_t i>
    auto& gridView_() const
    {
        return simulator_->template gridView<i>();
    }

    void initFirstIteration_()
    {
        firstIteration_ = false;
        // make sure all matrices have the same build mode
        setJacobianBuildMode_();
        // initialize the matices
        createMatrix_();
    }

    // Add the subdomain jacobian to the global jacobian
    template <std::size_t i, std::size_t j, class Jac>
    void addToJacobian(Dune::index_constant<i> I, Dune::index_constant<j> J, const Jac& localJac)
    {
        Dune::index_constant<0> _0;
        Dune::index_constant<1> _1;
        jacobian_[I][I] += localJac[_0][_0];
        jacobian_[I][J] += localJac[_0][_1];
        jacobian_[J][I] += localJac[_1][_0];
        jacobian_[J][J] += localJac[_1][_1];
    }

    // Add the subdomain residual to the global residual
    template <std::size_t i, std::size_t j, class Res>
    void addToResidual(Dune::index_constant<i> I, Dune::index_constant<j> J, const Res& localRes)
    {
        Dune::index_constant<0> _0;
        Dune::index_constant<1> _1;
        residual_[I] += localRes[_0];
        residual_[J] += localRes[_1];
    }

    // reset the global linear system of equations.
    void resetSystem_()
    {
        residual_ = 0.0;
        // zero all matrix entries
        jacobian_ = 0.0;
    }

    // Construct the BCRS matrix for the Jacobian of the residual function.
    // Sets all matrices to random build mode.
    // Throws a Dune::NotImplemented exception if the submatrices are not BCRS.
    void setJacobianBuildMode_()
    {
        using namespace Dune::Hybrid;
        forEach(jacobian_, [](auto& jacRow) {
            forEach(jacRow, [](auto& jacBlock) {
                using BlockType = std::decay_t<decltype(jacBlock)>;
                if (jacBlock.buildMode() == BlockType::BuildMode::unknown)
                    jacBlock.setBuildMode(BlockType::BuildMode::random);
                else if (jacBlock.buildMode() != BlockType::BuildMode::random)
                    DUNE_THROW(Dune::NotImplemented, "Only BCRS matrices with random build mode are supported at the moment");
            });
        });
    }

    // Reserve the sparcity pattern for a local matrix.
    template <std::size_t i, std::size_t j, class Set>
    void reserve_(const std::vector<Set>& sparsityPattern)
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
        for (size_t dofIdx = 0; dofIdx < rows_; ++dofIdx) {
            auto nIt = sparsityPattern[dofIdx].begin();
            auto nEndIt = sparsityPattern[dofIdx].end();
            for (; nIt != nEndIt; ++nIt)
                jacobian_[I][J].addindex(dofIdx, *nIt);
        }
        jacobian_[I][J].endindices();
    }

    // Initialize all local matrices
    void createMatrix_()
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(jacobian_)), [&](const auto domainI) {
            forEach(integralRange(Dune::Hybrid::size(jacobian_[domainI])), [&](const auto domainJ) {
                setJacobianPattern_(domainI, domainJ);
            });
        });
    }

    // Find and set the sparcity pattern of all matrices
    template <std::size_t i, std::size_t j>
    void setJacobianPattern_(Dune::index_constant<i> I, Dune::index_constant<j> J)
    {
        typedef std::set<unsigned> NeighborSet;
        std::vector<NeighborSet> sparsityPattern(model_<i>().numTotalDof());

        auto rows_ = model_<i>().numTotalDof();
        auto cols_ = model_<j>().numTotalDof();

        jacobian_[I][J].setSize(rows_, cols_);
        residual_[I].resize(rows_);
        if (i == j) {
            addDiagonalPattern_<i>(sparsityPattern);
        } else {
            addOffDiagonalPattern_<i, j>(sparsityPattern);
        }

        // allocate raw matrix

        // create matrix structure based on sparsity pattern
        reserve_<i, j>(sparsityPattern);
    }

    // add the sparcity pattern of the off diagonal matrices
    template <std::size_t i, std::size_t j, class Set>
    void addOffDiagonalPattern_(std::vector<Set>& sparsityPattern)
    {
        // chekc if there is a coupling between domain i and j
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

    // add the sparcity pattern of two coupled domains.
    // TODO: this assumes that the submodels and coupler uses an
    // ecfv discretization (tpfa).
    template <std::size_t k, class Set>
    void addCouplerPattern_(bool swap, Set& sparsityPattern)
    {
        const auto& model0 = coupler_<k>().template model<0>();
        const auto& model1 = coupler_<k>().template model<1>();
        const auto& mortarView = coupler_<k>().mortarView();
        const auto& map = coupler_<k>().projectionMapper();
        // find out the global indices of the neighboring degrees of
        // freedom of each primary degree of freedom
        typedef std::set<unsigned> NeighborSet;
        int numDofs;
        if (!swap)
            numDofs = model0.numTotalDof();
        else
            numDofs = model1.numTotalDof();

        const auto& idxSet = mortarView.indexSet();
        for (const auto& elem : elements(mortarView)) {
            const auto idx = idxSet.index(elem);

            const auto element0 = map.template toElement<0>(elem);
            unsigned elIdx0 = model0.gridView().indexSet().index(element0);

            const auto element1 = map.template toElement<1>(elem);
            unsigned elIdx1 = model1.gridView().indexSet().index(element1);
            if (!swap)
                sparsityPattern[elIdx0].insert(elIdx1);
            else // swap if model 1 is the main model
                sparsityPattern[elIdx1].insert(elIdx0);
        }
    }

    // Add the sparcity pattern of the diagonal matrices.
    // TODO: This assumes an ecfv discretization (tpfa). The sparcity pattern
    // should be taken from each submodel. Hoverer, the domain coupler also
    // only works for ecfv discretization, so keep this for now.
    template <std::size_t i, class Set>
    void addDiagonalPattern_(std::vector<Set>& sparsityPattern) const
    {
        // for the main model, find out the global indices of the neighboring degrees of
        // freedom of each primary degree of freedom
        Stencil<i> stencil(gridView_<i>(), model_<i>().dofMapper());

        auto elemIt = gridView_<i>().template begin<0>();
        const auto elemEndIt = gridView_<i>().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element<i>& elem = *elemIt;
            stencil.update(elem);

            for (unsigned primaryDofIdx = 0; primaryDofIdx < stencil.numPrimaryDof(); ++primaryDofIdx) {
                unsigned myIdx = stencil.globalSpaceIndex(primaryDofIdx);

                for (unsigned dofIdx = 0; dofIdx < stencil.numDof(); ++dofIdx) {
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

    Simulator* simulator_;
    Couplers* couplers_;
    GlobalEqVector residual_;
    JacobianMatrix jacobian_;

    bool firstIteration_ { true };
};

} // namespace Opm

#endif
