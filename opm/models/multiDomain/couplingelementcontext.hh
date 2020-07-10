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
 * \copydoc Opm::CouplingElementContext
 */
#ifndef EWOMS_COUPLING_ELEMENT_CONTEXT_HH
#define EWOMS_COUPLING_ELEMENT_CONTEXT_HH

#include <opm/models/utils/alignedallocator.hh>
#include <opm/models/discretization/common/fvbaseelementcontext.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/multiDomain/ecfvcouplingstencil.hh>

#include <opm/material/common/Exceptions.hpp>

#include <dune/common/fvector.hh>

#include <vector>

BEGIN_PROPERTIES
NEW_PROP_TAG(CouplingMapper);
NEW_PROP_TAG(SubTypeTag);
NEW_PROP_TAG(CouplingElementContext);
NEW_PROP_TAG(MortarView);

END_PROPERTIES

namespace Opm
{

/*!
 * \ingroup multiDomain
 *
 * \brief This class stores an array of IntensiveQuantities objects, one
 *        intensive quantities object for each of the two subdomains the mortar
 *        element couples.
 */
template <class TypeTag>
class CouplingElementContext // : public FvBaseElementContext<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, CouplingElementContext) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, CouplingMapper) CouplingMapper;
    typedef typename GET_PROP_TYPE(TypeTag, Stencil) Stencil;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) MortarView;
    template <std::size_t i>
    using Simulator = typename GET_PROP_TYPE(TypeTag, SubTypeTag)::template Simulator<i>;
    template <std::size_t i>
    using ElementContext = typename GET_PROP_TYPE(TypeTag, SubTypeTag)::template ElementContext<i>;
    template <std::size_t i>
    using IntensiveQuantities = typename GET_PROP_TYPE(TypeTag, SubTypeTag)::template IntensiveQuantities<i>;


    static const unsigned dimWorld = MortarView::dimensionworld;
    static const unsigned numEq = GET_PROP_VALUE(TypeTag, NumEq);

    enum
    {
        timeDiscHistorySize = GET_PROP_VALUE(TypeTag, TimeDiscHistorySize)
    };

    typedef typename MortarView::template Codim<0>::Entity Element;
    typedef std::tuple<Simulator<0> &, Simulator<1> &> Simulators;

public:
    /*!
     * \brief The constructor.
     */
    explicit CouplingElementContext(const MortarView &mortarView_, Simulators simulators, const CouplingMapper &dofMapper)
        : mortarView_(mortarView_), stencil_(mortarView_, dofMapper), subElemContext{std::get<0>(simulators), std::get<1>(simulators)}
    {
        // remember the simulator object
        enableStorageCache_ = true;
        stashedDofIdx_ = -1;
        focusDofIdx_ = -1;
    }

    /*!
     * \brief Construct all volume and extensive quantities of an element
     *        from scratch.
     *
     * \param elem The DUNE Codim<0> entity for which the volume
     *             variables ought to be calculated
     */
    void updateAll(const Element &elem)
    {
        asImp_().updateStencil(elem);
        asImp_().updateAllIntensiveQuantities();
        asImp_().updateAllExtensiveQuantities();
    }

    /*!
     * \brief Compute the finite volume geometry for an element.
     *
     * \param elem The grid element for which the finite volume geometry ought to be
     *             computed.
     */
    void updateStencil(const Element &elem)
    {
        // remember the current element
        elemPtr_ = &elem;

        // update the stencil. the center gradients are quite expensive to calculate and
        // most models don't need them, so that we only do this if the model explicitly
        // enables them
        stencil_.update(elem);
        std::get<0>(subElemContext).updateStencil(stencil_.template subElement<0>(0));
        std::get<1>(subElemContext).updateStencil(stencil_.template subElement<1>(0));
    }

    /*!
     * \brief Update the primary topological part of the stencil, but nothing else.
     *
     * \param elem The grid element for which the finite volume geometry ought to be
     *             computed.
     */
    void updatePrimaryStencil(const Element &elem)
    {
        // remember the current element
        elemPtr_ = &elem;

        // update the finite element geometry
        stencil_.updatePrimaryTopology(elem);
    }

    /*!
     * \brief Compute the intensive quantities of all sub-control volumes of the current
     *        element for all time indices.
     */
    void updateAllIntensiveQuantities()
    {
        if (!enableStorageCache_)
        {
            // if the storage cache is disabled, we need to calculate the storage term
            // from scratch, i.e. we need the intensive quantities of all of the history.
            for (unsigned timeIdx = 0; timeIdx < timeDiscHistorySize; ++timeIdx)
                asImp_().updateIntensiveQuantities(timeIdx);
        }
        else
            // if the storage cache is enabled, we only need to recalculate the storage
            // term for the most recent point of history (i.e., for the current iterative
            // solution)
            asImp_().updateIntensiveQuantities(/*timeIdx=*/0);
    }

    /*!
     * \brief Compute the intensive quantities of all sub-control volumes of the current
     *        element for a single time index.
     *
     * \param timeIdx The index of the solution vector used by the time discretization.
     */
    void updateIntensiveQuantities(unsigned timeIdx)
    {
        std::get<0>(subElemContext).updateIntensiveQuantities(timeIdx);
        std::get<1>(subElemContext).updateIntensiveQuantities(timeIdx);
    }

    void updatePrimaryIntensiveQuantities(unsigned timeIdx)
    {
        std::get<0>(subElemContext).updatePrimaryIntensiveQuantities(timeIdx);
        std::get<1>(subElemContext).updatePrimaryIntensiveQuantities(timeIdx);
    }

/*!
     * \brief Compute the extensive quantities of all sub-control volume
     *        faces of the current element for all time indices.
     */
    void updateAllExtensiveQuantities()
    { asImp_().updateExtensiveQuantities(/*timeIdx=*/0); }

    /*!
     * \brief Compute the extensive quantities of all sub-control volume
     *        faces of the current element for a single time index.
     *
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    void updateExtensiveQuantities(unsigned timeIdx)
    {
        std::get<0>(subElemContext).updateExtensiveQuantities(timeIdx);
        std::get<1>(subElemContext).updateExtensiveQuantities(timeIdx);
    }

    /*!
     * \brief Compute the intensive quantities of a single sub-control volume of the
     *        current element for a single time index.
     *
     * \param priVars The PrimaryVariables which should be used to calculate the
     *                intensive quantities.
     * \param dofIdx The local index in the current element of the sub-control volume
     *               which should be updated.
     * \param timeIdx The index of the solution vector used by the time discretization.
     */
    template<class PrimaryVariables>
    void updateIntensiveQuantities(const PrimaryVariables &priVars, unsigned dofIdx, unsigned timeIdx)
    {
        throw std::logic_error("NotImplemented: updateIntensiveQuantities() is not implemented for coupling element context");
    }

    /*!
     * \brief Return a reference to the intensive quantities of a
     *        sub-control volume of a subdomain at a given time.
     *
     * If the time step index is not given, return the volume
     * variables for the current time. Template parameter defines which subdomain
     * intensive quantities to be returned.
     *
     * \param dofIdx The local index of the degree of freedom in the current element.
     * \param timeIdx The index of the solution vector used by the time discretization.
     */
    template<std::size_t i>
    const IntensiveQuantities<i>& intensiveQuantities(unsigned dofIdx, unsigned timeIdx) const
    {
        const auto& subCtx = std::get<i>(subElemContext);
#ifndef NDEBUG
        assert(0 <= dofIdx && dofIdx < numDof(timeIdx));

        if (enableStorageCache_ && timeIdx != 0 && subCtx.problem().recycleFirstIterationStorage())
            throw std::logic_error("If caching of the storage term is enabled, only the intensive quantities "
                                   "for the most-recent substep (i.e. time index 0) are available!");
#endif

        return subCtx.intensiveQuantities(dofIdx, timeIdx);
    }

    /*!
     * \brief Sets the degree of freedom on which the simulator is currently "focused" on
     *
     * I.e., in the case of automatic differentiation, all derivatives are with regard to
     * the primary variables of that degree of freedom. Only "primary" DOFs can be
     * focused on.
     */
    void setFocusDofIndex(unsigned dofIdx)
    {
        focusDofIdx_ = dofIdx;
    }
    
    /*!
     * \brief Returns the degree of freedom on which the simulator is currently "focused" on
     *
     * \copydetails setFocusDof()
     */
    unsigned focusDofIndex() const
    { return focusDofIdx_; }

    /*!
     * \brief Return the current finite element geometry.
     *
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    const Stencil &stencil(unsigned timeIdx OPM_UNUSED) const
    {
        return stencil_;
    }

    /*!
     * \brief Return the number of sub-control volumes of the current element.
     */
    size_t numDof(unsigned timeIdx) const
    {
        return stencil(timeIdx).numDof();
    }

    /*!
     * \brief Return a reference to the problem of a subdomain.
     */
    template<std::size_t i>
    const auto& problem() const
    { return std::get<i>(subElemContext).problem();}

    /*!
     * \brief Return the global spatial index for a sub-control volume of a subdomain
     * 
     * The template parameter defines which subdomain the global index referes to.
     *
     * \param dofIdx The local index of the degree of freedom
     *               in the current element.
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    template<std::size_t i>
    unsigned globalSpaceIndex(unsigned dofIdx, unsigned timeIdx) const
    {
        return std::get<i>(subElemContext).globalSpaceIndex(dofIdx, timeIdx);
    }

    /*!
     * \brief Return the number of primary degrees of freedom of the current element.
     */
    size_t numPrimaryDof(unsigned timeIdx) const
    {
        return stencil(timeIdx).numPrimaryDof();
    }

private:
    Implementation &asImp_()
    {
        return *static_cast<Implementation *>(this);
    }

    const Implementation &asImp_() const
    {
        return *static_cast<const Implementation *>(this);
    }

protected:
    const Element *elemPtr_;
    const MortarView mortarView_;
    Stencil stencil_;
    std::tuple<ElementContext<0>, ElementContext<1>> subElemContext;

    int stashedDofIdx_;
    int focusDofIdx_;
    bool enableStorageCache_;
};

} // namespace Opm

#endif
