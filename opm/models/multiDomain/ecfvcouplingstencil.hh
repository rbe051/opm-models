
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
 */
#ifndef EWOMS_ECFV_COUPLING_STENCIL_HH
#define EWOMS_ECFV_COUPLING_STENCIL_HH

#include <opm/material/common/ConditionalStorage.hpp>

#include <dune/common/version.hh>
#include <vector>

namespace Opm {
/*!
 * \ingroup multiDomain
 *
 * \brief Represents the stencil (finite volume geometry) of a single
 *        mortar element in the ECFV discretization.
 *
 * This stencil couple two domains of the same dimension. For a mixed
 * dimensional coupling, see EcfvMixedDimStencil.
 * 
 * The ECFV discretization is a element centered finite volume
 * approach. This means that each element corresponds to a control
 * volume. It is assumed that the discretization stencil of the
 * two subdomains are also ecfv.
 */
template <class Scalar,
    class GridView,
    class Mapper,
    class SubTypeTag,
    bool needFaceIntegrationPos = true,
    bool needFaceNormal = true>
class EcfvCouplingStencil {
    enum {
        dimWorld = GridView::dimensionworld
    };
    typedef typename GridView::template Codim<0>::Entity MortarElement;
    template <std::size_t i>
    using SubGridView = typename SubTypeTag::template GridView<i>;
    template <std::size_t i>
    using CoordScalar = typename SubGridView<i>::ctype;
    template <std::size_t i>
    using SubElement = typename SubGridView<i>::template Codim<0>::Entity;
    template <std::size_t i>
    using Intersection = typename SubGridView<i>::Intersection;
    template <std::size_t i>
    using GlobalPosition = typename Dune::FieldVector<CoordScalar<i>, dimWorld>;
    template <std::size_t i>
    using WorldVector = typename Dune::FieldVector<Scalar, dimWorld>;

public:
    /*!
     * \brief Represents a sub-control volume.
     *
     * For element centered finite volumes, this is equivalent to the
     * element
     */
    template <int i>
    class SubControlVolume {
    public:
        // default construct an uninitialized object.
        // this is only here because std::vector needs it...
        SubControlVolume()
        {
        }

        SubControlVolume(const SubElement<i>& element)
            : element_(element)
        {
            update();
        }

        void update(const SubElement<i>& element)
        {
            element_ = element;
        }

        void update()
        {
            const auto& geometry = element_.geometry();
            centerPos_ = geometry.center();
            volume_ = geometry.volume();
        }

        /*!
         * \brief The global position associated with the sub-control volume
         */
        const GlobalPosition<i>& globalPos() const
        {
            return centerPos_;
        }

        /*!
         * \brief The center of the sub-control volume
         */
        const GlobalPosition<i>& center() const
        {
            return centerPos_;
        }

        /*!
         * \brief The volume [m^3] occupied by the sub-control volume
         */
        Scalar volume() const
        {
            return volume_;
        }

    private:
        GlobalPosition<i> centerPos_;
        Scalar volume_;
        SubElement<i> element_;
    };

    /*!
     * \brief Represents a sub-control volume Face.
     *
     * For element centered finite volumes, this is equivalent to the
     * intersection
     */
    template <bool needIntegrationPos, bool needNormal, int i>
    class CouplingSubControlVolumeFace {
    public:
        CouplingSubControlVolumeFace()
        {
        }

        CouplingSubControlVolumeFace(const Intersection<i>& intersection, unsigned localNeighborIdx)
        {
            exteriorIdx_ = static_cast<unsigned short>(localNeighborIdx);

            if (needNormal)
                (*normal_) = intersection.centerUnitOuterNormal();

            const auto& geometry = intersection.geometry();
            if (needIntegrationPos)
                (*integrationPos_) = geometry.center();
            area_ = geometry.volume();
        }

        /*!
         * \brief Returns the local index of the degree of freedom to
         *        the face's interior.
         */
        unsigned short interiorIndex() const
        {
            // The local index of the control volume in the interior
            // of a face of the stencil in the element centered finite
            // volume discretization is always the "central"
            // element. In this implementation this element always has
            // index 0....
            return 0;
        }

        /*!
         * \brief Returns the local index of the degree of freedom to
         *        the face's outside.
         */
        unsigned short exteriorIndex() const
        {
            return exteriorIdx_;
        }

        /*!
         * \brief Returns the global position of the face's
         *        integration point.
         */
        const GlobalPosition<i>& integrationPos() const
        {
            return *integrationPos_;
        }

        /*!
         * \brief Returns the outer unit normal at the face's
         *        integration point.
         */
        const WorldVector<i>& normal() const
        {
            return *normal_;
        }

        /*!
         * \brief Returns the area [m^2] of the face
         */
        Scalar area() const
        {
            return area_;
        }

    private:
        Opm::ConditionalStorage<needIntegrationPos, GlobalPosition<i>> integrationPos_;
        Opm::ConditionalStorage<needNormal, WorldVector<i>> normal_;
        Scalar area_;

        unsigned short exteriorIdx_;
    };

    template <std::size_t i>
    using SubControlVolumeFace = CouplingSubControlVolumeFace<needFaceIntegrationPos, needFaceNormal, i>;

    EcfvCouplingStencil(const GridView& gridView, const Mapper& mapper)
        : gridView_(gridView)
        , elementMapper_(mapper)
    {
        // try to ensure that the mapper passed indeed maps elements
        assert(int(gridView.size(/*codim=*/0)) == int(elementMapper_.size()));
    }

    /*!
     * \brief update the topology from a mortar element
     * 
     * The topology consist of the two subdomain elements that the mortar element couples.
     */
    void updateTopology(const MortarElement& element)
    {
        mortarElement_ = element;

        const auto& intersection0 = elementMapper_.template toIntersection<0>(element);
        std::get<0>(elements_) = intersection0.inside();
        std::get<0>(subControlVolumes_).reset(new SubControlVolume<0>(intersection0.inside()));
        std::get<0>(couplingFaces_).reset(new SubControlVolumeFace<0>(intersection0, 0));

        const auto& intersection1 = elementMapper_.template toIntersection<1>(element);
        std::get<1>(elements_) = intersection1.inside();
        std::get<1>(subControlVolumes_).reset(new SubControlVolume<1>(intersection1.inside()));
        std::get<1>(couplingFaces_).reset(new SubControlVolumeFace<1>(intersection1, 1));
    }

    void update(const MortarElement& element)
    {
        updateTopology(element);
    }

    void updatePrimaryTopology(const MortarElement& element)
    {
        const auto& intersection = elementMapper_.template toIntersection<0>(element);
        std::get<0>(elements_) = intersection.inside();
        std::get<0>(subControlVolumes_).reset(new SubControlVolume<0>(subControlVolumes_.inside()));
    }

    /*!
     * \brief Returns the number of degrees of freedom which the
     *        current element interacts with.
     * 
     * We assume that the mortar stencil couples exactly two dofs
     */
    size_t numDof() const
    {
        return 2;
    }

    /*!
     * \brief Returns the face object belonging to a given face index
     *        in the interior of the domain.
     */
    template <std::size_t i>
    const SubControlVolumeFace<i>& interiorFace(unsigned faceIdx) const
    {
        return *(std::get<i>(couplingFaces_));
    }
    /*!
     * \brief Returns the sub-control volume object belonging to a
     *        given degree of freedom.
     */
    template <std::size_t i>
    const SubControlVolume<i>& subControlVolume(unsigned dofIdx) const
    {
        return *(std::get<i>(subControlVolumes_));
    }

    /*!
     * \brief Return the element given the index of a degree of
     *        freedom.
     */
    const MortarElement& element(unsigned dofIdx) const
    {
        assert(0 == dofIdx);

        return mortarElement_;
    }
    /*!
     * \brief Return the element given the index of a degree of
     *        freedom.
     */
    template <std::size_t i>
    const SubElement<i>& subElement(unsigned dofIdx) const
    {
        assert(0 <= dofIdx && dofIdx < numDof());

        return std::get<i>(elements_);
    }
    /*!
     * \brief Returns the number of degrees of freedom which are contained
     *        by within the current element.
     *
     * Primary DOFs are always expected to have a lower index than
     * "secondary" DOFs.
     *
     * For element centered finite elements, this is only the two central DOFs.
     */
    size_t numPrimaryDof() const
    {
        return 2;
    }

protected:
    const GridView& gridView_;
    const Mapper& elementMapper_;

    MortarElement mortarElement_;
    std::tuple<SubElement<0>, SubElement<1>> elements_;
    std::tuple<std::unique_ptr<SubControlVolumeFace<0>>, std::unique_ptr<SubControlVolumeFace<1>>> couplingFaces_;
    std::tuple<std::unique_ptr<SubControlVolume<0>>, std::unique_ptr<SubControlVolume<1>>> subControlVolumes_;
};

/*!
 * \ingroup multiDomain
 *
 * \brief Represents the stencil (finite volume geometry) of a single
 *        mortar element in the ECFV discretization.
 *
 * This stencil couple two domains of different dimensions. The first
 * domain is assumed to have a dimension exactly one more than the 
 * mortar grid and second domain.
 * 
 * The ECFV discretization is a element centered finite volume
 * approach. This means that each element corresponds to a control
 * volume. It is assumed that the discretization stencil of the
 * two subdomains are also ecfv.
 */
template <class Scalar,
    class GridView,
    class Mapper,
    class SubTypeTag,
    bool needFaceIntegrationPos = true,
    bool needFaceNormal = true>
class EcfvMixedDimStencil : public EcfvCouplingStencil<Scalar, GridView, Mapper, SubTypeTag, needFaceIntegrationPos, needFaceNormal> {
    typedef EcfvCouplingStencil<Scalar, GridView, Mapper, SubTypeTag, needFaceIntegrationPos, needFaceNormal> ParentType;
    typedef typename GridView::template Codim<0>::Entity MortarElement;

public:
    using ParentType::ParentType;

    /*!
     * \brief update the topology from a mortar element
     * 
     * The topology consist of the two subdomain elements that the mortar element couples.
     */
    void updateTopology(const MortarElement& element)
    {
        this->mortarElement_ = element;

        const auto& intersection0 = this->elementMapper_.template toIntersection<0>(element);
        std::get<0>(this->elements_) = intersection0.inside();
        std::get<0>(this->subControlVolumes_).reset(new typename ParentType::template SubControlVolume<0>(intersection0.inside()));
        std::get<0>(this->couplingFaces_).reset(new typename ParentType::template SubControlVolumeFace<0>(intersection0, 0));

        const auto& element1 = this->elementMapper_.template toElement<1>(element);
        std::get<1>(this->elements_) = element1;
        std::get<1>(this->subControlVolumes_).reset(new typename ParentType::template SubControlVolume<1>(element1));
    }
    void update(const MortarElement& element)
    {
        updateTopology(element);
    }
};
} // namespace Opm

#endif
