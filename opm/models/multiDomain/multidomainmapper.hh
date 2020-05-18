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
 * \ingroup MultiDomain
 *
 * \brief Class that handles the mappings between domains
 */

namespace Opm {

/*!
 * \ingroup MultiDomainModel
 *
 * \brief Virtual class object that defines the base for the multidomain mappers
 *
 * This class defines the functions that must be overloaded by the domain
 * mappers
 *
 */
template <class TypeTag>
class BaseMapper {
    typedef typename GET_PROP_TYPE(TypeTag, CouplingMapper) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) MortarView;
    template <std::size_t i>
    using GridView = typename GET_PROP_TYPE(TypeTag,
        SubTypeTag)::template GridView<i>;
    template <std::size_t i>
    using Intersection = typename GridView<i>::Intersection;
    typedef typename MortarView::template Codim<0>::Entity MortarElement;
    template <std::size_t i>
    using EntitySeed = typename GridView<i>::Grid::template Codim<0>::EntitySeed;

    enum {
        dimWorld = MortarView::dimensionworld
    };

public:
    /*!
   * \brief The constructor. Read mapping from file.
   */
    BaseMapper(const std::string& file_name,
        const MortarView& mortarView,
        const GridView<0>& gridView0,
        const GridView<1>& gridView1)
        : mortarView_ { mortarView }
        , gridView0_ { gridView0 }
        , gridView1_ { gridView1 }
    {
    }

    /*!
   * \brief The constructor. Read find mapping by comparing geometry.
   */
    BaseMapper(const MortarView& mortarView,
        const GridView<0>& gridView1,
        const GridView<1>& gridView2)
        : mortarView_ { mortarView }
        , gridView0_ { gridView1 }
        , gridView1_ { gridView2 }
    {
    }

    /*!
   * \brief Returns the grid intersection that corresponds to the mortar element
   */
    template <std::size_t i>
    const Intersection<i> toIntersection(MortarElement e) const
    {
        throw std::logic_error(
            "Not implemented:toIntersection() method not implemented by the actual "
            "mapper");
    }

    /*!
   * \brief Returns the grid element that correspond to the the mortar element
   */
    template <std::size_t i>
    const typename GridView<i>::template Codim<0>::Entity toElement(
        MortarElement e) const
    {
        throw std::logic_error(
            "Not implemented:toElement() method not implemented by the actual "
            "mapper");
    }

    /*!
   * \brief Returns the number of mortar elements
   */
    unsigned size() const { return mortarView_.size(0); }

    /*!
   * \brief Returns the number of subdomain elements
   */
    unsigned size(int i) const
    {
        switch (i) {
        case 0:
            return gridView0_.size(0);
        case 1:
            return gridView1_.size(0);
        default:
            throw std::logic_error(
                "Logic error:size() method can only find the size of subdomain 0 "
                "or 1");
        }
    }

    /*!
   * \brief Read the projections from subgrids to mortar elements from file.
   */
    void setMapFromFile(const std::string& file_name)
    {
        throw std::logic_error(
            "Not implemented:setMapFromFile() method not implemented by the actual "
            "mapper");
    }

    /*!
   * \brief Calculate the projections from subgrids to mortar elements.
   */
    void calculateMap()
    {
        throw std::logic_error(
            "Not implemented:calculateMap() method not implemented by the actual "
            "mapper");
    }

protected:
    // Find the intersections of the gridView that matches the mortar elements
    template <std::size_t i>
    void calculateFaceMap_(const GridView<i>& gridView)
    {
        const auto& idxSetMortar = this->mortarView_.indexSet();
        const auto& idxSetGrid = gridView.indexSet();
        auto& map = asImp_().projectionMap();

        for (const auto& mortarEle : elements(this->mortarView_)) {
            const auto& mc = mortarEle.geometry().center();

            bool foundMap = false;
            for (const auto& ele : elements(gridView)) {
                for (const auto& intersection : intersections(gridView, ele)) {
                    const auto& fc = intersection.geometry().center();
                    bool equal = true;
                    for (int dim = 0; dim < dimWorld; ++dim) {
                        equal = equal && (std::abs(fc[dim] - mc[dim]) < 1e-10);
                    }
                    if (equal) {
                        std::get<i>(map)[idxSetMortar.index(mortarEle)] = intersection;
                        foundMap = true;
                        break;
                    }
                }
            }
            if (!foundMap)
                throw std::runtime_error(
                    "Could not find map. Are you sure the grids match?");
        }
    }
    //! get an element from an index i
    template <std::size_t I>
    const auto element(std::size_t i) const
    {
        const auto& seed = std::get<I>(seeds_)[i];
        const auto entity = gridView1_.grid().entity(seed);
        return entity;
    }
    //! get an element from an index i
    template <std::size_t I>
    const auto intersection(std::size_t i) const
    {
        return std::get<I>(intersection_)[i];
    }

    // Find the elements of the gridView that matches the mortar elements. The
    // gridView should have the same dimension as the mortar grid.
    template <std::size_t i>
    void calculateElementMap_(const GridView<i>& gridView)
    {
        const auto& idxSetMortar = this->mortarView_.indexSet();
        const auto& idxSetGrid = gridView.indexSet();
        auto& map = asImp_().projectionMap();

        for (const auto& mortarEle : elements(this->mortarView_)) {
            const auto& mc = mortarEle.geometry().center();
            bool foundMap = false;
            for (const auto& ele : elements(gridView)) {
                const auto& fc = ele.geometry().center();
                bool equal = true;
                for (int dim = 0; dim < dimWorld; ++dim) {
                    equal = equal && (std::abs(fc[dim] - mc[dim]) < 1e-10);
                }
                if (equal) {
                    // map[idxSetMortar.index(e)] = idxSetGrid.index(facet);
                    std::get<i>(map)[idxSetMortar.index(mortarEle)] = ele;
                    foundMap = true;
                    break;
                }
            }
            if (!foundMap) {
                throw std::runtime_error(
                    "Could not find map. Are you sure the grids match?");
            }
        }
    }

    void createElementToIndexMapper()
    {
        const auto& indxSet0 = gridView0_.indexSet();
        std::get<0>(seeds_).resize(gridView0_.size(0));
        for (const auto& ele : elements(gridView0_))
            std::get<0>(seeds_)[indxSet0.index(ele)] = ele.seed();

        const auto& indxSet1 = gridView1_.indexSet();
        std::get<1>(seeds_).resize(gridView1_.size(0));
        for (const auto& ele : elements(gridView1_))
            std::get<1>(seeds_)[indxSet1.index(ele)] = ele.seed();

        std::get<0>(intersection_).resize(gridView0_.size(1));
        for (const auto& ele : elements(gridView0_)) {
            for (const auto& intersection : intersections(gridView0_, ele)) {
                const int subIdx = intersection.indexInInside();
                const int faceIdx = indxSet0.subIndex(ele, subIdx, 1);
                std::get<0>(intersection_)[faceIdx] = intersection;
            }
        }
        std::get<1>(intersection_).resize(gridView1_.size(1));
        for (const auto& ele : elements(gridView1_)) {
            for (const auto& intersection : intersections(gridView1_, ele)) {
                const int subIdx = intersection.indexInInside();
                const int faceIdx = indxSet1.subIndex(ele, subIdx, 1);
                std::get<1>(intersection_)[faceIdx] = intersection;
            }
        }
    }

    // This method should be called after the projection map has been allocated
    void finalizeInit_(std::string file_name)
    {
        createElementToIndexMapper();
        std::get<0>(asImp_().projectionMap()).resize(mortarView_.size(0));
        std::get<1>(asImp_().projectionMap()).resize(mortarView_.size(0));
        asImp_().setMapFromFile(file_name);
    }

    // This method should be called after the projection map has been allocated
    void finalizeInit_()
    {
        std::get<0>(asImp_().projectionMap()).resize(mortarView_.size(0));
        std::get<1>(asImp_().projectionMap()).resize(mortarView_.size(0));
        asImp_().calculateMap();
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation& asImp_() { return *static_cast<Implementation*>(this); }
    //! \copydoc asImp_()
    const Implementation& asImp_() const
    {
        return *static_cast<const Implementation*>(this);
    }

protected:
    const MortarView& mortarView_;
    const GridView<0>& gridView0_;
    const GridView<1>& gridView1_;

    std::tuple<std::vector<EntitySeed<0>>, std::vector<EntitySeed<1>>> seeds_;
    std::tuple<std::vector<Intersection<0>>, std::vector<Intersection<1>>>
        intersection_;
};

/*!
 * \ingroup MultiDomainModel
 *
 * \brief Class object to identify face to face projections between domains.
 *
 * This class finds and stores the projections from the faces of the subdomains
 * to the cells of the mortar grid.
 * ------- | -------
 * |grid0| | |grid1 |
 * |_____| | |______|
 *         ^
 *     mortarGrid
 *
 */
template <class TypeTag>
class FaceFaceMapper : public BaseMapper<TypeTag> {
    typedef BaseMapper<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, CouplingMapper) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) MortarView;
    template <std::size_t i>
    using GridView = typename GET_PROP_TYPE(TypeTag,
        SubTypeTag)::template GridView<i>;
    template <std::size_t i>
    using Intersection = typename GridView<i>::Intersection;
    typedef std::tuple<std::vector<Intersection<0>>, std::vector<Intersection<1>>>
        ElementMap;
    typedef typename MortarView::template Codim<0>::Entity MortarElement;

    enum {
        dimWorld = MortarView::dimensionworld
    };

public:
    /*!
   * \copydoc BaseMapper(const std::string &, const MortarView &, const
   * GridView<0> &, const GridView<1> &)
   */
    FaceFaceMapper(const std::string& file_name,
        const MortarView& mortarView,
        const GridView<0>& gridView0,
        const GridView<1>& gridView1)
        : ParentType(file_name, mortarView, gridView0, gridView1)
    {
        this->finalizeInit_(file_name);
    }

    /*!
   * \copydoc BaseMapper(const MortarView &, const GridView<0> &, const
   * GridView<1> &)
   */
    FaceFaceMapper(const MortarView& mortarView,
        const GridView<0>& gridView0,
        const GridView<1>& gridView1)
        : ParentType(mortarView, gridView0, gridView1)
    {
        this->finalizeInit_();
    }

    /*!
   * \copydoc BaseMapper::toIntersection
   */
    template <std::size_t i>
    const Intersection<i> toIntersection(MortarElement e) const
    {
        const auto mortarIdx = this->mortarView_.indexSet().index(e);
        return std::get<i>(map_)[mortarIdx];
    }

    /*!
   * \copydoc BaseMapper::toElement
   */
    template <std::size_t i>
    const typename GridView<i>::template Codim<0>::Entity toElement(
        MortarElement e) const
    {
        return toIntersection<i>(e).inside();
    }

    //! Returns the projection
    ElementMap& projectionMap() { return map_; }

    /*!
   * \brief Calculate the projections from intersections to mortar elements.
   *
   * \copydoc BaseMapper::calculateMap
   *
   * The projections are calculated by comparing the intersection centers with
   * the cell centers. This assumes that the subdomain grids and mortar grid are
   * matching.
   */
    void calculateMap()
    {
        this->template calculateFaceMap_<0>(this->gridView0_);
        this->template calculateFaceMap_<1>(this->gridView1_);
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation& asImp_() { return *static_cast<Implementation*>(this); }
    //! \copydoc asImp_()
    const Implementation& asImp_() const
    {
        return *static_cast<const Implementation*>(this);
    }

protected:
    ElementMap map_;
};

/*!
 * \ingroup MultiDomainModel
 *
 * \brief Class object to identify face to element projections between domains.
 *
 * This class finds and stores the projections from the faces of the first
 * subdomain to the cells of the mortar grid, and the cells of the second
 * subdomain to the cells of the mortar grid
 * ------- | |
 * |grid0| | | < grid1
 * |_____| | |
 *         ^
 *     mortarGrid
 *
 */
template <class TypeTag>
class FaceElementMapper : public BaseMapper<TypeTag> {
    typedef BaseMapper<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, CouplingMapper) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) MortarView;
    template <std::size_t i>
    using GridView = typename GET_PROP_TYPE(TypeTag,
        SubTypeTag)::template GridView<i>;
    template <std::size_t i>
    using Intersection = typename GridView<i>::Intersection;
    template <std::size_t i>
    using Element = typename GridView<i>::template Codim<0>::Entity;
    typedef std::tuple<std::vector<Intersection<0>>, std::vector<Element<1>>>
        ElementMap;
    typedef typename MortarView::template Codim<0>::Entity MortarElement;
    enum {
        dimWorld = MortarView::dimensionworld
    };

public:
    /*!
   * \copydoc BaseMapper(const std::string &, const MortarView &, const
   * GridView<0> &, const GridView<1> &)
   */
    FaceElementMapper(const std::string& file_name,
        const MortarView& mortarView,
        const GridView<0>& gridView0,
        const GridView<1>& gridView1)
        : ParentType(file_name, mortarView, gridView0, gridView1)
    {
        this->finalizeInit_(file_name);
    }

    /*!
   * \copydoc BaseMapper(const MortarView &, const GridView<0> &, const
   * GridView<1> &)
   */
    FaceElementMapper(const MortarView& mortarView,
        const GridView<0>& gridView0,
        const GridView<1>& gridView1)
        : ParentType(mortarView, gridView0, gridView1)
    {
        this->finalizeInit_();
    }

    /*!
   * \copydoc BaseMapper::toIntersection
   */
    template <std::size_t i, typename std::enable_if_t<(i == 0), int> = 0>
    const Intersection<i> toIntersection(MortarElement e) const
    {
        const auto mortarIdx = this->mortarView_.indexSet().index(e);
        return std::get<i>(map_)[mortarIdx];
    }

    /*!
   * \copydoc BaseMapper::toElement
   */
    template <std::size_t i, typename std::enable_if_t<(i == 0), int> = 0>
    const typename GridView<i>::template Codim<0>::Entity toElement(
        MortarElement e) const
    {
        return toIntersection<0>(e).inside();
    }

    /*!
   * \copydoc BaseMapper::toElement
   */
    template <std::size_t i, typename std::enable_if_t<(i == 1), int> = 0>
    const typename GridView<i>::template Codim<0>::Entity toElement(
        MortarElement e) const
    {
        const auto mortarIdx = this->mortarView_.indexSet().index(e);
        const Element<1> E = std::get<1>(map_)[mortarIdx];
        return E;
    }

    //! Returns the projection between domains
    ElementMap& projectionMap() { return map_; }

    /*!
   * \copydoc BaseMapper::setMapFromFile
   */
    void setMapFromFile(const std::string& file_name)
    {
        std::array<std::vector<int>, 4> indexMap;

        readIndicesFromFile(file_name, indexMap);
        setElementMapFromIndices(indexMap);
        setFaceMapFromIndices(indexMap);
    }

    /*!
   * \copydoc BaseMapper::calculateMap
   *
   * The projections are calculated by comparing the intersection centers of the
   * first subdomain with the mortar cell centers, and the cell centers of the
   * second subdomain with the mortar cell centers. This assumes that the
   * subdomain grids and mortar grid are matching.
   */
    void calculateMap()
    {
        this->template calculateFaceMap_<0>(this->gridView0_);
        this->template calculateElementMap_<1>(this->gridView1_);
    }

protected:
    // Reads the intersection and element indices from file
    void readIndicesFromFile(const std::string& file_name,
        std::array<std::vector<int>, 4>& indexMap) const
    {
        std::string line;
        std::ifstream myfile(file_name);
        if (myfile.is_open()) {
            int i = 0;
            while (getline(myfile, line)) {
                int idx;
                std::istringstream stream(line);
                while (stream >> idx)
                    indexMap[i].push_back(idx);
                i++;
            }
            myfile.close();
        } else
            throw std::runtime_error("Could not read file");
    }

    // Assigns intersections to mapper from index map
    void setFaceMapFromIndices(const std::array<std::vector<int>, 4>& indexMap)
    {
        const auto& idxSetMortar = this->mortarView_.indexSet();
        const auto& idxSetGrid = this->gridView0_.indexSet();
        auto& map = asImp_().projectionMap();
        for (const auto& mortarEle : elements(this->mortarView_)) {
            const int mortarIdx = idxSetMortar.index(mortarEle);
            bool foundMap = false;
            if (indexMap[0][mortarIdx + 1] - indexMap[0][mortarIdx] != 1)
                throw std::runtime_error("Mortar must map to exactly 1 face");

            const auto intersection = this->template intersection<0>(indexMap[1][mortarIdx]);
            std::get<0>(map_)[mortarIdx] = intersection;
        }
    }

    // Assigns intersections to mapper from index map
    void setElementMapFromIndices(const std::array<std::vector<int>, 4>& indexMap)
    {
        const auto& idxSetMortar = this->mortarView_.indexSet();
        const auto& idxSetGrid = this->gridView1_.indexSet();

        for (const auto& mortarEle : elements(this->mortarView_)) {
            const auto mortarIdx = idxSetMortar.index(mortarEle);
            if (indexMap[2][mortarIdx + 1] - indexMap[2][mortarIdx] != 1)
                throw std::runtime_error("Mortar must map to exactly 1 face");

            const auto ele = this->template element<1>(indexMap[3][mortarIdx]);
            int idx = idxSetGrid.index(ele);
            if (indexMap[3][mortarIdx] == idxSetGrid.index(ele))
                std::get<1>(map_)[mortarIdx] = ele;
            else
                throw std::runtime_error(
                    "Could not find map. Are you sure the projection is correct?");
        }
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation& asImp_() { return *static_cast<Implementation*>(this); }
    //! \copydoc asImp_()
    const Implementation& asImp_() const
    {
        return *static_cast<const Implementation*>(this);
    }

protected:
    ElementMap map_;
};

} // namespace Opm
