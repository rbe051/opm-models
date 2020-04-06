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
 * \brief 
 */

namespace Opm
{

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
class FaceFaceMapper
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) MortarView;
    template <std::size_t i>
    using GridView = typename GET_PROP_TYPE(TypeTag, SubTypeTag)::template GridView<i>;
    template <std::size_t i>
    using Intersection = typename GridView<i>::Intersection;
    typedef std::tuple<std::vector<Intersection<0>>, std::vector<Intersection<1>>> ElementMap;
    typedef typename MortarView::template Codim<0>::Entity MortarElement;

    enum
    {
        dimWorld = MortarView::dimensionworld
    };

public:
    /*!
     * \brief The constructor. Read mapping from file.
     */
    FaceFaceMapper(const std::string &file_name, const MortarView &mortarView, const GridView<0> &gridView0, const GridView<1> &gridView1)
        : mortarView_{mortarView}, gridView0_{gridView0}, gridView1_{gridView1}
    {
        throw std::runtime_error("Not implemented: Can not read face face mapping from file");
    }

    /*!
     * \brief The constructor. Read find mapping by comparing geometry.
     */
    FaceFaceMapper(const MortarView &mortarView, const GridView<0> &gridView1, const GridView<1> &gridView2)
        : mortarView_{mortarView}, gridView0_{gridView1}, gridView1_{gridView2}
    {
        std::get<0>(map_).resize(mortarView_.size(0));
        std::get<1>(map_).resize(mortarView_.size(0));
        findMap<0>(gridView0_);
        findMap<1>(gridView1_, 1.0);
    }

    /*!
     * \brief Returns the grid intersection that corresponds to the mortar element
     */
    template <std::size_t i>
    const Intersection<i> toIntersection(MortarElement e) const
    {
        const auto mortarIdx = mortarView_.indexSet().index(e);
        return std::get<i>(map_)[mortarIdx];
    }

    /*!
     * \brief Returns the grid element that correspond to the the mortar element
     */
    template <std::size_t i>
    const typename GridView<i>::template Codim<0>::Entity toElement(MortarElement e) const
    {
        return toIntersection<i>(e).inside();
    }

    /*!
     * \brief Returns the number of mortar elements
     */
    unsigned size() const
    {
        return mortarView_.size(0);
    }

    /*!
     * \brief Returns the number of subdomain elements
     */
    unsigned size(int i) const
    {
        switch (i)
        {
        case 0:
            return gridView0_.size(0);
        case 1:
            return gridView1_.size(0);
        default:
            throw std::runtime_error("FaceFaceMapper can only find the size of subdomain 0 or 1");
        }
    }

protected:
    /*!
     * \brief Calculate the projections from intersections to mortar elements.
     * 
     * The projections are calculated by comparing the intersection centers with
     * the cell centers. This assumes that the subdomain grids and mortar grid are
     * matching.
     */
    template <std::size_t i>
    void findMap(const GridView<i> &gridView, double offset = 0)
    {
        const auto &idxSetMortar = mortarView_.indexSet();
        const auto &idxSetGrid = gridView.indexSet();
        for (const auto &mortarEle : elements(mortarView_))
        {
            const auto &mc = mortarEle.geometry().center();

            bool foundMap = false;
            for (const auto &ele : elements(gridView))
            {
                for (const auto &intersection : intersections(gridView, ele))
                {
                    const auto &fc = intersection.geometry().center();
                    bool equal = true;
                    for (int dim = 0; dim < dimWorld; ++dim)
                    {
                        equal = equal && (std::abs(fc[dim] - mc[dim]) < 1e-10);
                    }
                    if (equal)
                    {
                        std::get<i>(map_)[idxSetMortar.index(mortarEle)] = intersection;
                        foundMap = true;
                        break;
                    }
                }
            }
            if (!foundMap)
                throw std::runtime_error("Could not find map. Are you sure the grids match?");
        }
    }

public:
    ElementMap map_;
    const MortarView &mortarView_;
    const GridView<0> &gridView0_;
    const GridView<1> &gridView1_;
};

/*!
 * \ingroup MultiDomainModel
 *
 * \brief Class object to identify face to element projections between domains.
 * 
 * This class finds and stores the projections from the faces of the first subdomain
 * to the cells of the mortar grid, and the cells of the second subdomain to the
 * cells of the mortar grid
 * ------- | |
 * |grid0| | | < grid1
 * |_____| | |
 *         ^       
 *     mortarGrid
 * 
*/
template <class TypeTag>
class FaceElementMapper
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) MortarView;
    template <std::size_t i>
    using GridView = typename GET_PROP_TYPE(TypeTag, SubTypeTag)::template GridView<i>;
    template <std::size_t i>
    using Intersection = typename GridView<i>::Intersection;
    template <std::size_t i>
    using Element = typename GridView<i>::template Codim<0>::Entity;
    typedef std::tuple<std::vector<Intersection<0>>, std::vector<Element<1>>> ElementMap;
    typedef typename MortarView::template Codim<0>::Entity MortarElement;
    enum
    {
        dimWorld = MortarView::dimensionworld
    };

public:
    /*!
     * \brief The constructor. Read mapping from file.
     */
    FaceElementMapper(const std::string &file_name, const MortarView &mortarView, const GridView<0> &gridView0, const GridView<1> &gridView1)
        : mortarView_{mortarView}, gridView0_{gridView0}, gridView1_{gridView1}
    {
        std::get<0>(map_).resize(mortarView_.size(0));
        std::get<1>(map_).resize(mortarView_.size(0));
        setMapFromFile(file_name);
    }

    /*!
     * \brief The constructor. Read find mapping by comparing geometry.
     */
    FaceElementMapper(const MortarView &mortarView, const GridView<0> &gridView0, const GridView<1> &gridView1)
        : mortarView_{mortarView}, gridView0_{gridView0}, gridView1_{gridView1}
    {
        std::get<0>(map_).resize(mortarView_.size(0));
        std::get<1>(map_).resize(mortarView_.size(0));
        findMap();
    }

    template <std::size_t i, typename std::enable_if_t<(i == 0), int> = 0>
    const Intersection<i> toIntersection(MortarElement e) const
    {
        const auto mortarIdx = mortarView_.indexSet().index(e);
        return std::get<i>(map_)[mortarIdx];
    }

    template <std::size_t i, typename std::enable_if_t<(i == 0), int> = 0>
    const typename GridView<i>::template Codim<0>::Entity toElement(MortarElement e) const
    {
        return toIntersection<0>(e).inside();
    }

    template <std::size_t i, typename std::enable_if_t<(i == 1), int> = 0>
    const typename GridView<i>::template Codim<0>::Entity toElement(MortarElement e) const
    {
        const auto mortarIdx = this->mortarView_.indexSet().index(e);
        const Element<1> E = std::get<1>(map_)[mortarIdx];
        return E;
    }

    unsigned size() const
    {
        return mortarView_.size(0);
    }
    unsigned size(int i) const
    {
        switch (i)
        {
        case 0:
            return gridView0_.size(0);
        case 1:
            return gridView1_.size(0);
        default:
            assert(false);
        }
    }

protected:

    /*!
     * \brief Reads the projections from file and assign it to the mapper
     */
    void setMapFromFile(const std::string &file_name)
    {
        std::array<std::vector<int>, 4> indexMap;

        readIndicesFromFile(file_name, indexMap);
        setFaceMapFromIndices(indexMap);
        setElementMapFromIndices(indexMap);
    }


    void findMap()
    {
        findFaceMap();
        findElementMap();
    }

    void findFaceMap(double offset = 0)
    {
        const auto &idxSetMortar = mortarView_.indexSet();
        const auto &idxSetGrid = gridView0_.indexSet();
        for (const auto &mortarEle : elements(mortarView_))
        {
            const auto &mc = mortarEle.geometry().center();
            std::cout << "mc: " << mc[0] << ",  " << mc[1] << std::endl;
            bool foundMap = false;
            for (const auto &ele : elements(gridView0_))
            {
                for (const auto &intersection : intersections(gridView0_, ele))
                {
                    const auto &fc = intersection.geometry().center();
                    bool equal = true;
                    for (int dim = 0; dim < dimWorld; ++dim)
                    {
                        equal = equal && (std::abs(fc[dim] - mc[dim]) < 1e-10);
                    }
                    if (equal)
                    {
                        //map[idxSetMortar.index(e)] = idxSetGrid.index(facet);
                        std::get<0>(map_)[idxSetMortar.index(mortarEle)] = intersection;
                        foundMap = true;
                        break;
                    }
                }
                if (foundMap)
                    break;
            }
            if (!foundMap)
            {
                assert(false), "Could not find map. Are you sure the grids match?";
            }
        }
    }

    void findElementMap(double offset = 0)
    {
        const auto &idxSetMortar = mortarView_.indexSet();
        const auto &idxSetGrid = gridView1_.indexSet();
        for (const auto &mortarEle : elements(mortarView_))
        {
            const auto &mc = mortarEle.geometry().center();
            bool foundMap = false;
            for (const auto &ele : elements(gridView1_))
            {
                const auto &fc = ele.geometry().center();
                bool equal = true;
                for (int dim = 0; dim < dimWorld; ++dim)
                {
                    equal = equal && (std::abs(fc[dim] - mc[dim]) < 1e-10);
                }
                if (equal)
                {
                    //map[idxSetMortar.index(e)] = idxSetGrid.index(facet);
                    std::get<1>(map_)[idxSetMortar.index(mortarEle)] = ele;
                    foundMap = true;
                    break;
                }
            }
            if (!foundMap)
            {
                assert(false), "Could not find map. Are you sure the grids match?";
            }
        }
    }

    /*!
     * \brief Reads the intersection and element indices from file
     */
    void readIndicesFromFile(const std::string &file_name, std::array<std::vector<int>, 4> &indexMap) const
    {
        std::string line;
        std::ifstream myfile(file_name);
        if (myfile.is_open())
        {
            int i = 0;
            while (getline(myfile, line))
            {
                int idx;
                std::istringstream stream(line);
                while (stream >> idx)
                    indexMap[i].push_back(idx);
                i++;
            }
            myfile.close();
        }
        else
            throw std::runtime_error("Could not read file");
    }

    /*!
     * \brief Assigns intersections to mapper from index map
     */
    void setFaceMapFromIndices(const std::array<std::vector<int>, 4> &indexMap)
    {
        const auto &idxSetMortar = mortarView_.indexSet();
        const auto &idxSetGrid = gridView0_.indexSet();
        for (const auto &mortarEle : elements(mortarView_))
        {
            const int mortarIdx = idxSetMortar.index(mortarEle);
            bool foundMap = false;
            if (indexMap[0][mortarIdx + 1] - indexMap[0][mortarIdx] != 1)
                throw std::runtime_error("Mortar must map to exactly 1 face");
            for (const auto &ele : elements(gridView0_))
            {
                for (const auto &intersection : intersections(gridView0_, ele))
                {
                    const int subIdx = intersection.indexInInside();
                    const int faceIdx = idxSetGrid.subIndex(ele, subIdx, 1);
                    if (indexMap[1][mortarIdx] == faceIdx)
                    {
                        //map[idxSetMortar.index(e)] = idxSetGrid.index(facet);
                        std::get<0>(map_)[mortarIdx] = intersection;
                        foundMap = true;
                        const auto &fc = intersection.geometry().center();
                        const auto &mc = mortarEle.geometry().center();
                        bool equal = true;
                        for (int dim = 0; dim < dimWorld; ++dim)
                        {
                            equal = equal && (std::abs(fc[dim] - mc[dim]) < 1e-10);
                        }
                        if (!equal)
                            assert(false);
                    }
                }
                if (foundMap)
                    break;
            }
            if (!foundMap)
            {
                throw std::runtime_error("Could not find map. Are you sure the projection is correct?");
            }
        }
    }
    void setElementMapFromIndices(const std::array<std::vector<int>, 4> &indexMap)
    {
        const auto &idxSetMortar = mortarView_.indexSet();
        const auto &idxSetGrid = gridView1_.indexSet();
        for (const auto &mortarEle : elements(mortarView_))
        {
            const auto mortarIdx = idxSetMortar.index(mortarEle);
            bool foundMap = false;
            if (indexMap[2][mortarIdx + 1] - indexMap[2][mortarIdx] != 1)
                throw std::runtime_error("Mortar must map to exactly 1 face");
            for (const auto &ele : elements(gridView1_))
            {
                int idx = idxSetGrid.index(ele);
                if (indexMap[3][mortarIdx] == idxSetGrid.index(ele))
                {
                    //map[idxSetMortar.index(e)] = idxSetGrid.index(facet);
                    std::get<1>(map_)[mortarIdx] = ele;
                    foundMap = true;
                    break;
                }
            }
            if (!foundMap)
            {
                throw std::runtime_error("Could not find map. Are you sure the projection is correct?");
            }
        }
    }

    ElementMap map_;
    const MortarView &mortarView_;
    const GridView<0> &gridView0_;
    const GridView<1> &gridView1_;
};

} // namespace Opm
