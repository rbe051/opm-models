

namespace Opm
{

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
    FaceFaceMapper(const std::string &file_name, const MortarView &mortarView, const GridView<0> &gridView0, const GridView<1> &gridView1)
        : mortarView_{mortarView}, gridView0_{gridView0}, gridView1_{gridView1}
    {
        throw std::runtime_error("Not implemented: Can not read face face mapping");
    }

    FaceFaceMapper(const MortarView &mortarView, const GridView<0> &gridView1, const GridView<1> &gridView2)
        : mortarView_{mortarView}, gridView0_{gridView1}, gridView1_{gridView2}
    {
        std::get<0>(map_).resize(mortarView_.size(0));
        std::get<1>(map_).resize(mortarView_.size(0));
        findMap<0>(gridView0_);
        findMap<1>(gridView1_, 1.0);
    }

    template <std::size_t i>
    const Intersection<i> toIntersection(MortarElement e) const
    {
        const auto mortarIdx = mortarView_.indexSet().index(e);
        return std::get<i>(map_)[mortarIdx];
    }

    template <std::size_t i>
    const typename GridView<i>::template Codim<0>::Entity toElement(MortarElement e) const
    {
        return toIntersection<i>(e).inside();
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
                        //map[idxSetMortar.index(e)] = idxSetGrid.index(facet);
                        std::get<i>(map_)[idxSetMortar.index(mortarEle)] = intersection;
                        foundMap = true;
                        break;
                    }
                }
            }
            if (!foundMap)
            {
                assert(false), "Could not find map. Are you sure the grids match?";
            }
        }
    }

public:
    ElementMap map_;
    const MortarView &mortarView_;
    const GridView<0> &gridView0_;
    const GridView<1> &gridView1_;
};

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
    FaceElementMapper(const std::string &file_name, const MortarView &mortarView, const GridView<0> &gridView0, const GridView<1> &gridView1)
        : mortarView_{mortarView}, gridView0_{gridView0}, gridView1_{gridView1}
    {
        std::get<0>(map_).resize(mortarView_.size(0));
        std::get<1>(map_).resize(mortarView_.size(0));
        setMapFromFile(file_name);
    }

    FaceElementMapper(const MortarView &mortarView, const GridView<0> &gridView0, const GridView<1> &gridView1)
        : mortarView_{mortarView}, gridView0_{gridView0}, gridView1_{gridView1}
    {
        std::get<0>(map_).resize(mortarView_.size(0));
        std::get<1>(map_).resize(mortarView_.size(0));
        findMap();
    }

    void setMapFromFile(const std::string &file_name)
    {
        std::array<std::vector<int>, 4> indexMap;

        readIndicesFromFile(file_name, indexMap);
        setFaceMapFromIndices(indexMap);
        setElementMapFromIndices(indexMap);
    }

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

public:
    ElementMap map_;
    const MortarView &mortarView_;
    const GridView<0> &gridView0_;
    const GridView<1> &gridView1_;
}; // namespace Opm

} // namespace Opm
