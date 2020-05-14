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
 * \copydoc Opm::MultiDomainVanguard
 */
#ifndef EWOMS_MULTI_DOMAIN_VANGUARD_HH
#define EWOMS_MULTI_DOMAIN_VANGUARD_HH

#include <opm/models/io/basevanguard.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>

#include "opm/grid/polyhedralgrid.hh"
#include "opm/grid/cart_grid.h"

#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <vector>
#include <memory>

namespace Opm
{

template <class TypeTag>
class MultiDomainVanguard;

} // namespace Opm

BEGIN_PROPERTIES

NEW_TYPE_TAG(MultiDomainVanguard);

// declare the properties required by the for the structured grid simulator vanguard
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(GridFile);
NEW_PROP_TAG(GridDim);
NEW_PROP_TAG(WorldDim);
NEW_PROP_TAG(GridGlobalRefinements);

// set the Grid and Vanguard properties
SET_TYPE_PROP(MultiDomainVanguard, Grid, Dune::PolyhedralGrid<GET_PROP_VALUE(TypeTag, GridDim), GET_PROP_VALUE(TypeTag, WorldDim)>);

SET_TYPE_PROP(MultiDomainVanguard, Vanguard, Opm::MultiDomainVanguard<TypeTag>);

END_PROPERTIES

namespace Opm
{

/*!
 * \ingroup TestProblems
 *
 * \brief Helper class for grid instantiation of the lens problem.
 */
template <class TypeTag>
class MultiDomainVanguard : public BaseVanguard<TypeTag>
{
    typedef BaseVanguard<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridFile) GridFile;

    typedef std::unique_ptr<Grid> GridPointer;

    static const int dim = Grid::dimension;

public:
    /*!
     * \brief Register all run-time parameters for the structured grid simulator vanguard.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, GridGlobalRefinements,
                             "The number of global refinements of the grid "
                             "executed after it was loaded");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, GridFile,
                             "The file name of the file to load");
        EWOMS_REGISTER_PARAM(TypeTag, int, GridDim,
                             "The dimension of the grid");
        EWOMS_REGISTER_PARAM(TypeTag, int, WorldDim,
                             "The dimension of the world");
    }

    /*!
     * \brief Create the grid for the lens problem
     */
    MultiDomainVanguard(Simulator &simulator)
        : ParentType(simulator), verbose_{false}
    {
        
        const std::string gridFileName = EWOMS_GET_PARAM(TypeTag, std::string, GridFile);
        unsigned numRefinments = EWOMS_GET_PARAM(TypeTag, unsigned, GridGlobalRefinements);

        const char *c_str = gridFileName.c_str();
        if (verbose_)
        {
            std::cout << "reading grid " << std::endl
                      << std::flush;
        }
        UnstructuredGrid *grid = read_grid(c_str);
        if (grid==nullptr)
            throw std::runtime_error("RuntimeError: MultiDomainVanguard could not read grid file: " + gridFileName+ ". Are you sure the filename is correct?");
            
        if (verbose_)
        {
            std::cout << "finished " << std::endl
                      << std::flush;
            std::cout << "printing grid" << std::endl;
            print_grid(grid);
        }

        Grid polyGrid(*grid);
        GridPointer polygrid(new Grid(*grid));
        gridPtr_ = std::move(polygrid);
        if (numRefinments > 0)
            gridPtr_->globalRefine(static_cast<int>(numRefinments));

        this->finalizeInit_();
    }

    /*!
     * \brief Return a reference to the grid object.
     */
    Grid &grid()
    {
        return *gridPtr_;
    }

    /*!
     * \brief Return a constant reference to the grid object.
     */
    const Grid &grid() const
    {
        return *gridPtr_;
    }

private:
    GridPointer gridPtr_;
    bool verbose_;
};

} // namespace Opm

#endif