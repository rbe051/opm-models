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
 * \brief Test for the immisicible VCVF discretization with only a single phase
 */
#include "config.h"
#include "problems/benchmark1problem.hh"

#include <dune/istl/solvers.hh>
#include <dune/istl/umfpack.hh>


int main(int argc, char **argv)
{
    typedef TTAG(DupletModel) MixedDimModelTypeTag;
    typedef typename GET_PROP_TYPE(MixedDimModelTypeTag, SubTypeTag) SubTypes;
    typedef typename GET_PROP_TYPE(MixedDimModelTypeTag, CouplerTypeTag) CouplerTypes;
    typedef typename GET_PROP_TYPE(MixedDimModelTypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(MixedDimModelTypeTag, SubTypeTag)::JacobianMatrix JacobianMatrix;
    typedef typename GET_PROP_TYPE(MixedDimModelTypeTag, SubTypeTag)::GlobalEqVector GlobalEqVector;

    // Register all parameters
    Model::registerParameters();

    Dune::index_constant<SubTypes::numSubDomains> numDomains;
    using namespace Dune::Hybrid;
    forEach(integralRange(numDomains), [&](const auto typeI) {
        EWOMS_END_PARAM_REGISTRATION(SubTypes::template TypeTag<typeI>);
    });
    Dune::index_constant<CouplerTypes::numSubCouplers> numCoupers;
    using namespace Dune::Hybrid;
    forEach(integralRange(numCoupers), [&](const auto typeI) {
        EWOMS_END_PARAM_REGISTRATION(CouplerTypes::template TypeTag<typeI>);
    });

    // Initiate model
    Model model;
    model.applyInitialSolution();
    auto &solution = model.solution(0);
    // Calculate jacobian and residual
    model.linearizer().linearize();
    auto &jac = model.linearizer().jacobian();
    auto &res = model.linearizer().residual();

    // Set up linear solver
    const auto &bcrs = Opm::MatrixConverter<JacobianMatrix>::multiTypeToBCRSMatrix(jac);
    const auto &blockVec = Opm::VectorConverter<GlobalEqVector>::multiTypeToBlockVector(res);

    using MatrixBlock = typename Dune::FieldMatrix<double, 1, 1>;
    using BCRSMatrix = typename Dune::BCRSMatrix<MatrixBlock>;
    using VectorBlock = typename Dune::FieldVector<double, 1>;
    using SolutionVector = typename Dune::BlockVector<VectorBlock>;

    SolutionVector result(blockVec.size());
    SolutionVector x(blockVec);

    Dune::MatrixAdapter<BCRSMatrix, SolutionVector, SolutionVector> linearOperator(bcrs);
    Dune::SeqSOR<BCRSMatrix, SolutionVector, SolutionVector> preconditioner(bcrs, 1.0, 0.1);
    Dune::BiCGSTABSolver<SolutionVector> cg(linearOperator, preconditioner, 1e-9, 1050, 2);
    Dune::InverseOperatorResult statistics;
    //cg.apply(result, x, statistics);
    Dune::UMFPack<BCRSMatrix> solver(bcrs);
    solver.apply(result, x, statistics);

    GlobalEqVector blockUpdate(res);
    Opm::VectorConverter<GlobalEqVector>::retrieveValues(blockUpdate, result);
    // intialize the vtk output module
    // Update solution
    using namespace Dune::Hybrid;
    forEach(integralRange(Dune::Hybrid::size(solution)), [&](const auto domainI) {
        // update the DOFs of the auxiliary equations
        size_t numDof = model.model<domainI>().numTotalDof();
        for (size_t dofIdx = 0; dofIdx < numDof; ++dofIdx)
        {
            for (size_t phaseIdx = 0; phaseIdx < 1; ++phaseIdx)
            {
                solution[domainI][dofIdx][phaseIdx] -= blockUpdate[domainI][dofIdx][phaseIdx];
            }
        }
    });
    model.setSolution(solution);
    model.writeOutput();
}