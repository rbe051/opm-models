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
 * \copydoc Opm::MultiDomainBaseModel
 */
#ifndef OPM_Multi_DOMAIN_MODEL
#define OPM_Multi_DOMAIN_MODEL

#include <opm/models/utils/basicproperties.hh>
#include <opm/models/multiDomain/multidomainproperties.hh>
#include <opm/models/multiDomain/multidomainsimulator.hh>

#include <chrono>

namespace Opm
{
template <class TypeTag>
class MultiDomainBaseModel;
}

BEGIN_PROPERTIES
NEW_TYPE_TAG(MultiDomainBaseModel, INHERITS_FROM(MultiDomain));
NEW_PROP_TAG(Model);
NEW_PROP_TAG(CouplerTypeTag);
NEW_PROP_TAG(Simulator);

SET_TYPE_PROP(MultiDomainBaseModel, Model,
              Opm::MultiDomainBaseModel<TypeTag>);
SET_TYPE_PROP(MultiDomainBaseModel, Simulator,
              Opm::MultiDomainSimulator<TypeTag>);
END_PROPERTIES

namespace Opm
{
/*!
 * \ingroup MultiDomainModel
 *
 * \brief The base class for the  multidomain model.
 * 
 * A multidomain model consist of several subdomains which are initiated as
 * standard OPM models. This class is a container to initiate, couple, and
 * access all submodels.
 * 
*/
template <class TypeTag>
class MultiDomainBaseModel
{
    typedef typename GET_PROP_TYPE(TypeTag, Linearizer) Linearizer;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;

    using SubTypes = typename GET_PROP_TYPE(TypeTag, SubTypeTag);
    using CouplerTypes = typename GET_PROP_TYPE(TypeTag, CouplerTypeTag);

    template <std::size_t i>
    using SubSimulator = typename SubTypes::template Simulator<i>;
    template <std::size_t i>
    using SubModel = typename SubTypes::template Model<i>;
    template <std::size_t i>
    using CouplerSubDomain = typename CouplerTypes::template SubDomain<i>;
    template <std::size_t i>
    using Coupler = typename CouplerTypes::template Coupler<i>;

    using ModelsPtr = typename SubTypes::template TupleOfPtr<SubModel>;
    using CopulerPtr = typename CouplerTypes::template TupleOfSharedPtr<Coupler>;
    using GlobalEqVector = typename SubTypes::GlobalEqVector;
    using GlobalSolutionVector = typename SubTypes::SolutionVector;

public:
    MultiDomainBaseModel(Simulator& simulator) : linearizer_(new Linearizer()), simulator_{&simulator}
    {
        auto start_domain = std::chrono::high_resolution_clock::now();

        using namespace Dune::Hybrid;
        auto start_coupling = std::chrono::high_resolution_clock::now();
        forEach(integralRange(Dune::Hybrid::size(couplers_)), [&](const auto couplerK) {
            typename CouplerSubDomain<couplerK>::IndexI I;
            typename CouplerSubDomain<couplerK>::IndexJ J;

            std::get<couplerK>(couplers_).reset(new Coupler<couplerK>(simulator_->template subSimulator<I>(),
                                                                      simulator_->template subSimulator<J>()));
        });

        forEach(integralRange(numDomains), [&](const auto domainI) {
            std::get<domainI>(models_) = &simulator_->template model<domainI>();
        });

        auto start_linearizer = std::chrono::high_resolution_clock::now();
        linearizer_->init(*simulator_, couplers_);
        auto end_linearizer = std::chrono::high_resolution_clock::now();

        std::cout << "Time model initiation:\n"
                  << "   Domain: " << std::chrono::duration_cast<std::chrono::milliseconds>(start_coupling - start_domain).count() << '\n'
                  << "   Coupling: " << std::chrono::duration_cast<std::chrono::milliseconds>(start_linearizer - start_coupling).count() << '\n'
                  << "   Linearizer: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_linearizer - start_linearizer).count() << '\n';
    }

    /*!
     * \brief Returns the operator linearizer for the global jacobian of
     *        the problem.
     */
    Linearizer &linearizer()
    {
        return *linearizer_;
    }

    /*!
     * \brief Applies the initial solution for all degrees of freedom to which the model
     *        applies.
     */
    void applyInitialSolution()
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(numDomains), [&](const auto domainI) {
            std::get<domainI>(models_)->applyInitialSolution();
        });
    }

    /*!
     * \brief Reference to the solution at a given history index as a multiple type block vector.
     *
     * \param timeIdx The index of the solution used by the time discretization.
     */

    GlobalSolutionVector &solution(unsigned timeIdx)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(numDomains), [&](const auto domainI) {
            solution_[domainI] = std::get<domainI>(models_)->solution(timeIdx);
        });
        return solution_;
    }

    /*!
     * \brief Returns true if the cache for intensive quantities is enabled
     */
    bool storeIntensiveQuantities() const
    {
        bool store = false;
        using namespace Dune::Hybrid;
        forEach(integralRange(numDomains), [&](const auto domainI) {
            store = store || std::get<domainI>(models_)->storeIntensiveQuantities();
        });
        return store;
    }

    /*!
     * \brief Invalidate the whole intensive quantity cache for time index.
     *
     * \param timeIdx The index used by the time discretization.
     */
    void invalidateIntensiveQuantitiesCache(unsigned timeIdx) const
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(numDomains), [&](const auto domainI) {
            std::get<domainI>(models_)->invalidateIntensiveQuantitiesCache(timeIdx);
        });
    }

    /*!
     * \brief Called by the problem if a time integration was
     *        successful, post processing of the solution is done and
     *        the result has been written to disk.
     *
     * This should prepare the model for the next time integration.
     */
    void advanceTimeLevel()
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(numDomains), [&](const auto domainI) {
            std::get<domainI>(models_)->advanceTimeLevel();
        });
    }


    /*!
     * \brief Move the intensive quantities for a given time index to the back.
     *
     * This method should only be called by the time discretization.
     *
     * \param numSlots The number of time step slots for which the
     *                 hints should be shifted.
     */

    void shiftIntensiveQuantityCache(unsigned numSlots = 1)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(numDomains), [&](const auto domainI) {
            std::get<domainI>(models_)->shiftIntensiveQuantityCache(numSlots);
        });
    }

    /*!
     * \brief Register all run-time parameters for the model.
     */
    static void registerParameters()
    {        
        Dune::index_constant<SubTypes::numSubDomains> tempNum;
        using namespace Dune::Hybrid;
        Dune::index_constant<CouplerTypes::numSubCouplers> couplerNum;
        forEach(integralRange(couplerNum), [&](const auto domainI) {
            Coupler<domainI>::registerParameters();
        });
    }

    /*!
     * \brief Get the simulator
     */
    Simulator &simulator()
    {
        return *simulator_;
    }

    /*!
     * \brief Get the model of subdomain i.
     */
    template <std::size_t i>
    SubModel<i> &model()
    {
        return *std::get<i>(models_);
    }

    /*!
     * \brief Assign the subdomain solutions from a global solution vector.
     */
    void setSolution(GlobalSolutionVector & vec)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(numDomains), [&](const auto domainI) {
            auto &localSolution = std::get<domainI>(models_)->solution(/*timeIdx=*/0);
            for (int i = 0; i < std::get<domainI>(models_)->numTotalDof(); ++i)
            {
                for (size_t phaseIdx = 0; phaseIdx < 2; ++phaseIdx)
                {
                    localSolution[i][phaseIdx] = vec[domainI][i][phaseIdx];
                }
            }
        });
    }

    /*!
     * \brief Write the relevant secondary variables of the current
     *        solution into an VTK output file.
     */
    void writeOutput()
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(numDomains), [&](const auto domainI) {
            simulator_->template subSimulator<domainI>().problem().writeOutput(domainI);
        });
    }

private:
    Dune::index_constant<SubTypes::numSubDomains> numDomains;

    Simulator* simulator_;
    ModelsPtr models_;
    CopulerPtr couplers_;
    std::unique_ptr<Linearizer> linearizer_;
    GlobalSolutionVector solution_;
};
} // namespace Opm

#endif
