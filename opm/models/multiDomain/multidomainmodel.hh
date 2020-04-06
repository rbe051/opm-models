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

namespace Opm
{
template <class TypeTag>
class MultiDomainBaseModel;
}

BEGIN_PROPERTIES
NEW_TYPE_TAG(MultiDomainBaseModel, INHERITS_FROM(MultiDomain));
NEW_PROP_TAG(Model);
NEW_PROP_TAG(CouplerTypeTag);

SET_TYPE_PROP(MultiDomainBaseModel, Model,
              Opm::MultiDomainBaseModel<TypeTag>);
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
    using SubTypes = typename GET_PROP_TYPE(TypeTag, SubTypeTag);
    using CouplerTypes = typename GET_PROP_TYPE(TypeTag, CouplerTypeTag);

    template <std::size_t i>
    using Simulator = typename SubTypes::template Simulator<i>;
    template <std::size_t i>
    using Model = typename SubTypes::template Model<i>;
    template <std::size_t i>
    using CouplerSubDomain = typename CouplerTypes::template SubDomain<i>;
    template <std::size_t i>
    using Coupler = typename CouplerTypes::template Coupler<i>;

    using Simulators = typename SubTypes::template TupleOfSharedPtr<Simulator>;
    using Couplers = typename CouplerTypes::template TupleOfSharedPtr<Coupler>;
    using GlobalEqVector = typename SubTypes::GlobalEqVector;
    using GlobalSolutionVector = typename SubTypes::SolutionVector;

public:
    MultiDomainBaseModel() : linearizer_(new Linearizer())
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(numDomains), [&](const auto domainI) {
            std::get<domainI>(simulators_).reset(new Simulator<domainI>());
        });

        forEach(integralRange(Dune::Hybrid::size(couplers_)), [&](const auto couplerK) {
            typename CouplerSubDomain<couplerK>::IndexI I;
            typename CouplerSubDomain<couplerK>::IndexJ J;

            std::get<couplerK>(couplers_).reset(new Coupler<couplerK>(*std::get<I>(simulators_),
                                                                      *std::get<J>(simulators_)));
        });

        linearizer_->init(simulators_, couplers_);
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
            std::get<domainI>(simulators_)->model().applyInitialSolution();
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
            solution_[domainI] = std::get<domainI>(simulators_)->model().solution(timeIdx);
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
            store = store || model<domainI>().storeIntensiveQuantities();
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
            std::get<domainI>(simulators_)->model().invalidateIntensiveQuantitiesCache(timeIdx);
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
            std::get<domainI>(simulators_)->model().advanceTimeLevel();
        });
    }

    /*!
     * \brief Set the current time step size to a given value.
     *
     * If the step size would exceed the length of the current
     * episode, the timeStep() method will take care that the step
     * size won't exceed the episode or the end of the simulation,
     * though.
     *
     * \param timeStepSize The new value for the time step size \f$\mathrm{[s]}\f$
     */
    template <class Scalar>
    void setTimeStepSize(Scalar value)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(numDomains), [&](const auto domainI) {
            std::get<domainI>(simulators_)->setTimeStepSize(value);
        });
    }

    /*!
     * \brief Set the current time step index to a given value.
     *
     * \param timeStepIndex The new value for the time step index
     */
    void setTimeStepIndex(unsigned value)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(numDomains), [&](const auto domainI) {
            std::get<domainI>(simulators_)->setTimeStepIndex(value);
        });
    }

    /*!
     * \brief Set the current simulated time, don't change the current
     *        time step index.
     *
     * \param t The time \f$\mathrm{[s]}\f$ which should be jumped to
     */
    template <class Scalar>
    void setTime(Scalar t)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(numDomains), [&](const auto domainI) {
            std::get<domainI>(simulators_)->setTime(t);
        });
    }
    /*!
     * \brief Set the current simulated time and the time step index.
     *
     * \param t The time \f$\mathrm{[s]}\f$ which should be jumped to
     * \param stepIdx The new time step index
     */
    template <class Scalar>
    void setTime(Scalar t, unsigned stepIdx)
    {
        setTime(t);
        setTimeStepIndex(stepIdx);
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
            std::get<domainI>(simulators_)->model().shiftIntensiveQuantityCache(numSlots);
        });
    }

    /*!
     * \brief Register all run-time parameters for the model.
     */
    static void registerParameters()
    {
        Dune::index_constant<SubTypes::numSubDomains> tempNum;
        using namespace Dune::Hybrid;
        forEach(integralRange(tempNum), [&](const auto domainI) {
            Simulator<domainI>::registerParameters();
        });
        Dune::index_constant<CouplerTypes::numSubCouplers> couplerNum;
        forEach(integralRange(couplerNum), [&](const auto domainI) {
            Coupler<domainI>::registerParameters();
        });
    }

    /*!
     * \brief Get the simulator of subdomain i.
     */
    template <std::size_t i>
    Simulator<i> &simulator()
    {
        return *std::get<i>(simulators_);
    }

    /*!
     * \brief Get the model of subdomain i.
     */
    template <std::size_t i>
    Model<i> &model()
    {
        return std::get<i>(simulators_)->model();
    }

    /*!
     * \brief Assign the subdomain solutions from a global solution vector.
     */
    void setSolution(GlobalSolutionVector vec)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(numDomains), [&](const auto domainI) {
            auto &localSolution = std::get<domainI>(simulators_)->model().solution(/*timeIdx=*/0);
            for (int i = 0; i < model<domainI>().numTotalDof(); ++i)
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
            std::get<domainI>(simulators_)->problem().writeOutput(domainI);
        });
    }

private:
    Dune::index_constant<SubTypes::numSubDomains> numDomains;

    Simulators simulators_;
    Couplers couplers_;
    std::unique_ptr<Linearizer> linearizer_;
    GlobalSolutionVector solution_;
};
} // namespace Opm

#endif
