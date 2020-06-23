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
 * \copydoc Opm::MultiDomainSimulator
 */
#ifndef OPM_MULTI_DOMAIN_SIMULATOR
#define OPM_MULTI_DOMAIN_SIMULATOR

#include <opm/models/multiDomain/multidomainproperties.hh>
#include <opm/models/utils/basicproperties.hh>

#include <chrono>

namespace Opm {
template <class TypeTag>
class MultiDomainSimulator;
}

BEGIN_PROPERTIES
NEW_TYPE_TAG(MultiDomainSimulator, INHERITS_FROM(MultiDomain));
NEW_PROP_TAG(Simulator);

SET_TYPE_PROP(MultiDomainSimulator, Model,
    Opm::MultiDomainSimulator<TypeTag>);
END_PROPERTIES

namespace Opm {
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
class MultiDomainSimulator {
    using SubTypes = typename GET_PROP_TYPE(TypeTag, SubTypeTag);
    using Model = typename GET_PROP_TYPE(TypeTag, Model);
    template <std::size_t i>
    using SubSimulator = typename SubTypes::template Simulator<i>;
    template <std::size_t i>
    using SubModel = typename SubTypes::template Model<i>;
    using SimulatorPtr = typename SubTypes::template TupleOfSharedPtr<SubSimulator>;

public:
    MultiDomainSimulator()
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(numDomains), [&](const auto domainI) {
            std::get<domainI>(simulators_).reset(new SubSimulator<domainI>());
        });

        model_.reset(new Model(*this));
    }
    /*!
     * \brief Register all run-time parameters for the Simulator.
     */
    static void registerParameters()
    {
        Dune::index_constant<SubTypes::numSubDomains> tempNum;
        using namespace Dune::Hybrid;
        forEach(integralRange(tempNum), [&](const auto domainI) {
            SubSimulator<domainI>::registerParameters();
        });

        Model::registerParameters();
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
     * \brief Get reference to the model
     */
    Model& model()
    {
        return *model_;
    }
    /*!
     * \brief Get the model of subdomain i.
     */
    template <std::size_t i>
    SubModel<i>& model()
    {
        return std::get<i>(simulators_)->model();
    }
    /*!
     * \brief Get the GridView of subdomain i.
     */
    template <std::size_t i>
    auto& gridView()
    {
        return std::get<i>(simulators_)->gridView();
    }
    /*!
     * \brief Get the simulator of subdomain i.
     */
    template <std::size_t i>
    SubSimulator<i>& subSimulator()
    {
        return *std::get<i>(simulators_);
    }

private:
    Dune::index_constant<SubTypes::numSubDomains> numDomains;
    SimulatorPtr simulators_;
    std::unique_ptr<Model> model_;
};
} //namespace Opm

#endif