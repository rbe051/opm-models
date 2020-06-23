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
 * \brief Declare the properties used by the infrastructure code of
 *        the multidomain models.
 */
#ifndef EWOMS_MULTI_DOMAIN_PROPERTIES_HH
#define EWOMS_MULTI_DOMAIN_PROPERTIES_HH

#include "opm/models/utils/propertysystem.hh"
#include "opm/models/utils/parametersystem.hh"
#include <opm/models/multiDomain/multidomainlinearizer.hh>
#include <opm/models/io/multidomainvanguard.hh>

#include <dune/istl/bcrsmatrix.hh>

namespace Opm
{
template <typename... SubDomainTypeTags>
struct MultiDomainProperties;

} // namespace Opm

BEGIN_PROPERTIES
NEW_TYPE_TAG(MultiDomain, INHERITS_FROM(ParameterSystem, MultiDomainVanguard));
NEW_PROP_TAG(SubTypeTag);
NEW_PROP_TAG(CouplerTypeTag);
NEW_PROP_TAG(DomainI);
NEW_PROP_TAG(DomainJ);
NEW_PROP_TAG(Coupler);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(JacobianMatrix);
NEW_PROP_TAG(SolutionVector);
NEW_PROP_TAG(Linearizer);
NEW_PROP_TAG(Vanguard);

SET_INT_PROP(MultiDomain, GridGlobalRefinements, 0);
SET_TYPE_PROP(MultiDomain, Linearizer, Opm::MultiDomainLinearizer< TypeTag >);
SET_TYPE_PROP(MultiDomain, GridView, typename GET_PROP_TYPE(TypeTag, Grid)::LeafGridView);

END_PROPERTIES


namespace Opm
{
template <template <typename... Args> class Variadic, template <std::size_t> class Indexed, class U>
struct makeFromIndexedType;

template <template <typename... Args> class Variadic, template <std::size_t> class Indexed, std::size_t... IndexSeq>
struct makeFromIndexedType<Variadic, Indexed, std::index_sequence<IndexSeq...>>
{
    using type = Variadic<Indexed<IndexSeq>...>;
};

//! Helper type to determine whether a given type is a Dune::BCRSMatrix
template <class T>
struct isBCRSMatrix : public std::false_type
{
};

//! Helper type to determine whether a given type is a Dune::BCRSMatrix
template <class T>
struct isBCRSMatrix<Dune::BCRSMatrix<T>> : public std::true_type
{
};

//! Helper class to create a multitype matrix given the diagonal matrix blocks
template <class Scalar, class... JacobianBlocks>
class createMultiTypeBlockMatrixType
{
    //! TODO: replace by std::conjuction in C++17
    template <bool...>
    struct boolPack;
    template <bool... bools>
    using all_true = std::is_same<boolPack<bools..., true>, boolPack<true, bools...>>;

    static_assert(all_true<isBCRSMatrix<JacobianBlocks>::value...>::value, "Jacobian blocks have to be BCRSMatrices!");

    template <std::size_t id>
    using JacobianDiagBlock = typename std::tuple_element_t<id, std::tuple<JacobianBlocks...>>;

    template <std::size_t id>
    static constexpr decltype(auto) numEq()
    {
        return JacobianDiagBlock<id>::block_type::rows;
    }

    template <std::size_t id, class I>
    struct makeRow;

    template <std::size_t id, std::size_t... Is>
    struct makeRow<id, std::index_sequence<Is...>>
    {
        using type = Dune::MultiTypeBlockVector<Dune::BCRSMatrix<Opm::MatrixBlock<Scalar, numEq<id>(), numEq<Is>()>>...>;
    };

    template <class I>
    struct makeMatrix;

    template <std::size_t... Is>
    struct makeMatrix<std::index_sequence<Is...>>
    {
        using type = Dune::MultiTypeBlockMatrix<typename makeRow<Is, std::index_sequence<Is...>>::type...>;
    };

    using Indices = std::index_sequence_for<JacobianBlocks...>;

public:
    using type = typename makeMatrix<Indices>::type;
};

//! helper alias to create a tuple of ptr<...> from an indexed type
template <template <std::size_t> class T, class Indices>
struct MultiDomainTuplePtr
{
    template <std::size_t i>
    using PtrType = T<i>*;

    using type = typename makeFromIndexedType<std::tuple, PtrType, Indices>::type;
};

//! helper alias to create a tuple of shared_ptr<...> from an indexed type
template <template <std::size_t> class T, class Indices>
struct MultiDomainTupleSharedPtr
{
    template <std::size_t i>
    using PtrType = std::shared_ptr<T<i>>;

    using type = typename makeFromIndexedType<std::tuple, PtrType, Indices>::type;
};

//! helper alias to create a tuple of shared_ptr<const ...> from an indexed type
template <template <std::size_t> class T, class Indices>
struct MultiDomainTupleSharedPtrConst
{
    template <std::size_t i>
    using PtrType = std::shared_ptr<const T<i>>;

    using type = typename makeFromIndexedType<std::tuple, PtrType, Indices>::type;
};

//! helper alias to create the JacobianMatrix type
template <template <std::size_t> class SubDomainDiagBlocks, class Indices, class Scalar>
struct MultiDomainMatrixType
{
    template <typename... MatrixBlocks>
    using M = typename createMultiTypeBlockMatrixType<Scalar, MatrixBlocks...>::type::type;

    using type = typename makeFromIndexedType<M, SubDomainDiagBlocks, Indices>::type;
};

// Helper alias to access the subtypes of all subdomains
template <typename... SubDomainTypeTags>
struct MultiDomainProperties
{
    //! the number of subdomains
    static constexpr std::size_t numSubDomains = sizeof...(SubDomainTypeTags);

private:
    //! the type tag of a sub domain problem
    template <std::size_t id>
    using SubDomainTypeTag = typename std::tuple_element_t<id, std::tuple<SubDomainTypeTags...>>;

    //! helper alias to construct derived multidomain types like tuples
    using Indices = std::make_index_sequence<numSubDomains>;

    //! the scalar type of each sub domain
    template <std::size_t id>
    using SubDomainScalar = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Scalar);

    //! the solution type of each sub domain
    template <std::size_t id>
    using SubDomainSolutionVector = typename GET_PROP_TYPE(SubDomainTypeTag<id>, SolutionVector);
    //! the solution type of each sub domain
    template <std::size_t id>
    using SubDomainGlobalEqVector = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GlobalEqVector);

    //! the jacobian type of each sub domain
    template <std::size_t id>
    using SubDomainJacobianMatrix = typename GET_PROP_TYPE(SubDomainTypeTag<id>, SparseMatrixAdapter)::IstlMatrix;

public:
    /*
     * \brief sub domain types
     */
    //\{
    template <std::size_t id>
    using TypeTag = SubDomainTypeTag<id>;
    template <std::size_t id>
    using Simulator = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Simulator);
    template <std::size_t id>
    using Model = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Model);
    template <std::size_t id>
    using GridView = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridView);
    template <std::size_t id>
    using ElementContext = typename GET_PROP_TYPE(SubDomainTypeTag<id>, ElementContext);
    template <std::size_t id>
    using Stencil = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Stencil);
    template <std::size_t id>
    using IntensiveQuantities = typename GET_PROP_TYPE(SubDomainTypeTag<id>, IntensiveQuantities);

    template <std::size_t id>
    struct SubDomain
    {
        using Index = Dune::index_constant<id>;
        using TypeTag = SubDomainTypeTag<id>;
        using Problem = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Problem);
        using Model = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Model);

        using SolutionVector = typename GET_PROP_TYPE(SubDomainTypeTag<id>, SolutionVector);
    };

    //\}

    /*
     * \brief multi domain types
     */
    //\{

    //! the scalar type
    using Scalar = typename makeFromIndexedType<std::common_type_t, SubDomainScalar, Indices>::type;

    //!the solution vector type
    using SolutionVector = typename makeFromIndexedType<Dune::MultiTypeBlockVector, SubDomainSolutionVector, Indices>::type;

    //! the jacobian type
    using JacobianMatrix = typename MultiDomainMatrixType<SubDomainJacobianMatrix, Indices, Scalar>::type;

    //!the solution vector type
    using GlobalEqVector = typename makeFromIndexedType<Dune::MultiTypeBlockVector, SubDomainGlobalEqVector, Indices>::type;
    //!the solution vector type
    using GlobalSolutionVector = typename makeFromIndexedType<Dune::MultiTypeBlockVector, SubDomainSolutionVector, Indices>::type;

    //! helper alias to create tuple<...> from indexed type
    template <template <std::size_t> class T>
    using Tuple = typename makeFromIndexedType<std::tuple, T, Indices>::type;


    //! helper alias to create tuple<ptr<...>> from indexed type
    template <template <std::size_t> class T>
    using TupleOfPtr = typename MultiDomainTuplePtr<T, Indices>::type;
    //! helper alias to create tuple<std::shared_ptr<...>> from indexed type
    template <template <std::size_t> class T>
    using TupleOfSharedPtr = typename MultiDomainTupleSharedPtr<T, Indices>::type;

    //! helper alias to create tuple<std::shared_ptr<const ...>> from indexed type
    template <template <std::size_t> class T>
    using TupleOfSharedPtrConst = typename MultiDomainTupleSharedPtrConst<T, Indices>::type;
};

// Helper alias to access the subtypes of all subdomain couplings
template <typename... SubCouplerTypeTags>
struct MultiCouplerProperties
{
    //! the number of subdomains
    static constexpr std::size_t numSubCouplers = sizeof...(SubCouplerTypeTags);

private:
    //! helper alias to construct derived multidomain types like tuples
    using Indices = std::make_index_sequence<numSubCouplers>;
    //! the type tag of a sub domain problem
    template <std::size_t id>
    using SubCouplerTypeTag = typename std::tuple_element_t<id, std::tuple<SubCouplerTypeTags...>>;

public:
    /*
     * \brief sub domain types
     */
    //\{
    template <std::size_t id>
    struct SubDomain
    {
        using IndexI = typename GET_PROP_TYPE(SubCouplerTypeTag<id>, DomainI);
        using IndexJ = typename GET_PROP_TYPE(SubCouplerTypeTag<id>, DomainJ);
    };
    template <std::size_t id>
    using Coupler = typename GET_PROP_TYPE(SubCouplerTypeTag<id>, Coupler);

    template <std::size_t id>
    using TypeTag = SubCouplerTypeTag<id>;

    template <template <std::size_t> class T>
    using TupleOfSharedPtr = typename MultiDomainTupleSharedPtr<T, Indices>::type;
};
} //end namespace Opm

#endif
