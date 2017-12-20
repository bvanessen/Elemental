#ifndef EL_TRAITS_HPP_
#define EL_TRAITS_HPP_

#include "El/core/Element/Complex/decl.hpp"

namespace El
{

template<typename S,typename T>
using IsSame = std::is_same<S,T>;

template<typename Condition,class T=void>
using EnableIf = typename std::enable_if<Condition::value,T>::type;
template<typename Condition,class T=void>
using DisableIf = typename std::enable_if<!Condition::value,T>::type;

template<typename T>
struct IsIntegral { static const bool value = std::is_integral<T>::value; };
#ifdef HYDROGEN_HAVE_MPC
template<>
struct IsIntegral<BigInt> { static const bool value = true; };
#endif

// For querying whether an element's type is a scalar
// --------------------------------------------------
template<typename T> struct IsScalar
{ static const bool value=false; };
template<> struct IsScalar<unsigned>
{ static const bool value=true; };
template<> struct IsScalar<int>
{ static const bool value=true; };
template<> struct IsScalar<unsigned long>
{ static const bool value=true; };
template<> struct IsScalar<long int>
{ static const bool value=true; };
template<> struct IsScalar<unsigned long long>
{ static const bool value=true; };
template<> struct IsScalar<long long int>
{ static const bool value=true; };
template<> struct IsScalar<float>
{ static const bool value=true; };
template<> struct IsScalar<double>
{ static const bool value=true; };
template<> struct IsScalar<long double>
{ static const bool value=true; };
#ifdef HYDROGEN_HAVE_QD
template<> struct IsScalar<DoubleDouble>
{ static const bool value=true; };
template<> struct IsScalar<QuadDouble>
{ static const bool value=true; };
#endif
#ifdef HYDROGEN_HAVE_QUADMATH
template<> struct IsScalar<Quad>
{ static const bool value=true; };
#endif
#ifdef HYDROGEN_HAVE_MPC
template<> struct IsScalar<BigInt>
{ static const bool value=true; };
template<> struct IsScalar<BigFloat>
{ static const bool value=true; };
#endif

// For querying whether an element's type is a field
// -------------------------------------------------
template<typename T> struct IsField
{ static const bool value=false; };
template<> struct IsField<float>
{ static const bool value=true; };
template<> struct IsField<double>
{ static const bool value=true; };
template<> struct IsField<long double>
{ static const bool value=true; };
#ifdef HYDROGEN_HAVE_QD
template<> struct IsField<DoubleDouble>
{ static const bool value=true; };
template<> struct IsField<QuadDouble>
{ static const bool value=true; };
#endif
#ifdef HYDROGEN_HAVE_QUADMATH
template<> struct IsField<Quad>
{ static const bool value=true; };
#endif
#ifdef HYDROGEN_HAVE_MPC
template<> struct IsField<BigFloat>
{ static const bool value=true; };
#endif

// For querying whether an element's type is supported by the STL's math
// ---------------------------------------------------------------------
template<typename T> struct IsStdScalar
{ static const bool value=false; };
template<> struct IsStdScalar<unsigned>
{ static const bool value=true; };
template<> struct IsStdScalar<int>
{ static const bool value=true; };
template<> struct IsStdScalar<unsigned long>
{ static const bool value=true; };
template<> struct IsStdScalar<long int>
{ static const bool value=true; };
template<> struct IsStdScalar<unsigned long long>
{ static const bool value=true; };
template<> struct IsStdScalar<long long int>
{ static const bool value=true; };
template<> struct IsStdScalar<float>
{ static const bool value=true; };
template<> struct IsStdScalar<double>
{ static const bool value=true; };
template<> struct IsStdScalar<long double>
{ static const bool value=true; };
#ifdef HYDROGEN_HAVE_QUADMATH
template<> struct IsStdScalar<Quad>
{ static const bool value=true; };
#endif

// For querying whether an element's type is a field supported by STL
// ------------------------------------------------------------------
template<typename T> struct IsStdField
{ static const bool value=false; };
template<> struct IsStdField<float>
{ static const bool value=true; };
template<> struct IsStdField<double>
{ static const bool value=true; };
template<> struct IsStdField<long double>
{ static const bool value=true; };
#ifdef HYDROGEN_HAVE_QUADMATH
template<> struct IsStdField<Quad>
{ static const bool value=true; };
#endif

} // namespace El

#endif // EL_TRAITS_HPP_
