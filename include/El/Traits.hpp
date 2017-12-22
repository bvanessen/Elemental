#ifndef EL_TRAITS_HPP_
#define EL_TRAITS_HPP_

#include <type_traits>

#include "El/Types/Complex_decl.hpp"

namespace El
{

template<typename S,typename T>
using IsSame = std::is_same<S,T>;

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
template<typename T> struct IsScalar<Complex<T>>
{ static const bool value=IsScalar<T>::value; };

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
template<typename T> struct IsField<Complex<T>>
{ static const bool value=IsField<T>::value; };

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
template<typename T> struct IsStdScalar<Complex<T>>
{ static const bool value=IsStdScalar<T>::value; };

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
template<typename T> struct IsStdField<Complex<T>>
{ static const bool value=IsStdField<T>::value; };

// For querying whether or not an element's type is complex
// --------------------------------------------------------
// NOTE: This does not guarantee that the type is a field
// NOTE: IsReal is not the negation of IsComplex
template<typename Real> struct IsReal
{ static const bool value=IsScalar<Real>::value; };
template<typename Real> struct IsReal<Complex<Real>>
{ static const bool value=false; };

template<typename Real> struct IsComplex
{ static const bool value=false; };
template<typename Real> struct IsComplex<Complex<Real>>
{ static const bool value=true; };

// Returning the underlying, or "base", real field
// -----------------------------------------------
// Note: The following is for internal usage only; please use Base
template<typename Real> struct BaseHelper                { typedef Real type; };
template<typename Real> struct BaseHelper<Complex<Real>> { typedef Real type; };

template<typename F> using Base = typename BaseHelper<F>::type;

} // namespace El

#endif // EL_TRAITS_HPP_
