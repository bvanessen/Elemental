#ifndef EL_TYPES_COMPLEX_UTILS_HPP_
#define EL_TYPES_COMPLEX_UTILS_HPP_

#include "El/core/limits.hpp"
#include "El/Meta.hpp"
#include "El/Types/Complex_decl.hpp"
#include "El/Types/Complex_impl.hpp"
#include "El/Traits.hpp"

namespace El
{

// Public function declarations
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real NaiveDiv( const Real& a, const Real& b );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Complex<Real> NaiveDiv( const Real& a, const Complex<Real>& b );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Complex<Real> NaiveDiv( const Complex<Real>& a, const Real& b );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Complex<Real> NaiveDiv( const Complex<Real>& a, const Complex<Real>& b );

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real SmithDiv( const Real& a, const Real& b );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Complex<Real> SmithDiv( const Real& a, const Complex<Real>& b );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Complex<Real> SmithDiv( const Complex<Real>& a, const Real& b );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Complex<Real> SmithDiv( const Complex<Real>& a, const Complex<Real>& b );

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real SafeDiv( const Real& a, const Real& b );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Complex<Real> SafeDiv( const Real& a, const Complex<Real>& b );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Complex<Real> SafeDiv( const Complex<Real>& a, const Real& b );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Complex<Real> SafeDiv( const Complex<Real>& a, const Complex<Real>& b );

// c := a / b using the textbook algorithm
template<typename Real,typename=EnableIf<IsReal<Real>>>
void NaiveDiv
( const Real& aReal, const Real& aImag,
  const Real& bReal, const Real& bImag,
  Real& cReal,       Real& cImag );

// c := a / b using Smith's algorithm
// See Fig. 3 from Baudin and Smith
template<typename Real,typename=EnableIf<IsReal<Real>>>
void SmithDiv
( const Real& aReal, const Real& aImag,
  const Real& bReal, const Real& bImag,
  Real& cReal,       Real& cImag );

template<typename Real,typename=EnableIf<IsReal<Real>>>
void SafeDiv
( const Real& aReal, const Real& aImag,
  const Real& bReal, const Real& bImag,
  Real& cReal,       Real& cImag );

// Public function implementations

template<typename Real,typename>
void NaiveDiv
( const Real& aReal, const Real& aImag,
  const Real& bReal, const Real& bImag,
        Real& cReal,       Real& cImag )
{
    const Real den = bReal*bReal + bImag*bImag;
    cReal = (aReal*bReal + aImag*bImag) / den;
    cImag = (aImag*bReal - aReal*bImag) / den;
}

template<typename Real,typename>
void SmithDiv
( const Real& aReal, const Real& aImag,
  const Real& bReal, const Real& bImag,
        Real& cReal,       Real& cImag )
{
    if( Abs(bImag) <= Abs(bReal) )
    {
        const Real r = bImag/bReal;
        const Real den = bReal + bImag*r;
        cReal = (aReal + aImag*r) / den;
        cImag = (aImag - aReal*r) / den;
    }
    else
    {
        const Real r = bReal/bImag;
        const Real den = bReal*r + bImag;
        cReal = (aReal*r + aImag) / den;
        cImag = (aImag*r - aReal) / den;
    }
}

// SafeDiv helpers
namespace safe_div
{

template<typename Real,typename=EnableIf<IsReal<Real>>>
void InternalRealPart
( const Real& aReal, const Real& aImag,
  const Real& bReal, const Real& bImag,
  const Real& r,
  const Real& t,
        Real& result )
{
    const Real zero = Real(0);
    if( r != zero )
    {
        Real br = aImag*r;
        if( br != zero )
        {
            result = (aReal + br)*t;
        }
        else
        {
            result = aReal*t + (aImag*t)*r;
        }
    }
    else
    {
        result = (aReal + bImag*(aImag/bReal))*t;
    }
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void Subinternal
( const Real& aReal, const Real& aImag,
  const Real& bReal, const Real& bImag,
        Real& cReal,       Real& cImag )
{
    Real r = bImag/bReal;
    Real t = 1/(bReal + bImag*r);
    safe_div::InternalRealPart(aReal,aImag,bReal,bImag,r,t,cReal);
    safe_div::InternalRealPart(aImag,-aReal,bReal,bImag,r,t,cImag);
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void Internal
( const Real& aReal, const Real& aImag,
  const Real& bReal, const Real& bImag,
        Real& cReal,       Real& cImag )
{
    if( Abs(bImag) <= Abs(bReal) )
    {
        safe_div::Subinternal(aReal,aImag,bReal,bImag,cReal,cImag);
    }
    else
    {
        safe_div::Subinternal(aImag,aReal,bImag,bReal,cReal,cImag);
        cImag = -cImag;
    }
}

} // namespace safe_div

template<typename Real,typename>
void SafeDiv
( const Real& aReal, const Real& aImag,
  const Real& bReal, const Real& bImag,
        Real& cReal,       Real& cImag )
{
    const Real aMax = Max( Abs(aReal), Abs(aImag) );
    const Real bMax = Max( Abs(bReal), Abs(bImag) );
    const Real beta = 2;
    const Real overflow = limits::Max<Real>();
    const Real underflow = limits::SafeMin<Real>();
    const Real eps = limits::Epsilon<Real>();
    const Real betaEpsSq = beta / (eps*eps);

    Real sigma=1;
    Real aRealScaled=aReal, aImagScaled=aImag,
         bRealScaled=bReal, bImagScaled=bImag;
    if( aMax >= overflow/2 )
    {
        aRealScaled /= 2;
        aImagScaled /= 2;
        sigma *= 2;
    }
    if( bMax >= overflow/2 )
    {
        bRealScaled /= 2;
        bImagScaled /= 2;
        sigma /= 2;
    }
    if( aMax <= underflow*beta/eps )
    {
        aRealScaled *= betaEpsSq;
        aImagScaled *= betaEpsSq;
        sigma /= betaEpsSq;
    }
    if( bMax <= underflow*beta/eps )
    {
        bRealScaled *= betaEpsSq;
        bImagScaled *= betaEpsSq;
        sigma *= betaEpsSq;
    }

    safe_div::Internal
    ( aRealScaled, aImagScaled,
      bRealScaled, bImagScaled,
      cReal,       cImag );
    cReal *= sigma;
    cImag *= sigma;
}

template<typename Real,typename>
Real NaiveDiv( const Real& a, const Real& b )
{ return a / b; }
template<typename Real,typename>
Complex<Real> NaiveDiv
( const Real& a,
  const Complex<Real>& b )
{
    Real cReal, cImag;
    NaiveDiv( a, Real(0), b.real(), b.imag(), cReal, cImag );
    return Complex<Real>(cReal,cImag);
}
template<typename Real,typename>
Complex<Real> NaiveDiv
( const Complex<Real>& a,
  const Real& b )
{ return a / b; }
template<typename Real,typename>
Complex<Real> NaiveDiv
( const Complex<Real>& a,
  const Complex<Real>& b )
{
    Real cReal, cImag;
    NaiveDiv( a.real(), a.imag(), b.real(), b.imag(), cReal, cImag );
    return Complex<Real>(cReal,cImag);
}

template<typename Real,typename>
Real SmithDiv( const Real& a, const Real& b )
{ return a / b; }
template<typename Real,typename>
Complex<Real> SmithDiv
( const Real& a,
  const Complex<Real>& b )
{
    Real cReal, cImag;
    SmithDiv( a, Real(0), b.real(), b.imag(), cReal, cImag );
    return Complex<Real>(cReal,cImag);
}
template<typename Real,typename>
Complex<Real> SmithDiv
( const Complex<Real>& a,
  const Real& b )
{ return a / b; }
template<typename Real,typename>
Complex<Real> SmithDiv
( const Complex<Real>& a,
  const Complex<Real>& b )
{
    Real cReal, cImag;
    SmithDiv( a.real(), a.imag(), b.real(), b.imag(), cReal, cImag );
    return Complex<Real>(cReal,cImag);
}

template<typename Real,typename>
Real SafeDiv( const Real& a, const Real& b )
{ return a / b; }
template<typename Real,typename>
Complex<Real> SafeDiv
( const Real& a,
  const Complex<Real>& b )
{
    Real cReal, cImag;
    SafeDiv( a, Real(0), b.real(), b.imag(), cReal, cImag );
    return Complex<Real>(cReal,cImag);
}
template<typename Real,typename>
Complex<Real> SafeDiv
( const Complex<Real>& a,
  const Real& b )
{ return a / b; }
template<typename Real,typename>
Complex<Real> SafeDiv
( const Complex<Real>& a,
  const Complex<Real>& b )
{
    Real cReal, cImag;
    SafeDiv( a.real(), a.imag(), b.real(), b.imag(), cReal, cImag );
    return Complex<Real>(cReal,cImag);
}

}// namespace El
#endif // EL_TYPES_COMPLEX_UTILS_HPP_
