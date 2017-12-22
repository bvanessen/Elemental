#ifndef EL_TYPES_ENTRY_HPP_
#define EL_TYPES_ENTRY_HPP_

#include "El/Typedefs.hpp"

namespace El
{


template<typename Real>
struct Entry
{
    Int i, j;
    Real value;

    static bool Lesser( const Entry<Real>& a, const Entry<Real>& b )
    {
        return a.value < b.value;
    }

    static bool Greater( const Entry<Real>& a, const Entry<Real>& b )
    {
        return a.value > b.value;
    }
};// struct Entry

template<typename Real>
struct Entry<Complex<Real>>
{
    Int i, j;
    Complex<Real> value;

    static bool Lesser(
        const Entry<Complex<Real>>& a, const Entry<Complex<Real>>& b)
    {
        return ((RealPart(a.value) < RealPart(b.value)) ||
                ((RealPart(a.value) == RealPart(b.value)) &&
                 (ImagPart(a.value) < ImagPart(b.value))));
    }
    static bool Greater(
        const Entry<Complex<Real>>& a, const Entry<Complex<Real>>& b)
    {
        return ((RealPart(a.value) > RealPart(b.value)) ||
                ((RealPart(a.value) == RealPart(b.value)) &&
                 (ImagPart(a.value) > ImagPart(b.value))));
    }
};// struct Entry<Complex<Real>>

}// namespace El

#endif // EL_TYPES_ENTRY_HPP_
