#ifndef EL_TYPES_VALUEINT_HPP_
#define EL_TYPES_VALUEINT_HPP_

#include "El/Types/Complex_decl.hpp"
#include "El/Types/Complex_impl.hpp"

namespace El
{

template<typename Real>
struct ValueInt
{
    Real value;
    Int index;

    static bool Lesser(const ValueInt<Real>& a, const ValueInt<Real>& b)
    {
        return a.value < b.value;
    }

    static bool Greater(const ValueInt<Real>& a, const ValueInt<Real>& b)
    {
        return a.value > b.value;
    }
};// struct ValueInt

template<typename Real>
struct ValueInt<Complex<Real>>
{
    Complex<Real> value;
    Int index;

    static bool Lesser(
        const ValueInt<Complex<Real>>& a, const ValueInt<Complex<Real>>& b)
    {
        return RealPart(a.value) < RealPart(b.value) ||
               (RealPart(a.value) == RealPart(b.value) &&
                ImagPart(a.value) < ImagPart(b.value));
    }
    static bool Greater(
        const ValueInt<Complex<Real>>& a, const ValueInt<Complex<Real>>& b)
    {
        return RealPart(a.value) > RealPart(b.value) ||
               (RealPart(a.value) == RealPart(b.value) &&
                ImagPart(a.value) > ImagPart(b.value));
    }
};// struct ValueInt<Complex<Real>>

}// namespace El
#endif // EL_TYPES_VALUEINT_HPP_
