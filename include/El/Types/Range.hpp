#ifndef EL_TYPES_RANGE_HPP_
#define EL_TYPES_RANGE_HPP_

#include "El/Typedefs.hpp"

namespace El
{

template<typename T>
struct Range
{
    T beg, end;
    Range() : beg(0), end(0) { }
    Range( T begArg, T endArg ) : beg(begArg), end(endArg) { }

    Range<T> operator+( T shift ) const
    { return Range<T>(beg+shift,end+shift); }

    Range<T> operator-( T shift ) const
    { return Range<T>(beg-shift,end-shift); }

    Range<T> operator*( T scale ) const
    { return Range<T>(beg*scale,end*scale); }
};

static const Int END = -100;

template<>
struct Range<Int>
{
    Int beg, end;
    Range() : beg(0), end(0) { }
    Range( Int index ) : beg(index), end(index+1) { }
    Range( Int begArg, Int endArg ) : beg(begArg), end(endArg) { }

    Range<Int> operator+( Int shift ) const
    {
        if( end == END )
            throw std::logic_error("Unsupported shift");
        return Range<Int>(beg+shift,end+shift);
    }

    Range<Int> operator-( Int shift ) const
    {
        if( end == END )
            throw std::logic_error("Unsupported shift");
        return Range<Int>(beg-shift,end-shift);
    }

    Range<Int> operator*( Int scale ) const
    {
        if( end == END )
            throw std::logic_error("Unsupported scale");
        return Range<Int>(beg*scale,end*scale);
    }
};
typedef Range<Int> IR;

static const IR ALL(0,END);

template<typename T>
inline bool operator==( const Range<T>& a, const Range<T>& b )
{
    return a.beg == b.beg && a.end == b.end;
}

}// namespace El
#endif // EL_TYPES_RANGE_HPP_
