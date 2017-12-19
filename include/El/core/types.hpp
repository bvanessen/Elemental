/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TYPES_HPP
#define EL_TYPES_HPP

#include "El/core/Element/Complex/decl.hpp"

namespace El
{

// For the safe computation of products. The result is given by
//   product = rho * exp(kappa*n)
// where rho lies in (usually on) the unit circle and kappa is real-valued.
template<typename F>
struct SafeProduct
{
    F rho;
    Base<F> kappa;
    Int n;

    SafeProduct( Int numEntries );
};

// The basic eigenvalue structure of a Hermitian matrix
struct InertiaType
{
    Int numPositive, numNegative, numZero;
};

// At some point this enum will be a member of each matrix class and will
// be used to simplify the syntax of various operations
namespace MatrixClassNS {
enum MatrixClass
{
    UNSPECIFIED,
    GENERAL,
    HERMITIAN,
    HERMITIAN_LOWER,
    HERMITIAN_UPPER,
    SYMMETRIC,
    SYMMETRIC_LOWER,
    SYMMETRIC_UPPER,
    SKEW_SYMMETRIC,
    SKEW_SYMMETRIC_LOWER,
    SKEW_SYMMETRIC_UPPER,
    UNITARY,
    TRIANGULAR_LOWER,
    TRIANGULAR_LOWER_UNIT,
    TRIANGULAR_UPPER,
    TRIANGULAR_UPPER_UNIT,
    HESSENBERG_LOWER,
    HESSENBERG_UPPER,
    PERMUTATION,
    PERMUTATION_VECTOR,
    PIVOT_SEQUENCE,
    // Packed factorizations/decompositions
    LU_PACKED,
    QR_PACKED,
    RQ_PACKED,
    LQ_PACKED,
    QL_PACKED,
    BIDIAG_PACKED,
    TRIDIAG_LOWER_PACKED,
    TRIDIAG_UPPER_PACKED,
    HESSENBERG_LOWER_PACKED,
    HESSENBERG_UPPER_PACKED
};
}
using namespace MatrixClassNS;

namespace ConjugationNS {
enum Conjugation
{
    UNCONJUGATED,
    CONJUGATED
};
}
using namespace ConjugationNS;

namespace DistNS {
enum Dist
{
    MC,   // Col of a matrix distribution
    MD,   // Diagonal of a matrix distribution
    MR,   // Row of a matrix distribution
    VC,   // Col-major vector distribution
    VR,   // Row-major vector distribution
    STAR, // Give to every process
    CIRC  // Give to a single process
};
std::string DistToString( Dist distribution );
Dist StringToDist( std::string s );
}
using namespace DistNS;
typedef Dist Distribution;

namespace DistWrapNS {
enum DistWrap {
    ELEMENT,
    BLOCK
};
}
using namespace DistWrapNS;

// Return either the row or column piece of the implied diagonal distribution
// ==========================================================================

// Compile-time
// ------------
template<Dist U,Dist V> constexpr Dist DiagCol() { return ( U==STAR ? V : U ); }
template<Dist U,Dist V> constexpr Dist DiagRow() { return ( U==STAR ? U : V ); }
template<> constexpr Dist DiagCol<MC,Dist::MR>() { return MD; }
template<> constexpr Dist DiagRow<MC,Dist::MR>() { return STAR; }
template<> constexpr Dist DiagCol<MR,Dist::MC>() { return MD; }
template<> constexpr Dist DiagRow<MR,Dist::MC>() { return STAR; }

// Runtime
// -------
inline Dist DiagCol( Dist U, Dist V ) EL_NO_EXCEPT
{
    if( U == Dist::MC && V == Dist::MR )
        return MD;
    else if( U == Dist::MR && V == Dist::MC )
        return MD;
    else if( U == Dist::STAR )
        return V;
    else
        return U;
}
inline Dist DiagRow( Dist U, Dist V ) EL_NO_EXCEPT
{
    if( U == Dist::MC && V == Dist::MR )
        return STAR;
    else if( U == Dist::MR && V == Dist::MC )
        return STAR;
    else if( U == Dist::STAR )
        return U;
    else
        return V;
}

// Return a piece of a distribution which induces the given diagonal dist
// ======================================================================

// Compile-time
// ------------
template<Dist U,Dist V>
constexpr Dist DiagInvCol() { return ( U==STAR ? V : U ); }
template<Dist U,Dist V>
constexpr Dist DiagInvRow() { return ( U==STAR ? U : V ); }
template<> constexpr Dist DiagInvCol<MD,Dist::STAR>() { return MC; }
template<> constexpr Dist DiagInvRow<MD,Dist::STAR>() { return MR; }
template<> constexpr Dist DiagInvCol<STAR,Dist::MD>() { return MC; }
template<> constexpr Dist DiagInvRow<STAR,Dist::MD>() { return MR; }

// TODO: Runtime version?

// Union the distribution over its corresponding communicator
// ==========================================================
// Compile-time
// ------------
template<Dist U> constexpr Dist Collect()       { return STAR; }
template<>       constexpr Dist Collect<CIRC>() { return CIRC; }
// Run-time
// --------
inline Dist Collect( Dist U ) EL_NO_EXCEPT { return ( U==CIRC ? Dist::CIRC : Dist::STAR ); }

// Union the distribution over its corresponding partial communicator
// ==================================================================
// Compile-time
// ------------
template<Dist U> constexpr Dist Partial() { return U; }
template<>       constexpr Dist Partial<VC>() { return MC; }
template<>       constexpr Dist Partial<VR>() { return MR; }
// Run-time
// --------
inline Dist Partial( Dist U ) EL_NO_EXCEPT
{
    if( U == Dist::VC )
        return MC;
    else if( U == Dist::VR )
        return MR;
    else
        return U;
}

// Return the partial distribution that would be used for a partial union
// ======================================================================
// Compile-time
// ------------
template<Dist U,Dist V> constexpr Dist PartialUnionRow()          { return V;  }
template<>              constexpr Dist PartialUnionRow<VC,Dist::STAR>() { return MR; }
template<>              constexpr Dist PartialUnionRow<VR,Dist::STAR>() { return MC; }template<Dist U,Dist V> constexpr Dist PartialUnionCol()
{ return PartialUnionRow<V,U>(); }
// Run-time
// --------
inline Dist PartialUnionRow( Dist U, Dist V ) EL_NO_EXCEPT
{
    if( U == Dist::VC )
        return MR;
    else if( U == Dist::VR )
        return MC;
    else
        return V;
}
inline Dist PartialUnionCol( Dist U, Dist V ) EL_NO_EXCEPT
{ return PartialUnionRow( V, U ); }

// Return the product of two distributions
// =======================================
// Compile-time
// ------------
template<Dist U,Dist V> constexpr Dist ProductDist() { return CIRC; }
template<>              constexpr Dist ProductDist<MC,  Dist::MR  >() { return VC;   }
template<>              constexpr Dist ProductDist<MC,  Dist::STAR>() { return MC;   }
template<>              constexpr Dist ProductDist<MD,  Dist::STAR>() { return MD;   }
template<>              constexpr Dist ProductDist<MR,  Dist::MC  >() { return VR;   }
template<>              constexpr Dist ProductDist<MR,  Dist::STAR>() { return MR;   }
template<>              constexpr Dist ProductDist<STAR,Dist::MC  >() { return MC;   }
template<>              constexpr Dist ProductDist<STAR,Dist::MD  >() { return MD;   }
template<>              constexpr Dist ProductDist<STAR,Dist::MR  >() { return MR;   }
template<>              constexpr Dist ProductDist<STAR,Dist::STAR>() { return STAR; }
template<>              constexpr Dist ProductDist<STAR,Dist::VC  >() { return VC;   }
template<>              constexpr Dist ProductDist<STAR,Dist::VR  >() { return VR;   }
template<>              constexpr Dist ProductDist<VC,  Dist::STAR>() { return VC;   }
template<>              constexpr Dist ProductDist<VR,  Dist::STAR>() { return VR;   }
template<Dist U,Dist V>
constexpr Dist ProductDistPartner() { return STAR; }
template<>
constexpr Dist ProductDistPartner<CIRC,Dist::CIRC>() { return CIRC; }
// Runtime
// -------
inline Dist ProductDist( Dist U, Dist V ) EL_NO_EXCEPT
{
    if(      U == Dist::STAR ) return V;
    else if( V == Dist::STAR ) return U;
    else if( U == Dist::MC   && V == Dist::MR   ) return VC;
    else if( U == Dist::MR   && V == Dist::MC   ) return VR;
    else if( U == CIRC && V == CIRC ) return CIRC;
    else { return STAR; } // NOTE: This branch should not be possible
}
inline Dist ProductDistPartner( Dist U, Dist V ) EL_NO_EXCEPT
{
    if( U == CIRC && V == CIRC ) return CIRC;
    else                         return STAR;
}


namespace ForwardOrBackwardNS {
enum ForwardOrBackward
{
    FORWARD,
    BACKWARD
};
}
using namespace ForwardOrBackwardNS;

namespace GridOrderNS {
enum GridOrder
{
    ROW_MAJOR,
    COLUMN_MAJOR
};
}
using namespace GridOrderNS;

namespace LeftOrRightNS {
enum LeftOrRight
{
    LeftOrRight::LEFT,
    LeftOrRight::RIGHT
};
char LeftOrRightToChar( LeftOrRight side );
LeftOrRight CharToLeftOrRight( char c );
}
using namespace LeftOrRightNS;

namespace SortTypeNS {
enum SortType
{
    UNSORTED,
    DESCENDING,
    ASCENDING
};
}
using namespace SortTypeNS;

namespace NormTypeNS {
enum NormType
{
    ONE_NORM,           // Operator one norm
    INFINITY_NORM,      // Operator infinity norm
    ENTRYWISE_ONE_NORM, // One-norm of vectorized matrix
    MAX_NORM,           // Maximum entry-wise magnitude
    NUCLEAR_NORM,       // One-norm of the singular values
    FROBENIUS_NORM,     // Two-norm of the singular values
    TWO_NORM            // Infinity-norm of the singular values
};
}
using namespace NormTypeNS;

namespace OrientationNS {
enum Orientation
{
    Orientation::NORMAL,
    Orientation::TRANSPOSE,
    :ADJOINT
};
char OrientationToChar( Orientation orientation );
Orientation CharToOrientation( char c );
}
using namespace OrientationNS;

namespace UnitOrNonUnitNS {
enum UnitOrNonUnit
{
    NON_UNIT,
    UNIT
};
char UnitOrNonUnitToChar( UnitOrNonUnit diag );
UnitOrNonUnit CharToUnitOrNonUnit( char c );
}
using namespace UnitOrNonUnitNS;

namespace UpperOrLowerNS {
enum UpperOrLower
{
    LOWER,
    UPPER
};
char UpperOrLowerToChar( UpperOrLower uplo );
UpperOrLower CharToUpperOrLower( char c );
}
using namespace UpperOrLowerNS;

namespace VerticalOrHorizontalNS {
enum VerticalOrHorizontal
{
    VERTICAL,
    HORIZONTAL
};
}
using namespace VerticalOrHorizontalNS;

// TODO: Distributed file formats?
namespace FileFormatNS {
enum FileFormat
{
    AUTO, // Automatically detect from file extension
    ASCII,
    ASCII_MATLAB,
    BINARY,
    BINARY_FLAT,
    BMP,
    JPG,
    JPEG,
    MATRIX_MARKET,
    PNG,
    PPM,
    XBM,
    XPM,
    FileFormat_MAX // For detecting number of entries in enum
};
}
using namespace FileFormatNS;

} // namespace El

#endif // ifndef EL_TYPES_HPP
