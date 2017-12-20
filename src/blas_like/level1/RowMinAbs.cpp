/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El/blas_like/level1.hpp"

namespace El {

template<typename Ring>
void RowMinAbs( const Matrix<Ring>& A, Matrix<Base<Ring>>& mins )
{
    EL_DEBUG_CSE
    typedef Base<Ring> RealRing;
    const Int m = A.Height();
    const Int n = A.Width();
    mins.Resize( m, 1 );
    Zero( mins );
    for( Int i=0; i<m; ++i )
    {
        RealRing rowMin = limits::Max<RealRing>();
        for( Int j=0; j<n; ++j )
            rowMin = Min(rowMin,Abs(A(i,j)));
        mins(i) = rowMin;
    }
}

template<typename Ring>
void RowMinAbsNonzero
( const Matrix<Ring>& A,
  const Matrix<Base<Ring>>& upperBounds,
        Matrix<Base<Ring>>& mins )
{
    EL_DEBUG_CSE
    typedef Base<Ring> RealRing;
    const Int m = A.Height();
    const Int n = A.Width();
    mins.Resize( m, 1 );
    Zero( mins );
    for( Int i=0; i<m; ++i )
    {
        RealRing rowMin = upperBounds(i);
        for( Int j=0; j<n; ++j )
        {
            const RealRing absVal = Abs(A(i,j));
            if( absVal > RealRing(0) )
                rowMin = Min(rowMin,absVal);
        }
        mins(i) = rowMin;
    }
}

template<typename Ring,Dist U,Dist V>
void RowMinAbs
( const DistMatrix<Ring,U,V>& A, DistMatrix<Base<Ring>,U,Dist::STAR>& mins )
{
    EL_DEBUG_CSE
    mins.AlignWith( A );
    mins.Resize( A.Height(), 1 );
    RowMinAbs( A.LockedMatrix(), mins.Matrix() );
    AllReduce( mins, A.RowComm(), mpi::MIN );
}

template<typename Ring,Dist U,Dist V>
void RowMinAbsNonzero
( const DistMatrix<Ring,U,V>& A,
  const DistMatrix<Base<Ring>,U,Dist::STAR>& upperBounds,
        DistMatrix<Base<Ring>,U,Dist::STAR>& mins )
{
    EL_DEBUG_CSE
    if( upperBounds.ColAlign() != A.ColAlign() )
        LogicError("upperBounds was not aligned with A");
    mins.AlignWith( A );
    mins.Resize( A.Height(), 1 );
    RowMinAbsNonzero
    ( A.LockedMatrix(), upperBounds.LockedMatrix(), mins.Matrix() );
    AllReduce( mins, A.RowComm(), mpi::MIN );
}


#define PROTO_DIST(Ring,U,V) \
  template void RowMinAbs \
  ( const DistMatrix<Ring,U,V>& X, DistMatrix<Base<Ring>,U,Dist::STAR>& mins ); \
  template void RowMinAbsNonzero \
  ( const DistMatrix<Ring,U,V>& X, \
    const DistMatrix<Base<Ring>,U,Dist::STAR>& upperBounds, \
          DistMatrix<Base<Ring>,U,Dist::STAR>& mins );

#define PROTO(Ring) \
  template void RowMinAbs \
  ( const Matrix<Ring>& X, Matrix<Base<Ring>>& mins ); \
  template void RowMinAbsNonzero \
  ( const Matrix<Ring>& X, \
    const Matrix<Base<Ring>>& upperBounds, \
          Matrix<Base<Ring>>& mins ); \
  PROTO_DIST(Ring,Dist::MC,  Dist::MR  ) \
  PROTO_DIST(Ring,Dist::MC,  Dist::STAR) \
  PROTO_DIST(Ring,Dist::MD,  Dist::STAR) \
  PROTO_DIST(Ring,Dist::MR,  Dist::MC  ) \
  PROTO_DIST(Ring,Dist::MR,  Dist::STAR) \
  PROTO_DIST(Ring,Dist::STAR,Dist::MC  ) \
  PROTO_DIST(Ring,Dist::STAR,Dist::MD  ) \
  PROTO_DIST(Ring,Dist::STAR,Dist::MR  ) \
  PROTO_DIST(Ring,Dist::STAR,Dist::STAR) \
  PROTO_DIST(Ring,Dist::STAR,Dist::VC  ) \
  PROTO_DIST(Ring,Dist::STAR,Dist::VR  ) \
  PROTO_DIST(Ring,Dist::VC,  Dist::STAR) \
  PROTO_DIST(Ring,Dist::VR,  Dist::STAR)

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
