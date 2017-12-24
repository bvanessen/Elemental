/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "El/core/DistMatrix/Abstract.hpp"
#include "El/core/Matrix.hpp"
#include "El/core/Permutation.hpp"
#include "El/Types/Enums.hpp"
#include "El/lapack_like/factor.hpp"

namespace El
{

template<typename Field>
InertiaType Inertia
( UpperOrLower uplo,
  Matrix<Field>& A,
  const LDLPivotCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    if( uplo == UpperOrLower::UPPER )
        LogicError("This option not yet supported");
    Permutation p;
    Matrix<Field> dSub;
    LDL( A, dSub, p, true, ctrl );
    return ldl::Inertia( GetRealPartOfDiagonal(A), dSub );
}

template<typename Field>
InertiaType Inertia
( UpperOrLower uplo,
  AbstractDistMatrix<Field>& APre,
  const LDLPivotCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    if( uplo == UpperOrLower::UPPER )
        LogicError("This option not yet supported");

    DistMatrixReadProxy<Field,Field,Dist::MC,Dist::MR> AProx( APre );
    auto& A = AProx.Get();

    DistPermutation p( A.Grid() );
    DistMatrix<Field,Dist::MD,Dist::STAR> dSub( A.Grid() );
    LDL( A, dSub, p, true, ctrl );
    return ldl::Inertia( GetRealPartOfDiagonal(A), dSub );
}

#define PROTO(Field) \
  template InertiaType Inertia \
  ( UpperOrLower uplo, \
    Matrix<Field>& A, \
    const LDLPivotCtrl<Base<Field>>& ctrl ); \
  template InertiaType Inertia \
  ( UpperOrLower uplo, \
    AbstractDistMatrix<Field>& A, \
    const LDLPivotCtrl<Base<Field>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
