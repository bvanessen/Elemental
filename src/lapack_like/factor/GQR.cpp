/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "El/lapack_like/factor.hpp"

#include "El/core/DistMatrix/Abstract.hpp"
#include "El/core/Matrix/decl.hpp"
#include "El/Types/Enums.hpp"

namespace El
{

template<typename F>
void GQR
( Matrix<F>& A,
  Matrix<F>& householderScalarsA,
  Matrix<Base<F>>& signatureA,
  Matrix<F>& B,
  Matrix<F>& householderScalarsB,
  Matrix<Base<F>>& signatureB )
{
    EL_DEBUG_CSE
    QR( A, householderScalarsA, signatureA );
    qr::ApplyQ( LeftOrRight::LEFT, Orientation::ADJOINT, A, householderScalarsA, signatureA, B );
    RQ( B, householderScalarsB, signatureB );
}

template<typename F>
void GQR
( AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<F>& householderScalarsA,
  AbstractDistMatrix<Base<F>>& signatureA,
  AbstractDistMatrix<F>& BPre,
  AbstractDistMatrix<F>& householderScalarsB,
  AbstractDistMatrix<Base<F>>& signatureB )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<F,F,Dist::MC,Dist::MR> AProx( APre ), BProx( BPre );
    auto& A = AProx.Get();
    auto& B = BProx.Get();

    QR( A, householderScalarsA, signatureA );
    qr::ApplyQ( LeftOrRight::LEFT, Orientation::ADJOINT, A, householderScalarsA, signatureA, B );
    RQ( B, householderScalarsB, signatureB );
}

namespace gqr {

template<typename F>
void ExplicitTriang( Matrix<F>& A, Matrix<F>& B )
{
    EL_DEBUG_CSE
    Matrix<F> householderScalarsA;
    Matrix<Base<F>> signatureA;
    QR( A, householderScalarsA, signatureA );
    qr::ApplyQ( LeftOrRight::LEFT, Orientation::ADJOINT, A, householderScalarsA, signatureA, B );
    MakeTrapezoidal( UpperOrLower::UPPER, A );
    rq::ExplicitTriang( B );
}

template<typename F>
void ExplicitTriang( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& BPre )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<F,F,Dist::MC,Dist::MR> AProx( APre ), BProx( BPre );
    auto& A = AProx.Get();
    auto& B = BProx.Get();

    const Grid& g = A.Grid();
    DistMatrix<F,Dist::MD,Dist::STAR> householderScalarsA(g);
    DistMatrix<Base<F>,Dist::MD,Dist::STAR> signatureA(g);
    QR( A, householderScalarsA, signatureA );
    qr::ApplyQ( LeftOrRight::LEFT, Orientation::ADJOINT, A, householderScalarsA, signatureA, B );
    MakeTrapezoidal( UpperOrLower::UPPER, A );
    rq::ExplicitTriang( B );
}

} // namespace gqr

#define PROTO(F) \
  template void GQR \
  ( Matrix<F>& A, \
    Matrix<F>& householderScalarsA, \
    Matrix<Base<F>>& signatureA, \
    Matrix<F>& B, \
    Matrix<F>& householderScalarsB, \
    Matrix<Base<F>>& signatureB ); \
  template void GQR \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& householderScalarsA, \
    AbstractDistMatrix<Base<F>>& signatureA, \
    AbstractDistMatrix<F>& B, \
    AbstractDistMatrix<F>& householderScalarsB, \
    AbstractDistMatrix<Base<F>>& signatureB ); \
  template void gqr::ExplicitTriang \
  ( Matrix<F>& A, Matrix<F>& B ); \
  template void gqr::ExplicitTriang \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
