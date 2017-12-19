/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_APPLYPACKEDREFLECTORS_RUHF_HPP
#define EL_APPLYPACKEDREFLECTORS_RUHF_HPP

namespace El {
namespace apply_packed_reflectors {

//
// Since applying Householder transforms from vectors stored top-to-bottom
// implies that we will be forming a generalization of
//
//  (I - tau_0 v_0^T conj(v_0)) (I - tau_1 v_1^T conj(v_1)) =
//  I - [ v_0^T, v_1^T ] [ tau_0, -tau_0 tau_1 conj(v_0) v_1^T ] [ conj(v_0) ]
//                       [ 0,      tau_1                       ] [ conj(v_1) ],
//
// which has an upper-triangular center matrix, say S, we will form S as
// the inverse of a matrix T, which can easily be formed as
//
//   triu(T,1) = triu( conj(V V^H) ),
//   diag(T) = 1/householderScalars or 1/conj(householderScalars),
//
// where V is the matrix of Householder vectors and householderScalars is the
// vector of Householder reflection coefficients.
//
// V is stored row-wise in the matrix.
//

template<typename F>
void
RUHFUnblocked
( Conjugation conjugation,
  Int offset,
  const Matrix<F>& H,
  const Matrix<F>& householderScalars,
        Matrix<F>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Width() != H.Width() )
          LogicError("H and A must have the same width");
    )
    const Int nA = A.Width();
    const Int diagLength = H.DiagonalLength(offset);
    EL_DEBUG_ONLY(
      if( householderScalars.Height() != diagLength )
          LogicError
          ("householderScalars must be the same length as H's offset diag");
    )
    Matrix<F> hPanCopy, z;

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    for( Int k=0; k<diagLength; ++k )
    {
        const Int ki = k+iOff;
        const Int kj = k+jOff;

        auto hPan = H( IR(ki), IR(kj,nA) );
        auto ARight = A( ALL, IR(kj,nA) );
        const F tau = householderScalars(k);
        const F gamma = ( conjugation == CONJUGATED ? Conj(tau) : tau );

        // Convert to an explicit (scaled) Householder vector
        hPanCopy = hPan;
        hPanCopy(0) = F(1);

        // z := ARight hPan^T
        Gemv( Orientation::NORMAL, F(1), ARight, hPanCopy, z );
        // ARight := ARight (I - gamma hPan^T conj(hPan))
        Ger( -gamma, z, hPanCopy, ARight );
    }
}

template<typename F>
void
RUHFBlocked
( Conjugation conjugation,
  Int offset,
  const Matrix<F>& H,
  const Matrix<F>& householderScalars,
        Matrix<F>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Width() != H.Width() )
          LogicError("H and A must have the same width");
    )
    const Int nA = A.Width();
    const Int diagLength = H.DiagonalLength(offset);
    EL_DEBUG_ONLY(
      if( householderScalars.Height() != diagLength )
          LogicError
          ("householderScalars must be the same length as H's offset diag");
    )
    Matrix<F> HPanConj, SInv, Z;

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    const Int bsize = Blocksize();
    for( Int k=0; k<diagLength; k+=bsize )
    {
        const Int nb = Min(bsize,diagLength-k);
        const Int ki = k+iOff;
        const Int kj = k+jOff;

        auto HPan   = H( IR(ki,ki+nb), IR(kj,nA) );
        auto ARight = A( ALL,          IR(kj,nA) );
        auto householderScalars1 = householderScalars( IR(k,k+nb), ALL );

        // Convert to an explicit matrix of (scaled) Householder vectors
        Conjugate( HPan, HPanConj );
        MakeTrapezoidal( UpperOrLower::UPPER, HPanConj );
        FillDiagonal( HPanConj, F(1) );

        // Form the small triangular matrix needed for the UT transform
        Herk( UpperOrLower::UPPER, Orientation::NORMAL, Base<F>(1), HPanConj, SInv );
        FixDiagonal( conjugation, householderScalars1, SInv );

        // Z := ARight HPan^T
        Gemm( Orientation::NORMAL, Orientation::ADJOINT, F(1), ARight, HPanConj, Z );
        // Z := ARight HPan^T inv(SInv)
        Trsm( LeftOrRight::RIGHT, UpperOrLower::UPPER, Orientation::NORMAL, UnitOrNonUnit::NON_UNIT, F(1), SInv, Z );
        // ARight := ARight (I - HPan^T inv(SInv) conj(HPan))
        Gemm( Orientation::NORMAL, Orientation::NORMAL, F(-1), Z, HPanConj, F(1), ARight );
    }
}

template<typename F>
void
RUHF
( Conjugation conjugation,
  Int offset,
  const Matrix<F>& H,
  const Matrix<F>& householderScalars,
        Matrix<F>& A )
{
    EL_DEBUG_CSE
    const Int numLHS = A.Height();
    const Int blocksize = Blocksize();
    if( numLHS < blocksize )
    {
        RUHFUnblocked( conjugation, offset, H, householderScalars, A );
    }
    else
    {
        RUHFBlocked( conjugation, offset, H, householderScalars, A );
    }
}

template<typename F>
void
RUHFUnblocked
( Conjugation conjugation,
  Int offset,
  const AbstractDistMatrix<F>& H,
  const AbstractDistMatrix<F>& householderScalarsPre,
        AbstractDistMatrix<F>& APre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( APre.Width() != H.Width() )
          LogicError("A and H must have the same width");
      AssertSameGrids( H, householderScalarsPre, APre );
    )

    // We gather the entire set of Householder scalars at the start rather than
    // continually paying the latency cost of the broadcasts in a 'Get' call
    DistMatrixReadProxy<F,F,Dist::STAR,Dist::STAR>
      householderScalarsProx( householderScalarsPre );
    auto& householderScalars = householderScalarsProx.GetLocked();

    DistMatrixReadWriteProxy<F,F,Dist::MC,Dist::MR> AProx( APre );
    auto& A = AProx.Get();

    const Int nA = A.Width();
    const Int diagLength = H.DiagonalLength(offset);
    EL_DEBUG_ONLY(
      if( householderScalars.Height() != diagLength )
          LogicError
          ("householderScalars must be the same length as H's offset diag");
    )
    const Grid& g = H.Grid();
    auto hPan = unique_ptr<AbstractDistMatrix<F>>( H.Construct(g,H.Root()) );
    DistMatrix<F,Dist::STAR,Dist::MR> hPan_STAR_MR(g);
    DistMatrix<F,Dist::MC,Dist::STAR> z_MC_STAR(g);

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    for( Int k=0; k<diagLength; ++k )
    {
        const Int ki = k+iOff;
        const Int kj = k+jOff;

        auto ARight = A( ALL, IR(kj,nA) );
        const F tau = householderScalars.GetLocal( k, 0 );
        const F gamma = ( conjugation == CONJUGATED ? Conj(tau) : tau );

        // Convert to an explicit (scaled) Householder vector
        LockedView( *hPan, H, IR(ki), IR(kj,nA) );
        hPan_STAR_MR.AlignWith( ARight );
        Copy( *hPan, hPan_STAR_MR );
        hPan_STAR_MR.Set( 0, 0, F(1) );

        // z := ARight hPan^T
        z_MC_STAR.AlignWith( ARight );
        Zeros( z_MC_STAR, ARight.Height(), 1 );
        LocalGemv( Orientation::NORMAL, F(1), ARight, hPan_STAR_MR, F(0), z_MC_STAR );
        El::AllReduce( z_MC_STAR, ARight.RowComm() );

        // ARight := ARight (I - gamma hPan^T conj(hPan))
        LocalGer( -gamma, z_MC_STAR, hPan_STAR_MR, ARight );
    }
}

template<typename F>
void
RUHFBlocked
( Conjugation conjugation,
  Int offset,
  const AbstractDistMatrix<F>& H,
  const AbstractDistMatrix<F>& householderScalarsPre,
        AbstractDistMatrix<F>& APre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( APre.Width() != H.Width() )
          LogicError("A and H must have the same width");
      AssertSameGrids( H, householderScalarsPre, APre );
    )

    DistMatrixReadProxy<F,F,Dist::MC,Dist::STAR>
      householderScalarsProx( householderScalarsPre );
    auto& householderScalars = householderScalarsProx.GetLocked();

    DistMatrixReadWriteProxy<F,F,Dist::MC,Dist::MR> AProx( APre );
    auto& A = AProx.Get();

    const Int nA = A.Width();
    const Int diagLength = H.DiagonalLength(offset);
    EL_DEBUG_ONLY(
      if( householderScalars.Height() != diagLength )
          LogicError
          ("householderScalars must be the same length as H's offset diag");
    )
    const Grid& g = H.Grid();
    auto HPan = unique_ptr<AbstractDistMatrix<F>>( H.Construct(g,H.Root()) );
    DistMatrix<F> HPanConj(g);
    DistMatrix<F,Dist::STAR,Dist::VR  > HPan_STAR_VR(g);
    DistMatrix<F,Dist::STAR,Dist::MR  > HPan_STAR_MR(g);
    DistMatrix<F,Dist::STAR,Dist::STAR> householderScalars1_STAR_STAR(g);
    DistMatrix<F,Dist::STAR,Dist::STAR> SInv_STAR_STAR(g);
    DistMatrix<F,Dist::STAR,Dist::MC  > ZAdj_STAR_MC(g);
    DistMatrix<F,Dist::STAR,Dist::VC  > ZAdj_STAR_VC(g);

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    const Int bsize = Blocksize();
    for( Int k=0; k<diagLength; k+=bsize )
    {
        const Int nb = Min(bsize,diagLength-k);
        const Int ki = k+iOff;
        const Int kj = k+jOff;

        auto ARight = A( ALL, IR(kj,nA) );
        auto householderScalars1 = householderScalars( IR(k,k+nb), ALL );

        // Convert to an explicit matrix of (scaled) Householder vectors
        LockedView( *HPan, H, IR(ki,ki+nb), IR(kj,nA) );
        Conjugate( *HPan, HPanConj );
        MakeTrapezoidal( UpperOrLower::UPPER, HPanConj );
        FillDiagonal( HPanConj, F(1) );

        // Form the small triangular matrix needed for the UT transform
        HPan_STAR_VR = HPanConj;
        Zeros( SInv_STAR_STAR, nb, nb );
        Herk
        ( UpperOrLower::UPPER, Orientation::NORMAL,
          Base<F>(1), HPan_STAR_VR.LockedMatrix(),
          Base<F>(0), SInv_STAR_STAR.Matrix() );
        El::AllReduce( SInv_STAR_STAR, HPan_STAR_VR.RowComm() );
        householderScalars1_STAR_STAR = householderScalars1;
        FixDiagonal
        ( conjugation, householderScalars1_STAR_STAR, SInv_STAR_STAR );

        // Z := ARight HPan^T
        HPan_STAR_MR.AlignWith( ARight );
        HPan_STAR_MR = HPan_STAR_VR;
        ZAdj_STAR_MC.AlignWith( ARight );
        LocalGemm( Orientation::NORMAL, Orientation::ADJOINT, F(1), HPan_STAR_MR, ARight, ZAdj_STAR_MC );
        ZAdj_STAR_VC.AlignWith( ARight );
        Contract( ZAdj_STAR_MC, ZAdj_STAR_VC );

        // Z := ARight HPan^T inv(SInv)
        LocalTrsm
        ( LeftOrRight::LEFT, UpperOrLower::UPPER, Orientation::ADJOINT, UnitOrNonUnit::NON_UNIT, F(1), SInv_STAR_STAR, ZAdj_STAR_VC );

        // ARight := ARight (I - HPan^T inv(SInv) conj(HPan))
        ZAdj_STAR_MC = ZAdj_STAR_VC;
        LocalGemm
        ( Orientation::ADJOINT, Orientation::NORMAL, F(-1), ZAdj_STAR_MC, HPan_STAR_MR, F(1), ARight );
    }
}

template<typename F>
void
RUHF
( Conjugation conjugation,
  Int offset,
  const AbstractDistMatrix<F>& H,
  const AbstractDistMatrix<F>& householderScalars,
        AbstractDistMatrix<F>& A )
{
    EL_DEBUG_CSE
    const Int numLHS = A.Height();
    const Int blocksize = Blocksize();
    if( numLHS < blocksize )
    {
        RUHFUnblocked( conjugation, offset, H, householderScalars, A );
    }
    else
    {
        RUHFBlocked( conjugation, offset, H, householderScalars, A );
    }
}

} // namespace apply_packed_reflectors
} // namespace El

#endif // ifndef EL_APPLYPACKEDREFLECTORS_RUHF_HPP
