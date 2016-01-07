/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SVD_GOLUBREINSCH_HPP
#define EL_SVD_GOLUBREINSCH_HPP

#include "./Util.hpp"

namespace El {
namespace svd {

template<typename F>
inline void
GolubReinsch
( DistMatrix<F>& A,
  ElementalMatrix<Base<F>>& s, 
  DistMatrix<F>& V )
{
    DEBUG_ONLY(CSE cse("svd::GolubReinsch"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = Min( m, n );
    const Int offdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );

    // Bidiagonalize A
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,STAR> tP(g), tQ(g);
    Bidiag( A, tP, tQ );

    // Grab copies of the diagonal and sub/super-diagonal of A
    auto d_MD_STAR = GetRealPartOfDiagonal(A);
    auto e_MD_STAR = GetRealPartOfDiagonal(A,offdiagonal);

    // NOTE: lapack::BidiagQRAlg expects e to be of length k
    typedef Base<F> Real;
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                               eHat_STAR_STAR( k, 1, g );
    auto e_STAR_STAR = eHat_STAR_STAR( IR(0,k-1), ALL );
    e_STAR_STAR = e_MD_STAR;

    // Initialize U and VAdj to the appropriate identity matrices
    DistMatrix<F,VC,STAR> U_VC_STAR( g );
    U_VC_STAR.AlignWith( A );
    Identity( U_VC_STAR, m, k );
    DistMatrix<F,STAR,VC> VAdj_STAR_VC( g );
    VAdj_STAR_VC.AlignWith( V );
    Identity( VAdj_STAR_VC, k, n );

    // Compute the SVD of the bidiagonal matrix and accumulate the Givens
    // rotations into our local portion of U and VAdj
    Matrix<F>& ULoc = U_VC_STAR.Matrix();
    Matrix<F>& VAdjLoc = VAdj_STAR_VC.Matrix();
    lapack::BidiagQRAlg
    ( uplo, k, VAdjLoc.Width(), ULoc.Height(),
      d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(), 
      VAdjLoc.Buffer(), VAdjLoc.LDim(), 
      ULoc.Buffer(), ULoc.LDim() );

    // Make a copy of A (for the Householder vectors) and pull the necessary 
    // portions of U and VAdj into a standard matrix dist.
    auto B( A );
    if( m >= n )
    {
        auto AT = A( IR(0,n  ), ALL );
        auto AB = A( IR(n,END), ALL );
        auto UT_VC_STAR = U_VC_STAR( IR(0,n), ALL );
        AT = UT_VC_STAR;
        Zero( AB );
        Adjoint( VAdj_STAR_VC, V );
    }
    else
    {
        auto VAdjL_STAR_VC = VAdj_STAR_VC( IR(0,k), IR(0,m) );
        auto VT = V( IR(0,m  ), ALL );
        auto VB = V( IR(m,END), ALL );
        Adjoint( VAdjL_STAR_VC, VT );
        Zero( VB );
    }

    // Backtransform U and V
    bidiag::ApplyQ( LEFT, NORMAL, B, tQ, A );
    bidiag::ApplyP( LEFT, NORMAL, B, tP, V );

    // Copy out the appropriate subset of the singular values
    Copy( d_STAR_STAR, s );
}

template<typename F>
inline void
GolubReinsch
( ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& s, 
  ElementalMatrix<F>& VPre )
{
    DEBUG_ONLY(CSE cse("svd::GolubReinsch"))
    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> VProx( VPre );
    auto& A = AProx.Get();
    auto& V = VProx.Get();
    GolubReinsch( A, s, V );
}

#ifdef EL_HAVE_FLA_BSVD
template<typename F>
inline void
GolubReinschFlame
( DistMatrix<F>& A,
  ElementalMatrix<Base<F>>& s, 
  DistMatrix<F>& V )
{
    DEBUG_ONLY(CSE cse("svd::GolubReinschFlame"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = Min( m, n );
    const Int offdiagonal = ( m>=n ? 1 : -1 );

    // Bidiagonalize A
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,STAR> tP(g), tQ(g);
    Bidiag( A, tP, tQ );

    // Grab copies of the diagonal and sub/super-diagonal of A
    auto d_MD_STAR = GetRealPartOfDiagonal(A);
    auto e_MD_STAR = GetRealPartOfDiagonal(A,offdiagonal);

    // In order to use serial QR kernels, we need the full bidiagonal matrix
    // on each process
    typedef Base<F> Real;
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                               e_STAR_STAR( e_MD_STAR );

    // Initialize U and VAdj to the appropriate identity matrices
    DistMatrix<F,VC,STAR> U_VC_STAR(g), V_VC_STAR(g);
    U_VC_STAR.AlignWith( A );
    V_VC_STAR.AlignWith( V );
    Identity( U_VC_STAR, m, k );
    Identity( V_VC_STAR, n, k );

    // Since libFLAME, to the best of my current knowledge, only supports the
    // upper-bidiagonal case, we may instead work with the adjoint in the 
    // lower-bidiagonal case.
    if( m >= n )
    {
        flame::BidiagSVD
        ( k, U_VC_STAR.LocalHeight(), V_VC_STAR.LocalHeight(),
          d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          U_VC_STAR.Buffer(), U_VC_STAR.LDim(),
          V_VC_STAR.Buffer(), V_VC_STAR.LDim() );
    }
    else
    {
        flame::BidiagSVD
        ( k, V_VC_STAR.LocalHeight(), U_VC_STAR.LocalHeight(),
          d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          V_VC_STAR.Buffer(), V_VC_STAR.LDim(),
          U_VC_STAR.Buffer(), U_VC_STAR.LDim() );
    }

    // Make a copy of A (for the Householder vectors) and pull the necessary 
    // portions of U and V into a standard matrix dist.
    auto B( A );
    if( m >= n )
    {
        auto UT_VC_STAR = U_VC_STAR( IR(0,n), IR(0,k) );
        auto AT = A( IR(0,n), ALL );
        auto AB = A( IR(n,END), ALL );
        AT = UT_VC_STAR;
        Zero( AB );
        V = V_VC_STAR;
    }
    else
    {
        auto VT_VC_STAR = V_VC_STAR( IR(0,m), IR(0,k) );
        auto VT = V( IR(0,m  ), ALL );
        auto VB = V( IR(m,END), ALL ); 
        VT = VT_VC_STAR;
        Zero( VB );
    }

    // Backtransform U and V
    bidiag::ApplyQ( LEFT, NORMAL, B, tQ, A );
    bidiag::ApplyP( LEFT, NORMAL, B, tP, V );

    // Copy out the appropriate subset of the singular values
    Copy( d_STAR_STAR, s );
}

template<typename F>
inline void
GolubReinschFlame
( ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& s, 
  ElementalMatrix<F>& VPre )
{
    DEBUG_ONLY(CSE cse("svd::GolubReinschFlame"))
    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> VProx( VPre );
    auto& A = AProx.Get();
    auto& V = VProx.Get();
    GolubReinschFlame( A, s, V );
}

template<>
inline void
GolubReinsch
( ElementalMatrix<double>& A,
  ElementalMatrix<double>& s, 
  ElementalMatrix<double>& V )
{
    DEBUG_ONLY(CSE cse("svd::GolubReinsch"))
    GolubReinschFlame( A, s, V );
}

template<>
inline void
GolubReinsch
( ElementalMatrix<Complex<double>>& A,
  ElementalMatrix<double>& s, 
  ElementalMatrix<Complex<double>>& V )
{
    DEBUG_ONLY(CSE cse("svd::GolubReinsch"))
    GolubReinschFlame( A, s, V );
}
#endif // EL_HAVE_FLA_BSVD

template<typename F>
inline void
GolubReinsch( DistMatrix<F>& A, ElementalMatrix<Base<F>>& s )
{
    DEBUG_ONLY(CSE cse("svd::GolubReinsch"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = Min( m, n );
    const Int offdiagonal = ( m>=n ? 1 : -1 );

    // Bidiagonalize A
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,STAR> tP(g), tQ(g);
    Bidiag( A, tP, tQ );

    // Grab copies of the diagonal and sub/super-diagonal of A
    auto d_MD_STAR = GetRealPartOfDiagonal(A);
    auto e_MD_STAR = GetRealPartOfDiagonal(A,offdiagonal);

    // In order to use serial DQDS kernels, we need the full bidiagonal matrix
    // on each process
    //
    // NOTE: lapack::BidiagDQDS expects e to be of length k
    typedef Base<F> Real;
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                               eHat_STAR_STAR( k, 1, g );
    auto e_STAR_STAR = eHat_STAR_STAR( IR(0,k-1), ALL );
    e_STAR_STAR = e_MD_STAR;

    // Compute the singular values of the bidiagonal matrix via DQDS
    lapack::BidiagDQDS( k, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer() );

    // Copy out the appropriate subset of the singular values
    Copy( d_STAR_STAR, s );
}

template<typename F>
inline void
GolubReinsch( ElementalMatrix<F>& APre, ElementalMatrix<Base<F>>& s )
{
    DEBUG_ONLY(CSE cse("svd::GolubReinsch"))
    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();
    GolubReinsch( A, s );
}

} // namespace svd
} // namespace El

#endif // ifndef EL_SVD_GOLUBREINSCH_HPP
