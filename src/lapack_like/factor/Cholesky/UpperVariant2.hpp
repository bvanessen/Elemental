/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CHOLESKY_UPPER_VARIANT2_HPP
#define EL_CHOLESKY_UPPER_VARIANT2_HPP

#include "./UVar3.hpp"

// TODO: Reverse variants

namespace El {
namespace cholesky {

template<typename F>
void UpperVariant2Blocked( Matrix<F>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("Can only compute Cholesky factor of square matrices");
    )
    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A01 = A( ind0, ind1 );
        auto A02 = A( ind0, ind2 );
        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );

        Herk( UpperOrLower::UPPER, Orientation::ADJOINT, F(-1), A01, F(1), A11 );
        cholesky::UpperVariant3Unblocked( A11 );
        Gemm( Orientation::ADJOINT, Orientation::NORMAL, F(-1), A02, A01, F(1), A12 );
        Trsm( LeftOrRight::LEFT, UpperOrLower::UPPER, Orientation::ADJOINT, UnitOrNonUnit::NON_UNIT, F(1), A11, A12 );
    }
}

template<typename F>
void UpperVariant2Blocked( AbstractDistMatrix<F>& APre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( APre.Height() != APre.Width() )
          LogicError("Can only compute Cholesky factor of square matrices");
    )
    const Grid& grid = APre.Grid();

    DistMatrixReadWriteProxy<F,F,Dist::MC,Dist::MR> AProx( APre );
    auto& A = AProx.Get();

    DistMatrix<F,Dist::MC,  Dist::STAR> A01_MC_STAR(grid);
    DistMatrix<F,Dist::STAR,Dist::STAR> A11_STAR_STAR(grid);
    DistMatrix<F,Dist::STAR,Dist::VR  > A12_STAR_VR(grid);
    DistMatrix<F,Dist::MR,  Dist::STAR> X11Adj_MR_STAR(grid);
    DistMatrix<F,Dist::MR,  Dist::MC  > X11Adj_MR_MC(grid);
    DistMatrix<F,Dist::MR,  Dist::STAR> X12Adj_MR_STAR(grid);
    DistMatrix<F,Dist::MR,  Dist::MC  > X12Adj_MR_MC(grid);
    DistMatrix<F> X11(grid), X12(grid);

    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A01 = A( ind0, ind1 );
        auto A02 = A( ind0, ind2 );
        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );

        A01_MC_STAR.AlignWith( A01 );
        A01_MC_STAR = A01;
        X11Adj_MR_STAR.AlignWith( A01 );
        LocalGemm( Orientation::ADJOINT, Orientation::NORMAL, F(1), A01, A01_MC_STAR, X11Adj_MR_STAR );
        X11Adj_MR_MC.AlignWith( A11 );
        Contract( X11Adj_MR_STAR, X11Adj_MR_MC );
        X11.AlignWith( A11 );
        Adjoint( X11Adj_MR_MC, X11 );
        A11 -= X11;

        A11_STAR_STAR = A11;
        Cholesky( UpperOrLower::UPPER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        X12Adj_MR_STAR.AlignWith( A02 );
        LocalGemm( Orientation::ADJOINT, Orientation::NORMAL, F(1), A02, A01_MC_STAR, X12Adj_MR_STAR );
        X12Adj_MR_MC.AlignWith( A12 );
        Contract( X12Adj_MR_STAR, X12Adj_MR_MC );
        X12.AlignWith( A12 );
        Adjoint( X12Adj_MR_MC, X12 );
        A12 -= X12;

        A12_STAR_VR = A12;
        LocalTrsm
        ( LeftOrRight::LEFT, UpperOrLower::UPPER, Orientation::ADJOINT, UnitOrNonUnit::NON_UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );
        A12 = A12_STAR_VR;
    }
}

} // namespace cholesky
} // namespace El

#endif // ifndef EL_CHOLESKY_UPPER_VARIANT2_HPP
