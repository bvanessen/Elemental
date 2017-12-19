/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_INVERSE_HPD_CHOLESKYUVAR2_HPP
#define EL_INVERSE_HPD_CHOLESKYUVAR2_HPP

namespace El {
namespace hpd_inv {

// This approach is based upon the reordered Variant 2 algorithm from Fig. 9 in
// Bientinesi et al.'s "Families of Algorithms Related to the Inversion of
// a Symmetric Positive Definite Matrix".

template<typename Field>
void
CholeskyUVar2( Matrix<Field>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("Nonsquare matrices cannot be triangular");
    )

    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(n-k,bsize);
        const Range<Int> ind0( 0, k ), ind1( k, k+nb ), ind2( k+nb, n );

        auto A00 = A( ind0, ind0 );
        auto A01 = A( ind0, ind1 );
        auto A02 = A( ind0, ind2 );
        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );
        auto A22 = A( ind2, ind2 );

        Cholesky( UpperOrLower::UPPER, A11 );
        Trsm( LeftOrRight::RIGHT, UpperOrLower::UPPER, Orientation::NORMAL, UnitOrNonUnit::NON_UNIT, Field(1), A11, A01 );
        Trsm( LeftOrRight::LEFT, UpperOrLower::UPPER, Orientation::ADJOINT, UnitOrNonUnit::NON_UNIT, Field(1), A11, A12 );
        Herk( UpperOrLower::UPPER, Orientation::NORMAL, Base<Field>(1), A01, Base<Field>(1), A00 );
        Gemm( Orientation::NORMAL, Orientation::NORMAL, Field(-1), A01, A12, Field(1), A02 );
        Herk( UpperOrLower::UPPER, Orientation::ADJOINT, Base<Field>(-1), A12, Base<Field>(1), A22 );
        Trsm( LeftOrRight::RIGHT, UpperOrLower::UPPER, Orientation::ADJOINT, UnitOrNonUnit::NON_UNIT, Field(1), A11, A01 );
        Trsm( LeftOrRight::LEFT, UpperOrLower::UPPER, Orientation::NORMAL, UnitOrNonUnit::NON_UNIT, Field(-1), A11, A12 );
        TriangularInverse( UpperOrLower::UPPER, UnitOrNonUnit::NON_UNIT, A11 );
        Trtrmm( UpperOrLower::UPPER, A11, true );
    }
}

template<typename Field>
void
CholeskyUVar2( AbstractDistMatrix<Field>& APre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( APre.Height() != APre.Width() )
          LogicError("Nonsquare matrices cannot be triangular");
    )

    DistMatrixReadWriteProxy<Field,Field,Dist::MC,Dist::MR> AProx( APre );
    auto& A = AProx.Get();

    const Grid& g = A.Grid();
    DistMatrix<Field,Dist::STAR,Dist::STAR> A11_STAR_STAR(g);
    DistMatrix<Field,Dist::VC,  Dist::STAR> A01_VC_STAR(g);
    DistMatrix<Field,Dist::VR,  Dist::STAR> A01_VR_STAR(g);
    DistMatrix<Field,Dist::STAR,Dist::VR  > A12_STAR_VR(g);
    DistMatrix<Field,Dist::STAR,Dist::MC  > A01Trans_STAR_MC(g);
    DistMatrix<Field,Dist::MR,  Dist::STAR> A01_MR_STAR(g);
    DistMatrix<Field,Dist::STAR,Dist::MR  > A01Adj_STAR_MR(g);
    DistMatrix<Field,Dist::STAR,Dist::MR  > A12_STAR_MR(g);
    DistMatrix<Field,Dist::STAR,Dist::MC  > A12_STAR_MC(g);

    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(n-k,bsize);
        const Range<Int> ind0( 0, k ), ind1( k, k+nb ), ind2( k+nb, n );

        auto A00 = A( ind0, ind0 );
        auto A01 = A( ind0, ind1 );
        auto A02 = A( ind0, ind2 );
        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );
        auto A22 = A( ind2, ind2 );

        A11_STAR_STAR = A11;
        Cholesky( UpperOrLower::UPPER, A11_STAR_STAR );

        A01_VC_STAR.AlignWith( A00 );
        A01_VC_STAR = A01;
        LocalTrsm
        ( LeftOrRight::RIGHT, UpperOrLower::UPPER, Orientation::NORMAL, UnitOrNonUnit::NON_UNIT,
          Field(1), A11_STAR_STAR, A01_VC_STAR );

        A12_STAR_VR.AlignWith( A02 );
        A12_STAR_VR = A12;
        LocalTrsm
        ( LeftOrRight::LEFT, UpperOrLower::UPPER, Orientation::ADJOINT, UnitOrNonUnit::NON_UNIT,
          Field(1), A11_STAR_STAR, A12_STAR_VR );

        A01Trans_STAR_MC.AlignWith( A00 );
        Transpose( A01_VC_STAR, A01Trans_STAR_MC );
        A01_VR_STAR.AlignWith( A00 );
        A01_VR_STAR = A01_VC_STAR;
        A01Adj_STAR_MR.AlignWith( A00 );
        Adjoint( A01_VR_STAR, A01Adj_STAR_MR );
        LocalTrrk
        ( UpperOrLower::UPPER, Orientation::TRANSPOSE,
          Field(1), A01Trans_STAR_MC, A01Adj_STAR_MR, Field(1), A00 );

        A12_STAR_MR.AlignWith( A02 );
        A12_STAR_MR = A12_STAR_VR;
        LocalGemm
        ( Orientation::TRANSPOSE, Orientation::NORMAL,
          Field(-1), A01Trans_STAR_MC, A12_STAR_MR, Field(1), A02 );

        A12_STAR_MC.AlignWith( A22 );
        A12_STAR_MC = A12_STAR_VR;
        LocalTrrk
        ( UpperOrLower::UPPER, Orientation::ADJOINT,
          Field(-1), A12_STAR_MC, A12_STAR_MR, Field(1), A22 );

        LocalTrsm
        ( LeftOrRight::RIGHT, UpperOrLower::UPPER, Orientation::ADJOINT, UnitOrNonUnit::NON_UNIT,
          Field(1), A11_STAR_STAR, A01_VC_STAR );
        A01 = A01_VC_STAR;

        LocalTrsm
        ( LeftOrRight::LEFT, UpperOrLower::UPPER, Orientation::NORMAL, UnitOrNonUnit::NON_UNIT,
          Field(-1), A11_STAR_STAR, A12_STAR_VR );
        A12 = A12_STAR_VR;

        LocalTriangularInverse( UpperOrLower::UPPER, UnitOrNonUnit::NON_UNIT, A11_STAR_STAR );

        Trtrmm( UpperOrLower::UPPER, A11_STAR_STAR, true );
        A11 = A11_STAR_STAR;
    }
}

} // namespace hpd_inv
} // namespace El

#endif // ifndef EL_INVERSE_HPD_CHOLESKYUVAR2_HPP
