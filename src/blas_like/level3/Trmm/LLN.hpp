/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2013, The University of Texas at Austin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRMM_LLN_HPP
#define EL_TRMM_LLN_HPP

namespace El {
namespace trmm {

template<typename T>
void LocalAccumulateLLN
( Orientation orientation,
  UnitOrNonUnit diag,
  T alpha,
  const DistMatrix<T,Dist::MC,  Dist::MR  >& L,
  const DistMatrix<T,Dist::STAR,Dist::MR  >& XTrans,
        DistMatrix<T,Dist::MC,  Dist::STAR>& Z )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( L, XTrans, Z );
      if( L.Height() != L.Width() ||
          L.Height() != XTrans.Width() ||
          L.Height() != Z.Height() ||
          XTrans.Height() != Z.Width() )
          LogicError
          ("Nonconformal: \n",DimsString(L,"L"),"\n",
           DimsString(XTrans,"X'"),"\n",DimsString(Z,"Z"));
      if( XTrans.RowAlign() != L.RowAlign() ||
          Z.ColAlign() != L.ColAlign() )
          LogicError("Partial matrix distributions are misaligned");
    )
    const Int m = Z.Height();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<T> D11(g);
    const Int ratio = Max( g.Height(), g.Width() );
    for( Int k=0; k<m; k+=ratio*bsize )
    {
        const Int nb = Min(ratio*bsize,m-k);

        auto L11 = L( IR(k,k+nb), IR(k,k+nb) );
        auto L21 = L( IR(k+nb,m), IR(k,k+nb) );

        auto X1Trans = XTrans( ALL, IR(k,k+nb) );

        auto Z1 = Z( IR(k,k+nb), ALL );
        auto Z2 = Z( IR(k+nb,m), ALL );

        D11.AlignWith( L11 );
        D11 = L11;
        MakeTrapezoidal( UpperOrLower::LOWER, D11 );
        if( diag == UnitOrNonUnit::UNIT )
            FillDiagonal( D11, T(1) );
        LocalGemm( Orientation::NORMAL, orientation, alpha, D11, X1Trans, T(1), Z1 );
        LocalGemm( Orientation::NORMAL, orientation, alpha, L21, X1Trans, T(1), Z2 );
    }
}

template<typename T>
void LLNA
( UnitOrNonUnit diag,
  const AbstractDistMatrix<T>& LPre,
        AbstractDistMatrix<T>& XPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( LPre, XPre );
      if( LPre.Height() != LPre.Width() || LPre.Width() != XPre.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(LPre,"L"),"\n",DimsString(XPre,"X"))
    )
    const Int m = XPre.Height();
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<T,T,Dist::MC,Dist::MR> LProx( LPre );
    DistMatrixReadWriteProxy<T,T,Dist::MC,Dist::MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,Dist::VR,  Dist::STAR> X1_VR_STAR(g);
    DistMatrix<T,Dist::STAR,Dist::MR  > X1Trans_STAR_MR(g);
    DistMatrix<T,Dist::MC,  Dist::STAR> Z1_MC_STAR(g);

    X1_VR_STAR.AlignWith( L );
    X1Trans_STAR_MR.AlignWith( L );
    Z1_MC_STAR.AlignWith( L );

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto X1 = X( ALL, IR(k,k+nb) );

        X1_VR_STAR = X1;
        Transpose( X1_VR_STAR, X1Trans_STAR_MR );
        Z1_MC_STAR.Resize( m, nb );
        Zero( Z1_MC_STAR );
        LocalAccumulateLLN
        ( Orientation::TRANSPOSE, diag, T(1), L, X1Trans_STAR_MR, Z1_MC_STAR );
        Contract( Z1_MC_STAR, X1 );
    }
}

template<typename T>
void LLNCOld
( UnitOrNonUnit diag,
  const AbstractDistMatrix<T>& LPre,
        AbstractDistMatrix<T>& XPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( LPre, XPre );
      if( LPre.Height() != LPre.Width() || LPre.Width() != XPre.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(LPre,"L"),"\n",DimsString(XPre,"X"))
    )
    const Int m = XPre.Height();
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<T,T,Dist::MC,Dist::MR> LProx( LPre );
    DistMatrixReadWriteProxy<T,T,Dist::MC,Dist::MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,Dist::STAR,Dist::MC  > L10_STAR_MC(g);
    DistMatrix<T,Dist::STAR,Dist::STAR> L11_STAR_STAR(g);
    DistMatrix<T,Dist::STAR,Dist::VR  > X1_STAR_VR(g);
    DistMatrix<T,Dist::MR,  Dist::STAR> D1Trans_MR_STAR(g);
    DistMatrix<T,Dist::MR,  Dist::MC  > D1Trans_MR_MC(g);
    DistMatrix<T,Dist::MC,  Dist::MR  > D1(g);

    const Int kLast = LastOffset( m, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto L10 = L( IR(k,k+nb), IR(0,k)    );
        auto L11 = L( IR(k,k+nb), IR(k,k+nb) );

        auto X0 = X( IR(0,k),    ALL );
        auto X1 = X( IR(k,k+nb), ALL );

        L11_STAR_STAR = L11;
        X1_STAR_VR = X1;
        LocalTrmm( LeftOrRight::LEFT, UpperOrLower::LOWER, Orientation::NORMAL, diag, T(1), L11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;

        L10_STAR_MC.AlignWith( X0 );
        L10_STAR_MC = L10;
        D1Trans_MR_STAR.AlignWith( X1 );
        LocalGemm
        ( Orientation::TRANSPOSE, Orientation::TRANSPOSE, T(1), X0, L10_STAR_MC, D1Trans_MR_STAR );
        D1Trans_MR_MC.AlignWith( X1 );
        Contract( D1Trans_MR_STAR, D1Trans_MR_MC );
        D1.AlignWith( X1 );
        D1.Resize( nb, n );
        Zero( D1 );
        Transpose( D1Trans_MR_MC.Matrix(), D1.Matrix() );
        Axpy( T(1), D1, X1 );
    }
}

template<typename T>
void LLNC
( UnitOrNonUnit diag,
  const AbstractDistMatrix<T>& LPre,
        AbstractDistMatrix<T>& XPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( LPre, XPre );
      if( LPre.Height() != LPre.Width() || LPre.Width() != XPre.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(LPre,"L"),"\n",DimsString(XPre,"X"))
    )
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<T,T,Dist::MC,Dist::MR> LProx( LPre );
    DistMatrixReadWriteProxy<T,T,Dist::MC,Dist::MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,Dist::MC,  Dist::STAR> L21_MC_STAR(g);
    DistMatrix<T,Dist::STAR,Dist::STAR> L11_STAR_STAR(g);
    DistMatrix<T,Dist::STAR,Dist::VR  > X1_STAR_VR(g);
    DistMatrix<T,Dist::MR,  Dist::STAR> X1Trans_MR_STAR(g);

    const Int kLast = LastOffset( m, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto L11 = L( IR(k,k+nb), IR(k,k+nb) );
        auto L21 = L( IR(k+nb,m), IR(k,k+nb) );

        auto X1 = X( IR(k,k+nb), ALL );
        auto X2 = X( IR(k+nb,m), ALL );

        L21_MC_STAR.AlignWith( X2 );
        L21_MC_STAR = L21;
        X1Trans_MR_STAR.AlignWith( X2 );
        Transpose( X1, X1Trans_MR_STAR );
        LocalGemm
        ( Orientation::NORMAL, Orientation::TRANSPOSE, T(1), L21_MC_STAR, X1Trans_MR_STAR, T(1), X2 );

        L11_STAR_STAR = L11;
        X1_STAR_VR.AlignWith( X1 );
        Transpose( X1Trans_MR_STAR, X1_STAR_VR );
        LocalTrmm( LeftOrRight::LEFT, UpperOrLower::LOWER, Orientation::NORMAL, diag, T(1), L11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
    }
}

// Left Lower Normal (Non)Unit Trmm
//   X := tril(L)  X, or
//   X := trilu(L) X
template<typename T>
void LLN
( UnitOrNonUnit diag,
  const AbstractDistMatrix<T>& L,
        AbstractDistMatrix<T>& X )
{
    EL_DEBUG_CSE
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Width() )
        LLNA( diag, L, X );
    else
        LLNC( diag, L, X );
}

} // namespace trmm
} // namespace El

#endif // ifndef EL_TRMM_LLN_HPP
