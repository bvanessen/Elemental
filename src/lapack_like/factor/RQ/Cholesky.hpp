/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_RQ_CHOLESKY_HPP
#define EL_RQ_CHOLESKY_HPP

namespace El {
namespace rq {

// NOTE: This version is designed for short-fat matrices and is much less
//       numerically stable than Householder-based RQ factorizations
//

template<typename F>
void Cholesky( Matrix<F>& A, Matrix<F>& R )
{
    EL_DEBUG_CSE
    if( A.Height() > A.Width() )
        LogicError("A A^H will be singular");
    Herk( UpperOrLower::UPPER, Orientation::NORMAL, Base<F>(1), A, R );
    El::ReverseCholesky( UpperOrLower::UPPER, R );
    Trsm( LeftOrRight::LEFT, UpperOrLower::UPPER, Orientation::NORMAL, UnitOrNonUnit::NON_UNIT, F(1), R, A );
}

template<typename F>
void Cholesky( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& RPre )
{
    EL_DEBUG_CSE
    const Int m = APre.Height();
    const Int n = APre.Width();
    if( m > n )
        LogicError("A A^H will be singular");

    DistMatrixReadWriteProxy<F,F,Dist::VR,Dist::STAR> AProx( APre );
    DistMatrixWriteProxy<F,F,Dist::STAR,Dist::STAR> RProx( RPre );
    auto& A = AProx.Get();
    auto& R = RProx.Get();

    Zeros( R, m, m );
    Herk( UpperOrLower::UPPER, Orientation::NORMAL, Base<F>(1), A.Matrix(), Base<F>(0), R.Matrix() );
    El::AllReduce( R, A.RowComm() );
    El::ReverseCholesky( UpperOrLower::UPPER, R.Matrix() );
    Trsm( LeftOrRight::LEFT, UpperOrLower::UPPER, Orientation::NORMAL, UnitOrNonUnit::NON_UNIT, F(1), R.Matrix(), A.Matrix() );
}

} // namespace rq
} // namespace El

#endif // ifndef EL_RQ_CHOLESKY_HPP
