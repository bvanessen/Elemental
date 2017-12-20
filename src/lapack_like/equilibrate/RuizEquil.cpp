/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {

template<typename Real>
Real DampScaling( const Real& alpha )
{
    const Real tol = Pow(limits::Epsilon<Real>(),Real(0.33));
    if( alpha == Real(0) )
        return 1;
    else
        return Max(alpha,tol);
}

template<typename Field>
void RuizEquil
( Matrix<Field>& A,
  Matrix<Base<Field>>& dRow,
  Matrix<Base<Field>>& dCol,
  bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Ones( dRow, m, 1 );
    Ones( dCol, n, 1 );

    // TODO(poulson): Expose these as control parameters
    // For now, simply hard-code the number of iterations
    const Int maxIter = 4;

    Matrix<Real> rowScale, colScale;
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Rescale the columns
        // -------------------
        ColumnMaxNorms( A, colScale );
        EntrywiseMap( colScale, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LeftOrRight::LEFT, Orientation::NORMAL, colScale, dCol );
        DiagonalSolve( LeftOrRight::RIGHT, Orientation::NORMAL, colScale, A );

        // Rescale the rows
        // ----------------
        RowMaxNorms( A, rowScale );
        EntrywiseMap( rowScale, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LeftOrRight::LEFT, Orientation::NORMAL, rowScale, dRow );
        DiagonalSolve( LeftOrRight::LEFT, Orientation::NORMAL, rowScale, A );
    }
    SetIndent( indent );
}

template<typename Field>
void RuizEquil
( AbstractDistMatrix<Field>& APre,
  AbstractDistMatrix<Base<Field>>& dRowPre,
  AbstractDistMatrix<Base<Field>>& dColPre,
  bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    ElementalProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;

    DistMatrixReadWriteProxy<Field,Field,Dist::MC,Dist::MR> AProx( APre, control );
    DistMatrixWriteProxy<Real,Real,Dist::MC,Dist::STAR> dRowProx( dRowPre, control );
    DistMatrixWriteProxy<Real,Real,Dist::MR,Dist::STAR> dColProx( dColPre, control );
    auto& A = AProx.Get();
    auto& dRow = dRowProx.Get();
    auto& dCol = dColProx.Get();

    const Int m = A.Height();
    const Int n = A.Width();
    Ones( dRow, m, 1 );
    Ones( dCol, n, 1 );

    // TODO(poulson): Expose these as control parameters
    // For now, simply hard-code the number of iterations
    const Int maxIter = 4;

    DistMatrix<Real,Dist::MC,Dist::STAR> rowScale(A.Grid());
    DistMatrix<Real,Dist::MR,Dist::STAR> colScale(A.Grid());
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Rescale the columns
        // -------------------
        ColumnMaxNorms( A, colScale );
        EntrywiseMap( colScale, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LeftOrRight::LEFT, Orientation::NORMAL, colScale, dCol );
        DiagonalSolve( LeftOrRight::RIGHT, Orientation::NORMAL, colScale, A );

        // Rescale the rows
        // ----------------
        RowMaxNorms( A, rowScale );
        EntrywiseMap( rowScale, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LeftOrRight::LEFT, Orientation::NORMAL, rowScale, dRow );
        DiagonalSolve( LeftOrRight::LEFT, Orientation::NORMAL, rowScale, A );
    }
    SetIndent( indent );
}

template<typename Field>
void StackedRuizEquil
( Matrix<Field>& A,
  Matrix<Field>& B,
  Matrix<Base<Field>>& dRowA,
  Matrix<Base<Field>>& dRowB,
  Matrix<Base<Field>>& dCol,
  bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();
    Ones( dRowA, mA, 1 );
    Ones( dRowB, mB, 1 );
    Ones( dCol, n, 1 );

    // TODO(poulson): Expose these as control parameters
    // For now, simply hard-code the number of iterations
    const Int maxIter = 4;

    Matrix<Real> rowScale, colScale, colScaleB;
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Rescale the columns
        // -------------------
        ColumnMaxNorms( A, colScale );
        ColumnMaxNorms( B, colScaleB );
        for( Int j=0; j<n; ++j )
            colScale(j) = Max(colScale(j),colScaleB(j));
        EntrywiseMap( colScale, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LeftOrRight::LEFT, Orientation::NORMAL, colScale, dCol );
        DiagonalSolve( LeftOrRight::RIGHT, Orientation::NORMAL, colScale, A );
        DiagonalSolve( LeftOrRight::RIGHT, Orientation::NORMAL, colScale, B );

        // Rescale the rows
        // ----------------
        RowMaxNorms( A, rowScale );
        EntrywiseMap( rowScale, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LeftOrRight::LEFT, Orientation::NORMAL, rowScale, dRowA );
        DiagonalSolve( LeftOrRight::LEFT, Orientation::NORMAL, rowScale, A );

        RowMaxNorms( B, rowScale );
        EntrywiseMap( rowScale, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LeftOrRight::LEFT, Orientation::NORMAL, rowScale, dRowB );
        DiagonalSolve( LeftOrRight::LEFT, Orientation::NORMAL, rowScale, B );
    }
    SetIndent( indent );
}

template<typename Field>
void StackedRuizEquil
( AbstractDistMatrix<Field>& APre,
  AbstractDistMatrix<Field>& BPre,
  AbstractDistMatrix<Base<Field>>& dRowAPre,
  AbstractDistMatrix<Base<Field>>& dRowBPre,
  AbstractDistMatrix<Base<Field>>& dColPre,
  bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    ElementalProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;

    DistMatrixReadWriteProxy<Field,Field,Dist::MC,Dist::MR> AProx( APre, control );
    DistMatrixReadWriteProxy<Field,Field,Dist::MC,Dist::MR> BProx( BPre, control );
    DistMatrixWriteProxy<Real,Real,Dist::MC,Dist::STAR> dRowAProx( dRowAPre, control );
    DistMatrixWriteProxy<Real,Real,Dist::MC,Dist::STAR> dRowBProx( dRowBPre, control );
    DistMatrixWriteProxy<Real,Real,Dist::MR,Dist::STAR> dColProx( dColPre, control );
    auto& A = AProx.Get();
    auto& B = BProx.Get();
    auto& dRowA = dRowAProx.Get();
    auto& dRowB = dRowBProx.Get();
    auto& dCol = dColProx.Get();

    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();
    const Int nLocal = A.LocalWidth();
    Ones( dRowA, mA, 1 );
    Ones( dRowB, mB, 1 );
    Ones( dCol, n, 1 );

    // TODO(poulson): Expose these as control parameters
    // For now, simply hard-code the number of iterations
    const Int maxIter = 4;

    DistMatrix<Real,Dist::MC,Dist::STAR> rowScale(A.Grid());
    DistMatrix<Real,Dist::MR,Dist::STAR> colScale(A.Grid()), colScaleB(B.Grid());
    auto& colScaleLoc = colScale.Matrix();
    auto& colScaleBLoc = colScaleB.Matrix();
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Rescale the columns
        // -------------------
        ColumnMaxNorms( A, colScale );
        ColumnMaxNorms( B, colScaleB );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            colScaleLoc(jLoc) =
              Max(colScaleLoc(jLoc),colScaleBLoc(jLoc));
        EntrywiseMap( colScale, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LeftOrRight::LEFT, Orientation::NORMAL, colScale, dCol );
        DiagonalSolve( LeftOrRight::RIGHT, Orientation::NORMAL, colScale, A );
        DiagonalSolve( LeftOrRight::RIGHT, Orientation::NORMAL, colScale, B );

        // Rescale the rows
        // ----------------
        RowMaxNorms( A, rowScale );
        EntrywiseMap( rowScale, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LeftOrRight::LEFT, Orientation::NORMAL, rowScale, dRowA );
        DiagonalSolve( LeftOrRight::LEFT, Orientation::NORMAL, rowScale, A );

        RowMaxNorms( B, rowScale );
        EntrywiseMap( rowScale, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LeftOrRight::LEFT, Orientation::NORMAL, rowScale, dRowB );
        DiagonalSolve( LeftOrRight::LEFT, Orientation::NORMAL, rowScale, B );
    }
    SetIndent( indent );
}


#define PROTO(Field) \
  template void RuizEquil \
  ( Matrix<Field>& A, \
    Matrix<Base<Field>>& dRow, \
    Matrix<Base<Field>>& dCol, \
    bool progress ); \
  template void RuizEquil \
  ( AbstractDistMatrix<Field>& A, \
    AbstractDistMatrix<Base<Field>>& dRow, \
    AbstractDistMatrix<Base<Field>>& dCol, \
    bool progress ); \
  template void StackedRuizEquil \
  ( Matrix<Field>& A, \
    Matrix<Field>& B, \
    Matrix<Base<Field>>& dRowA, \
    Matrix<Base<Field>>& dRowB, \
    Matrix<Base<Field>>& dCol, \
    bool progress ); \
  template void StackedRuizEquil \
  ( AbstractDistMatrix<Field>& A, \
    AbstractDistMatrix<Field>& B, \
    AbstractDistMatrix<Base<Field>>& dRowA, \
    AbstractDistMatrix<Base<Field>>& dRowB, \
    AbstractDistMatrix<Base<Field>>& dCol, \
    bool progress );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
