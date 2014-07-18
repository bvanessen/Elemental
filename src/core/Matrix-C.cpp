/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El-C.h"
using namespace El;

extern "C" {

#define MATRIX_CONSTRUCT(SIG,SIGBASE,T) \
  /* Matrix<T>::Matrix() */ \
  ElError ElMatrixCreate_ ## SIG ( ElMatrix_ ## SIG * A ) \
  { EL_TRY( *A = Reinterpret( new Matrix<T> ) ) } \
  /* Matrix<T>::~Matrix() */ \
  ElError ElMatrixDestroy_ ## SIG ( ElConstMatrix_ ## SIG A ) \
  { EL_TRY( delete Reinterpret(A) ) }

#define MATRIX_RECONFIG(SIG,SIGBASE,T) \
  /* void Matrix<T>::Empty() */ \
  ElError ElMatrixEmpty_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( Reinterpret(A)->Empty() ) } \
  /* void Matrix<T>::Resize( Int height, Int width ) */ \
  ElError ElMatrixResize_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt height, ElInt width ) \
  { EL_TRY( Reinterpret(A)->Resize(height,width) ) } \
  /* void Matrix<T>::Resize( Int height, Int width, Int ldim ) */ \
  ElError ElMatrixResizeWithLDim_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt height, ElInt width, ElInt ldim ) \
  { EL_TRY( Reinterpret(A)->Resize(height,width,ldim) ) } \
  /* void Matrix<T>::Attach( Int height, Int width, T* buffer, Int ldim ) */ \
  ElError ElMatrixAttach_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElInt height, ElInt width, CREFLECT(T)* buffer, ElInt ldim ) \
  { EL_TRY\
    ( Reinterpret(A)->Attach(height,width,Reinterpret(buffer),ldim) ) } \
  /* void Matrix<T>::LockedAttach
     ( Int height, Int width, const T* buffer, Int ldim ) */ \
  ElError ElMatrixLockedAttach_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElInt height, ElInt width, const CREFLECT(T)* buffer, ElInt ldim ) \
  { EL_TRY \
    ( Reinterpret(A)->LockedAttach \
      (height,width,Reinterpret(buffer),ldim) ) } \
  /* void Matrix<T>::Control( Int height, Int width, T* buffer, Int ldim ) */ \
  ElError ElMatrixControl_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt height, ElInt width, CREFLECT(T)* buffer, \
    ElInt ldim ) \
  { EL_TRY \
    ( Reinterpret(A)->Control(height,width,Reinterpret(buffer),ldim) ) } \
  /* B = A */ \
  ElError ElMatrixCopy_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( *Reinterpret(B) = *Reinterpret(A) ) }

#define MATRIX_BASIC(SIG,SIGBASE,T) \
  /* Int Matrix<T>::Height() const */ \
  ElError ElMatrixHeight_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElInt* height ) \
  { EL_TRY( *height = Reinterpret(A)->Height() ) } \
  /* Int Matrix<T>::Width() const */ \
  ElError ElMatrixWidth_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElInt* width ) \
  { EL_TRY( *width = Reinterpret(A)->Width() ) } \
  /* Int Matrix<T>::LDim() const */ \
  ElError ElMatrixLDim_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElInt* ldim ) \
  { EL_TRY( *ldim = Reinterpret(A)->LDim() ) } \
  /* Int Matrix<T>::MemorySize() const */ \
  ElError ElMatrixMemorySize_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElInt* mem ) \
  { EL_TRY( *mem = Reinterpret(A)->MemorySize() ) } \
  /* Int Matrix<T>::DiagonalLength( Int offset ) const */ \
  ElError ElMatrixDiagonalLength_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElInt offset, ElInt* length ) \
  { EL_TRY( *length = Reinterpret(A)->DiagonalLength(offset) ) } \
  /* T* Matrix<T>::Buffer() */ \
  ElError ElMatrixBuffer_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T)** buffer ) \
  { EL_TRY( *buffer = Reinterpret(Reinterpret(A)->Buffer()) ) } \
  /* const T* Matrix<T>::LockedBuffer() const */ \
  ElError ElMatrixLockedBuffer_ ## SIG \
  ( ElConstMatrix_ ## SIG A, const CREFLECT(T)** buffer ) \
  { EL_TRY( *buffer = Reinterpret(Reinterpret(A)->LockedBuffer()) ) } \
  /* bool Matrix<T>::Viewing() const */ \
  ElError ElMatrixViewing_ ## SIG \
  ( ElConstMatrix_ ## SIG A, bool* viewing ) \
  { EL_TRY( *viewing = Reinterpret(A)->Viewing() ) } \
  /* bool Matrix<T>::FixedSize() const */ \
  ElError ElMatrixFixedSize_ ## SIG \
  ( ElConstMatrix_ ## SIG A, bool* fixedSize ) \
  { EL_TRY( *fixedSize = Reinterpret(A)->FixedSize() ) } \
  /* bool Matrix<T>::Locked() const */ \
  ElError ElMatrixLocked_ ## SIG \
  ( ElConstMatrix_ ## SIG A, bool* locked ) \
  { EL_TRY( *locked = Reinterpret(A)->Locked() ) }

#define MATRIX_SINGLEENTRY(SIG,SIGBASE,T) \
  /* T Matrix<T>::Get( Int i, Int j ) const */ \
  ElError ElMatrixGet_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElInt i, ElInt j, CREFLECT(T)* val ) \
  { EL_TRY( *val = Reinterpret(Reinterpret(A)->Get(i,j)) ) } \
  /* void Matrix<T>::Set( Int i, Int j, T alpha ) */ \
  ElError ElMatrixSet_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt i, ElInt j, CREFLECT(T) alpha ) \
  { EL_TRY( Reinterpret(A)->Set(i,j,Reinterpret(alpha)) ) } \
  /* void Matrix<T>::Update( Int i, Int j, T alpha ) */ \
  ElError ElMatrixUpdate_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt i, ElInt j, CREFLECT(T) alpha ) \
  { EL_TRY( Reinterpret(A)->Update(i,j,Reinterpret(alpha)) ) }

#define MATRIX_SINGLEENTRY_COMPLEX(SIG,SIGBASE,T) \
  /* Base<T> Matrix<T>::GetRealPart( Int i, Int j ) const */ \
  ElError ElMatrixGetRealPart_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElInt i, ElInt j, Base<T>* val ) \
  { EL_TRY( *val = Reinterpret(A)->GetRealPart(i,j) ) } \
  /* Base<T> Matrix<T>::GetImagPart( Int i, Int j ) const */ \
  ElError ElMatrixGetImagPart_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElInt i, ElInt j, Base<T>* val ) \
  { EL_TRY( *val = Reinterpret(A)->GetImagPart(i,j) ) } \
  /* void Matrix<T>::SetRealPart( Int i, Int j, Base<T> alpha ) */ \
  ElError ElMatrixSetRealPart_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt i, ElInt j, Base<T> alpha ) \
  { EL_TRY( Reinterpret(A)->SetRealPart(i,j,alpha) ) } \
  /* void Matrix<T>::SetImagPart( Int i, Int j, Base<T> alpha ) */ \
  ElError ElMatrixSetImagPart_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt i, ElInt j, Base<T> alpha ) \
  { EL_TRY( Reinterpret(A)->SetImagPart(i,j,alpha) ) } \
  /* void Matrix<T>::UpdateRealPart( Int i, Int j, Base<T> alpha ) */ \
  ElError ElMatrixUpdateRealPart_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt i, ElInt j, Base<T> alpha ) \
  { EL_TRY( Reinterpret(A)->UpdateRealPart(i,j,alpha) ) } \
  /* void Matrix<T>::UpdateImagPart( Int i, Int j, Base<T> alpha ) */ \
  ElError ElMatrixUpdateImagPart_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt i, ElInt j, Base<T> alpha ) \
  { EL_TRY( Reinterpret(A)->UpdateImagPart(i,j,alpha) ) } \
  /* void Matrix<T>::MakeReal( Int i, Int j ) */ \
  ElError ElMatrixMakeReal_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt i, ElInt j ) \
  { EL_TRY( Reinterpret(A)->MakeReal(i,j) ) } \
  /* void Matrix<T>::Conjugate( Int i, Int j ) */ \
  ElError ElMatrixConjugate_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt i, ElInt j ) \
  { EL_TRY( Reinterpret(A)->Conjugate(i,j) ) }

#define MATRIX_DIAGONAL(SIG,SIGBASE,T) \
  /* Matrix<T> Matrix<T>::GetDiagonal( Int offset ) const */ \
  ElError ElMatrixGetDiagonal_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElInt offset, ElMatrix_ ## SIG *d ) \
  { EL_TRY( auto dPtr = new Matrix<T>; \
            Reinterpret(A)->GetDiagonal( *dPtr, offset ); \
            *d = Reinterpret(dPtr) ) } \
  /* void Matrix<T>::SetDiagonal( const Matrix<T>& d, Int offset ) */ \
  ElError ElMatrixSetDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_ ## SIG d, ElInt offset ) \
  { EL_TRY \
    ( Reinterpret(A)->SetDiagonal( Reinterpret(d), offset ) ) } \
  /* void Matrix<T>::UpdateDiagonal( const Matrix<T>& d, Int offset ) */ \
  ElError ElMatrixUpdateDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_ ## SIG d, ElInt offset ) \
  { EL_TRY \
    ( Reinterpret(A)->UpdateDiagonal( Reinterpret(d), offset ) ) }

#define MATRIX_DIAGONAL_COMPLEX(SIG,SIGBASE,T) \
  /* Matrix<Base<T>> Matrix<T>::GetRealPartOfDiagonal( Int offset ) const */ \
  ElError ElMatrixGetRealPartOfDiagonal_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElInt offset, \
    ElMatrix_ ## SIGBASE *d ) \
  { EL_TRY( auto dPtr = new Matrix<Base<T>>; \
            Reinterpret(A)->GetRealPartOfDiagonal( *dPtr, offset ); \
            *d = Reinterpret(dPtr) ) } \
  /* Matrix<Base<T>> Matrix<T>::GetImagPartOfDiagonal( Int offset ) const */ \
  ElError ElMatrixGetImagPartOfDiagonal_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElInt offset, \
    ElMatrix_ ## SIGBASE *d ) \
  { EL_TRY( auto dPtr = new Matrix<Base<T>>; \
            Reinterpret(A)->GetImagPartOfDiagonal( *dPtr, offset ); \
            *d = Reinterpret(dPtr) ) } \
  /* void Matrix<T>::SetRealPartOfDiagonal \
     ( const Matrix<Base<T>>& d, Int offset ) */ \
  ElError ElMatrixSetRealPartOfDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIGBASE d, ElInt offset ) \
  { EL_TRY( Reinterpret(A)->SetRealPartOfDiagonal \
            (Reinterpret(d),offset) ) } \
  /* void Matrix<T>::SetImagPartOfDiagonal \
     ( const Matrix<Base<T>>& d, Int offset ) */ \
  ElError ElMatrixSetImagPartOfDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIGBASE d, ElInt offset ) \
  { EL_TRY( Reinterpret(A)->SetImagPartOfDiagonal \
            (Reinterpret(d),offset) ) } \
  /* void Matrix<T>::UpdateRealPartOfDiagonal \
     ( const Matrix<Base<T>>& d, Int offset ) */ \
  ElError ElMatrixUpdateRealPartOfDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIGBASE d, ElInt offset ) \
  { EL_TRY( Reinterpret(A)->UpdateRealPartOfDiagonal \
            (Reinterpret(d),offset) ) } \
  /* void Matrix<T>::UpdateImagPartOfDiagonal \
     ( const Matrix<Base<T>>& d, Int offset ) */ \
  ElError ElMatrixUpdateImagPartOfDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIGBASE d, ElInt offset ) \
  { EL_TRY( Reinterpret(A)->UpdateImagPartOfDiagonal \
            (Reinterpret(d),offset) ) } \
  /* void Matrix<T>::MakeDiagonalReal( Int offset ) */ \
  ElError ElMatrixMakeDiagonalReal_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( Reinterpret(A)->MakeDiagonalReal(offset) ) } \
  /* void Matrix<T>::ConjugateDiagonal( Int offset ) */ \
  ElError ElMatrixConjugateDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( Reinterpret(A)->ConjugateDiagonal(offset) ) }

#define MATRIX_SUBMATRIX(SIG,SIGBASE,T) \
  /* Matrix<T> Matrix<T>::GetSubmatrix
     ( const std::vector<Int>& rowInds, \
       const std::vector<Int>& colInds ) const */ \
  ElError ElMatrixGetSubmatrix_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElInt numRowInds, const ElInt* rowInds, \
    ElInt numColInds, const ElInt* colInds, ElMatrix_ ## SIG *ASub ) \
  { EL_TRY( std::vector<Int> rowIndVec( rowInds, rowInds+numRowInds ); \
            std::vector<Int> colIndVec( colInds, colInds+numColInds ); \
            auto ASubPtr = new Matrix<T>; \
            Reinterpret(A)->GetSubmatrix( rowIndVec, colIndVec, *ASubPtr ); \
            *ASub = Reinterpret(ASubPtr) ) }

#define C_PROTO(SIG,SIGBASE,T) \
  MATRIX_CONSTRUCT(SIG,SIGBASE,T) \
  MATRIX_RECONFIG(SIG,SIGBASE,T) \
  MATRIX_BASIC(SIG,SIGBASE,T) \
  MATRIX_SINGLEENTRY(SIG,SIGBASE,T) \
  MATRIX_DIAGONAL(SIG,SIGBASE,T) \
  MATRIX_SUBMATRIX(SIG,SIGBASE,T)

#define C_PROTO_COMPLEX(SIG,SIGBASE,T) \
  C_PROTO(SIG,SIGBASE,T) \
  MATRIX_SINGLEENTRY_COMPLEX(SIG,SIGBASE,T) \
  MATRIX_DIAGONAL_COMPLEX(SIG,SIGBASE,T)

#include "El/macros/CInstantiate.h"

} // extern "C"
