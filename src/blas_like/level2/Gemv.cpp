/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El/blas_like/level1/decl.hpp"
#include "El/blas_like/level2.hpp"
#include "El/core/Element/Complex/decl.hpp"
#include "El/core/Element/Complex/impl.hpp"
#include "El/core/imports/blas.hpp"
#include "El/core/imports/scalapack.hpp"
#include "El/core/Matrix/decl.hpp"

#include "./Gemv/Normal.hpp"
#include "./Gemv/Transpose.hpp"

namespace El
{

template<typename T>
void Gemv
(Orientation orientation,
  T alpha, const Matrix<T>& A,
           const Matrix<T>& x,
  T beta,        Matrix<T>& y)
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if((x.Height() != 1 && x.Width() != 1) ||
          (y.Height() != 1 && y.Width() != 1))
          LogicError
          ("Nonconformal: \n",DimsString(x,"x"),"\n",DimsString(y,"y"));
      const Int xLength = (x.Width()==1 ? x.Height() : x.Width());
      const Int yLength = (y.Width()==1 ? y.Height() : y.Width());
      if(orientation == Orientation::NORMAL)
      {
          if(A.Height() != yLength || A.Width() != xLength)
              LogicError
              ("Nonconformal: \n",DimsString(A,"A"),"\n",
               DimsString(x,"x"),"\n",DimsString(y,"y"));
      }
      else
      {
          if(A.Width() != yLength || A.Height() != xLength)
              LogicError
              ("Nonconformal: \n",DimsString(A,"A"),"\n",
               DimsString(x,"x"),"\n",DimsString(y,"y"));
      }
   )
    const char transChar = OrientationToChar(orientation);
    const Int m = A.Height();
    const Int n = A.Width();
    const Int xDim = (transChar == 'N' ? n : m);
    const Int yDim = (transChar == 'N' ? m : n);
    const Int incx = (x.Width()==1 ? 1 : x.LDim());
    const Int incy = (y.Width()==1 ? 1 : y.LDim());
    if(xDim != 0)
    {
        if(yDim != 0)
        {
            blas::Gemv
            (transChar, m, n,
              alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx,
              beta,  y.Buffer(), incy);
        }
    }
    else
    {
        y *= beta;
    }
}

template<typename T>
void Gemv
(Orientation orientation,
  T alpha, const Matrix<T>& A,
           const Matrix<T>& x,
                 Matrix<T>& y)
{
    EL_DEBUG_CSE
    if(orientation == Orientation::NORMAL)
        y.Resize(A.Height(), 1);
    else
        y.Resize(A.Width(), 1);
    Zero(y);
    Gemv(orientation, alpha, A, x, T(0), y);
}

template<typename T>
void Gemv
(Orientation orientation,
  T alpha, const AbstractDistMatrix<T>& A,
           const AbstractDistMatrix<T>& x,
  T beta,        AbstractDistMatrix<T>& y)
{
    EL_DEBUG_CSE
    if(orientation == Orientation::NORMAL)
        gemv::Normal(alpha, A, x, beta, y);
    else
        gemv::Transpose(orientation, alpha, A, x, beta, y);
}

template<typename T>
void Gemv
(Orientation orientation,
  T alpha, const AbstractDistMatrix<T>& A,
           const AbstractDistMatrix<T>& x,
                 AbstractDistMatrix<T>& y)
{
    EL_DEBUG_CSE
    y.AlignWith(A);
    if(orientation == Orientation::NORMAL)
        y.Resize(A.Height(), 1);
    else
        y.Resize(A.Width(), 1);
    Zero(y);
    Gemv(orientation, alpha, A, x, T(0), y);
}

template<typename T>
void LocalGemv
(Orientation orientation,
  T alpha, const AbstractDistMatrix<T>& A,
           const AbstractDistMatrix<T>& x,
  T beta,        AbstractDistMatrix<T>& y)
{
    EL_DEBUG_CSE
    // TODO(poulson): Add error checking here
    Gemv
    (orientation ,
      alpha, A.LockedMatrix(), x.LockedMatrix(),
      beta,                    y.Matrix());
}

namespace gemv {

template<typename T,typename=EnableIf<IsBlasScalar<T>>>
void ScaLAPACKHelper
(Orientation orientation,
  T alpha, const DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& A,
           const DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& x,
  T beta,        DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& y)
{
    AssertScaLAPACKSupport();
#ifdef EL_HAVE_SCALAPACK
    const Int m = A.Height();
    const Int n = A.Width();
    const char orientChar = OrientationToChar(orientation);

    auto descA = FillDesc(A);
    auto descx = FillDesc(x);
    auto descy = FillDesc(y);
    pblas::Gemv
    (orientChar, m, n,
      alpha,
      A.LockedBuffer(), descA.data(),
      x.LockedBuffer(), descx.data(), 1,
      beta,
      y.Buffer(),       descy.data(), 1);
#endif
}

template<typename T,typename=DisableIf<IsBlasScalar<T>>,typename=void>
void ScaLAPACKHelper
(Orientation orientation,
  T alpha, const DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& A,
           const DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& x,
  T beta,        DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& y)
{
    LogicError("ScaLAPACK does not support this datatype");
}

} // namespace gemv

template<typename T>
void Gemv
(Orientation orientation,
  T alpha, const DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& A,
           const DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& x,
  T beta,        DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& y)
{
    EL_DEBUG_CSE
    gemv::ScaLAPACKHelper(orientation, alpha, A, x, beta, y);
}

template<>
void Gemv
(Orientation orientation,
  Int alpha, const DistMatrix<Int,Dist::MC,Dist::MR,DistWrap::BLOCK>& A,
             const DistMatrix<Int,Dist::MC,Dist::MR,DistWrap::BLOCK>& x,
  Int beta,        DistMatrix<Int,Dist::MC,Dist::MR,DistWrap::BLOCK>& y)
{
    EL_DEBUG_CSE
    LogicError("ScaLAPACK does not support integer data");
}
#define EL_NO_INT_PROTO

#ifdef HYDROGEN_HAVE_QUADMATH
template<>
void Gemv
(Orientation orientation,
  Quad alpha, const DistMatrix<Quad,Dist::MC,Dist::MR,DistWrap::BLOCK>& A,
              const DistMatrix<Quad,Dist::MC,Dist::MR,DistWrap::BLOCK>& x,
  Quad beta,        DistMatrix<Quad,Dist::MC,Dist::MR,DistWrap::BLOCK>& y)
{
    EL_DEBUG_CSE
    LogicError("ScaLAPACK does not support quad-precision data");
}

template<>
void Gemv
(Orientation orientation,
  Complex<Quad> alpha, const DistMatrix<Complex<Quad>,Dist::MC,Dist::MR,DistWrap::BLOCK>& A,
                       const DistMatrix<Complex<Quad>,Dist::MC,Dist::MR,DistWrap::BLOCK>& x,
  Complex<Quad> beta,        DistMatrix<Complex<Quad>,Dist::MC,Dist::MR,DistWrap::BLOCK>& y)
{
    EL_DEBUG_CSE
    LogicError("ScaLAPACK does not support quad-precision data");
}
#endif // ifdef HYDROGEN_HAVE_QUADMATH

template<typename T>
void Gemv
(Orientation orientation,
  T alpha, const DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& A,
           const DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& x,
                 DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& y)
{
    EL_DEBUG_CSE
    y.AlignWith(A);
    if(orientation == Orientation::NORMAL)
        y.Resize(A.Height(), 1);
    else
        y.Resize(A.Width(), 1);
    Zero(y);
    Gemv(orientation, alpha, A, x, T(0), y);
}

#define PROTO(T) \
  template void Gemv \
  (Orientation orientation, \
    T alpha, const Matrix<T>& A, \
             const Matrix<T>& x, \
    T beta,        Matrix<T>& y); \
  template void Gemv \
  (Orientation orientation, \
    T alpha, const Matrix<T>& A, \
             const Matrix<T>& x, \
                   Matrix<T>& y); \
  template void Gemv \
  (Orientation orientation, \
    T alpha, const AbstractDistMatrix<T>& A, \
             const AbstractDistMatrix<T>& x, \
    T beta,        AbstractDistMatrix<T>& y); \
  template void Gemv \
  (Orientation orientation, \
    T alpha, const AbstractDistMatrix<T>& A, \
             const AbstractDistMatrix<T>& x, \
                   AbstractDistMatrix<T>& y); \
  template void Gemv \
  (Orientation orientation, \
    T alpha, const DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& A, \
             const DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& x, \
    T beta,        DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& y); \
  template void Gemv \
  (Orientation orientation, \
    T alpha, const DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& A, \
             const DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& x, \
                   DistMatrix<T,Dist::MC,Dist::MR,DistWrap::BLOCK>& y); \
  template void LocalGemv \
  (Orientation orientation, \
    T alpha, const AbstractDistMatrix<T>& A, \
             const AbstractDistMatrix<T>& x, \
    T beta,        AbstractDistMatrix<T>& y);

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
