/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_ONES_HPP
#define EL_ONES_HPP

namespace El {

template<typename T> 
void MakeOnes( Matrix<T>& A );

template<typename T,Dist U,Dist V>
void MakeOnes( DistMatrix<T,U,V>& A );

template<typename T,Dist U,Dist V>
void MakeOnes( BlockDistMatrix<T,U,V>& A );

template<typename T>
void Ones( Matrix<T>& A, Int m, Int n );

template<typename T,Dist U,Dist V>
void Ones( DistMatrix<T,U,V>& A, Int m, Int n );

template<typename T,Dist U,Dist V>
void Ones( BlockDistMatrix<T,U,V>& A, Int m, Int n );

template<typename T>
void Ones( AbstractDistMatrix<T>& A, Int m, Int n );

template<typename T>
void Ones( AbstractBlockDistMatrix<T>& A, Int m, Int n );

} // namespace El

#endif // ifndef EL_ONES_HPP