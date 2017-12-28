/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#define EL_INSTANTIATE_CORE

#include "El/config.h"

#include "El/core/imports/valgrind.hpp"
#include "El/core/imports/omp.hpp"
#include "El/core/imports/qd.hpp"
#include "El/core/imports/mpfr.hpp"
#include "El/core/imports/qt5.hpp"

#include "El/core/Element/decl.hpp"
#include "El/core/Serialize.hpp"
#include "El/core/imports/mpi.hpp"
#include "El/core/imports/choice.hpp"
#include "El/core/imports/mpi_choice.hpp"
#include "El/core/environment/decl.hpp"

#include "El/core/Timer.hpp"
#include "El/core/indexing/decl.hpp"
#include "El/core/imports/blas.hpp"
#include "El/core/imports/lapack.hpp"
#include "El/core/imports/flame.hpp"
#include "El/core/imports/mkl.hpp"
#include "El/core/imports/openblas.hpp"
#include "El/core/imports/pmrrr.hpp"
#include "El/core/imports/scalapack.hpp"

#include "El/core/limits.hpp"

#include "El/core/Memory.hpp"

#include "El/core/Matrix/decl.hpp"
#include "El/core/DistMap/decl.hpp"
#include "El/core/View/decl.hpp"
#include "El/blas_like/level1/decl.hpp"

#include "El/core/Matrix/impl.hpp"
#include "El/core/Grid.hpp"
#include "El/core/DistMatrix.hpp"
#include "El/core/Proxy.hpp"

// Implement the intertwined parts of the library
#include "El/core/Element/impl.hpp"
#include "El/core/environment/impl.hpp"
#include "El/core/indexing/impl.hpp"

// Declare and implement the decoupled parts of the core of the library
// (perhaps these should be moved into their own directory?)
#include "El/core/View/impl.hpp"
#include "El/core/FlamePart.hpp"
#include "El/core/random/decl.hpp"
#include "El/core/random/impl.hpp"

// TODO: Sequential map
//#include "El/core/Map.hpp"

#include "El/core/DistMap.hpp"

#include "El/core/Permutation.hpp"
#include "El/core/DistPermutation.hpp"
