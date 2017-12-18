#ifndef EL_TYPEDEFS_H_
#define EL_TYPEDEFS_H_

#include "El/config.h"

namespace El
{

using byte = unsigned char;

// If these are changes, you must make sure that they have
// existing MPI datatypes. This is only sometimes true for 'long long'
#ifdef EL_USE_64BIT_INTS
typedef long long int Int;
typedef long long unsigned Unsigned;
#else
typedef int Int;
typedef unsigned Unsigned;
#endif

#ifdef HYDROGEN_HAVE_QUADMATH
typedef __float128 Quad;
#endif

// Forward declarations
// --------------------
#ifdef HYDROGEN_HAVE_QD
struct DoubleDouble;
struct QuadDouble;
#endif
#ifdef HYDROGEN_HAVE_MPC
class BigInt;
class BigFloat;
#endif
template<typename Real>
class Complex;

// Convert CMake configuration into a typedef and an enum
typedef EL_FORT_LOGICAL FortranLogical;
enum FortranLogicalEnum
{
  FORTRAN_TRUE=EL_FORT_TRUE,
  FORTRAN_FALSE=EL_FORT_FALSE
};// enum FortranLogicalEnum

}// namespace El
#endif // EL_TYPEDEFS_H_
