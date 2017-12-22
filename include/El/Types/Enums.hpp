#ifndef EL_TYPES_ENUMS_HPP_
#define EL_TYPES_ENUMS_HPP_

// FIXME (trb 12/17/17): Not sure where these belong, so they live
// here now.
//
// They probably belong in a refactored BLAS header...
namespace El
{

enum class UpperOrLower
{
    LOWER,
    UPPER
};
char UpperOrLowerToChar(UpperOrLower uplo);
UpperOrLower CharToUpperOrLower(char c);

enum class Orientation
{
    NORMAL,
    TRANSPOSE,
    ADJOINT
};
char OrientationToChar(Orientation orientation);
Orientation CharToOrientation(char c);

enum class LeftOrRight
{
    LEFT,
    RIGHT
};
char LeftOrRightToChar(LeftOrRight side);
LeftOrRight CharToLeftOrRight(char c);

enum class UnitOrNonUnit
{
    NON_UNIT,
    UNIT
};
char UnitOrNonUnitToChar( UnitOrNonUnit diag );
UnitOrNonUnit CharToUnitOrNonUnit( char c );

enum class ForwardOrBackward
{
    FORWARD,
    BACKWARD
};

enum class VerticalOrHorizontal
{
    VERTICAL,
    HORIZONTAL
};

enum class SortType
{
    UNSORTED,
    DESCENDING,
    ASCENDING
};

enum class Conjugation
{
    UNCONJUGATED,
    CONJUGATED
};

enum class NormType
{
    ONE_NORM,           // Operator one norm
    INFINITY_NORM,      // Operator infinity norm
    ENTRYWISE_ONE_NORM, // One-norm of vectorized matrix
    MAX_NORM,           // Maximum entry-wise magnitude
    NUCLEAR_NORM,       // One-norm of the singular values
    FROBENIUS_NORM,     // Two-norm of the singular values
    TWO_NORM            // Infinity-norm of the singular values
};

enum class FileFormat
{
    AUTO, // Automatically detect from file extension
    ASCII,
    ASCII_MATLAB,
    BINARY,
    BINARY_FLAT,
    BMP,
    JPG,
    JPEG,
    MATRIX_MARKET,
    PNG,
    PPM,
    XBM,
    XPM,
    FileFormat_MAX // For detecting number of entries in enum
};

}// namespace El

#endif // EL_TYPES_ENUMS_HPP_
