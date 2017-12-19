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

}// namespace El
#endif // EL_TYPES_ENUMS_HPP_
