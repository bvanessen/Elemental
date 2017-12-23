#include <string>

#include "El/core/Dist.hpp"
#include "El/core/environment/decl.hpp"
#include "El/Types/Complex.hpp"
#include "El/Types/Enums.hpp"
#include "El/Types/SafeProduct.hpp"

namespace El
{

template<typename F>
SafeProduct<F>::SafeProduct( Int numEntries )
: rho(1), kappa(0), n(numEntries)
{ }

std::string DistToString( Dist dist )
{
    std::string distString;
    switch( dist )
    {
        case Dist::CIRC: distString = "o "; break;
        case Dist::MC:   distString = "MC";   break;
        case Dist::MD:   distString = "MD";   break;
        case Dist::MR:   distString = "MR";   break;
        case Dist::VC:   distString = "VC";   break;
        case Dist::VR:   distString = "VR";   break;
        default:   distString = "* ";   break;
    }
    return distString;
}

Dist StringToDist( std::string s )
{
    // Most compilers' logic for detecting potentially uninitialized variables
    // is horrendously bad.
    Dist dist=Dist::MC;
    if( s == "MC" )
        dist = Dist::MC;
    else if( s == "MD" )
        dist = Dist::MD;
    else if( s == "MR" )
        dist = Dist::MR;
    else if( s == "VC" )
        dist = Dist::VC;
    else if( s == "VR" )
        dist = Dist::VR;
    else if( s == "* " || s == " *" || s == "*" )
        dist = Dist::STAR;
    else if( s == "o " || s == " o" || s == "o" )
        dist = Dist::CIRC;
    else
        LogicError
        ("StringToDist expects string in "
         "{\"MC\",\"MD\",\"MR\",\"VC\",\"VR\",\"* \",\" *\",\"*\"}");
    return dist;
}

char LeftOrRightToChar( LeftOrRight side )
{
    char sideChar;
    switch( side )
    {
        case LeftOrRight::LEFT:  sideChar = 'L'; break;
        default:    sideChar = 'R'; break;
    }
    return sideChar;
}

LeftOrRight CharToLeftOrRight( char c )
{
    // Most compilers' logic for detecting potentially uninitialized variables
    // is horrendously bad.
    LeftOrRight side=LeftOrRight::LEFT;
    switch( c )
    {
        case 'L': side = LeftOrRight::LEFT;  break;
        case 'R': side = LeftOrRight::RIGHT; break;
        default:
            LogicError("CharToLeftOrRight expects char in {L,R}");
    }
    return side;
}

char OrientationToChar( Orientation orientation )
{
    char orientationChar;
    switch( orientation )
    {
        case Orientation::NORMAL:    orientationChar = 'N'; break;
        case Orientation::TRANSPOSE: orientationChar = 'T'; break;
        default:        orientationChar = 'C'; break;
    }
    return orientationChar;
}

Orientation CharToOrientation( char c )
{
    // Most compilers' logic for detecting potentially uninitialized variables
    // is horrendously bad.
    Orientation orientation=Orientation::NORMAL;
    switch( c )
    {
        case 'N': orientation = Orientation::NORMAL;    break;
        case 'T': orientation = Orientation::TRANSPOSE; break;
        case 'C': orientation = Orientation::ADJOINT;   break;
        default:
            LogicError
            ("CharToOrientation expects char in {N,T,C}");
    }
    return orientation;
}

char UnitOrNonUnitToChar( UnitOrNonUnit diag )
{
    char diagChar;
    switch( diag )
    {
        case UnitOrNonUnit::NON_UNIT: diagChar = 'N'; break;
        default:       diagChar = 'U'; break;
    }
    return diagChar;
}

UnitOrNonUnit CharToUnitOrNonUnit( char c )
{
    // Most compilers' logic for detecting potentially uninitialized variables
    // is horrendously bad.
    UnitOrNonUnit diag=UnitOrNonUnit::NON_UNIT;
    switch( c )
    {
        case 'N': diag = UnitOrNonUnit::NON_UNIT; break;
        case 'U': diag = UnitOrNonUnit::UNIT;     break;
        default:
            LogicError("CharToUnitOrNonUnit expects char in {N,U}");
    }
    return diag;
}

char UpperOrLowerToChar( UpperOrLower uplo )
{
    char uploChar;
    switch( uplo )
    {
        case UpperOrLower::LOWER: uploChar = 'L'; break;
        default:    uploChar = 'U'; break;
    }
    return uploChar;
}

UpperOrLower CharToUpperOrLower( char c )
{
    // Most compilers' logic for detecting potentially uninitialized variables
    // is horrendously bad.
    UpperOrLower uplo=UpperOrLower::LOWER;
    switch( c )
    {
        case 'L': uplo = UpperOrLower::LOWER; break;
        case 'U': uplo = UpperOrLower::UPPER; break;
        default:
            LogicError("CharToUpperOrLower expects char in {L,U}");
    }
    return uplo;
}

template struct SafeProduct<float>;
template struct SafeProduct<double>;
template struct SafeProduct<Complex<float>>;
template struct SafeProduct<Complex<double>>;
#ifdef HYDROGEN_HAVE_QD
template struct SafeProduct<DoubleDouble>;
template struct SafeProduct<QuadDouble>;
template struct SafeProduct<Complex<DoubleDouble>>;
template struct SafeProduct<Complex<QuadDouble>>;
#endif
#ifdef HYDROGEN_HAVE_QUADMATH
template struct SafeProduct<Quad>;
template struct SafeProduct<Complex<Quad>>;
#endif
#ifdef HYDROGEN_HAVE_MPC
template struct SafeProduct<BigFloat>;
template struct SafeProduct<Complex<BigFloat>>;
#endif

}// namespace El
