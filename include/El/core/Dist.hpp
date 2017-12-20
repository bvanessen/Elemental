#ifndef EL_CORE_DIST_HPP_
#define EL_CORE_DIST_HPP_

#include "El/macros.h"

namespace El
{

/** \brief Types of distributions used to describe matrix/vector
 *      distributions.
 */
enum class Dist
{
    MC,   // Col of a matrix distribution
    MD,   // Diagonal of a matrix distribution
    MR,   // Row of a matrix distribution
    VC,   // Col-major vector distribution
    VR,   // Row-major vector distribution
    STAR, // Give to every process
    CIRC  // Give to a single process
};

std::string DistToString(Dist distribution);

Dist StringToDist(std::string s);

// Handy typedef? IDEK...
using Distribution = Dist;

template<Dist U,Dist V> constexpr Dist DiagCol() { return (U==Dist::STAR ? V : U); }
template<Dist U,Dist V> constexpr Dist DiagRow() { return (U==Dist::STAR ? U : V); }
template<> constexpr Dist DiagCol<Dist::MC,Dist::MR>() { return Dist::MD; }
template<> constexpr Dist DiagRow<Dist::MC,Dist::MR>() { return Dist::STAR; }
template<> constexpr Dist DiagCol<Dist::MR,Dist::MC>() { return Dist::MD; }
template<> constexpr Dist DiagRow<Dist::MR,Dist::MC>() { return Dist::STAR; }

// Runtime
// -------
inline Dist DiagCol(Dist U, Dist V) EL_NO_EXCEPT
{
    if (U == Dist::MC && V == Dist::MR)
        return Dist::MD;
    else if (U == Dist::MR && V == Dist::MC)
        return Dist::MD;
    else if (U == Dist::STAR)
        return V;
    else
        return U;
}
inline Dist DiagRow(Dist U, Dist V) EL_NO_EXCEPT
{
    if (U == Dist::MC && V == Dist::MR)
        return Dist::STAR;
    else if (U == Dist::MR && V == Dist::MC)
        return Dist::STAR;
    else if (U == Dist::STAR)
        return U;
    else
        return V;
}

// Return a piece of a distribution which induces the given diagonal dist
// ======================================================================

// Compile-time
// ------------
template<Dist U,Dist V>
constexpr Dist DiagInvCol() { return (U==Dist::STAR ? V : U); }
template<Dist U,Dist V>
constexpr Dist DiagInvRow() { return (U==Dist::STAR ? U : V); }
template<> constexpr Dist DiagInvCol<Dist::MD,Dist::STAR>() { return Dist::MC; }
template<> constexpr Dist DiagInvRow<Dist::MD,Dist::STAR>() { return Dist::MR; }
template<> constexpr Dist DiagInvCol<Dist::STAR,Dist::MD>() { return Dist::MC; }
template<> constexpr Dist DiagInvRow<Dist::STAR,Dist::MD>() { return Dist::MR; }

// TODO: Runtime version?

// Union the distribution over its corresponding communicator
// ==========================================================
// Compile-time
// ------------
template<Dist U> constexpr Dist Collect() { return Dist::STAR; }
template<> constexpr Dist Collect<Dist::CIRC>() { return Dist::CIRC; }
// Run-time
// --------
inline Dist Collect(Dist U) EL_NO_EXCEPT { return (U==Dist::CIRC ? Dist::CIRC : Dist::STAR); }

// Union the distribution over its corresponding partial communicator
// ==================================================================
// Compile-time
// ------------
template<Dist U> constexpr Dist Partial() { return U; }
template<> constexpr Dist Partial<Dist::VC>() { return Dist::MC; }
template<> constexpr Dist Partial<Dist::VR>() { return Dist::MR; }
// Run-time
// --------
inline Dist Partial(Dist U) EL_NO_EXCEPT
{
    if (U == Dist::VC)
        return Dist::MC;
    else if (U == Dist::VR)
        return Dist::MR;
    else
        return U;
}

// Return the partial distribution that would be used for a partial union
// ======================================================================
// Compile-time
// ------------
template<Dist U,Dist V> constexpr Dist PartialUnionRow() { return V; }
template<> constexpr Dist PartialUnionRow<Dist::VC,Dist::STAR>() { return Dist::MR; }
template<> constexpr Dist PartialUnionRow<Dist::VR,Dist::STAR>() { return Dist::MC; }template<Dist U,Dist V> constexpr Dist PartialUnionCol()
{ return PartialUnionRow<V,U>(); }
// Run-time
// --------
inline Dist PartialUnionRow(Dist U, Dist V) EL_NO_EXCEPT
{
    if (U == Dist::VC)
        return Dist::MR;
    else if (U == Dist::VR)
        return Dist::MC;
    else
        return V;
}
inline Dist PartialUnionCol(Dist U, Dist V) EL_NO_EXCEPT
{ return PartialUnionRow(V, U); }

// Return the product of two distributions
// =======================================
// Compile-time
// ------------
template<Dist U,Dist V> constexpr Dist ProductDist() { return Dist::CIRC; }
template<> constexpr Dist ProductDist<Dist::MC, Dist::MR>() { return Dist::VC; }
template<> constexpr Dist ProductDist<Dist::MC, Dist::STAR>() { return Dist::MC; }
template<> constexpr Dist ProductDist<Dist::MD, Dist::STAR>() { return Dist::MD; }
template<> constexpr Dist ProductDist<Dist::MR, Dist::MC>() { return Dist::VR; }
template<> constexpr Dist ProductDist<Dist::MR, Dist::STAR>() { return Dist::MR; }
template<> constexpr Dist ProductDist<Dist::STAR,Dist::MC>() { return Dist::MC; }
template<> constexpr Dist ProductDist<Dist::STAR,Dist::MD>() { return Dist::MD; }
template<> constexpr Dist ProductDist<Dist::STAR,Dist::MR>() { return Dist::MR; }
template<> constexpr Dist ProductDist<Dist::STAR,Dist::STAR>() { return Dist::STAR; }
template<> constexpr Dist ProductDist<Dist::STAR,Dist::VC>() { return Dist::VC; }
template<> constexpr Dist ProductDist<Dist::STAR,Dist::VR>() { return Dist::VR; }
template<> constexpr Dist ProductDist<Dist::VC, Dist::STAR>() { return Dist::VC; }
template<> constexpr Dist ProductDist<Dist::VR, Dist::STAR>() { return Dist::VR; }
template<Dist U,Dist V>
constexpr Dist ProductDistPartner() { return Dist::STAR; }
template<>
constexpr Dist ProductDistPartner<Dist::CIRC,Dist::CIRC>() { return Dist::CIRC; }
// Runtime
// -------
inline Dist ProductDist(Dist U, Dist V) EL_NO_EXCEPT
{
    if (U == Dist::STAR) return V;
    else if (V == Dist::STAR) return U;
    else if (U == Dist::MC && V == Dist::MR) return Dist::VC;
    else if (U == Dist::MR && V == Dist::MC) return Dist::VR;
    else if (U == Dist::CIRC && V == Dist::CIRC) return Dist::CIRC;
    else { return Dist::STAR; } // NOTE: This branch should not be possible
}
inline Dist ProductDistPartner(Dist U, Dist V) EL_NO_EXCEPT
{
    if (U == Dist::CIRC && V == Dist::CIRC)
        return Dist::CIRC;
    else
        return Dist::STAR;
}

}// namespace El

#endif // EL_CORE_DIST_HPP_
