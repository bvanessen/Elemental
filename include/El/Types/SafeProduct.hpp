#ifndef EL_TYPES_SAFEPRODUCT_HPP_
#define EL_TYPES_SAFEPRODUCT_HPP_

namespace El
{

// For the safe computation of products. The result is given by
//   product = rho * exp(kappa*n)
// where rho lies in (usually on) the unit circle and kappa is real-valued.
template<typename F>
struct SafeProduct
{
    F rho;
    Base<F> kappa;
    Int n;

    SafeProduct( Int numEntries );
};
}// namespace El
#endif // EL_TYPES_SAFEPRODUCT_HPP_
