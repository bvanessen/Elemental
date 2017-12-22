#ifndef EL_TYPES_COMPLEX_DECL_HPP_
#define EL_TYPES_COMPLEX_DECL_HPP_

#include <complex>

#include "El/config.h"

#ifdef HYDROGEN_HAVE_MPC
#include <mpc.h>
#endif // HYDROGEN_HAVE_MPC

#ifdef HYDROGEN_HAVE_QUADMATH
#include <quadmath.h>
#endif

namespace El
{

template<typename Real>
class Complex
{
public:
    Real realPart, imagPart;
};

template<>
class Complex<float> : public std::complex<float>
{
public:
    typedef float realType;
    // TODO: Extend operators to other types?
    using std::complex<realType>::operator=;
    using std::complex<realType>::operator-=;
    using std::complex<realType>::operator+=;
    using std::complex<realType>::operator*=;
    using std::complex<realType>::operator/=;

    template<typename S>
    inline Complex( const S& a );
    template<typename S>
    inline Complex( const Complex<S>& a );

    template<typename S,typename T>
    inline Complex( const S& a, const T& b );

    inline Complex();
    inline Complex( const std::complex<realType>& a );
};

template<>
class Complex<double> : public std::complex<double>
{
public:
    typedef double realType;
    // TODO: Extend operators to other types?
    using std::complex<realType>::operator=;
    using std::complex<realType>::operator-=;
    using std::complex<realType>::operator+=;
    using std::complex<realType>::operator*=;
    using std::complex<realType>::operator/=;

    template<typename S>
    inline Complex( const S& a );
    template<typename S>
    inline Complex( const Complex<S>& a );

    template<typename S,typename T>
    inline Complex( const S& a, const T& b );

    inline Complex();
    inline Complex( const std::complex<realType>& a );
};

#ifdef HYDROGEN_HAVE_QUADMATH
template<>
class Complex<Quad> : public std::complex<Quad>
{
public:
    typedef Quad realType;
    // TODO: Extend operators to other types?
    using std::complex<realType>::operator=;
    using std::complex<realType>::operator-=;
    using std::complex<realType>::operator+=;
    using std::complex<realType>::operator*=;
    using std::complex<realType>::operator/=;

    template<typename S>
    inline Complex( const S& a );
    template<typename S>
    inline Complex( const Complex<S>& a );

    template<typename S,typename T>
    inline Complex( const S& a, const T& b );

    inline Complex();
    inline Complex( const std::complex<realType>& a );
};
#endif

#ifdef HYDROGEN_HAVE_QD
template<>
class Complex<DoubleDouble>
{
public:
    typedef DoubleDouble realType;
    realType realPart, imagPart;

    template<typename S>
    inline Complex( const S& a );
    template<typename S>
    inline Complex( const Complex<S>& a );

    template<typename S,typename T>
    inline Complex( const S& a, const T& b );

    inline Complex();
    inline Complex( const Complex<realType>& a );

    inline ~Complex();

    inline realType real() const;
    inline realType imag() const;

    inline void real( const realType& newReal );
    inline void imag( const realType& newImag );

    template<typename S>
    inline Complex<realType>& operator=( const S& a );
    template<typename S>
    inline Complex<realType>& operator=( const Complex<S>& a );
    inline Complex<realType>& operator=( const Complex<realType>& a );

    template<typename S>
    inline Complex<realType>& operator+=( const S& a );
    template<typename S>
    inline Complex<realType>& operator+=( const Complex<S>& a );
    inline Complex<realType>& operator+=( const Complex<realType>& a );

    template<typename S>
    inline Complex<realType>& operator-=( const S& a );
    template<typename S>
    inline Complex<realType>& operator-=( const Complex<S>& a );
    inline Complex<realType>& operator-=( const Complex<realType>& a );

    template<typename S>
    inline Complex<realType>& operator*=( const S& a );
    template<typename S>
    inline Complex<realType>& operator*=( const Complex<S>& a );
    inline Complex<realType>& operator*=( const Complex<realType>& a );

    template<typename S>
    inline Complex<realType>& operator/=( const S& a );
    template<typename S>
    inline Complex<realType>& operator/=( const Complex<S>& a );
    inline Complex<realType>& operator/=( const Complex<realType>& a );
};

template<>
class Complex<QuadDouble>
{
public:
    typedef QuadDouble realType;
    realType realPart, imagPart;

    template<typename S>
    inline Complex( const S& a );
    template<typename S>
    inline Complex( const Complex<S>& a );

    template<typename S,typename T>
    inline Complex( const S& a, const T& b );

    inline Complex();
    inline Complex( const Complex<realType>& a );

    inline ~Complex();

    inline realType real() const;
    inline realType imag() const;

    inline void real( const realType& newReal );
    inline void imag( const realType& newImag );

    template<typename S>
    inline Complex<realType>& operator=( const S& a );
    template<typename S>
    inline Complex<realType>& operator=( const Complex<S>& a );
    inline Complex<realType>& operator=( const Complex<realType>& a );

    template<typename S>
    inline Complex<realType>& operator+=( const S& a );
    template<typename S>
    inline Complex<realType>& operator+=( const Complex<S>& a );
    inline Complex<realType>& operator+=( const Complex<realType>& a );

    template<typename S>
    inline Complex<realType>& operator-=( const S& a );
    template<typename S>
    inline Complex<realType>& operator-=( const Complex<S>& a );
    inline Complex<realType>& operator-=( const Complex<realType>& a );

    template<typename S>
    inline Complex<realType>& operator*=( const S& a );
    template<typename S>
    inline Complex<realType>& operator*=( const Complex<S>& a );
    inline Complex<realType>& operator*=( const Complex<realType>& a );

    template<typename S>
    inline Complex<realType>& operator/=( const S& a );
    template<typename S>
    inline Complex<realType>& operator/=( const Complex<S>& a );
    inline Complex<realType>& operator/=( const Complex<realType>& a );
};
#endif

#ifdef HYDROGEN_HAVE_MPC
// We force the precision of the real and imaginary components to be the same...
// Note that we are requiring that the 'j' functions accept long long integers,
// which implies that intmax_t is at least as large as long long int
template<>
class Complex<BigFloat>
{
private:
    mpc_t mpcFloat_;
    size_t numLimbs_;

    inline void SetNumLimbs( mpfr_prec_t prec );
    inline void Init( mpfr_prec_t prec=mpfr::Precision() );

public:
    typedef BigFloat realType;

    inline mpc_ptr Pointer();
    inline mpc_srcptr LockedPointer() const;

    inline mpfr_ptr RealPointer();
    inline mpfr_ptr ImagPointer();
    inline mpfr_srcptr LockedRealPointer() const;
    inline mpfr_srcptr LockedImagPointer() const;

    inline mpfr_prec_t Precision() const;
    inline void SetPrecision( mpfr_prec_t );
    inline size_t NumLimbs() const;

    inline realType real() const;
    inline realType imag() const;

    inline void real( const realType& newReal );
    inline void imag( const realType& newImag );

    inline Complex();

    template<typename S>
    inline Complex( const S& a, mpfr_prec_t prec=mpfr::Precision() );
    template<typename S>
    inline Complex( const Complex<S>& a, mpfr_prec_t prec=mpfr::Precision() );

    inline Complex
    ( const unsigned& a, mpfr_prec_t prec=mpfr::Precision() );
    inline Complex
    ( const unsigned long long& a, mpfr_prec_t prec=mpfr::Precision() );
    inline Complex
    ( const int& a, mpfr_prec_t prec=mpfr::Precision() );
    inline Complex
    ( const long long int& a, mpfr_prec_t prec=mpfr::Precision() );
    inline Complex
    ( const BigInt& a, mpfr_prec_t prec=mpfr::Precision() );
    inline Complex
    ( const BigInt& a,
      const BigInt& b,
            mpfr_prec_t prec=mpfr::Precision() );
    inline Complex
    ( const float& a, mpfr_prec_t prec=mpfr::Precision() );
    inline Complex
    ( const std::complex<float>& a, mpfr_prec_t prec=mpfr::Precision() );
    inline Complex
    ( const double& a, mpfr_prec_t prec=mpfr::Precision() );
    inline Complex
    ( const std::complex<double>& a, mpfr_prec_t prec=mpfr::Precision() );
    inline Complex
    ( const realType& a, mpfr_prec_t prec=mpfr::Precision() );
    inline Complex
    ( const realType& a,
      const realType& b,
            mpfr_prec_t prec=mpfr::Precision() );
    inline Complex
    ( const Complex<realType>& a,
            mpfr_prec_t prec=mpfr::Precision() );
    inline Complex( Complex<realType>&& a );
    inline ~Complex();

    template<typename S>
    inline Complex<realType>& operator=( const Complex<S>& a );
    template<typename S>
    inline Complex<realType>& operator=( const S& a );

    inline Complex<realType>& operator=( Complex<realType>&& a );
    inline Complex<realType>& operator=( const Complex<realType>& a );
    inline Complex<realType>& operator=( const realType& a );
    inline Complex<realType>& operator=( const BigInt& a );
    inline Complex<realType>& operator=( const Complex<double>& a );
    inline Complex<realType>& operator=( const double& a );
    inline Complex<realType>& operator=( const Complex<float>& a );
    inline Complex<realType>& operator=( const float& a );
    inline Complex<realType>& operator=( const long long int& a );
    inline Complex<realType>& operator=( const long int& a );
    inline Complex<realType>& operator=( const int& a );
    inline Complex<realType>& operator=( const unsigned long long& a );
    inline Complex<realType>& operator=( const unsigned long& a );
    inline Complex<realType>& operator=( const unsigned& a );

    template<typename S>
    inline Complex<realType>& operator+=( const S& a );
    template<typename S>
    inline Complex<realType>& operator+=( const Complex<S>& a );

    inline Complex<realType>& operator+=( const Complex<realType>& a );
    inline Complex<realType>& operator+=( const realType& a );
    inline Complex<realType>& operator+=( const Complex<double>& a );
    inline Complex<realType>& operator+=( const double& a );
    inline Complex<realType>& operator+=( const Complex<float>& a );
    inline Complex<realType>& operator+=( const float& a );
    inline Complex<realType>& operator+=( const long long int& a );
    inline Complex<realType>& operator+=( const long int& a );
    inline Complex<realType>& operator+=( const int& a );
    inline Complex<realType>& operator+=( const unsigned long long& a );
    inline Complex<realType>& operator+=( const unsigned long& a );
    inline Complex<realType>& operator+=( const unsigned& a );

    template<typename S>
    inline Complex<realType>& operator-=( const S& a );
    template<typename S>
    inline Complex<realType>& operator-=( const Complex<S>& a );

    inline Complex<realType>& operator-=( const Complex<realType>& a );
    inline Complex<realType>& operator-=( const realType& a );
    inline Complex<realType>& operator-=( const Complex<double>& a );
    inline Complex<realType>& operator-=( const double& a );
    inline Complex<realType>& operator-=( const Complex<float>& a );
    inline Complex<realType>& operator-=( const float& a );
    inline Complex<realType>& operator-=( const long long int& a );
    inline Complex<realType>& operator-=( const long int& a );
    inline Complex<realType>& operator-=( const int& a );
    inline Complex<realType>& operator-=( const unsigned long long& a );
    inline Complex<realType>& operator-=( const unsigned long& a );
    inline Complex<realType>& operator-=( const unsigned& a );

    template<typename S>
    inline Complex<realType>& operator*=( const S& a );
    template<typename S>
    inline Complex<realType>& operator*=( const Complex<S>& a );

    inline Complex<realType>& operator*=( const Complex<realType>& a );
    inline Complex<realType>& operator*=( const realType& a );
    inline Complex<realType>& operator*=( const Complex<double>& a );
    inline Complex<realType>& operator*=( const double& a );
    inline Complex<realType>& operator*=( const Complex<float>& a );
    inline Complex<realType>& operator*=( const float& a );
    inline Complex<realType>& operator*=( const long long int& a );
    inline Complex<realType>& operator*=( const long int& a );
    inline Complex<realType>& operator*=( const int& a );
    inline Complex<realType>& operator*=( const unsigned long long& a );
    inline Complex<realType>& operator*=( const unsigned long& a );
    inline Complex<realType>& operator*=( const unsigned& a );

    template<typename S>
    inline Complex<realType>& operator/=( const S& a );
    template<typename S>
    inline Complex<realType>& operator/=( const Complex<S>& a );

    inline Complex<realType>& operator/=( const Complex<realType>& a );
    inline Complex<realType>& operator/=( const realType& a );
    inline Complex<realType>& operator/=( const Complex<double>& a );
    inline Complex<realType>& operator/=( const double& a );
    inline Complex<realType>& operator/=( const Complex<float>& a );
    inline Complex<realType>& operator/=( const float& a );
    inline Complex<realType>& operator/=( const long long int& a );
    inline Complex<realType>& operator/=( const long int& a );
    inline Complex<realType>& operator/=( const int& a );
    inline Complex<realType>& operator/=( const unsigned long long& a );
    inline Complex<realType>& operator/=( const unsigned long& a );
    inline Complex<realType>& operator/=( const unsigned& a );

    inline void Zero();

    inline size_t SerializedSize() const;

    inline       byte* Serialize( byte* buf ) const;
    inline const byte* Deserialize( const byte* buf );
    inline       byte* Deserialize( byte* buf );
};
#endif // HYDROGEN_HAVE_MPC

#ifdef HYDROGEN_HAVE_QD
inline bool operator==
( const Complex<DoubleDouble>& a, const Complex<DoubleDouble>& b );
inline bool operator!=
( const Complex<DoubleDouble>& a, const Complex<DoubleDouble>& b );
inline bool operator==
( const Complex<QuadDouble>& a, const Complex<QuadDouble>& b );
inline bool operator!=
( const Complex<QuadDouble>& a, const Complex<QuadDouble>& b );

inline bool operator==
( const Complex<DoubleDouble>& a, const DoubleDouble& b );
inline bool operator!=
( const Complex<DoubleDouble>& a, const DoubleDouble& b );
inline bool operator==
( const Complex<QuadDouble>& a, const QuadDouble& b );
inline bool operator!=
( const Complex<QuadDouble>& a, const QuadDouble& b );

inline bool operator==
( const DoubleDouble& a, const Complex<DoubleDouble>& b );
inline bool operator!=
( const DoubleDouble& a, const Complex<DoubleDouble>& b );
inline bool operator==
( const QuadDouble& a, const Complex<QuadDouble>& b );
inline bool operator!=
( const QuadDouble& a, const Complex<QuadDouble>& b );
#endif
#ifdef HYDROGEN_HAVE_MPC
inline bool operator==
( const Complex<BigFloat>& a, const Complex<BigFloat>& b );
inline bool operator!=
( const Complex<BigFloat>& a, const Complex<BigFloat>& b );

inline bool operator==
( const Complex<BigFloat>& a, const BigFloat& b );
inline bool operator!=
( const Complex<BigFloat>& a, const BigFloat& b );

inline bool operator==
( const BigFloat& a, const Complex<BigFloat>& b );
inline bool operator!=
( const BigFloat& a, const Complex<BigFloat>& b );
#endif

template<typename Real>
inline Complex<Real> operator-( const Complex<Real>& a );
#ifdef HYDROGEN_HAVE_QD
inline Complex<DoubleDouble> operator-( const Complex<DoubleDouble>& a );
inline Complex<QuadDouble> operator-( const Complex<QuadDouble>& a );
#endif
#ifdef HYDROGEN_HAVE_MPC
inline Complex<BigFloat> operator-( const Complex<BigFloat>& a );
#endif

template<typename Real>
inline Complex<Real>
operator+( const Complex<Real>& a, const Complex<Real>& b );
template<typename Real>
inline Complex<Real>
operator-( const Complex<Real>& a, const Complex<Real>& b );
template<typename Real>
inline Complex<Real>
operator*( const Complex<Real>& a, const Complex<Real>& b );
template<typename Real>
inline Complex<Real>
operator/( const Complex<Real>& a, const Complex<Real>& b );
#ifdef HYDROGEN_HAVE_QD
inline Complex<DoubleDouble>
operator+( const Complex<DoubleDouble>& a, const Complex<DoubleDouble>& b );
inline Complex<DoubleDouble>
operator-( const Complex<DoubleDouble>& a, const Complex<DoubleDouble>& b );
inline Complex<DoubleDouble>
operator*( const Complex<DoubleDouble>& a, const Complex<DoubleDouble>& b );
inline Complex<DoubleDouble>
operator/( const Complex<DoubleDouble>& a, const Complex<DoubleDouble>& b );
inline Complex<QuadDouble>
operator+( const Complex<QuadDouble>& a, const Complex<QuadDouble>& b );
inline Complex<QuadDouble>
operator-( const Complex<QuadDouble>& a, const Complex<QuadDouble>& b );
inline Complex<QuadDouble>
operator*( const Complex<QuadDouble>& a, const Complex<QuadDouble>& b );
inline Complex<QuadDouble>
operator/( const Complex<QuadDouble>& a, const Complex<QuadDouble>& b );
#endif
#ifdef HYDROGEN_HAVE_MPC
inline Complex<BigFloat>
operator+( const Complex<BigFloat>& a, const Complex<BigFloat>& b );
inline Complex<BigFloat>
operator-( const Complex<BigFloat>& a, const Complex<BigFloat>& b );
inline Complex<BigFloat>
operator*( const Complex<BigFloat>& a, const Complex<BigFloat>& b );
inline Complex<BigFloat>
operator/( const Complex<BigFloat>& a, const Complex<BigFloat>& b );
#endif

template<typename Real>
inline Complex<Real>
operator+( const Complex<Real>& a, const Real& b );
template<typename Real>
inline Complex<Real>
operator-( const Complex<Real>& a, const Real& b );
template<typename Real>
inline Complex<Real>
operator*( const Complex<Real>& a, const Real& b );
template<typename Real>
inline Complex<Real>
operator/( const Complex<Real>& a, const Real& b );
#ifdef HYDROGEN_HAVE_QD
inline Complex<DoubleDouble>
operator+( const Complex<DoubleDouble>& a, const DoubleDouble& b );
inline Complex<DoubleDouble>
operator-( const Complex<DoubleDouble>& a, const DoubleDouble& b );
inline Complex<DoubleDouble>
operator*( const Complex<DoubleDouble>& a, const DoubleDouble& b );
inline Complex<DoubleDouble>
operator/( const Complex<DoubleDouble>& a, const DoubleDouble& b );
inline Complex<QuadDouble>
operator+( const Complex<QuadDouble>& a, const QuadDouble& b );
inline Complex<QuadDouble>
operator-( const Complex<QuadDouble>& a, const QuadDouble& b );
inline Complex<QuadDouble>
operator*( const Complex<QuadDouble>& a, const QuadDouble& b );
inline Complex<QuadDouble>
operator/( const Complex<QuadDouble>& a, const QuadDouble& b );
#endif
#ifdef HYDROGEN_HAVE_MPC
inline Complex<BigFloat>
operator+( const Complex<BigFloat>& a, const BigFloat& b );
inline Complex<BigFloat>
operator-( const Complex<BigFloat>& a, const BigFloat& b );
inline Complex<BigFloat>
operator*( const Complex<BigFloat>& a, const BigFloat& b );
inline Complex<BigFloat>
operator/( const Complex<BigFloat>& a, const BigFloat& b );
#endif

template<typename Real>
inline Complex<Real>
operator+( const Real& a, const Complex<Real>& b );
template<typename Real>
inline Complex<Real>
operator-( const Real& a, const Complex<Real>& b );
template<typename Real>
inline Complex<Real>
operator*( const Real& a, const Complex<Real>& b );
template<typename Real>
inline Complex<Real>
operator/( const Real& a, const Complex<Real>& b );
#ifdef HYDROGEN_HAVE_QD
inline Complex<DoubleDouble>
operator+( const DoubleDouble& a, const Complex<DoubleDouble>& b );
inline Complex<DoubleDouble>
operator-( const DoubleDouble& a, const Complex<DoubleDouble>& b );
inline Complex<DoubleDouble>
operator*( const DoubleDouble& a, const Complex<DoubleDouble>& b );
inline Complex<DoubleDouble>
operator/( const DoubleDouble& a, const Complex<DoubleDouble>& b );
inline Complex<QuadDouble>
operator+( const QuadDouble& a, const Complex<QuadDouble>& b );
inline Complex<QuadDouble>
operator-( const QuadDouble& a, const Complex<QuadDouble>& b );
inline Complex<QuadDouble>
operator*( const QuadDouble& a, const Complex<QuadDouble>& b );
inline Complex<QuadDouble>
operator/( const QuadDouble& a, const Complex<QuadDouble>& b );
#endif
#ifdef HYDROGEN_HAVE_MPC
inline Complex<BigFloat>
operator+( const BigFloat& a, const Complex<BigFloat>& b );
inline Complex<BigFloat>
operator-( const BigFloat& a, const Complex<BigFloat>& b );
inline Complex<BigFloat>
operator*( const BigFloat& a, const Complex<BigFloat>& b );
inline Complex<BigFloat>
operator/( const BigFloat& a, const Complex<BigFloat>& b );
#endif

using scomplex = Complex<float>;
using dcomplex = Complex<double>;
#ifdef HYDROGEN_HAVE_QUADMATH
using qcomplex = Complex<Quad>;
#endif
#ifdef HYDROGEN_HAVE_QD
using ddcomplex = Complex<DoubleDouble>;
using qdcomplex = Complex<QuadDouble>;
#endif
#ifdef HYDROGEN_HAVE_MPC
using acomplex = Complex<BigFloat>;
#endif
}// namespace El

#endif // EL_TYPES_COMPLEX_DECL_HPP_
