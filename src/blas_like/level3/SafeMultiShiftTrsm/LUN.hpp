/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace safemstrsm {

/* Determine machine dependent parameters to control overflow
 *   Note: LAPACK uses more complicated parameters to handle 
 *   issues that can happen on Cray machines.
 */
template<typename Real>
inline void
OverflowParameters( Real& smlnum, Real& bignum )
{
    const Real unfl = lapack::MachineSafeMin<Real>();
    const Real ovfl = lapack::MachineOverflowThreshold<Real>();
    const Real ulp  = lapack::MachinePrecision<Real>();
    smlnum = std::max( unfl/ulp, 1/(ovfl*ulp) );
    bignum = 1/smlnum;
}
  
template<typename F>
void
LUNBlock
(       Matrix<F>& U,
  const Matrix<F>& shifts,
        Matrix<F>& X,
	Matrix<F>& scales )
{

    typedef Base<F> Real;
  
    DEBUG_ONLY(
      CSE cse("safemstrsm::LUNBlock");
      if( shifts.Height() != X.Width() )
          LogicError("Incompatible number of shifts");
    )
    auto diag = GetDiagonal(U);
    const Int n = U.Height();
    const Int ldim = U.LDim();
    const Int numShifts = shifts.Height();

    Real smlnum, bignum;
    OverflowParameters<Real>( smlnum, bignum );
    
    // Default scale is 1
    Ones( scales, numShifts, 1 );

    // Compute infinity norms of columns of U (excluding diagonal)
    // TODO: scale cnorm if an entry is bigger than bignum
    Matrix<Real> cnorm( n, 1 );
    cnorm.Set( 0, 0, Real(0) );
    for( Int j=1; j<n; ++j )
    {
        cnorm.Set( j, 0, MaxNorm( U(IR(0,j),IR(j)) ) );
    }

    // Iterate through RHS's
    for( Int j=0; j<numShifts; ++j )
    {

        // Initialize triangular system
        ShiftDiagonal( U, -shifts.Get(j,0) );
	auto xj = X( ALL, IR(j) );

	// Determine largest entry of RHS
	Real xjMax = MaxNorm( xj );
	if( xjMax >= bignum )
	{
	    const Real s = 0.5*bignum/xjMax;
	    xj *= s;
	    xjMax *= s;
	    scales.Set( j, 0, s*scales.Get(j,0) );
	}
	xjMax = std::max( xjMax, 2*smlnum );

	// Estimate growth of entries in triangular solve
	//   Note: See "Robust Triangular Solves for Use in Condition
	//   Estimation" by Edward Anderson for explanation of bounds.
	Real invGi = 1/xjMax;
	Real invMi = invGi;
	for( Int i=n-1; i>=0; --i )
	{
	    const Real absUii = SafeAbs( U.Get(i,i) );
	    if( invGi<=smlnum || invMi<=smlnum || absUii<=smlnum )
	    {
	        invGi = 0;
	        break;
	    }
	    invMi = std::min( invMi, absUii*invGi );
	    invGi *= absUii/(absUii+cnorm.Get(i,0));
	}
	invGi = std::min( invGi, invMi );

	// Call TRSV if estimated growth is not too large
	if( invGi > smlnum )
	{
	  blas::Trsv
	  ( 'U', 'N', 'N', n,
	    U.LockedBuffer(), ldim, X.Buffer(0,j), 1 );
	}

	// Perform backward substitution if estimated growth is large
	else
	{
	    for( Int i=n-1; i>=0; --i )
	    {

		// Perform division and check for overflow
		const Real absUii = SafeAbs( U.Get(i,i) );
		Real absXij = SafeAbs( xj.Get(i,0) );
		if( absUii > smlnum )
		{
		    if( absUii<=1 && absXij>=absUii*bignum )
		    {
			// Set overflowing entry to 0.5/U[i,i]
		        const Real s = 0.5/absXij;
			scales.Set( j, 0, s*scales.Get(j,0) );
			xj *= s;
			xjMax *= s;
		    }
		    xj.Set( i, 0, xj.Get(i,0)/U.Get(i,i) );
		}
		else if( absUii > 0 )
		{
		    if( absXij >= absUii*bignum )
		    {
			// Set overflowing entry to bignum/2
		        const Real s = 0.5*absUii*bignum/absXij;
			scales.Set( j, 0, s*scales.Get(j,0) );
			xj *= s;
			xjMax *= s;
		    }
		    xj.Set( i, 0, xj.Get(i,0)/U.Get(i,i) );
		}
		else
		{
		    // TODO: maybe this tolerance should be loosened to
		    //   | Xij | >= || A || * eps
		    if( absXij >= smlnum )
		    {
		        scales.Set( j, 0, F(0) );
			Zero( xj );
			xjMax = 0;
		    }
		    xj.Set( i, 0, F(1) );
		}

		if( i > 0 )
	        {

		    // Check for possible overflows in AXPY
		    // Note: G(i+1) <= G(i) + | Xij | * cnorm(i)
		    absXij = SafeAbs( xj.Get(i,0) );
		    if( absXij >= 1 &&
			cnorm.Get(i,0) >= (bignum-xjMax)/absXij )
		    {
		        const Real s = 0.25/absXij;
			scales.Set( j, 0, s*scales.Get(j,0) );
			xj *= s;
			xjMax *= s;
			absXij *= s;
		    }
		    else if( absXij < 1 &&
			     absXij*cnorm.Get(i,0) >= bignum-xjMax )
		    {
		        const Real s = 0.25;
		        scales.Set( j, 0, s*scales.Get(j,0) );
			xj *= s;
			xjMax *= s;
			absXij *= s;
		    }
		    xjMax += absXij * cnorm.Get(i,0);

		    // AXPY
		    auto U01 = U( IR(0,i), IR(i) );
		    auto X1  = X( IR(0,i), IR(j) );
		    Axpy( -xj.Get(i,0), U01, X1 );

		}

	    }

	}

	// Reset matrix diagonal
        SetDiagonal( U, diag );
    }
}

template<typename F>
void
LUN( Matrix<F>& U, const Matrix<F>& shifts,
     Matrix<F>& X, Matrix<F>& scales ) 
{
    typedef Base<F> Real;

    DEBUG_ONLY(CSE cse("safemstrsm::LUN"))
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( m, bsize );

    Real smlnum, bignum;
    OverflowParameters<Real>( smlnum, bignum );

    Ones( scales, n, 1 );
    Matrix<F> scalesUpdate( n, 1 );

    // Determine largest entry of each RHS
    Matrix<Real> XMax( n, 1 );
    for( Int j=0; j<n; ++j )
    {
        auto xj = X( ALL, IR(j) );
	Real xjMax = MaxNorm( xj );
	if( xjMax >= bignum )
	{
	    const Real s = 0.5*bignum/xjMax;
	    xj *= s;
	    xjMax *= s;
	    scales.Set( j, 0, s*scales.Get(j,0) );
	}
	xjMax = std::max( xjMax, 2*smlnum );
        XMax.Set( j, 0, xjMax );
    }
	
    // Perform block triangular solve
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        const Range<Int> ind0( 0,    k    ),
	                 ind1( k,    k+nb ),
	                 ind2( k+nb, END  );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );
	auto X2 = X( ind2, ALL );

	// Perform triangular solve on diagonal block
	LUNBlock( U11, shifts, X1, scalesUpdate );

	// Apply scalings on RHS
	for( Int j=0; j<n; ++j )
	{
	    const Real sj = scalesUpdate.GetRealPart(j,0);
	    if( sj < 1 )
	    {
	        scales.Set( j, 0, sj*scales.Get(j,0) );
	        auto X0j = X0( ALL, IR(j) );
	        auto X2j = X2( ALL, IR(j) );
		X0j *= sj;
		X2j *= sj;
		XMax.Set( j, 0, sj*XMax.Get(j,0) );
	    }
	}

	if( k > 0 )
	{

	    // Compute infinity norms of columns in U01
	    // Note: nb*cnorm is the sum of infinity norms
	    // TODO: scale cnorm if an entry is bigger than bignum
	    Real cnorm = 0;
	    for( Int j=0; j<nb; ++j )
	    {
	        cnorm += MaxNorm( U01(ALL,IR(j)) ) / nb;
	    }

	    // Check for possible overflows in GEMM
	    // Note: G(i+1) <= G(i) + nb*cnorm*|| X1[:,j] ||_infty
	    for( Int j=0; j<n; ++j )
	    {
		auto xj = X( ALL, IR(j) );
	        Real xjMax = XMax.Get(j,0);
		Real X1Max = MaxNorm( X1(ALL,IR(j)) );
		if( X1Max >= 1 &&
		    cnorm >= (bignum-xjMax)/X1Max/nb )
		{
		    const Real s = 0.5/X1Max/nb;
		    scales.Set( j, 0, s*scales.Get(j,0) );
		    xj *= s;
		    xjMax *= s;
		    X1Max *= s;
		}
		else if( X1Max < 1 &&
			 cnorm*X1Max >= (bignum-xjMax)/nb )
		{
		    const Real s = 0.5/nb;
		    scales.Set( j, 0, s*scales.Get(j,0) );
		    xj *= s;
		    xjMax *= s;
		    X1Max *= s;
		}
		xjMax += nb*cnorm*X1Max;
		XMax.Set( j, 0, xjMax );
	    }

	    // Update RHS with GEMM
	    Gemm( NORMAL, NORMAL, F(-1), U01, X1, F(1), X0 );
	
	}
    }
}

#if 0
  // TODO: template specialization for elemental matrices
template<typename F>
inline void
LUN
( const ElementalMatrix<F>& UPre, 
  const ElementalMatrix<F>& shiftsPre,
        ElementalMatrix<F>& XPre ) 
{
    DEBUG_ONLY(CSE cse("mstrsm::LUN"))

    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixReadProxy<F,F,VR,STAR> shiftsProx( shiftsPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& shifts = shiftsProx.GetLocked();
    auto& X = XProx.Get();

    const Grid& g = U.Grid();
    DistMatrix<F,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    const Int m = X.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( m, bsize );

    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );

        // X1[* ,VR] := U11^-1[* ,* ] X1[* ,VR]
        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1_STAR_VR.AlignWith( shifts );
        X1_STAR_VR = X1; // X1[* ,VR] <- X1[MC,MR]
        LUN
        ( U11_STAR_STAR.Matrix(), shifts.LockedMatrix(), 
          X1_STAR_VR.Matrix() );

        X1_STAR_MR.AlignWith( X0 );
        X1_STAR_MR = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1 = X1_STAR_MR; // X1[MC,MR] <- X1[* ,MR]

        // X0[MC,MR] -= U01[MC,* ] X1[* ,MR]
        U01_MC_STAR.AlignWith( X0 );
        U01_MC_STAR = U01; // U01[MC,* ] <- U01[MC,MR]
        LocalGemm( NORMAL, NORMAL, F(-1), U01_MC_STAR, X1_STAR_MR, F(1), X0 );
    }
}
#endif
  
} // namespace safemstrsm
} // namespace El
