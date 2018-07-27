//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hlle.cpp
//  \brief HLLE Riemann solver for collisionless stellar dynamics 
//
//  Computes 1D fluxes using the Harten-Lax-van Leer (HLL) Riemann solver.  This flux is
//  very diffusive, especially for contacts, and so it is not recommended for use in
//  applications.  However, as shown by Einfeldt et al.(1991), it is positively
//  conservative (cannot return negative densities or pressure), so it is a useful
//  option when other approximate solvers fail and/or when extra dissipation is needed.
//
// REFERENCES:
// - E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics", 2nd ed.,
//   Springer-Verlag, Berlin, (1999) chpt. 10.
// - Einfeldt et al., "On Godunov-type methods near low densities", JCP, 92, 273 (1991)
// - A. Harten, P. D. Lax and B. van Leer, "On upstream differencing and Godunov-type
//   schemes for hyperbolic conservation laws", SIAM Review 25, 35-61 (1983).
// - N.L. Mitchell, E.I. Vorobyov, and G. Hensler, "Collisionless stellar hydrodynamics as
//	 an efficient alternative to N-body methods", MNRAS 428, 2674-2687 (2013) 

// C/C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena++ headers
#include "../cless.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../eos/eos.hpp"

//----------------------------------------------------------------------------------------
//! \func

void Cless::RiemannSolverCL(const int kl, const int ku, const int jl, const int ju,
  const int il, const int iu, const int ivx, const int ip12,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr, AthenaArray<Real> &flx) {
	// get proper directional dependence 
  int ivy  = IVX  + ((ivx -IVX )+1)%3;
  int ivz  = IVX  + ((ivx -IVX )+2)%3;
	// diagnol comp 
	int ip11 = IP11 + ((ivx -IVX )  )%3; 
	int ip22 = IP11 + ((ivx -IVX )+1)%3;
	int ip33 = IP11 + ((ivx -IVX )+2)%3; 
	// off-diag comp 
	int ip13 = IP12 + ((ip12-IP12)+1)%3;
	int ip23 = IP12 + ((ip12-IP12)+2)%3; 

  Real wli[(NCLESS)],wri[(NCLESS)],wroe[(NCLESS)];
  Real fl[(NCLESS)],fr[(NCLESS)],flxi[(NCLESS)];

  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {

//--- Step 1.  Load L/R states into local variables

    wli[IDN ]=wl(IDN ,k,j,i);
    wli[IVX ]=wl(ivx ,k,j,i);
    wli[IVY ]=wl(ivy ,k,j,i);
    wli[IVZ ]=wl(ivz ,k,j,i);
		wli[IP11]=wl(ip11,k,j,i);
		wli[IP22]=wl(ip22,k,j,i);
		wli[IP33]=wl(ip33,k,j,i);
		wli[IP12]=wl(ip12,k,j,i);
		wli[IP13]=wl(ip13,k,j,i);
		wli[IP23]=wl(ip23,k,j,i);

    wri[IDN ]=wr(IDN ,k,j,i);
    wri[IVX ]=wr(ivx ,k,j,i);
    wri[IVY ]=wr(ivy ,k,j,i);
    wri[IVZ ]=wr(ivz ,k,j,i);
		wri[IP11]=wr(ip11,k,j,i);
		wri[IP22]=wr(ip22,k,j,i);
		wri[IP33]=wr(ip33,k,j,i);
		wri[IP12]=wr(ip12,k,j,i);
		wri[IP13]=wr(ip13,k,j,i);
		wri[IP23]=wr(ip23,k,j,i);
		

//--- Step 2.  Compute Roe-averaged state

    Real sqrtdl = std::sqrt(wli[IDN]);
    Real sqrtdr = std::sqrt(wri[IDN]);
		Real isqrtdl = 1.0/sqrtdl;
		Real isqrtdr = 1.0/sqrtdr; 
    Real isdlpdr	  = 1.0/(sqrtdl + sqrtdr);
		Real isdlpdrSQR = SQR(isdlpdr);  

		// Stress tensor components are Sij = Pij/d, so Sij sqrt(d) = Pij/sqrt(d)
		// c.f. Eqn. 6.21 of Brown (1996) PhD thesis 

    wroe[IDN ] = sqrtdl*sqrtdr;
    wroe[IVX ] = ( sqrtdl*wli[IVX ] +  sqrtdr*wri[IVX ])*isdlpdr;
    wroe[IVY ] = ( sqrtdl*wli[IVY ] +  sqrtdr*wri[IVY ])*isdlpdr;
    wroe[IVZ ] = ( sqrtdl*wli[IVZ ] +  sqrtdr*wri[IVZ ])*isdlpdr;
		wroe[IP11] = (isqrtdl*wli[IP11] + isqrtdr*wri[IP11])*isdlpdr;
		wroe[IP22] = (isqrtdl*wli[IP22] + isqrtdr*wri[IP22])*isdlpdr;
		wroe[IP33] = (isqrtdl*wli[IP33] + isqrtdr*wri[IP33])*isdlpdr;
		wroe[IP12] = (isqrtdl*wli[IP12] + isqrtdr*wri[IP12])*isdlpdr;
		wroe[IP13] = (isqrtdl*wli[IP13] + isqrtdr*wri[IP13])*isdlpdr;
		wroe[IP23] = (isqrtdl*wli[IP23] + isqrtdr*wri[IP23])*isdlpdr;

		// add ui uj terms of Brown (1996) to Sij, this is necessary b/c the system 
		// is solved in terms of the 'specific enthalpy Hij = 3 Sij + ui uj, but written 
		// in terms of Sij

		wroe[IP11]+=ONE_3RD*sqrtdl*sqrtdr*(wri[IVX]-wli[IVX])*(wri[IVX]-wli[IVX])*isdlpdrSQR;
		wroe[IP22]+=ONE_3RD*sqrtdl*sqrtdr*(wri[IVY]-wli[IVY])*(wri[IVY]-wli[IVY])*isdlpdrSQR;
		wroe[IP33]+=ONE_3RD*sqrtdl*sqrtdr*(wri[IVZ]-wli[IVZ])*(wri[IVZ]-wli[IVZ])*isdlpdrSQR;
		wroe[IP12]+=ONE_3RD*sqrtdl*sqrtdr*(wri[IVX]-wli[IVX])*(wri[IVY]-wli[IVY])*isdlpdrSQR;
		wroe[IP13]+=ONE_3RD*sqrtdl*sqrtdr*(wri[IVX]-wli[IVX])*(wri[IVZ]-wli[IVZ])*isdlpdrSQR;
		wroe[IP23]+=ONE_3RD*sqrtdl*sqrtdr*(wri[IVY]-wli[IVY])*(wri[IVZ]-wli[IVZ])*isdlpdrSQR;


//--- Step 3.  Compute sound speed in L,R, and Roe-averaged states

    Real cl = std::sqrt(3.0*wli[IP11]/wli[IDN]);
    Real cr = std::sqrt(3.0*wri[IP11]/wri[IDN]);
    Real a  = std::sqrt(3.0*wroe[IP11]);  

//--- Step 4. Compute the max/min wave speeds based on L/R and Roe-averaged values

    Real al = std::min((wroe[IVX] - a),(wli[IVX] - cl));
    Real ar = std::max((wroe[IVX] + a),(wri[IVX] + cr));

    Real bp = ar > 0.0 ? ar : 0.0;
    Real bm = al < 0.0 ? al : 0.0;

//-- Step 5. Compute L/R fluxes along the lines bm/bp: F_L - (S_L)U_L; F_R - (S_R)U_R

    Real vxl = wli[IVX] - bm;
    Real vxr = wri[IVX] - bp;

    fl[IDN ] = wli[IDN]*vxl;
    fr[IDN ] = wri[IDN]*vxr;

    fl[IVX ] = wli[IDN]*wli[IVX]*vxl + wli[IP11];
    fr[IVX ] = wri[IDN]*wri[IVX]*vxr + wri[IP11];

    fl[IVY ] = wli[IDN]*wli[IVY]*vxl + wli[IP12];
    fr[IVY ] = wri[IDN]*wri[IVY]*vxr + wri[IP12];

    fl[IVZ ] = wli[IDN]*wli[IVZ]*vxl + wli[IP13];
    fr[IVZ ] = wri[IDN]*wri[IVZ]*vxr + wri[IP13];

		fl[IP11] = (wli[IP11] + wli[IDN]*wli[IVX]*wli[IVX])*vxl + 2.0*wli[IP11]*wli[IVX]; 
		fr[IP11] = (wri[IP11] + wri[IDN]*wri[IVX]*wri[IVX])*vxr + 2.0*wri[IP11]*wri[IVX]; 

		fl[IP22] = (wli[IP22] + wli[IDN]*wli[IVY]*wli[IVY])*vxl + 2.0*wli[IP12]*wli[IVY]; 
		fr[IP22] = (wri[IP22] + wri[IDN]*wri[IVY]*wri[IVY])*vxr + 2.0*wri[IP12]*wri[IVY]; 
		
		fl[IP33] = (wli[IP33] + wli[IDN]*wli[IVZ]*wli[IVZ])*vxl + 2.0*wli[IP13]*wli[IVZ]; 
		fr[IP33] = (wri[IP33] + wri[IDN]*wri[IVZ]*wri[IVZ])*vxr + 2.0*wri[IP13]*wri[IVZ]; 

		fl[IP12] = (wli[IP12] + wli[IDN]*wli[IVX]*wli[IVY])*vxl + ( wli[IP12]*wli[IVX] 
																									          +   wli[IP11]*wli[IVY] );
		fr[IP12] = (wri[IP12] + wri[IDN]*wri[IVX]*wri[IVY])*vxr + ( wri[IP12]*wri[IVX] 
																									          +   wri[IP11]*wri[IVY] );

		fl[IP13] = (wli[IP13] + wli[IDN]*wli[IVX]*wli[IVZ])*vxl + ( wli[IP13]*wli[IVX] 
																									          +   wli[IP11]*wli[IVZ] );
		fr[IP13] = (wri[IP13] + wri[IDN]*wri[IVX]*wri[IVZ])*vxr + ( wri[IP13]*wri[IVX] 
																									          +   wri[IP11]*wri[IVZ] );

		fl[IP23] = (wli[IP23] + wli[IDN]*wli[IVY]*wli[IVZ])*vxl + ( wli[IP13]*wli[IVY] 
																									          +   wli[IP12]*wli[IVZ] );
		fr[IP23] = (wri[IP23] + wri[IDN]*wri[IVY]*wri[IVZ])*vxr + ( wri[IP13]*wri[IVY] 
																									          +   wri[IP12]*wri[IVZ] );


//--- Step 6. Compute the HLLE flux at interface.

    Real tmp=0.0;
    if (bp != bm) tmp = 0.5*(bp + bm)/(bp - bm);

    flxi[IDN ] = 0.5*(fl[IDN ]+fr[IDN ]) + (fl[IDN ]-fr[IDN ])*tmp;
    flxi[IVX ] = 0.5*(fl[IVX ]+fr[IVX ]) + (fl[IVX ]-fr[IVX ])*tmp;
    flxi[IVY ] = 0.5*(fl[IVY ]+fr[IVY ]) + (fl[IVY ]-fr[IVY ])*tmp;
    flxi[IVZ ] = 0.5*(fl[IVZ ]+fr[IVZ ]) + (fl[IVZ ]-fr[IVZ ])*tmp;
		flxi[IP11] = 0.5*(fl[IP11]+fr[IP11]) + (fl[IP11]-fr[IP11])*tmp;
		flxi[IP22] = 0.5*(fl[IP22]+fr[IP22]) + (fl[IP22]-fr[IP22])*tmp;
		flxi[IP33] = 0.5*(fl[IP33]+fr[IP33]) + (fl[IP33]-fr[IP33])*tmp;
		flxi[IP12] = 0.5*(fl[IP12]+fr[IP12]) + (fl[IP12]-fr[IP12])*tmp;
		flxi[IP13] = 0.5*(fl[IP13]+fr[IP13]) + (fl[IP13]-fr[IP13])*tmp;
		flxi[IP23] = 0.5*(fl[IP23]+fr[IP23]) + (fl[IP23]-fr[IP23])*tmp;

    flx(IDN ,k,j,i) = flxi[IDN ];
    flx(ivx ,k,j,i) = flxi[IVX ];
    flx(ivy ,k,j,i) = flxi[IVY ];
    flx(ivz ,k,j,i) = flxi[IVZ ];
		flx(ip11,k,j,i) = flxi[IP11];
		flx(ip22,k,j,i) = flxi[IP22];
		flx(ip33,k,j,i) = flxi[IP33];
		flx(ip12,k,j,i) = flxi[IP12];
		flx(ip13,k,j,i) = flxi[IP13];
		flx(ip23,k,j,i) = flxi[IP23];

  }
  }}

  return;
}
