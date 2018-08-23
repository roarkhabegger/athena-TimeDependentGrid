//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file  roe_cl.cpp
//  \brief Roe's linearized Riemann solver for the collisionless system.
//
// Computes 1D fluxes using Roe's linearization.  When Roe's method fails because of
// negative density or pressure in the intermediate states, LLF fluxes are used instead.
//
// REFERENCES:
// - P. Roe, "Approximate Riemann solvers, parameter vectors, and difference schemes",
//   JCP, 43, 357 (1981).

// C/C++ headers
#include <algorithm>  // max()
#include <cmath>      // sqrt()

// Athena++ headers
#include "../cless.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../eos/eos.hpp"

// prototype for function to compute eigenvalues and eigenvectors of Roe's matrix A
inline void RoeEigensystemCL(const Real wroe[], Real eigenvalues[],
  Real right_eigenmatrix[][(NWAVECL)], Real left_eigenmatrix[][(NWAVECL)]);
// set to 1 to test intermediate states (currently not working)
#define TEST_INTERMEDIATE_STATES 0 

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
  
	Real wli[NWAVECL],wri[NWAVECL];
	Real wroe[NWAVECL],fl[NWAVECL],fr[NWAVECL],flxi[NWAVECL];

  Real coeff[NWAVECL];
  Real ev[NWAVECL],rem[NWAVECL][NWAVECL],lem[NWAVECL][NWAVECL];
  Real du[NWAVECL],a[NWAVECL],u[NWAVECL];

  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
#pragma omg simd
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

		// get some necessary conserved variables
		Real e11l = wli[IP11] + wli[IDN]*wli[IVX]*wli[IVX];
		Real e11r = wri[IP11] + wri[IDN]*wri[IVX]*wri[IVX];
		Real e22l = wli[IP22] + wli[IDN]*wli[IVY]*wli[IVY];
		Real e22r = wri[IP22] + wri[IDN]*wri[IVY]*wri[IVY];
		Real e33l = wli[IP33] + wli[IDN]*wli[IVZ]*wli[IVZ];
		Real e33r = wri[IP33] + wri[IDN]*wri[IVZ]*wri[IVZ];
		Real e12l = wli[IP12] + wli[IDN]*wli[IVX]*wli[IVY];
		Real e12r = wri[IP12] + wri[IDN]*wri[IVX]*wri[IVY];
		Real e13l = wli[IP13] + wli[IDN]*wli[IVX]*wli[IVZ];
		Real e13r = wri[IP13] + wri[IDN]*wri[IVX]*wri[IVZ];
		Real e23l = wli[IP23] + wli[IDN]*wli[IVY]*wli[IVZ];
		Real e23r = wri[IP23] + wri[IDN]*wri[IVY]*wri[IVZ];

//--- Step 2.  Compute Roe-averaged data from left- and right-states

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

		//wroe[IP11]+=ONE_3RD*sqrtdl*sqrtdr*(wri[IVX]-wli[IVX])*(wri[IVX]-wli[IVX])*isdlpdrSQR;
		//wroe[IP22]+=ONE_3RD*sqrtdl*sqrtdr*(wri[IVY]-wli[IVY])*(wri[IVY]-wli[IVY])*isdlpdrSQR;
		//wroe[IP33]+=ONE_3RD*sqrtdl*sqrtdr*(wri[IVZ]-wli[IVZ])*(wri[IVZ]-wli[IVZ])*isdlpdrSQR;
		//wroe[IP12]+=ONE_3RD*sqrtdl*sqrtdr*(wri[IVX]-wli[IVX])*(wri[IVY]-wli[IVY])*isdlpdrSQR;
		//wroe[IP13]+=ONE_3RD*sqrtdl*sqrtdr*(wri[IVX]-wli[IVX])*(wri[IVZ]-wli[IVZ])*isdlpdrSQR;
		//wroe[IP23]+=ONE_3RD*sqrtdl*sqrtdr*(wri[IVY]-wli[IVY])*(wri[IVZ]-wli[IVZ])*isdlpdrSQR;

//--- Step 3.  Compute eigenvalues and eigenmatrices using Roe-averaged values

    RoeEigensystemCL(wroe,ev,rem,lem);

//--- Step 4.  Compute L/R fluxes

		Real vxl = wli[IVX];
		Real vxr = wri[IVX]; 
    
		fl[IDN ] = wli[IDN]*vxl;
    fr[IDN ] = wli[IDN]*vxr;

    fl[IVX ] = wli[IDN]*wli[IVX]*vxl + wli[IP11];
    fr[IVX ] = wri[IDN]*wri[IVX]*vxr + wri[IP11];

    fl[IVY ] = wli[IDN]*wli[IVY]*vxl + wli[IP12];
    fr[IVY ] = wri[IDN]*wri[IVY]*vxr + wri[IP12];

    fl[IVZ ] = wli[IDN]*wli[IVZ]*vxl + wli[IP13];
    fr[IVZ ] = wri[IDN]*wri[IVZ]*vxr + wri[IP13];

		fl[IP11] = e11l*vxl + 2.0*wli[IP11]*wli[IVX]; 
		fr[IP11] = e11r*vxr + 2.0*wri[IP11]*wri[IVX]; 

		fl[IP22] = e22l*vxl + 2.0*wli[IP12]*wli[IVY]; 
		fr[IP22] = e22r*vxr + 2.0*wri[IP12]*wri[IVY]; 
		
		fl[IP33] = e33l*vxl + 2.0*wli[IP13]*wli[IVZ]; 
		fr[IP33] = e33r*vxr + 2.0*wri[IP13]*wri[IVZ]; 

		fl[IP12] = e12l*vxl + ( wli[IP12]*wli[IVX] + wli[IP11]*wli[IVY] );
		fr[IP12] = e12r*vxr + ( wri[IP12]*wri[IVX] + wri[IP11]*wri[IVY] );

		fl[IP13] = e13l*vxl + ( wli[IP13]*wli[IVX] + wli[IP11]*wli[IVZ] );
		fr[IP13] = e13r*vxr + ( wri[IP13]*wri[IVX] + wri[IP11]*wri[IVZ] );

		fl[IP23] = e23l*vxl + ( wli[IP13]*wli[IVY] + wli[IP12]*wli[IVZ] );
		fr[IP23] = e23r*vxr + ( wri[IP13]*wri[IVY] + wri[IP12]*wri[IVZ] );


//--- Step 5.  Compute projection of dU onto L eigenvectors ("vector A")

		du[IDN ] = wri[IDN]          - wli[IDN];
		du[IVX ] = wri[IDN]*wri[IVX] - wli[IDN]*wli[IVX]; 
		du[IVY ] = wri[IDN]*wri[IVY] - wli[IDN]*wli[IVY]; 
		du[IVZ ] = wri[IDN]*wri[IVZ] - wli[IDN]*wli[IVZ]; 
		du[IP11] = e11r - e11l;  
		du[IP22] = e22r - e22l; 
		du[IP33] = e33r - e33l; 
		du[IP12] = e12r - e12l;  
		du[IP13] = e13r - e13l; 
		du[IP23] = e23r - e23l; 
		
		a[IDN ]  = lem[IDN ][IDN ]*du[IDN ];
		a[IDN ] += lem[IDN ][IVX ]*du[IVX ];
		a[IDN ] += lem[IDN ][IVY ]*du[IVY ];
		a[IDN ] += lem[IDN ][IVZ ]*du[IVZ ];
		a[IDN ] += lem[IDN ][IP11]*du[IP11];
		a[IDN ] += lem[IDN ][IP22]*du[IP22];
		a[IDN ] += lem[IDN ][IP33]*du[IP33];
		a[IDN ] += lem[IDN ][IP12]*du[IP12]; 
		a[IDN ] += lem[IDN ][IP13]*du[IP13];
		a[IDN ] += lem[IDN ][IP23]*du[IP23];

		a[IVX ]  = lem[IVX ][IDN ]*du[IDN ];
		a[IVX ] += lem[IVX ][IVX ]*du[IVX ];
		a[IVX ] += lem[IVX ][IVY ]*du[IVY ];
		a[IVX ] += lem[IVX ][IVZ ]*du[IVZ ];
		a[IVX ] += lem[IVX ][IP11]*du[IP11];
		a[IVX ] += lem[IVX ][IP22]*du[IP22];
		a[IVX ] += lem[IVX ][IP33]*du[IP33];
		a[IVX ] += lem[IVX ][IP12]*du[IP12]; 
		a[IVX ] += lem[IVX ][IP13]*du[IP13];
		a[IVX ] += lem[IVX ][IP23]*du[IP23];

		a[IVY ]  = lem[IVY ][IDN ]*du[IDN ];
		a[IVY ] += lem[IVY ][IVX ]*du[IVX ];
		a[IVY ] += lem[IVY ][IVY ]*du[IVY ];
		a[IVY ] += lem[IVY ][IVZ ]*du[IVZ ];
		a[IVY ] += lem[IVY ][IP11]*du[IP11];
		a[IVY ] += lem[IVY ][IP22]*du[IP22];
		a[IVY ] += lem[IVY ][IP33]*du[IP33];
		a[IVY ] += lem[IVY ][IP12]*du[IP12]; 
		a[IVY ] += lem[IVY ][IP13]*du[IP13];
		a[IVY ] += lem[IVY ][IP23]*du[IP23];

		a[IVZ ]  = lem[IVZ ][IDN ]*du[IDN ];
		a[IVZ ] += lem[IVZ ][IVX ]*du[IVX ];
		a[IVZ ] += lem[IVZ ][IVY ]*du[IVY ];
		a[IVZ ] += lem[IVZ ][IVZ ]*du[IVZ ];
		a[IVZ ] += lem[IVZ ][IP11]*du[IP11];
		a[IVZ ] += lem[IVZ ][IP22]*du[IP22];
		a[IVZ ] += lem[IVZ ][IP33]*du[IP33];
		a[IVZ ] += lem[IVZ ][IP12]*du[IP12]; 
		a[IVZ ] += lem[IVZ ][IP13]*du[IP13];
		a[IVZ ] += lem[IVZ ][IP23]*du[IP23];

		a[IP11]  = lem[IP11][IDN ]*du[IDN ];
		a[IP11] += lem[IP11][IVX ]*du[IVX ];
		a[IP11] += lem[IP11][IVY ]*du[IVY ];
		a[IP11] += lem[IP11][IVZ ]*du[IVZ ];
		a[IP11] += lem[IP11][IP11]*du[IP11];
		a[IP11] += lem[IP11][IP22]*du[IP22];
		a[IP11] += lem[IP11][IP33]*du[IP33];
		a[IP11] += lem[IP11][IP12]*du[IP12]; 
		a[IP11] += lem[IP11][IP13]*du[IP13];
		a[IP11] += lem[IP11][IP23]*du[IP23];

		a[IP22]  = lem[IP22][IDN ]*du[IDN ];
		a[IP22] += lem[IP22][IVX ]*du[IVX ];
		a[IP22] += lem[IP22][IVY ]*du[IVY ];
		a[IP22] += lem[IP22][IVZ ]*du[IVZ ];
		a[IP22] += lem[IP22][IP11]*du[IP11];
		a[IP22] += lem[IP22][IP22]*du[IP22];
		a[IP22] += lem[IP22][IP33]*du[IP33];
		a[IP22] += lem[IP22][IP12]*du[IP12]; 
		a[IP22] += lem[IP22][IP13]*du[IP13];
		a[IP22] += lem[IP22][IP23]*du[IP23];

		a[IP33]  = lem[IP33][IDN ]*du[IDN ];
		a[IP33] += lem[IP33][IVX ]*du[IVX ];
		a[IP33] += lem[IP33][IVY ]*du[IVY ];
		a[IP33] += lem[IP33][IVZ ]*du[IVZ ];
		a[IP33] += lem[IP33][IP11]*du[IP11];
		a[IP33] += lem[IP33][IP22]*du[IP22];
		a[IP33] += lem[IP33][IP33]*du[IP33];
		a[IP33] += lem[IP33][IP12]*du[IP12]; 
		a[IP33] += lem[IP33][IP13]*du[IP13];
		a[IP33] += lem[IP33][IP23]*du[IP23];

		a[IP12]  = lem[IP12][IDN ]*du[IDN ];
		a[IP12] += lem[IP12][IVX ]*du[IVX ];
		a[IP12] += lem[IP12][IVY ]*du[IVY ];
		a[IP12] += lem[IP12][IVZ ]*du[IVZ ];
		a[IP12] += lem[IP12][IP11]*du[IP11];
		a[IP12] += lem[IP12][IP22]*du[IP22];
		a[IP12] += lem[IP12][IP33]*du[IP33];
		a[IP12] += lem[IP12][IP12]*du[IP12]; 
		a[IP12] += lem[IP12][IP13]*du[IP13];
		a[IP12] += lem[IP12][IP23]*du[IP23];

		a[IP13]  = lem[IP13][IDN ]*du[IDN ];
		a[IP13] += lem[IP13][IVX ]*du[IVX ];
		a[IP13] += lem[IP13][IVY ]*du[IVY ];
		a[IP13] += lem[IP13][IVZ ]*du[IVZ ];
		a[IP13] += lem[IP13][IP11]*du[IP11];
		a[IP13] += lem[IP13][IP22]*du[IP22];
		a[IP13] += lem[IP13][IP33]*du[IP33];
		a[IP13] += lem[IP13][IP12]*du[IP12]; 
		a[IP13] += lem[IP13][IP13]*du[IP13];
		a[IP13] += lem[IP13][IP23]*du[IP23];

		a[IP23]  = lem[IP23][IDN ]*du[IDN ];
		a[IP23] += lem[IP23][IVX ]*du[IVX ];
		a[IP23] += lem[IP23][IVY ]*du[IVY ];
		a[IP23] += lem[IP23][IVZ ]*du[IVZ ];
		a[IP23] += lem[IP23][IP11]*du[IP11];
		a[IP23] += lem[IP23][IP22]*du[IP22];
		a[IP23] += lem[IP23][IP33]*du[IP33];
		a[IP23] += lem[IP23][IP12]*du[IP12]; 
		a[IP23] += lem[IP23][IP13]*du[IP13];
		a[IP23] += lem[IP23][IP23]*du[IP23];

//--- Step 6.  Check that the density and pressure in the intermediate states are
// positive.  If not, set a flag that will be checked below.

    int llf_flag = 0;
		if (TEST_INTERMEDIATE_STATES) {	
			u[IDN ] = wli[IDN]; 
			u[IVX ] = wli[IDN]*wli[IVX]; 
			u[IVY ] = wli[IDN]*wli[IVY];
			u[IVZ ] = wli[IDN]*wli[IVZ]; 
			u[IE11] = e11l;
			u[IE22] = e22l;
			u[IE33] = e33l;
			u[IE12] = e12l;
			u[IE13] = e13l;
			u[IE23] = e23l;

			// jump across wave[0] 
			u[IDN ] += a[0]*rem[IDN ][0];
			if (u[IDN] < 0.0) llf_flag=1; 
			u[IVX ] += a[0]*rem[IVX ][0];
			u[IVY ] += a[0]*rem[IVY ][0];
			u[IVZ ] += a[0]*rem[IVZ ][0];
			u[IE11] += a[0]*rem[IP11][0];
			u[IE22] += a[0]*rem[IP22][0];
			u[IE33] += a[0]*rem[IP33][0];
			Real p11 = u[IE11] - SQR(u[IVX])/u[IDN];
			Real p22 = u[IE22] - SQR(u[IVY])/u[IDN];
			Real p33 = u[IE33] - SQR(u[IVZ])/u[IDN];

			if (p11 < 0.0) llf_flag=2;
			if (p22 < 0.0) llf_flag=3;
			if (p33 < 0.0) llf_flag=4; 

			// jump across wave[1] 
			u[IDN ] += a[1]*rem[IDN ][1];
			if (u[IDN] < 0.0) llf_flag=1; 
			u[IVX ] += a[1]*rem[IVX ][1];
			u[IVY ] += a[1]*rem[IVY ][1];
			u[IVZ ] += a[1]*rem[IVZ ][1];
			u[IE11] += a[1]*rem[IP11][1];
			u[IE22] += a[1]*rem[IP22][1];
			u[IE33] += a[1]*rem[IP33][1];
			p11 = u[IE11] - SQR(u[IVX])/u[IDN];
			p22 = u[IE22] - SQR(u[IVY])/u[IDN];
			p33 = u[IE33] - SQR(u[IVZ])/u[IDN];

			if (p11 < 0.0) llf_flag=2;
			if (p22 < 0.0) llf_flag=3;
			if (p33 < 0.0) llf_flag=4; 

			// jump across wave[2] 
			u[IDN ] += a[2]*rem[IDN ][2];
			if (u[IDN] < 0.0) llf_flag=1; 
			u[IVX ] += a[2]*rem[IVX ][2];
			u[IVY ] += a[2]*rem[IVY ][2];
			u[IVZ ] += a[2]*rem[IVZ ][2];
			u[IE11] += a[2]*rem[IP11][2];
			u[IE22] += a[2]*rem[IP22][2];
			u[IE33] += a[2]*rem[IP33][2];
			p11 = u[IE11] - SQR(u[IVX])/u[IDN];
			p22 = u[IE22] - SQR(u[IVY])/u[IDN];
			p33 = u[IE33] - SQR(u[IVZ])/u[IDN];

			if (p11 < 0.0) llf_flag=2;
			if (p22 < 0.0) llf_flag=3;
			if (p33 < 0.0) llf_flag=4; 

			// jump across wave[3] 
			u[IDN ] += a[3]*rem[IDN ][3];
			if (u[IDN] < 0.0) llf_flag=1; 
			u[IVX ] += a[3]*rem[IVX ][3];
			u[IVY ] += a[3]*rem[IVY ][3];
			u[IVZ ] += a[3]*rem[IVZ ][3];
			u[IE11] += a[3]*rem[IP11][3];
			u[IE22] += a[3]*rem[IP22][3];
			u[IE33] += a[3]*rem[IP33][3];
			p11 = u[IE11] - SQR(u[IVX])/u[IDN];
			p22 = u[IE22] - SQR(u[IVY])/u[IDN];
			p33 = u[IE33] - SQR(u[IVZ])/u[IDN];

			if (p11 < 0.0) llf_flag=2;
			if (p22 < 0.0) llf_flag=3;
			if (p33 < 0.0) llf_flag=4; 
		}

//--- Step 7.  Compute Roe flux
		coeff[IDN ] = 0.5*fabs(ev[IDN ])*a[IDN ];
		coeff[IVX ] = 0.5*fabs(ev[IVX ])*a[IVX ];
		coeff[IVY ] = 0.5*fabs(ev[IVY ])*a[IVY ];
		coeff[IVZ ] = 0.5*fabs(ev[IVZ ])*a[IVZ ];
		coeff[IP11] = 0.5*fabs(ev[IP11])*a[IP11];
		coeff[IP22] = 0.5*fabs(ev[IP22])*a[IP22];
		coeff[IP33] = 0.5*fabs(ev[IP33])*a[IP33];
		coeff[IP12] = 0.5*fabs(ev[IP12])*a[IP12];
		coeff[IP13] = 0.5*fabs(ev[IP13])*a[IP13];
		coeff[IP23] = 0.5*fabs(ev[IP23])*a[IP23];

		flxi[IDN ]  = 0.5*(fl[IDN] + fr[IDN]);
		flxi[IDN ] -= rem[IDN ][IDN ]*coeff[IDN ];
		flxi[IDN ] -= rem[IDN ][IVX ]*coeff[IVX ];
		flxi[IDN ] -= rem[IDN ][IVY ]*coeff[IVY ];
		flxi[IDN ] -= rem[IDN ][IVZ ]*coeff[IVZ ];
		flxi[IDN ] -= rem[IDN ][IP11]*coeff[IP11];
		flxi[IDN ] -= rem[IDN ][IP22]*coeff[IP22];
		flxi[IDN ] -= rem[IDN ][IP33]*coeff[IP33];
		flxi[IDN ] -= rem[IDN ][IP12]*coeff[IP12]; 
		flxi[IDN ] -= rem[IDN ][IP13]*coeff[IP13];
		flxi[IDN ] -= rem[IDN ][IP23]*coeff[IP23];

		flxi[IVX ]  = 0.5*(fl[IVX] + fr[IVX]); 
		flxi[IVX ] -= rem[IVX ][IDN ]*coeff[IDN ];
		flxi[IVX ] -= rem[IVX ][IVX ]*coeff[IVX ];
		flxi[IVX ] -= rem[IVX ][IVY ]*coeff[IVY ];
		flxi[IVX ] -= rem[IVX ][IVZ ]*coeff[IVZ ];
		flxi[IVX ] -= rem[IVX ][IP11]*coeff[IP11];
		flxi[IVX ] -= rem[IVX ][IP22]*coeff[IP22];
		flxi[IVX ] -= rem[IVX ][IP33]*coeff[IP33];
		flxi[IVX ] -= rem[IVX ][IP12]*coeff[IP12]; 
		flxi[IVX ] -= rem[IVX ][IP13]*coeff[IP13];
		flxi[IVX ] -= rem[IVX ][IP23]*coeff[IP23];

		flxi[IVY ]  = 0.5*(fl[IVY] + fr[IVY]);
		flxi[IVY ] -= rem[IVY ][IDN ]*coeff[IDN ];
		flxi[IVY ] -= rem[IVY ][IVX ]*coeff[IVX ];
		flxi[IVY ] -= rem[IVY ][IVY ]*coeff[IVY ];
		flxi[IVY ] -= rem[IVY ][IVZ ]*coeff[IVZ ];
		flxi[IVY ] -= rem[IVY ][IP11]*coeff[IP11];
		flxi[IVY ] -= rem[IVY ][IP22]*coeff[IP22];
		flxi[IVY ] -= rem[IVY ][IP33]*coeff[IP33];
		flxi[IVY ] -= rem[IVY ][IP12]*coeff[IP12]; 
		flxi[IVY ] -= rem[IVY ][IP13]*coeff[IP13];
		flxi[IVY ] -= rem[IVY ][IP23]*coeff[IP23];

		flxi[IVZ ]  = 0.5*(fl[IVZ] + fr[IVZ]);
		flxi[IVZ ] -= rem[IVZ ][IDN ]*coeff[IDN ];
		flxi[IVZ ] -= rem[IVZ ][IVX ]*coeff[IVX ];
		flxi[IVZ ] -= rem[IVZ ][IVY ]*coeff[IVY ];
		flxi[IVZ ] -= rem[IVZ ][IVZ ]*coeff[IVZ ];
		flxi[IVZ ] -= rem[IVZ ][IP11]*coeff[IP11];
		flxi[IVZ ] -= rem[IVZ ][IP22]*coeff[IP22];
		flxi[IVZ ] -= rem[IVZ ][IP33]*coeff[IP33];
		flxi[IVZ ] -= rem[IVZ ][IP12]*coeff[IP12]; 
		flxi[IVZ ] -= rem[IVZ ][IP13]*coeff[IP13];
		flxi[IVZ ] -= rem[IVZ ][IP23]*coeff[IP23];

		flxi[IP11]  = 0.5*(fl[IP11] + fr[IP11]);
		flxi[IP11] -= rem[IP11][IDN ]*coeff[IDN ];
		flxi[IP11] -= rem[IP11][IVX ]*coeff[IVX ];
		flxi[IP11] -= rem[IP11][IVY ]*coeff[IVY ];
		flxi[IP11] -= rem[IP11][IVZ ]*coeff[IVZ ];
		flxi[IP11] -= rem[IP11][IP11]*coeff[IP11];
		flxi[IP11] -= rem[IP11][IP22]*coeff[IP22];
		flxi[IP11] -= rem[IP11][IP33]*coeff[IP33];
		flxi[IP11] -= rem[IP11][IP12]*coeff[IP12]; 
		flxi[IP11] -= rem[IP11][IP13]*coeff[IP13];
		flxi[IP11] -= rem[IP11][IP23]*coeff[IP23];

		flxi[IP22]  = 0.5*(fl[IP22] + fr[IP22]);
		flxi[IP22] -= rem[IP22][IDN ]*coeff[IDN ];
		flxi[IP22] -= rem[IP22][IVX ]*coeff[IVX ];
		flxi[IP22] -= rem[IP22][IVY ]*coeff[IVY ];
		flxi[IP22] -= rem[IP22][IVZ ]*coeff[IVZ ];
		flxi[IP22] -= rem[IP22][IP11]*coeff[IP11];
		flxi[IP22] -= rem[IP22][IP22]*coeff[IP22];
		flxi[IP22] -= rem[IP22][IP33]*coeff[IP33];
		flxi[IP22] -= rem[IP22][IP12]*coeff[IP12]; 
		flxi[IP22] -= rem[IP22][IP13]*coeff[IP13];
		flxi[IP22] -= rem[IP22][IP23]*coeff[IP23];

		flxi[IP33]  = 0.5*(fl[IP33] + fr[IP33]);
		flxi[IP33] -= rem[IP33][IDN ]*coeff[IDN ];
		flxi[IP33] -= rem[IP33][IVX ]*coeff[IVX ];
		flxi[IP33] -= rem[IP33][IVY ]*coeff[IVY ];
		flxi[IP33] -= rem[IP33][IVZ ]*coeff[IVZ ];
		flxi[IP33] -= rem[IP33][IP11]*coeff[IP11];
		flxi[IP33] -= rem[IP33][IP22]*coeff[IP22];
		flxi[IP33] -= rem[IP33][IP33]*coeff[IP33];
		flxi[IP33] -= rem[IP33][IP12]*coeff[IP12]; 
		flxi[IP33] -= rem[IP33][IP13]*coeff[IP13];
		flxi[IP33] -= rem[IP33][IP23]*coeff[IP23];

		flxi[IP12]  = 0.5*(fl[IP12] + fr[IP12]);
		flxi[IP12] -= rem[IP12][IDN ]*coeff[IDN ];
		flxi[IP12] -= rem[IP12][IVX ]*coeff[IVX ];
		flxi[IP12] -= rem[IP12][IVY ]*coeff[IVY ];
		flxi[IP12] -= rem[IP12][IVZ ]*coeff[IVZ ];
		flxi[IP12] -= rem[IP12][IP11]*coeff[IP11];
		flxi[IP12] -= rem[IP12][IP22]*coeff[IP22];
		flxi[IP12] -= rem[IP12][IP33]*coeff[IP33];
		flxi[IP12] -= rem[IP12][IP12]*coeff[IP12]; 
		flxi[IP12] -= rem[IP12][IP13]*coeff[IP13];
		flxi[IP12] -= rem[IP12][IP23]*coeff[IP23];

		flxi[IP13]  = 0.5*(fl[IP13] + fr[IP13]);
		flxi[IP13] -= rem[IP13][IDN ]*coeff[IDN ];
		flxi[IP13] -= rem[IP13][IVX ]*coeff[IVX ];
		flxi[IP13] -= rem[IP13][IVY ]*coeff[IVY ];
		flxi[IP13] -= rem[IP13][IVZ ]*coeff[IVZ ];
		flxi[IP13] -= rem[IP13][IP11]*coeff[IP11];
		flxi[IP13] -= rem[IP13][IP22]*coeff[IP22];
		flxi[IP13] -= rem[IP13][IP33]*coeff[IP33];
		flxi[IP13] -= rem[IP13][IP12]*coeff[IP12]; 
		flxi[IP13] -= rem[IP13][IP13]*coeff[IP13];
		flxi[IP13] -= rem[IP13][IP23]*coeff[IP23];

		flxi[IP23]  = 0.5*(fl[IP23] + fr[IP23]);
		flxi[IP23] -= rem[IP23][IDN ]*coeff[IDN ];
		flxi[IP23] -= rem[IP23][IVX ]*coeff[IVX ];
		flxi[IP23] -= rem[IP23][IVY ]*coeff[IVY ];
		flxi[IP23] -= rem[IP23][IVZ ]*coeff[IVZ ];
		flxi[IP23] -= rem[IP23][IP11]*coeff[IP11];
		flxi[IP23] -= rem[IP23][IP22]*coeff[IP22];
		flxi[IP23] -= rem[IP23][IP33]*coeff[IP33];
		flxi[IP23] -= rem[IP23][IP12]*coeff[IP12]; 
		flxi[IP23] -= rem[IP23][IP13]*coeff[IP13];
		flxi[IP23] -= rem[IP23][IP23]*coeff[IP23];


//--- Step 8.  Overwrite with upwind flux if flow is supersonic
		if (ev[0] >= 0.0) {
			flxi[IDN ] = fl[IDN ];
			flxi[IVX ] = fl[IVX ];  
			flxi[IVY ] = fl[IVY ];
			flxi[IVZ ] = fl[IVZ ];
			flxi[IP11] = fl[IP11];
			flxi[IP22] = fl[IP22];
			flxi[IP33] = fl[IP33];
			flxi[IP12] = fl[IP12];
			flxi[IP13] = fl[IP13];
			flxi[IP23] = fl[IP23];
		}
		if (ev[NWAVECL-1] <= 0.0) {
			flxi[IDN ] = fr[IDN ];
			flxi[IVX ] = fr[IVX ];  
			flxi[IVY ] = fr[IVY ];
			flxi[IVZ ] = fr[IVZ ];
			flxi[IP11] = fr[IP11];
			flxi[IP22] = fr[IP22];
			flxi[IP33] = fr[IP33];
			flxi[IP12] = fr[IP12];
			flxi[IP13] = fr[IP13];
			flxi[IP23] = fr[IP23];
		}

//--- Step 9.  Overwrite with LLF flux if any of intermediate states are negative
		if (llf_flag != 0) {
			Real a = std::max(fabs(ev[0]), fabs(ev[NWAVECL-1]));

			flxi[IDN ] = 0.5*(fl[IDN ] + fr[IDN ]) - a*du[IDN ];
			flxi[IVX ] = 0.5*(fl[IVX ] + fr[IVX ]) - a*du[IVX ];
			flxi[IVY ] = 0.5*(fl[IVY ] + fr[IVY ]) - a*du[IVY ];
			flxi[IVZ ] = 0.5*(fl[IVZ ] + fr[IVZ ]) - a*du[IVZ ];
			flxi[IP11] = 0.5*(fl[IP11] + fr[IP11]) - a*du[IP11];
			flxi[IP22] = 0.5*(fl[IP22] + fr[IP22]) - a*du[IP22];
			flxi[IP33] = 0.5*(fl[IP33] + fr[IP33]) - a*du[IP33];
			flxi[IP12] = 0.5*(fl[IP12] + fr[IP12]) - a*du[IP12];
			flxi[IP13] = 0.5*(fl[IP13] + fr[IP13]) - a*du[IP13];
			flxi[IP23] = 0.5*(fl[IP23] + fr[IP23]) - a*du[IP23];
		}

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

//----------------------------------------------------------------------------------------
// \!fn RoeEigensystemCL()
// \brief computes eigenvalues and eigenvectors for hydrodynamics
//
// PURPOSE: Functions to evaluate the eigenvalues, and left- and right-eigenvectors of
// "Roe's matrix A" for the linearized system in the CONSERVED variables, i.e.
// U,t = AU,x, where U=(d,d*vx,d*vy,d*vz,[E],[By,Bz]). The eigenvalues are returned
// through the argument list as a vector of length NWAVECL.  The eigenvectors are returned
// as matrices of size (NWAVECL)x(NWAVECL), with right-eigenvectors stored as COLUMNS
// (so R_i = right_eigenmatrix[*][i]), and left-eigenvectors stored as ROWS
// (so L_i = left_eigenmatrix[i][*]).
//     - Input: v1,v2,v3,h = Roe averaged velocities and enthalpy
//     - Output: eigenvalues[], right_eigenmatrix[][], left_eigenmatrix[][];
//
// REFERENCES:
// - P. Cargo & G. Gallice, "Roe matrices for ideal MHD and systematic construction of
//   Roe matrices for systems of conservation laws", JCP, 136, 446 (1997)
//
// - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new code for
//   astrophysical MHD", ApJS, (2008), Appendix B  Equation numbers refer to this paper.

inline void RoeEigensystemCL(const Real wroe[], Real eigenvalues[],
  Real right_eigenmatrix[][(NWAVECL)], Real left_eigenmatrix[][(NWAVECL)]) {
	// get all variables
	Real d   = wroe[IDN ];
	Real v1  = wroe[IVX ];
	Real v2  = wroe[IVY ];
	Real v3  = wroe[IVZ ];
	Real S11 = wroe[IP11];
	Real S22 = wroe[IP22];
	Real S33 = wroe[IP33];
	Real S12 = wroe[IP12];
	Real S13 = wroe[IP13];
	Real S23 = wroe[IP23];
	
	Real a		 = std::sqrt(S11);
	Real sqt3  = std::sqrt(3.0);
	Real sqt3a = sqt3*a;
	Real a2		 = S11;
	Real inva  = 1.0/a;
	Real inva2 = 1.0/a2; 
 
	// compute eigenvalues 
	eigenvalues[0] = v1 - sqt3a;
	eigenvalues[1] = v1 - a;
	eigenvalues[2] = v1 - a;
	eigenvalues[3] = v1;
	eigenvalues[4] = v1;
	eigenvalues[5] = v1;
	eigenvalues[6] = v1;
	eigenvalues[7] = v1 + a;
	eigenvalues[8] = v1 + a;
	eigenvalues[9] = v1 + sqt3a;  

	// Right-eigenvectors, stored as COLUMNS 
  right_eigenmatrix[0][0] = 1.0;
  right_eigenmatrix[1][0] = v1 - sqt3a;
  right_eigenmatrix[2][0] = v2 - sqt3*S12*inva;
  right_eigenmatrix[3][0] = v3 - sqt3*S13*inva;
  right_eigenmatrix[4][0] = SQR(v1 - sqt3a);
  right_eigenmatrix[5][0] = SQR(v2) - 2.0*sqt3*v2*S12*inva + S22 + 2.0*SQR(S12*inva);
  right_eigenmatrix[6][0] = SQR(v3) - 2.0*sqt3*v3*S13*inva + S33 + 2.0*SQR(S13*inva);
  right_eigenmatrix[7][0] = (v2-sqt3*S12*inva)*(v1-sqt3a);
  right_eigenmatrix[8][0] = (v3-sqt3*S13*inva)*(v1-sqt3a);
  right_eigenmatrix[9][0] = v2*v3 + 2.0*S12*S13*inva2 + S23 - sqt3*(S13*v2+S12*v3)*inva;

	right_eigenmatrix[0][1] = 0.0; 
	right_eigenmatrix[1][1] = 0.0;
  right_eigenmatrix[2][1] = inva;
	right_eigenmatrix[3][1] = 0.0; 
	right_eigenmatrix[4][1] = 0.0; 
  right_eigenmatrix[5][1] = 2.0*(a*v2-S12)*inva2;
	right_eigenmatrix[6][1] = 0.0; 
  right_eigenmatrix[7][1] = v1*inva - 1.0;
	right_eigenmatrix[8][1] = 0.0; 
  right_eigenmatrix[9][1] = (a*v3 - S13)*inva2;

	right_eigenmatrix[0][2] = 0.0; 
	right_eigenmatrix[1][2] = 0.0; 
	right_eigenmatrix[2][2] = 0.0; 
  right_eigenmatrix[3][2] = -inva;
	right_eigenmatrix[4][2] = 0.0; 
	right_eigenmatrix[5][2] = 0.0; 
  right_eigenmatrix[6][2] = 2.0*(S13-a*v3)*inva2;
	right_eigenmatrix[7][2] = 0.0; 
  right_eigenmatrix[8][2] = 1.0-v1*inva;
  right_eigenmatrix[9][2] = (S12-a*v2)*inva2;

	right_eigenmatrix[0][3] = 0.0; 
	right_eigenmatrix[1][3] = 0.0;
	right_eigenmatrix[2][3] = 0.0; 
	right_eigenmatrix[3][3] = 0.0; 
	right_eigenmatrix[4][3] = 0.0; 
  right_eigenmatrix[5][3] = 1.0;
	right_eigenmatrix[6][3] = 0.0; 
	right_eigenmatrix[7][3] = 0.0; 
	right_eigenmatrix[8][3] = 0.0; 
	right_eigenmatrix[9][3] = 0.0; 

	right_eigenmatrix[0][4] = 0.0; 
	right_eigenmatrix[1][4] = 0.0; 
	right_eigenmatrix[2][4] = 0.0; 
	right_eigenmatrix[3][4] = 0.0; 
	right_eigenmatrix[4][4] = 0.0; 
	right_eigenmatrix[5][4] = 0.0; 
  right_eigenmatrix[6][4] = 1.0;
	right_eigenmatrix[7][4] = 0.0; 
	right_eigenmatrix[8][4] = 0.0; 
	right_eigenmatrix[9][4] = 0.0; 

  right_eigenmatrix[0][5] = 1.0;
  right_eigenmatrix[1][5] = v1;
  right_eigenmatrix[2][5] = v2;
  right_eigenmatrix[3][5] = v3;
  right_eigenmatrix[4][5] = SQR(v1);
  right_eigenmatrix[5][5] = SQR(v2) + S22 - 4.0*SQR(S12*inva);
  right_eigenmatrix[6][5] = SQR(v3) + S33 - 4.0*SQR(S13*inva);
  right_eigenmatrix[7][5] = v1*v2;
  right_eigenmatrix[8][5] = v1*v3;
  right_eigenmatrix[9][5] = v2*v3 + S23 - 4.0*S12*S13*inva2;

	right_eigenmatrix[0][6] = 0.0; 
	right_eigenmatrix[1][6] = 0.0; 
	right_eigenmatrix[2][6] = 0.0; 
	right_eigenmatrix[3][6] = 0.0; 
	right_eigenmatrix[4][6] = 0.0; 
	right_eigenmatrix[5][6] = 0.0; 
	right_eigenmatrix[6][6] = 0.0; 
	right_eigenmatrix[7][6] = 0.0; 
	right_eigenmatrix[8][6] = 0.0; 
  right_eigenmatrix[9][6] = 1.0;

	right_eigenmatrix[0][7] = 0.0; 
	right_eigenmatrix[1][7] = 0.0; 
	right_eigenmatrix[2][7] = 0.0; 
  right_eigenmatrix[3][7] = inva;
	right_eigenmatrix[4][7] = 0.0; 
	right_eigenmatrix[5][7] = 0.0; 
  right_eigenmatrix[6][7] = 2.0*(S13 + a*v3)*inva2;
	right_eigenmatrix[7][7] = 0.0; 
	right_eigenmatrix[8][7] = 1.0+v1*inva;
  right_eigenmatrix[9][7] = (S12 + a*v2)*inva2;

	right_eigenmatrix[0][8] = 0.0; 
	right_eigenmatrix[1][8] = 0.0; 
  right_eigenmatrix[2][8] = inva;
	right_eigenmatrix[3][8] = 0.0; 
	right_eigenmatrix[4][8] = 0.0; 
  right_eigenmatrix[5][8] = 2.0*(a*v2 + S12)*inva2;
	right_eigenmatrix[6][8] = 0.0; 
  right_eigenmatrix[7][8] = 1.0+v1*inva;
	right_eigenmatrix[8][8] = 0.0; 
  right_eigenmatrix[9][8] = (a*v3 + S13)*inva2;

  right_eigenmatrix[0][9] = 1.0;
  right_eigenmatrix[1][9] = v1 + sqt3a;
  right_eigenmatrix[2][9] = v2 + sqt3*S12*inva;
  right_eigenmatrix[3][9] = v3 + sqt3*S13*inva;
  right_eigenmatrix[4][9] = SQR(v1 + sqt3a);
  right_eigenmatrix[5][9] = SQR(v2) + 2.0*sqt3*v2*S12*inva + S22 + 2.0*SQR(S12*inva);
  right_eigenmatrix[6][9] = SQR(v3) + 2.0*sqt3*v3*S13*inva + S33 + 2.0*SQR(S13*inva);
  right_eigenmatrix[7][9] = (v2+sqt3*S12*inva)*(v1+sqt3a);
  right_eigenmatrix[8][9] = (v3+sqt3*S13*inva)*(v1+sqt3a);
  right_eigenmatrix[9][9] = v2*v3 + 2.0*S12*S13*inva2 + S23 + sqt3*(S13*v2+S12*v3)*inva;

	// Left-eigenvectors, stored as ROWS 
  left_eigenmatrix[0][0] = v1*(v1 + sqt3a)*inva2/6.0;
  left_eigenmatrix[0][1] = -(2.0*v1 + sqt3a)*inva2/6.0;
	left_eigenmatrix[0][2] = 0.0; 
	left_eigenmatrix[0][3] = 0.0; 
  left_eigenmatrix[0][4] = inva2/6.0;
	left_eigenmatrix[0][5] = 0.0; 
	left_eigenmatrix[0][6] = 0.0; 
	left_eigenmatrix[0][7] = 0.0; 
	left_eigenmatrix[0][8] = 0.0; 
	left_eigenmatrix[0][9] = 0.0; 

  left_eigenmatrix[1][0] = -0.5*(v1 + a)*(a2*v2 - S12*v1)*inva2;
  left_eigenmatrix[1][1] = 0.5*(v2 - S12*(2*v1 + a)*inva2);
  left_eigenmatrix[1][2] = 0.5*(v1+a);
  left_eigenmatrix[1][3] = 0.0; 
  left_eigenmatrix[1][4] = 0.5*S12*inva2;
  left_eigenmatrix[1][5] = 0.0; 
  left_eigenmatrix[1][6] = 0.0; 
  left_eigenmatrix[1][7] = -0.5;
  left_eigenmatrix[1][8] = 0.0; 
  left_eigenmatrix[1][9] = 0.0; 

  left_eigenmatrix[2][0] = 0.5*(v1 + a)*(a2*v3 - S13*v1)*inva2;
  left_eigenmatrix[2][1] = -0.5*(v3 - S13*(2.0*v1 + a)*inva2);
  left_eigenmatrix[2][2] = 0.0; 
  left_eigenmatrix[2][3] = -0.5*(v1+a);
  left_eigenmatrix[2][4] = -0.5*S13*inva2;
  left_eigenmatrix[2][5] = 0.0; 
  left_eigenmatrix[2][6] = 0.0; 
  left_eigenmatrix[2][7] = 0.0; 
  left_eigenmatrix[2][8] = 0.5;
  left_eigenmatrix[2][9] = 0.0; 

  left_eigenmatrix[3][0] = SQR(v2) + 2.0*S12*(2.0*S12 - v1*v2)*inva2 - S22;
  left_eigenmatrix[3][1] = 2.0*S12*v2*inva2;
  left_eigenmatrix[3][2] = 2.0*(S12*v1*inva2 - v2);
  left_eigenmatrix[3][3] = 0.0; 
	left_eigenmatrix[3][4] = 0.0; 
  left_eigenmatrix[3][5] = 1.0;
	left_eigenmatrix[3][6] = 0.0;
  left_eigenmatrix[3][7] = -2.0*S12*inva2;
	left_eigenmatrix[3][8] = 0.0; 
	left_eigenmatrix[3][9] = 0.0; 

  left_eigenmatrix[4][0] = SQR(v3) + 2.0*S13*(2.0*S13 - v1*v3)*inva2 - S33;
  left_eigenmatrix[4][1] = 2.0*S13*v3*inva2;
	left_eigenmatrix[4][2] = 0.0; 
  left_eigenmatrix[4][3] = 2.0*(S13*v1*inva2 - v3);
	left_eigenmatrix[4][4] = 0.0; 
	left_eigenmatrix[4][5] = 0.0; 
  left_eigenmatrix[4][6] = 1.0;
	left_eigenmatrix[4][7] = 0.0; 
  left_eigenmatrix[4][8] = -2.0*S13*inva2;
	left_eigenmatrix[4][9] = 0.0;

  left_eigenmatrix[5][0] = 1.0 - SQR(v1)*inva2/3.0;
  left_eigenmatrix[5][1] = 2.0*v1*inva2/3.0;
	left_eigenmatrix[5][2] = 0.0; 
	left_eigenmatrix[5][3] = 0.0; 
  left_eigenmatrix[5][4] = -inva2/3.0;
	left_eigenmatrix[5][5] = 0.0; 
	left_eigenmatrix[5][6] = 0.0; 
	left_eigenmatrix[5][7] = 0.0; 
	left_eigenmatrix[5][8] = 0.0; 
	left_eigenmatrix[5][9] = 0.0; 

  left_eigenmatrix[6][0] = -(S13*v1*v2 + S12*(v1*v3 - 4.0*S13) + a2*(S23 - v2*v3))*inva2;
  left_eigenmatrix[6][1] = (S13*v2 + S12*v3)*inva2;
  left_eigenmatrix[6][2] = S13*v1*inva2 - v3;
  left_eigenmatrix[6][3] = S12*v1*inva2 - v2;
	left_eigenmatrix[6][4] = 0.0; 
	left_eigenmatrix[6][5] = 0.0; 
	left_eigenmatrix[6][6] = 0.0; 
  left_eigenmatrix[6][7] = -S13*inva2;
  left_eigenmatrix[6][8] = -S12*inva2;
  left_eigenmatrix[6][9] = 1.0;

  left_eigenmatrix[7][0] = -0.5*(a - v1)*(a2*v3 - S13*v1)*inva2;
  left_eigenmatrix[7][1] = -0.5*(v3 - S13*(2.0*v1 - a)*inva2);
	left_eigenmatrix[7][2] = 0.0; 
  left_eigenmatrix[7][3] = -0.5*(v1 - a);
  left_eigenmatrix[7][4] = -0.5*S13*inva2;
	left_eigenmatrix[7][5] = 0.0; 
	left_eigenmatrix[7][6] = 0.0; 
	left_eigenmatrix[7][7] = 0.0; 
  left_eigenmatrix[7][8] = 0.5;
	left_eigenmatrix[7][9] = 0.0; 

  left_eigenmatrix[8][0] = -0.5*(a - v1)*(a2*v2 - S12*v1)*inva2;
  left_eigenmatrix[8][1] = -0.5*(v2 - S12*(2.0*v1 - a)*inva2);
  left_eigenmatrix[8][2] = -0.5*(v1 - a);
	left_eigenmatrix[8][3] = 0.0; 
  left_eigenmatrix[8][4] = -0.5*S12*inva2;
	left_eigenmatrix[8][5] = 0.0; 
	left_eigenmatrix[8][6] = 0.0; 
  left_eigenmatrix[8][7] = 0.5;
	left_eigenmatrix[8][8] = 0.0; 
	left_eigenmatrix[8][9] = 0.0; 

  left_eigenmatrix[9][0] = v1*(v1 - sqt3a)*inva2/6.0;
  left_eigenmatrix[9][1] = -(2.0*v1 - sqt3a)*inva2/6.0;
	left_eigenmatrix[9][2] = 0.0; 
	left_eigenmatrix[9][3] = 0.0; 
  left_eigenmatrix[9][4] = inva2/6.0;
	left_eigenmatrix[9][5] = 0.0; 
	left_eigenmatrix[9][6] = 0.0; 
	left_eigenmatrix[9][7] = 0.0; 
	left_eigenmatrix[9][8] = 0.0; 
	left_eigenmatrix[9][9] = 0.0; 
}
