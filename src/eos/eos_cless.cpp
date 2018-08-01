
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file eos_cless.cpp
//  \brief implements functions in class EquationOfState for collisionless system 

// C/C++ headers
#include <cmath>   // sqrt()
#include <cfloat>  // FLT_MIN

// Athena++ headers
#include "eos.hpp"
#include "../hydro/hydro.hpp"
#include "../cless/cless.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::ConsclToPrimcl(AthenaArray<Real> &cons,
//           const AthenaArray<Real> &prim_old,
//           AthenaArray<Real> &prim, Coordinates *pco,
//           int il, int iu, int jl, int ju, int kl, int ku)
// \brief Converts conserved into primitive variables for CLESS vars  

void EquationOfState::ConsclToPrimcl(AthenaArray<Real> &cons,
  const AthenaArray<Real> &prim_old, AthenaArray<Real> &prim,
  Coordinates *pco, int il,int iu, int jl,int ju, int kl,int ku) {

  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real& u_d   = cons(IDN ,k,j,i);
      Real& u_m1  = cons(IM1 ,k,j,i);
      Real& u_m2  = cons(IM2 ,k,j,i);
      Real& u_m3  = cons(IM3 ,k,j,i);
      Real& u_e11 = cons(IE11,k,j,i);
			Real& u_e22 = cons(IE22,k,j,i);
			Real& u_e33 = cons(IE33,k,j,i);
			Real& u_e12 = cons(IE12,k,j,i);
			Real& u_e13 = cons(IE13,k,j,i);
			Real& u_e23 = cons(IE23,k,j,i);

      Real& w_d   = prim(IDN ,k,j,i);
      Real& w_vx  = prim(IVX ,k,j,i);
      Real& w_vy  = prim(IVY ,k,j,i);
      Real& w_vz  = prim(IVZ ,k,j,i);
      Real& w_p11 = prim(IP11,k,j,i);
			Real& w_p22 = prim(IP22,k,j,i);
			Real& w_p33 = prim(IP33,k,j,i);
      Real& w_p12 = prim(IP12,k,j,i);
			Real& w_p13 = prim(IP13,k,j,i);
			Real& w_p23 = prim(IP23,k,j,i);
			

      // apply density floor, without changing momentum or energy
      u_d = (u_d > density_floor_) ?  u_d : density_floor_;
      w_d = u_d;

      Real di = 1.0/u_d;
      w_vx = u_m1*di;
      w_vy = u_m2*di;
      w_vz  = u_m3*di;

      Real ke11 = di*( u_m1*u_m1 );
			Real ke22 = di*( u_m2*u_m2 );
			Real ke33 = di*( u_m3*u_m3 );
      Real ke12 = di*( u_m1*u_m2 );
			Real ke13 = di*( u_m1*u_m3 );
			Real ke23 = di*( u_m2*u_m3 );

			w_p11 = u_e11 - ke11; 
			w_p22 = u_e22 - ke22; 
			w_p33 = u_e33 - ke33; 
			w_p12 = u_e12 - ke12; 
			w_p13 = u_e13 - ke13; 
			w_p23 = u_e23 - ke23; 
			
			// Apply pressure floors to diagnoal components of 'presure' tensor
			u_e11 = (w_p11 > pressure_floor_) ?  u_e11 : ((pressure_floor_) + ke11);
			u_e22 = (w_p22 > pressure_floor_) ?  u_e22 : ((pressure_floor_) + ke22);
			u_e33 = (w_p33 > pressure_floor_) ?  u_e33 : ((pressure_floor_) + ke33);

			w_p11 = (w_p11 > pressure_floor_) ?  w_p11 : pressure_floor_;
			w_p22 = (w_p22 > pressure_floor_) ?  w_p22 : pressure_floor_;
			w_p33 = (w_p33 > pressure_floor_) ?  w_p33 : pressure_floor_;
    }
	}}

  return;
}


//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::PrimclToConscl(const AthenaArray<Real> &prim,
//           AthenaArray<Real> &cons, Coordinates *pco,
//           int il, int iu, int jl, int ju, int kl, int ku);
// \brief Converts primitive variables into conservative variables

void EquationOfState::PrimclToConscl(const AthenaArray<Real> &prim,
     AthenaArray<Real> &cons, Coordinates *pco,
     int il, int iu, int jl, int ju, int kl, int ku) {

  // Force outer-loop vectorization
#pragma omp simd
  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
    //#pragma omp simd
#pragma novector
    for (int i=il; i<=iu; ++i) {
      Real& u_d   = cons(IDN ,k,j,i);
      Real& u_m1  = cons(IM1 ,k,j,i);
      Real& u_m2  = cons(IM2 ,k,j,i);
      Real& u_m3  = cons(IM3 ,k,j,i);
      Real& u_e11 = cons(IE11,k,j,i);
			Real& u_e22 = cons(IE22,k,j,i);
			Real& u_e33 = cons(IE33,k,j,i);
			Real& u_e12 = cons(IE12,k,j,i);
			Real& u_e13 = cons(IE13,k,j,i);
			Real& u_e23 = cons(IE23,k,j,i);

      const Real& w_d   = prim(IDN ,k,j,i);
      const Real& w_vx  = prim(IVX ,k,j,i);
      const Real& w_vy  = prim(IVY ,k,j,i);
      const Real& w_vz  = prim(IVZ ,k,j,i);
      const Real& w_p11 = prim(IP11,k,j,i);
			const Real& w_p22 = prim(IP22,k,j,i);
			const Real& w_p33 = prim(IP33,k,j,i);
      const Real& w_p12 = prim(IP12,k,j,i);
			const Real& w_p13 = prim(IP13,k,j,i);
			const Real& w_p23 = prim(IP23,k,j,i);

      u_d   = w_d;
      u_m1  = w_vx*w_d;
      u_m2  = w_vy*w_d;
      u_m3  = w_vz*w_d;
      u_e11 = w_p11 + w_d*w_vx*w_vx;
			u_e22 = w_p22 + w_d*w_vy*w_vy;
      u_e33 = w_p11 + w_d*w_vz*w_vz;
			u_e12 = w_p12 + w_d*w_vx*w_vy;
      u_e13 = w_p13 + w_d*w_vx*w_vz;
			u_e23 = w_p23 + w_d*w_vy*w_vz;
    }
  }}
  return;
}


//----------------------------------------------------------------------------------------
// \!fn Real EquationOfState::SoundSpeedsCL(Real prim[NCLESS])
// \brief returns sound speed in each direction for CLESS variables
void EquationOfState::SoundSpeedsCL(const Real prim[NCLESS], 
																		Real *c11, Real *c22, Real *c33) {
	*c11 = std::sqrt(3*prim[IP11]/prim[IDN]);
	*c22 = std::sqrt(3*prim[IP22]/prim[IDN]);
	*c33 = std::sqrt(3*prim[IP33]/prim[IDN]);
	return;
}

//---------------------------------------------------------------------------------------
// \!fn void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim,
//           int k, int j, int i)
// \brief Apply density and pressure floors to reconstructed L/R cell interface states

void EquationOfState::ApplyPrimitiveFloorsCL(AthenaArray<Real> &prim, int k, int j, int i) {
  Real& w_d   = prim(IDN ,k,j,i);
  Real& w_p11 = prim(IP11,k,j,i);
	Real& w_p22 = prim(IP22,k,j,i);
	Real& w_p33 = prim(IP33,k,j,i);

  // apply density floor
  w_d = (w_d > density_floor_) ?  w_d : density_floor_;
  // apply pressure floor
  w_p11 = (w_p11 > pressure_floor_) ?  w_p11 : pressure_floor_;
  w_p22 = (w_p22 > pressure_floor_) ?  w_p22 : pressure_floor_;
  w_p33 = (w_p33 > pressure_floor_) ?  w_p33 : pressure_floor_;

  return;
}
