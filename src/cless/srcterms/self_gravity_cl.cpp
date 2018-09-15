//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief source terms due to self-gravity

// Athena++ headers
#include "cless_srcterms.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../cless.hpp"
#include "../../gravity/gravity.hpp"

//----------------------------------------------------------------------------------------
//! \fn void ClessSourceTerms::SelfGravity
//  \brief Adds source terms for self-gravitational acceleration to conserved variables
//				 Momentum terms are added directly to the flux (cf. gravity_fluxes_cl.cpp)
//				 Energy source terms are added in this file 

void ClessSourceTerms::SelfGravityCL(const Real dt,const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_cless_->pmy_block;
  if (SELF_GRAVITY_ENABLED) {
    Gravity *pgrav = pmb->pgrav;
    Real four_pi_G = pgrav->four_pi_G;
    Real grav_mean_rho = pgrav->grav_mean_rho;
    Real phic, phir, phil;
// acceleration in 1-direction
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real dx1 = pmb->pcoord->dx1v(i);
          Real dx2 = pmb->pcoord->dx2v(j);
          Real dx3 = pmb->pcoord->dx3v(k);
          Real dtodx1 = dt/dx1;
          phic = pgrav->phi(k,j,i);
          phil = 0.5*(pgrav->phi(k,j,i-1)+pgrav->phi(k,j,i  ));
          phir = 0.5*(pgrav->phi(k,j,i  )+pgrav->phi(k,j,i+1));
          // Update Eij with d/dx1 terms
					// For CLESS use the primitive variables to update the conserved variables
          cons(IE11,k,j,i) -= 2.0*dtodx1*(phir-phil)*prim(IDN,k,j,i)*prim(IVX,k,j,i);
					cons(IE12,k,j,i) -=     dtodx1*(phir-phil)*prim(IDN,k,j,i)*prim(IVY,k,j,i);
					cons(IE13,k,j,i) -=			dtodx1*(phir-phil)*prim(IDN,k,j,i)*prim(IVZ,k,j,i);
        }
      }
    }

    if (pmb->block_size.nx2 > 1) {
      // acceleration in 2-direction
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real dx1 = pmb->pcoord->dx1v(i);
            Real dx2 = pmb->pcoord->dx2v(j);
            Real dx3 = pmb->pcoord->dx3v(k);
            Real dtodx2 = dt/dx2;
            phic = pgrav->phi(k,j,i);
            phil = 0.5*(pgrav->phi(k,j-1,i)+pgrav->phi(k,j  ,i));
            phir = 0.5*(pgrav->phi(k,j  ,i)+pgrav->phi(k,j+1,i));
						// Update Eij with d/dx2 terms
						// For CLESS use the primitive variables to update the conserved variables
						cons(IE22,k,j,i) -= 2.0*dtodx2*(phir-phil)*prim(IDN,k,j,i)*prim(IVY,k,j,i);
						cons(IE12,k,j,i) -=     dtodx2*(phir-phil)*prim(IDN,k,j,i)*prim(IVX,k,j,i);
						cons(IE23,k,j,i) -=     dtodx2*(phir-phil)*prim(IDN,k,j,i)*prim(IVZ,k,j,i);
          }
        }
      }
    }

    if (pmb->block_size.nx3 > 1) {
      // acceleration in 3-direction
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real dx1 = pmb->pcoord->dx1v(i);
            Real dx2 = pmb->pcoord->dx2v(j);
            Real dx3 = pmb->pcoord->dx3v(k);
            Real dtodx3 = dt/dx3;
            phic = pgrav->phi(k,j,i);
            phil = 0.5*(pgrav->phi(k-1,j,i)+pgrav->phi(k  ,j,i));
            phir = 0.5*(pgrav->phi(k  ,j,i)+pgrav->phi(k+1,j,i));
						// Update Eij with d/dx3 terms
						// For CLESS use the primitive variables to update the conserved variables
						cons(IE33,k,j,i) -= 2.0*dtodx3*(phir-phil)*prim(IDN,k,j,i)*prim(IVZ,k,j,i);
						cons(IE13,k,j,i) -=     dtodx3*(phir-phil)*prim(IDN,k,j,i)*prim(IVX,k,j,i);
						cons(IE23,k,j,i) -=     dtodx3*(phir-phil)*prim(IDN,k,j,i)*prim(IVY,k,j,i);
          }
        }
      }
    }
  }
  return;
}

