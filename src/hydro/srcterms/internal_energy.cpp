//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief Adds source terms for internal energy, when using dual-energy formalism 

// Athena++ headers
#include "hydro_srcterms.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
#include "../../eos/eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::InternalEnergy
//  \brief Adds source terms for internal energy, when using dual-energy formalism 
//				 de/dt = -P div(v) 
//				 The calls to VolCenterLength are necessary for taking the gradient
//				 in arbitrary coordinate systems. 

void HydroSourceTerms::InternalEnergy(const Real dt, const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;
	if (DUAL_ENERGY && NON_BAROTROPIC_EOS) {
		Real Pc, vr, vl; 
		Real gm1 = pmb->peos->GetGamma() - 1.0;
		AthenaArray<Real> dx1, dx2m1, dx2p1, dx3m1, dx3p1; 
		
		dx1.InitWithShallowCopy(pmb->pcoord->dx1v); 
		if (pmb->block_size.nx2 > 1) {
			dx2m1.InitWithShallowCopy(pmb->pcoord->dx1v);
			dx2p1.InitWithShallowCopy(pmb->pcoord->dx1v);
		}
		if (pmb->block_size.nx3 > 1) {
			dx3m1.InitWithShallowCopy(pmb->pcoord->dx1v);
			dx3p1.InitWithShallowCopy(pmb->pcoord->dx1v); 
		}
		
		for (int k=pmb->ks; k<=pmb->ke; ++k) {
			for (int j=pmb->js; j<=pmb->je; ++j) {

				// calculate x1-grad
				pmb->pcoord->VolCenter1Length(k, j, pmb->is-1, pmb->ie+1, dx1); 
#pragma omp simd
				for (int i=pmb->is; i<=pmb->ie; ++i) {
					Real dtodx1 = dt/(dx1(i-1) + dx1(i)); 
					// cell-centered pressure
					Pc	= cons(IIE,k,j,i)*gm1;
					// vr, vl
					vr  = prim(IVX,k,j,i+1);
					vl  = prim(IVX,k,j,i-1); 
					// for the time being, this assumes uniform spacing  
					cons(IIE,k,j,i) -= dtodx1*Pc*(vr-vl); 
				}

				if (pmb->block_size.nx2 > 1) {
					// calculate x2-grad
					pmb->pcoord->VolCenter2Length(k, j-1, pmb->is, pmb->ie, dx2m1);
					pmb->pcoord->VolCenter2Length(k, j  , pmb->is, pmb->ie, dx2p1); 
#pragma omp simd
					for (int i=pmb->is; i<=pmb->ie; ++i) {
						Real dtodx2 = dt/((dx2p1(i) + dx2m1(i))); 
						// cell-centered pressure
						Pc	= cons(IIE,k,j,i)*gm1;
						// vr, vl
						vr  = prim(IVY,k,j+1,i);
						vl  = prim(IVY,k,j-1,i); 
						// for the time being, this assumes uniform spacing
						cons(IIE,k,j,i) -= dtodx2*Pc*(vr-vl); 
					}
				}
				if (pmb->block_size.nx3 > 1) {
					// calculate x3-grad
					pmb->pcoord->VolCenter3Length(k-1, j, pmb->is, pmb->ie, dx3m1);
					pmb->pcoord->VolCenter3Length(k  , j, pmb->is, pmb->ie, dx3p1);
#pragma omp simd
					for (int i=pmb->is; i<=pmb->ie; ++i) {
						Real dtodx3 = dt/(dx3p1(i) + dx3m1(i)); 
						// cell-centered pressure
						Pc	= cons(IIE,k,j,i)*gm1;
						// vr, vl
						vr  = prim(IVZ,k+1,j,i);
						vl  = prim(IVZ,k-1,j,i); 
						// for the time being, this assumes uniform spacing
						cons(IIE,k,j,i) -= dtodx3*Pc*(vr-vl); 
					}
				}
			}
		}
	}

  return;
}
