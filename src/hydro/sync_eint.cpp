//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file sync_eint.cpp
//  \brief Syncs internal with total energy

// Athena++ headers
#include "hydro.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../eos/eos.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::SyncEint
//  \brief Sync internal energy with total energy, if necessary 

void Hydro::SyncEint(AthenaArray<Real> &u) {
  MeshBlock *pmb=pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
	Real i2 = pmb->peos->GetIeta2(); 
	Real eint, emax;  
	int dim = 1; 

	// Get dimension
	if (pmb->block_size.nx2 > 1) dim = 2;
	if (pmb->block_size.nx3 > 1) dim = 3;  

	// 1D-case
	if ( dim == 1 ) {
#pragma omp simd
		for (int i=is; i<=ie; ++i) {
			eint = u(IEN,ks,js,i) - 0.5*(SQR(u(IVX,ks,js,i)))/u(IDN,ks,js,i);
			emax = u(IEN,ks,js,i); 
			for (int ii=std::max(is,i-1); ii<=std::min(ie,i+1); ii++) {
				emax = std::max(emax,u(IEN,ks,js,ii));
			}
			if ((eint/emax > i2) && (eint > 0.0)) {
				u(IIE,ks,js,i) = eint; 
			}
		}
	}
	// 2D-case
	else if ( dim == 2 ) {
		for (int j=js; j<=je; ++j) {
#pragma omp simd
			for (int i=is; i<=ie; ++i) {
				eint = u(IEN,ks,j,i) - 0.5*(SQR(u(IVX,ks,j,i))
															+     SQR(u(IVY,ks,j,i)))/u(IDN,ks,j,i);
				emax = u(IEN,ks,j,i); 
				for (int jj=std::max(js,j-1); jj<=std::min(je,j+1); jj++) {
					for (int ii=std::max(is,i-1); ii<=std::min(ie,i+1); ii++) {
						emax = std::max(emax,u(IEN,ks,jj,ii));
					}
				}
				if ((eint/emax > i2) && (eint > 0.0)) {
					//std::cout << "[SyncEint]: k " << ks << " j " << j  << " i "  << i << " eint=" << eint
					//					<< " IE: " << u(IIE,ks,j,i) << " ieta2: " << i2 << std::endl;
					u(IIE,ks,j,i) = eint; 

				}
			}
		}
	}
	// 3D-case
	else {
		for (int k=ks; k<=ke; ++k) {
			for (int j=js; j<=je; ++j) {
#pragma omp simd
				for (int i=is; i<=ie; ++i) {
					eint = u(IEN,k,j,i) - 0.5*(SQR(u(IVX,k,j,i))
																+    SQR(u(IVY,k,j,i)) 
																+    SQR(u(IVZ,k,j,i)))/u(IDN,k,j,i);
					emax = u(IEN,k,j,i); 
					//std::cout << "[SyncEint]: k" << k << " j " << j  << " i "  << i << " eint=" << eint
					//					<< " IE: " << u(IIE,k,j,i) << std::endl;
					for (int kk=std::max(ks,k-1); kk<=std::min(ke,k+1); kk++) {
						for (int jj=std::max(js,j-1); jj<=std::min(je,j+1); jj++) {
							for (int ii=std::max(is,i-1); ii<=std::min(ie,i+1); ii++) {
								emax = std::max(emax,u(IEN,kk,jj,ii));
							}
						}
					}
					if ((eint/emax > i2) && (eint > 0.0)) {
						u(IIE,k,j,i) = eint; 
					}
				}
			}
		}
	}
	return;
}

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CheckEint
//  \brief Check if internal energy is necessary at start of each substep  

void Hydro::CheckEint(AthenaArray<Real> &u) {
  MeshBlock *pmb=pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
	Real i1 = pmb->peos->GetIeta1(); 
	Real eint, etot, kine;  
	int dim = 1; 

	// Get dimension
	if (pmb->block_size.nx2 > 1) dim = 2;
	if (pmb->block_size.nx3 > 1) dim = 3;  

	// 1D-case
	if ( dim == 1 ) {
#pragma omp simd
		for (int i=is-NGHOST; i<=ie+NGHOST; ++i) {
			etot = u(IEN,ks,js,i);
			kine = 0.5*(SQR(u(IVX,ks,js,i)))/u(IDN,ks,js,i);
			eint = etot - kine;
			if ( ((eint/etot) < i1) || (etot <= 0.0) || (eint <= 0.0) ) {
				u(IEN,ks,js,i)  = u(IIE,ks,js,i);
				u(IEN,ks,js,i) += kine; 
			}
		}
	}
	// 2D-case
	else if ( dim == 2 ) {
		for (int j=js-NGHOST; j<=je+NGHOST; ++j) {
#pragma omp simd
			for (int i=is-NGHOST; i<=ie+NGHOST; ++i) {
				etot = u(IEN,ks,j,i);
				kine = 0.5*( SQR(u(IVX,ks,j,i))
									 + SQR(u(IVY,ks,j,i)) )/u(IDN,ks,j,i);
				eint = etot - kine;
				if ( ((eint/etot) < i1) || (etot <= 0.0) || (eint <= 0.0) ) {
					u(IEN,ks,j,i)  = u(IIE,ks,j,i);
					u(IEN,ks,j,i) += kine; 
				}
			}
	}	}
	// 3D-case
	else {
		for (int k=ks-NGHOST; k<=ke+NGHOST; ++k) {
			for (int j=js-NGHOST; j<=je+NGHOST; ++j) {
#pragma omp simd
				for (int i=is-NGHOST; i<=ie+NGHOST; ++i) {
					etot = u(IEN,k,j,i); 
					kine = 0.5*( SQR(u(IVX,k,j,i))
								 	   + SQR(u(IVY,k,j,i)) 
										 + SQR(u(IVZ,k,j,i)) )/u(IDN,k,j,i);
					eint = etot - kine; 
					if ( ((eint/etot) < i1) || (etot <= 0.0) || (eint <= 0.0) ) {
						u(IEN,k,j,i)  = u(IIE,k,j,i);
						u(IEN,k,j,i) += kine;
					}
				}
			}
		}
	}
	return;
}
