//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file calculate_fluxes.cpp
//  \brief Calculate hydro/MHD fluxes

// C/C++ headers
#include <algorithm>   // min,max

// Athena++ headers
#include "cless.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../bvals/bvals.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "../gravity/gravity.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CalculateFluxes
//  \brief Calculate Hydrodynamic Fluxes using the Riemann solver

void Cless::CalculateFluxes(AthenaArray<Real> &w, int order) {
  MeshBlock *pmb=pmy_block;
  AthenaArray<Real> &x1flux=flux[X1DIR];
  AthenaArray<Real> &x2flux=flux[X2DIR];
  AthenaArray<Real> &x3flux=flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;

  AthenaArray<Real> w_x1f,w_x2f,w_x3f;

  AthenaArray<Real> wl, wr, dxw;
  wl.InitWithShallowCopy(wl_);
  wr.InitWithShallowCopy(wr_);
  dxw.InitWithShallowCopy(dxw_);

//----------------------------------------------------------------------------------------
// i-direction

  // set the loop limits
  jl=js, ju=je, kl=ks, ku=ke;

  // reconstruct L/R states
  if (order == 1) {
    pmb->precon->DonorCellCLX1(pmb,kl,ku,jl,ju,is,ie+1,w,wl,wr);
  } else if (order == 2) {
    pmb->precon->PiecewiseLinearCLX1(pmb,kl,ku,jl,ju,is,ie+1,w,wl,wr);
  } else {
    pmb->precon->PiecewiseParabolicCLX1(pmb,kl,ku,jl,ju,is,ie+1,w,wl,wr);
  }

  // compute fluxes, store directly into 3D arrays
  RiemannSolver(kl,ku,jl,ju,is,ie+1,IVX,IP12,wl,wr,x1flux);


//----------------------------------------------------------------------------------------
// j-direction

  if (pmb->block_size.nx2 > 1) {

    // set the loop limits
    il=is, iu=ie, kl=ks, ku=ke;
    
		// reconstruct L/R states at j
    if (order == 1) {
      pmb->precon->DonorCellCLX2(pmb,kl,ku,js,je+1,il,iu,w,bcc,wl,wr);
    } else if (order == 2) {
      pmb->precon->PiecewiseLinearCLX2(pmb,kl,ku,js,je+1,il,iu,w,wl,wr);
    } else {
      pmb->precon->PiecewiseParabolicCLX2(pmb,kl,ku,js,je+1,il,iu,w,wl,wr);
    }

    // compute fluxes, store directly into 3D arrays
    RiemannSolver(kl,ku,js,je+1,il,iu,IVY,IP23,wl,wr,x2flux);

    // compute weights for GS07 CT algorithm
  }

//----------------------------------------------------------------------------------------
// k-direction

  if (pmb->block_size.nx3 > 1) {

    // set the loop limits
    il=is, iu=ie, jl=js, ju=je;

    // reconstruct L/R states at k
    if (order == 1) {
      pmb->precon->DonorCellCLX3(pmb,ks,ke+1,jl,ju,il,iu,w,wl,wr);
    } else if (order == 2) {
      pmb->precon->PiecewiseLinearCLX3(pmb,ks,ke+1,jl,ju,il,iu,w,wl,wr);
    } else {
      pmb->precon->PiecewiseParabolicCLX3(pmb,ks,ke+1,jl,ju,il,iu,w,wl,wr);
    }

    // compute fluxes, store directly into 3D arrays
    RiemannSolver(ks,ke+1,jl,ju,il,iu,IVZ,IP13,wl,wr,x3flux);

    // compute weights for GS07 CT algorithm
  }

  if (SELF_GRAVITY_ENABLED) AddGravityFlux(); // add gravity flux directly
  
	return;
}
