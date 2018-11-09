//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief source terms due to self-gravity

// Athena++ headers
#include "hydro_srcterms.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
//#include "../../gravity/gravity.hpp"

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::StaticGravity
//  \brief Adds source terms for static gravitational accelerations to energy

void HydroSourceTerms::StaticGravity(const Real dt,const Real time, const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  Real phic, phir, phil;
  // acceleration in 1-direction
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    Real x3v = pmb->pcoord->x3v(k);
    for (int j=pmb->js; j<=pmb->je; ++j) {
      Real x2v = pmb->pcoord->x2v(j);
#pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real x1v         = pmb->pcoord->x1v(i);
        Real dx1v        = pmb->pcoord->dx1v(i);
        Real dtodx1      = dt/dx1v;
        phic             = StaticGravPot(x1v         ,x2v,x3v,time);
        phir             = StaticGravPot(x1v+0.5*dx1v,x2v,x3v,time);
        phil             = StaticGravPot(x1v-0.5*dx1v,x2v,x3v,time);
        cons(IM1,k,j,i) -= dtodx1*cons(IDN,k,j,i)*(phir-phil);
        if (NON_BAROTROPIC_EOS)
          cons(IEN,k,j,i) -= dtodx1*( flux[X1DIR](IDN,k,j,i  )*(phic-phil)
                                     +flux[X1DIR](IDN,k,j,i+1)*(phir-phic));
      }
    }
  }

  if (pmb->block_size.nx2 > 1) {
    // acceleration in 2-direction
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      Real x3v = pmb->pcoord->x3v(k);
      for (int j=pmb->js; j<=pmb->je; ++j) {
        Real x2v    = pmb->pcoord->x2v(j);
        Real dx2v   = pmb->pcoord->dx2v(j);
        Real dtodx2 = dt/dx2v;
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real x1v         = pmb->pcoord->x1v(i);
          phic             = StaticGravPot(x1v,x2v         ,x3v,time);
          phir             = StaticGravPot(x1v,x2v+0.5*dx2v,x3v,time);
          phil             = StaticGravPot(x1v,x2v-0.5*dx2v,x3v,time);
          cons(IM2,k,j,i) -= dtodx2*cons(IDN,k,j,i)*(phir-phil);
          if (NON_BAROTROPIC_EOS)
            cons(IEN,k,j,i) -= dtodx2*( flux[X2DIR](IDN,k,j  ,i)*(phic - phil)
                                       +flux[X2DIR](IDN,k,j+1,i)*(phir - phic));
        }
      }
    }
  }

  if (pmb->block_size.nx3 > 1) {
    // acceleration in 3-direction
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      Real x3v    = pmb->pcoord->x3v(k);
      Real dx3v   = pmb->pcoord->dx3v(k);
      Real dtodx3 = dt/dx3v;
      for (int j=pmb->js; j<=pmb->je; ++j) {
        Real x2v = pmb->pcoord->x2v(j);
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real x1v         = pmb->pcoord->x1v(i);
          phic             = StaticGravPot(x1v,x2v,x3v         ,time);
          phir             = StaticGravPot(x1v,x2v,x3v+0.5*dx3v,time);     
          phil             = StaticGravPot(x1v,x2v,x3v-0.5*dx3v,time);
          cons(IM3,k,j,i) -= dtodx3*cons(IDN,k,j,i)*(phir-phil);
          if (NON_BAROTROPIC_EOS)
            cons(IEN,k,j,i) -= dtodx3*( flux[X3DIR](IDN,k  ,j,i)*(phic - phil) 
                                       +flux[X3DIR](IDN,k+1,j,i)*(phir - phic));
        }
      }
    }
  }
  return;
}

