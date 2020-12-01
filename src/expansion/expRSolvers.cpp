//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ppmEX.cpp
//  \brief piecewise parabolic reconstruction with modified McCorquodale/Colella limiter
//         for a uniform Cartesian mesh, Mignone limiter for nonuniform mesh, Accounting
//         for expanding grid
//
// REFERENCES:
// (CW) P. Colella & P. Woodward, "The Piecewise Parabolic Method (PPM) for Gas-Dynamical
// Simulations", JCP, 54, 174 (1984)
//
// (CS) P. Colella & M. Sekora, "A limiter for PPM that preserves accuracy at smooth
// extrema", JCP, 227, 7069 (2008)
//
// (MC) P. McCorquodale & P. Colella,  "A high-order finite-volume method for conservation
// laws on locally refined grids", CAMCoS, 6, 1 (2011)
//
// (CD) P. Colella, M.R. Dorr, J. Hittinger, D. Martin, "High-order, finite-volume methods
// in mapped coordinates", JCP, 230, 2952 (2011)
//
// (Mignone) A. Mignone, "High-order conservative reconstruction schemes for finite volume
// methods in cylindrical and spherical coordinates", JCP, 270, 784 (2014)
//========================================================================================

// C++ headers
#include <algorithm>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "expansion.hpp"


void Expansion::FluxSolver(const int kl, const int ku, const int jl, const int ju,
    const int il, const int iu, const int ivx, const AthenaArray<Real> &bx,
    AthenaArray<Real> &wL, AthenaArray<Real> &wR, 
    AthenaArray<Real> &flx,
    AthenaArray<Real> &e1, AthenaArray<Real> &e2, AthenaArray<Real> &vArr){

  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  Real wi[(NHYDRO)];
  Real flxi[(NHYDRO)],f[(NHYDRO)];
  Real gm1 = pmy_block->peos->GetGamma() - 1.0;
  Real igm1 = 1.0/gm1;
  Real wallV;

  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
#pragma distribute_point
#pragma omp simd private(wli,wri,wroe,flxi,fl,fr)
  for (int i=il; i<=iu; ++i) {
 
//--- Step 1.  Load L/R states into local variables

    wallV = vArr(i);
    Real midV = 0.5*(wR(ivx,k,j,i) + wL(ivx,k,j,i));
    if (wallV > 0.0) {    
      wi[IDN]=wR(IDN,k,j,i);
      wi[IVX]=wR(ivx,k,j,i);
      wi[IVY]=wR(ivy,k,j,i);
      wi[IVZ]=wR(ivz,k,j,i);
      wi[IPR]=wR(IPR,k,j,i);
      if (DUAL_ENERGY) wi[IGE]=wR(IGE,k,j,i);
    } else if (wallV < 0.0) {    
      wi[IDN]=wL(IDN,k,j,i);
      wi[IVX]=wL(ivx,k,j,i);
      wi[IVY]=wL(ivy,k,j,i);
      wi[IVZ]=wL(ivz,k,j,i);
      wi[IPR]=wL(IPR,k,j,i);
      if (DUAL_ENERGY) wi[IGE]=wL(IGE,k,j,i);
    } else if (wallV == 0){
      wi[IDN]=0.0;
      wi[IVX]=0.0;
      wi[IVY]=0.0;
      wi[IVZ]=0.0;
      wi[IPR]=0.0;
      if (DUAL_ENERGY) wi[IGE]=0.0;
    }

    Real e; 
    e = wi[IPR]*igm1 + 0.5*wi[IDN]*(SQR(wi[IVX]) + SQR(wi[IVY]) + SQR(wi[IVZ]));


//--- Step 2.  Compute L/R fluxes for each state

    f[IDN] = wi[IDN]*wallV;

    f[IVX] = wi[IDN]*wi[IVX]*wallV;

    f[IVY] = wi[IDN]*wi[IVY]*wallV;

    f[IVZ] = wi[IDN]*wi[IVZ]*wallV;

    f[IEN] = e*wallV;

//--- Step 3.  Compute the  flux at interface
    flxi[IDN] = f[IDN];
    flxi[IVX] = f[IVX];
    flxi[IVY] = f[IVY];
    flxi[IVZ] = f[IVZ];
    flxi[IEN] = f[IEN];
   
    flx(IDN,k,j,i) = flxi[IDN];
    flx(ivx,k,j,i) = flxi[IVX];
    flx(ivy,k,j,i) = flxi[IVY];
    flx(ivz,k,j,i) = flxi[IVZ];
    flx(IEN,k,j,i) = flxi[IEN];
    if (DUAL_ENERGY) {
      if (flxi[IDN]  >= 0) {
        flx(IIE,k,j,i) = flxi[IDN]*wi[IGE];
      }
      else {
        flx(IIE,k,j,i) = flxi[IDN]*wi[IGE];
      }
    }
  }
  }}
 
  return;
}

