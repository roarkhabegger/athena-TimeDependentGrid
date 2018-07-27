//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file dc_cl.cpp
//  \brief piecewise constant (donor cell) reconstruction for cless vars

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellCLX1()
//  \brief

void Reconstruction::DonorCellCLX1(MeshBlock *pmb,
  const int kl, const int ku, const int jl, const int ju, const int il, const int iu,
  const AthenaArray<Real> &w,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  // compute L/R states for each variable
  for (int n=0; n<(NCLESS); ++n) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        wl(n,k,j,i) = w(n,k,j,i-1);
        wr(n,k,j,i) = w(n,k,j,i  );
      }
    }}
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellCLX2()
//  \brief

void Reconstruction::DonorCellCLX2(MeshBlock *pmb,
  const int kl, const int ku, const int jl, const int ju, const int il, const int iu,
  const AthenaArray<Real> &w,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  // compute L/R states for each variable
  for (int n=0; n<(NCLESS); ++n) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        wl(n,k,j,i) = w(n,k,j-1,i);
        wr(n,k,j,i) = w(n,k,j  ,i);
      }
    }}
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellCLX3()
//  \brief

void Reconstruction::DonorCellCLX3(MeshBlock *pmb,
  const int kl, const int ku, const int jl, const int ju, const int il, const int iu,
  const AthenaArray<Real> &w,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  // compute L/R states for each variable
  for (int n=0; n<(NCLESS); ++n) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        wl(n,k,j,i) = w(n,k-1,j,i);
        wr(n,k,j,i) = w(n,k  ,j,i);
      }
    }}
  }
  return;
}
