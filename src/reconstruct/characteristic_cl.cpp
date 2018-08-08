//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file characteristic.cpp
//  \brief Functions to transform vectors between primitive and characteristic variables

// C++ headers
#include <cmath>

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../eos/eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::LeftEigenmatrixDotVectorCL()
//  \brief Computes inner-product of left-eigenmatrix of Roe's matrix A in the primitive
//  variables and an input vector.  This operation converts primitive to characteristic
//  variables.  The result is returned in the input vector, with the components of the
//  characteristic field stored such that vect(1,i) is in the direction of the sweep.
//
//  The order of the components in the input vector should be:
//     (IDN,IVX,IVY,IVZ,[IPR],[IBY,IBZ])
//  and these are permuted according to the direction specified by the input flag "ivx".
//
// REFERENCES:
// - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new code for
//   astrophysical MHD", ApJS, (2008), Appendix A.  Equation numbers refer to this paper.

void Reconstruction::LeftEigenmatrixDotVectorCL(MeshBlock *pmb, const int ivx,
  const int il, const int iu, const AthenaArray<Real> &w,
  AthenaArray<Real> &vect) {
  // permute components of input primitive vector depending on direction
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::RightEigenmatrixDotVectorCL()
//  \brief Computes inner-product of right-eigenmatrix of Roe's matrix A in the primitive
//  variables and an input vector.  This operation converts characteristic to primitive
//  variables.  The result is returned in the input vector.
//
//  The order of the components in the input vector (characteristic fields) should be:
//     (IDN,ivx,ivy,ivz,[IPR],[IBY,IBZ])
//  where the lower-case indices indicate that the characteristic field in the direction
//  of the sweep (designated by the input flag "ivx") is stored first.  On output, the
//  components of velocity are in the standard order used for primitive variables.
//
// REFERENCES:
// - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new code for
//   astrophysical MHD", ApJS, (2008), Appendix A.  Equation numbers refer to this paper.

void Reconstruction::RightEigenmatrixDotVectorCL(MeshBlock *pmb, const int ivx,
  const int il, const int iu, const AthenaArray<Real> &w,
  AthenaArray<Real> &vect) {
  // permute components of output primitive vector depending on direction
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

   return;
}
