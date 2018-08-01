//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief Class to implement source terms in the hydro equations

// C/C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "cless_srcterms.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../cless.hpp"
#include "../../parameter_input.hpp"

// ClessSourceTerms constructor

ClessSourceTerms::ClessSourceTerms(Cless *pcl, ParameterInput *pin) {
  pmy_cless_ = pcl;
  cless_sourceterms_defined = false;

  // read point mass or constant acceleration parameters from input block

  if (SELF_GRAVITY_ENABLED) cless_sourceterms_defined = true;

	if (DUAL_ENERGY) cless_sourceterms_defined = true; 

  UserSourceTerm = pcl->pmy_block->pmy_mesh->UserSourceTerm_;
  if (UserSourceTerm != NULL) cless_sourceterms_defined = true;
}

// destructor

ClessSourceTerms::~ClessSourceTerms() {
}

//----------------------------------------------------------------------------------------
//! \fn void ClessSourceTerms::AddClessSourceTerms
//  \brief Adds source terms to conserved variables

void ClessSourceTerms::AddClessSourceTerms(const Real time, const Real dt,
     const AthenaArray<Real> *flux, const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bcc, AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_cless_->pmy_block;

  // Add new source terms here
  if (SELF_GRAVITY_ENABLED) SelfGravity(dt, flux, prim, cons);
	
  // MyNewSourceTerms()
  //  user-defined source terms
  if (UserSourceTerm != NULL)
    UserSourceTerm(pmb, time,dt,prim,bcc,cons);

  return;
}

