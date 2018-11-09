//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file advect.cpp
//  \brief Simple advection problem, to test scalars
//
//

// C++ headers
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <iostream>
#include <iomanip>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp" 
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"


//====================================================================================
// Enroll user-specific functions
void Mesh::InitUserMeshData(ParameterInput *pin) {
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Spherical blast wave test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real v0      = 1.0;
  Real d0      = 1.0;
  Real d1      = 2.0;
  Real p0      = 1.0;
  Real gamma   = peos->GetGamma();
  Real gm1     = gamma - 1.0;

  // setup uniform ambient medium with spherical over-pressured region
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real x   = pcoord->x1v(i);
        Real y   = pcoord->x2v(j);
        Real z   = pcoord->x3v(k);
        Real prf = 0.25*(1.0+std::tanh((x+0.25)/0.05))*(1.0-std::tanh((x-0.25)/0.05));
        phydro->u(IDN,k,j,i) = d0+(d1-d0)*prf;
        phydro->u(IM1,k,j,i) = v0*phydro->u(IDN,k,j,i);
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = p0/gm1+0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          if (DUAL_ENERGY) {
            phydro->u(IIE,k,j,i) = p0/gm1;
          }
        }
        if (NSCALARS == 2) {
          phydro->u(NHYDRO-NSCALARS  ,k,j,i) = d1*prf;
          phydro->u(NHYDRO-NSCALARS+1,k,j,i) = d1*(1.0-prf);
        }
      }
    }
  }
}

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//========================================================================================
void MeshBlock::UserWorkInLoop(void) {
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief Check radius of sphere to make sure it is round
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  return;
}
