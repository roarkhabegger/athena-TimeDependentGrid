//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file shock_tube.cpp
//  \brief Problem generator for shock tube problems.
//
// Problem generator for shock tube (1-D Riemann) problems. Initializes plane-parallel
// shock along x1 (in 1D, 2D, 3D), along x2 (in 2D, 3D), and along x3 (in 3D).
//========================================================================================

// C headers
#include <stdio.h>

// C++ headers
#include <cmath>      // sqrt()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../cless/cless.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"


Real WallVel(Real xf, int i, Real time, Real dt, int dir, AthenaArray<Real> gridData);
void UpdateGridData(Mesh *pm);
void OuterX1_UniformMedium(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void InnerX1_UniformMedium(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

//====================================================================================
// IMplement Expanding Functions

//Wall Velocity. Depending on direction and GridData and time, return the velocity at xf, which
//will be the location of the ith cell wall at time time
Real WallVel(Real xf, int i, Real time, Real dt, int dir, AthenaArray<Real> gridData) {
  Real retVal = 0.0;
  Real x0l = gridData(3);
  Real x0r = gridData(5);
  Real myX = xf;
  //std::cout << gridData(0) << std::endl;
  if (dir != gridData(1)){
    retVal = 0.0;
  } else if (myX==gridData(0)){
    retVal = 0.0;
  } else if (myX < gridData(0)){
    if (gridData(2)==0.0) retVal = 0.0;
    else retVal = gridData(2) * (gridData(0)-myX)/(gridData(0)-x0l);
  } else if (myX > gridData(0)){ 
    if (gridData(4)==0.0) retVal = 0.0;
    else retVal = gridData(4) * (myX-gridData(0))/(x0r-gridData(0));
  } 
  return retVal; 
}

void UpdateGridData(Mesh *pm) {
   

    Real xMin;
    Real xMax;

    int shk_dir = pm->GridData(1);
    switch(shk_dir) {
      //--- shock in 1-direction
      case 1:
        xMin = pm->mesh_size.x1min;
        xMax = pm->mesh_size.x1max;
      break;
      //--- shock in 2-direction
      case 2:
        xMin = pm->mesh_size.x2min;
        xMax = pm->mesh_size.x2max;
      break;
      //--- shock in 3-direction
      case 3:   
        xMin = pm->mesh_size.x3min;
        xMax = pm->mesh_size.x3max;
      break;
    }

    pm->GridData(3) = xMin;
    pm->GridData(5) = xMax;

  return;
}

//Global Variables for OuterX1
Real outerDens;
Real outerVel;
Real outerPres;

void OuterX1_UniformMedium(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        prim(IVX,k,j,ie+i) = outerVel;
        prim(IDN,k,j,ie+i) = outerDens;
        prim(IPR,k,j,ie+i) = outerPres;  
        prim(IVY,k,j,ie+i) = 0.0;
        prim(IVZ,k,j,ie+i) = 0.0;
      }
    }
  }

  // no magnetic fields in ambient medium
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x1f(k,j,(ie+i)) = 0.0;  
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x2f(k,j,(ie+i)) = 0.0;
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x3f(k,j,(ie+i)) =  0.0;
        }
      }
    }
  }
  return;

}

//Global Variables for OuterX1
Real innerDens;
Real innerVel;
Real innerPres;

void InnerX1_UniformMedium(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        prim(IVX,k,j,is-i) = innerVel;
        prim(IDN,k,j,is-i) = innerDens;
        prim(IPR,k,j,is-i) = innerPres;  
        prim(IVY,k,j,is-i) = 0.0;
        prim(IVZ,k,j,is-i) = 0.0;
      }
    }
  }

  // no magnetic fields in ambient medium
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x1f(k,j,(is-i)) = 0.0;  
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x2f(k,j,(is-i)) = 0.0;
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x3f(k,j,(is-i)) =  0.0;
        }
      }
    }
  }
  return;

}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  
  if (EXPANDING) {
    outerDens = pin->GetReal("problem","dr");
    outerVel  = pin->GetReal("problem","vr");
    outerPres = pin->GetReal("problem","pr");
    innerDens = pin->GetReal("problem","dl");
    innerVel  = pin->GetReal("problem","vl");
    innerPres = pin->GetReal("problem","pl");
    SetGridData(6);
    EnrollGridDiffEq(WallVel);
    EnrollCalcGridData(UpdateGridData);
    
    if (mesh_bcs[OUTER_X1] == GetBoundaryFlag("user")) {
      EnrollUserBoundaryFunction(OUTER_X1,OuterX1_UniformMedium);
    }
    if (mesh_bcs[INNER_X1] == GetBoundaryFlag("user")) {
      EnrollUserBoundaryFunction(INNER_X1,InnerX1_UniformMedium);
    }

    int shk_dir = pin->GetInteger("problem","shock_dir");
    GridData(1) = pin->GetReal("problem","shock_dir");
    Real xMin;
    Real xMax;
    Real n;
    Real vGridL = 0.0;
    Real vGridR = 0.0;
    Real tLim = pin->GetReal("time","tlim");

    switch(shk_dir) {
      //--- shock in 1-direction
      case 1:
        xMin = pin->GetReal("mesh","x1min");
        xMax = pin->GetReal("mesh","x1max");
        n = (Real)(pin->GetInteger("mesh","nx1"));
      break;
      //--- shock in 2-direction
      case 2:
        xMin = pin->GetReal("mesh","x2min");
	n = pin->GetReal("mesh","nx2");
        xMax = pin->GetReal("mesh","x2max");
      break;
      //--- shock in 3-direction
      case 3:   
        xMin = pin->GetReal("mesh","x3min");
        xMax = pin->GetReal("mesh","x3max");
        n = pin->GetReal("mesh","nx3");
      break;
    }
    GridData(0) = 0.0; //(xMin+xMax)/n*0.5;
    Real Fudge = pin->GetOrAddReal("problem","Fudge",1.0);
    
    Real maxXR = Fudge*(pin->GetReal("problem","maxRight"));
    Real maxXL = Fudge*(pin->GetReal("problem","maxLeft"));
    vGridL = (maxXL-xMin)/tLim;
    vGridR = (maxXR-xMax)/tLim;
    GridData(2) = vGridL;
    GridData(3) = xMin; 
    GridData(4) = vGridR;
    GridData(5) = xMax;
    
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the shock tube tests
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  std::stringstream msg;

  // parse shock direction: {1,2,3} -> {x1,x2,x3}
  int shk_dir = pin->GetInteger("problem","shock_dir");

  // parse shock location (must be inside grid)
  Real xshock = pin->GetReal("problem","xshock");
  if (shk_dir == 1 && (xshock < pmy_mesh->mesh_size.x1min ||
                       xshock > pmy_mesh->mesh_size.x1max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x1 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (shk_dir == 2 && (xshock < pmy_mesh->mesh_size.x2min ||
                       xshock > pmy_mesh->mesh_size.x2max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x2 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (shk_dir == 3 && (xshock < pmy_mesh->mesh_size.x3min ||
                       xshock > pmy_mesh->mesh_size.x3max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x3 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // Parse left state read from input file: dl,ul,vl,wl,[pl]
  Real wl[NHYDRO+NFIELD];
  wl[IDN] = pin->GetReal("problem","dl");
  wl[IVX] = pin->GetReal("problem","ul");
  wl[IVY] = pin->GetReal("problem","vl");
  wl[IVZ] = pin->GetReal("problem","wl");
  if (NON_BAROTROPIC_EOS) wl[IPR] = pin->GetReal("problem","pl");
  if (MAGNETIC_FIELDS_ENABLED) {
    wl[NHYDRO  ] = pin->GetReal("problem","bxl");
    wl[NHYDRO+1] = pin->GetReal("problem","byl");
    wl[NHYDRO+2] = pin->GetReal("problem","bzl");
  }

  // Parse right state read from input file: dr,ur,vr,wr,[pr]
  Real wr[NHYDRO+NFIELD];
  wr[IDN] = pin->GetReal("problem","dr");
  wr[IVX] = pin->GetReal("problem","ur");
  wr[IVY] = pin->GetReal("problem","vr");
  wr[IVZ] = pin->GetReal("problem","wr");
  if (NON_BAROTROPIC_EOS) wr[IPR] = pin->GetReal("problem","pr");
  if (MAGNETIC_FIELDS_ENABLED) {
    wr[NHYDRO  ] = pin->GetReal("problem","bxr");
    wr[NHYDRO+1] = pin->GetReal("problem","byr");
    wr[NHYDRO+2] = pin->GetReal("problem","bzr");
  }

  
   
// Initialize the discontinuity in the Hydro variables ---------------------------------

  switch(shk_dir) {

//--- shock in 1-direction
  case 1:
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        if (pcoord->x1v(i) < xshock) {
          phydro->u(IDN,k,j,i) = wl[IDN];
          phydro->u(IM1,k,j,i) = wl[IVX]*wl[IDN];
          phydro->u(IM2,k,j,i) = wl[IVY]*wl[IDN];
          phydro->u(IM3,k,j,i) = wl[IVZ]*wl[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wl[IPR]/(peos->GetGamma() - 1.0)
            + 0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY] + wl[IVZ]*wl[IVZ]);
				
        } else {
          phydro->u(IDN,k,j,i) = wr[IDN];
          phydro->u(IM1,k,j,i) = wr[IVX]*wr[IDN];
          phydro->u(IM2,k,j,i) = wr[IVY]*wr[IDN];
          phydro->u(IM3,k,j,i) = wr[IVZ]*wr[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wr[IPR]/(peos->GetGamma() - 1.0)
            + 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY] + wr[IVZ]*wr[IVZ]);
					
        }
      }
    }}
    break;

//--- shock in 2-direction
  case 2:
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      if (pcoord->x2v(j) < xshock) {
        for (int i=is; i<=ie; ++i) {
          phydro->u(IDN,k,j,i) = wl[IDN];
          phydro->u(IM2,k,j,i) = wl[IVX]*wl[IDN];
          phydro->u(IM3,k,j,i) = wl[IVY]*wl[IDN];
          phydro->u(IM1,k,j,i) = wl[IVZ]*wl[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wl[IPR]/(peos->GetGamma() - 1.0)
            + 0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY] + wl[IVZ]*wl[IVZ]);
        }
      } else {
        for (int i=is; i<=ie; ++i) {
          phydro->u(IDN,k,j,i) = wr[IDN];
          phydro->u(IM2,k,j,i) = wr[IVX]*wr[IDN];
          phydro->u(IM3,k,j,i) = wr[IVY]*wr[IDN];
          phydro->u(IM1,k,j,i) = wr[IVZ]*wr[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wr[IPR]/(peos->GetGamma() - 1.0)
            + 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY] + wr[IVZ]*wr[IVZ]);
        }
      }
    }}
    break;

//--- shock in 3-direction

  case 3:
    for (int k=ks; k<=ke; ++k) {
      if (pcoord->x3v(k) < xshock) {
        for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          phydro->u(IDN,k,j,i) = wl[IDN];
          phydro->u(IM3,k,j,i) = wl[IVX]*wl[IDN];
          phydro->u(IM1,k,j,i) = wl[IVY]*wl[IDN];
          phydro->u(IM2,k,j,i) = wl[IVZ]*wl[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wl[IPR]/(peos->GetGamma() - 1.0)
            + 0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY] + wl[IVZ]*wl[IVZ]);
        }}
      } else {
        for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          phydro->u(IDN,k,j,i) = wr[IDN];
          phydro->u(IM3,k,j,i) = wr[IVX]*wr[IDN];
          phydro->u(IM1,k,j,i) = wr[IVY]*wr[IDN];
          phydro->u(IM2,k,j,i) = wr[IVZ]*wr[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wr[IPR]/(peos->GetGamma() - 1.0)
            + 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY] + wr[IVZ]*wr[IVZ]);
        }}
      }
    }
    break;

  default:
    msg << "### FATAL ERROR in Problem Generator" << std::endl
        << "shock_dir=" << shk_dir << " must be either 1,2, or 3" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// now set face-centered (interface) magnetic fields -----------------------------------

  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        if (shk_dir==1 && pcoord->x1v(i) < xshock) {
          pfield->b.x1f(k,j,i) = wl[NHYDRO  ];
          pfield->b.x2f(k,j,i) = wl[NHYDRO+1];
          pfield->b.x3f(k,j,i) = wl[NHYDRO+2];
        } else if (shk_dir==2 && pcoord->x2v(j) < xshock) {
          pfield->b.x1f(k,j,i) = wl[NHYDRO+2];
          pfield->b.x2f(k,j,i) = wl[NHYDRO  ];
          pfield->b.x3f(k,j,i) = wl[NHYDRO+1];
        } else if (shk_dir==3 && pcoord->x3v(k) < xshock) {
          pfield->b.x1f(k,j,i) = wl[NHYDRO+1];
          pfield->b.x2f(k,j,i) = wl[NHYDRO+2];
          pfield->b.x3f(k,j,i) = wl[NHYDRO];
        }

        if (shk_dir==1 && pcoord->x1v(i) >= xshock) {
          pfield->b.x1f(k,j,i) = wr[NHYDRO  ];
          pfield->b.x2f(k,j,i) = wr[NHYDRO+1];
          pfield->b.x3f(k,j,i) = wr[NHYDRO+2];
        } else if (shk_dir==2 && pcoord->x2v(j) >= xshock) {
          pfield->b.x1f(k,j,i) = wr[NHYDRO+2];
          pfield->b.x2f(k,j,i) = wr[NHYDRO  ];
          pfield->b.x3f(k,j,i) = wr[NHYDRO+1];
        } else if (shk_dir==3 && pcoord->x3v(k) >= xshock)  {
          pfield->b.x1f(k,j,i) = wr[NHYDRO+1];
          pfield->b.x2f(k,j,i) = wr[NHYDRO+2];
          pfield->b.x3f(k,j,i) = wr[NHYDRO];
        }
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) += 0.5*(SQR(pfield->b.x1f(k,j,i))
            + SQR(pfield->b.x2f(k,j,i)) + SQR(pfield->b.x3f(k,j,i)));
        }
      }
    }}

    // end by adding bi.x1 at ie+1, bi.x2 at je+1, and bi.x3 at ke+1

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      pfield->b.x1f(k,j,ie+1) = pfield->b.x1f(k,j,ie);
    }}
    for (int k=ks; k<=ke; ++k) {
    for (int i=is; i<=ie; ++i) {
      pfield->b.x2f(k,je+1,i) = pfield->b.x2f(k,je,i);
    }}
    for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      pfield->b.x3f(ke+1,j,i) = pfield->b.x3f(ke,j,i);
    }}
  }
  return;
}
