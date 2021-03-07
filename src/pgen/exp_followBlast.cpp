//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file exp_default_pgen.cpp
//  \brief Provides default (empty) versions of all functions in problem generator files
//  This means user does not have to implement these functions if they are not needed.
//
//

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <iostream>
// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"

#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
//========================================================================================
// Time Dependent Grid Functions
//  \brief Functions for time dependent grid, including two example boundary conditions
//========================================================================================
Real WallVel(Real xf, int i, Real time, Real dt, int dir, AthenaArray<Real> gridData);
void UpdateGridData(Mesh *pm);

//Global Variables for OuterX1
Real ambDens;
Real ambVel;
Real ambPres;
void OuterX1_UniformMedium(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void OuterX2_UniformMedium(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void InnerX1_UniformMedium(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void InnerX2_UniformMedium(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void ShockDetector(AthenaArray<Real> data, AthenaArray<Real> grid, int outArr[], Real eps);

//========================================================================================
//! \fn void WallVel(Real xf, int i, Real time, Real dt, int dir, AthenaArray<Real> gridData)
//  \brief Function that returns the velocity of cell wall i at location xf. Time, total
//  time step and direction are all given. Direction is one of 0,1,2, corresponding to x1,x2,x3
//  and gridData is an athena array that contains overall mesh data. gridData is updated
//  before every time sub-step by the UpdateGridData function. Some instances do not need
//  this data to be updated and the UpdateGridData function can be left blank. The gridData
//  array is supposed to carry all mesh-level information, i.e. the information used for 
//  multiple cell walls in the simulation.
//========================================================================================
Real WallVel(Real xf, int i, Real time, Real dt, int dir, AthenaArray<Real> gridData) {
  Real retVal = 0.0;
  
 
  Real myX = xf;
  if (COORDINATE_SYSTEM == "cartesian") {
    if (dir == gridData(1)){
      if ((myX > 0.0)&&(gridData(3)>0.0)){ 
        if (gridData(2)==0.0) retVal = 0.0;
        else retVal = gridData(2) * myX/gridData(3);
      } else if ((myX < 0.0)&&(gridData(0)<0.0)){ 
        if (gridData(2) == 0.0) retVal = 0.0;
        else retVal = -1.0*gridData(2) * myX/gridData(0);
      }
    } else if (dir == gridData(5)) {
      if ((myX > 0.0)&&(gridData(7)>0.0)){ 
        if (gridData(6)==0.0) retVal = 0.0;
        else retVal = gridData(6) * myX/gridData(7);
      } else if ((myX < 0.0)&&(gridData(4)<0.0)){ 
        if (gridData(6) == 0.0) retVal = 0.0;
        else retVal = -1.0*gridData(6) * myX/gridData(4);
      }
    }   
  } else if (COORDINATE_SYSTEM == "cylindrical") {
    if (dir != gridData(1)){
      retVal = 0.0;
    } else if (myX<=gridData(0)){
      retVal = 0.0;
    } else if (myX > gridData(0)){ 
      if (gridData(2)==0.0) retVal = 0.0;
      else retVal = gridData(2) * (myX-gridData(0))/(gridData(3)-gridData(0));
    } 
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    if (dir != gridData(1)){
      retVal = 0.0;
    } else if (myX<=gridData(0)){
      retVal = 0.0;
    } else if (myX > gridData(0)){ 
      if (gridData(2)==0.0) retVal = 0.0;
      else retVal = gridData(2) * (myX-gridData(0))/(gridData(3)-gridData(0));
    } 
  }

  return retVal; 
}

//========================================================================================
//! \fn void UpdateGridData(Mesh *pm)
//  \brief Function which can edit and calculate any terms in gridData, which is used 
//  in the WallVel function. The object in mesh is GridData(i) and i can range over the
//  integers, limited by SetGridData argument in InitMeshUserData. See exp_blast for an 
//  example use of this function.
//========================================================================================
void UpdateGridData(Mesh *pm) {
  Real xMax;
  Real xMin;
  if (COORDINATE_SYSTEM == "cartesian") {
    xMax = pm->mesh_size.x1max;
    xMin = pm->mesh_size.x1min;
    pm->GridData(3) = xMax;
    pm->GridData(0) = xMin;

    xMax = pm->mesh_size.x2max;  
    xMin = pm->mesh_size.x2min;
    pm->GridData(7) = xMax;
    pm->GridData(4) = xMin;
    MeshBlock *pmb = pm->pblock;
    Real myVel = 0.0;
    Real pos  = 0.0;
    Real posUp = pmb->pcoord->x1v(((pmb->ie)-3));
    Real posLow = pmb->pcoord->x1v(((pmb->ie)-15));
    Real mom1 = 0.0;
    Real vol = 0.0;   
    //while (pmb != NULL) {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          pos = pow( pow( pmb->pcoord->x1v(i),2.0) + pow(pmb->pcoord->x2v(j),2.0),0.5);
          if ((pos<=posUp) && (pos>=posLow)) {
            mom1 += 5.0*pow( pow(pmb->phydro->u(IM1,k,j,i),2.0) 
                           + pow(pmb->phydro->u(IM2,k,j,i),2.0),0.5)
                      /  pmb->phydro->u(IDN,k,j,i)*pmb->pcoord->GetCellVolume(k,j,i);
            vol += pmb->pcoord->GetCellVolume(k,j,i);
          }                  
        }
      }
    }
    myVel = mom1/vol*(pm->GridData(3)/(posLow));
    //std::cout << myVel << std::endl;
    if ((myVel <=0.0)) {
      myVel = 0.0;
    }
    pm->GridData(2) = myVel;
    pm->GridData(6) = myVel;
  } else {
    xMax = pm->mesh_size.x1max;
    pm->GridData(3) = xMax;
    MeshBlock *pmb = pm->pblock;
    Real myVel = 0.0;
    Real pos  = 0.0;
    Real posUp = pmb->pcoord->x1v(((pmb->ie)-3));
    Real posLow = pmb->pcoord->x1v(((pmb->ie)-15));
    Real mom1 = 0.0;
    Real vol = 0.0;   
    //while (pmb != NULL) {
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            pos = pmb->pcoord->x1v(i);
            if ((pos<=posUp) && (pos>=posLow)) {
              mom1 += 5.0*pmb->phydro->u(IM1,k,j,i)/pmb->phydro->u(IDN,k,j,i)*pmb->pcoord->GetCellVolume(k,j,i);
              vol += pmb->pcoord->GetCellVolume(k,j,i);
            }
                        
          }
        }
      }
    myVel = mom1/vol*(pm->GridData(3)/(posLow));
    //std::cout << myVel << std::endl;
    if ((myVel <=0.0)) {
      myVel = 0.0;
    }
    pm->GridData(2) = myVel;
  }

   

  return;
}


//Harten Van Leer Shock detection algorithm Out data should be a 1 dimensional array,
// with the same length as indata and grid. indata is the array of Real values where
// we look for the shocks. eps is the slope magnitude limiter, i.e. if the slope is
// above eps, then the location has a shock.
void ShockDetector(AthenaArray<Real> data, AthenaArray<Real> grid, int outArr[], Real eps ) {
  int n, loc;
  Real a, b, c;
  n = data.GetDim1();
  AthenaArray<Real> shockData;
  shockData.NewAthenaArray(n-1);
  loc = 0;
  for (int i=1; i<(n-1); ++i) {
    a = 0;
    b = 0;
    a = std::abs(data(i)-data(i-1));
    b = std::abs(data(i+1)-data(i));
    c = a+b;
    shockData(i-1) = std::pow(a-b,2);

    if ( c <= eps) { 
      shockData(i-1) = 0.0;
    } else { 
      shockData(i-1) *= std::pow(a+b,-2);
    }  
  }  
  int k=0;
  for (int i=0; i< (n-1); ++i) {
    if (shockData(i) >= 0.95){
      outArr[k] = i;
      k+=1;    
    }

  }
   
 
  shockData.DeleteAthenaArray();
  return;
}

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in Mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  //========================================================================================
  //! \brief For a time dependent grid, make sure to use SetGridData, EnrollGridDiffEq, and
  //   EnrollCalcGridData here. The boundary conditions are of course optional. Reflecting 
  //   is a good boundary function if a wall of the simulation is static. But if there is
  //   any expansion of the grid, it is recommended that you use the UniformMedium condition
  //   for the expanding boundary. Otherwise, reconstruction might fail because the data is
  //   inaccurate (for example, periodic boundary conditions do not make sense
  //   for an expanding grid).
  //========================================================================================
  if (EXPANDING) {
    EnrollGridDiffEq(WallVel);
      
    if (COORDINATE_SYSTEM == "cartesian") {
      SetGridData(8);
    } else {
      SetGridData(4);
    }
    EnrollCalcGridData(UpdateGridData);
    
    ambDens = pin->GetReal("problem","damb");
    ambVel  = 0.0;
    ambPres = pin->GetReal("problem","pamb");
    if (mesh_bcs[OUTER_X1] == GetBoundaryFlag("user")) {
      EnrollUserBoundaryFunction(OUTER_X1,OuterX1_UniformMedium);
    }
    if (mesh_bcs[OUTER_X2] == GetBoundaryFlag("user")) {
      EnrollUserBoundaryFunction(OUTER_X2,OuterX2_UniformMedium);
    }
    if (mesh_bcs[INNER_X1] == GetBoundaryFlag("user")) {
      EnrollUserBoundaryFunction(INNER_X1,InnerX1_UniformMedium);
    }
    if (mesh_bcs[INNER_X2] == GetBoundaryFlag("user")) {
      EnrollUserBoundaryFunction(INNER_X2,InnerX2_UniformMedium);
    }
    
    Real rout = pin->GetReal("problem","radius");
    Real rin  = rout - pin->GetOrAddReal("problem","ramp",0.0);
    Real vs   = pin->GetOrAddReal("problem","vel",0.0);

    if (COORDINATE_SYSTEM == "cartesian") {
      GridData(0) = mesh_size.x1min;
      GridData(1) = 1; 
      GridData(2) = 0.0;
      GridData(3) = mesh_size.x1max; 
      GridData(4) = mesh_size.x2min;
      GridData(5) = 2; 
      GridData(6) = 0.0;
      GridData(7) = mesh_size.x2max; 
    } else {
      GridData(0) = mesh_size.x1min;
      GridData(1) = 1; 
      GridData(2) = 0.0;
      GridData(3) = mesh_size.x1max; 
    }
  }

  return;
}

//========================================================================================
//! \fn void OuterX1_UniformMedium(MeshBlock *pmb, Coordinates *pco, 
//                                 AthenaArray<Real> &prim,FaceField &b, Real time,
//                                 Real dt, int is, int ie, int js, int je,
//                                 int ks, int ke, int ngh) {
//  \brief Function for outer boundary being a uniform medium with density, velocity,
//   and pressure given by the global variables listed at the beginning of the file.
//========================================================================================
void OuterX1_UniformMedium(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        prim(IVX,k,j,ie+i) = ambVel;
        prim(IDN,k,j,ie+i) = ambDens;
        prim(IPR,k,j,ie+i) = ambPres;  
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


//========================================================================================
//! \fn void OuterX2_UniformMedium(MeshBlock *pmb, Coordinates *pco, 
//                                 AthenaArray<Real> &prim,FaceField &b, Real time,
//                                 Real dt, int is, int ie, int js, int je,
//                                 int ks, int ke, int ngh) {
//  \brief Function for outer boundary being a uniform medium with density, velocity,
//   and pressure given by the global variables listed at the beginning of the file.
//========================================================================================
void OuterX2_UniformMedium(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  
  for (int k=ks; k<=ke; ++k) {
    for (int i=is; i<=ie; ++i) {
#pragma omp simd
      for (int j=1; j<=ngh; ++j) {
        prim(IVX,k,je+j,i) = ambVel;
        prim(IDN,k,je+j,i) = ambDens;
        prim(IPR,k,je+j,i) = ambPres;  
        prim(IVY,k,je+j,i) = 0.0;
        prim(IVZ,k,je+j,i) = 0.0;
      }
    }
  }

  // no magnetic fields in ambient medium
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int i=is; i<=ie+1; ++i) {
#pragma omp simd
        for (int j=1; j<=ngh; ++j) {
          b.x1f(k,je+j,i) = 0.0;  
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int i=is; i<=ie; ++i) {
#pragma omp simd
        for (int j=1; j<=ngh; ++j) {
          b.x2f(k,je+j,i) = 0.0;
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int i=is; i<=ie; ++i) {
#pragma omp simd
        for (int j=1; j<=ngh; ++j) {
          b.x3f(k,je+j,i) =  0.0;
        }
      }
    }
  }
  return;

}

//========================================================================================
//! \fn void InnerX1_UniformMedium(MeshBlock *pmb, Coordinates *pco, 
//                                 AthenaArray<Real> &prim,FaceField &b, Real time,
//                                 Real dt, int is, int ie, int js, int je,
//                                 int ks, int ke, int ngh) {
//  \brief Function for inner boundary being a uniform medium with density, velocity,
//   and pressure given by the global variables listed at the beginning of the file.
//========================================================================================

void InnerX1_UniformMedium(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        prim(IVX,k,j,is-i) = ambVel;
        prim(IDN,k,j,is-i) = ambDens;
        prim(IPR,k,j,is-i) = ambPres;  
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

//========================================================================================
//! \fn void InnerX2_UniformMedium(MeshBlock *pmb, Coordinates *pco, 
//                                 AthenaArray<Real> &prim,FaceField &b, Real time,
//                                 Real dt, int is, int ie, int js, int je,
//                                 int ks, int ke, int ngh) {
//  \brief Function for inner boundary being a uniform medium with density, velocity,
//   and pressure given by the global variables listed at the beginning of the file.
//========================================================================================

void InnerX2_UniformMedium(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        prim(IVX,k,js-j,i) = ambVel;
        prim(IDN,k,js-j,i) = ambDens;
        prim(IPR,k,js-j,i) = ambPres;  
        prim(IVY,k,js-j,i) = 0.0;
        prim(IVZ,k,js-j,i) = 0.0;
      }
    }
  }

  // no magnetic fields in ambient medium
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          b.x1f(k,js-j,i) = 0.0;  
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x2f(k,js-j,i) = 0.0;
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

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Should be used to set initial conditions.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // In practice, this function should *always* be replaced by a version
  // that sets the initial conditions for the problem of interest.
  Real rout = pin->GetReal("problem","radius");
  Real rin  = rout - pin->GetOrAddReal("problem","ramp",0.0);
  Real pa   = pin->GetOrAddReal("problem","pamb",1.0);
  Real da   = pin->GetOrAddReal("problem","damb",1.0);
  Real prat = pin->GetReal("problem","prat");
  Real drat = pin->GetOrAddReal("problem","drat",1.0);
  Real vSh   = pin->GetOrAddReal("problem","vel",0.0);
  Real b0,angle;
  if (MAGNETIC_FIELDS_ENABLED) {
    b0 = pin->GetReal("problem","b0");
    angle = (PI/180.0)*pin->GetReal("problem","angle");
  }
  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;

  // get coordinates of center of blast, and convert to Cartesian if necessary
  Real x1_0   = pin->GetOrAddReal("problem","x1_0",0.0);
  Real x2_0   = pin->GetOrAddReal("problem","x2_0",0.0);
  Real x3_0   = pin->GetOrAddReal("problem","x3_0",0.0);
  Real x0,y0,z0;
  if (COORDINATE_SYSTEM == "cartesian") {
    x0 = x1_0;
    y0 = x2_0;
    z0 = x3_0;
  } else if (COORDINATE_SYSTEM == "cylindrical") {
    x0 = x1_0*std::cos(x2_0);
    y0 = x1_0*std::sin(x2_0);
    z0 = x3_0;
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    x0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
    y0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
    z0 = x1_0*std::cos(x2_0);
  } else {
    // Only check legality of COORDINATE_SYSTEM once in this function
    std::cout << "### FATAL ERROR in blast.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM= " << COORDINATE_SYSTEM << std::endl;
  }

  // setup uniform ambient medium with spherical over-pressured region
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    Real rad;
    if (COORDINATE_SYSTEM == "cartesian") {
      Real x = pcoord->x1v(i);
      Real y = pcoord->x2v(j);
      Real z = pcoord->x3v(k);
      rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
    } else if (COORDINATE_SYSTEM == "cylindrical") {
      Real x = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
      Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
      Real z = pcoord->x3v(k);
      rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
    } else { // if (COORDINATE_SYSTEM == "spherical_polar")
      Real x = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
      Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
      Real z = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
      rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
    }
    Real den = da;
    Real v1  = 0.0;
    if (rad < rout) {
      v1 = vSh;
      if (rad < rin) {
        den = drat*da;
      } else {   // add smooth ramp in density
        Real f = (rad-rin) / (rout-rin);
        Real log_den = (1.0-f) * std::log(drat*da) + f * std::log(da);
        den = std::exp(log_den);
      }
    }

    phydro->u(IDN,k,j,i) = den;
    phydro->u(IM1,k,j,i) = den*v1;
    phydro->u(IM2,k,j,i) = 0.0;
    phydro->u(IM3,k,j,i) = 0.0;
    if (NON_BAROTROPIC_EOS) {
      Real pres = pa;
      if (rad < rout) {
        if (rad < rin) {
          pres = prat*pa;
        } else {  // add smooth ramp in pressure
          Real f = (rad-rin) / (rout-rin);
          Real log_pres = (1.0-f) * std::log(prat*pa) + f * std::log(pa);
          pres = std::exp(log_pres);
        }
      }
      phydro->u(IEN,k,j,i) = 0.5*den*pow(v1,2.0)+pres/gm1;
    }

    if (NSCALARS > 0) {
      for (int n=NHYDRO-NSCALARS; n<NHYDRO; ++n) {
        phydro->u(n,k,j,i) = den;
      }
    }
		
  }}}

  return;
}


