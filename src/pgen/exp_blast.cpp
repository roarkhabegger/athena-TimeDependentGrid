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
Real outerDens;
Real outerVel;
Real outerPres;
void OuterX1_UniformMedium(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
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
  Real shockR = gridData(0);
  Real shockV = gridData(1);
  
  if (dir==0) {
    retVal = shockV;
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
  for (int i=1;i < n-1; ++i) {
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
  int k=0
  for (int i=0; i< n-1; ++i) {
    if (shockData(i) >= 0.95){
      outArr[k] = i;
      k+=1    
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
    SetGridData(2);
    EnrollCalcGridData(UpdateGridData);
    
    outerDens = pin->GetReal("problem","damb");
    outerVel  = 0.0;
    outerPres = pin->GetReal("problem","pamb");
    if (mesh_bcs[OUTER_X1] == GetBoundaryFlag("user")) {
      EnrollUserBoundaryFunction(OUTER_X1,OuterX1_UniformMedium);
    }
    
    Real rout = pin->GetReal("problem","radius");
    Real rin  = rout - pin->GetOrAddReal("problem","ramp",0.0);
    Real vs   = pin->GetOrAddReal("problem","vel",0.0);
    GridData(0) = rin;
    GridData(1) = vs; 

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
  Real vs   = pin->GetOrAddReal("problem","vel",0.0);
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
    Real vel = 0.0;
    Real den = da;
    if (rad < rout) {
      if (rad < rin) {
        den = drat*da;
        vel = vs;
      } else {   // add smooth ramp in density
        Real f = (rad-rin) / (rout-rin);
        Real log_den = (1.0-f) * std::log(drat*da) + f * std::log(da);
        Real log_vel = (1.0-f) * std::log(vs) + f * std::log(vs*pow(10.0,-13));
        den = std::exp(log_den);
        vel = std::exp(log_vel);
      }
    }

    phydro->u(IDN,k,j,i) = den;
    phydro->u(IM1,k,j,i) = den*vel;
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
      phydro->u(IEN,k,j,i) = pres/gm1 + 0.5*den*pow(vel,2.0);
      if (DUAL_ENERGY) {
	phydro->u(IIE,k,j,i) = pres/gm1; 
      }
    }

    if (NSCALARS > 0) {
      for (int n=NHYDRO-NSCALARS; n<NHYDRO; ++n) {
        phydro->u(n,k,j,i) = den;
      }
    }
		
  }}}


  // initialize interface B and total energy
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie+1; ++i) {
          if (COORDINATE_SYSTEM == "cartesian") {
            pfield->b.x1f(k,j,i) = b0 * std::cos(angle);
          } else if (COORDINATE_SYSTEM == "cylindrical") {
            Real phi = pcoord->x2v(j);
            pfield->b.x1f(k,j,i) =
                b0 * (std::cos(angle) * std::cos(phi) + std::sin(angle) * std::sin(phi));
          } else { //if (COORDINATE_SYSTEM == "spherical_polar") {
            Real theta = pcoord->x2v(j);
            Real phi = pcoord->x3v(k);
            pfield->b.x1f(k,j,i) = b0 * std::abs(std::sin(theta))
                * (std::cos(angle) * std::cos(phi) + std::sin(angle) * std::sin(phi));
          }
        }
      }
    }
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je+1; ++j) {
        for (int i = is; i <= ie; ++i) {
          if (COORDINATE_SYSTEM == "cartesian") {
            pfield->b.x2f(k,j,i) = b0 * std::sin(angle);
          } else if (COORDINATE_SYSTEM == "cylindrical") {
            Real phi = pcoord->x2v(j);
            pfield->b.x2f(k,j,i) =
                b0 * (std::sin(angle) * std::cos(phi) - std::cos(angle) * std::sin(phi));
          } else { //if (COORDINATE_SYSTEM == "spherical_polar") {
            Real theta = pcoord->x2v(j);
            Real phi = pcoord->x3v(k);
            pfield->b.x2f(k,j,i) = b0 * std::cos(theta)
                * (std::cos(angle) * std::cos(phi) + std::sin(angle) * std::sin(phi));
            if (std::sin(theta) < 0.0)
              pfield->b.x2f(k,j,i) *= -1.0;
          }
        }
      }
    }
    for (int k = ks; k <= ke+1; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          if (COORDINATE_SYSTEM == "cartesian" || COORDINATE_SYSTEM == "cylindrical") {
            pfield->b.x3f(k,j,i) = 0.0;
          } else { //if (COORDINATE_SYSTEM == "spherical_polar") {
            Real phi = pcoord->x3v(k);
            pfield->b.x3f(k,j,i) =
                b0 * (std::sin(angle) * std::cos(phi) - std::cos(angle) * std::sin(phi));
          }
        }
      }
    }
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          phydro->u(IEN,k,j,i) += 0.5*b0*b0;
        }
      }
    }
  }
  return;
}


