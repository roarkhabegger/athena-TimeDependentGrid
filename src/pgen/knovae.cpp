//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file knovae.cpp
//  \brief Problem generator for spherical knova  problem.  
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
#include "../bvals/bvals.hpp"
#include "../expansion/expansion.hpp"
//====================================================================================
// global variables

#ifdef MPI_PARALLEL
typedef struct MPI_Comm_Sub {
  MPI_Group gsub, gworld;
  MPI_Comm  comsub;
} MPI_Comm_Sub;

MPI_Comm_Sub comm_x1;
MPI_Comm_Sub comm_slab;
#endif

//====================================================================================
// local functions
Real LogMeshSpacingX1(Real x, RegionSize rs);

void ReflectInnerX1_nonuniform(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

Real WallVel(Real xf, int i, Real time, Real dt, int dir, AthenaArray<Real> gridData);
void UpdateGridData(Mesh *pm);
int ShockDetector(AthenaArray<Real> data, AthenaArray<Real> grid, Real eps);


void ExpandingOuterX1_UniformMedium(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

//====================================================================================
// IMplement Expanding Functions

//New face coord. Depending on direction and GridAData and time, return where the cell
//wall at xf needs to move to. I.e. given xf_n return xf_{n+1} the position at the 
//new time step
Real WallVel(Real xf, int i, Real time, Real dt, int dir, AthenaArray<Real> gridData) {
  Real retVal = 0.0;
  //std::cout << xf << std::endl;
  // 0 -> x0, 1 -> shock Location, 2 -> Shock Velocity, 3 -> half of comoving width     

  if (dir != 0){
    retVal += 0.0;
  } else if (xf==gridData(0)){
    retVal += 0.0;
  } else if (gridData(1) == gridData(0)){
    retVal += 0.0;
  //} else if ((xf>=gridData(1)-gridData(3))&(xf<=gridData(1)+gridData(3))) {
 //   retVal = dt*gridData(2);
  //} else if (fabs(xf-gridData(1)) <= gridData(3)) {
  //  retVal += gridData(2)*dt;
  } else { 
    //Real x0 = LockData(0);
    //Real alpha = LockData(1);
    //std::cout << "alpha " << alpha << std::endl;
    //Real delAlpha = (dT)*(LockData(2));
    //if in x-direction
    retVal = gridData(2)*(xf-gridData(0)) / (gridData(1)-gridData(0));
    //std::cout << Vel  <<std::endl;
    //std::cout << dT << std::endl;
  } 
  //std::cout << retVal - xf << std::endl; 
  //if (retVal < 0) std::cout << "Got a negative" << std::endl;  
  return retVal; 
}

//Update grid data. Called once per stage to take in entire mesh and grid data. Main
//use is shock location for following a shock through a simulation. This data is the 
//first expansion function called during a given time step and the data is used 
//throughout the process of calcualting source terms and moving the grid
//NOTE: You have access to Mesh here as well as the ExpGridData array in Mesh which
//is where you should store it
void UpdateGridData(Mesh *pm) {
  //This gets called once per mesh dt
  //std::cout << "Updating Grid Data" << std::endl; 

  //Real xu = pm->mesh_size.x1min;
  //Real xl = pm->mesh_size.x1max;
  int n1 =  pm->mesh_size.nx1;
  int ng = (NGHOST);
  AthenaArray<Real> myData;
  AthenaArray<Real> myGrid;
  AthenaArray<Real> myVel;
  myData.NewAthenaArray(n1);
  myGrid.NewAthenaArray(n1);
  myVel.NewAthenaArray(n1);

  MeshBlock *pmb = pm->pblock;  
  //Should be on first meshBlock at this point.
  //Now go through meshblocks and average over other dimensions.
  int myI = 0;
  Real del = 100.0;
  Real Dye  = 0.0;
  Real DyeMoment = 0.0;

  for (int i=pmb->is;i<=pmb->ie;++i) {
    //std::cout << i << " "<< pmb->phydro->u(IDN,0,0,i) <<std::endl;
    myData(i-ng) = pmb->phydro->u(IDN,pmb->ks,pmb->js,i);
    myGrid(i-ng) = pmb->pcoord->x1v(i);
    myVel(i-ng)  = pmb->phydro->u(IVX,pmb->ks,pmb->js,i)/myData(i-ng);
    //Dye += pmb->pcoord->dx1f(i)*myData(i);
    //DyeMoment += pmb->pcoord->dx1f(i)*myData(i)*myGrid(i);
  }
//  while(pmb->next !=NULL) {
  //  std::cout << "More than one MeshBlock" << std::endl;
  //  pmb = pmb->next;  
  //  for (int i=pmb->is;i<=pmb->ie;++i) {
      //myData(i) += pmb->phydro->u(IDN,pmb->ks,pmb->js,i);
      //myGrid(i) = pmb->pcoord->x1v(i);
      //myVel(i)  = pmb->phydro->u(IVX,pmb->ks,pmb->js,i);
  //  }
  //}
  myI = ShockDetector(myData, myGrid, 1.0);
  Real vel = 0.0;
  for (int i=myI-ng;i<=myI+ng;++i){
    vel += myVel(i);
  }
  vel *= 1.0/(2.0*ng);
  //Real CoM = DyeMoment/Dye;
  //for (int i=pmb->is;i<=pmb->ie;++i) {
    //DyeMoment += pmb->pcoord->dx1f(i)*myData(i);
    //if (DyeMoment/Dye > 0.9){
      //myI = i;
      //break;
    //}
    //if (del > fabs(myGrid(i) - CoM)) {
    //  del = fabs(myGrid(i) - CoM);
    //  myI = i;
    //}
  //}
  //pm->GridData(1) = myGrid(myI);
  //pm->GridData(2) = vel; 
  std::cout << vel << std::endl;
  myVel.DeleteAthenaArray();
  myData.DeleteAthenaArray();
  myGrid.DeleteAthenaArray();

  return; 
  //Real dt = (stage_wghts[(stage-1)].beta)*(pm->dt);
  //pm->ExpMeshSize(0) = pmb->pmy_mesh->NewXf_(pm->mesh_size.x1min,pm->time,dt,0,pm->ExpGridData) - pm->mesh_size.x1min;

//  std::cout << myI << std::endl;
//  std::cout << pm->ExpGridData(1) << std::endl;     
//  std::cout << pm->ExpGridData(2) << std::endl;     
}

//Harten Van Leer Chock detection algorithm Out data should be a 1 dimensional array,
// with the same length as indata and grid. indata is the array of Real values where
// we look for the shocks. eps is the slope magnitude limiter, i.e. if the slope is
// above eps, then the location has a shock.
int ShockDetector(AthenaArray<Real> data, AthenaArray<Real> grid, Real eps ) {
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
  
  for (int i=0; i< n-1; ++i) {
    if (shockData(i) >= 0.95){
      loc = i;    
    }

  }
  
 
  shockData.DeleteAthenaArray();
  return loc;
}

//====================================================================================
// "Logarithmic" (power-law) mesh
// Note that the grid setup asks for the **local** x1min, x1max etc. 
// x is the "logical" position in the grid, with the logical grid runing
// from 0 to 1, i.e. x = i/Nx
Real LogMeshSpacingX1(Real x, RegionSize rs) {
  Real xf, xrat;
  xrat   = pow(rs.x1max/rs.x1min,1.0/((Real) rs.nx1)); // Only valid for fixed grid, no MPI
  xf     = rs.x1min*pow(xrat,x*rs.nx1); // x = i/Nx
  return xf;
}

//Global Variables for OuterX1
Real outerDens;
Real outerVel;
Real outerPres;
Real outerS;

void ExpandingOuterX1_UniformMedium(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
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
        if (NSCALARS ==1) {
          prim(IS0,k,j,ie+i) = outerS;
        }
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
// Reflecting inner X1 boundary conditions for radially non-uniform grids

void ReflectInnerX1_nonuniform(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  // copy hydro variables into ghost zones, reflecting v1
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVX)) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=1; i<=ngh; ++i) {
            prim(IVX,k,j,is-i) = -prim(IVX,k,j,(is+i-1));  // reflect 1-velocity
          }
        }
      }
    } else {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=1; i<=ngh; ++i) {
            prim(n,k,j,is-i) = prim(n,k,j,(is+i-1));
          }
        }
      }
    }
  }
  // copy face-centered magnetic fields into ghost zones, reflecting b1
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x1f(k,j,(is-i)) = -b.x1f(k,j,(is+i  ));  // reflect 1-field
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x2f(k,j,(is-i)) =  b.x2f(k,j,(is+i-1));
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x3f(k,j,(is-i)) =  b.x3f(k,j,(is+i-1));
        }
      }
    }
  }
  return;
}

//====================================================================================
// Enroll user-specific functions
void Mesh::InitUserMeshData(ParameterInput *pin) {
  Real x1rat = pin->GetOrAddReal("mesh","x1rat",0.0);
  
  if (x1rat < 0.0) {
    EnrollUserMeshGenerator(X1DIR, LogMeshSpacingX1);
  }
  if (mesh_bcs[INNER_X1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(INNER_X1, ReflectInnerX1_nonuniform);
  }
  
  if (EXPANDING) {
    SetGridData(4);
    EnrollGridDiffEq(WallVel);
    EnrollCalcGridData(UpdateGridData);
    
    if (mesh_bcs[OUTER_X1] == GetBoundaryFlag("user")) {
      EnrollUserBoundaryFunction(OUTER_X1,ExpandingOuterX1_UniformMedium);
    }
    //EnrollUserTimeStepFunction(ExpGridTimeStep);

    GridData(0) = pin->GetReal("mesh","x1min");
    GridData(1) = pin->GetReal("problem","rGrid");
    GridData(2) = pin->GetReal("problem","vGrid");
    GridData(3) = pin->GetReal("problem","comWidth");
    
  }

  return;
}

//Make sure grid does not move too far
//========================================================================================
//! \fn void MeshBlock::InitOTFOutput(ParameterInput *pin)
//  \brief Sets data structures etc for on-the-fly analysis.
//========================================================================================

void MeshBlock::InitOTFOutput(ParameterInput *pin) {
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Spherical blast wave test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real pi         = 4.0*std::atan(1.0);
  Real rej       = pin->GetReal("problem","rej"); // radius of initial shock
  Real dr         = pin->GetReal("problem","dr"); // width of transition
  Real vej       = pin->GetReal("problem","vej"); // velocity of ejecta
  Real dej       = pin->GetReal("problem","dej"); // mass of ejecta
  Real pej       = pin->GetReal("problem","pej"); // pressure of ejecta
  Real x1min      = pin->GetReal("mesh","x1min");

  Real p0         = pin->GetReal("problem","p0"); // ambient pressure
  Real d0         = pin->GetReal("problem","d0"); // ambient density

 
  outerDens = d0;
  outerPres = p0;
  outerVel  = 0.0;
  outerS = 1.0;

  Real gamma   = peos->GetGamma();
  Real gm1     = gamma - 1.0;

  //Real Tej     = gm1*Eej/mej; 
  Real eej     = pej/gm1;
  Real e0      = p0/gm1;

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
    std::stringstream msg;
    msg << "### FATAL ERROR in knovae.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM= " << COORDINATE_SYSTEM << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // setup uniform ambient medium with spherical over-pressured region
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real rad,thet,x;
        if (COORDINATE_SYSTEM == "cartesian") {
          x   = pcoord->x1v(i);
          Real y   = pcoord->x2v(j);
          Real z   = pcoord->x3v(k);
          rad      = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
          thet   = std::acos(z/rad);
        } else if (COORDINATE_SYSTEM == "cylindrical") {
          x = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
          Real z = pcoord->x3v(k);
          rad    = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
          thet   = std::acos(z/rad);
        } else { // if (COORDINATE_SYSTEM == "spherical_polar")
          x = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
          Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
          Real z = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          rad    = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
          thet   = pcoord->x2v(j);
        }
        
        //Real ejecta = (rej - x)/dr +0.5;
        Real ejecta =  0.5 * (1.0-std::tanh((rad-rej)/dr));
        //Real ambient =1.0;
        Real ambient = 0.5 * (1.0-std::tanh((rej-rad)/dr));       
        //if (ejecta >= 1.0) {
        //  ejecta = 1.0;
        //} else if (ejecta <=0.0) {
        //  ejecta = 0.0;
        //}
        phydro->u(IDN,k,j,i) = dej*ejecta + d0*ambient; 
        phydro->u(IM1,k,j,i) = dej*vej*ejecta; 
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = e0*ambient+ (0.5*SQR(vej)*dej + eej)*ejecta; 
        }
        if (DUAL_ENERGY) {
         // phydro->u(IIE,k,j,i) = e0;
        }
        if (NSCALARS == 1) {
          //if (rad >.95* r0 && rad <1.04* r0) phydro->u(NHYDRO-NSCALARS, k,j,i) = 1.0 ;
          phydro->u(IS0,k,j,i) = (2.0*ejecta+outerS*ambient)*phydro->u(IDN,k,j,i);
        }
        if (EXPANDING) {
          //std::cout << "dx1 = " << pex->delx1f(i) << std::endl;
          //phydro->u(IM1,k,j,i) += 0.5*(pex->delx1f(i)+pex->delx1f(i+1))*phydro->u(IDN,k,j,i);
        }
 

      }
    }
  }
  return;
}

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief otf diagnostics (shell position, sphericity etc)
//========================================================================================

void MeshBlock::UserWorkInLoop(void) {
  return;
}

//========================================================================================
//! \fn void MeshBlock::OTFWorkBeforeOutput(void)
//  \brief otf diagnostics (shell position, sphericity etc)
//========================================================================================

void MeshBlock::OTFWorkBeforeOutput(ParameterInput *pin) {
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief Any clean-up etc.
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  return;
}
