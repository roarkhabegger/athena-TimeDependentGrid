//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file expansion.cpp
//  \brief implementation of functions in class Expansion

// C/C++ headers
#include <string>
#include <algorithm>  // min()
#include <cfloat>     // FLT_MAX
#include <cmath>      // fabs(), sqrt()

// Athena++ headers
#include "expansion.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../reconstruct/reconstruction.hpp"
// constructor, initializes data structures and parameters

Expansion::Expansion(MeshBlock *pmb, ParameterInput *pin) {
  bool coarse_flag=pmb->pcoord->coarse_flag;
  pmy_block = pmb;
  if (coarse_flag==true) {
    is = pmb->cis; js = pmb->cjs; ks = pmb->cks;
    ie = pmb->cie; je = pmb->cje; ke = pmb->cke;
    ng=pmb->cnghost;
  } else {
    is = pmb->is; js = pmb->js; ks = pmb->ks;
    ie = pmb->ie; je = pmb->je; ke = pmb->ke;
    ng=NGHOST;
  }
  // Allocate memory for mesh dimensions
  int ncells1 = (ie-is+1) + 2*ng;
  il = is-ng;
  iu = ie+ng;
  int ncells2 = 1, ncells3 = 1;
  jl = js;
  ju = je;
  kl = ks;
  ku = ke;
  if (pmb->block_size.nx2 > 1) {ncells2 = (je-js+1) + 2*ng; jl = js-ng; ju = je+ng;}
  if (pmb->block_size.nx3 > 1) {ncells3 = (ke-ks+1) + 2*ng; kl = ks-ng; ku = ke+ng;}

   //Allocate velocities and grid data
  v1f.NewAthenaArray((ncells1+1));
  v2f.NewAthenaArray((ncells2+1));
  v3f.NewAthenaArray((ncells3+1));
  x1_0.NewAthenaArray((ncells1+1));
  x2_0.NewAthenaArray((ncells2+1));
  x3_0.NewAthenaArray((ncells3+1));
  x1_1.NewAthenaArray((ncells1+1));
  x2_1.NewAthenaArray((ncells2+1));
  x3_1.NewAthenaArray((ncells3+1));
  x1_2.NewAthenaArray((ncells1+1));
  x2_2.NewAthenaArray((ncells2+1));
  x3_2.NewAthenaArray((ncells3+1));
  vol.NewAthenaArray(ncells3,ncells2,ncells1);
 
  Expwl.NewAthenaArray((NWAVE+NINT+NSCALARS),ncells3,ncells2,ncells1);
  Expwr.NewAthenaArray((NWAVE+NINT+NSCALARS),ncells3,ncells2,ncells1);
  expFlux[X1DIR].NewAthenaArray(NHYDRO,ncells3,ncells2,ncells1+1);
  if (pmy_block->block_size.nx2 > 1)
    expFlux[X2DIR].NewAthenaArray(NHYDRO,ncells3,ncells2+1,ncells1);
  if (pmy_block->block_size.nx3 > 1)
    expFlux[X3DIR].NewAthenaArray(NHYDRO,ncells3+1,ncells2,ncells1);
  for (int i=il; i<=iu+1;++i) {
    v1f(i) = 0.0;
    x1_0(i) = pmb->pcoord->x1f(i);
 }

  for (int j=jl; j<=ju+1;++j) {
    v2f(j) = 0.0;
    x2_0(j) = pmb->pcoord->x2f(j);
  }
  
  for (int k=kl; k<=ku+1;++k) {
    v3f(k) = 0.0;
    x3_0(k) = pmb->pcoord->x3f(k);
  }
  
  for (int k=kl; k<=ku;++k) {
    for (int j=jl; j<=ju;++j) {
      for (int i=il; i<=iu;++i) {
        vol(k,j,i) = pmb->pcoord->GetCellVolume(k,j,i);
      }
    }
  }
  mydt = (FLT_MAX);
}
// destructor

Expansion::~Expansion() {
  v1f.DeleteAthenaArray();
  v2f.DeleteAthenaArray();
  v3f.DeleteAthenaArray(); 

  x1_0.DeleteAthenaArray();
  x2_0.DeleteAthenaArray();
  x3_0.DeleteAthenaArray();
  x1_1.DeleteAthenaArray();
  x2_1.DeleteAthenaArray();
  x3_1.DeleteAthenaArray();
  x1_2.DeleteAthenaArray();
  x2_2.DeleteAthenaArray();
  x3_2.DeleteAthenaArray();

  Expwl.DeleteAthenaArray();
  Expwr.DeleteAthenaArray();
  expFlux[X1DIR].DeleteAthenaArray();
  if (pmy_block->block_size.nx2 > 1)
    expFlux[X2DIR].DeleteAthenaArray();
  if (pmy_block->block_size.nx3 > 1)
    expFlux[X3DIR].DeleteAthenaArray();
}

void Expansion::WeightedAveX(const int low, const int up, AthenaArray<Real> &x_out, AthenaArray<Real> &x_in1, AthenaArray<Real> &x_in2, const Real wght[3]){

  if (wght[2] != 0.0) {
#pragma omp simd
    for (int i=low; i<=up; ++i) {
      x_out(i) = wght[0]*x_out(i) + wght[1]*x_in1(i) + wght[2]*x_in2(i);
    }     
  } else { // do not dereference u_in2
    if (wght[1] != 0.0) {
#pragma omp simd
      for (int i=low; i<=up; ++i) {
        x_out(i) = wght[0]*x_out(i) + wght[1]*x_in1(i);
      }     
    } else { // do not dereference u_in1
      if (wght[0] != 0.0) {
#pragma omp simd
        for (int i=low; i<=up; ++i) {
          x_out(i) = wght[0]*x_out(i);
        }
      } else { // directly initialize u_out to 0
#pragma omp simd
        for (int i=low; i<=up; ++i) {
          x_out(i) = 0.0;
        }
      }
    }
  }
  return;
}


void Expansion::IntegrateWalls(Real dt){

  for (int i=il; i<=iu+1; ++i) {
    x1_0(i) += dt*v1f(i);
  }
  for (int j=jl; j<=ju+1; ++j) {
    x2_0(j) += dt*v2f(j);
  }
  for (int k=kl; k<=ku+1; ++k) {
    x3_0(k) += dt*v3f(k);
  }
  //Save current volumes
  Coordinates *pc = pmy_block->pcoord;
  for (int k=kl; k<=ku;++k) {
    for (int j=jl; j<=ju;++j) {
      for (int i=il; i<=iu;++i) {
        Real volume;
        Real dx1, dx2, dx3;
        dx1 = pc->dx1f(i)+v1f(i+1)*dt - v1f(i)*dt ;          
        dx2 = pc->dx2f(j)+v2f(j+1)*dt - v2f(j)*dt ;          
        dx3 = pc->dx3f(k)+v3f(k+1)*dt - v3f(k)*dt ;          

        if (COORDINATE_SYSTEM == "cartesian") {
          volume = dx1*dx2*dx3;          
        } else if (COORDINATE_SYSTEM == "cylindrical") {


        } else if (COORDINATE_SYSTEM == "spherical_polar") {


        } else {
          std::cout << "Unsupported cordinate system in IntegrateWalls. " << std::endl;
        }
        vol(k,j,i) = volume;
      }
    }
  }

  return;
}
//Source Term Funciton
void Expansion::AddWallFluxDivergence(Real dt, AthenaArray<Real> &prim, AthenaArray<Real> &cons){
  Real vol,A1,A2 = 0.0;
  Real flxL,flxR = 0.0;
  Real divF = 0.0;
  Real qLi,qRi,qLp1,qRp1 = 0.0;
  Real gm1 = pmy_block->peos->GetGamma()-1.0;
  int order = 3;
  MeshBlock *pmb = pmy_block;
  AthenaArray<Real> &x1flux=expFlux[X1DIR];
  AthenaArray<Real> &x2flux=expFlux[X2DIR];
  AthenaArray<Real> &x3flux=expFlux[X3DIR];
  AthenaArray<Real> b1,b2,b3,w_x1f,w_x2f,w_x3f,e2x1,e3x1,e1x2,e3x2,e1x3,e2x3;
  if (MAGNETIC_FIELDS_ENABLED) {
    b1.InitWithShallowCopy(pmb->pfield->b.x1f);
    b2.InitWithShallowCopy(pmb->pfield->b.x2f);
    b3.InitWithShallowCopy(pmb->pfield->b.x3f);
    w_x1f.InitWithShallowCopy(pmb->pfield->wght.x1f);
    w_x2f.InitWithShallowCopy(pmb->pfield->wght.x2f);
    w_x3f.InitWithShallowCopy(pmb->pfield->wght.x3f);
    e2x1.InitWithShallowCopy(pmb->pfield->e2_x1f);
    e3x1.InitWithShallowCopy(pmb->pfield->e3_x1f);
    e1x2.InitWithShallowCopy(pmb->pfield->e1_x2f);
    e3x2.InitWithShallowCopy(pmb->pfield->e3_x2f);
    e1x3.InitWithShallowCopy(pmb->pfield->e1_x3f);
    e2x3.InitWithShallowCopy(pmb->pfield->e2_x3f);
  }

  if (order ==3) {
    PiecewiseParabolicOffsetX1(pmb,ks,ke,js,je,is,ie+1,prim,pmb->pfield->bcc,Expwl,Expwr,dt);
  } else {
    PiecewiseLinearOffsetX1(pmb,ks,ke,js,je,is,ie+1,prim,pmb->pfield->bcc,Expwl,Expwr,dt);
  }

  FluxSolver(ks,ke,js,je,is,ie+1,IVX,b1,Expwl,Expwr,x1flux,e3x1,e2x1,v1f);


  for (int k = ks; k<=ke;++k) {
    for (int j = js; j<=je;++j) {      
      for (int i = is; i<=ie;++i) {
        A1 = pmb->pcoord->GetFace1Area(k,j,i);
        A2 = pmb->pcoord->GetFace1Area(k,j,i+1);
        vol = pmb->pcoord->GetCellVolume(k,j,i);
        for (int n=0; n<NHYDRO;++n) {
          divF = x1flux(n,k,j,i+1)*A2 - x1flux(n,k,j,i) * A1;
          cons(n,k,j,i) += dt/vol*divF;        
        }
      }
    }
  }
  return;
}

void Expansion::ExpansionSourceTerms(const Real dt, const AthenaArray<Real> *flux, const AthenaArray<Real> &prim, AthenaArray<Real> &cons) {
  Real vm = 0.0;
  Coordinates *pc = pmy_block->pcoord;
  for (int k = ks; k<=ke;++k) {
    for (int j = js; j<=je;++j) {      
      for (int i = is; i<=ie;++i) {
        for (int n = 0; n<(NHYDRO+NSCALARS); ++n) {
          Real newVol = vol(k,j,i);
          Real oldVol = pc->GetCellVolume(k,j,i);
          //std::cout << oldVol-newVol << std::endl;
          cons(n,k,j,i) *= oldVol/newVol;        
        }
      }
    }
  }

  return;
}


void Expansion::UpdateVelData(MeshBlock *pmb ,Real time, Real dt){
  mydt = dt;
  if (pmb->pmy_mesh->GridDiffEq_ != NULL) {
    //Edit each delx1f, delx2f, delx3f before source terms and editing of grid      
    for (int k = kl; k<=ku+1;++k){
      v3f(k) = pmb->pmy_mesh->GridDiffEq_(pmb->pcoord->x3f(k),k,pmb->pmy_mesh->time,dt,2,pmb->pmy_mesh->GridData);
    }
    for (int j = jl;j<=ju+1;++j){      
      v2f(j) = pmb->pmy_mesh->GridDiffEq_(pmb->pcoord->x2f(j),j,pmb->pmy_mesh->time,dt,1,pmb->pmy_mesh->GridData);
    }
    for (int i = il; i<=iu+1;++i){
      v1f(i) = pmb->pmy_mesh->GridDiffEq_(pmb->pcoord->x1f(i),i,pmb->pmy_mesh->time,dt,0,pmb->pmy_mesh->GridData); 
    }
  } 
  return;
}

void Expansion::GridEdit(MeshBlock *pmb){  
  //FACE CENTERED
  //x1
  for (int i=il; i<=iu+1; ++i){
    pmb->pcoord->x1f(i) = x1_0(i);  
  }
  for (int i=il; i<=iu; ++i) {
    pmb->pcoord->dx1f(i) = pmb->pcoord->x1f(i+1)-pmb->pcoord->x1f(i);
  }
  
  //x2  
  if (pmb->block_size.nx2 !=1){
    for (int j=jl; j<=ju+1; ++j){
      pmb->pcoord->x2f(j) = x2_0(j);
    }
    for (int j=jl; j<=ju; ++j) {
      pmb->pcoord->dx2f(j) = pmb->pcoord->x2f(j+1)- pmb->pcoord->x2f(j);
    }
  }                      
  //x3
  if (pmb->block_size.nx3 !=1){
    for (int k=kl; k<=ku+1; ++k){
      pmb->pcoord->x3f(k) = x3_0(k);
    }
    for (int k=kl; k<=ku; ++k) {
      pmb->pcoord->dx3f(k) = pmb->pcoord->x3f(k+1)- pmb->pcoord->x3f(k);
    }
  }

  //Set Reconstruction Coefficients
  for (int i=il+1; i<=iu-1; ++i) {
    Real& dx_im1 = pmb->pcoord->dx1f(i-1);
    Real& dx_i   = pmb->pcoord->dx1f(i  );
    Real& dx_ip1 = pmb->pcoord->dx1f(i+1);
    Real qe = dx_i/(dx_im1 + dx_i + dx_ip1);       // Outermost coeff in CW eq 1.7
    pmb->precon->c1i(i) = qe*(2.0*dx_im1+dx_i)/(dx_ip1 + dx_i); // First term in CW eq 1.7
    pmb->precon->c2i(i) = qe*(2.0*dx_ip1+dx_i)/(dx_im1 + dx_i); // Second term in CW eq 1.7
    if (i > il+1) {  // c3-c6 are not computed in first iteration
      Real& dx_im2 = pmb->pcoord->dx1f(i-2);
      Real qa = dx_im2 + dx_im1 + dx_i + dx_ip1;
      Real qb = dx_im1/(dx_im1 + dx_i);
      Real qc = (dx_im2 + dx_im1)/(2.0*dx_im1 + dx_i);
      Real qd = (dx_ip1 + dx_i)/(2.0*dx_i + dx_im1);
      qb = qb + 2.0*dx_i*qb/qa*(qc-qd);
      pmb->precon->c3i(i) = 1.0 - qb;
      pmb->precon->c4i(i) = qb;
      pmb->precon->c5i(i) = dx_i/qa*qd;
      pmb->precon->c6i(i) = -dx_im1/qa*qc;
    }
  }
  // Compute curvilinear geometric factors for limiter (Mignone eq 48)
  for (int i=il+1; i<=iu-1; ++i) {
    if ((COORDINATE_SYSTEM == "cylindrical") ||
        (COORDINATE_SYSTEM == "spherical_polar")) {
      Real h_plus, h_minus;
      Real& dx_i   = pmb->pcoord->dx1f(i);
      Real& xv_i   = pmb->pcoord->x1v(i);
      if (COORDINATE_SYSTEM == "cylindrical") {
        // cylindrical radial coordinate
        h_plus = 3.0 + dx_i/(2.0*xv_i);
        h_minus = 3.0 - dx_i/(2.0*xv_i);
      } else {
        // spherical radial coordinate
        h_plus = 3.0 + (2.0*dx_i*(10.0*xv_i + dx_i))/(20.0*SQR(xv_i) + SQR(dx_i));
        h_minus = 3.0 + (2.0*dx_i*(-10.0*xv_i + dx_i))/(20.0*SQR(xv_i) + SQR(dx_i));
      }
      pmb->precon->hplus_ratio_i(i) = (h_plus + 1.0)/(h_minus - 1.0);
      pmb->precon->hminus_ratio_i(i) = (h_minus + 1.0)/(h_plus - 1.0);
    } else { // Cartesian, SR, GR
        // Ratios are = 2 for Cartesian coords, as in original PPM overshoot limiter
        pmb->precon->hplus_ratio_i(i) = 2.0;
        pmb->precon->hminus_ratio_i(i) = 2.0;
    }
  }
  //VOLUME BASED QUANTITIES 
  if (COORDINATE_SYSTEM == "cartesian") {
    //Cartesian
    // initialize volume-averaged coordinates and spacing
    // x1-direction: x1v = dx/2
    for (int i=il; i<=iu; ++i) {
      pmb->pcoord->x1v(i) = 0.5*(pmb->pcoord->x1f(i+1) + pmb->pcoord->x1f(i));
    }
    for (int i=il; i<=iu-1; ++i) {
      if (pmb->block_size.x1rat != 1.0) {
        pmb->pcoord->dx1v(i) = pmb->pcoord->x1v(i+1) - pmb->pcoord->x1v(i);
      } else {
        // dx1v = dx1f constant for uniform mesh; may disagree with x1v(i+1) - x1v(i)
        pmb->pcoord->dx1v(i) = pmb->pcoord->dx1f(i);
      }
    }
  
    // x2-direction: x2v = dy/2
    if (pmb->block_size.nx2 == 1) {
      pmb->pcoord->x2v(jl) = 0.5*(pmb->pcoord->x2f(jl+1) + pmb->pcoord->x2f(jl));
      pmb->pcoord->dx2v(jl) = pmb->pcoord->dx2f(jl);
    } else {
      for (int j=jl; j<=ju; ++j) {
        pmb->pcoord->x2v(j) = 0.5*(pmb->pcoord->x2f(j+1) + pmb->pcoord->x2f(j));
      }
      for (int j=jl; j<=ju-1; ++j) {
        if (pmb->block_size.x2rat != 1.0) {
          pmb->pcoord->dx2v(j) = pmb->pcoord->x2v(j+1) - pmb->pcoord->x2v(j);
        } else {
          // dx2v = dx2f constant for uniform mesh; may disagree with x2v(j+1) - x2v(j)
          pmb->pcoord->dx2v(j) = pmb->pcoord->dx2f(j);
        }
      }
    }
  
    // x3-direction: x3v = dz/2
    if (pmb->block_size.nx3 == 1) {
      pmb->pcoord->x3v(kl) = 0.5*(pmb->pcoord->x3f(kl+1) + pmb->pcoord->x3f(kl));
      pmb->pcoord->dx3v(kl) = pmb->pcoord->dx3f(kl);
    } else {
      for (int k=kl; k<=ku; ++k) {
        pmb->pcoord->x3v(k) = 0.5*(pmb->pcoord->x3f(k+1) + pmb->pcoord->x3f(k));
      }
      for (int k=kl; k<=ku-1; ++k) {
        if (pmb->block_size.x3rat != 1.0) {
          pmb->pcoord->dx3v(k) = pmb->pcoord->x3v(k+1) - pmb->pcoord->x3v(k);
        } else {
          // dxkv = dx3f constant for uniform mesh; may disagree with x3v(k+1) - x3v(k)
          pmb->pcoord->dx3v(k) = pmb->pcoord->dx3f(k);
        }
      }
    }


    // initialize area-averaged coordinates used with MHD AMR
    if ((pmb->pmy_mesh->multilevel==true) && MAGNETIC_FIELDS_ENABLED) {
      for (int i=il; i<=iu; ++i) {
        pmb->pcoord->x1s2(i) = pmb->pcoord->x1s3(i) = pmb->pcoord->x1v(i);
      }
      if (pmb->block_size.nx2 == 1) {
        pmb->pcoord->x2s1(jl) = pmb->pcoord->x2s3(jl) = pmb->pcoord->x2v(jl);
      } else {
        for (int j=jl; j<=ju; ++j) {
          pmb->pcoord->x2s1(j) = pmb->pcoord->x2s3(j) = pmb->pcoord->x2v(j);
        }
      }
      if (pmb->block_size.nx3 == 1) {
        pmb->pcoord->x3s1(kl) = pmb->pcoord->x3s2(kl) = pmb->pcoord->x3v(kl);
      } else {
        for (int k=kl; k<=ku; ++k) {
          pmb->pcoord->x3s1(k) = pmb->pcoord->x3s2(k) = pmb->pcoord->x3v(k);
        }
      }
    }

    //Reset Reconstruction coefficients.


  } else if (COORDINATE_SYSTEM == "cylindrical") { 
    //Cylindrical
    for (int i=il; i<=iu; ++i) {
      pmb->pcoord->x1v(i) = (TWO_3RD)*(pow(pmb->pcoord->x1f(i+1),3)-pow(pmb->pcoord->x1f(i),3))/(pow(pmb->pcoord->x1f(i+1),2) - pow(pmb->pcoord->x1f(i),2));
    }
    for (int i=il; i<=iu-1; ++i) {
      pmb->pcoord->dx1v(i) = pmb->pcoord->x1v(i+1) - pmb->pcoord->x1v(i);
    }
    
    // x2-direction: x2v = (\int phi dV / \int dV) = dphi/2
    if (pmb->block_size.nx2 == 1) { 
      pmb->pcoord->x2v(jl) = 0.5*(pmb->pcoord->x2f(jl+1) + pmb->pcoord->x2f(jl));
      pmb->pcoord->dx2v(jl) = pmb->pcoord->dx2f(jl);
    } else {
      for (int j=jl; j<=ju; ++j) {
        pmb->pcoord->x2v(j) = 0.5*(pmb->pcoord->x2f(j+1) + pmb->pcoord->x2f(j));
      }
      for (int j=jl; j<=ju-1; ++j) {
        pmb->pcoord->dx2v(j) = pmb->pcoord->x2v(j+1) - pmb->pcoord->x2v(j);
      }
    }
    
    // x3-direction: x3v = (\int z dV / \int dV) = dz/2
    if (pmb->block_size.nx3 == 1) {
      pmb->pcoord->x3v(kl) = 0.5*(pmb->pcoord->x3f(kl+1) + pmb->pcoord->x3f(kl));
      pmb->pcoord->dx3v(kl) = pmb->pcoord->dx3f(kl);
    } else {
      for (int k=kl; k<=ku; ++k) {
        pmb->pcoord->x3v(k) = 0.5*(pmb->pcoord->x3f(k+1) + pmb->pcoord->x3f(k));
      }
      for (int k=kl; k<=ku-1; ++k) {
        pmb->pcoord->dx3v(k) = pmb->pcoord->x3v(k+1) - pmb->pcoord->x3v(k);
      }
    }
     
    //Geometry Coefficients
    for (int i=il; i<=iu; ++i){
      pmb->pcoord->h2v(i) = pmb->pcoord->x1v(i);
      pmb->pcoord->h2f(i) = pmb->pcoord->x1f(i); 
    }
    
    //Area averaged coordinates
    if ((pmb->pmy_mesh->multilevel==true) && MAGNETIC_FIELDS_ENABLED) {
      for (int i=il; i<=iu; ++i) {
        pmb->pcoord->x1s2(i) = pmb->pcoord->x1s3(i) = pmb->pcoord->x1v(i);
      }
      if (pmb->block_size.nx2 == 1) {
        pmb->pcoord->x2s1(jl) = pmb->pcoord->x2s3(jl) = pmb->pcoord->x2v(jl);
      } else {
        for (int j=jl; j<=ju; ++j) {
          pmb->pcoord->x2s1(j) = pmb->pcoord->x2s3(j) = pmb->pcoord->x2v(j);
        }
      }
      if (pmb->block_size.nx3 == 1) {
        pmb->pcoord->x3s1(kl) = pmb->pcoord->x3s2(kl) = pmb->pcoord->x3v(kl);
      } else {
        for (int k=kl; k<=ku; ++k) {
          pmb->pcoord->x3s1(k) = pmb->pcoord->x3s2(k) = pmb->pcoord->x3v(k);
        }
      }
    }

    //Edit Scratch Arrays with coefficients
    if (pmb->pcoord->coarse_flag==false) {
      // Compute and store constant coefficients needed for face-areas, cell-volumes, etc.
      // This helps improve performance
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        Real rm = pmb->pcoord->x1f(i  );
        Real rp = pmb->pcoord->x1f(i+1);
        // dV = 0.5*(R_{i+1}^2 - R_{i}^2)
        pmb->pcoord->coord_area3_i_(i)= 0.5*(rp*rp - rm*rm);
        // dV = 0.5*(R_{i+1}^2 - R_{i}^2)
        pmb->pcoord->coord_vol_i_(i) = pmb->pcoord->coord_area3_i_(i);
        // (A1^{+} - A1^{-})/dV
        pmb->pcoord->coord_src1_i_(i) = pmb->pcoord->dx1f(i)/pmb->pcoord->coord_vol_i_(i);
        // (dR/2)/(R_c dV)
        pmb->pcoord->coord_src2_i_(i) = pmb->pcoord->dx1f(i)/((rm + rp)*pmb->pcoord->coord_vol_i_(i));
        // Rf_{i}/R_{i}/Rf_{i}^2
        pmb->pcoord->phy_src1_i_(i) = 1.0/(pmb->pcoord->x1v(i)*pmb->pcoord->x1f(i));
      }
#pragma omp simd
      for (int i=il; i<=iu-1; ++i) {
        // Rf_{i+1}/R_{i}/Rf_{i+1}^2
        pmb->pcoord->phy_src2_i_(i) = 1.0/(pmb->pcoord->x1v(i)*pmb->pcoord->x1f(i+1));
        // dV = 0.5*(R_{i+1}^2 - R_{i}^2)
        pmb->pcoord->coord_area3vc_i_(i)= 0.5*(SQR(pmb->pcoord->x1v(i+1)) - SQR(pmb->pcoord->x1v(i)));
      }
    }

  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    //Spherical
    //x1 deltas, volumes  
    // x1-direction: x1v = (\int r dV / \int dV) = d(r^4/4)/d(r^3/3)
    for (int i=il; i<iu; ++i) {
      pmb->pcoord->x1v(i) = 0.75*(pow(pmb->pcoord->x1f(i+1),4) - pow(pmb->pcoord->x1f(i),4))/(pow(pmb->pcoord->x1f(i+1),3) - pow(pmb->pcoord->x1f(i),3));
    }
    for (int i=il; i<=iu-1; ++i) {
      pmb->pcoord->dx1v(i) = pmb->pcoord->x1v(i+1) - pmb->pcoord->x1v(i);
    }
  
    //x2 Deltas and volumes
    // x2-direction: x2v = (\int sin[theta] theta dV / \int dV) =
    //  d(sin[theta] - theta cos[theta])/d(-cos[theta])
    if (pmb->block_size.nx2 == 1) {
      pmb->pcoord->x2v(jl) = 0.5*(pmb->pcoord->x2f(jl+1) + pmb->pcoord->x2f(jl));
      pmb->pcoord->dx2v(jl) = pmb->pcoord->dx2f(jl);
    } else {
      for (int j=jl; j<=ju; ++j) {
        pmb->pcoord->x2v(j) = ((sin(pmb->pcoord->x2f(j+1)) - pmb->pcoord->x2f(j+1)*cos(pmb->pcoord->x2f(j+1))) -
                  (sin(pmb->pcoord->x2f(j  )) - pmb->pcoord->x2f(j  )*cos(pmb->pcoord->x2f(j  ))))/
                  (cos(pmb->pcoord->x2f(j  )) - cos(pmb->pcoord->x2f(j+1)));
      }
      for (int j=jl; j<=ju-1; ++j) {
        pmb->pcoord->dx2v(j) = pmb->pcoord->x2v(j+1) - pmb->pcoord->x2v(j);
      }
    }
      
    //x3 Deltas and volumes
    // x3-direction: x3v = (\int phi dV / \int dV) = dphi/2 
    if (pmb->block_size.nx3 == 1) {
      pmb->pcoord->x3v(kl) = 0.5*(pmb->pcoord->x3f(kl+1) + pmb->pcoord->x3f(kl));
      pmb->pcoord->dx3v(kl) = pmb->pcoord->dx3f(kl);
    } else {
      for (int k=kl; k<=ku; ++k) {
        pmb->pcoord->x3v(k) = 0.5*(pmb->pcoord->x3f(k+1) + pmb->pcoord->x3f(k));
      }
      for (int k=kl; k<=ku-1; ++k) {
        pmb->pcoord->dx3v(k) = pmb->pcoord->x3v(k+1) - pmb->pcoord->x3v(k);
      }
    }

    //Geometry Coefficients
    for (int i=il; i<iu; ++i) {
      pmb->pcoord->h2v(i) = pmb->pcoord->x1v(i);
      pmb->pcoord->h2f(i) = pmb->pcoord->x1f(i);
      pmb->pcoord->h31v(i) = pmb->pcoord->x1v(i);
      pmb->pcoord->h31f(i) = pmb->pcoord->x1f(i);
    }
    // x2-direction
    if (pmb->block_size.nx2 == 1) {
      pmb->pcoord->h32v(jl) = sin(pmb->pcoord->x2v(jl));
      pmb->pcoord->h32f(jl) = sin(pmb->pcoord->x2f(jl));
      pmb->pcoord->dh32vd2(jl) = cos(pmb->pcoord->x2v(jl));
      pmb->pcoord->dh32fd2(jl) = cos(pmb->pcoord->x2f(jl));
    } else {
      for (int j=jl; j<=ju; ++j) {
        pmb->pcoord->h32v(j) = sin(pmb->pcoord->x2v(j));
        pmb->pcoord->h32f(j) = sin(pmb->pcoord->x2f(j));
        pmb->pcoord->dh32vd2(j) = cos(pmb->pcoord->x2v(j));
        pmb->pcoord->dh32fd2(j) = cos(pmb->pcoord->x2f(j));
      }
    }
    if ((pmb->pmy_mesh->multilevel==true) && MAGNETIC_FIELDS_ENABLED) {
      for (int i=il; i<=iu; ++i) {
        pmb->pcoord->x1s2(i) = pmb->pcoord->x1s3(i) = (2.0/3.0)*(pow(pmb->pcoord->x1f(i+1),3) - pow(pmb->pcoord->x1f(i),3))
                            /(SQR(pmb->pcoord->x1f(i+1)) - SQR(pmb->pcoord->x1f(i)));
      }
      if (pmb->block_size.nx2 == 1) {
        pmb->pcoord->x2s1(jl) = pmb->pcoord->x2s3(jl) = pmb->pcoord->x2v(jl);
      } else {
        for (int j=jl; j<=ju; ++j) {
          pmb->pcoord->x2s1(j) = pmb->pcoord->x2s3(j) = pmb->pcoord->x2v(j);
        }
      }
      if (pmb->block_size.nx3 == 1) {
        pmb->pcoord->x3s1(kl) = pmb->pcoord->x3s2(kl) = pmb->pcoord->x3v(kl);
      } else {
        for (int k=kl; k<=ku; ++k) {
          pmb->pcoord->x3s1(k) = pmb->pcoord->x3s2(k) = pmb->pcoord->x3v(k);
        }
      }
    }

    if (pmb->pcoord->coarse_flag==false){
#pragma omp simd
     //Fix internal scratch arrays for good measure
     for (int i=il; i<=iu; ++i) {
         Real rm = pmb->pcoord->x1f(i  );
         Real rp = pmb->pcoord->x1f(i+1);
         // R^2
         pmb->pcoord->coord_area1_i_(i) = rm*rm;
         // 0.5*(R_{i+1}^2 - R_{i}^2)
         pmb->pcoord->coord_area2_i_(i) = 0.5*(rp*rp - rm*rm);
         // 0.5*(R_{i+1}^2 - R_{i}^2)
         pmb->pcoord->coord_area3_i_(i) = pmb->pcoord->coord_area2_i_(i);
         // dV = (R_{i+1}^3 - R_{i}^3)/3
         pmb->pcoord->coord_vol_i_(i) = (ONE_3RD)*(rp*rp*rp - rm*rm*rm);
         // (A1^{+} - A1^{-})/dV
         pmb->pcoord->coord_src1_i_(i) = pmb->pcoord->coord_area2_i_(i)/pmb->pcoord->coord_vol_i_(i);
         // (dR/2)/(R_c dV)
         pmb->pcoord->coord_src2_i_(i) = pmb->pcoord->dx1f(i)/((rm + rp)*pmb->pcoord->coord_vol_i_(i));
         // Rf_{i}^2/R_{i}^2/Rf_{i}^2
         pmb->pcoord->phy_src1_i_(i) = 1.0/SQR(pmb->pcoord->x1v(i));
         // Rf_{i+1}^2/R_{i}^2/Rf_{i+1}^2
         pmb->pcoord->phy_src2_i_(i) = pmb->pcoord->phy_src1_i_(i);
         // R^2 at the volume center for non-ideal MHD
         pmb->pcoord->coord_area1vc_i_(i) = SQR(pmb->pcoord->x1v(i));
     }
     pmb->pcoord->coord_area1_i_(iu+ng+1) = pmb->pcoord->x1f(iu+ng+1)*pmb->pcoord->x1f(iu+ng+1);
#pragma omp simd
      for (int i=il; i<=iu-1; ++i) {//non-ideal MHD
        // 0.5*(R_{i+1}^2 - R_{i}^2)
        pmb->pcoord->coord_area2vc_i_(i)= 0.5*(SQR(pmb->pcoord->x1v(i+1))-SQR(pmb->pcoord->x1v(i)));
        // 0.5*(R_{i+1}^2 - R_{i}^2)
        pmb->pcoord->coord_area3vc_i_(i)= pmb->pcoord->coord_area2vc_i_(i);
      }

      if (pmb->block_size.nx2 > 1) {
#pragma omp simd
        for (int j=jl; j<=ju; ++j) {
          Real sm = fabs(sin(pmb->pcoord->x2f(j  )));
          Real sp = fabs(sin(pmb->pcoord->x2f(j+1)));
          Real cm = cos(pmb->pcoord->x2f(j  ));
          Real cp = cos(pmb->pcoord->x2f(j+1));
          // d(sin theta) = d(-cos theta)
          pmb->pcoord->coord_area1_j_(j) = fabs(cm - cp);
          // sin theta
          pmb->pcoord->coord_area2_j_(j) = sm;
          // d(sin theta) = d(-cos theta)
          pmb->pcoord->coord_vol_j_(j) = pmb->pcoord->coord_area1_j_(j);
          // (A2^{+} - A2^{-})/dV
          pmb->pcoord->coord_src1_j_(j) = (sp - sm)/pmb->pcoord->coord_vol_j_(j);
          // (dS/2)/(S_c dV)
          pmb->pcoord->coord_src2_j_(j) = (sp - sm)/((sm + sp)*pmb->pcoord->coord_vol_j_(j));
          // < cot theta > = (|sin th_p| - |sin th_m|) / |cos th_m - cos th_p|
          pmb->pcoord->coord_src3_j_(j) = (sp - sm)/pmb->pcoord->coord_vol_j_(j);
          // d(sin theta) = d(-cos theta) at the volume center for non-ideal MHD
          pmb->pcoord->coord_area1vc_j_(j)= fabs(cos(pmb->pcoord->x2v(j))-cos(pmb->pcoord->x2v(j+1)));
          // sin theta at the volume center for non-ideal MHD
          pmb->pcoord->coord_area2vc_j_(j)= fabs(sin(pmb->pcoord->x2v(j)));
        }
        pmb->pcoord->coord_area2_j_(ju+ng+1) = fabs(sin(pmb->pcoord->x2f(ju+ng+1)));
        if (pmb->pcoord->IsPole(jl))   // inner polar boundary
          pmb->pcoord->coord_area1vc_j_(jl-1)= 2.0-cos(pmb->pcoord->x2v(jl-1))-cos(pmb->pcoord->x2v(jl));
        if (pmb->pcoord->IsPole(ju))   // outer polar boundary
          pmb->pcoord->coord_area1vc_j_(ju)  = 2.0+cos(pmb->pcoord->x2v(ju))+cos(pmb->pcoord->x2v(ju+1));
      } else {
        Real sm = fabs(sin(pmb->pcoord->x2f(jl  )));
        Real sp = fabs(sin(pmb->pcoord->x2f(jl+1)));
        Real cm = cos(pmb->pcoord->x2f(jl  ));
        Real cp = cos(pmb->pcoord->x2f(jl+1));
        pmb->pcoord->coord_area1_j_(jl) = fabs(cm - cp);
        pmb->pcoord->coord_area2_j_(jl) = sm;
        pmb->pcoord->coord_area1vc_j_(jl)= pmb->pcoord->coord_area1_j_(jl);
        pmb->pcoord->coord_area2vc_j_(jl)= sin(pmb->pcoord->x2v(jl));
        pmb->pcoord->coord_vol_j_(jl) = pmb->pcoord->coord_area1_j_(jl);
        pmb->pcoord->coord_src1_j_(jl) = (sp - sm)/pmb->pcoord->coord_vol_j_(jl);
        pmb->pcoord->coord_src2_j_(jl) = (sp - sm)/((sm + sp)*pmb->pcoord->coord_vol_j_(jl));
        pmb->pcoord->coord_src3_j_(jl) = (sp - sm)/pmb->pcoord->coord_vol_j_(jl);
        pmb->pcoord->coord_area2_j_(jl+1) = sp;
      }
    }
  }

  //Check Mesh_Size object and reset bounds if necessary
  //x1
  if (pmb->block_size.x1min == pmb->pmy_mesh->mesh_size.x1min) {
    pmb->pmy_mesh->mesh_size.x1min = x1_0(is);
  }  
  if (pmb->block_size.x1max == pmb->pmy_mesh->mesh_size.x1max) {
    pmb->pmy_mesh->mesh_size.x1max = x1_0(ie);
  } 
  //x2 
  if (pmb->block_size.x2min == pmb->pmy_mesh->mesh_size.x2min) {
    pmb->pmy_mesh->mesh_size.x2min = x2_0(js);
  }  
  if (pmb->block_size.x2max == pmb->pmy_mesh->mesh_size.x2max) {
    pmb->pmy_mesh->mesh_size.x2max = x2_0(je);
  }  
  //x3
  if (pmb->block_size.x3min == pmb->pmy_mesh->mesh_size.x3min) {
    pmb->pmy_mesh->mesh_size.x3min = x3_0(ks);
  }  
  if (pmb->block_size.x3max == pmb->pmy_mesh->mesh_size.x3max) {
    pmb->pmy_mesh->mesh_size.x3max = x3_0(ke);
  } 

  //Reset block_size object 
  pmb->block_size.x1min = x1_0(is);
  pmb->block_size.x2min = x2_0(js);
  pmb->block_size.x3min = x3_0(ks);
  
  pmb->block_size.x1max = x1_0(ie);
  pmb->block_size.x2max = x2_0(je);
  pmb->block_size.x3max = x3_0(ke);
   
  return;
}


Real Expansion::GridTimeStep(MeshBlock *pmb){ 
  
  Real nextPosDelta, minCellSize;//, dt;
  
  Mesh *pmesh = pmb->pmy_mesh;
  int is, ie, js, je, ks, ke, ng;
  is = pmb->is; js = pmb->js; ks = pmb->ks;
  ie = pmb->ie; je = pmb->je; ke = pmb->ke;
  ng = NGHOST;

  Real mydt = pmesh->dt;

  Real min_dt = mydt*100000;
  for (int i = is; i<=ie+1;i++){
    nextPosDelta = pmesh->GridDiffEq_(pmb->pcoord->x1f(i),i,pmesh->time,mydt,0,pmesh->GridData)*mydt;

    if (nextPosDelta < 0 && i != ie-ng){ 
      minCellSize = pmb->pcoord->dx1f(i-1);
    } else {
      minCellSize = pmb->pcoord->dx1f(i);
    }

    minCellSize *= pmesh->cfl_number;
    nextPosDelta = fabs(nextPosDelta);

    if (nextPosDelta != 0.0 && minCellSize != 0.0){
      Real overStep = minCellSize - nextPosDelta;
      int count = 0;
      while (overStep <= 0){
        mydt *= 0.9;
        nextPosDelta = pmesh->GridDiffEq_(pmb->pcoord->x1f(i),i,pmesh->time,mydt,0,pmesh->GridData)*mydt;
        overStep = minCellSize - nextPosDelta;      
        count++;
        if (count >= 100) {
          mydt*=0.5;
        }
        
      }
      Real dtEx = mydt;//*pmesh->cfl_number;//fabs(minCellSize* 0.5/nextPosDelta * (pmesh->dt));
      min_dt = std::min(min_dt, dtEx);
    } else {
      Real dtEx = pmesh->dt;
      min_dt = std::min(min_dt,dtEx);
    }
  }
  return min_dt;
}


