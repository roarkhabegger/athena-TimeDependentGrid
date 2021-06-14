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
#include "../globals.hpp"
#include "expansion.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../reconstruct/reconstruction.hpp"
// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif
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


  if (COORDINATE_SYSTEM == "cartesian") {
    x1Move = true;
    x2Move = true;
    x3Move = true;
    if (ncells2 == 1) x2Move = false;
    if (ncells3 == 1) x3Move = false;
  } else if (COORDINATE_SYSTEM == "cylindrical") {
    x1Move = true;
    x2Move = false;
    x3Move = false;
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    x1Move = true;
    x2Move = false;
    x3Move = false;
  }




   //Allocate velocities and grid data
  //v1f.NewAthenaArray((ncells1+1));
  //v2f.NewAthenaArray((ncells2+1));
  //v3f.NewAthenaArray((ncells3+1));
  x1_0.NewAthenaArray((ncells1+1));
  x2_0.NewAthenaArray((ncells2+1));
  x3_0.NewAthenaArray((ncells3+1));
  x1_1.NewAthenaArray((ncells1+1));
  x2_1.NewAthenaArray((ncells2+1));
  x3_1.NewAthenaArray((ncells3+1));
  x1_2.NewAthenaArray((ncells1+1));
  x2_2.NewAthenaArray((ncells2+1));
  x3_2.NewAthenaArray((ncells3+1));
  vf[X1DIR].NewAthenaArray(ncells1+1);
  vf[X2DIR].NewAthenaArray(ncells2+1);
  vf[X3DIR].NewAthenaArray(ncells3+1);
 
  ExpwL.NewAthenaArray((NWAVE+NINT+NSCALARS),ncells3,ncells2,ncells1);
  ExpwR.NewAthenaArray((NWAVE+NINT+NSCALARS),ncells3,ncells2,ncells1);
  if (x1Move){ 
     expFlux[X1DIR].NewAthenaArray(NHYDRO,ncells3,ncells2,ncells1+1);
  }
  if (x2Move) {
    expFlux[X2DIR].NewAthenaArray(NHYDRO,ncells3,ncells2+1,ncells1);
  }
  if (x3Move){
    expFlux[X3DIR].NewAthenaArray(NHYDRO,ncells3+1,ncells2,ncells1);
  }
  xArr.NewAthenaArray(5);
  dxArr.NewAthenaArray(4);
  intPnts.NewAthenaArray(5);
  prims.NewAthenaArray(4);
  coeffL.NewAthenaArray(4);
  coeffR.NewAthenaArray(4);
  

  wi.NewAthenaArray(NHYDRO);
  
  AthenaArray<Real> &v1f = vf[X1DIR];
  AthenaArray<Real> &v2f = vf[X2DIR];
  AthenaArray<Real> &v3f = vf[X3DIR];

  if (x1Move) {
#pragma omp simd
    for (int i=il; i<=iu+1;++i) {
      v1f(i) = 0.0;
      x1_0(i) = pmb->pcoord->x1f(i);
    }
  }
  if (x2Move) {
#pragma omp simd
    for (int j=jl; j<=ju+1;++j) {
      v2f(j) = 0.0;
      x2_0(j) = pmb->pcoord->x2f(j);
    }
  }
  if (x3Move){
#pragma omp simd
    for (int k=kl; k<=ku+1;++k) {
      v3f(k) = 0.0;
      x3_0(k) = pmb->pcoord->x3f(k);
    }
  }
  mydt = (FLT_MAX);
}
// destructor

Expansion::~Expansion() {

  //v1f.DeleteAthenaArray();
  //v2f.DeleteAthenaArray();
  //v3f.DeleteAthenaArray(); 
  x1_0.DeleteAthenaArray();
  x2_0.DeleteAthenaArray();
  x3_0.DeleteAthenaArray();
  x1_1.DeleteAthenaArray();
  x2_1.DeleteAthenaArray();
  x3_1.DeleteAthenaArray();
  x1_2.DeleteAthenaArray();
  x2_2.DeleteAthenaArray();
  x3_2.DeleteAthenaArray();

  ExpwL.DeleteAthenaArray();
  ExpwR.DeleteAthenaArray();
  if (x1Move) {
    expFlux[X1DIR].DeleteAthenaArray();
    vf[X1DIR].DeleteAthenaArray();
  }
  if (x2Move){
    expFlux[X2DIR].DeleteAthenaArray();
    vf[X2DIR].DeleteAthenaArray();
  }
  if (x3Move) {
    expFlux[X3DIR].DeleteAthenaArray();
    vf[X3DIR].DeleteAthenaArray();
  }
  dxArr.DeleteAthenaArray();
  xArr.DeleteAthenaArray();
  intPnts.DeleteAthenaArray();
  prims.DeleteAthenaArray();
  coeffL.DeleteAthenaArray();
  coeffR.DeleteAthenaArray();
  wi.DeleteAthenaArray();

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
 
  AthenaArray<Real> &v1f = vf[X1DIR];
  AthenaArray<Real> &v2f = vf[X2DIR];
  AthenaArray<Real> &v3f = vf[X3DIR];
  if (x1Move) {
#pragma omp simd
    for (int i=il; i<=iu+1; ++i) {
      x1_0(i) += dt*v1f(i);
    }
  }
  if (x2Move) {
#pragma omp simd
    for (int j=jl; j<=ju+1; ++j) {
      x2_0(j) += dt*v2f(j);
    }
  }
  if (x3Move) {
#pragma omp simd
    for (int k=kl; k<=ku+1; ++k) {
      x3_0(k) += dt*v3f(k);
    }
  }
  
  return;
}



void Expansion::ExpansionSourceTerms(const Real dt, const AthenaArray<Real> *flux, const AthenaArray<Real> &prim, AthenaArray<Real> &cons) {
  Real vm = 0.0;
  Coordinates *pc = pmy_block->pcoord;
  Real dx1 = 0.0;
  Real dx2 = 0.0;
  Real dx3 = 0.0;
  Real volume = 0.0;
  Real newVol = 0.0;
  Real oldVol = 0.0;
  Real vol = 0.0;
  Real A1 = 0.0;
  Real A2 = 0.0;
  Real divF1 = 0.0;
  Real divF2 = 0.0; 
  Real divF3 = 0.0;
  Real gm1 = pmy_block->peos->GetGamma()-1.0;
  MeshBlock *pmb = pmy_block;
  AthenaArray<Real> &x1flux=expFlux[X1DIR];
  AthenaArray<Real> &x2flux=expFlux[X2DIR];
  AthenaArray<Real> &x3flux=expFlux[X3DIR];

  AthenaArray<Real> &v1f = vf[X1DIR];
  AthenaArray<Real> &v2f = vf[X2DIR];
  AthenaArray<Real> &v3f = vf[X3DIR];
  
  for (int k = ks; k<=ke;++k) {
    for (int j = js; j<=je;++j) {      
#pragma omp simd
      for (int i = is; i<=ie;++i) {
	if (COORDINATE_SYSTEM == "cartesian") {
	  dx1 = pc->dx1f(i)+v1f(i+1)*dt - v1f(i)*dt ;          
	  dx2 = pc->dx2f(j)+v2f(j+1)*dt - v2f(j)*dt ;          
	  dx3 = pc->dx3f(k)+v3f(k+1)*dt - v3f(k)*dt ;          
	  volume = dx1*dx2*dx3;          
	} else if (COORDINATE_SYSTEM == "cylindrical") {
	  dx1 = pc->coord_vol_i_(i) + dt*(v1f(i+1)*pc->x1f(i+1)-v1f(i)*pc->x1f(i)) 
				    + 0.5*pow(dt,2.0)*(pow(v1f(i+1),2.0) - pow(v1f(i),2.0));          
	  dx2 = pc->dx2f(j)+v2f(j+1)*dt - v2f(j)*dt ;          
	  dx3 = pc->dx3f(k)+v3f(k+1)*dt - v3f(k)*dt ;          
	  volume = dx1*dx2*dx3;         
	} else if (COORDINATE_SYSTEM == "spherical_polar") {
	  dx1 = pc->coord_vol_i_(i) + dt*(v1f(i+1)*pow(pc->x1f(i+1),2.0)-v1f(i)*pow(pc->x1f(i),2.0)) 
				    + pow(dt,2.0)*(pow(v1f(i+1),2.0)*pc->x1f(i+1) - pow(v1f(i),2.0)*pc->x1f(i))
				    + 1.0/3.0*pow(dt,3.0)*(pow(v1f(i+1),3.0) - pow(v1f(i),3.0));          
	  dx2 = fabs(cos(pc->x2f(j)+v2f(j)*dt) - cos(pc->x2f(j+1)+v2f(j+1)*dt));          
	  dx3 = pc->dx3f(k)+v3f(k+1)*dt - v3f(k)*dt ;          
	  volume = dx1*dx2*dx3;         
	} 
	newVol = volume;
	oldVol = pc->GetCellVolume(k,j,i);
 

        for (int n = 0; n<(NHYDRO+NSCALARS); ++n) { 
          vol = pmb->pcoord->GetCellVolume(k,j,i);
          divF1 = 0.0;
          divF2 = 0.0;
          divF3 = 0.0;
          A1 = pmb->pcoord->GetFace1Area(k,j,i);
          A2 = pmb->pcoord->GetFace1Area(k,j,i+1);
          if (x1Move) {
            divF1 = x1flux(n,k,j,i+1)*A2 - x1flux(n,k,j,i) * A1;
          }
               
          if (x2Move){
            A1 = pmb->pcoord->GetFace2Area(k,j,i);
            A2 = pmb->pcoord->GetFace2Area(k,j+1,i);
            divF2 = x2flux(n,k,j+1,i)*A2 - x2flux(n,k,j,i) * A1;
          }
          if (x3Move){
            A1 = pmb->pcoord->GetFace3Area(k,j,i);
            A2 = pmb->pcoord->GetFace3Area(k+1,j,i);
            divF3 = x3flux(n,k+1,j,i)*A2 - x3flux(n,k,j,i) * A1;
          }
          cons(n,k,j,i) += dt/vol*(divF1 + divF2 + divF3);        
  
          cons(n,k,j,i) *= oldVol/newVol;        
        }
      }
    }
  }

  return;
}


void Expansion::UpdateVelData(MeshBlock *pmb ,Real time, Real dt){
  AthenaArray<Real> &v1f = vf[X1DIR];
  AthenaArray<Real> &v2f = vf[X2DIR];
  AthenaArray<Real> &v3f = vf[X3DIR];
  mydt = dt;
  if (pmb->pmy_mesh->GridDiffEq_ != NULL) {
    //Edit each delx1f, delx2f, delx3f before source terms and editing of grid      
    if (x3Move) {
      for (int k = kl; k<=ku+1;++k){
        v3f(k) = pmb->pmy_mesh->GridDiffEq_(pmb->pcoord->x3f(k),k,pmb->pmy_mesh->time,dt,3,pmb->pmy_mesh->GridData);
      }
    }
    if (x2Move) {
      for (int j = jl;j<=ju+1;++j){      
        v2f(j) = pmb->pmy_mesh->GridDiffEq_(pmb->pcoord->x2f(j),j,pmb->pmy_mesh->time,dt,2,pmb->pmy_mesh->GridData);
      }
    }
    if (x1Move) {
      for (int i = il; i<=iu+1;++i){
        v1f(i) = pmb->pmy_mesh->GridDiffEq_(pmb->pcoord->x1f(i),i,pmb->pmy_mesh->time,dt,1,pmb->pmy_mesh->GridData); 
      }
    }
  } 
  return;
}

void Expansion::GridEdit(MeshBlock *pmb,bool lastStage){  
  //FACE CENTERED
  //x1
  if (x1Move) {
    for (int i=il; i<=iu+1; ++i){
      pmb->pcoord->x1f(i) = x1_0(i);  
    }
    for (int i=il; i<=iu; ++i) {
      pmb->pcoord->dx1f(i) = pmb->pcoord->x1f(i+1)-pmb->pcoord->x1f(i);
    }
  }
  
  //x2  
  if (x2Move){
    for (int j=jl; j<=ju+1; ++j){
      pmb->pcoord->x2f(j) = x2_0(j);
    }
    for (int j=jl; j<=ju; ++j) {
      pmb->pcoord->dx2f(j) = pmb->pcoord->x2f(j+1)- pmb->pcoord->x2f(j);
    }
  }                      
  //x3
  if (x3Move){
    for (int k=kl; k<=ku+1; ++k){
      pmb->pcoord->x3f(k) = x3_0(k);
    }
    for (int k=kl; k<=ku; ++k) {
      pmb->pcoord->dx3f(k) = pmb->pcoord->x3f(k+1)- pmb->pcoord->x3f(k);
    }
  }

  //Set Reconstruction Coefficients x1
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

  if (pmb->block_size.nx2 !=1){
    for (int j=jl+1; j<=ju-1; ++j) {
      Real& dx_jm1 = pmb->pcoord->dx2f(j-1);
      Real& dx_j   = pmb->pcoord->dx2f(j  );
      Real& dx_jp1 = pmb->pcoord->dx2f(j+1);
      Real qe = dx_j/(dx_jm1 + dx_j + dx_jp1);       // Outermost coeff in CW eq 1.7
      pmb->precon->c1j(j) = qe*(2.0*dx_jm1+dx_j)/(dx_jp1 + dx_j); // First term in CW eq 1.7
      pmb->precon->c2j(j) = qe*(2.0*dx_jp1+dx_j)/(dx_jm1 + dx_j); // Second term in CW eq 1.7
      if (j > jl+1) {  // c3-c6 are not computed in first iteration
        Real& dx_jm2 = pmb->pcoord->dx2f(j-2);
        Real qa = dx_jm2 + dx_jm1 + dx_j + dx_jp1;
        Real qb = dx_jm1/(dx_jm1 + dx_j);
        Real qc = (dx_jm2 + dx_jm1)/(2.0*dx_jm1 + dx_j);
        Real qd = (dx_jp1 + dx_j)/(2.0*dx_j + dx_jm1);
        qb = qb + 2.0*dx_j*qb/qa*(qc-qd);
        pmb->precon->c3j(j) = 1.0 - qb;
        pmb->precon->c4j(j) = qb;
        pmb->precon->c5j(j) = dx_j/qa*qd;
        pmb->precon->c6j(j) = -dx_jm1/qa*qc;
      }
    }
    // Compute curvilinear geometric factors for limiter (Mignone eq 48)
    for (int j=jl+1; j<=je-1; ++j) {
      // corrections to PPMx2 only for spherical polar coordinates
      if (COORDINATE_SYSTEM == "spherical_polar") {
        // x2 = theta polar coordinate adjustment
        Real h_plus, h_minus;
        Real& dx_j   = pmb->pcoord->dx2f(j);
        Real& xf_j   = pmb->pcoord->x2f(j);
        Real& xf_jp1   = pmb->pcoord->x2f(j+1);
        Real dmu = cos(xf_j) - cos(xf_jp1);
        Real dmu_tilde = sin(xf_j) - sin(xf_jp1);
        h_plus = (dx_j*(dmu_tilde + dx_j*cos(xf_jp1)))/(
                  dx_j*(sin(xf_j) + sin(xf_jp1)) - 2.0*dmu);
        h_minus = -(dx_j*(dmu_tilde + dx_j*cos(xf_j)))/(
                    dx_j*(sin(xf_j) + sin(xf_jp1)) - 2.0*dmu);
        pmb->precon->hplus_ratio_j(j) = (h_plus + 1.0)/(h_minus - 1.0);
        pmb->precon->hminus_ratio_j(j) = (h_minus + 1.0)/(h_plus - 1.0);
      } else {
        // h_plus = 3.0;
        // h_minus = 3.0;
        // Ratios are both = 2, as in orig PPM overshoot limiter
        pmb->precon->hplus_ratio_j(j) = 2.0;
        pmb->precon->hminus_ratio_j(j) = 2.0;
      }
    }
  }


  if (pmb->block_size.nx3 !=1){
    for (int k=kl+1; k<=ku-1; ++k) {
      Real& dx_km1 = pmb->pcoord->dx3f(k-1);
      Real& dx_k   = pmb->pcoord->dx3f(k  );
      Real& dx_kp1 = pmb->pcoord->dx3f(k+1);
      Real qe = dx_k/(dx_km1 + dx_k + dx_kp1);       // Outermost coeff in CW eq 1.7
      pmb->precon->c1k(k) = qe*(2.0*dx_km1+dx_k)/(dx_kp1 + dx_k); // First term in CW eq 1.7
      pmb->precon->c2k(k) = qe*(2.0*dx_kp1+dx_k)/(dx_km1 + dx_k); // Second term in CW eq 1.7
  
      if (k > kl+1) {  // c3-c6 are not computed in first iteration
        Real& dx_km2 = pmb->pcoord->dx3f(k-2);
        Real qa = dx_km2 + dx_km1 + dx_k + dx_kp1;
        Real qb = dx_km1/(dx_km1 + dx_k);
        Real qc = (dx_km2 + dx_km1)/(2.0*dx_km1 + dx_k);
        Real qd = (dx_kp1 + dx_k)/(2.0*dx_k + dx_km1);
        qb = qb + 2.0*dx_k*qb/qa*(qc-qd);
        pmb->precon->c3k(k) = 1.0 - qb;
        pmb->precon->c4k(k) = qb;
        pmb->precon->c5k(k) = dx_k/qa*qd;
        pmb->precon->c6k(k) = -dx_km1/qa*qc;
      }
    }
    // Compute curvilinear geometric factors for limiter (Mignone eq 48)
    // No corrections in x3 for the built-in Newtonian coordinate systems
    for (int k=kl+1; k<=ku-1; ++k) {
      // h_plus = 3.0;
      // h_minus = 3.0;
      // Ratios are both = 2 for Cartesian and all curviliniear coords
      pmb->precon->hplus_ratio_k(k) = 2.0;
      pmb->precon->hminus_ratio_k(k) = 2.0;
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
    } else if (x2Move) {
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
    } else if (x3Move) {
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
  if (x1Move) {
    //x1
    Real inner = 0.0;
    Real outer = 0.0;
    inner = x1_0(is);
    outer = x1_0(ie+1);
    pmb->block_size.x1min = inner;
    pmb->block_size.x1max = outer;
  }
  if (x2Move) {
    //x2 
    Real inner = 0.0;
    Real outer = 0.0;
    inner = x2_0(js);
    outer = x2_0(je+1);
    pmb->block_size.x2min = inner;
    pmb->block_size.x2max = outer;
  }
  if (x3Move) {
    //x3
   
    Real inner = 0.0;
    Real outer = 0.0;
    inner = x3_0(ks);
    outer = x3_0(ke+1);
    pmb->block_size.x3min = inner;
    pmb->block_size.x3max = outer;
  }
   
  return;
}


Real Expansion::GridTimeStep(MeshBlock *pmb){ 
  
  Real nextPosDelta, minCellSize;//, dt;
  
  Mesh *pmesh = pmb->pmy_mesh;

  Real mydt = pmesh->dt;

  Real min_dt = (FLT_MAX);
  Real dtEx = min_dt;
  int count = 0;
  Real overStep = 0.0;
  if (x1Move) { 
    for (int i = il; i<=iu+1;i++){
      nextPosDelta = pmesh->GridDiffEq_(pmb->pcoord->x1f(i),i,pmesh->time,mydt,1,pmesh->GridData)*mydt;
  
      if (nextPosDelta < 0 && i != il){ 
        minCellSize = pmb->pcoord->dx1f(i-1);
      } else {
        minCellSize = pmb->pcoord->dx1f(i);
      }
  
      nextPosDelta = fabs(nextPosDelta);
  
      if (nextPosDelta != 0.0 && minCellSize != 0.0){
        overStep = minCellSize - nextPosDelta;
        count = 0;
        while (overStep <= 0){
          mydt *= 0.9;
          nextPosDelta = pmesh->GridDiffEq_(pmb->pcoord->x1f(i),i,pmesh->time,mydt,1,pmesh->GridData)*mydt;
          nextPosDelta = fabs(nextPosDelta);
          overStep = minCellSize - nextPosDelta;      
          count++;
          if (count >=20) {
            mydt *= 0.5;
            count = 0.0;
          } 
        }
        dtEx = mydt;//*pmesh->cfl_number;//fabs(minCellSize* 0.5/nextPosDelta * (pmesh->dt));
        min_dt = std::min(min_dt, dtEx);
      } else {
        dtEx = pmesh->dt;
        min_dt = std::min(min_dt,dtEx);
      }
    }
  }
  
  if (x2Move) {
    for (int j = jl; j<=ju+1; j++){
      nextPosDelta = pmesh->GridDiffEq_(pmb->pcoord->x2f(j),j,pmesh->time,mydt,2,pmesh->GridData)*mydt;
  
      if (nextPosDelta < 0 && j != jl){ 
        minCellSize = pmb->pcoord->dx2f(j-1);
      } else {
        minCellSize = pmb->pcoord->dx2f(j);
      }
  
      nextPosDelta = fabs(nextPosDelta);
  
      if (nextPosDelta != 0.0 && minCellSize != 0.0){
        overStep = minCellSize - nextPosDelta;
        count = 0;
        while (overStep <= 0){
          mydt *= 0.9;
          nextPosDelta = pmesh->GridDiffEq_(pmb->pcoord->x2f(j),j,pmesh->time,mydt,2,pmesh->GridData)*mydt;
          nextPosDelta = fabs(nextPosDelta);
          overStep = minCellSize - nextPosDelta;      
          count++;
          if (count >=20) {
            mydt *= 0.5;
            count = 0.0;
          } 
        }
        dtEx = mydt;//*pmesh->cfl_number;//fabs(minCellSize* 0.5/nextPosDelta * (pmesh->dt));
        min_dt = std::min(min_dt, dtEx);
      } else {
        dtEx = pmesh->dt;
        min_dt = std::min(min_dt,dtEx);
      }
    }
  }

  if (x3Move) {
    for (int k = kl; k<=ku+1; k++){
      nextPosDelta = pmesh->GridDiffEq_(pmb->pcoord->x3f(k),k,pmesh->time,mydt,3,pmesh->GridData)*mydt;
  
      if (nextPosDelta < 0 && k != kl){ 
        minCellSize = pmb->pcoord->dx3f(k-1);
      } else {
        minCellSize = pmb->pcoord->dx3f(k);
      }
  
      nextPosDelta = fabs(nextPosDelta);
  
      if (nextPosDelta != 0.0 && minCellSize != 0.0){
        overStep = minCellSize - nextPosDelta;
        count = 0;
        while (overStep <= 0){
          mydt *= 0.9;
          nextPosDelta = pmesh->GridDiffEq_(pmb->pcoord->x3f(k),k,pmesh->time,mydt,3,pmesh->GridData)*mydt;
          nextPosDelta = fabs(nextPosDelta);
          overStep = minCellSize - nextPosDelta;      
          count++;
          if (count >=20) {
            mydt *= 0.5;
            count = 0.0;
          } 
        }
        dtEx = mydt;//*pmesh->cfl_number;//fabs(minCellSize* 0.5/nextPosDelta * (pmesh->dt));
        min_dt = std::min(min_dt, dtEx);
      } else {
        dtEx = pmesh->dt;
        min_dt = std::min(min_dt,dtEx);
      }
    }
  }

  return min_dt;
}


