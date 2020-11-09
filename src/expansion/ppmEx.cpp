//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ppmEX.cpp
//  \brief piecewise parabolic reconstruction with modified McCorquodale/Colella limiter
//         for a uniform Cartesian mesh, Mignone limiter for nonuniform mesh, Accounting
//         for expanding grid
//
// REFERENCES:
// (CW) P. Colella & P. Woodward, "The Piecewise Parabolic Method (PPM) for Gas-Dynamical
// Simulations", JCP, 54, 174 (1984)
//
// (CS) P. Colella & M. Sekora, "A limiter for PPM that preserves accuracy at smooth
// extrema", JCP, 227, 7069 (2008)
//
// (MC) P. McCorquodale & P. Colella,  "A high-order finite-volume method for conservation
// laws on locally refined grids", CAMCoS, 6, 1 (2011)
//
// (CD) P. Colella, M.R. Dorr, J. Hittinger, D. Martin, "High-order, finite-volume methods
// in mapped coordinates", JCP, 230, 2952 (2011)
//
// (Mignone) A. Mignone, "High-order conservative reconstruction schemes for finite volume
// methods in cylindrical and spherical coordinates", JCP, 270, 784 (2014)
//========================================================================================

// C++ headers
#include <algorithm>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "expansion.hpp"


Real WallIntegration(double xArr[],double prims[],Real myX1, Real myX2, int xArrSize);
Real CenterIntegration(double xArr[],double prims[], Real myX, int xArrSize);

Real WallIntegration(double xArr[], double prims[], Real myX1, Real myX2, int xArrSize) {
  // Prims are cell centered, xArr has cell wall locations 
  // Prims should be length one less than xArr
  int n  = xArrSize;
  double dxArr[n-1];
  double intPnts[n];
  int m,p;
  double val;
  Real a1, a2, a3, a4;
  Real b11, b12, b13, b14;
  Real b21, b22, b23, b24;
  Real b31, b32, b33, b34;
  Real g4, g3, g2, g1;
 
  for (m=0;m<=(n-2);++m) { 
     dxArr[m]=xArr[m+1]-xArr[m];
  }

  for (m=0;m<=(n-1);++m) {
    intPnts[m] = 0.0;
  }


  for (m=1;m<=(n-1);++m) {
    for (p=m;p<=(n-1);++p){
      intPnts[p] += dxArr[m-1]*prims[m-1];
    }
  }
  a1 = -1.0*intPnts[1]
	/((dxArr[0])*(dxArr[1])*(dxArr[1]+dxArr[2])*(dxArr[1]+dxArr[2]+dxArr[3]));
  a2 = 1.0*intPnts[2]
	/((dxArr[0]+dxArr[1])*(dxArr[1])*(dxArr[2])*(dxArr[2]+dxArr[3]));
  a3 = -1.0*intPnts[3]
        /((dxArr[0]+dxArr[1]+dxArr[2])*(dxArr[1]+dxArr[2])*(dxArr[2])*(dxArr[3]));
  a4 = 1.0*intPnts[4]
	/((dxArr[0]+dxArr[1]+dxArr[2]+dxArr[3])*(dxArr[1]+dxArr[2]+dxArr[3])*(dxArr[2]+dxArr[3])*(dxArr[3]));


  b11 = xArr[0]*xArr[2]*xArr[3] + xArr[0]*xArr[2]*xArr[4] + xArr[0]*xArr[3]*xArr[4] + xArr[2]*xArr[3]*xArr[4];
  b12 = xArr[0]*xArr[1]*xArr[3] + xArr[0]*xArr[1]*xArr[4] + xArr[0]*xArr[3]*xArr[4] + xArr[1]*xArr[3]*xArr[4];
  b13 = xArr[0]*xArr[1]*xArr[2] + xArr[0]*xArr[1]*xArr[4] + xArr[0]*xArr[2]*xArr[4] + xArr[1]*xArr[2]*xArr[4];
  b14 = xArr[0]*xArr[1]*xArr[2] + xArr[0]*xArr[1]*xArr[3] + xArr[0]*xArr[2]*xArr[3] + xArr[1]*xArr[2]*xArr[3];

  b21 = xArr[0]*xArr[2] + xArr[0]*xArr[3] + xArr[2]*xArr[3] + xArr[0]*xArr[4] + xArr[2]*xArr[4] + xArr[3]*xArr[4];
  b22 = xArr[0]*xArr[1] + xArr[0]*xArr[3] + xArr[1]*xArr[3] + xArr[0]*xArr[4] + xArr[1]*xArr[4] + xArr[3]*xArr[4];
  b23 = xArr[0]*xArr[1] + xArr[0]*xArr[2] + xArr[1]*xArr[2] + xArr[0]*xArr[4] + xArr[1]*xArr[4] + xArr[2]*xArr[4];
  b24 = xArr[0]*xArr[1] + xArr[0]*xArr[2] + xArr[1]*xArr[2] + xArr[0]*xArr[3] + xArr[1]*xArr[3] + xArr[2]*xArr[3];

  b31 = xArr[0] + xArr[2] + xArr[3] + xArr[4];
  b32 = xArr[0] + xArr[1] + xArr[3] + xArr[4];
  b33 = xArr[0] + xArr[1] + xArr[2] + xArr[4];
  b34 = xArr[0] + xArr[1] + xArr[2] + xArr[3];

  g4 = a1 + a2 + a3 + a4;  
  g3 = -1.0*(a1*b31 + a2*b32 + a3*b33 + a4*b34);
  g2 =  1.0*(a1*b21 + a2*b22 + a3*b23 + a4*b24);
  g1 = -1.0*(a1*b11 + a2*b12 + a3*b13 + a4*b14);

  if (myX2!=myX1){
    val = (g4*(pow(myX2,4.0)-pow(myX1,4.0)) + g3*(pow(myX2,3.0)-pow(myX1,3.0)) 
	 + g2*(pow(myX2,2.0)-pow(myX1,2.0)) + g1*(pow(myX2,1.0)-pow(myX1,1.0)))/(myX2-myX1);
  } else {
    val = 4.0*g4*pow(myX1,3.0) + 3.0*g3*pow(myX1,2.0) + 2.0*g2*pow(myX1,1.0) + g1;
  }

  return val;
}

Real CenterIntegration(double xArr[],double prims[], Real myX, int xArrSize) {
  int n = xArrSize;
  double termsC[5];
  double Myx;
  int m,l;
  Real valCi = 0.0;
  for (m=0;m<=4;++m){
    termsC[m] = 1.0;
    for (l=0;l<=4;++l){
      if (l != m) {
        termsC[m] *= (myX-xArr[l])/(xArr[m]-xArr[l]);
      }
    }
    valCi += prims[m]*termsC[m];
  }       
  return valCi; 

}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseParabolicX1()
//  \brief Returns L/R interface values in X1-dir constructed using fourth-order PPM and
//         Colella-Sekora or Mignone limiting over [kl,ku][jl,ju][il,iu]
void Expansion::PiecewiseLinearOffsetX1(MeshBlock *pmb,
  const int kl, const int ku, const int jl, const int ju, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wL, AthenaArray<Real> &wR, const Real dt) {


  return;
}

void Expansion::PiecewiseParabolicOffsetX1(MeshBlock *pmb,
  const int kl, const int ku, const int jl, const int ju, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wL, AthenaArray<Real> &wR, const Real dt) {

  double xArrL[6], xArrR[6];
  double primsL[5], primsR[5];
  double MyxL, MyxR;
  double valLi, valRi;

  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
    // cache the x1-sliced primitive states for eigensystem calculation
    for (int n=0; n<(NHYDRO); ++n) {

      //Deal with Left Boundary
      xArrL[0] = pmb->pcoord->x1f(il-3);
      xArrL[1] = pmb->pcoord->x1f(il-2);
      xArrL[2] = pmb->pcoord->x1f(il-1);
      xArrL[3] = pmb->pcoord->x1f(il  );
      xArrL[4] = pmb->pcoord->x1f(il+1);
      
      primsL[0] = w(n,k,j,il-3);
      primsL[1] = w(n,k,j,il-2);
      primsL[2] = w(n,k,j,il-1);
      primsL[3] = w(n,k,j,il);
  
      MyxL = xArrL[2]+v1f(il-1)*dt;	
      valLi =  WallIntegration(xArrL,primsL,xArrL[2],MyxL,5);
 
      if (valLi > std::max(primsL[2],primsL[1])) valLi = std::max(primsL[2],primsL[1]);
      if (valLi < std::min(primsL[2],primsL[1])) valLi = std::min(primsL[2],primsL[1]);
    
      xArrR[0] = pmb->pcoord->x1f(il-2);
      xArrR[1] = pmb->pcoord->x1f(il-1);
      xArrR[2] = pmb->pcoord->x1f(il  );
      xArrR[3] = pmb->pcoord->x1f(il+1);
      xArrR[4] = pmb->pcoord->x1f(il+2);
          
      primsR[0] = w(n,k,j,il-2);
      primsR[1] = w(n,k,j,il-1);
      primsR[2] = w(n,k,j,il  );
      primsR[3] = w(n,k,j,il+1);
 
      MyxR = xArrR[2]+v1f(il)*dt;	
      valRi =  WallIntegration(xArrR,primsR,xArrR[2],MyxR,5);
  
      if (valRi > std::max(primsR[2],primsR[1])) valRi = std::max(primsR[2],primsR[1]);
      if (valRi < std::min(primsR[2],primsR[1])) valRi = std::min(primsR[2],primsR[1]);
  
      //Extrema Detector
      Real dqm = w(n,k,j,il-1) - valLi;
      Real dqp = valRi-w(n,k,j,il-1);
    
      if (dqp*dqm <= 0.0){
        valLi = w(n,k,j,il-1);
        valRi = w(n,k,j,il-1);
      } else {
        // Overshoot i-1/2,R / i,(-) state
        if (fabs(dqm) >= pmb->precon->hplus_ratio_i(il-1)*fabs(dqp)) {
          valLi = w(n,k,j,il-1) - pmb->precon->hplus_ratio_i(il-1)*dqp;
        }
        // Overshoot i+1/2,L / i,(+) state
        if (fabs(dqp) >= pmb->precon->hminus_ratio_i(il-1)*fabs(dqm)) {
          valRi = w(n,k,j,il-1) + pmb->precon->hminus_ratio_i(il-1)*dqm;
        }
      }
      wL(n,k,j,il) = valRi;
      wR(n,k,j,il-1) = valLi;

      //Deal with Right Boundary
      xArrL[0] = pmb->pcoord->x1f(iu-2);
      xArrL[1] = pmb->pcoord->x1f(iu-1);
      xArrL[2] = pmb->pcoord->x1f(iu  );
      xArrL[3] = pmb->pcoord->x1f(iu+1);
      xArrL[4] = pmb->pcoord->x1f(iu+2);
      
      primsL[0] = w(n,k,j,iu-2);
      primsL[1] = w(n,k,j,iu-1);
      primsL[2] = w(n,k,j,iu  );
      primsL[3] = w(n,k,j,iu+1);
  
      MyxL = xArrL[2]+v1f(iu)*dt;	
      valLi =  WallIntegration(xArrL,primsL,xArrL[2],MyxL,5);
 
      if (valLi > std::max(primsL[2],primsL[1])) valLi = std::max(primsL[2],primsL[1]);
      if (valLi < std::min(primsL[2],primsL[1])) valLi = std::min(primsL[2],primsL[1]);
    
      xArrR[0] = pmb->pcoord->x1f(iu-1);
      xArrR[1] = pmb->pcoord->x1f(iu);
      xArrR[2] = pmb->pcoord->x1f(iu+1);
      xArrR[3] = pmb->pcoord->x1f(iu+2);
      xArrR[4] = pmb->pcoord->x1f(iu+3);
          
      primsR[0] = w(n,k,j,iu-1);
      primsR[1] = w(n,k,j,iu  );
      primsR[2] = w(n,k,j,iu+1);
      primsR[3] = w(n,k,j,iu+2);
 
      MyxR = xArrR[2]+v1f(iu+1)*dt;	
      valRi =  WallIntegration(xArrR,primsR,xArrR[2],MyxR,5);
  
      if (valRi > std::max(primsR[2],primsR[1])) valRi = std::max(primsR[2],primsR[1]);
      if (valRi < std::min(primsR[2],primsR[1])) valRi = std::min(primsR[2],primsR[1]);
  
      //Extrema Detector
      dqm = w(n,k,j,iu) - valLi;
      dqp = valRi-w(n,k,j,iu);
    
      if (dqp*dqm <= 0.0){
        valLi = w(n,k,j,iu);
        valRi = w(n,k,j,iu);
      } else {
        // Overshoot i-1/2,R / i,(-) state
        if (fabs(dqm) >= pmb->precon->hplus_ratio_i(iu)*fabs(dqp)) {
          valLi = w(n,k,j,iu) - pmb->precon->hplus_ratio_i(iu)*dqp;
        }
        // Overshoot i+1/2,L / i,(+) state
        if (fabs(dqp) >= pmb->precon->hminus_ratio_i(iu)*fabs(dqm)) {
          valRi = w(n,k,j,iu) + pmb->precon->hminus_ratio_i(iu)*dqm;
        }
      }
      wL(n,k,j,iu+1) = valRi;
      wR(n,k,j,iu) = valLi;
#pragma omp simd

      for (int i=il; i<=iu-1; ++i) {
        if (v1f(i) > 0.0) {
          //Get i-1/2 average state
          xArrL[0] = pmb->pcoord->x1f(i-1);
          xArrL[1] = pmb->pcoord->x1f(i  );
          xArrL[2] = pmb->pcoord->x1f(i+1);
          xArrL[3] = pmb->pcoord->x1f(i+2);
          xArrL[4] = pmb->pcoord->x1f(i+3);
          
          primsL[0] = w(n,k,j,i-1);
          primsL[1] = w(n,k,j,i  );
          primsL[2] = w(n,k,j,i+1  );
          primsL[3] = w(n,k,j,i+2);
    
          MyxL = xArrL[1]+v1f(i)*dt;	
          valLi =  WallIntegration(xArrL,primsL,xArrL[1],MyxL,5);
      
          //Get i+1/2 average state
          xArrR[0] = pmb->pcoord->x1f(i);
          xArrR[1] = pmb->pcoord->x1f(i+1);
          xArrR[2] = pmb->pcoord->x1f(i+2);
          xArrR[3] = pmb->pcoord->x1f(i+3);
          xArrR[4] = pmb->pcoord->x1f(i+4);
            
          primsR[0] = w(n,k,j,i);
          primsR[1] = w(n,k,j,i+1);
          primsR[2] = w(n,k,j,i+2);
          primsR[3] = w(n,k,j,i+3);
   
          MyxR = xArrR[1]+v1f(i+1)*dt;	
          valRi =  WallIntegration(xArrR,primsR,xArrR[1],MyxR,5);

        } else if (v1f(i) < 0.0) {
          //Get i-1/2 average state
          xArrL[0] = pmb->pcoord->x1f(i-3);
          xArrL[1] = pmb->pcoord->x1f(i-2);
          xArrL[2] = pmb->pcoord->x1f(i-1);
          xArrL[3] = pmb->pcoord->x1f(i  );
          xArrL[4] = pmb->pcoord->x1f(i+1);
          
          primsL[0] = w(n,k,j,i-3);
          primsL[1] = w(n,k,j,i-2);
          primsL[2] = w(n,k,j,i-1);
          primsL[3] = w(n,k,j,i  );
    
          MyxL = xArrL[3]+v1f(i)*dt;	
          valLi =  WallIntegration(xArrL,primsL,xArrL[3],MyxL,5);
      
          //Get i+1/2 average state
          xArrR[0] = pmb->pcoord->x1f(i-2);
          xArrR[1] = pmb->pcoord->x1f(i-1);
          xArrR[2] = pmb->pcoord->x1f(i  );
          xArrR[3] = pmb->pcoord->x1f(i+1);
          xArrR[4] = pmb->pcoord->x1f(i+2);
            
          primsR[0] = w(n,k,j,i-2);
          primsR[1] = w(n,k,j,i-1);
          primsR[2] = w(n,k,j,i );
          primsR[3] = w(n,k,j,i+1);
   
          MyxR = xArrR[2]+v1f(i+1)*dt;	
          valRi =  WallIntegration(xArrR,primsR,xArrR[2],MyxR,5);

        } else if (v1f(i) == 0.0) {
          //Get i-1/2 average state
          xArrL[0] = pmb->pcoord->x1f(i-2);
          xArrL[1] = pmb->pcoord->x1f(i-1);
          xArrL[2] = pmb->pcoord->x1f(i);
          xArrL[3] = pmb->pcoord->x1f(i+1);
          xArrL[4] = pmb->pcoord->x1f(i+2);
          
          primsL[0] = w(n,k,j,i-2);
          primsL[1] = w(n,k,j,i-1);
          primsL[2] = w(n,k,j,i  );
          primsL[3] = w(n,k,j,i+1);
    
          MyxL = xArrL[2];//+v1f(i)*dt;	
          valLi =  WallIntegration(xArrL,primsL,xArrL[2],MyxL,5);
      
          //Get i+1/2 average state
          xArrR[0] = pmb->pcoord->x1f(i-1);
          xArrR[1] = pmb->pcoord->x1f(i);
          xArrR[2] = pmb->pcoord->x1f(i+1);
          xArrR[3] = pmb->pcoord->x1f(i+2);
          xArrR[4] = pmb->pcoord->x1f(i+3);
            
          primsR[0] = w(n,k,j,i-1);
          primsR[1] = w(n,k,j,i  );
          primsR[2] = w(n,k,j,i+1);
          primsR[3] = w(n,k,j,i+2);
   
          MyxR = xArrR[2];//+v1f(i+1)*dt;	
          valRi =  WallIntegration(xArrR,primsR,xArrR[2],MyxR,5);     
        }

        //SLOPE LIMITING i-1/2
        if (valLi > std::max(w(n,k,j,i),w(n,k,j,i-1))) valLi = std::max(w(n,k,j,i),w(n,k,j,i-1));
        if (valLi < std::min(w(n,k,j,i),w(n,k,j,i-1))) valLi = std::min(w(n,k,j,i),w(n,k,j,i-1));
        //SLOPE LIMITING i+1/2
        if (valRi > std::max(w(n,k,j,i+1),w(n,k,j,i))) valRi = std::max(w(n,k,j,i+1),w(n,k,j,i));
        if (valRi < std::min(w(n,k,j,i+1),w(n,k,j,i))) valRi = std::min(w(n,k,j,i+1),w(n,k,j,i));

        //Extrema Detector
        dqm = w(n,k,j,i) - valLi;
        dqp = valRi-w(n,k,j,i);
      
        if (dqp*dqm <= 0.0){
          valLi = w(n,k,j,i);
          valRi = w(n,k,j,i);
        } else {
          // Overshoot i-1/2,R / i,(-) state
          if (fabs(dqm) >= pmb->precon->hplus_ratio_i(i)*fabs(dqp)) {
            valLi = w(n,k,j,i) - pmb->precon->hplus_ratio_i(i)*dqp;
          }
          // Overshoot i+1/2,L / i,(+) state
          if (fabs(dqp) >= pmb->precon->hminus_ratio_i(i)*fabs(dqm)) {
            valRi = w(n,k,j,i) + pmb->precon->hminus_ratio_i(i)*dqm;
          }
        }
   
        wL(n,k,j,i+1) = valRi;
        wR(n,k,j,i  ) = valLi;
      }
    }
  }}

  return;
}







