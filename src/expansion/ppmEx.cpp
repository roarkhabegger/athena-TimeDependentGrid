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

void WallIntegration(double xArr[], double prims[],  int xArrSize, double coeff[]) {
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

  coeff[0] = g1;
  coeff[1] = g2;
  coeff[2] = g3;
  coeff[3] = g4;


  return;
}

Real CenterIntegration(double xArr[],double prims[], Real myX, int xArrSize) {
  int n = xArrSize;
  double dxArr[n-1];
  int m;

  Real a0, a1, a2, a3, a4;
  Real b00, b01, b02, b03, b04;
  Real b10, b11, b12, b13, b14;
  Real b20, b21, b22, b23, b24;
  Real b30, b31, b32, b33, b34;
  Real b40, b41, b42, b43, b44;
  Real g4, g3, g2, g1, g0;
  Real val;

  for (m=0;m<=(n-2);++m) { 
     dxArr[m]=xArr[m+1]-xArr[m];
  }

  a0 = 1.0*prims[0]
	/((dxArr[0])*(dxArr[0]+dxArr[1])*(dxArr[0]+dxArr[1]+dxArr[2])*(dxArr[0]+dxArr[1]+dxArr[2]+dxArr[3]));
  a1 = -1.0*prims[1]
	/((dxArr[0])*(dxArr[1])*(dxArr[1]+dxArr[2])*(dxArr[1]+dxArr[2]+dxArr[3]));
  a2 = 1.0*prims[2]
	/((dxArr[0]+dxArr[1])*(dxArr[1])*(dxArr[2])*(dxArr[2]+dxArr[3]));
  a3 = -1.0*prims[3]
        /((dxArr[0]+dxArr[1]+dxArr[2])*(dxArr[1]+dxArr[2])*(dxArr[2])*(dxArr[3]));
  a4 = 1.0*prims[4]
	/((dxArr[0]+dxArr[1]+dxArr[2]+dxArr[3])*(dxArr[1]+dxArr[2]+dxArr[3])*(dxArr[2]+dxArr[3])*(dxArr[3]));

  b00 = xArr[1]*xArr[2]*xArr[3]*xArr[4];
  b10 = xArr[0]*xArr[2]*xArr[3]*xArr[4];
  b20 = xArr[0]*xArr[1]*xArr[3]*xArr[4];
  b30 = xArr[0]*xArr[1]*xArr[2]*xArr[4];
  b40 = xArr[0]*xArr[1]*xArr[2]*xArr[3];

  b01 = -1.0*(xArr[1]*xArr[2]*xArr[3] + xArr[1]*xArr[2]*xArr[4] + xArr[1]*xArr[3]*xArr[4] + xArr[2]*xArr[3]*xArr[4]);
  b11 = -1.0*(xArr[0]*xArr[2]*xArr[3] + xArr[0]*xArr[2]*xArr[4] + xArr[0]*xArr[3]*xArr[4] + xArr[2]*xArr[3]*xArr[4]);
  b21 = -1.0*(xArr[0]*xArr[1]*xArr[3] + xArr[0]*xArr[1]*xArr[4] + xArr[0]*xArr[3]*xArr[4] + xArr[1]*xArr[3]*xArr[4]);
  b31 = -1.0*(xArr[0]*xArr[1]*xArr[2] + xArr[0]*xArr[1]*xArr[4] + xArr[0]*xArr[2]*xArr[4] + xArr[1]*xArr[2]*xArr[4]);
  b41 = -1.0*(xArr[0]*xArr[1]*xArr[2] + xArr[0]*xArr[1]*xArr[3] + xArr[0]*xArr[2]*xArr[3] + xArr[1]*xArr[2]*xArr[3]);

  b02 = xArr[1]*xArr[2] + xArr[1]*xArr[3] + xArr[2]*xArr[3] + xArr[1]*xArr[4] + xArr[2]*xArr[4] + xArr[3]*xArr[4];
  b12 = xArr[0]*xArr[2] + xArr[0]*xArr[3] + xArr[2]*xArr[3] + xArr[0]*xArr[4] + xArr[2]*xArr[4] + xArr[3]*xArr[4];
  b22 = xArr[0]*xArr[1] + xArr[0]*xArr[3] + xArr[1]*xArr[3] + xArr[0]*xArr[4] + xArr[1]*xArr[4] + xArr[3]*xArr[4];
  b32 = xArr[0]*xArr[1] + xArr[0]*xArr[2] + xArr[1]*xArr[2] + xArr[0]*xArr[4] + xArr[1]*xArr[4] + xArr[2]*xArr[4];
  b42 = xArr[0]*xArr[1] + xArr[0]*xArr[2] + xArr[1]*xArr[2] + xArr[0]*xArr[3] + xArr[1]*xArr[3] + xArr[2]*xArr[3];

  b03 = -1.0*(xArr[1] + xArr[2] + xArr[3] + xArr[4]);
  b13 = -1.0*(xArr[0] + xArr[2] + xArr[3] + xArr[4]);
  b23 = -1.0*(xArr[0] + xArr[1] + xArr[3] + xArr[4]);
  b33 = -1.0*(xArr[0] + xArr[1] + xArr[2] + xArr[4]);
  b43 = -1.0*(xArr[0] + xArr[1] + xArr[2] + xArr[3]);

  b04 = 1.0;
  b14 = 1.0;
  b24 = 1.0;
  b34 = 1.0;
  b44 = 1.0;

  g0 = a0*b00 + a1*b10 + a2*b20 + a3*b30 + a4*b40;
  g1 = a0*b01 + a1*b11 + a2*b21 + a3*b31 + a4*b41;
  g2 = a0*b02 + a1*b12 + a2*b22 + a3*b32 + a4*b42;
  g3 = a0*b03 + a1*b13 + a2*b23 + a3*b33 + a4*b43;
  g4 = a0*b04 + a1*b14 + a2*b24 + a3*b34 + a4*b44;

  val = g0 + g1*(pow(myX,1.0));
  val += g2*(pow(myX,2.0)) + g3*(pow(myX,3.0));
  val += g4*(pow(myX,4.0));

  Real OldX;
  if (myX < xArr[2]) OldX = xArr[1];
  if (myX > xArr[2]) OldX = xArr[3];
  
 


  val = g0*(myX-OldX) + g1/2.0*(pow(myX,2.0)-pow(OldX,2.0));
  val += g2/3.0*(pow(myX,3.0)-pow(OldX,3.0)) + g3/4.0*(pow(myX,4.0)-pow(OldX,4.0));
  val += g4/5.0*(pow(myX,5.0)-pow(OldX,5.0));
  val *= 1/(myX-OldX); 


  return val; 

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

  double xArrL[5], xArrR[5];
  double primsL[4], primsR[4];
  double primsC[5], xArrC[5];
  double MyxL, MyxR;
  double valLn, valRn, valLp1, valRp1;
  double coeffL[4], coeffR[4];

  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
    // cache the x1-sliced primitive states for eigensystem calculation
    for (int n=0; n<(NHYDRO); ++n) {
#pragma omp simd
      for (int i=il-1; i<=iu; ++i) {
        //get i-1/2 profile
        xArrL[0] = pmb->pcoord->x1f(i-2) - pmb->pcoord->x1f(i);
        xArrL[1] = pmb->pcoord->x1f(i-1) - pmb->pcoord->x1f(i);
        xArrL[2] = 0.0;
        xArrL[3] = pmb->pcoord->x1f(i+1) - pmb->pcoord->x1f(i);
        xArrL[4] = pmb->pcoord->x1f(i+2) - pmb->pcoord->x1f(i);
        
        primsL[0] = w(n,k,j,i-2);
        primsL[1] = w(n,k,j,i-1);
        primsL[2] = w(n,k,j,i  );
        primsL[3] = w(n,k,j,i+1);
        WallIntegration(xArrL,primsL,5,coeffL);
     
        //Get i+1/2 profile
        xArrR[0] = pmb->pcoord->x1f(i-1) - pmb->pcoord->x1f(i+1);
        xArrR[1] = pmb->pcoord->x1f(i)   - pmb->pcoord->x1f(i+1);
        xArrR[2] = 0.0;
        xArrR[3] = pmb->pcoord->x1f(i+2) - pmb->pcoord->x1f(i+1);
        xArrR[4] = pmb->pcoord->x1f(i+3) - pmb->pcoord->x1f(i+1);
            
        primsR[0] = w(n,k,j,i-1);
        primsR[1] = w(n,k,j,i  );
        primsR[2] = w(n,k,j,i+1);
        primsR[3] = w(n,k,j,i+2);
        WallIntegration(xArrR,primsR,5,coeffR);     

        // Get base line L R states
        valLn = coeffL[0];
        valRn = coeffR[0];

        if (valLn > std::max(primsL[2],primsL[1])) valLn = std::max(primsL[2],primsL[1]);
        if (valLn < std::min(primsL[2],primsL[1])) valLn = std::min(primsL[2],primsL[1]);

        if (valRn > std::max(primsR[2],primsR[1])) valRn = std::max(primsR[2],primsR[1]);
        if (valRn < std::min(primsR[2],primsR[1])) valRn = std::min(primsR[2],primsR[1]);

        //Extrema Detector
        Real dqm = w(n,k,j,i) - valLn;
        Real dqp = valRn - w(n,k,j,i);
      
        if (dqp*dqm <= 0.0){
          valLn = w(n,k,j,i);
          valRn = w(n,k,j,i);
        } else {
          // Overshoot i-1/2,R / i,(-) state
          if (fabs(dqm) >= pmb->precon->hplus_ratio_i(i)*fabs(dqp)) {
            valLn = w(n,k,j,i) - pmb->precon->hplus_ratio_i(i)*dqp;
          }
          // Overshoot i+1/2,L / i,(+) state
          if (fabs(dqp) >= pmb->precon->hminus_ratio_i(i)*fabs(dqm)) {
            valRn = w(n,k,j,i) + pmb->precon->hminus_ratio_i(i)*dqm;
          }
        }

        xArrC[1] = pmb->pcoord->x1f(i)   - pmb->pcoord->x1v(i);
        xArrC[2] = 0.0;
        xArrC[3] = pmb->pcoord->x1f(i+1) - pmb->pcoord->x1v(i);
        
        primsC[1] = valLn;
        primsC[2] = w(n,k,j,i);
        primsC[3] = valRn;

        Real A, B, C;
        A = primsC[1]/(xArrC[1]*(xArrC[1]-xArrC[3]));
        B = primsC[2]/(xArrC[1]*xArrC[3]);
        C = primsC[3]/(xArrC[3]*(xArrC[3]-xArrC[1]));

        Real c1, c2, c0;
        c2 = A + B + C;
        c1 = -1.0*(A*(xArrC[3])+B*(xArrC[1]+xArrC[3])+C*xArrC[1]);
        c0 = B*xArrC[1]*xArrC[3];
        

        // Get next time step states
        MyxL = xArrC[1] + abs(v1f(i))*dt;	
        MyxR = xArrC[3] - abs(v1f(i+1))*dt;	

        // Get base line L R states
        valLp1 = ( c2*(pow(MyxL,3.0)-pow(xArrC[1],3.0))/3.0 
                 + c1*(pow(MyxL,2.0)-pow(xArrC[1],2.0))/2.0 
                 + c0*(pow(MyxL,1.0)-pow(xArrC[1],1.0))) / (MyxL - xArrC[1]) ;
        valRp1 = ( c2*(pow(MyxR,3.0)-pow(xArrC[3],3.0))/3.0 
                 + c1*(pow(MyxR,2.0)-pow(xArrC[3],2.0))/2.0 
                 + c0*(pow(MyxR,1.0)-pow(xArrC[3],1.0))) / (MyxR - xArrC[3]) ;

        if (valLp1 > std::max(primsC[2],primsC[1])) valLp1 = std::max(primsC[2],primsC[1]);
        if (valLp1 < std::min(primsC[2],primsC[1])) valLp1 = std::min(primsC[2],primsC[1]);
        if (valRp1 > std::max(primsC[2],primsC[3])) valRp1 = std::max(primsC[2],primsC[3]);
        if (valRp1 < std::min(primsC[2],primsC[3])) valRp1 = std::min(primsC[2],primsC[3]);

        if (MyxR==xArrC[3]) valRp1 = valRn;
        if (MyxL==xArrC[1]) valLp1 = valLn;
        wL(n,k,j,i+1) = valRp1;
        wR(n,k,j,i  ) = valLp1;
      }
    }
  }}

  return;
}







