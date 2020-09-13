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


Real WallIntegration(double xArr[],double prims[],Real myX, int xArrSize);
Real CenterIntegration(double xArr[],double prims[], Real myX, int xArrSize);

Real WallIntegration(double xArr[], double prims[], Real myX, int xArrSize) {
  // Prims are cell centered, xArr has cell wall locations 
  // Prims should be length one less than xArr
  int n  = xArrSize;
  double dxArr[n-1];
  double intPnts[n];
  double terms[n];
  int m,p,l;
  double xm,xp,xl;
  double NewTerm;
  double val;
  
  for (m=0;m<=(n-2);++m) { 
     dxArr[m]=xArr[m+1]-xArr[m];
  }

  for (m=0;m<=(n-1);++m) {
    intPnts[m] = 0.0;
    terms[m]   = 0.0;
  }


  for (m=1;m<=(n-1);++m) {
    for (p=m;p<=(n-1);++p){
      intPnts[p] += dxArr[m-1]*prims[m-1];
    }
  }

  for (m=0;m<=(n-1);++m){
    xm = xArr[m];
    for (p=0;p<=(n-1);++p) {
      xp = xArr[p];
      NewTerm = 0.0;
      if (xp != xm) {
        NewTerm = 1.0/(xm-xp);
        for (l=0;l<=(n-1);++l) {
          xl = xArr[l];
          if ((xl != xm) && (xl!=xp)) {
            NewTerm *=(myX-xl)/(xm-xl);
          }
        }
      }     
      terms[m] += NewTerm;
    }
  } //Term Arr Complete!

  val = 0.0;
  for (m=0;m<=(n-1);++m) {
    val += terms[m]*intPnts[m];
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
void  Expansion::SrcTermDataX1( const int myn, const int myk, const int myj, const int myi, const AthenaArray<Real> grid, const AthenaArray<Real> gridDel, const AthenaArray<Real> data, AthenaArray<Real> &output) {

  double xArrL[7], xArrC[7], xArrR[7];
  double primsL[6], primsC[7], primsR[6];
  double MyxL, MyxR, MyxC;
  double valLi, valRi, valCi;

    
  xArrL[0] = grid(myi-2);
  xArrL[1] = grid(myi-1);
  xArrL[2] = grid(myi);
  xArrL[3] = grid(myi+1);
  xArrL[4] = grid(myi+2);
  //xArrL[5] = grid(myi+2);
  //xArrL[6] = grid(myi+3);  

  xArrR[0] = grid(myi-1);
  xArrR[1] = grid(myi);
  xArrR[2] = grid(myi+1);
  xArrR[3] = grid(myi+2);
  xArrR[4] = grid(myi+3);
  //xArrR[5] = grid(myi+3);
  //xArrR[6] = grid(myi+4);


  for (int v=0; v<=4;++v) {
    xArrC[v] = 0.5*(xArrL[v]+xArrR[v]);
  }

  primsL[0] = data(myn,myk,myj,myi-2);   
  primsL[1] = data(myn,myk,myj,myi-1);
  primsL[2] = data(myn,myk,myj,myi);
  primsL[3] = data(myn,myk,myj,myi+1);
  //primsL[4] = data(myn,myk,myj,myi+1);
  //primsL[5] = data(myn,myk,myj,myi+2);

  primsR[0] = data(myn,myk,myj,myi-1);
  primsR[1] = data(myn,myk,myj,myi); 
  primsR[2] = data(myn,myk,myj,myi+1);
  primsR[3] = data(myn,myk,myj,myi+2);        
  //primsR[4] = data(myn,myk,myj,myi+2);        
  //primsR[5] = data(myn,myk,myj,myi+3);        
 
  primsC[0] = data(myn,myk,myj,myi-2);
  primsC[1] = data(myn,myk,myj,myi-1);
  primsC[2] = data(myn,myk,myj,myi);
  primsC[3] = data(myn,myk,myj,myi+1);
  primsC[4] = data(myn,myk,myj,myi+2);
  //primsC[5] = data(myn,myk,myj,myi+2);
  //primsC[6] = data(myn,myk,myj,myi+3);
        
  MyxL = xArrL[2];	
  MyxR = xArrR[2];        
  MyxC = xArrC[2];
  valLi =  WallIntegration(xArrL,primsL,MyxL,5);
  valRi = WallIntegration(xArrR,primsR,MyxR,5);
  valCi = primsC[2];
  output(0) = valLi; // 0.5*(valLi+ppmLi);
  output(1) = valCi;
  output(2) = valRi;

  MyxL = xArrL[2]+gridDel(myi);	
  MyxR = xArrR[2]+gridDel(myi+1);    
  MyxC = xArrC[2] + 0.5*(gridDel(myi)+gridDel(myi+1));    
  valLi =  WallIntegration(xArrL,primsL,MyxL,5);
  valRi = WallIntegration(xArrR,primsR,MyxR,5);
  valCi = CenterIntegration(xArrC, primsC,MyxC,5);
  output(3) = valLi; // 0.5*(valLi+ppmLi);
  output(4) = valCi;
  output(5) = valRi; 
  return;
}


void  Expansion::FlxCenterX1( const int myn, const int myk, const int myj, const int myi,
                                         const AthenaArray<Real> grid, const AthenaArray<Real> gridDel,
                                         const AthenaArray<Real> data, Real &output, const int dir) {

  double xArrC[5];
  double primsC[5];
  double termsC[5];
  double Myx;
  int m,p,l;
  double valCi;

    xArrC[0] = 0.5*(grid(myi-2)+grid(myi-1));
    xArrC[1] = 0.5*(grid(myi-1)+grid(myi));
    xArrC[2] = 0.5*(grid(myi)+grid(myi+1));
    xArrC[3] = 0.5*(grid(myi+1)+grid(myi+2));
    xArrC[4] = 0.5*(grid(myi+2)+grid(myi+3));
    Myx = xArrC[2]+0.25*(gridDel(myi)+gridDel(myi+1));	

  //std::cout << xArrC << std::endl;
  primsC[0] = data(myn,myk,myj,myi-2);   
  primsC[1] = data(myn,myk,myj,myi-1);
  primsC[2] = data(myn,myk,myj,myi  );
  primsC[3] = data(myn,myk,myj,myi+1);
  primsC[4] = data(myn,myk,myj,myi+2);
  valCi = 0.0;
  valCi = CenterIntegration(xArrC, primsC,Myx,5);
  output = valCi; 
  return;
}

