#ifndef EXPANSION_EXPANSION_HPP_
#define EXPANSION_EXPANSION_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file expansion.hpp
//  \brief definitions for Expansion class

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
//#include "../task_list/task_list.hpp"

class MeshBlock;
class ParameterInput;
struct IntegratorWeight;

//! \class Expansion
//  \brief Expanding Grid information and edits

class Expansion {
//friend class Field;
friend class Hydro;
friend class Mesh;
friend class Reconstruction; 
public:
  Expansion(MeshBlock *pmb, ParameterInput *pin);
  ~Expansion();

  //Boolean Direction variables
  bool x1Move;
  bool x2Move;
  bool x3Move;


  // Expansion Data
  //AthenaArray<Real> v1f;
  //AthenaArray<Real> v2f;
  //AthenaArray<Real> v3f;
  AthenaArray<Real> vol;

  //Integration Registers
  AthenaArray<Real> x1_0, x2_0, x3_0;
  AthenaArray<Real> x1_1, x2_1, x3_1;
  AthenaArray<Real> x1_2, x2_2, x3_2;
  
  AthenaArray<Real> ExpwL, ExpwR;
  AthenaArray<Real> expFlux[3];  // face-averaged flux vector
  AthenaArray<Real> vf[3];  // face-averaged flux vector

  Real mydt;
  int il, iu, jl, ju, kl, ku, ng; //With Ghost cells
  int ie,is,je,js,ke,ks; //Without ghost cells

  void WeightedAveX(const int low, const int up, AthenaArray<Real> &x_out, AthenaArray<Real> &x_in1, AthenaArray<Real> &x_in2, const Real wght[3]);
  //void AddWallFluxDivergence( Real dt, AthenaArray<Real> &prim, AthenaArray<Real> &cons);
  void IntegrateWalls(Real dt);

  void ExpansionSourceTerms(const Real dt, const AthenaArray<Real> *flux, const AthenaArray<Real> &prim, AthenaArray<Real> &cons); 

  void GridEdit(MeshBlock *pmb, bool lastStage);
  void UpdateVelData(MeshBlock *pmb,Real time, Real dt);
  Real GridTimeStep(MeshBlock *pmb);
  //void UpdateMeshSize(MeshBlock *pmb);


  //void WallFlux(const int k, const int j, const int i, const int ivx,
  //  const AthenaArray<Real> &bx, AthenaArray<Real> &e1, AthenaArray<Real> &e2,
  //  AthenaArray<Real> &vArr, const int dir);
  //void AddWallFlux(const int k, const int j, const int i, const int dir, Real dt, AthenaArray<Real> &cons);
//  void PPMOffsetX1(MeshBlock *pmb, const int n,
//    const int k, const int j, const int i, 
//    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc);
//
  void PPMOffsetX2(MeshBlock *pmb, const int n,
    const int k, const int j, const int i, 
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc);

  void PPMOffsetX3(MeshBlock *pmb, const int n,
    const int k, const int j, const int i, 
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc);
private:
  MeshBlock* pmy_block;    // ptr to MeshBlock containing this Expansion

  //void WallIntegrationParabolic(AthenaArray<Real> &xArr, AthenaArray<Real> &prims, AthenaArray<Real> &coeff);
  //Real CenterIntegration(AthenaArray<Real> &xArr, AthenaArray<Real> &prims, Real myX, int xArrSize);

  //Reconstruction variables
  AthenaArray<Real> xArr, prims, coeffL, coeffR;
  AthenaArray<Real> dxArr, intPnts;

  Real a1, a2, a3, a4;
  Real b11, b12, b13, b14;
  Real b21, b22, b23, b24;
  Real b31, b32, b33, b34;
  Real g4, g3, g2, g1;

  //Flux Variables
  AthenaArray<Real> wi;


};
#endif // EXPANSION_EXPANSION_HPP_
