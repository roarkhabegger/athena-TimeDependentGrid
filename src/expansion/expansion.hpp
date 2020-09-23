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

  // Expansion Data
  AthenaArray<Real> v1f;
  AthenaArray<Real> v2f;
  AthenaArray<Real> v3f;
  //Integration Registers
  AthenaArray<Real> x1_0, x2_0, x3_0;
  AthenaArray<Real> x1_1, x2_1, x3_1;
  AthenaArray<Real> x1_2, x2_2, x3_2;

  Real mydt;
  int il, iu, jl, ju, kl, ku, ng; //With Ghost cells
  int ie,is,je,js,ke,ks; //Without ghost cells

  void WeightedAveX(const int low, const int up, AthenaArray<Real> &x_out, AthenaArray<Real> &x_in1, AthenaArray<Real> &x_in2, const Real wght[3]);
  void AddWallFluxDivergence( Real dt, AthenaArray<Real> &prim, AthenaArray<Real> &cons);
  void IntegrateWalls(Real dt);

  void ExpansionSourceTerms(const Real dt, const AthenaArray<Real> *flx, 
                      const AthenaArray<Real> &p, AthenaArray<Real> &c); 

  void GridEdit(MeshBlock *pmb);
  void UpdateVelData(MeshBlock *pmb,Real time, Real dt);
  Real GridTimeStep(MeshBlock *pmb);
  void UpdateMeshSize(MeshBlock *pmb);

  void InterpData( const int myn, const int myk, const int myj, const int myi, const double dt,
                   const AthenaArray<Real> grid, const AthenaArray<Real> gridVel,
                   const AthenaArray<Real> data, AthenaArray<Real> &output); 
private:
  MeshBlock* pmy_block;    // ptr to MeshBlock containing this Expansion

};
#endif // EXPANSION_EXPANSION_HPP_
