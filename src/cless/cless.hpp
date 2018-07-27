#ifndef CLESS_CLESS_HPP_
#define CLESS_CLESS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cless.hpp
//  \brief definitions for Hydro class

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
//#include "../task_list/task_list.hpp"

class MeshBlock;
class ParameterInput;
class ClessSourceTerms;
struct IntegratorWeight;

//! \class Cless 
//  \brief cless data and functions

class Cless {
//friend class Hydro;
public:
  Cless(MeshBlock *pmb, ParameterInput *pin);
  ~Cless();

  // data
  MeshBlock* pmy_block;    // ptr to MeshBlock containing this Cless 
  // conserved and primitive variables
  AthenaArray<Real> u,w;      // time-integrator memory register #1
  AthenaArray<Real> u1,w1;    // time-integrator memory register #2
  AthenaArray<Real> u2;       // time-integrator memory register #3
  AthenaArray<Real> flux[3];  // face-averaged flux vector

  ClessSourceTerms *psrc;
  
	// functions
	Real NewBlockTimeStepCL(void); // computes new timestep for CLESS_ONLY mode 
  void WeightedAveUCL(AthenaArray<Real> &u_out, AthenaArray<Real> &u_in1,
    AthenaArray<Real> &u_in2, const Real wght[3]);
  void AddFluxDivergenceToAverageCL(AthenaArray<Real> &w,
    const Real wght, AthenaArray<Real> &u_out);
  void CalculateFluxesCL(AthenaArray<Real> &w, int order);
  void RiemannSolverCL(const int kl, const int ku, const int jl, const int ju,
    const int il, const int iu, const int ivx, const int ip12, 
    AthenaArray<Real> &wl, AthenaArray<Real> &wr, AthenaArray<Real> &flx);

  void AddGravityFluxCL(void);
  void AddGravityFluxWithGflxCL(void);
  void CalculateGravityFluxCL(AthenaArray<Real> &phi_in);
  void CorrectGravityFluxCL(void);

private:
	AthenaArray<Real> dt1_,dt2_,dt3_; 
  // scratch space used to compute fluxes
  AthenaArray<Real> wl_, wr_;
  AthenaArray<Real> dxw_;
  AthenaArray<Real> x1face_area_, x2face_area_, x3face_area_;
  AthenaArray<Real> x2face_area_p1_, x3face_area_p1_;
  AthenaArray<Real> cell_volume_;
  AthenaArray<Real> dflx_;
  // self-gravity
  AthenaArray<Real> gflx[3], gflx_old[3]; // gravity tensor (old Athena style)
};
#endif // CLESS_CLESS_HPP_
