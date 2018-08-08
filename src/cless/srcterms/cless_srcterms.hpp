#ifndef CLESS_SRCTERMS_CLESS_SRCTERMS_HPP_
#define CLESS_SRCTERMS_CLESS_SRCTERMS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cless_srcterms.hpp
//  \brief defines class HydroSourceTerms
//  Contains data and functions that implement physical (not coordinate) source terms

// Athena headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"

// Forward declarations
class Cless;
class ParameterInput;

//! \class HydroSourceTerms
//  \brief data and functions for physical source terms in the hydro

class ClessSourceTerms {
public:
  ClessSourceTerms(Cless *pcl, ParameterInput *pin);
  ~ClessSourceTerms();

  // data
  bool cless_sourceterms_defined;

  // functions
  void AddClessSourceTerms(const Real time, const Real dt, const AthenaArray<Real> *flx,
    const AthenaArray<Real> &p, AthenaArray<Real> &c);
  void SelfGravityCL(const Real dt, const AthenaArray<Real> *flx,
    const AthenaArray<Real> &p, AthenaArray<Real> &c);
  //void EnrollSrcTermFunction(SrcTermFunc_t my_func);
  //SrcTermFunc_t UserSourceTermCL;

private:
  Cless *pmy_cless_;  // ptr to Cless containing this ClessSourceTerms
};
#endif // CLESS_SRCTERMS_CLESS_SRCTERMS_HPP_
