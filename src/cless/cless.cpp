//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cless.cpp
//  \brief implementation of functions in class Cless 

// C/C++ headers
#include <string>

// Athena++ headers
#include "cless.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "srcterms/cless_srcterms.hpp"

// constructor, initializes data structures and parameters

Cless::Cless(MeshBlock *pmb, ParameterInput *pin) {
  pmy_block = pmb;

  // Allocate memory for primitive/conserved variables
  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmy_block->block_size.nx2 > 1) ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if (pmy_block->block_size.nx3 > 1) ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);

  // Allocate memory registers for primitive/conserved variables for time-integrator
	if (CLESS_ENABLED) {
		u.NewAthenaArray(NCLESS,ncells3,ncells2,ncells1);
		w.NewAthenaArray(NCLESS,ncells3,ncells2,ncells1);
		u1.NewAthenaArray(NCLESS,ncells3,ncells2,ncells1);
		w1.NewAthenaArray(NCLESS,ncells3,ncells2,ncells1);
		// If user-requested time integrator is type 3S*, allocate additional memory registers
		std::string integrator = pin->GetOrAddString("time","integrator","vl2");
		if (integrator == "ssprk5_4")
			// future extension may add "int nregister" to Hydro class
			u2.NewAthenaArray(NCLESS,ncells3,ncells2,ncells1);

		flux[X1DIR].NewAthenaArray(NCLESS,ncells3,ncells2,ncells1+1);
		if (pmy_block->block_size.nx2 > 1)
			flux[X2DIR].NewAthenaArray(NCLESS,ncells3,ncells2+1,ncells1);
		if (pmy_block->block_size.nx3 > 1)
			flux[X3DIR].NewAthenaArray(NCLESS,ncells3+1,ncells2,ncells1);

		// Allocate memory for scratch arrays
		dt1_.NewAthenaArray(ncells1);
		dt2_.NewAthenaArray(ncells1);
		dt3_.NewAthenaArray(ncells1); 
		dxw_.NewAthenaArray(ncells1);
		wl_.NewAthenaArray((NWAVECL),ncells3,ncells2,ncells1);
		wr_.NewAthenaArray((NWAVECL),ncells3,ncells2,ncells1);
		x1face_area_.NewAthenaArray(ncells1+1);
		if (pmy_block->block_size.nx2 > 1) {
			x2face_area_.NewAthenaArray(ncells1);
			x2face_area_p1_.NewAthenaArray(ncells1);
		}
		if (pmy_block->block_size.nx3 > 1) {
			x3face_area_.NewAthenaArray(ncells1);
			x3face_area_p1_.NewAthenaArray(ncells1);
		}
		cell_volume_.NewAthenaArray(ncells1);
		dflx_.NewAthenaArray((NCLESS),ncells1);
		// for one-time potential calcuation and correction (old Athena)
		if (SELF_GRAVITY_ENABLED == 3) {
			gflx[X1DIR].NewAthenaArray(NCLESS,ncells3,ncells2,ncells1+1);
			if (pmy_block->block_size.nx2 > 1)
				gflx[X2DIR].NewAthenaArray(NCLESS,ncells3,ncells2+1,ncells1);
			if (pmy_block->block_size.nx3 > 1)
				gflx[X3DIR].NewAthenaArray(NCLESS,ncells3+1,ncells2,ncells1);

			gflx_old[X1DIR].NewAthenaArray(NCLESS,ncells3,ncells2,ncells1+1);
			if (pmy_block->block_size.nx2 > 1)
				gflx_old[X2DIR].NewAthenaArray(NCLESS,ncells3,ncells2+1,ncells1);
			if (pmy_block->block_size.nx3 > 1)
				gflx_old[X3DIR].NewAthenaArray(NCLESS,ncells3+1,ncells2,ncells1);
		}

		// Construct ptrs to objects of various classes needed to integrate hydro/MHD eqns
		psrc  = new ClessSourceTerms(this,pin);
  }
}

// destructor

Cless::~Cless() {
  u.DeleteAthenaArray();
  w.DeleteAthenaArray();
  u1.DeleteAthenaArray();
  w1.DeleteAthenaArray();
  // only allocated if integrator was 3S* integrator
  u2.DeleteAthenaArray();

  flux[X1DIR].DeleteAthenaArray();
  if (pmy_block->block_size.nx2 > 1) flux[X2DIR].DeleteAthenaArray();
  if (pmy_block->block_size.nx3 > 1) flux[X3DIR].DeleteAthenaArray();

	dt1_.DeleteAthenaArray();
	dt2_.DeleteAthenaArray();
	dt3_.DeleteAthenaArray();
  dxw_.DeleteAthenaArray();
  wl_.DeleteAthenaArray();
  wr_.DeleteAthenaArray();
  x1face_area_.DeleteAthenaArray();
  if (pmy_block->block_size.nx2 > 1) {
    x2face_area_.DeleteAthenaArray();
    x2face_area_p1_.DeleteAthenaArray();
  }
  if (pmy_block->block_size.nx3 > 1) {
    x3face_area_.DeleteAthenaArray();
    x3face_area_p1_.DeleteAthenaArray();
  }
  cell_volume_.DeleteAthenaArray();
  dflx_.DeleteAthenaArray();
  // for one-time potential calcuation and correction (old Athena)
  if (SELF_GRAVITY_ENABLED == 3) {
    gflx[X1DIR].DeleteAthenaArray();
    if (pmy_block->block_size.nx2 > 1) gflx[X2DIR].DeleteAthenaArray();
    if (pmy_block->block_size.nx3 > 1) gflx[X3DIR].DeleteAthenaArray();
    gflx_old[X1DIR].DeleteAthenaArray();
    if (pmy_block->block_size.nx2 > 1) gflx_old[X2DIR].DeleteAthenaArray();
    if (pmy_block->block_size.nx3 > 1) gflx_old[X3DIR].DeleteAthenaArray();
  }
  delete psrc;
}
