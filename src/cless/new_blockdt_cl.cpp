//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file new_blockdt.cpp
//  \brief computes timestep using CFL condition on a MEshBlock

// C/C++ headers
#include <algorithm>  // min()
#include <cfloat>     // FLT_MAX
#include <cmath>      // fabs(), sqrt()

// Athena++ headers
#include "cless.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../cless/cless.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
// \!fn Real Cless::NewBlockTimeStep(void)
// \brief calculate the minimum timestep within a MeshBlock
// only used in CLESS_ONLY mode 

Real Cless::NewBlockTimeStepCL(void) {
  MeshBlock *pmb=pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  AthenaArray<Real> w,bcc,b_x1f,b_x2f,b_x3f, wcl;

	if (CLESS_ENABLED) {
		wcl.InitWithShallowCopy(pmb->pcless->w); 
	}

  AthenaArray<Real> dt1, dt2, dt3;
  dt1.InitWithShallowCopy(dt1_);
  dt2.InitWithShallowCopy(dt2_);
  dt3.InitWithShallowCopy(dt3_);
  Real wi[(NWAVE+NINT)];
	Real wicl[(NWAVECL)]; 

  Real min_dt = (FLT_MAX);

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      pmb->pcoord->CenterWidth1(k,j,is,ie,dt1);
      pmb->pcoord->CenterWidth2(k,j,is,ie,dt2);
      pmb->pcoord->CenterWidth3(k,j,is,ie,dt3);
#pragma ivdep
			for (int i=is; i<=ie; ++i) {
				Real c1f, c2f, c3f; 
				wicl[IDN ] = wcl(IDN ,k,j,i);
				wicl[IVX ] = wcl(IVX ,k,j,i);
				wicl[IVY ] = wcl(IVY ,k,j,i);
				wicl[IVZ ] = wcl(IVZ ,k,j,i);
				wicl[IP11] = wcl(IP11,k,j,i);
				wicl[IP22] = wcl(IP22,k,j,i);
				wicl[IP33] = wcl(IP33,k,j,i);
				wicl[IP12] = wcl(IP12,k,j,i);
				wicl[IP13] = wcl(IP13,k,j,i);
				wicl[IP23] = wcl(IP23,k,j,i);

				pmb->peos->SoundSpeedsCL(wicl,&c1f,&c2f,&c3f); 

				dt1(i) /= fabs(wicl[IVX] + c1f);
				dt2(i) /= fabs(wicl[IVY] + c2f);
				dt3(i) /= fabs(wicl[IVZ] + c3f); 
			}

      // compute minimum of (v1 +/- C)
      for (int i=is; i<=ie; ++i) {
        Real& dt_1 = dt1(i);
        min_dt = std::min(min_dt,dt_1);
      }

      // if grid is 2D/3D, compute minimum of (v2 +/- C)
      if (pmb->block_size.nx2 > 1) {
        for (int i=is; i<=ie; ++i) {
          Real& dt_2 = dt2(i);
          min_dt = std::min(min_dt,dt_2);
        }
      }

      // if grid is 3D, compute minimum of (v3 +/- C)
      if (pmb->block_size.nx3 > 1) {
        for (int i=is; i<=ie; ++i) {
          Real& dt_3 = dt3(i);
          min_dt = std::min(min_dt,dt_3);
        }
      }

    }
  }

  min_dt *= pmb->pmy_mesh->cfl_number;

  //if (UserTimeStep_!=NULL) {
  //  min_dt = std::min(min_dt, UserTimeStep_(pmb));
  //}

  pmb->new_block_dt=min_dt;
  return min_dt;
}
