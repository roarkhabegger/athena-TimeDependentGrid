//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file knovae.cpp
//  \brief Problem generator for spherical knova  problem.  
//
//

// C++ headers
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <iostream>
#include <iomanip>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../bvals/bvals.hpp"

//====================================================================================
// global variables

#ifdef MPI_PARALLEL
typedef struct MPI_Comm_Sub {
  MPI_Group gsub, gworld;
  MPI_Comm  comsub;
} MPI_Comm_Sub;

MPI_Comm_Sub comm_x1;
MPI_Comm_Sub comm_slab;
#endif

//====================================================================================
// local functions
Real LogMeshSpacingX1(Real x, RegionSize rs);

void ReflectInnerX1_nonuniform(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

//====================================================================================
// Enroll user-specific functions
void Mesh::InitUserMeshData(ParameterInput *pin) {
  Real x1rat = pin->GetOrAddReal("mesh","x1rat",0.0);
  
  if (x1rat < 0.0) {
    EnrollUserMeshGenerator(X1DIR, LogMeshSpacingX1);
    if (mesh_bcs[INNER_X1] == GetBoundaryFlag("user")) {
      EnrollUserBoundaryFunction(INNER_X1, ReflectInnerX1_nonuniform);
    }
  }
  return;
}

//====================================================================================
// "Logarithmic" (power-law) mesh
// Note that the grid setup asks for the **local** x1min, x1max etc. 
// x is the "logical" position in the grid, with the logical grid runing
// from 0 to 1, i.e. x = i/Nx
Real LogMeshSpacingX1(Real x, RegionSize rs) {
  Real xf, xrat;
  xrat   = pow(rs.x1max/rs.x1min,1.0/((Real) rs.nx1)); // Only valid for fixed grid, no MPI
  xf     = rs.x1min*pow(xrat,x*rs.nx1); // x = i/Nx
  return xf;
}

//========================================================================================
// Reflecting inner X1 boundary conditions for radially non-uniform grids

void ReflectInnerX1_nonuniform(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  // copy hydro variables into ghost zones, reflecting v1
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVX)) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=1; i<=ngh; ++i) {
            prim(IVX,k,j,is-i) = -prim(IVX,k,j,(is+i-1));  // reflect 1-velocity
          }
        }
      }
    } else {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=1; i<=ngh; ++i) {
            prim(n,k,j,is-i) = prim(n,k,j,(is+i-1));
          }
        }
      }
    }
  }
  // copy face-centered magnetic fields into ghost zones, reflecting b1
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x1f(k,j,(is-i)) = -b.x1f(k,j,(is+i  ));  // reflect 1-field
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x2f(k,j,(is-i)) =  b.x2f(k,j,(is+i-1));
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x3f(k,j,(is-i)) =  b.x3f(k,j,(is+i-1));
        }
      }
    }
  }
  return;
}

//========================================================================================
//! \fn void MeshBlock::InitOTFOutput(ParameterInput *pin)
//  \brief Sets data structures etc for on-the-fly analysis.
//========================================================================================

void MeshBlock::InitOTFOutput(ParameterInput *pin) {
  otf_data.len = 8*(NSCALARS+1)+3+8*pmy_mesh->mesh_size.nx2*pmy_mesh->mesh_size.nx3;
  otf_data.data = new Real[otf_data.len];
  //std::cout << "rank=" << Globals::my_rank << " loc.level=" << loc.level << " lid=" << lid << std::endl;
  if (lid == 0) { // only defined on root level
    // First, we'll do the radial communicator for radial sums.
    // count ranks in local loc.lx2,loc.lx3 column. Complicated way to express nrbx1.
    int nbx1count=0;
    for (int nbt=0; nbt<pmy_mesh->nbtotal; nbt++) {
      if (  (loc.lx2   == pmy_mesh->loclist[nbt].lx2)
          &&(loc.lx3   == pmy_mesh->loclist[nbt].lx3)
          &&(loc.level == pmy_mesh->loclist[nbt].level)) {
        nbx1count++;
      }
    }
    int* bx1ranks = new int[nbx1count]; // rank indices along one column (duplicates possible)
    int* bx1block = new int[nbx1count]; // block indices along one column (no duplicates)
    // again, now we assign the ranks
    int ibx1count = 0; // use as index
    for (int nbt=0; nbt<pmy_mesh->nbtotal; nbt++) {
      if (  (loc.lx2   == pmy_mesh->loclist[nbt].lx2)
          &&(loc.lx3   == pmy_mesh->loclist[nbt].lx3)
          &&(loc.level == pmy_mesh->loclist[nbt].level)) {
        bx1ranks[ibx1count] = pmy_mesh->ranklist[nbt];
        bx1block[ibx1count] = nbt;
        ibx1count++;
      }
    }
    //for (int nbt=0; nbt<nbx1count; nbt++) {
    //  std::cout << "block=" << bx1block[nbt] << " rank=" << bx1ranks[nbt] << " lx2=" << loc.lx2 << " lx3=" << loc.lx3 << std::endl;
    //}  
    // At this point, each rank knows which blocks to sum over.
    // Now we need to clean the rank list (removing duplicates).
    int nrx1count = 0;
    int* rcount = new int[Globals::nranks];
    for (int nr=0; nr<Globals::nranks; nr++) {
      rcount[nr] = 0;
      for (int nbt=0; nbt<nbx1count; nbt++) {
        if (bx1ranks[nbt] == nr) {
          rcount[nr] = 1;
        }
      }
      nrx1count += rcount[nr];
    }
    int irx1count = 0;
    int* rx1ranks = new int[nrx1count];
    for (int nr=0; nr<Globals::nranks; nr++) {
      if (rcount[nr] == 1) { // if rank is in global list for this column...
        rx1ranks[irx1count] = nr;  // ... store in communicator list
        irx1count++; // cannot have more than one rank per block.
      }
    }
    //for (int irx1=0; irx1<nrx1count; irx1++) {
    //  std::cout << "cleaned rank list: rank=" << Globals::my_rank << " irx1=" << irx1 << " rx1ranks=" << rx1ranks[irx1] << std::endl;
    //}
    // Second, we do the slab communicator, to stitch maps together
    int nbslabcount=0;
    for (int nbt=0; nbt<pmy_mesh->nbtotal; nbt++) {
      if (pmy_mesh->loclist[nbt].lx1 == 0) { // only blocks at lowest radial coordinate
        nbslabcount++;
      }
    }
    int* bslabranks = new int[nbslabcount]; // rank indices in bottom slab (duplicates possible)
    int* bslabblock = new int[nbslabcount]; // block indices in bottom slab (no duplicates)
    // again, now we assign the ranks
    int ibslabcount = 0; // use as index
    for (int nbt=0; nbt<pmy_mesh->nbtotal; nbt++) {
      if (pmy_mesh->loclist[nbt].lx1 == 0) {
        bslabranks[ibslabcount] = pmy_mesh->ranklist[nbt];
        bslabblock[ibslabcount] = nbt;
        ibslabcount++;
      }
    }
    //for (int nbt=0; nbt<nbslabcount; nbt++) {
    //  std::cout << "block=" << bslabblock[nbt] << " rank=" << bslabranks[nbt] << " lx1=" << loc.lx1 << std::endl;
    //}
    int nrslabcount = 0;
    int* rcountslab = new int[Globals::nranks];
    for (int nr=0; nr<Globals::nranks; nr++) {
      rcountslab[nr] = 0;
      for (int nbt=0; nbt<nbslabcount; nbt++) {
        if (bslabranks[nbt] == nr) {
          rcountslab[nr] = 1;
        }
      }
      nrslabcount += rcountslab[nr];
    }
    int irslabcount = 0;
    int* rslabranks = new int[nrslabcount];
    for (int nr=0; nr<Globals::nranks; nr++) {
      if (rcountslab[nr] == 1) { // if rank is in global list for this column...
        rslabranks[irslabcount] = nr;  // ... store in communicator list
        irslabcount++; // cannot have more than one rank per block.
      }
    }
    //for (int irslab=0; irslab<nrslabcount; irslab++) {
    //  std::cout << "cleaned rank list: rank=" << Globals::my_rank << " irslab=" << irslab << " rslabranks=" << rslabranks[irslab] << std::endl;
    //}
#ifdef MPI_PARALLEL
    int mpierr;
    std::stringstream msg;
    mpierr = MPI_Comm_group(MPI_COMM_WORLD, &comm_x1.gworld);
    if (mpierr) {
      msg << "[knovae]: MPI_Comm_group error x1 = " << mpierr << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    mpierr = MPI_Group_incl(comm_x1.gworld,nrx1count,rx1ranks,&comm_x1.gsub);
    if (mpierr) {
      msg << "[knovae]: MPI_Group_incl error x1 = " << mpierr << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    mpierr = MPI_Comm_create(MPI_COMM_WORLD,comm_x1.gsub,&comm_x1.comsub);
    if (mpierr) {
      msg << "[knovae]: MPI_Comm_create error x1 = " << mpierr << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    mpierr = MPI_Comm_group(MPI_COMM_WORLD, &comm_slab.gworld);
    if (mpierr) {
      msg << "[knovae]: MPI_Comm_group error slab = " << mpierr << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    mpierr = MPI_Group_incl(comm_slab.gworld,nrslabcount,rslabranks,&comm_slab.gsub);
    if (mpierr) {
      msg << "[knovae]: MPI_Group_incl error slab = " << mpierr << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    mpierr = MPI_Comm_create(MPI_COMM_WORLD,comm_slab.gsub,&comm_slab.comsub);
    if (mpierr) {
      msg << "[knovae]: MPI_Comm_create error slab = " << mpierr << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    // now we have the communicator to sum along the x1 coordinate.
#endif /* MPI_PARALLEL */
    delete [] rx1ranks;
    delete [] rcount;
    delete [] bx1block;
    delete [] bx1ranks;
    delete [] rslabranks;
    delete [] rcountslab;
    delete [] bslabblock;
    delete [] bslabranks;
  }
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Spherical blast wave test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real b0=0.0, angle=0.0;
  Real pi         = 4.0*std::atan(1.0);
  Real r0         = pin->GetReal("problem","r0"); // radius of initial ejecta
  Real dr         = pin->GetReal("problem","dr"); // width of transition
  Real p0         = pin->GetReal("problem","p0"); // ambient pressure
  Real d0         = pin->GetReal("problem","d0"); // ambient density
  Real v0         = pin->GetReal("problem","v0"); // velocity of ejecta
  Real m0         = pin->GetReal("problem","m0"); // mass of ejecta
  Real E0         = pin->GetReal("problem","E0"); // total energy of ejecta
  int  ihomol     = pin->GetOrAddInteger("problem","ihomol",0); // initial homologous expansion (no "ring")
  Real thopen     = pin->GetOrAddReal("problem","thopen",60.0)*PI/180.0; // opening angle of disk ejecta in degrees. 0 means all tidal ejecta.
  Real rwindtidev = pin->GetOrAddReal("problem","rwindtidev",1.0); // ratio between wind and tidal ejecta speed.
  Real rwindtidem = pin->GetOrAddReal("problem","rwindtidem",1.0); // ratio between wind and tidal ejecta mass.
  int  ilog       = pin->GetOrAddInteger("mesh", "ilog", 0);
  Real x1min      = pin->GetReal("mesh","x1min");
  if ((!ihomol) && (x1min >= 0.5*r0)) { // for ring, make sure that x1min < r0/2
    std::stringstream msg;
    msg << "### FATAL ERROR in knovae.cpp ProblemGenerator" << std::endl
        << "x1min > 0.5*r0 " << x1min << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (thopen > 0.5*PI) { 
    std::stringstream msg;
    msg << "### FATAL ERROR in knovae.cpp ProblemGenerator" << std::endl
        << "thopen > 0.5*pi " << thopen << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    b0    = pin->GetReal("problem","b0");
    angle = (PI/180.0)*pin->GetReal("problem","angle");
  }
  Real gamma   = peos->GetGamma();
  Real gm1     = gamma - 1.0;
  Real T       = gm1*E0/m0; 
  Real mtide   = m0/(1.0+rwindtidem);
  Real mwind   = rwindtidem*mtide;
  Real voltide = 4.0*PI*pow(r0,3)*std::cos(thopen)/3.0;
  Real volwind = 4.0*PI*pow(r0,3)*(1.0-std::cos(thopen))/3.0;
  Real rhotide = mtide/voltide;
  Real rhowind = mwind/volwind;
  Real ethtide = rhotide*T/gm1;
  Real ethwind = rhowind*T/gm1;
  Real vtide   = v0;
  Real vwind   = v0*rwindtidev;
  Real rho0    = d0;
  Real e0      = p0/gm1;

  // get coordinates of center of blast, and convert to Cartesian if necessary
  Real x1_0   = pin->GetOrAddReal("problem","x1_0",0.0);
  Real x2_0   = pin->GetOrAddReal("problem","x2_0",0.0);
  Real x3_0   = pin->GetOrAddReal("problem","x3_0",0.0);
  Real x0,y0,z0;
  if (COORDINATE_SYSTEM == "cartesian") {
    x0 = x1_0;
    y0 = x2_0;
    z0 = x3_0;
  } else if (COORDINATE_SYSTEM == "cylindrical") {
    x0 = x1_0*std::cos(x2_0);
    y0 = x1_0*std::sin(x2_0);
    z0 = x3_0;
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    x0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
    y0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
    z0 = x1_0*std::cos(x2_0);
  } else {
    // Only check legality of COORDINATE_SYSTEM once in this function
    std::stringstream msg;
    msg << "### FATAL ERROR in knovae.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM= " << COORDINATE_SYSTEM << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // setup uniform ambient medium with spherical over-pressured region
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real rad,vwindrad,vtiderad,radv,thet;
        if (COORDINATE_SYSTEM == "cartesian") {
          Real x   = pcoord->x1v(i);
          Real y   = pcoord->x2v(j);
          Real z   = pcoord->x3v(k);
          rad      = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
          thet   = std::acos(z/rad);
        } else if (COORDINATE_SYSTEM == "cylindrical") {
          Real x = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
          Real z = pcoord->x3v(k);
          rad    = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
          thet   = std::acos(z/rad);
        } else { // if (COORDINATE_SYSTEM == "spherical_polar")
          Real x = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
          Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
          Real z = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          rad    = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
          thet   = pcoord->x2v(j);
        }

        // factors for radial and polar dependence
        Real fm1rad            = 0.5*(1.0-std::tanh((rad-r0)/dr)); // radial drop off
        Real fm1the            = 0.25*((1.0+std::tanh((thet-0.5*(pi-thopen))/0.01))*(1.0-std::tanh((thet-0.5*(pi+thopen))/0.01))); // polar drop off
        if (ilog) {
          radv = (rad-x1min)/(r0-x1min);
        } else {
          radv = rad;
        }
        if (ihomol) { // ramp to make ring
          vwindrad = vwind*rad/r0;
          vtiderad = vtide*rad/r0;
        } else {
          if (rad <= 0.5*r0) {
            vwindrad      = 2.0*vwind*rad/r0;
            vtiderad      = 2.0*vtide*rad/r0;
          } else {
            vwindrad      = vwind;
            vtiderad      = vtide;
          }
        }
        phydro->u(IDN,k,j,i) = rho0+(rhowind-rho0)*fm1rad+(rhotide-rhowind)*fm1the*fm1rad; 
        phydro->u(IM1,k,j,i) = (vwindrad+(vtiderad-vwindrad)*fm1the)*fm1rad*phydro->u(IDN,k,j,i);
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = e0+(ethwind-e0)*fm1rad+(ethtide-ethwind)*fm1the*fm1rad
                                +0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          if (RELATIVISTIC_DYNAMICS)  // this should only ever be SR with this file
            phydro->u(IEN,k,j,i) += phydro->u(IDN,k,j,i);
        }
        if (DUAL_ENERGY) {
          phydro->u(IIE,k,j,i) = e0+(ethwind-e0)*fm1rad+(ethtide-ethwind)*fm1the*fm1rad;
        }
        if (NSCALARS == 2) {
          // the first index corresponds to wind ejecta, the second to tidal ejecta
          phydro->u(NHYDRO-NSCALARS  ,k,j,i) = rhowind*fm1rad*(1.0-fm1the)+1e-4;       // wind ejecta
          phydro->u(NHYDRO-NSCALARS+1,k,j,i) = rhotide*fm1rad*fm1the+1e-4; // tidal ejecta
        }
      }
    }
  }

  // initialize interface B and total energy
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie+1; ++i) {
          if (COORDINATE_SYSTEM == "cartesian") {
            pfield->b.x1f(k,j,i) = b0 * std::cos(angle);
          } else if (COORDINATE_SYSTEM == "cylindrical") {
            Real phi = pcoord->x2v(j);
            pfield->b.x1f(k,j,i) =
                b0 * (std::cos(angle) * std::cos(phi) + std::sin(angle) * std::sin(phi));
          } else { //if (COORDINATE_SYSTEM == "spherical_polar") {
            Real theta = pcoord->x2v(j);
            Real phi = pcoord->x3v(k);
            pfield->b.x1f(k,j,i) = b0 * std::abs(std::sin(theta))
                * (std::cos(angle) * std::cos(phi) + std::sin(angle) * std::sin(phi));
          }
        }
      }
    }
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je+1; ++j) {
        for (int i = is; i <= ie; ++i) {
          if (COORDINATE_SYSTEM == "cartesian") {
            pfield->b.x2f(k,j,i) = b0 * std::sin(angle);
          } else if (COORDINATE_SYSTEM == "cylindrical") {
            Real phi = pcoord->x2v(j);
            pfield->b.x2f(k,j,i) =
                b0 * (std::sin(angle) * std::cos(phi) - std::cos(angle) * std::sin(phi));
          } else { //if (COORDINATE_SYSTEM == "spherical_polar") {
            Real theta = pcoord->x2v(j);
            Real phi = pcoord->x3v(k);
            pfield->b.x2f(k,j,i) = b0 * std::cos(theta)
                * (std::cos(angle) * std::cos(phi) + std::sin(angle) * std::sin(phi));
            if (std::sin(theta) < 0.0)
              pfield->b.x2f(k,j,i) *= -1.0;
          }
        }
      }
    }
    for (int k = ks; k <= ke+1; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          if (COORDINATE_SYSTEM == "cartesian" || COORDINATE_SYSTEM == "cylindrical") {
            pfield->b.x3f(k,j,i) = 0.0;
          } else { //if (COORDINATE_SYSTEM == "spherical_polar") {
            Real phi = pcoord->x3v(k);
            pfield->b.x3f(k,j,i) =
                b0 * (std::sin(angle) * std::cos(phi) - std::cos(angle) * std::sin(phi));
          }
        }
      }
    }
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          phydro->u(IEN,k,j,i) += 0.5*b0*b0;
        }
      }
    }
  }

#ifdef WURSTBROT
  //std::cout << "p=" << std::setw(3) << Globals::my_rank << " lid=" <<  std::setw(3) << MeshBlock::lid << " gid=" << std::setw(3) << MeshBlock::gid;
  //std::cout << " lx1=" << std::setw(3) << MeshBlock::loc.lx1 << " lx2=" << std::setw(3) << MeshBlock::loc.lx2 << " lx3=" << std::setw(3) << MeshBlock::loc.lx3 << " level=" << std::setw(3) << MeshBlock::loc.level;
  //std::cout << " nneighbor=" << nngb << " RS=" << std::setprecision(5) << MeshBlock::block_size.x1min << "," << MeshBlock::block_size.x1max << "," << MeshBlock::block_size.x2min << "," << MeshBlock::block_size.x2max << std::endl;
  
  //std::cout << "p=" << std::setw(3) << Globals::my_rank;
  //for (int n=0; n<nngb; n++) {
  //  std::cout << " |n=" << std::setw(3) << n << " rank=" << std::setw(3) << (MeshBlock::pbval->neighbor[n]).rank << " lid=" << std::setw(3) << (MeshBlock::pbval->neighbor[n]).lid << " gid=" << std::setw(3) << (MeshBlock::pbval->neighbor[n]).gid;
   // std::cout << " ox1,2,3=" << std::setw(2) << (MeshBlock::pbval->neighbor[n]).ox1 << "," << std::setw(2) << (MeshBlock::pbval->neighbor[n]).ox2 << "," << std::setw(2) << (MeshBlock::pbval->neighbor[n]).ox3;
  //}
  //std::cout << std::endl;

  // Ok. We need to set up the communication infrastructure. 
  // (1) We need a communicator (like we had in Athena) for summation along the radial direction.
  // We assume spherical coordinates for now. Grouping into communicators means we need to 
  // figure out all processors with a certain logical location. In our case, we sum along r, so
  // we need lx1=0...nrbx1-1 for all lx2. Also note that several MeshBlocks can live on the same
  // (processor) rank. 
  // (2) Sum over all MeshBlocks (one or more per rank) at level (lid) = 0. Since restriction is
  //     conservative, this should be good enough for mixing. 
  // (3) MPI_Allreduce over correct communicator.
  // For AMR, we would have to redo this occasionally ....

  //if (Globals::my_rank == 0) {
  //  for (int nbt=0; nbt<pmy_mesh->nbtotal; nbt++) {
  //    std::cout << "block=" << nbt << " rank=" << pmy_mesh->ranklist[nbt];
  //    std::cout <<  " level=" << pmy_mesh->loclist[nbt].level;
  //    std::cout << " loc=" << pmy_mesh->loclist[nbt].lx1 << " " << pmy_mesh->loclist[nbt].lx2 << " "<< pmy_mesh->loclist[nbt].lx3;
  //    std::cout << std::endl;
  //  }
  //}

  //std::cout << "rank=" << Globals::my_rank << " loc.level=" << loc.level << " lid=" << lid << std::endl;
  if (lid == 0) { // only defined on root level
    // count ranks in local loc.lx2,loc.lx3 column. Complicated way to express nrbx1.
    int nbx1count=0;
    for (int nbt=0; nbt<pmy_mesh->nbtotal; nbt++) {
      if (  (loc.lx2   == pmy_mesh->loclist[nbt].lx2)
          &&(loc.lx3   == pmy_mesh->loclist[nbt].lx3)
          &&(loc.level == pmy_mesh->loclist[nbt].level)) {
        nbx1count++;
      }
    }
    int* bx1ranks = new int[nbx1count]; // rank indices along one column (duplicates possible)
    int* bx1block = new int[nbx1count]; // block indices along one column (no duplicates)
    // again, now we assign the ranks
    int ibx1count = 0; // use as index
    for (int nbt=0; nbt<pmy_mesh->nbtotal; nbt++) {
      if (  (loc.lx2   == pmy_mesh->loclist[nbt].lx2)
          &&(loc.lx3   == pmy_mesh->loclist[nbt].lx3)
          &&(loc.level == pmy_mesh->loclist[nbt].level)) {
        bx1ranks[ibx1count] = pmy_mesh->ranklist[nbt];
        bx1block[ibx1count] = nbt;
        ibx1count++;
      }
    }
    //for (int nbt=0; nbt<nbx1count; nbt++) {
    //  std::cout << "block=" << bx1block[nbt] << " rank=" << bx1ranks[nbt] << " lx2=" << loc.lx2 << " lx3=" << loc.lx3 << std::endl;
    //}  
    // At this point, each rank knows which blocks to sum over.
    // Now we need to clean the rank list (removing duplicates).
    int nrx1count = 0;
    int* rcount = new int[Globals::nranks];
    for (int nr=0; nr<Globals::nranks; nr++) {
      rcount[nr] = 0;
      for (int nbt=0; nbt<nbx1count; nbt++) { 
        if (bx1ranks[nbt] == nr) {
          rcount[nr] = 1;
        }
      }
      nrx1count += rcount[nr];
    }
    int irx1count = 0;
    int* rx1ranks = new int[nrx1count];
    for (int nr=0; nr<Globals::nranks; nr++) {
      if (rcount[nr] == 1) { // if rank is in global list for this column...
        rx1ranks[irx1count] = nr;  // ... store in communicator list
        irx1count++; // cannot have more than one rank per block.
      }
    }
    for (int irx1=0; irx1<nrx1count; irx1++) {
      std::cout << "cleaned rank list: rank=" << Globals::my_rank << " irx1=" << irx1 << " rx1ranks=" << rx1ranks[irx1] << std::endl;
    }
#ifdef MPI_PARALLEL
    int mpierr;
    mpierr = MPI_Comm_group(MPI_COMM_WORLD, &comm_x1.gworld);
    mpierr = MPI_Group_incl(comm_x1.gworld,nrx1count,rx1ranks,&comm_x1.gsub);
    mpierr = MPI_Comm_create(MPI_COMM_WORLD,comm_x1.gsub,&comm_x1.comsub);
    // now we have the communicator to sum along the x1 coordinate.
#endif /* MPI_PARALLEL */
    delete [] rx1ranks;
    delete [] rcount;
    delete [] bx1block;
    delete [] bx1ranks;
  }
#endif // WURSTBROT

}

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief otf diagnostics (shell position, sphericity etc)
//========================================================================================

void MeshBlock::UserWorkInLoop(void) {
  return;
}

//========================================================================================
//! \fn void MeshBlock::OTFWorkBeforeOutput(void)
//  \brief otf diagnostics (shell position, sphericity etc)
//========================================================================================

void MeshBlock::OTFWorkBeforeOutput(ParameterInput *pin) {
#ifdef MPI_PARALLEL
  int mpierr;
#endif
  if (NSCALARS == 2) {
    int nscl = NSCALARS;
    int nx3  = ke-ks+1; // These are local!
    int nx2  = je-js+1;
    int ioff = 4*(nscl+1);
    Real rad, thet, x, y, z, dv, d, p, dc, c0, c1, v1, mixpar, mixmas;
    Real *avgparam = new Real[4*(nscl+1)+2*(nscl+1)+3]; // average rad, the and their mass weights, for c0, c1 and c0+c1
    Real *stdparam = new Real[2*(nscl+1)]; // standard deviation of rad, the for c0, c1 and c0+c1
    //Real *hydparam = new Real[2*(nscl+1)]; // density and temperature for c0, c1, and c0+c1
    Real *mixparam = new Real[2*nx2*nx3];  // mixing parameter and mass weight, averaged over radius
    Real *velparam = new Real[6*nx2*nx3];  // radial velocities and mass weights for c0, c1, c0+c1
    //Real *volparam = new Real[3];          // volume
    for (int s=0; s<2*(nscl+1); s++) {
      stdparam[s] = 0.0;
    }
    for (int s=0; s<6*(nscl+1)+3; s++) 
      avgparam[s] = 0.0;
    for (int s=0; s<2*nx2*nx3; s++)
      mixparam[s] = 0.0;
    for (int s=0; s<6*nx2*nx3; s++)
      velparam[s] = 0.0;

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          dv = pcoord->GetCellVolume(k,j,i);
          if (COORDINATE_SYSTEM == "cartesian") {
            x    = pcoord->x1v(i);
            y    = pcoord->x2v(j);
            z    = pcoord->x3v(k);
            rad  = std::sqrt(SQR(x) + SQR(y) + SQR(z));
            thet = std::acos(z/rad);
          } else if (COORDINATE_SYSTEM == "cylindrical") {
            x    = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
            y    = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
            z    = pcoord->x3v(k);
            rad  = std::sqrt(SQR(x) + SQR(y) + SQR(z));
            thet = std::acos(z/rad);
          } else { // if (COORDINATE_SYSTEM == "spherical_polar")
            x    = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
            y    = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
            z    = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
            rad  = std::sqrt(SQR(x) + SQR(y) + SQR(z));
            thet = pcoord->x2v(j);
          }
          d                                      = phydro->u(IDN,k,j,i);
          v1                                     = phydro->u(IM1,k,j,i)/d;
          c0                                     = phydro->u(NHYDRO-nscl  ,k,j,i);
          c1                                     = phydro->u(NHYDRO-nscl+1,k,j,i);
          dc                                     = c0+c1;
          p                                      = phydro->w(IPR,k,j,i)/d;
          // fastest index: rad, thet. middle index: scal. slowest index: quantity, mass
          avgparam[(nscl+1)*2*0 +2*0 + 0]       += rad *c0*dv; //positions
          avgparam[(nscl+1)*2*0 +2*0 + 1]       += thet*c0*dv;
          avgparam[(nscl+1)*2*0 +2*1 + 0]       += rad *c1*dv;
          avgparam[(nscl+1)*2*0 +2*1 + 1]       += thet*c1*dv;
          avgparam[(nscl+1)*2*0 +2*2 + 0]       += rad *dc*dv;
          avgparam[(nscl+1)*2*0 +2*2 + 1]       += thet*dc*dv;
          avgparam[(nscl+1)*2*1 +2*0 + 0]       += c0*dv; // masses
          avgparam[(nscl+1)*2*1 +2*0 + 1]       += c0*dv;
          avgparam[(nscl+1)*2*1 +2*1 + 0]       += c1*dv;
          avgparam[(nscl+1)*2*1 +2*1 + 1]       += c1*dv;
          avgparam[(nscl+1)*2*1 +2*2 + 0]       += dc*dv;
          avgparam[(nscl+1)*2*1 +2*2 + 1]       += dc*dv;
          avgparam[ioff+2*0+0]                  += p*c0*dv; // temperature * density...
          avgparam[ioff+2*0+1]                  += c0*dv; // density
          avgparam[ioff+2*1+0]                  += p*c1*dv;
          avgparam[ioff+2*1+1]                  += c1*dv; 
          avgparam[ioff+2*2+0]                  += p*dc*dv; 
          avgparam[ioff+2*2+1]                  += dc*dv;
          avgparam[ioff+2*(nscl+1)+0]           += c0/d*dv; // volume
          avgparam[ioff+2*(nscl+1)+1]           += c1/d*dv;
          avgparam[ioff+2*(nscl+1)+2]           += dc/d*dv;
          mixpar                                 = 0.5*(1.0+(c1-c0)/(c1+c0));
          mixmas                                 = dc*dv;
          mixparam[nx3*nx2*0+nx2*(k-ks)+(j-js)] += mixpar*mixmas;
          mixparam[nx3*nx2*1+nx2*(k-ks)+(j-js)] += mixmas;
          velparam[nx3*nx2*0+nx2*(k-ks)+(j-js)] += v1*dc*dv;
          velparam[nx3*nx2*1+nx2*(k-ks)+(j-js)] += dc*dv;
          velparam[nx3*nx2*2+nx2*(k-ks)+(j-js)] += v1*c0*dv;
          velparam[nx3*nx2*3+nx2*(k-ks)+(j-js)] += c0*dv;
          velparam[nx3*nx2*4+nx2*(k-ks)+(j-js)] += v1*c1*dv;
          velparam[nx3*nx2*5+nx2*(k-ks)+(j-js)] += c1*dv;
        }
      }
    }
#ifdef MPI_PARALLEL
    std::stringstream msg;
    mpierr = MPI_Allreduce(MPI_IN_PLACE,avgparam,6*(nscl+1)+3,MPI_ATHENA_REAL,MPI_SUM,MPI_COMM_WORLD); 
    if (mpierr) {
      msg << "[knovae]: MPI_Allreduce error avgparam = " << mpierr << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    if (comm_x1.comsub != MPI_COMM_NULL) {
      mpierr = MPI_Allreduce(MPI_IN_PLACE,mixparam,2*nx2*nx3 ,MPI_ATHENA_REAL,MPI_SUM,comm_x1.comsub);
      if (mpierr) {
        msg << "[knovae]: MPI_Allreduce error mixparam = " << mpierr << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      mpierr = MPI_Allreduce(MPI_IN_PLACE,velparam,6*nx2*nx3 ,MPI_ATHENA_REAL,MPI_SUM,comm_x1.comsub);
      if (mpierr) {
        msg << "[knovae]: MPI_Allreduce error velparam = " << mpierr << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }
#endif // MPI_PARALLEL
    for (int s=0; s<2*(nscl+1); s++) 
      avgparam[s] /= avgparam[s+2*(nscl+1)]; // these are rad and thet
    for (int s=0; s<nscl+1; s++) {
      avgparam[ioff+2*s+0] /= avgparam[ioff+2*s+1]; // temperature
    }
    for (int s=0; s<nscl+1; s++) {
      avgparam[ioff+2*s+1] /= avgparam[ioff+2*(nscl+1)+s]; // density
    }
    for (int s=0; s<nx2*nx3; s++) {
      mixparam[s]           /= mixparam[s+nx2*nx3];
      velparam[s+0*nx2*nx3] /= velparam[s+1*nx2*nx3];
      velparam[s+2*nx2*nx3] /= velparam[s+3*nx2*nx3];
      velparam[s+4*nx2*nx3] /= velparam[s+5*nx2*nx3];
    }
      
    // still need to get dispersion of color-averaged theta
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          dv = pcoord->GetCellVolume(k,j,i);
          if (COORDINATE_SYSTEM == "cartesian") {
            x    = pcoord->x1v(i);
            y    = pcoord->x2v(j);
            z    = pcoord->x3v(k);
            rad  = std::sqrt(SQR(x) + SQR(y) + SQR(z));
            thet = std::acos(z/rad);
          } else if (COORDINATE_SYSTEM == "cylindrical") {
            x    = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
            y    = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
            z    = pcoord->x3v(k);
            rad  = std::sqrt(SQR(x) + SQR(y) + SQR(z));
            thet = std::acos(z/rad);
          } else { // if (COORDINATE_SYSTEM == "spherical_polar")
            x    = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
            y    = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
            z    = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
            rad  = std::sqrt(SQR(x) + SQR(y) + SQR(z));
            thet = pcoord->x2v(j);
          }
          c0                                     = phydro->u(NHYDRO-nscl  ,k,j,i);
          c1                                     = phydro->u(NHYDRO-nscl+1,k,j,i);
          dc                                     = c0+c1;
          // fastest index: rad, thet. middle index: scal. slowest index: quantity, mass
          stdparam[2*0 + 0] += SQR(rad -avgparam[(nscl+1)*2*0 + 2*0 + 0]) *c0*dv;
          stdparam[2*0 + 1] += SQR(thet-avgparam[(nscl+1)*2*0 + 2*0 + 1]) *c0*dv;
          stdparam[2*1 + 0] += SQR(rad -avgparam[(nscl+1)*2*0 + 2*1 + 0]) *c1*dv;
          stdparam[2*1 + 1] += SQR(thet-avgparam[(nscl+1)*2*0 + 2*1 + 1]) *c1*dv;
          stdparam[2*2 + 0] += SQR(rad -avgparam[(nscl+1)*2*0 + 2*2 + 0]) *dc*dv;
          stdparam[2*2 + 1] += SQR(thet-avgparam[(nscl+1)*2*0 + 2*2 + 1]) *dc*dv;
        }
      }
    }
#ifdef MPI_PARALLEL
    mpierr = MPI_Allreduce(MPI_IN_PLACE,stdparam,2*(nscl+1),MPI_ATHENA_REAL,MPI_SUM,MPI_COMM_WORLD); 
    if (mpierr) {
      msg << "[knovae]: MPI_Allreduce error stdparam = " << mpierr << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
#endif // MPI_PARALLEL
    for (int s=0; s<=2*(nscl+1); s++) 
      stdparam[s] = sqrt(stdparam[s]/avgparam[s+2*(nscl+1)]); 

    // at this point, we have the averages over all ranks, and the mixing parameter for each column separately.
    // Now stitch together the columns.
    // The averages can be directly copied into otf_data.data.
    // The mixing parameters mixparam are actually 2D maps (theta,phi). All ranks along a column in r have
    // a copy of the local map. These need to be stitched together via Allgather. Only one set (at constant r)
    // of ranks has to communicate the results, though. We need a special communicator for that, too.
    // We have now the slab communicator for the bottom slab. 
    int off = 0;
    for (int s=0; s<6*(nscl+1)+3; s++) 
      otf_data.data[    s] = avgparam[s];
    off = 6*(nscl+1)+3;
    for (int s=0; s<2*(nscl+1); s++)
      otf_data.data[off+s] = stdparam[s];
    off = 8*(nscl+1)+3;
#ifdef MPI_PARALLEL
    if (comm_slab.comsub != MPI_COMM_NULL) {
      mpierr = MPI_Allgather(mixparam, 2*nx2*nx3, MPI_ATHENA_REAL, &(otf_data.data[off]), 2*nx2*nx3, MPI_ATHENA_REAL, comm_slab.comsub);
      if (mpierr) {
        msg << "[knovae]: MPI_Allgather error mixparam = " << mpierr << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      off = 8*(nscl+1)+3+2*pmy_mesh->mesh_size.nx2*pmy_mesh->mesh_size.nx3;
      mpierr = MPI_Allgather(velparam, 6*nx2*nx3, MPI_ATHENA_REAL, &(otf_data.data[off]), 6*nx2*nx3, MPI_ATHENA_REAL, comm_slab.comsub);
      if (mpierr) {
        msg << "[knovae]: MPI_Allgather error mixparam = " << mpierr << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }
#else
    for (int s=0; s<8*nx2*nx3; s++)
      otf_data.data[off+s] = mixparam[s];
#endif
    delete [] avgparam;
    delete [] stdparam;
    delete [] mixparam;
  } // if (NSCALARS == 2)

  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief Any clean-up etc.
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  MeshBlock *pmb = pblock;
  delete [] pmb->otf_data.data;
  return;
}
