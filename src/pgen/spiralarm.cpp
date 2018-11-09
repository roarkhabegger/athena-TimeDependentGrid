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

//====================================================================================
// local functions
Real find_dprof(AthenaArray<Real> &dprof,const int nx2, const Real x2min, const Real x2max,
                Real csound, Real surfn0, Real surfs0, Real sigs0, Real beta0);
void rk4drive (Real *y0, const Real x1, const Real x2, const int nstep, const Real b, const Real c, AthenaArray<Real> &y);
void derivs(Real x, Real *y, Real b, Real c, Real *dydx);
static Real gravpot_spirarm1(const Real x1, const Real x2, const Real x3, const Real time);
void SpiralInflowInnerX1(MeshBlock *pmb, Coordinates *pco,
                         AthenaArray<Real> &a,
                         FaceField &b, Real time, Real dt,
                         int is, int ie, int js, int je, int ks, int ke, int ngh);
void SpiralOutflowOuterX1(MeshBlock *pmb, Coordinates *pco,
                         AthenaArray<Real> &a,
                         FaceField &b, Real time, Real dt,
                         int is, int ie, int js, int je, int ks, int ke, int ngh);
Real GetMeanDensity(Coordinates *pcoord, Hydro *phydro,
                    const int is, const int ie, const int js, const int je, const int ks, const int ke);

static void stop_this();

//====================================================================================
// Enroll user-specific functions
void Mesh::InitUserMeshData(ParameterInput *pin) {
  // vertical external gravitational potential
  EnrollStaticGravPotFunction(gravpot_spirarm1);

  // enroll user-defined boundary conditions
  if (mesh_bcs[INNER_X1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(INNER_X1, SpiralInflowInnerX1);
  }
  if (mesh_bcs[OUTER_X1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(OUTER_X1, SpiralOutflowOuterX1);
  }

  Real four_pi_G = 48.0*std::atan(1.0);
  Real eps = 0.0;
  SetFourPiG(four_pi_G);
  SetGravityThreshold(eps);
  SetMeanDensity(0.0); // temporary -- will be set to mean density in initialization

  return;
}

//====================================================================================
// short for debugging interrupt
static void stop_this() {
  std::stringstream msg;
  msg << "stop" << std::endl;
  throw std::runtime_error(msg.str().c_str());
}

//====================================================================================
// Real GetMeanDensity(void)
Real GetMeanDensity(Coordinates *pcoord, Hydro *phydro,
                    const int is, const int ie, const int js, const int je, const int ks, const int ke) 
{
  std::stringstream msg;
#ifdef MPI_PARALLEL
  int  mpierr;
#endif
  Real mass[2], gmass[2];
  mass[0] = 0.0; mass[1] = 0.0;
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real dv  = pcoord->GetCellVolume(k,j,i);
        mass[0] += phydro->u(IDN,k,j,i)*dv;
        mass[1] += dv;
      }
    }
  }
#ifdef MPI_PARALLEL
  mpierr = MPI_Allreduce(&mass, &gmass, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (mpierr) {
    msg << "[GetMeanDensity]: MPI_Allreduce error = "
        << mpierr << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  for (int n=0; n<2; n++) mass[n] = gmass[n];
#endif // MPI_PARALLEL
  mass[0] /= mass[1];
  return mass[0];
}

//====================================================================================
// Global variables for boundaries and gravity
Real vx0_g, x1len_g, cphi0_g, phisp_g, ramp_g;
AthenaArray<Real> dprofloc;
Real vx0_bc,vy0_bc,rgc0_bc,omega0_bc,T0_bc,gm1_bc,bx0_bc,by0_bc,bz0_bc,csound_bc,n00_bc,beta0_bc,bpow_bc;
int imode_bc;

//====================================================================================
// Static gravitational potential for spiral arms
static Real gravpot_spirarm1(const Real x1, const Real x2, const Real x3, const Real time)
{
  Real gfrac  = std::min(ramp_g*time*vx0_g/x1len_g,1.0);
  Real grvpot = cphi0_g*x2*x2 - gfrac*phisp_g*cos(6.283185307179586*x1/x1len_g);
  return grvpot;
}

//====================================================================================
//
/*------------------------------------------------------------------*/
/* Function find_dprof: finds initial vertical density profile      */
/* see Kim & Ostriker 2002, ApJ 581, 1080                           */
/* Given a gas column density surfn0, find the corresponding        */
/* density profile via bracketing and integration.                  */

Real find_dprof(AthenaArray<Real> &dprof,const int nx2, const Real x2min, const Real x2max,
                Real csound, Real surfn0, Real surfs0, Real sigs0, Real beta0)
{
  int j, ic;
  Real pi = 4.0*std::atan(1.0);
  Real x2len, Hinf, H0, s0, n00, hhout, dx2, a, b, c, c1, c2, ha, hb, hh0, hh1, hhm;
  Real eps = 1e-7;
  Real y0[2];
  AthenaArray<Real> dprof2;

  x2len  = x2max-x2min;
  dx2    = 0.5*x2len/((Real) nx2);
  dprof2.NewAthenaArray(2,nx2);
/* prepare the coefficients */
  Hinf   = SQR(csound)/(pi*surfn0);
  s0     = SQR(sigs0*surfn0/(csound*surfs0));
  hb     = 1.0+0.5*beta0;
  ha     = 0.5*sqrt(pi*s0*(1.0+0.5*beta0));
  hh0    = std::min(ha,hb);
  hh1    = std::max(ha,hb);
  a      = 0.5*(1.0+0.5*beta0);
  c      = -1.0/s0;
  c2     = c/a;
  ic     = 0;
  y0[0] = 1.0;
  y0[1] = 0.0;
  //if (Globals::my_rank == 0) {
  //  std::cout << "[find_dprof]: nx2=" << nx2 << " x2len=" << x2len << " Hinf=" << Hinf << " c1=" << c1 << " c2=" << c2 << std::endl;
  //  std::cout << "[find_dprof]: dx2=" << dx2 << " csound=" << csound << " surfn0=" << surfn0 << " sigs0=" << sigs0 << " beta0=" << beta0 <<  std::endl;
  //}
  do {
    hhm = 0.5*(hh0+hh1);
    c1  = -1.0/(hhm*a);
    rk4drive(y0,0.0,0.5*x2len/Hinf,nx2,c1,c2,dprof2);

    //stop_this();

    /* test for surfn0 */
    hhout = 0.0;
    for (j=0;j<nx2;j++) {
      //if (Globals::my_rank == 0) std::cout << "j=" << j << " dprof2(0,j)=" << dprof2(0,j) << std::endl;
      hhout += dprof2(0,j);
    }
    hhout *= 2.0*dx2/Hinf;
    n00 = 0.5*surfn0/(hhout*Hinf);
    //if (Globals::my_rank == 0) {
    //  std::cout << "[find_dprof]: ic=" << std::setw(5) << ic 
    //            << " n00="   << std::scientific << std::setprecision(5) << n00
    //            << " hhout=" << std::scientific << std::setprecision(5) << hhout
    //            << " surf="  << std::scientific << std::setprecision(5) << 2.0*hhout*Hinf*n00
    //            << " hh0="   << std::scientific << std::setprecision(5) << hh0
    //            << " hh1="   << std::scientific << std::setprecision(5) << hh1
    //            << " Hinf="  << std::scientific << std::setprecision(5) << Hinf
    //            << std::endl;
    //}
    if (hhout > hhm) hh0 = hhm;
    if (hhout < hhm) hh1 = hhm;
    ic++;
  } while (fabs(hh0-hh1) > eps);
  for (j=0; j<nx2; j++) dprof(j) = dprof2(0,j)*n00;
  dprof2.DeleteAthenaArray();
  return n00;
}

/*------------------------------------------------------------------*/
/* Function rk4drive                                                */
/* Note that this assumes quantities at cell faces. Need to inter-  */
/* polate at the end.                                               */
/* Need two ODE here. First is actual result (dy/dx=z)              */
/* nstep only nx2/2                                                 */

void rk4drive (Real* y0, const Real x1, const Real x2, const int nstep, const Real b, const Real c, AthenaArray<Real> &y) {
  int i,k,nvar=2;
  Real x,h;
  Real yt[2], dydx[2], k1[2], k2[2], k3[2], k4[2];

  for (i=0;i<nvar;i++) {
    yt[i]    = y0[i];
    y(i,0)   = y0[i];
  }
  x    = x1;
  h    = (x2-x1)/nstep;
  for (k=0; k<nstep; k++) {
    derivs(x,yt,b,c,dydx); 
    for (i=0;i<nvar;i++) {
      k1[i] = h*dydx[i];
      yt[i] = y(i,k)+0.5*k1[i];
    }
    derivs(x+0.5*h,yt,b,c,dydx);
    for (i=0;i<nvar;i++) {
      k2[i] = h*dydx[i];
      yt[i] = y(i,k)+0.5*k2[i];
    }
    derivs(x+0.5*h,yt,b,c,dydx);
    for (i=0;i<nvar;i++) {
      k3[i] = h*dydx[i];
      yt[i] = y(i,k)+k3[i];
    }
    derivs(x+h,yt,b,c,dydx);
    for (i=0;i<nvar;i++) {
      k4[i]     = h*dydx[i];
      y(i,k+1)  = y(i,k)+(k1[i]+0.5*(k2[i]+k3[i])+k4[i])/6.0;
    }
    //if (Globals::my_rank == 0) {
    //  std::cout << "[rk4drive]: k=" << std::setw(5) << k 
    //            << " x="   << std::scientific << std::setprecision(5) << x 
    //            << " y(k)=" << std::scientific << std::setprecision(5) << y(0,k)
    //            << " y(k+1)="  << std::scientific << std::setprecision(5) << y(0,k+1)
    //            << " k1="   << std::scientific << std::setprecision(5) << k1[0]
    //            << " k2="   << std::scientific << std::setprecision(5) << k2[0]
    //            << " k3="  << std::scientific << std::setprecision(5) << k3[0]
    //            << " k4="  << std::scientific << std::setprecision(5) << k4[0]
    //            << std::endl;
    //}
    x = x+h;
  }
  return;
}

/*------------------------------------------------------------------*/
/* Function derivs                                                  */

void derivs(Real x, Real *y, Real b, Real c, Real *dydx) {
  dydx[0] = y[1];
  dydx[1] = y[1]*y[1]/y[0] + b*y[0]*y[0] + c*y[0];
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Spherical blast wave test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real T0, csound,d0,vx0,vy0,vz0,bx0,by0,bz0,bpow,qtoom0,beta0,kap0,surfn0;
  AthenaArray<Real> dprof, dprofcmp;
  Real pi       = 4.0*std::atan(1.0);
  int  nx1      = pin->GetInteger("mesh","nx1");
  int  nx2      = pin->GetInteger("mesh","nx2");
  int  nx3      = pin->GetInteger("mesh","nx3");
  int  nz, nzloc, iz, izloc, izs, IZ, IY;
  Real x1min    = pin->GetReal("mesh","x1min");
  Real x1max    = pin->GetReal("mesh","x1max");
  Real x2min    = pin->GetReal("mesh","x2min");
  Real x2max    = pin->GetReal("mesh","x2max");
  Real x3min    = pin->GetReal("mesh","x3min");
  Real x3max    = pin->GetReal("mesh","x3max");
  Real x1len    = x1max-x1min;
  Real x2len    = x2max-x2min;
  Real x3len    = x3max-x3min;
  Real zmin,zmax,zloc0;
  if (nx3 > 1) {
    izs   = ks;
    nz    = nx3;
    nzloc = ke-ks+1;
    zmin  = x3min; 
    zmax  = x3max;
    zloc0 = pcoord->x3v(izs);
    IY    = IM2;
    IZ    = IM3;
  } else {
    izs   = js;
    nz    = nx2; 
    nzloc = je-js+1;
    zmin  = x2min;
    zmax  = x2max;
    zloc0 = pcoord->x2v(izs);
    IY    = IM3;
    IZ    = IM2;
  }
  Real zlen     = zmax-zmin;
  int  imode    = pin->GetInteger("problem","imode");  // 0 physical variables, 1: scaled variables
  Real sigs0    = pin->GetReal("problem","sigs0");     // stellar velocity dispersion in z
  Real surfs0   = pin->GetReal("problem","surfs0");    // stellar surface density
  Real gamma    = peos->GetGamma();
  Real gm1      = gamma-1.0;
  if (NON_BAROTROPIC_EOS) {
    T0          = pin->GetReal("problem","T0");        // temperature
    csound      = sqrt(gamma*T0);
  } else {
    csound      = peos->GetIsoSoundSpeed();
  }
  Real rgc0     = pin->GetReal("problem","rgc0");      // galacto-centric radius 
  Real omega0   = pin->GetReal("problem","omega0");    // orbital period at rgc0
  Real sini0    = pin->GetReal("problem","sini0");     // sin i
  Real m0       = pin->GetReal("problem","m0");        // number of spiral arms
  Real frac0    = pin->GetReal("problem","frac0");     // strength of potential 
  Real ramp     = pin->GetOrAddReal("problem","ramp",0.1); // ramp for potential
  Real sinim0   = sini0/m0;
  Real cphi0    = SQR((PI*surfs0/sigs0));
  Real phisp    = frac0*sinim0*SQR(rgc0*omega0);

  switch (imode)
  {
    case 0 : /* explicit numbers */
    {
      d0       = pin->GetReal("problem","d0"  );      /* density. Linked to surface density... */
      vx0      = pin->GetReal("problem","vx0" );      /* flow velocity */
      vy0      = pin->GetReal("problem","vy0" );      /* transverse 1 (usually 0) */
      vz0      = pin->GetReal("problem","vz0" );      /* transverse 2 (usually 0) */
      if (MAGNETIC_FIELDS_ENABLED) {
        bx0      = pin->GetReal("problem","bx0");       /* field along inflow */
        by0      = pin->GetReal("problem","by0");       /* field perp to inflow 1 */
        bz0      = pin->GetReal("problem","bz0");       /* field perp to inflow 2 */
        bpow     = pin->GetReal("problem","bpow");      /* B propto n^bpow */
      }
    }; break;
    case 1 : /* dimensionless initialization (Kim \& Ostriker) */
    {
      qtoom0   = pin->GetReal("problem","qtoom0");    /* Toomre Q */
      beta0    = pin->GetReal("problem","beta0");     /* inv of plasma beta */
      kap0     = sqrt(2.0)*omega0;
      surfn0   = kap0*csound/(pi*qtoom0);       /* now reconstruct field variables */
      vx0      = 0.5*rgc0*omega0*sini0; /* see Kim & Ostriker 2002, footnote 2 */
      vy0      = 0.0; /* set to (0.5*rgc0*omega0-omega0*x1) in loop below */;
      if (MAGNETIC_FIELDS_ENABLED) {
        bx0      = 0.0;
        by0      = 0.0;
        bz0      = 0.0; /* set to sqrt(csound*csound*d0*beta0) in loop below */
        bpow     = pin->GetReal("problem","bpow");      /* B propto n^bpow */
      }
    }; break;
    default: 
    {
      std::stringstream msg;
      msg << "### FATAL ERROR in spiralarm.cpp ProblemGenerator" << std::endl
          << "invalid imode: " << imode << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    qtoom0 = kap0*csound/(pi*surfn0);
    frac0  = phisp/(SQR(rgc0*omega0)*sinim0);
  }
  // set global variables for hydrostatic potential etc.
  cphi0_g  = cphi0;
  phisp_g  = phisp;
  x1len_g  = x1len;
  vx0_g    = vx0;
  ramp_g   = ramp;

  std::cout << "cphi0_g=" << cphi0_g << " phisp_g=" << phisp_g << " x1len_g=" << x1len_g << " vx0_g=" << vx0_g << std::endl;

  // Get hydrostatic equilibrium
  // and copy global profile into local grid. All processes do full 
  // integration (doesn't cost much).
  std::cout << "nz/2+1 = " << nz/2+1 << " nzloc = " << nzloc << std::endl; 
  dprof.NewAthenaArray(nz/2+1);
  dprofcmp.NewAthenaArray(nz+1);
  dprofloc.NewAthenaArray(nzloc);
  Real n00 = find_dprof(dprof,nz/2+1,zmin,zmax,csound,surfn0,surfs0,sigs0,beta0);
  for (iz=0; iz<=nz/2; iz++) {
    dprofcmp(nz/2+iz) = dprof(iz);
    dprofcmp(nz/2-iz) = dprof(iz);
  }

  izloc = (int) (zloc0-zmin)/(zmax-zmin)*((Real) nz); // starting position
  for (iz=0; iz<nzloc; iz++) {
    dprofloc(iz) = 0.5*(dprofcmp(iz+izloc)+dprofcmp(iz+izloc+1));
  }

  // set global variables for boundary conditions
  imode_bc = imode;
  vx0_bc   = vx0;
  vy0_bc   = vy0;
  bx0_bc   = bx0;
  by0_bc   = by0;
  bz0_bc   = bz0;
  omega0_bc= omega0;
  rgc0_bc  = rgc0;
  T0_bc    = T0;
  gm1_bc   = gm1;
  n00_bc   = n00;
  csound_bc= csound;
  beta0_bc = beta0;
  bpow_bc  = bpow;

  // fill main grid: vertically stratified disk with spiral arm potential
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        if (nx3 > 1) {
          iz = k;
        } else {
          iz = j;
        }
        phydro->u(IDN,k,j,i)  = dprofloc(iz-izs);
        switch (imode)
        {
          case 0:
          {
            phydro->u(IM1,k,j,i) = dprofloc(iz-izs)*vx0;
            phydro->u(IM2,k,j,i) = dprofloc(iz-izs)*vy0;
            phydro->u(IM3,k,j,i) = 0.0;
          } break;
          case 1:
          {
            phydro->u(IM1,k,j,i) = dprofloc(iz-izs)*vx0;
            phydro->u(IY ,k,j,i) = dprofloc(iz-izs)*(0.5*rgc0*omega0-omega0*pcoord->x1v(i));
            phydro->u(IZ ,k,j,i)  = 0.0;
          } break;
        } 
        if (NON_BAROTROPIC_EOS) {
          //phydro->u(IEN,k,j,i) = p0/gm1 + 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) = T0*dprofloc(iz-izs)/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))+SQR(phydro->u(IM3,k,j,i)))/dprofloc(iz-izs);
        }
        if (DUAL_ENERGY) {
          phydro->u(IIE,k,j,i) = T0*dprofloc(iz-izs)/gm1;
        }
      }
    }
  }

  // initialize interface B and total energy
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie+1; ++i) {
          pfield->b.x1f(k,j,i) = bx0;
        }
      }
    }
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je+1; ++j) {
        for (int i = is; i <= ie; ++i) {
          if (nx3 > 1) {
            iz = k;
            pfield->b.x2f(k,j,i) = by0*sqrt(csound*csound*n00*beta0)*pow(dprofloc(iz-izs)/n00,bpow);
          } else {
            iz = j;
            pfield->b.x2f(k,j,i) = by0;
          }
        }
      }
    }
    for (int k = ks; k <= ke+1; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          if (nx3 > 1) {
            iz = k;
            pfield->b.x3f(k,j,i) = bz0;
          } else {
            iz = j;
            pfield->b.x3f(k,j,i) = bz0*sqrt(csound*csound*n00*beta0)*pow(dprofloc(iz-izs)/n00,bpow);
          }
        }
      }
    }
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          phydro->u(IEN,k,j,i) += 0.5*(SQR(pfield->b.x1f(k,j,i))+SQR(pfield->b.x2f(k,j,i))+SQR(pfield->b.x3f(k,j,i)));
        }
      }
    }
  }
  dprofcmp.DeleteAthenaArray();
  dprof.DeleteAthenaArray();

  // set the mean density
  Real dmean = GetMeanDensity(pcoord, phydro, is, ie, js, je, ks, je);
  pmy_mesh->SetMeanDensity(dmean);

  return;
}

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief User-defined actions. No source terms!!
//========================================================================================
//

void MeshBlock::UserWorkInLoop(void) {
  Real d0 = GetMeanDensity(pcoord, phydro, is, ie, js, je, ks, je);
  pmy_mesh->SetMeanDensity(d0);
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief Check radius of sphere to make sure it is round
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  return;
}

//========================================================================================
//! \fn void SpiralInnerX1(...)
//  \brief Inner (inflow) X1 boundary conditions for spiral arm
//========================================================================================

void SpiralInflowInnerX1(MeshBlock *pmb, Coordinates *pco,
    AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js,
    int je, int ks, int ke, int ngh) {

  int iz,izs,IY,IZ;
  int nx3 = ke-ks+1;
  if (nx3 > 1) {
    izs = ks;
    IY  = IM2;
    IZ  = IM3;
  } else {
    izs = js;
    IY  = IM3;
    IZ  = IM2;
  }

  // Assign fields 
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=0; i<ngh; i++) {
          b.x1f(k,j,i) = bx0_bc;
        }
      }
    }
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=0; i<ngh; i++) {
          if (nx3 > 1) {
            iz = k;
            b.x2f(k,j,i) = by0_bc*sqrt(csound_bc*csound_bc*n00_bc*beta0_bc)*pow(dprofloc(iz-izs)/n00_bc,bpow_bc);
          } else {
            iz = j;
            b.x2f(k,j,i) = by0_bc;
          }
        }
      }
    }
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=0; i<ngh; i++) {
          if (nx3 > 1) {
            iz = k;
            b.x3f(k,j,i) = bz0_bc;
          } else {
            iz = j;
            b.x3f(k,j,i) = bz0_bc*sqrt(csound_bc*csound_bc*n00_bc*beta0_bc)*pow(dprofloc(iz-izs)/n00_bc,bpow_bc);
          }
        }
      }
    }
  }
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=0; i<ngh; i++) {
        if (nx3 > 1) {
          iz = k;
        } else {
          iz = j;
        }
        pmb->phydro->u(IDN,k,j,i)  = dprofloc(iz-izs);
        switch (imode_bc)
        {
          case 0:
          {
            pmb->phydro->u(IM1,k,j,i) = dprofloc(iz-izs)*vx0_bc;
            pmb->phydro->u(IY ,k,j,i) = dprofloc(iz-izs)*vy0_bc;
            pmb->phydro->u(IZ ,k,j,i) = 0.0;
          } break;
          case 1:
          {
            pmb->phydro->u(IM1,k,j,i) = dprofloc(iz-izs)*vx0_bc;
            pmb->phydro->u(IY ,k,j,i) = dprofloc(iz-izs)*(0.5*rgc0_bc*omega0_bc-omega0_bc*pmb->pcoord->x1v(i));
            pmb->phydro->u(IZ ,k,j,i) = 0.0;
          } break;
        }
        if (NON_BAROTROPIC_EOS) {
          pmb->phydro->u(IEN,k,j,i) =  T0_bc*dprofloc(iz-izs)/gm1_bc 
                                     + 0.5*( SQR(pmb->phydro->u(IM1,k,j,i))
                                            +SQR(pmb->phydro->u(IM2,k,j,i))
                                            +SQR(pmb->phydro->u(IM3,k,j,i)))/dprofloc(iz-izs);
        }
        if (DUAL_ENERGY) {
          pmb->phydro->u(IIE,k,j,i) = T0_bc*dprofloc(iz-izs)/gm1_bc;
        }
        if (MAGNETIC_FIELDS_ENABLED) {
          pmb->phydro->u(IEN,k,j,i) += 0.5*( SQR(b.x1f(k,j,i))
                                            +SQR(b.x2f(k,j,i))
                                            +SQR(b.x3f(k,j,i)));
        }
      }
    }
  }
  return;
}

//========================================================================================
//! \fn void SpiralOutflowOuterX1(...)
//  \brief Outer (outflow) X1 boundary conditions for spiral arm
//========================================================================================

void SpiralOutflowOuterX1(MeshBlock *pmb, Coordinates *pco,
    AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js,
    int je, int ks, int ke, int ngh) {

  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=1; i<=ngh; i++) {
          b.x1f(k,j,ie+i) = b.x1f(k,j,ie);
        }
      }
    }
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=1; i<=ngh; i++) {
          b.x2f(k,j,ie+i) = b.x2f(k,j,ie);
        }
      }
    }     
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=1; i<=ngh; i++) {
          b.x3f(k,j,ie+i) = b.x3f(k,j,ie);
        }        
      }
    }
  }
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=1; i<=ngh; i++) {
      }
    }
  }
  return;
}


