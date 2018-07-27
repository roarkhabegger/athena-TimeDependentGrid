//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file characteristic.cpp
//  \brief Functions to transform vectors between primitive and characteristic variables

// C++ headers
#include <cmath>

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../eos/eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::LeftEigenmatrixDotVectorCL()
//  \brief Computes inner-product of left-eigenmatrix of Roe's matrix A in the primitive
//  variables and an input vector.  This operation converts primitive to characteristic
//  variables.  The result is returned in the input vector, with the components of the
//  characteristic field stored such that vect(1,i) is in the direction of the sweep.
//
//  The order of the components in the input vector should be:
//     (IDN,IVX,IVY,IVZ,IP11,IP22,IP33,IP12,IP13,IP23)
//  and these are permuted according to the direction specified by the input flag "ivx".
//
// REFERENCES:
// - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new code for
//   astrophysical MHD", ApJS, (2008), Appendix A.  Equation numbers refer to this paper.
// - S. Brown, (1996), PhD thesis
// - N.L. Mitchell, E.I. Vorobyov, and G. Hensler, "Collisionless stellar hydrodynamics as
//	 an efficient alternative to N-body methods", MNRAS 428, 2674-2687 (2013) 

void Reconstruction::LeftEigenmatrixDotVectorCL(MeshBlock *pmb, const int ivx, 
  const int ip12, const int il, const int iu, const AthenaArray<Real> &w,
  AthenaArray<Real> &vect) {
  // permute components of input primitive vector depending on direction
  int ivy  = IVX  + ((ivx -IVX )+1)%3;
  int ivz  = IVX  + ((ivx -IVX )+2)%3;
	// diagnol comp 
	int ip11 = IP11 + ((ivx -IVX )  )%3; 
	int ip22 = IP11 + ((ivx -IVX )+1)%3;
	int ip33 = IP11 + ((ivx -IVX )+2)%3; 
	// off-diag comp 
	int ip13 = IP12 + ((ip12-IP12)+1)%3;
	int ip23 = IP12 + ((ip12-IP12)+2)%3; 

#pragma omp simd
	for (int i=il; i<=iu; ++i) {
		Real asq    = w(ip11,i)/w(IDN,i);
		Real d			= w(IDN,i); 
		Real a	    = std::sqrt(asq);
		Real sqt3   = std::sqrt(3.0);
		Real sqt3a  = sqt3*a; 
		Real inva   = 1.0/a;
		Real invasq = 1.0/asq; 
		Real dsq		= SQR(d);
		Real invd   = 1.0/d;
		Real invdsq = 1.0/dsq; 

		Real P11    = w(ip11, i);
		Real P22		= w(ip22, i);
		Real P33		= w(ip33, i);
		Real P12		= w(ip12, i);
		Real P13		= w(ip13, i);
		Real P23		= w(ip23, i);
		
		// Multiply row of L-eigenmatrix with vector using matrix elements 
		Real v_0 = -0.5*(invd/sqt3a)*vect(ivx,i) + (invasq*invdsq/6.0) * vect(ip11,i); 
		Real v_1 =  0.5*P12*inva*invd*vect(ivx,i) - 0.5*a*vect(ivy,i) 
						 -  0.5*P12*invasq*invdsq*vect(ip11,i) + 0.5*invd*vect(ip12,i);
		Real v_2 =  0.5*P13*inva*invd*vect(ivx,i) - 0.5*a*vect(ivz,0) 
						 -  0.5*P13*invasq*invdsq*vect(ip11,i) + 0.5*invd*vect(ip13,i);
		Real v_3 =  invdsq*vect(IDN,i) - (invasq*invdsq/3.0) * vect(ip11,i); 
		Real v_4 =  (4.0*SQR(P12)-asq*P22*d)*invdsq*vect(IDN,i) + asq*vect(ip22,i)
						 -  2.0*P12*invd*vect(ip12,i); 
		Real v_5 =  (4.0*SQR(P13)-asq*P33*d)*invdsq*vect(IDN,i) + asq*vect(ip33,i)
						 -	2.0*P13*invd*vect(ip13,i); 
		Real v_6 =  (4.0*P12*P13-asq*P23*d)*invdsq*vect(IDN,i) - P13*invd*vect(ip12,i)
						 -  P13*invd*vect(ip13,i) + asq*vect(ip23,i); 
		Real v_7 = -0.5*P13*inva*invd*vect(ivx,i) + 0.5*a*vect(ivz,i) 
						 -  0.5*P13*invasq*invdsq*vect(ip11,i) + 0.5*invd*vect(ip13,i); 		
		Real v_8 = -0.5*P12*inva*invd*vect(ivx,i) + 0.5*a*vect(2,i) 
						 -  0.5*P12*invasq*invdsq*vect(ip11,i) + 0.5*invd*vect(ip12,i); 
		Real v_9 =  0.5*(invd/sqt3a)*vect(ivx,i) + (invasq*invdsq/6.0) * vect(ip11,i);


		// Permute components back into standard order for primitives on output
		vect(0,i) = v_0;
		vect(1,i) = v_1;
		vect(2,i) = v_2;
		vect(3,i) = v_3;
		vect(4,i) = v_4;
		vect(5,i) = v_5;
		vect(7,i) = v_6;
		vect(7,i) = v_7;
		vect(8,i) = v_8;
		vect(9,i) = v_9;
	}
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::RightEigenmatrixDotVectorCL()
//  \brief Computes inner-product of right-eigenmatrix of Roe's matrix A in the primitive
//  variables and an input vector.  This operation converts characteristic to primitive
//  variables.  The result is returned in the input vector. This is for CLESS variables
//
//  The order of the components in the input vector (characteristic fields) should be:
//     (IDN,ivx,ivy,ivz,ip11,ip22,ip33,ip12,ip13,ip23)
//  where the lower-case indices indicate that the characteristic field in the direction
//  of the sweep (designated by the input flag "ivx") is stored first.  On output, the
//  components of velocity are in the standard order used for primitive variables.
//
// REFERENCES:
// - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new code for
//   astrophysical MHD", ApJS, (2008), Appendix A.  Equation numbers refer to this paper.
// - S. Brown, (1996), PhD thesis. 

void Reconstruction::RightEigenmatrixDotVectorCL(MeshBlock *pmb, const int ivx,
  const int ip12, const int il, const int iu, const AthenaArray<Real> &w,
  AthenaArray<Real> &vect) {
  // permute components of output primitive vector depending on direction
  int ivy  = IVX  + ((ivx -IVX )+1)%3;
  int ivz  = IVX  + ((ivx -IVX )+2)%3;
	// diagnol comp 
	int ip11 = IP11 + ((ivx -IVX )  )%3; 
	int ip22 = IP11 + ((ivx -IVX )+1)%3;
	int ip33 = IP11 + ((ivx -IVX )+2)%3; 
	// off-diag comp 
	int ip13 = IP12 + ((ip12-IP12)+1)%3;
	int ip23 = IP12 + ((ip12-IP12)+2)%3; 

#pragma omp simd
	for (int i=il; i<=iu; ++i) {
		Real asq    = w(ip11,i)/w(IDN,i);
		Real d			= w(IDN,i); 
		Real a	    = std::sqrt(asq);
		Real sqt3   = std::sqrt(3.0);
		Real sqt3a  = sqt3*a; 
		Real inva   = 1.0/a;
		Real invasq = 1.0/asq; 
		Real dsq		= SQR(d);
		Real invd   = 1.0/d;
		Real invdsq = 1.0/dsq; 

		Real P11    = w(ip11, i);
		Real P22		= w(ip22, i);
		Real P33		= w(ip33, i);
		Real P12		= w(ip12, i);
		Real P13		= w(ip13, i);
		Real P23		= w(ip23, i);
		
		// Multiply row of R-eigenmatrix with vector using matrix elements 
		// Components of vect() are addressed directly as they are input in permuted order
		Real v_0 = dsq*(vect(0,i) + vect(3,i) + vect(9,i));	
		Real v_1 = -sqt3a*d*vect(0,i) + sqt3a*d*vect(9,i);
		Real v_2 = -sqt3*P12*inva*vect(0,i) - inva*vect(1,i) + inva*vect(8,i)
						 + sqt3*P12*inva*vect(9,i); 	
		Real v_3 = -sqt3*P13*inva*vect(0,i) - inva*vect(2,i) + inva*vect(7,i) 
						 + sqt3*P13*inva*vect(9,i);  
		Real v_4 = 3.0*asq*dsq*(vect(0,i) + vect(9,i));  
		Real v_5 = (2.0*SQR(P12)+asq*P22*d)*invasq*vect(0,i) + 2.0*P12*invasq*vect(1,i)
						 + (asq*P22*d-4.0*SQR(P12))*invasq*vect(3,i) + invasq*vect(4,i) 
						 + 2.0*P12*invasq*vect(8,i) + (2.0*SQR(P12)+asq*P22*d)*invasq*vect(9,i); 
		Real v_6 = (2.0*SQR(P13)+asq*P33*d)*invasq*vect(0,i) + 2.0*P13*invasq*vect(2,i)
						 + (asq*P33*d-4.0*SQR(P13))*invasq*vect(3,i) + invasq*vect(5,i)
						 + 2.0*P13*invasq*vect(7,i) + (2.0*SQR(P13)+asq*P33*d)*invasq*vect(9,i);   
		Real v_7 = 3.0*P12*d*vect(0,i) + d*(vect(1,i) + vect(8,i)) + 3.0*P12*d*vect(9,i);  		
		Real v_8 = 3.0*P13*d*vect(0,i) + d*(vect(2,i) + vect(7,i)) + 3.0*P13*d*vect(9,i);  
		Real v_9 = (2.0*P12*P13+asq*P23*d)*invasq*vect(0,i) + P13*invasq*vect(1,i) 
						 + P12*invasq*vect(2,i) + (asq*P23*d-4.0*P12*P13)*invasq*vect(3,i) 
						 + invasq*vect(6,i) + P12*invasq*vect(7,i) + P13*invasq*vect(8,i) 
						 + (2.0*P12*P13+asq*P23*d)*invasq*vect(9,i);   

		// Permute components back into standard order for primitives on output
		vect(IDN ,i) = v_0;
		vect(ivx ,i) = v_1;
		vect(ivy ,i) = v_2;
		vect(ivz ,i) = v_3;
		vect(ip11,i) = v_4;
		vect(ip22,i) = v_5;
		vect(ip33,i) = v_6;
		vect(ip12,i) = v_7;
		vect(ip13,i) = v_8;
		vect(ip23,i) = v_9;
	}
  return;
}
