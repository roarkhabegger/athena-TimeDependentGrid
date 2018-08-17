//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reflect.cpp
//  \brief implementation of reflecting BCs in each dimension

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "bvals.hpp"

//----------------------------------------------------------------------------------------
//! \fn void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, const Real time, const Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, inner x1 boundary

void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy hydro variables into ghost zones, reflecting v1
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVX)) {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=(NGHOST); ++i) {
          prim(IVX,k,j,is-i) = -prim(IVX,k,j,(is+i-1));  // reflect 1-velocity
        }
      }}
    } else {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=(NGHOST); ++i) {
          prim(n,k,j,is-i) = prim(n,k,j,(is+i-1));
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b1
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x1f(k,j,(is-i)) = -b.x1f(k,j,(is+i  ));  // reflect 1-field
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,(is-i)) =  b.x2f(k,j,(is+i-1));
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x3f(k,j,(is-i)) =  b.x3f(k,j,(is+i-1));
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, const Real time, const Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, outer x1 boundary

void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy hydro variables into ghost zones, reflecting v1
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVX)) {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=(NGHOST); ++i) {
          prim(IVX,k,j,ie+i) = -prim(IVX,k,j,(ie-i+1));  // reflect 1-velocity
        }
      }}
    } else {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=(NGHOST); ++i) {
          prim(n,k,j,ie+i) = prim(n,k,j,(ie-i+1));
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b1
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x1f(k,j,(ie+i+1)) = -b.x1f(k,j,(ie-i+1));  // reflect 1-field
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,(ie+i  )) =  b.x2f(k,j,(ie-i+1));
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x3f(k,j,(ie+i  )) =  b.x3f(k,j,(ie-i+1));
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReflecInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                         FaceField &b, const Real time, const Real dt,
//                         int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, inner x2 boundary

void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy hydro variables into ghost zones, reflecting v2
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVY)) {
      for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(IVY,k,js-j,i) = -prim(IVY,k,js+j-1,i);  // reflect 2-velocity
        }
      }}
    } else {
      for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(n,k,js-j,i) = prim(n,k,js+j-1,i);
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(js-j),i) =  b.x1f(k,(js+j-1),i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(js-j),i) = -b.x2f(k,(js+j  ),i);  // reflect 2-field
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(js-j),i) =  b.x3f(k,(js+j-1),i);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, const Real time, const Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, outer x2 boundary

void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy hydro variables into ghost zones, reflecting v2
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVY)) {
      for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(IVY,k,je+j,i) = -prim(IVY,k,je-j+1,i);  // reflect 2-velocity
        }
      }}
    } else {
      for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(n,k,je+j,i) = prim(n,k,je-j+1,i);
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(je+j  ),i) =  b.x1f(k,(je-j+1),i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(je+j+1),i) = -b.x2f(k,(je-j+1),i);  // reflect 2-field
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(je+j  ),i) =  b.x3f(k,(je-j+1),i);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, const Real time, const Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, inner x3 boundary

void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy hydro variables into ghost zones, reflecting v3
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVZ)) {
      for (int k=1; k<=(NGHOST); ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(IVZ,ks-k,j,i) = -prim(IVZ,ks+k-1,j,i);  // reflect 3-velocity
        }
      }}
    } else {
      for (int k=1; k<=(NGHOST); ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(n,ks-k,j,i) = prim(n,ks+k-1,j,i);
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b3
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f((ks-k),j,i) =  b.x1f((ks+k-1),j,i);
      }
    }}

    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x2f((ks-k),j,i) =  b.x2f((ks+k-1),j,i);
      }
    }}

    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x3f((ks-k),j,i) = -b.x3f((ks+k  ),j,i);  // reflect 3-field
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, const Real time, const Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, outer x3 boundary

void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy hydro variables into ghost zones, reflecting v3
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVZ)) {
      for (int k=1; k<=(NGHOST); ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(IVZ,ke+k,j,i) = -prim(IVZ,ke-k+1,j,i);  // reflect 3-velocity
        }
      }}
    } else {
      for (int k=1; k<=(NGHOST); ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(n,ke+k,j,i) = prim(n,ke-k+1,j,i);
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b3
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f((ke+k  ),j,i) =  b.x1f((ke-k+1),j,i);
      }
    }}

    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
        b.x2f((ke+k  ),j,i) =  b.x2f((ke-k+1),j,i);
      }
    }}

    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x3f((ke+k+1),j,i) = -b.x3f((ke-k+1),j,i);  // reflect 3-field
      }
    }}
  }

  return;
}

// Cless versions of reflecting BC  
//----------------------------------------------------------------------------------------
//! \fn void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, Real time, Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, inner x1 boundary

void ReflectInnerCLX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy cless variables into ghost zones
  for (int n=0; n<(NCLESS); ++n) {
		if (n==(IVX)) {
	    for (int k=ks; k<=ke; ++k) {
  	  for (int j=js; j<=je; ++j) {
#pragma omp simd
	      for (int i=1; i<=(NGHOST); ++i) {
		      prim(IVX,k,j,is-i) = -prim(IVX,k,j,(is+i-1));
			  }
		 }}
		}
		else if (n==(IP12)) {
	    for (int k=ks; k<=ke; ++k) {
  	  for (int j=js; j<=je; ++j) {
#pragma omp simd
	      for (int i=1; i<=(NGHOST); ++i) {
		      prim(IP12,k,j,is-i) = -prim(IP12,k,j,(is+i-1));
			  }
		 }}
		}
		else if (n==(IP13)) {
	    for (int k=ks; k<=ke; ++k) {
  	  for (int j=js; j<=je; ++j) {
#pragma omp simd
	      for (int i=1; i<=(NGHOST); ++i) {
		      prim(IP13,k,j,is-i) = -prim(IP13,k,j,(is+i-1));
			  }
		 }}
		}
		else {
	    for (int k=ks; k<=ke; ++k) {
  	  for (int j=js; j<=je; ++j) {
#pragma omp simd
	      for (int i=1; i<=(NGHOST); ++i) {
		      prim(n,k,j,is-i) = prim(n,k,j,(is+i-1));
			  }
		 }}
		}
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                         FaceField &b, Real time, Real dt,
//                         int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, outer x1 boundary

void ReflectOuterCLX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy cless variables into ghost zones
  for (int n=0; n<(NCLESS); ++n) {
		if (n==(IVX)) {
	    for (int k=ks; k<=ke; ++k) {
  	  for (int j=js; j<=je; ++j) {
#pragma omp simd
	      for (int i=1; i<=(NGHOST); ++i) {
		      prim(IVX,k,j,ie+i) = -prim(IVX,k,j,(ie-i+1));
			  }
		 }}
		}
		else if (n==(IP12)) {
	    for (int k=ks; k<=ke; ++k) {
  	  for (int j=js; j<=je; ++j) {
#pragma omp simd
	      for (int i=1; i<=(NGHOST); ++i) {
		      prim(IP12,k,j,ie+i) = -prim(IP12,k,j,(ie-i+1));
			  }
		 }}
		}
		else if (n==(IP13)) {
	    for (int k=ks; k<=ke; ++k) {
  	  for (int j=js; j<=je; ++j) {
#pragma omp simd
	      for (int i=1; i<=(NGHOST); ++i) {
		      prim(IP13,k,j,ie+i) = -prim(IP13,k,j,(ie-i+1));
			  }
		 }}
		}
		else {
	    for (int k=ks; k<=ke; ++k) {
  	  for (int j=js; j<=je; ++j) {
#pragma omp simd
	      for (int i=1; i<=(NGHOST); ++i) {
		      prim(n,k,j,ie+i) = prim(n,k,j,(ie-i+1));
			  }
		 }}
		}
  }
	return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, Real time, Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, inner x2 boundary

void ReflectInnerCLX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy cless variables into ghost zones
  for (int n=0; n<(NCLESS); ++n) {
		if (n==(IVY)) {
	    for (int k=ks; k<=ke; ++k) {
  	  for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
	      for (int i=is; i<=ie; ++i) {
		      prim(IVY,k,js-j,i) = -prim(IVY,k,(js+j-1),i);
			  }
		 }}
		}
		else if (n==(IP12)) {
	    for (int k=ks; k<=ke; ++k) {
  	  for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
	      for (int i=is; i<=ie; ++i) {
		      prim(IP12,k,js-j,i) = -prim(IP12,k,(js+j-1),i);
			  }
		 }}
		}
		else if (n==(IP23)) {
	    for (int k=ks; k<=ke; ++k) {
  	  for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
	      for (int i=is; i<=ie; ++i) {
		      prim(IP23,k,js-j,i) = -prim(IP23,k,(js+j-1),i);
			  }
		 }}
		}
		else {
	    for (int k=ks; k<=ke; ++k) {
  	  for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
	      for (int i=is; i<=ie; ++i) {
		      prim(n,k,js-j,i) = prim(n,k,(js+j-1),i);
			  }
		 }}
		}
	}
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, Real time, Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, outer x2 boundary

void ReflectOuterCLX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy cless variables into ghost zones
  for (int n=0; n<(NCLESS); ++n) {
		if (n==(IVY)) {
	    for (int k=ks; k<=ke; ++k) {
  	  for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
	      for (int i=is; i<=ie; ++i) {
		      prim(IVY,k,je+j,i) = -prim(IVY,k,(je-j+1),i);
			  }
		 }}
		}
		else if (n==(IP12)) {
	    for (int k=ks; k<=ke; ++k) {
  	  for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
	      for (int i=is; i<=ie; ++i) {
		      prim(IP12,k,je+j,i) = -prim(IP12,k,(je-j+1),i);
			  }
		 }}
		}
		else if (n==(IP23)) {
	    for (int k=ks; k<=ke; ++k) {
  	  for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
	      for (int i=is; i<=ie; ++i) {
		      prim(IP23,k,je+j,i) = -prim(IP23,k,(je-j+1),i);
			  }
		 }}
		}
		else {
	    for (int k=ks; k<=ke; ++k) {
  	  for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
	      for (int i=is; i<=ie; ++i) {
		      prim(n,k,je+j,i) = prim(n,k,(je-j+1),i);
			  }
		 }}
		}
	}
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, Real time, Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, inner x3 boundary

void ReflectInnerCLX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy cless variables into ghost zones
  for (int n=0; n<(NCLESS); ++n) {
		if (n==(IVZ)) {
	    for (int k=1; k<=(NGHOST); ++k) {
  	  for (int j=js; j<=je; ++j) {
#pragma omp simd
	      for (int i=is; i<=ie; ++i) {
		      prim(IVZ,ks-k,j,i) = -prim(IVZ,(ks+k-1),j,i);
			  }
		 }}
		}
		else if (n==(IP13)) {
	    for (int k=1; k<=(NGHOST); ++k) {
  	  for (int j=js; j<=je; ++j) {
#pragma omp simd
	      for (int i=is; i<=ie; ++i) {
		      prim(IP13,ks-k,j,i) = -prim(IP13,(ks+k-1),j,i);
			  }
		 }}
		}
		else if (n==(IP23)) {
	    for (int k=1; k<=(NGHOST); ++k) {
  	  for (int j=js; j<=je; ++j) {
#pragma omp simd
	      for (int i=is; i<=ie; ++i) {
		      prim(IP23,ks-k,j,i) = -prim(IP23,(ks+k-1),j,i);
			  }
		 }}
		}
		else {
	    for (int k=1; k<=(NGHOST); ++k) {
  	  for (int j=js; j<=je; ++j) {
#pragma omp simd
	      for (int i=is; i<=ie; ++i) {
		      prim(n,ks-k,j,i) = prim(n,(ks+k-1),j,i);
			  }
		 }}
		}
	}  
	return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, Real time, Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, outer x3 boundary

void ReflectOuterCLX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy cless variables into ghost zones
  for (int n=0; n<(NCLESS); ++n) {
		if (n==(IVZ)) {
	    for (int k=1; k<=(NGHOST); ++k) {
  	  for (int j=js; j<=je; ++j) {
#pragma omp simd
	      for (int i=is; i<=ie; ++i) {
		      prim(IVZ,ke+k,j,i) = -prim(IVZ,(ke-k+1),j,i);
			  }
		 }}
		}
		else if (n==(IP13)) {
	    for (int k=1; k<=(NGHOST); ++k) {
  	  for (int j=js; j<=je; ++j) {
#pragma omp simd
	      for (int i=is; i<=ie; ++i) {
		      prim(IP13,ke+k,j,i) = -prim(IP13,(ke-k+1),j,i);
			  }
		 }}
		}
		else if (n==(IP23)) {
	    for (int k=1; k<=(NGHOST); ++k) {
  	  for (int j=js; j<=je; ++j) {
#pragma omp simd
	      for (int i=is; i<=ie; ++i) {
		      prim(IP23,ke+k,j,i) = -prim(IP23,(ke-k+1),j,i);
			  }
		 }}
		}
		else {
	    for (int k=1; k<=(NGHOST); ++k) {
  	  for (int j=js; j<=je; ++j) {
#pragma omp simd
	      for (int i=is; i<=ie; ++i) {
		      prim(n,ke+k,j,i) = prim(n,(ke-k+1),j,i);
			  }
		 }}
		}
	} 
	return;
}
