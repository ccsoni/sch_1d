#include "sch_1d.h"
#include <mkl.h>

void evolve_5pnt(complexd *psi, run_param & tr, double dt)
{
  // (I-T)*\Psi^{n+1) = (I+T)*\Psi^n
  // (I+MT)*\Psi^(n+1) = (I-MT)*\Psi^n with MT=-T
  
  static bool first_call = true;
  static complexd *MT;
  static complexd *psi_new, *rhs, *rhs_;
  static lapack_int *ipiv;

  complexd coeff(0.0, tr.rho*tr.hbar/48.0);

  // diagonal and first and second off-diagonal elements of MT
  complexd diag(0.0, 30.0*tr.hbar*tr.rho/48.0);
  complexd off_diag1(0.0, -16.0*tr.hbar*tr.rho/48.0);
  complexd off_diag2(0.0, tr.hbar*tr.rho/48.0);

  if(first_call) {
    MT = static_cast<complexd*>(std::aligned_alloc(64, sizeof(complexd)*SQR(tr.nmesh_x)));
    psi_new = static_cast<complexd*>(std::aligned_alloc(64, sizeof(complexd)*tr.nmesh_x));
    rhs = static_cast<complexd*>(std::aligned_alloc(64, sizeof(complexd)*tr.nmesh_x));
    rhs_ = static_cast<complexd*>(std::aligned_alloc(64, sizeof(complexd)*tr.nmesh_x));

    ipiv = static_cast<lapack_int*>(std::aligned_alloc(64, sizeof(lapack_int)*tr.nmesh_x));

    first_call = false;
  }

  // zero out
#pragma omp parallel for schedule(auto)
  for(int32_t i=0;i<SQR(tr.nmesh_x);i++) MT[i] = 0.0;

  // setup the matrix MT
#pragma omp parallel for schedule(auto)
  for(int32_t i=0;i<tr.nmesh_x;i++) {
    int32_t diag_indx = i + tr.nmesh_x*i;

    MT[diag_indx] = diag;

    int32_t ip1 = (i+1) % tr.nmesh_x;
    int32_t im1 = (i-1+tr.nmesh_x) % tr.nmesh_x;
    MT[ip1 + tr.nmesh_x*i] = off_diag1;
    MT[im1 + tr.nmesh_x*i] = off_diag1;

    int32_t ip2 = (i+2) % tr.nmesh_x;
    int32_t im2 = (i-2+tr.nmesh_x) % tr.nmesh_x;
    MT[ip2 + tr.nmesh_x*i] = off_diag2;
    MT[im2 + tr.nmesh_x*i] = off_diag2;
  }

  // copy current wavefunc to rhs
#pragma omp parallel for schedule(auto)  
  for(int32_t i=0;i<tr.nmesh_x;i++) {
    rhs[i] = psi[i];
    rhs_[i] = psi[i];
  }

  // (I-MT)*rhs = -MT*rhs + rhs
  complexd alpha(-1.0, 0.0);
  complexd beta(1.0, 0.0);
  cblas_zgemv(CblasRowMajor, CblasNoTrans, tr.nmesh_x, tr.nmesh_x,
	      &alpha, MT, tr.nmesh_x, rhs_, 1,
	      &beta, rhs, 1);

  // compute (I+MT) matrix
#pragma omp parallel for schedule(auto)  
  for(int32_t i=0;i<tr.nmesh_x;i++) {
    int32_t diag_indx = i + tr.nmesh_x*i;
    MT[diag_indx] += 1.0;
  }

  // solve the linear equation 
  LAPACKE_zgesv(LAPACK_ROW_MAJOR, tr.nmesh_x, 1,
		reinterpret_cast<MKL_Complex16*>(MT), tr.nmesh_x,
		ipiv, reinterpret_cast<MKL_Complex16*>(rhs), 1);

  //copy the solution into the array psi[*]
#pragma omp parallel for schedule(auto)  
  for(int32_t i=0;i<tr.nmesh_x;i++) psi[i] = rhs[i];

}
