#include "sch_1d.h"
#include <mkl.h>
#include <fstream>
#include <cstring>




// Index helper
// ab is 1D array, COL_MAJOR, ldab=7
// range: max(0,j-ku) <= i <= min(N-1,j+kl)
static inline int ab_idx(int i, int j, int ldab, int kl, int ku) {
    return (kl + ku + i - j) + j * ldab;
}


// ============================================================
//  evolve_5pnt_band: Crank-Nicolson method
//  kl=2, ku=2 (5 points)
//  LAPACK band matrix:  ab[kl+ku+i-j + j*ldab]  (COL_MAJOR)
//  ldab = 2*kl + ku + 1 = 7
// ============================================================
void evolve_5pnt(complexd *psi, double *pot, run_param &tr, double dt)
{
  // (I-T)*\Psi^{n+1) = (I+T)*\Psi^n
  // (I+MT)*\Psi^(n+1) = (I-MT)*\Psi^n with MT=-T

  static bool   first_call = true;
  static complexd *ab;       // band matrix (I+MT)
  static complexd *rhs;      // 
  static lapack_int *ipiv;   // LU pivot

  const int N    = tr.nmesh_x;
  const int kl   = 2;
  const int ku   = 2;
  const int ldab = 2 * kl + ku + 1;   // = 7

  // diagonal and first and second off-diagonal elements of MT
  const complexd diag  (0.0,  30.0 * tr.hbar * tr.rho / 48.0);
  const complexd off1  (0.0, -16.0 * tr.hbar * tr.rho / 48.0);
  const complexd off2  (0.0,         tr.hbar * tr.rho / 48.0);

  if (first_call) {

    ab   = static_cast<complexd*>(
               std::aligned_alloc(64, sizeof(complexd) * ldab * N));
    rhs  = static_cast<complexd*>(
               std::aligned_alloc(64, sizeof(complexd) * N));
    ipiv = static_cast<lapack_int*>(
               std::aligned_alloc(64, sizeof(lapack_int) * N));

    // initialize ab
    std::memset(ab, 0, sizeof(complexd) * ldab * N);

    // construct ( I + MT ) in band matrix form.
#pragma omp parallel for schedule(static)
    for (int j = 0; j < N; j++) {
      complexd pot_term(0.0, 0.5 * pot[j] * dt / tr.hbar);
      ab[ab_idx(j, j, ldab, kl, ku)] = 1.0 + diag + pot_term;
    }

    // 1st off diag (upper band ku=1, lower band kl=1 )
#pragma omp parallel for schedule(static)
    for (int j = 0; j < N - 1; j++) {
      ab[ab_idx(j,   j+1, ldab, kl, ku)] = off1; // i=j,  col=j+1
      ab[ab_idx(j+1, j,   ldab, kl, ku)] = off1; // i=j+1, col=j
    }

    // 2nd off diag (upper band ku=2, lower band kl=2 )
#pragma omp parallel for schedule(static)
    for (int j = 0; j < N - 2; j++) {
      ab[ab_idx(j,   j+2, ldab, kl, ku)] = off2;
      ab[ab_idx(j+2, j,   ldab, kl, ku)] = off2;
    }

    // LU decomposition
    lapack_int info = LAPACKE_zgbtrf(
      LAPACK_COL_MAJOR,
      N, N,
      kl, ku,
      reinterpret_cast<lapack_complex_double*>(ab),
      ldab,
      ipiv
    );
    if (info != 0) {
      fprintf(stderr, "[evolve_5pnt] zgbtrf failed: info=%d\n", info);
      std::abort();
    }

    first_call = false;
  }

  // rhs = (I - MT) * psi
#pragma omp parallel for schedule(static)
  for (int i = 0; i < N; i++) {
    // diag: (1 - diag - i*0.5*dt*pot/hbar) * psi[i]
    complexd pot_term(0.0, 0.5 * pot[i] * dt / tr.hbar);
    complexd val = (1.0 - diag - pot_term) * psi[i];

    // 1st
    if (i > 0)     val -= off1 * psi[i-1];
    if (i < N-1)   val -= off1 * psi[i+1];

    // 2nd
    if (i > 1)     val -= off2 * psi[i-2];
    if (i < N-2)   val -= off2 * psi[i+2];

    rhs[i] = val;
  }

  // LU solution
  lapack_int info = LAPACKE_zgbtrs(
    LAPACK_COL_MAJOR,
    'N',
    N, kl, ku,
    1, 
    reinterpret_cast<const lapack_complex_double*>(ab),
    ldab,
    ipiv,
    reinterpret_cast<lapack_complex_double*>(rhs),
    N
  );
  if (info != 0) {
    fprintf(stderr, "[evolve_5pnt] zgbtrs failed: info=%d\n", info);
    std::abort();
  }

  // ---- update psi
  cblas_zcopy(N, reinterpret_cast<const void*>(rhs), 1,
                 reinterpret_cast<void*>(psi), 1);
}




// This function usea a general matrix library for LU decomposition,
// but it is no longer used because a new fanction was created
// that employs a faster library optimized for band matrices.
void evolve_5pnt_old_v2(complexd *psi, double *pot, run_param & tr, double dt)
{
  // (I-T)*\Psi^{n+1) = (I+T)*\Psi^n
  // (I+MT)*\Psi^(n+1) = (I-MT)*\Psi^n with MT=-T
  
  static bool first_call = true;
  static complexd *A, *MT;
  static complexd *rhs, *rhs_;
  static lapack_int *ipiv;

  complexd coeff(0.0, tr.rho*tr.hbar/48.0);

  // diagonal and first and second off-diagonal elements of MT
  complexd diag(0.0, 30.0*tr.hbar*tr.rho/48.0);
  complexd off_diag1(0.0, -16.0*tr.hbar*tr.rho/48.0);
  complexd off_diag2(0.0, tr.hbar*tr.rho/48.0);

  int32_t info;

  if(first_call) {
    A = static_cast<complexd*>(std::aligned_alloc(64, sizeof(complexd)*SQR(tr.nmesh_x)));
    MT = static_cast<complexd*>(std::aligned_alloc(64, sizeof(complexd)*SQR(tr.nmesh_x)));
    rhs = static_cast<complexd*>(std::aligned_alloc(64, sizeof(complexd)*tr.nmesh_x));
    rhs_ = static_cast<complexd*>(std::aligned_alloc(64, sizeof(complexd)*tr.nmesh_x));

    ipiv = static_cast<lapack_int*>(std::aligned_alloc(64, sizeof(lapack_int)*tr.nmesh_x));

    // zero out
#pragma omp parallel for schedule(auto)
    for(int32_t i=0;i<SQR(tr.nmesh_x);i++) MT[i] = 0.0;

  // setup the matrix MT
#ifdef  __PERIODIC__
#pragma omp parallel for schedule(auto)
    for(int32_t i=0;i<tr.nmesh_x;i++) {
      int32_t diag_indx = i + tr.nmesh_x*i;

      complexd pot_diag(0.0, 0.5*pot[i]*dt/tr.hbar);

      MT[diag_indx] = diag + pot_diag;

      int32_t ip1 = (i+1) % tr.nmesh_x;
      int32_t im1 = (i-1+tr.nmesh_x) % tr.nmesh_x;
      MT[ip1 + tr.nmesh_x*i] = off_diag1;
      MT[im1 + tr.nmesh_x*i] = off_diag1;

      int32_t ip2 = (i+2) % tr.nmesh_x;
      int32_t im2 = (i-2+tr.nmesh_x) % tr.nmesh_x;
      MT[ip2 + tr.nmesh_x*i] = off_diag2;
      MT[im2 + tr.nmesh_x*i] = off_diag2;
    }
#else // !__PERIODIC__
    // diagonal component
#pragma omp parallel for schedule(auto)
    for(int32_t i=0;i<tr.nmesh_x;i++) {
      complexd pot_diag(0.0, 0.5*pot[i]*dt/tr.hbar);
      MT[i + tr.nmesh_x*i] = diag + pot_diag;
    }

    // first off-diagonal component
#pragma omp parallel for schedule(auto)
    for(int32_t i=0;i<tr.nmesh_x-1;i++) MT[i + 1 + tr.nmesh_x*i] = off_diag1;
#pragma omp parallel for schedule(auto)
    for(int32_t i=1;i<tr.nmesh_x;i++) MT[i - 1 + tr.nmesh_x*i] = off_diag1;

    // second off-diagonal component
#pragma omp parallel for schedule(auto)
    for(int32_t i=0;i<tr.nmesh_x-2;i++) MT[i + 2 + tr.nmesh_x*i] = off_diag2;
#pragma omp parallel for schedule(auto)
    for(int32_t i=2;i<tr.nmesh_x;i++) MT[i - 2 + tr.nmesh_x*i] = off_diag2;
#endif
    
    // compute (I+MT) matrix
    //  for(int32_t i=0;i<SQR(tr.nmesh_x);i++) A[i] = MT[i];
    cblas_zcopy(tr.nmesh_x*tr.nmesh_x, MT, 1, A, 1);
#pragma omp parallel for schedule(auto)  
    for(int32_t i=0;i<tr.nmesh_x;i++) {
      int32_t diag_indx = i + tr.nmesh_x*i;
      A[diag_indx] += 1.0;
    }

    // LU-decomposition
    info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, tr.nmesh_x, tr.nmesh_x,
			  reinterpret_cast<MKL_Complex16*>(A), tr.nmesh_x, ipiv);

    first_call = false;
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


  info = LAPACKE_zgetrs(LAPACK_ROW_MAJOR, 'N', tr.nmesh_x, 1,
			reinterpret_cast<MKL_Complex16*>(A), tr.nmesh_x, ipiv,
			reinterpret_cast<MKL_Complex16*>(rhs), 1);
#if 0
  // solve the linear equation 
  LAPACKE_zgesv(LAPACK_ROW_MAJOR, tr.nmesh_x, 1,
		reinterpret_cast<MKL_Complex16*>(A), tr.nmesh_x,
		ipiv, reinterpret_cast<MKL_Complex16*>(rhs), 1);
#endif

  //copy the solution into the array psi[*]
#pragma omp parallel for schedule(auto)
  for(int32_t i=0;i<tr.nmesh_x;i++) psi[i] = rhs[i];

}



void evolve_5pnt_old(complexd *psi, run_param & tr, double dt)
{
  // (I-T)*\Psi^{n+1) = (I+T)*\Psi^n
  // (I+MT)*\Psi^(n+1) = (I-MT)*\Psi^n with MT=-T

  static bool first_call = true;
  static complexd *MT;
  static complexd *rhs, *rhs_;
  static lapack_int *ipiv;

  complexd coeff(0.0, tr.rho*tr.hbar/48.0);

  // diagonal and first and second off-diagonal elements of MT
  complexd diag(0.0, 30.0*tr.hbar*tr.rho/48.0);
  complexd off_diag1(0.0, -16.0*tr.hbar*tr.rho/48.0);
  complexd off_diag2(0.0, tr.hbar*tr.rho/48.0);

  if(first_call) {
    MT = static_cast<complexd*>(std::aligned_alloc(64, sizeof(complexd)*SQR(tr.nmesh_x)));
    rhs = static_cast<complexd*>(std::aligned_alloc(64, sizeof(complexd)*tr.nmesh_x));
    rhs_ = static_cast<complexd*>(std::aligned_alloc(64, sizeof(complexd)*tr.nmesh_x));

    ipiv = static_cast<lapack_int*>(std::aligned_alloc(64, sizeof(lapack_int)*tr.nmesh_x));

    first_call = false;
  }

  // zero out
#pragma omp parallel for schedule(auto)
  for(int32_t i=0;i<SQR(tr.nmesh_x);i++) MT[i] = 0.0;

  // setup the matrix MT
#ifdef  __PERIODIC__
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
#else // !__PERIODIC__
  // diagonal component
#pragma omp parallel for schedule(auto)
  for(int32_t i=0;i<tr.nmesh_x;i++) MT[i + tr.nmesh_x*i] = diag;

  // first off-diagonal component
#pragma omp parallel for schedule(auto)
  for(int32_t i=0;i<tr.nmesh_x-1;i++) MT[i + 1 + tr.nmesh_x*i] = off_diag1;
#pragma omp parallel for schedule(auto)
  for(int32_t i=1;i<tr.nmesh_x;i++) MT[i - 1 + tr.nmesh_x*i] = off_diag1;

  // second off-diagonal component
#pragma omp parallel for schedule(auto)
  for(int32_t i=0;i<tr.nmesh_x-2;i++) MT[i + 2 + tr.nmesh_x*i] = off_diag2;
#pragma omp parallel for schedule(auto)
  for(int32_t i=2;i<tr.nmesh_x;i++) MT[i - 2 + tr.nmesh_x*i] = off_diag2;
#endif

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
