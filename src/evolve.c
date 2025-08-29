#include "sch_1d.h"

void print_matrix(gsl_matrix_complex *M, struct run_param *tr)
{
  for(int32_t i=0;i<tr->nmesh_x;i++) {
    for(int32_t j=0;j<tr->nmesh_x;j++) {
      gsl_complex tmp = gsl_matrix_complex_get(M, i, j);
      printf("%5.3f+(%+5.3f)*I ", GSL_REAL(tmp), GSL_IMAG(tmp));
    }
    printf("\n");
  }
}

void print_vector(gsl_vector_complex *v, struct run_param *tr)
{
  for(int32_t i=0;i<tr->nmesh_x;i++) {
    gsl_complex tmp = gsl_vector_complex_get(v,i);
    printf("%5.3f+(%+5.3f)*I\n", GSL_REAL(tmp), GSL_IMAG(tmp));
  }    
}

// time evolution of the schroedinger equation with 5pt approximation
// of the Laplace operator
void evolve_5pnt(double complex *psi, double *pot, struct run_param *tr,
		 double dt)
{
  // (I-T)*\Psi^{n+1) = (I+T)*\Psi^n
  // (I+MT)*\Psi^(n+1) = (I-MT)*\Psi^n with MT=-T

  static bool first_call = true;
  
  double rho = dt/SQR(tr->delta_x);
  // diag and off-diag component of the matrix T
  double diag_imag = -30.0*tr->hbar*rho/48.0;
  double off_diag1_imag = 16.0*tr->hbar*rho/48.0;
  double off_diag2_imag = -1.0*tr->hbar*rho/48.0;

  static gsl_matrix_complex *A, *MT;
  static gsl_vector_complex *rhs, *rhs_, *psi_new;
  static gsl_permutation *perm;

  if(first_call) {
    A  = gsl_matrix_complex_alloc(tr->nmesh_x, tr->nmesh_x);
    MT = gsl_matrix_complex_alloc(tr->nmesh_x, tr->nmesh_x);
    
    psi_new = gsl_vector_complex_alloc(tr->nmesh_x);
    rhs = gsl_vector_complex_alloc(tr->nmesh_x);

    perm = gsl_permutation_alloc(tr->nmesh_x);

    rhs_ = gsl_vector_complex_alloc(tr->nmesh_x);

    gsl_permutation_init(perm);
    gsl_matrix_complex_set_zero(MT);

    for(int32_t i=0;i<tr->nmesh_x;i++) {
      double pot_imag = 0.5*pot[i]*dt/tr->hbar;
      gsl_matrix_complex_set(MT, i, i, gsl_complex_rect(0.0, -diag_imag+pot_imag));

      int32_t ip1=(i+1) % tr->nmesh_x;
      int32_t im1=(i-1+tr->nmesh_x) % tr->nmesh_x;
      gsl_matrix_complex_set(MT, i, ip1, gsl_complex_rect(0.0, -off_diag1_imag));
      gsl_matrix_complex_set(MT, i, im1, gsl_complex_rect(0.0, -off_diag1_imag));
      
      int32_t ip2=(i+2) % tr->nmesh_x;
      int32_t im2=(i-2+tr->nmesh_x) % tr->nmesh_x;
      gsl_matrix_complex_set(MT, i, ip2, gsl_complex_rect(0.0, -off_diag2_imag));
      gsl_matrix_complex_set(MT, i, im2, gsl_complex_rect(0.0, -off_diag2_imag));
    }

    // Compute I-T = I+MT matrix
    gsl_matrix_complex_memcpy(A, MT);
    for(int32_t i=0;i<tr->nmesh_x;i++) {
      gsl_complex Aii = gsl_matrix_complex_get(A, i, i);
      Aii = gsl_complex_add(Aii, gsl_complex_rect(1.0, 0.0));
      gsl_matrix_complex_set(A, i, i, Aii);
    }

    int signum;
    gsl_linalg_complex_LU_decomp(A, perm, &signum);

    first_call = false;
  }


  // Copy wave function to rsh
#pragma omp parallel for schedule(auto)  
  for(int32_t i=0;i<tr->nmesh_x;i++) {
    gsl_vector_complex_set(rhs_, i,
			   gsl_complex_rect(creal(psi[i]), cimag(psi[i])));
    gsl_vector_complex_set(rhs, i,
			   gsl_complex_rect(creal(psi[i]), cimag(psi[i])));
  }

  // (I+T)*rhs = (I-MT)*rhs = -MT*rhs + rhs
  gsl_blas_zgemv(CblasNoTrans,
		 gsl_complex_rect(-1.0, 0.0), MT, rhs_,
		 gsl_complex_rect(1.0, 0.0), rhs);

  gsl_linalg_complex_LU_solve(A, perm, rhs, psi_new);

  // Copy the solution into the double complex array.
#pragma omp parallel for schedule(auto)
  for(int32_t i=0;i<tr->nmesh_x;i++) {
    gsl_complex psi_i = gsl_vector_complex_get(psi_new, i);
    psi[i] = GSL_REAL(psi_i) + _Complex_I*GSL_IMAG(psi_i);
  }
}


// time evolution of the schroedinger equation with 5pt approximation
// of the Laplace operator
void evolve_5pnt_free(double complex *psi, struct run_param *tr, double dt)
{
  // (I-T)*\Psi^{n+1) = (I+T)*\Psi^n
  // (I+MT)*\Psi^(n+1) = (I-MT)*\Psi^n with MT=-T

  static bool first_call = true;
  
  double rho = dt/SQR(tr->delta_x);
  // diag and off-diag component of the matrix T
  double diag_imag = -30.0*tr->hbar*rho/48.0;
  double off_diag1_imag = 16.0*tr->hbar*rho/48.0;
  double off_diag2_imag = -1.0*tr->hbar*rho/48.0;

  static gsl_matrix_complex *A, *MT;
  static gsl_vector_complex *rhs, *rhs_, *psi_new;
  static gsl_permutation *perm;

  if(first_call) {
    A  = gsl_matrix_complex_alloc(tr->nmesh_x, tr->nmesh_x);
    MT = gsl_matrix_complex_alloc(tr->nmesh_x, tr->nmesh_x);
    
    psi_new = gsl_vector_complex_alloc(tr->nmesh_x);
    rhs = gsl_vector_complex_alloc(tr->nmesh_x);

    perm = gsl_permutation_alloc(tr->nmesh_x);

    rhs_ = gsl_vector_complex_alloc(tr->nmesh_x);

    gsl_permutation_init(perm);
    gsl_matrix_complex_set_zero(MT);

    for(int32_t i=0;i<tr->nmesh_x;i++) {
      gsl_matrix_complex_set(MT, i, i, gsl_complex_rect(0.0, -diag_imag));

      int32_t ip1=(i+1) % tr->nmesh_x;
      int32_t im1=(i-1+tr->nmesh_x) % tr->nmesh_x;
      gsl_matrix_complex_set(MT, i, ip1, gsl_complex_rect(0.0, -off_diag1_imag));
      gsl_matrix_complex_set(MT, i, im1, gsl_complex_rect(0.0, -off_diag1_imag));
      
      int32_t ip2=(i+2) % tr->nmesh_x;
      int32_t im2=(i-2+tr->nmesh_x) % tr->nmesh_x;
      gsl_matrix_complex_set(MT, i, ip2, gsl_complex_rect(0.0, -off_diag2_imag));
      gsl_matrix_complex_set(MT, i, im2, gsl_complex_rect(0.0, -off_diag2_imag));
    }

    // Compute I-T = I+MT matrix
    gsl_matrix_complex_memcpy(A, MT);
    for(int32_t i=0;i<tr->nmesh_x;i++) {
      gsl_complex Aii = gsl_matrix_complex_get(A, i, i);
      Aii = gsl_complex_add(Aii, gsl_complex_rect(1.0, 0.0));
      gsl_matrix_complex_set(A, i, i, Aii);
    }

    int signum;
    gsl_linalg_complex_LU_decomp(A, perm, &signum);

    first_call = false;
  }


  // Copy wave function to rsh
#pragma omp parallel for schedule(auto)  
  for(int32_t i=0;i<tr->nmesh_x;i++) {
    gsl_vector_complex_set(rhs_, i,
			   gsl_complex_rect(creal(psi[i]), cimag(psi[i])));
    gsl_vector_complex_set(rhs, i,
			   gsl_complex_rect(creal(psi[i]), cimag(psi[i])));
  }

  // (I+T)*rhs = (I-MT)*rhs = -MT*rhs + rhs
  gsl_blas_zgemv(CblasNoTrans,
		 gsl_complex_rect(-1.0, 0.0), MT, rhs_,
		 gsl_complex_rect(1.0, 0.0), rhs);

  gsl_linalg_complex_LU_solve(A, perm, rhs, psi_new);

  // Copy the solution into the double complex array.
#pragma omp parallel for schedule(auto)
  for(int32_t i=0;i<tr->nmesh_x;i++) {
    gsl_complex psi_i = gsl_vector_complex_get(psi_new, i);
    psi[i] = GSL_REAL(psi_i) + _Complex_I*GSL_IMAG(psi_i);
  }
}

// time evolution of the schroedinger equation with 3pt approximation
// of the Laplace operator
void evolve_3pnt(double complex *psi, double *pot, struct run_param *tr, double dt)
{
  // (I-T)*\Psi^{n+1) = (I+T)*\Psi^n
  // (I+MT)*\Psi^(n+1) = (I-MT)*\Psi^n with MT=-T

  static bool first_call = true;
  
  double rho = dt/SQR(tr->delta_x);
  // diag and off-diag component of the matrix T
  double diag_imag = -2.0*tr->hbar*rho/4.0;
  double off_diag_imag = tr->hbar*rho/4.0;

  static gsl_matrix_complex *A, *MT;
  static gsl_vector_complex *rhs, *rhs_, *psi_new;
  static gsl_permutation *perm;

  if(first_call) {
    A  = gsl_matrix_complex_alloc(tr->nmesh_x, tr->nmesh_x);
    MT = gsl_matrix_complex_alloc(tr->nmesh_x, tr->nmesh_x);
    
    psi_new = gsl_vector_complex_alloc(tr->nmesh_x);
    rhs = gsl_vector_complex_alloc(tr->nmesh_x);

    perm = gsl_permutation_alloc(tr->nmesh_x);

    rhs_ = gsl_vector_complex_alloc(tr->nmesh_x);

    gsl_permutation_init(perm);
    gsl_matrix_complex_set_zero(MT);

    for(int32_t i=0;i<tr->nmesh_x;i++) {
      double pot_imag = 0.5*pot[i]*dt/tr->hbar;
      gsl_matrix_complex_set(MT, i, i, gsl_complex_rect(0.0, -diag_imag+pot_imag));

      int32_t ip1=(i+1) % tr->nmesh_x;
      int32_t im1=(i-1+tr->nmesh_x) % tr->nmesh_x;
      gsl_matrix_complex_set(MT, i, ip1, gsl_complex_rect(0.0, -off_diag_imag));
      gsl_matrix_complex_set(MT, i, im1, gsl_complex_rect(0.0, -off_diag_imag));
    }
    
    // compute I+MT matrix
    gsl_matrix_complex_memcpy(A, MT);
    for(int32_t i=0;i<tr->nmesh_x;i++) {
      gsl_complex Aii = gsl_matrix_complex_get(A, i, i);
      Aii = gsl_complex_add(Aii, gsl_complex_rect(1.0, 0.0));
      gsl_matrix_complex_set(A, i, i, Aii);
    }

    int signum;
    gsl_linalg_complex_LU_decomp(A, perm, &signum);

    first_call = false;
  }


  // Copy wave function to rsh
  for(int32_t i=0;i<tr->nmesh_x;i++) {
    gsl_vector_complex_set(rhs_, i,
			   gsl_complex_rect(creal(psi[i]), cimag(psi[i])));
    gsl_vector_complex_set(rhs, i,
			   gsl_complex_rect(creal(psi[i]), cimag(psi[i])));
  }

  // (I-MT)*rhs = -MT*rhs + rhs
  gsl_blas_zgemv(CblasNoTrans,
		 gsl_complex_rect(-1.0, 0.0), MT, rhs_,
		 gsl_complex_rect(1.0, 0.0), rhs);

  gsl_linalg_complex_LU_solve(A, perm, rhs, psi_new);

  // copy the solution into the double complex array.
  for(int32_t i=0;i<tr->nmesh_x;i++) {
    gsl_complex psi_i = gsl_vector_complex_get(psi_new, i);
    psi[i] = GSL_REAL(psi_i) + _Complex_I*GSL_IMAG(psi_i);
  }
  
}


// time evolution of the schroedinger equation with 3pt approximation
// of the Laplace operator
void evolve_3pnt_free(double complex *psi, struct run_param *tr, double dt)
{
  // (I-T)*\Psi^{n+1) = (I+T)*\Psi^n
  // (I+MT)*\Psi^(n+1) = (I-MT)*\Psi^n with MT=-T

  static bool first_call = true;
  
  double rho = dt/SQR(tr->delta_x);
  // diag and off-diag component of the matrix T
  double diag_imag = -2.0*tr->hbar*rho/4.0;
  double off_diag_imag = tr->hbar*rho/4.0;

  static gsl_matrix_complex *A, *MT;
  static gsl_vector_complex *rhs, *rhs_, *psi_new;
  static gsl_permutation *perm;

  if(first_call) {
    A  = gsl_matrix_complex_alloc(tr->nmesh_x, tr->nmesh_x);
    MT = gsl_matrix_complex_alloc(tr->nmesh_x, tr->nmesh_x);
    
    psi_new = gsl_vector_complex_alloc(tr->nmesh_x);
    rhs = gsl_vector_complex_alloc(tr->nmesh_x);

    perm = gsl_permutation_alloc(tr->nmesh_x);

    rhs_ = gsl_vector_complex_alloc(tr->nmesh_x);

    gsl_permutation_init(perm);
    gsl_matrix_complex_set_zero(MT);

    for(int32_t i=0;i<tr->nmesh_x;i++) {
      gsl_matrix_complex_set(MT, i, i, gsl_complex_rect(0.0, -diag_imag));

      int32_t ip1=(i+1) % tr->nmesh_x;
      int32_t im1=(i-1+tr->nmesh_x) % tr->nmesh_x;
      gsl_matrix_complex_set(MT, i, ip1, gsl_complex_rect(0.0, -off_diag_imag));
      gsl_matrix_complex_set(MT, i, im1, gsl_complex_rect(0.0, -off_diag_imag));
    }
    
    // compute I+MT matrix
    gsl_matrix_complex_memcpy(A, MT);
    for(int32_t i=0;i<tr->nmesh_x;i++) {
      gsl_complex Aii = gsl_matrix_complex_get(A, i, i);
      Aii = gsl_complex_add(Aii, gsl_complex_rect(1.0, 0.0));
      gsl_matrix_complex_set(A, i, i, Aii);
    }

    int signum;
    gsl_linalg_complex_LU_decomp(A, perm, &signum);

    first_call = false;
  }


  // Copy wave function to rsh
  for(int32_t i=0;i<tr->nmesh_x;i++) {
    gsl_vector_complex_set(rhs_, i,
			   gsl_complex_rect(creal(psi[i]), cimag(psi[i])));
    gsl_vector_complex_set(rhs, i,
			   gsl_complex_rect(creal(psi[i]), cimag(psi[i])));
  }

  // (I-MT)*rhs = -MT*rhs + rhs
  gsl_blas_zgemv(CblasNoTrans,
		 gsl_complex_rect(-1.0, 0.0), MT, rhs_,
		 gsl_complex_rect(1.0, 0.0), rhs);

  gsl_linalg_complex_LU_solve(A, perm, rhs, psi_new);

  // copy the solution into the double complex array.
  for(int32_t i=0;i<tr->nmesh_x;i++) {
    gsl_complex psi_i = gsl_vector_complex_get(psi_new, i);
    psi[i] = GSL_REAL(psi_i) + _Complex_I*GSL_IMAG(psi_i);
  }
  
}

