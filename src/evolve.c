#include "sch_1d.h"

void print_matrix(gsl_matrix_complex *M, struct run_param *this_run)
{
  for(int32_t i=0;i<this_run->nmesh_x;i++) {
    for(int32_t j=0;j<this_run->nmesh_x;j++) {
      gsl_complex tmp = gsl_matrix_complex_get(M, i, j);
      printf("%5.3f+(%+5.3f)*I ", GSL_REAL(tmp), GSL_IMAG(tmp));
    }
    printf("\n");
  }
}

void print_vector(gsl_vector_complex *v, struct run_param *this_run)
{
  for(int32_t i=0;i<this_run->nmesh_x;i++) {
    gsl_complex tmp = gsl_vector_complex_get(v,i);
    printf("%5.3f+(%+5.3f)*I\n", GSL_REAL(tmp), GSL_IMAG(tmp));
  }    
}

void evolve_5pnt(double complex *psi, struct run_param *this_run, double dt)
{
  // (I-T)*\Psi^{n+1) = (I+T)*\Psi^n
  // (I+MT)*\Psi^(n+1) = (I-MT)*\Psi^n with MT=-T

  static bool first_call = true;
  
  double rho = dt/SQR(this_run->delta_x);
  // diag and off-diag component of the matrix T
  double diag_imag = -30.0*this_run->hbar*rho/48.0;
  double off_diag1_imag = 16.0*this_run->hbar*rho/48.0;
  double off_diag2_imag = -1.0*this_run->hbar*rho/48.0;

  static gsl_matrix_complex *MT;
  static gsl_vector_complex *rhs, *rhs_, *lhs, *psi_new;
  static gsl_permutation *perm;

  if(first_call) {
    MT = gsl_matrix_complex_alloc(this_run->nmesh_x, this_run->nmesh_x);
    //    gsl_matrix_complex_set_zero(MT);
    
    psi_new = gsl_vector_complex_alloc(this_run->nmesh_x);
    rhs = gsl_vector_complex_alloc(this_run->nmesh_x);

    perm = gsl_permutation_alloc(this_run->nmesh_x);

    //    MT_ = gsl_matrix_complex_alloc(this_run->nmesh_x, this_run->nmesh_x);
    rhs_ = gsl_vector_complex_alloc(this_run->nmesh_x);
    //    lhs = gsl_vector_complex_alloc(this_run->nmesh_x);    

    first_call = false;
  }

  gsl_permutation_init(perm);
  gsl_matrix_complex_set_zero(MT);

  for(int32_t i=0;i<this_run->nmesh_x;i++) {
    gsl_matrix_complex_set(MT, i, i, gsl_complex_rect(0.0, -diag_imag));

    int32_t ip1=(i+1) % this_run->nmesh_x;
    int32_t im1=(i-1+this_run->nmesh_x) % this_run->nmesh_x;
    gsl_matrix_complex_set(MT, i, ip1, gsl_complex_rect(0.0, -off_diag1_imag));
    gsl_matrix_complex_set(MT, i, im1, gsl_complex_rect(0.0, -off_diag1_imag));

    int32_t ip2=(i+2) % this_run->nmesh_x;
    int32_t im2=(i-2+this_run->nmesh_x) % this_run->nmesh_x;
    gsl_matrix_complex_set(MT, i, ip2, gsl_complex_rect(0.0, -off_diag2_imag));
    gsl_matrix_complex_set(MT, i, im2, gsl_complex_rect(0.0, -off_diag2_imag));
  }

  // Copy wave function to rsh
  for(int32_t i=0;i<this_run->nmesh_x;i++) {
    gsl_vector_complex_set(rhs_, i,
			   gsl_complex_rect(creal(psi[i]), cimag(psi[i])));
    gsl_vector_complex_set(rhs, i,
			   gsl_complex_rect(creal(psi[i]), cimag(psi[i])));
  }

  // (I-MT)*rhs = -MT*rhs + rhs
  gsl_blas_zgemv(CblasNoTrans,
		 gsl_complex_rect(-1.0, 0.0), MT, rhs_,
		 gsl_complex_rect(1.0, 0.0), rhs);

  // Compute I+MT matrix
  for(int32_t i=0;i<this_run->nmesh_x;i++) {
    gsl_complex MTii = gsl_matrix_complex_get(MT, i, i);
    MTii = gsl_complex_add(MTii, gsl_complex_rect(1.0, 0.0));
    gsl_matrix_complex_set(MT, i, i, MTii);
  }

  int signum;
  gsl_linalg_complex_LU_decomp(MT, perm, &signum);
  gsl_linalg_complex_LU_solve(MT, perm, rhs, psi_new);

  // Copy the solution into the double complex array.
  for(int32_t i=0;i<this_run->nmesh_x;i++) {
    gsl_complex psi_i = gsl_vector_complex_get(psi_new, i);
    psi[i] = GSL_REAL(psi_i) + _Complex_I*GSL_IMAG(psi_i);
  }
}

// time evolution of the schroedinger equation with 3pt approximation
// of the Laplace operator
void evolve_3pnt(double complex *psi, struct run_param *this_run, double dt)
{
  // (I-T)*\Psi^{n+1) = (I+T)*\Psi^n
  // (I+MT)*\Psi^(n+1) = (I-MT)*\Psi^n with MT=-T

  static bool first_call = true;
  
  double rho = dt/SQR(this_run->delta_x);
  // diag and off-diag component of the matrix T
  double diag_imag = -2.0*this_run->hbar*rho/4.0;
  double off_diag_imag = this_run->hbar*rho/4.0;

  static gsl_matrix_complex *MT, *MT_;
  static gsl_vector_complex *rhs, *rhs_, *lhs, *psi_new;
  static gsl_permutation *perm;

  if(first_call) {
    MT = gsl_matrix_complex_alloc(this_run->nmesh_x, this_run->nmesh_x);
    //    gsl_matrix_complex_set_zero(MT);
    
    psi_new = gsl_vector_complex_alloc(this_run->nmesh_x);
    rhs = gsl_vector_complex_alloc(this_run->nmesh_x);

    perm = gsl_permutation_alloc(this_run->nmesh_x);

    //    MT_ = gsl_matrix_complex_alloc(this_run->nmesh_x, this_run->nmesh_x);
    rhs_ = gsl_vector_complex_alloc(this_run->nmesh_x);
    //    lhs = gsl_vector_complex_alloc(this_run->nmesh_x);    

    first_call = false;
  }

  gsl_permutation_init(perm);
  gsl_matrix_complex_set_zero(MT);

  for(int32_t i=0;i<this_run->nmesh_x;i++) {
    gsl_matrix_complex_set(MT, i, i, gsl_complex_rect(0.0, -diag_imag));

    int32_t ip1=(i+1) % this_run->nmesh_x;
    int32_t im1=(i-1+this_run->nmesh_x) % this_run->nmesh_x;
    gsl_matrix_complex_set(MT, i, ip1, gsl_complex_rect(0.0, -off_diag_imag));
    gsl_matrix_complex_set(MT, i, im1, gsl_complex_rect(0.0, -off_diag_imag));
  }

  // Copy wave function to rsh
  for(int32_t i=0;i<this_run->nmesh_x;i++) {
    gsl_vector_complex_set(rhs_, i,
			   gsl_complex_rect(creal(psi[i]), cimag(psi[i])));
    gsl_vector_complex_set(rhs, i,
			   gsl_complex_rect(creal(psi[i]), cimag(psi[i])));
  }

  // (I-MT)*rhs = -MT*rhs + rhs
  gsl_blas_zgemv(CblasNoTrans,
		 gsl_complex_rect(-1.0, 0.0), MT, rhs_,
		 gsl_complex_rect(1.0, 0.0), rhs);

  // compute I+MT matrix
  for(int32_t i=0;i<this_run->nmesh_x;i++) {
    gsl_complex MTii = gsl_matrix_complex_get(MT, i, i);
    MTii = gsl_complex_add(MTii, gsl_complex_rect(1.0, 0.0));
    gsl_matrix_complex_set(MT, i, i, MTii);
  }

#if 0
  print_matrix(MT, this_run);
#endif

#if 0
  gsl_matrix_complex_memcpy(MT_, MT);
  gsl_vector_complex_memcpy(rhs_, rhs);
#endif

  int signum;
  gsl_linalg_complex_LU_decomp(MT, perm, &signum);
  gsl_linalg_complex_LU_solve(MT, perm, rhs, psi_new);

#if 0
  gsl_blas_zgemv(CblasNoTrans,
		 gsl_complex_rect(1.0, 0.0), MT_, psi_new,
		 gsl_complex_rect(0.0, 0.0), lhs);

  for(int32_t i=0;i<this_run->nmesh_x;i++) {
    gsl_complex rhs_i = gsl_vector_complex_get(rhs_, i);
    gsl_complex lhs_i = gsl_vector_complex_get(lhs, i);
    printf("%d %12.4e+I*(%12.4e), %12.4e+I*(%12.4e)\n",
	   i, GSL_REAL(rhs_i), GSL_IMAG(rhs_i),
	   GSL_REAL(lhs_i), GSL_IMAG(lhs_i));


  }
#endif

  // copy the solution into the double complex array.
  for(int32_t i=0;i<this_run->nmesh_x;i++) {
    gsl_complex psi_i = gsl_vector_complex_get(psi_new, i);
    psi[i] = GSL_REAL(psi_i) + _Complex_I*GSL_IMAG(psi_i);
  }
  
}

