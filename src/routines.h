#ifndef ROUTINES_H
#define ROUTINES_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <assert.h>

#include "../cx_sparse/Include/cs.h"
#include "stdbool.h"
#include "ac_spec.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ABS(x)    ((x) < 0.0 ? (-x) : (x))

#define COMPLEX_ZERO(z)    ((GSL_REAL(z) == 0.0 && GSL_IMAG(z) == 0.0))
#define CS_COMPLEX_NEG(z)  __cs_complex_neg(z)
#define CS_COMPLEX_CONJ(z) __cs_complex_conj(z)

double dot_product(double *x, double *y, int n);
gsl_complex complex_dot_product(gsl_vector_complex *x, gsl_vector_complex *y, int n);
double norm2(double *x, int n);
double complex_norm2(gsl_vector_complex *x, int n);
void axpy(double *dest, double a, double *x, double *y, int n);
void complex_axpy(gsl_vector_complex *dest, gsl_complex a, gsl_vector_complex *x, gsl_vector_complex *y, int n);
void mat_vec_mul(double *Ax, double **A, double *x, int n);
void mat_vec_mul_trans(double *Ax, double **A, double *x, int n);
void cs_mat_vec_mul(double *dest, cs *A, double *x);
void cs_mat_vec_mul_trans(double *dest, cs *A, double *x);
void complex_cs_mat_vec_mul(gsl_vector_complex *Ax, cs_ci *A, gsl_vector_complex *x);
void complex_cs_mat_vec_mul_herm(gsl_vector_complex *Ax, cs_ci *A, gsl_vector_complex *x);
void jacobi_precond(double *M, double **A, cs *C, int n, bool SPARSE);
void complex_jacobi_precond(gsl_vector_complex *M, gsl_matrix_complex *A, cs_ci *C, int n, bool SPARSE);
void precond_solve(double *M_fin, double *M, double *x, int n);
void complex_precond_solve();
void sub_vector(double *dest, double *x, double *y, int n);
void add_vector(double *dest, double *x, double *y, int n);
void zero_out_vector(double *x, int dimension);
void zero_out_matrix(double **matrix, int row, int col);
void set_vec_val(double *x, double val, int dimension);
ac_t rect_to_polar(gsl_complex z);
double to_degrees(double radians);
double complex_abs(gsl_complex z);
void real_to_gsl_complex_vector(gsl_vector_complex *x_complex, double *x, int dimension);
void real_to_cs_complex_vector(cs_complex_t *x_complex, double *x, int dimension);
void vector_conjugate(gsl_vector_complex *M_conj, gsl_vector_complex *M, int dimension);
double **init_array(int row, int col);
double *init_vector(int row);
gsl_matrix_complex *init_gsl_complex_array(int row, int col);
gsl_vector_complex *init_gsl_complex_vector(int row);
double *init_val_vector(int row, double val);
gsl_permutation *init_permutation(int dimension);
void gsl_to_cs_complex(cs_complex_t *dst, gsl_vector_complex *src, int dimension);
void cs_complex_to_gsl(gsl_vector_complex *dst, cs_complex_t *src, int dimension);
cs_complex_t pol_to_rect(double magnitude, double phase);
gsl_complex __gsl_complex_conj(gsl_complex x);
gsl_complex __gsl_complex_neg(gsl_complex x);
cs_complex_t __cs_complex_neg(cs_complex_t x);
cs_complex_t __cs_complex_conj(cs_complex_t x);

#endif