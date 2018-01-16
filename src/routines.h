#ifndef ROUTINES_H
#define ROUTINES_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex_math.h>

#include "../cx_sparse/Include/cs.h"
#include "stdbool.h"
#include "ac_spec.h"

double dot_product(double *x, double *y, int n);
gsl_complex complex_dot_product(gsl_vector_complex *x, gsl_vector_complex *y, int n);
double norm2(double *x, int n);
void axpy(double *dest, double a, double *x, double *y, int n);
void mat_vec_mul(double *Ax, double **A, double *x, int n);
void mat_vec_mul_trans(double *Ax, double **A, double *x, int n);
void cs_mat_vec_mul(double *dest, cs *A, double *x);
void cs_mat_vec_mul_trans(double *dest, cs *A, double *x);
void jacobi_precond(double *M, double **A, cs *C, int n, bool SPARSE);
void precond_solve(double *M_fin, double *M, double *x, int n);
void sub_vector(double *dest, double *x, double *y, int n);
void add_vector(double *dest, double *x, double *y, int n);
void zero_out_vec(double *x, int dimension);
void set_vec_val(double *x, double val, int dimension);
ac_t rect_to_polar(gsl_complex z);
double to_degrees(double radians);

#endif