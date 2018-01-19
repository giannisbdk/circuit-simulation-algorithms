#ifndef ITER_H
#define ITER_H

#include <stdbool.h>
#include <gsl/gsl_blas.h>

#include "routines.h"
#include "../cx_sparse/Include/cs.h"

#define EPSILON				1e-16
#define MAX_ITER_THRESHOLD 	20

int conj_grad(double **A, cs *C, double *x, double *b, double *M, int dimension, double itol, int maxiter, bool SPARSE);
int bi_conj_grad(double **A, cs *C, double *x, double *b, double *M, int dimension, double itol, int maxiter, bool SPARSE);

int complex_conj_grad(gsl_matrix_complex *A, gsl_vector_complex *x, gsl_vector_complex *b, gsl_vector_complex *M,
                      int dimension, double itol, int maxiter, bool SPARSE);

int complex_bi_conj_grad(gsl_matrix_complex *A, gsl_vector_complex *x, gsl_vector_complex *b, gsl_vector_complex *M,
                         gsl_vector_complex *M_conj, int dimension, double itol, int maxiter, bool SPARSE);

#endif