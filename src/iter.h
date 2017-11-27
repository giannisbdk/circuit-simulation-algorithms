#ifndef ITER_H
#define ITER_H

#include <stdbool.h>

#include "routines.h"
#include "csparse.h"

#define EPSILON 1e-16

int conj_grad(double **A, cs *C, double *x, double *b, double *M, int dimension, double itol, int maxiter, bool SPARSE);
int bi_conj_grad(double **A, cs *C, double *x, double *b, double *M, int dimension, double itol, int maxiter, bool SPARSE);

#endif