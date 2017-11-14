#ifndef ITER_H
#define ITER_H

#define EPSILON 1e-16

int conj_grad(double **A, double *x, double *b, double *M, int dimension, double itol, int maxiter);
int bi_conj_grad(double **A, double *x, double *b, double **A_trans, double *M, double *M_trans, int dimension, double itol, int maxiter);

#endif