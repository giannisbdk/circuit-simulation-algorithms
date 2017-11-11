#ifndef ITER_H
#define ITER_H

void axpy(double *dest, double a, double *x, double *y, int n);
double ddot(double *x, double *y, int n);
void matvec(double *Ax, double **Adata, double *xvect, int n);
void jacobi_precond(double *M, double **Adata, int n);
void prec_solve(double *Minvx, double *Mdata, double *x, int n);
double norm2(double *x, int n);
void sub_vector(double *dest, double *x, double *y, int n);
void add_vector(double *dest, double *x, double *y, int n);
int cg(double *x, double itol, int dimension, int maxiter, double *b, double **A);
void bi_cg();

#endif