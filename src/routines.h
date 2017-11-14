#ifndef ROUTINES_H
#define ROUTINES_H

double dot_product(double *x, double *y, int n);
double norm2(double *x, int n);
void axpy(double *dest, double a, double *x, double *y, int n);
void mat_vec_mul(double *Ax, double **A, double *x, int n);
void jacobi_precond(double *M, double **A, int n);
void precond_solve(double *M_fin, double *M, double *x, int n);
void trans_matrix(double **A_trans, double **A, int n);
void sub_vector(double *dest, double *x, double *y, int n);
void add_vector(double *dest, double *x, double *y, int n);
void zero_out_vec(double *x, int dimension);

#endif