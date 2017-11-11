#include "iter.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

void axpy(double *dest, double a, double *x, double *y, int n) {
    for (int i = 0; i < n; ++i)
        dest[i] = a * x[i] + y[i];
}

void sub_vector(double *dest, double *x, double *y, int n) {

	for(int i = 0; i< n; i++) {
		dest[i] = x[i] - y[i];
	}
}

void add_vector(double *dest, double *x, double *y, int n) {

	for(int i = 0; i< n; i++) {
		dest[i] = x[i] + y[i];
	}
}

double ddot(double *x, double *y, int n) {
    double final_sum = 0;

    for (int i = 0; i < n; ++i)
        final_sum += x[i] * y[i];

    return final_sum;
}

double norm2(double *x, int n) {
	return sqrt(ddot(x, x, n));
}

void matvec(double *Ax, double **Adata, double *xvect, int n) {
    for (int i = 0; i < n; ++i) {
        Ax[i] = 0;
        for (int j = 0; j < n; ++j) {
            Ax[i] += Adata[i][j]*xvect[j];
        }
    }
}

/* Creation of a Jacobi preconditioner - Zeros are not stored * 
 * Mij = Aij, i=j                                             *
 * Mij = 0, i<>j                                              */                           
void jacobi_precond(double *M, double **Adata, int n) {
    for (int i = 0; i < n; ++i) {
    	if(Adata[i][i] == 0) {
    		M[i] = 1;
    	}
    	else {
			M[i] = 1/Adata[i][i];
		}
    }
}

/* Apply Jacobi preconditioner */
void prec_solve(double *Minvx, double *Mdata, double *x, int n) {
    for (int i = 0; i < n; i++) {
		Minvx[i] = Mdata[i]*x[i];
    }
}

/* Solve the system with conjugate gradient method and return the number of iterations */
int cg(double *x, double itol, int dimension, int maxiter, double *b, double **A) {

	/* Residual vector r */
	double *r = (double *)malloc(dimension * sizeof(double));
	/* Preconditioner M */
	double *M = (double *)malloc(dimension * sizeof(double));
	double *Ax = (double *)malloc(dimension * sizeof(double));
	/* Alocate z vector: solution of preconditioner */
	double *z = (double *)malloc(dimension * sizeof(double));
	double *p = (double *)malloc(dimension * sizeof(double));
	double *q = (double *)malloc(dimension * sizeof(double));
	double r_norm, b_norm, rho, rho1, alpha, beta;
	int iter = 0;

	/* Initialize M preconditioner */
	jacobi_precond(M, A, dimension);
	/* Initialize vector r with vector b */
	memcpy(r, b, dimension*sizeof(double));
	/* Compute A*x and store it to Ax */
	matvec(Ax, A, x, dimension);
	/* Compute r = b - Ax */
	sub_vector(r, r, Ax, dimension);

	/* Initialize norm2 of b and r */
	r_norm = norm2(r, dimension);
	b_norm = norm2(b, dimension);
	/* Set b_norm = 1 in case it's zero to avoid seg fault */
	b_norm = b_norm == 0 ? 1 : b_norm;
	while (iter < maxiter && (r_norm / b_norm) > itol) {
		iter++;
		/* Solution of the preconditioner */
		prec_solve(z, M, r, dimension);
		rho = ddot(r, z, dimension);
		if(iter == 1) {
			p = z;
		}
		else {
			beta = rho / rho1;
			axpy(p, beta, p, z, dimension);
		}
		rho1 = rho;
		matvec(q, A, p, dimension);
		alpha = rho / ddot(p, q, dimension);
		axpy(x, alpha, p, x, dimension);
		axpy(r, -alpha, q, r, dimension);
		r_norm = norm2(r, dimension);
	}
	return iter;
}

void bi_cg() {

}