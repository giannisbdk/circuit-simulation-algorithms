#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "iter.h"

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

/* Computes the dest = a*x + y , x and y are vectors and a is a constant */
void axpy(double *dest, double a, double *x, double *y, int n) {
    for (int i = 0; i < n; i++)
        dest[i] = a * x[i] + y[i];
}

/* Substracts two vectors and stores the result in dest */
void sub_vector(double *dest, double *x, double *y, int n) {
	for (int i = 0; i < n; i++) {
		dest[i] = x[i] - y[i];
	}
}

/* Adds two vectors and stores the result in dest */
void add_vector(double *dest, double *x, double *y, int n) {
	for (int i = 0; i < n; i++) {
		dest[i] = x[i] + y[i];
	}
}

/* Computes the dot product of vectors x and y */
double dot_product(double *x, double *y, int n) {
    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += x[i] * y[i];
    }
    return sum;
}

/* Computes the euclidean norm of vector x */
double norm2(double *x, int n) {
	return sqrt(dot_product(x, x, n));
}

/* Multiplies matrix A with vector x and stores the reuslt in supplied vector Ax */
void mat_vec_mul(double *Ax, double **A, double *x, int n) {
    for (int i = 0; i < n; i++) {
        Ax[i] = 0;
        for (int j = 0; j < n; j++) {
            Ax[i] += A[i][j] * x[j];
        }
    }
}

/* Creation of a Jacobi preconditioner and stores it in supplied vector M, zeros are not stored */
void jacobi_precond(double *M, double **A, int n) {
    for (int i = 0; i < n; i++) {
    	/* We don't want to add zeros, we replace them with 1 instead */
    	//TODO perhaps instead of 0 we check with EPSILON? CAUSE DOUBLE??
    	if(A[i][i] == 0) {
    		M[i] = 1;
    	}
    	else {
			M[i] = 1/A[i][i];
		}
    }
}

/* Apply Jacobi preconditioner and store it in vector M_fin */
void precond_solve(double *M_fin, double *M, double *x, int n) {
    for (int i = 0; i < n; i++) {
		M_fin[i] = M[i] * x[i];
    }
}

/* Transposes the given matrix A to the destination matrix A_trans */
void trans_matrix(double *A_trans, double **A, int n) {
	//TODO This has to be done with an efficient way
}

/* Solve the SPD system with the iterative conjugate gradient method
 * stores the result in vector x and also returns the number of iterations
 */ 
int conj_grad(double **A, double *x, double *b, int dimension, double itol, int maxiter) {
	double *Ax = (double *)malloc(dimension * sizeof(double));
	/* Preconditioner M */
	double *M  = (double *)malloc(dimension * sizeof(double));
	/* Residual vector r */
	double *r  = (double *)malloc(dimension * sizeof(double));
	/* Alocate z vector: solution of preconditioner */
	double *z  = (double *)malloc(dimension * sizeof(double));
	double *p  = (double *)malloc(dimension * sizeof(double));
	double *q  = (double *)malloc(dimension * sizeof(double));
	double r_norm, b_norm, rho, rho1, alpha, beta;
	int iter = 0;

	/* Initialize M preconditioner */
	jacobi_precond(M, A, dimension);
	/* Compute A*x and store it to Ax */
	mat_vec_mul(Ax, A, x, dimension);
	/* Compute r = b - Ax */
	sub_vector(r, b, Ax, dimension);

	/* Initialize norm2 of b and r vectors */
	r_norm = norm2(r, dimension);
	b_norm = norm2(b, dimension);
	/* Set b_norm = 1 in case it's zero to avoid seg fault */
	b_norm = b_norm == 0 ? 1 : b_norm;

	while (iter < maxiter && (r_norm / b_norm) > itol) {
		iter++;
		/* Solution of the preconditioner Mz = r */
		precond_solve(z, M, r, dimension);
		/* rho = r*z */
		rho = dot_product(r, z, dimension);
		if (iter == 1) {
			/* Set p = z */
			memcpy(p, z, dimension * sizeof(double));
		}
		else {
			beta = rho / rho1;
			/* p = z + beta*p */
			axpy(p, beta, p, z, dimension);
		}
		rho1 = rho;
		/* q = A*p */
		mat_vec_mul(q, A, p, dimension);
		/* a = rho / p*q */
		alpha = rho / dot_product(p, q, dimension);
		/* x = x + alpha*p */
		axpy(x,  alpha, p, x, dimension);
		/* r = r - alpha*q */
		axpy(r, -alpha, q, r, dimension);
		r_norm = norm2(r, dimension);
	}
	/* Free all the memory we allocated */
	free(Ax);
	free(M);
	free(r);
	free(z);
	free(p);
	free(q);
	return iter;
}

/* Solve the system with bi-conjugate gradient method and return the number of iterations */
int bi_conj_grad(double **A, double *x, double *b, int dimension, double itol, int maxiter) {
	/* Preconditioner M */
	double *M       = (double *)malloc(dimension * sizeof(double));
	/* Transpose of Preconditioner M */
	double *M_trans = (double *)malloc(dimension * sizeof(double));
	/* Transpose of Matrix A */
	double *A_trans = (double *)malloc(dimension * sizeof(double));
	/* Vector to store A*x */
	double *Ax = (double *)malloc(dimension * sizeof(double));
	/* Residual vector r */
	double *r  = (double *)malloc(dimension * sizeof(double));
	/* Alocate z vector: solution of preconditioner */
	double *z  = (double *)malloc(dimension * sizeof(double));
	double *p  = (double *)malloc(dimension * sizeof(double));
	double *q  = (double *)malloc(dimension * sizeof(double));
	/* Residual vector r_tilde */
	double *r_tilde = (double *)malloc(dimension * sizeof(double));
	/* Alocate z_tilde vector: solution of preconditioner */
	double *z_tilde = (double *)malloc(dimension * sizeof(double));
	double *p_tilde = (double *)malloc(dimension * sizeof(double));
	double *q_tilde = (double *)malloc(dimension * sizeof(double));
	double r_norm, b_norm, rho, rho1, alpha, beta, omega;
	int iter = 0;

	/* Initialize M preconditioner */
	jacobi_precond(M, A, dimension);
	/* Copy M preconditioner to M_trans they are both equal */
	memcpy(M_trans, M, dimension * sizeof(double));
	/* Compute A*x and store it to Ax */
	mat_vec_mul(Ax, A, x, dimension);
	/* Compute r = b - Ax */
	sub_vector(r, b, Ax, dimension);
	/* Compute r_tilde = b - Ax = r,               *
	 * sub_vector(r_tilde, b, Ax, dimension);      */
	memcpy(r_tilde, r, dimension * sizeof(double));
	/* Get the transpose matrix of A */
	trans_matrix(A_trans, A, dimension);

	/* Initialize norm2 of b and r vectors */
	r_norm = norm2(r, dimension);
	b_norm = norm2(b, dimension);
	/* Set b_norm = 1 in case it's zero to avoid seg fault */
	b_norm = b_norm == 0 ? 1 : b_norm;

	while (iter < maxiter && (r_norm / b_norm) > itol) {
		iter++;
		/* Solution of the preconditioner Mz = r */
		precond_solve(z, M, r, dimension);
		/* Solution of the preconditioner M'z_tilde = r_tilde, *
		 * M' = M because is a diagonal matrix                 */
		precond_solve(z_tilde, M, r_tilde, dimension);
		/* rho = r_tilde*z */
		rho = dot_product(r_tilde, z, dimension);
		/* Check for Algorithm Failure */
		if (fabs(rho) < EPSILON) {
			return -1;
		}
		if (iter == 1) {
			/* Set p = z */
			memcpy(p, z, dimension * sizeof(double));
			/* Set p_tilde = z_tilde */
			memcpy(p_tilde, z_tilde, dimension * sizeof(double));
		}
		else {
			beta = rho / rho1;
			/* p = z + beta*p */
			axpy(p, beta, p, z, dimension);
			/* p_tilde = z_tilde + beta*p_tilde */
			axpy(p_tilde, beta, p_tilde, z_tilde, dimension);
		}
		rho1 = rho;
		/* q = A*p */
		mat_vec_mul(q, A, p, dimension);
		/* q_tilde = A'*p_tilde, A' = A */
		//TODO Use A_Trans instead of A
		mat_vec_mul(q_tilde, A, p_tilde, dimension);
		/* omega = p_tilde * q */
		omega = dot_product(p_tilde, q, dimension);
		/* Check for Algorithm Failure */
		if (fabs(omega) < EPSILON) {
			return -1;
		}
		alpha = rho / omega;
		/* x = x + alpha*p */
		axpy(x,  alpha, p, x, dimension);
		/* r = r - alpha*q */
		axpy(r, -alpha, q, r, dimension);
		/* r_tilde = r_tilde - alpha*q_tilde */
		axpy(r_tilde, -alpha, q_tilde, r_tilde, dimension);
		r_norm = norm2(r, dimension);
	}
	/* Free all the memory we allocated */
	free(Ax);
	free(M);
	free(r);
	free(z);
	free(p);
	free(q);
	free(r_tilde);
	free(z_tilde);
	free(p_tilde);
	free(q_tilde);
	return iter;
}