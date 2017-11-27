#include <math.h>

#include "routines.h"

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
        Ax[i] = 0.0;
        for (int j = 0; j < n; j++) {
            Ax[i] += A[i][j] * x[j];
        }
    }
}

void mat_vec_mul_trans(double *Ax, double **A, double *x, int n) {
    for (int i = 0; i < n; i++) {
        Ax[i] = 0.0;
        for (int j = 0; j < n; j++) {
            Ax[i] += A[j][i] * x[j];
        }
    }
}

/* Multiplies matrix A with vector x and stores the reuslt in supplied vector Ax */
void cs_mat_vec_mul(double *Ax, cs *A, double *x) {
    zero_out_vec(Ax, A->n);
    for (int j = 0; j < A->n; j++) {
        for (int p = A->p[j]; p < A->p[j+1]; p++) {
            Ax[A->i[p]] += A->x[p] * x[j];
        }
    }
}

/* Multiplies matrix A with vector x and stores the reuslt in supplied vector Ax */
void cs_mat_vec_mul_trans(double *Ax, cs *A, double *x) {
    zero_out_vec(Ax, A->n);
    for (int j = 0; j < A->n; j++) {
        for (int p = A->p[j]; p < A->p[j+1]; p++) {
            Ax[j] += A->x[p] * x[A->i[p]];
        }
    }
}

/* Creation of a Jacobi preconditioner and stores it in supplied vector M, zeros are not stored */
void jacobi_precond(double *M, double **A, cs *C, int n, bool SPARSE) {
    if (!SPARSE) {
        for (int i = 0; i < n; i++) {
            /* We don't want to add zeros, we replace them with 1 instead */
            //TODO perhaps instead of 0 we check with EPSILON? CAUSE DOUBLE??
            if (A[i][i] == 0) {
                M[i] = 1.0;
            }
            else {
                M[i] = 1 / A[i][i];
            }
        }
    }
    else {
        int num_nodes;
        for (int j = 0; j < C->n; j++) {
            for (int p = C->p[j]; p < C->p[j+1]; p++) {
                if (C->i[p] == j) {
                    M[j] = 1 / C->x[p];
                    /* Save the last position of the non zero diagonal element */
                    num_nodes = j;
                }
            }
        }
        /* Handle the g2 elements which have zero on diagonal of A */
        for (int i = num_nodes + 1; i < C->n; i++) {
            M[i] = 1.0;
        }
    }
}

/* Apply Jacobi preconditioner and store it in vector M_fin */
void precond_solve(double *M_fin, double *M, double *x, int n) {
    for (int i = 0; i < n; i++) {
		M_fin[i] = M[i] * x[i];
    }
}

/* Zero outs the supplied vector x */
void zero_out_vec(double *x, int dimension) {
    for (int i = 0; i < dimension; i++) {
        x[i] = 0.0;
    }
}
