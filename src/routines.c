#include <math.h>

#include "routines.h"

/* Computes the dest = a*x + y , x and y are vectors and a is a constant */
void axpy(double *dest, double a, double *x, double *y, int n) {
	for (int i = 0; i < n; i++)
		dest[i] = a * x[i] + y[i];
}

/* Computes the dest = a*x + y , x and y are vectors and a is a constant */
void complex_axpy(gsl_vector_complex *dest, gsl_complex a, gsl_vector_complex *x, gsl_vector_complex *y, int n) {
	gsl_complex xi, yi, ax, axy;
	for (int i = 0; i < n; i++) {
		/* Get the current complex numbers from the vectors */
		xi = gsl_vector_complex_get(x, i);
		yi = gsl_vector_complex_get(y, i);
		/* Multiply the first with alpha */
		ax = gsl_complex_mul(a, xi);
		/* Add them together */
		axy = gsl_complex_add(ax, yi);
		/* Save it to dest */
		gsl_vector_complex_set(dest, i, axy);
	}
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
	double sum = 0.0;
	for (int i = 0; i < n; i++) {
		sum += x[i] * y[i];
	}
	return sum;
}

/* Computes the complex dot product of complex vectors x and y */
gsl_complex complex_dot_product(gsl_vector_complex *x, gsl_vector_complex *y, int n) {
	gsl_complex z, xy, xi, yi;
	/* Initialize to z = 0.0 + i0.0 */
	GSL_SET_COMPLEX(&z, 0.0, 0.0);
	for (int i = 0; i < n; i++) {
		/* Get the current complex number from the vectors */
		xi = gsl_vector_complex_get(x, i);
		yi = gsl_vector_complex_get(y, i);
		/* Get the conjugate of the first */
		xi = gsl_complex_conjugate(xi);
		/* Get their multiplication */
		xy = gsl_complex_mul(xi, yi);
		/* Get the sum */
		z  = gsl_complex_add(z, xy);
	}
	return z;
}

/* Computes the euclidean norm of vector x */
double norm2(double *x, int n) {
	return sqrt(dot_product(x, x, n));
}

/* Computes the euclidean norm of complex vector x */
double complex_norm2(gsl_vector_complex *x, int n) {
	/* 1. ||x|| = sqrt(Σ(xr^2 + xi^2)) for i...n */
	/* 2. ||x|| = sqrt(Σ(z*z_conj)) for i...n */
	/* Method 1 is used below */
	gsl_complex z;
	double z_real, z_imag, sum = 0.0;
	for (int i = 0; i < n; i++) {
		/* Get current complex number from vector */
		z = gsl_vector_complex_get(x, i);
		z_real = GSL_REAL(z) * GSL_REAL(z);
		z_imag = GSL_IMAG(z) * GSL_IMAG(z);
		sum += z_real + z_imag;
	}
	return sqrt(sum);
}

/* Returns the absolute value of the complex vector */
double complex_abs(gsl_complex z) {
	return sqrt(GSL_REAL(z) * GSL_REAL(z) + GSL_IMAG(z) * GSL_IMAG(z));
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
			if (A[i][i] == 0.0) {
				M[i] = 1.0;
			}
			else {
				M[i] = 1 / A[i][i];
			}
		}
	}
	else {
		for (int j = 0; j < C->n; j++) {
			for (int p = C->p[j]; p < C->p[j+1]; p++) {
				if (C->i[p] == j) {
					M[j] =  1 / C->x[p];
				}
			}
		}
	}
}

/* Creation of a complex Jacobi preconditioner and stores it in supplied vector M, zeros are not stored */
void complex_jacobi_precond(gsl_vector_complex *M, gsl_matrix_complex *A, int n) {
	// TODO add sparse struct and flag
	for (int i = 0; i < n; i++) {
		/* We don't want to add zeros, we replace them with 1 + 0i instead */
		if (IS_COMPLEX_ZERO(gsl_matrix_complex_get(A, i, i))) {
			gsl_vector_complex_set(M, i, GSL_COMPLEX_ONE);
		}
		else {
			/* We set M to 1/A[i][i] */
			gsl_vector_complex_set(M, i, gsl_complex_div(GSL_COMPLEX_ONE, gsl_matrix_complex_get(A, i, i)));
		}
	}
}

/* Apply Jacobi preconditioner and store it in vector M_fin */
void precond_solve(double *M_fin, double *M, double *x, int n) {
	for (int i = 0; i < n; i++) {
		M_fin[i] = M[i] * x[i];
	}
}

/* Apply complex Jacobi preconditioner and store it in vector M_fin */
void complex_precond_solve(gsl_vector_complex *M_fin, gsl_vector_complex *M, gsl_vector_complex *x, int n) {
	gsl_complex z, mi, xi;
	for (int i = 0; i < n; i++) {
		/* Get the current complex number from the vectors */
		mi = gsl_vector_complex_get(M, i);
		xi = gsl_vector_complex_get(x, i);
		/* Get their multiplication */
		z = gsl_complex_mul(mi, xi);
		/* Set the result to vector M_fin at i */
		gsl_vector_complex_set(M_fin, i, z);
	}
}

/* Zero outs the supplied vector x */
void zero_out_vec(double *x, int dimension) {
	for (int i = 0; i < dimension; i++) {
		x[i] = 0.0;
	}
}

/* Sets the value <val> to the supplied vector x */
void set_vec_val(double *x, double val, int dimension) {
	for (int i = 0; i < dimension; i++) {
		x[i] = val;
	}
}

/* Converts the complex z into polar form */
ac_t rect_to_polar(gsl_complex z) {
	ac_t ac;
	double real = GSL_REAL(z);
	double imag = GSL_IMAG(z);
	ac.magnitude = sqrt(real * real + imag * imag);
	ac.phase = to_degrees(atan2(imag, real));
	return ac;
}

/* Converts radians to degrees in the range of 0-360 */
double to_degrees(double radians) {
	double degrees = radians * (180.0 / M_PI);
	degrees = fmod(degrees + 360.0, 360.0);
	return degrees;
}

/* Converts a real vector x to a complex one with 0 imaginary parts */
void real_to_complex_vector(gsl_vector_complex *x_complex, double *x, int dimension) {
	/* Set a zero complex number 0 + 0i */
	gsl_complex z = GSL_COMPLEX_ZERO;
	for (int i = 0; i < dimension; i++) {
		/* Set the real part of z to x[i] */
		GSL_SET_REAL(&z, x[i]);
		/* Set z to v[i] */
		gsl_vector_complex_set(x_complex, i, z);
	}
}

/* Computes the conjugate of a vector M and stores it to vector M_conj */
void vector_conjugate(gsl_vector_complex *M_conj, gsl_vector_complex *M, int dimension) {
	gsl_complex z;
	for (int i = 0; i < dimension; i++) {
		/* Get the current complex number from the vector */
		z = gsl_vector_complex_get(M, i);
		/* Save its conjugate to vector M_conj */
		gsl_vector_complex_set(M_conj, i, gsl_complex_conjugate(z));
	}
}

/* Allocate memory for the matrix of the MNA system */
double **init_array(int row, int col) {
	double **array = (double **)malloc(row * sizeof(double *));
	assert(array != NULL);
	array[0] = (double *)calloc(row * col, sizeof(double));
	assert(array[0] != NULL);
	for (int i = 1; i < row; i++) {
		array[i] = array[i-1] + col;
	}
	return array;
}

/* Allocate memory for the vector of the MNA system */
double *init_vector(int row) {
	double *vector = (double *)calloc(row, sizeof(double));;
	assert(vector != NULL);
	return vector;
}

/* Allocate memory for the complex gsl matrix of the MNA system */
gsl_matrix_complex *init_gsl_complex_array(int row, int col) {
	gsl_matrix_complex *matrix = gsl_matrix_complex_alloc(row, col);
	assert(matrix != NULL);
	return matrix;
}

/* Allocate memory for the complex gsl vector of the MNA system */
gsl_vector_complex *init_gsl_complex_vector(int row) {
	gsl_vector_complex *vector = gsl_vector_complex_alloc(row);
	assert(vector != NULL);
	return vector;
}

/* Allocate memory for the new vector and set all values to val */
double *init_val_vector(int row, double val) {
	double *vector = (double *)malloc(row * sizeof(double));
	assert(vector != NULL);
	for (int i = 0; i < row ; i++) {
		vector[i] = val;
	}
	return vector;
}

/* Allocate memory for the GSL permutation matrix of the MNA system */
gsl_permutation *init_permutation(int dimension) {
	gsl_permutation *P = gsl_permutation_calloc(dimension);
	return P;
}