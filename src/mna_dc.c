#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "mna_dc.h"
#include "iter.h"

/* Allocate memory for the MNA system */
mna_system_t *init_mna_system(int num_nodes, int num_g2_elem) {
	int dimension = num_nodes + num_g2_elem;
	mna_system_t *mna = (mna_system_t *)malloc(sizeof(mna_system_t));
	assert(mna != NULL);
	mna->A = init_array(dimension, dimension);
	mna->b = init_vector(dimension);
	mna->P = init_permutation(dimension);
	mna->is_decomp = false;
	mna->num_nodes = num_nodes;
	mna->num_g2_elem = num_g2_elem;
	mna->g2_indx = (g2_indx_t *)malloc(num_g2_elem * sizeof(g2_indx_t));
	assert(mna->g2_indx != NULL);
	return mna;
}

/* Allocate memory for the matrix of the MNA system */
double **init_array(int row, int col) {
	double **array = (double **)malloc(row * sizeof(double *));
	assert(array !=NULL);
	array[0] = (double *)calloc(row * col, sizeof(double));
	assert(array[0] !=NULL);
	for (int i = 1; i < row; i++) {
		array[i] = array[i-1] + col;
	}
	return array;
}

/* Allocate memory for the vector of the MNA system */
double *init_vector(int row) {
	double *vector = (double *) calloc(row, sizeof(double));;
	assert(vector != NULL);
	return vector;
}

/* Allocate memory for the GSL permutation matrix of the MNA system */
gsl_permutation *init_permutation(int dimension) {
	gsl_permutation *P;
	P = gsl_permutation_calloc(dimension);
	return P;
}

/* Create the MNA system array and vector */
void create_mna_system(mna_system_t *mna, index_t *index, hash_table_t *hash_table, int offset) {
	list1_t *curr;
	double value;
	int volt_sources_cnt = 0;
	for (curr = index->head1; curr != NULL; curr = curr->next) {
		int probe1_id = ht_get_id(hash_table, curr->probe1);
		int probe2_id = ht_get_id(hash_table, curr->probe2);
		int i = probe1_id - 1;
		int j = probe2_id - 1;
		if (curr->type == 'C' || curr->type == 'c') {
			continue;
		}
		else if (curr->type == 'R' || curr->type == 'r') {
			value = 1 / curr->value;
			if (probe1_id == 0) {
				mna->A[j][j] += value; 
			}
			else if (probe2_id == 0) {
				mna->A[i][i] += value; 
			}
			else {
				mna->A[i][i] += value;
				mna->A[j][j] += value;
				mna->A[i][j] -= value;
				mna->A[j][i] -= value;
			}
		}
		else if (curr->type == 'I' || curr->type == 'i') {
			value = curr->value;
			if (probe1_id == 0) {
				mna->b[j] += value;
			}
			else if (probe2_id == 0) {
				mna->b[i] -= value;
			}
			else {
				mna->b[i] -= value;
				mna->b[j] += value;
			}
		}
		else {
			if (curr->type == 'L' || curr->type == 'l') {
				/* We treat the Inductor like a voltage source with value 0 */
				value = 0.0;
			}
			else if (curr->type == 'V' || curr->type == 'v') {
				value = curr->value;
			}
			/* Save the g2 element source you find, keep indexing */
			mna->g2_indx[volt_sources_cnt].element = (char *)malloc(strlen(curr->element) * sizeof(char));
			assert(mna->g2_indx[volt_sources_cnt].element != NULL);
			strcpy(mna->g2_indx[volt_sources_cnt].element, curr->element);
			if (probe1_id == 0) {
				mna->A[j][offset + volt_sources_cnt] = -1.0;
				mna->A[offset + volt_sources_cnt][j] = -1.0;
				mna->b[offset + volt_sources_cnt] += value; 
			}
			else if (probe2_id == 0) {
				mna->A[i][offset + volt_sources_cnt] = 1.0;
				mna->A[offset + volt_sources_cnt][i] = 1.0;
				mna->b[offset + volt_sources_cnt] += value;
			}
			else {
				mna->A[i][offset + volt_sources_cnt] =  1.0;
				mna->A[j][offset + volt_sources_cnt] = -1.0;
				mna->A[offset + volt_sources_cnt][i] =  1.0;
				mna->A[offset + volt_sources_cnt][j] = -1.0;
				mna->b[offset + volt_sources_cnt] += value;
			}
			/* Keep track of how many voltage sources or inductors (which are treated like voltages with 0), we have already found */
			volt_sources_cnt++;
		}
	}
}

/* Searches through the g2_indx of the mna system and returns the index of the argument element */
int g2_elem_indx(g2_indx_t *g2_indx, int num_nodes, int num_g2_elem, char *element) {
	for (int i = 0; i < num_g2_elem; i++) {
		if (strcmp(element, g2_indx[i].element) == 0) {
			/* Return the offset index */
			return num_nodes + i;
		}
	}
	return 0;
}

/* Solves the mna system according to the specified method provided by options argument
 * (LU, Cholesky, Iterative Conj_Grad / Bi-Conj_Grad) and stores the solution on the supplied vector x
 */
void solve_mna_system(mna_system_t *mna, double **x, options_t *options) {
	int dimension = mna->num_nodes + mna->num_g2_elem;
	gsl_vector_view view_x = gsl_vector_view_array(*x, dimension);
	if (options->ITER) {
		if (options->SPD) {
			int iterations = conj_grad(mna->A, *x, mna->b, dimension, options->itol, dimension);
			printf("Conjugate gradient method did %d iterations.\n", iterations);
		}
		else {
			// int iterations = bi_conj_grad(mna->A, *x, mna->b, dimension, options->itol, dimension);
		}
	}
	else {
		if (options->SPD) {
			solve_cholesky(mna, view_x);
		}
		else {
			solve_lu(mna, view_x);
		}
	}
}

/* Solve the MNA system using LU decomposition */
void solve_lu(mna_system_t *mna, gsl_vector_view x) {
	/* The sign of the permutation matrix */
	int signum;
	/* Allocate memory for the solution vector */
	int dimension = mna->num_nodes + mna->num_g2_elem;
	gsl_matrix_view view_A = gsl_matrix_view_array(mna->A[0], dimension, dimension);
	gsl_vector_view view_b = gsl_vector_view_array(mna->b, dimension);
	if (!mna->is_decomp) {
		/* LU decomposition on A, PA = LU */
		gsl_linalg_LU_decomp(&view_A.matrix, mna->P, &signum);
		mna->is_decomp = true;
		printf("LU Matrix:\n\n");
		print_array(mna->A, dimension);
		printf("Permutation Vector:\n\n");
		print_permutation(mna->P);
	}
	/* Solve the LU system */
	gsl_linalg_LU_solve(&view_A.matrix, mna->P, &view_b.vector, &x.vector);
}

/* Solve the MNA system using cholesky decomposition */
void solve_cholesky(mna_system_t *mna, gsl_vector_view x) {
	int dimension = mna->num_nodes + mna->num_g2_elem;
	/* Allocate memory for the solution vector */
	gsl_matrix_view view_A = gsl_matrix_view_array(mna->A[0], dimension, dimension);
	gsl_vector_view view_b = gsl_vector_view_array(mna->b, dimension);
	if (!mna->is_decomp) {
		/* Cholesky decomposition A = LL^T*/
		gsl_linalg_cholesky_decomp(&view_A.matrix);
		mna->is_decomp = true;
		printf("Cholesky Matrix:\n\n");
		print_array(mna->A, dimension);
		printf("\n\n");
	}
	/* Solve the cholesky system */
	gsl_linalg_cholesky_solve(&view_A.matrix, &view_b.vector, &x.vector);
}

/* Print the MNA system */
void print_mna_system(mna_system_t *mna) {
	int dimension = mna->num_nodes + mna->num_g2_elem;
	printf("MNA A array:\n\n");
	print_array(mna->A, dimension);
	printf("MNA b vector:\n\n");
	print_vector(mna->b, dimension);
}

/* Print the array */
void print_array(double **A, int dimension) {
	int rows = dimension;
	int cols = dimension;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			printf("%-15.4lf", A[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

/* Print the vector */
void print_vector(double *b, int dimension) {
	for (int i = 0; i< dimension; i++) {
		printf("%lf\n", b[i]);
	}
    printf("\n");
}

/* Print the permutation */
void print_permutation(gsl_permutation *P) {
	int dimension = P->size;
	for (int i = 0; i < dimension; i++) {
		printf("%zu ", gsl_permutation_get(P, i));
	}
	printf("\n\n");
}

/* Free all the memory allocated for the MNA system */
void free_mna_system(mna_system_t **mna) {
	/* Free everything we allocated for the MNA and GSL */
	free((*mna)->A[0]);
	free((*mna)->A);
	free((*mna)->b);
	gsl_permutation_free((*mna)->P);
	/* Free every string allocated for the group2 elements */
	for (int i = 0; i < (*mna)->num_g2_elem; i++) {
		free((*mna)->g2_indx[i].element);
	}
	/* Free the array */
	free((*mna)->g2_indx);
	free(*mna);
	/* Set mna to NULL to limit further acesses */
	*mna = NULL;
}
