#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "mna_dc.h"

/* Allocate memory for the MNA system */
mna_system_t *init_mna_system(int num_nodes, int num_g2_elem) {
	int dimension = num_nodes + num_g2_elem;
	mna_system_t *mna = (mna_system_t *)malloc(sizeof(mna_system_t));
	assert(mna != NULL);
	mna->A = init_array(dimension, dimension);
	mna->b = init_vector(dimension);
	mna->P = init_permutation(dimension);
	mna->is_decomp = 0;
	mna->num_nodes = num_nodes;
	mna->num_g2_elem = num_g2_elem;
	mna->g2_indx = (g2_indx_t *)malloc(num_g2_elem * sizeof(g2_indx_t));
	assert(mna->g2_indx != NULL);
	return mna;
}

/* Allocate memory for the GSL matrix of the MNA system */
gsl_matrix *init_array(int row, int col) {
    gsl_matrix *array;
    array = gsl_matrix_calloc(row, col);
    return array;
}

/* Allocate memory for the GSL vector of the MNA system */
gsl_vector *init_vector(int row) {
	gsl_vector *vector;
	vector = gsl_vector_calloc(row);
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
				gsl_matrix_set(mna->A, j, j, gsl_matrix_get(mna->A, j, j) + value);
			}
			else if (probe2_id == 0) {
				gsl_matrix_set(mna->A, i, i, gsl_matrix_get(mna->A, i, i) + value);
			}
			else {
				gsl_matrix_set(mna->A, i, i, gsl_matrix_get(mna->A, i, i) + value);
				gsl_matrix_set(mna->A, j, j, gsl_matrix_get(mna->A, j, j) + value);
				gsl_matrix_set(mna->A, i, j, gsl_matrix_get(mna->A, i, j) - value);
				gsl_matrix_set(mna->A, j, i, gsl_matrix_get(mna->A, j, i) - value);
			}
		}
		else if (curr->type == 'I' || curr->type == 'i') {
			value = curr->value;
			if (probe1_id == 0) {
				gsl_vector_set(mna->b, j, gsl_vector_get(mna->b, j) + value);
			}
			else if (probe2_id == 0) {
				gsl_vector_set(mna->b, i, gsl_vector_get(mna->b, i) + value);
			}
			else {
				gsl_vector_set(mna->b, i, gsl_vector_get(mna->b, i) + value);
				gsl_vector_set(mna->b, j, gsl_vector_get(mna->b, j) + value);
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
				gsl_matrix_set(mna->A, j, offset + volt_sources_cnt, 1.0);
				gsl_matrix_set(mna->A, offset + volt_sources_cnt, j, 1.0);
				gsl_vector_set(mna->b, offset + volt_sources_cnt, gsl_vector_get(mna->b, offset + volt_sources_cnt) + value);
			}
			else if (probe2_id == 0) {
				gsl_matrix_set(mna->A, i, offset + volt_sources_cnt, 1.0);
				gsl_matrix_set(mna->A, offset + volt_sources_cnt, i, 1.0);
				gsl_vector_set(mna->b, offset + volt_sources_cnt, gsl_vector_get(mna->b, offset + volt_sources_cnt) + value);
			}
			else {
				gsl_matrix_set(mna->A, i, offset + volt_sources_cnt, 1.0);
				gsl_matrix_set(mna->A, j, offset + volt_sources_cnt, 1.0);
				gsl_matrix_set(mna->A, offset + volt_sources_cnt, i, 1.0);
				gsl_matrix_set(mna->A, offset + volt_sources_cnt, j, 1.0);
				gsl_vector_set(mna->b, offset + volt_sources_cnt, gsl_vector_get(mna->b, offset + volt_sources_cnt) + value);
			}
			/* Keep track of how many voltage sources or inductors (which are treated like voltages with 0), we have already found */
			volt_sources_cnt++;
		}
	}
}

/* Searches through the g2_indx of the mna system and returns the index of the argument element */
int g2_elem_indx(g2_indx_t *g2_indx, int num_nodes, int num_g2_elem, char *element) {
	for (int i = 0; i < num_g2_elem; i++) {
		if(strcmp(element, g2_indx[i].element) == 0) {
			/* Return the offset index */
			return num_nodes + i;
		}
	}
	return FAILURE;
}

/* LU or Cholesky decomposition and solution of the MNA system Ax=b and returns the solution vector x */
gsl_vector *solve_mna_system(mna_system_t *mna, bool SPD) {
	if (SPD) {
		return solve_cholesky(mna);
	}
	return solve_lu(mna);
}

/* Solve the MNA system using LU decomposition */
gsl_vector *solve_lu(mna_system_t *mna) {
	/* The sign of the permutation matrix */
	int signum;
	int dimension = mna->num_nodes + mna->num_g2_elem;
	/* Allocate memory for the solution vector */
	gsl_vector *x = gsl_vector_calloc(dimension);
	if (!mna->is_decomp) {
		/* LU decomposition on A, PA = LU */
		gsl_linalg_LU_decomp(mna->A, mna->P, &signum);
		mna->is_decomp = true;
	}
	/* Solve the LU system */
	gsl_linalg_LU_solve(mna->A, mna->P, mna->b, x);
	return x;
}

/* Solve the MNA system using cholesky decomposition */
gsl_vector *solve_cholesky(mna_system_t *mna) {
	int dimension = mna->num_nodes + mna->num_g2_elem;
	/* Allocate memory for the solution vector */
	gsl_vector *x = gsl_vector_calloc(dimension);
	if (!mna->is_decomp) {
		/* Cholesky decomposition A = LL^T*/
		gsl_linalg_cholesky_decomp(mna->A);
		mna->is_decomp = true;
	}
	/* Solve the cholesky system */
	gsl_linalg_cholesky_solve(mna->A, mna->b, x);
	return x;
}

/* Print the MNA system */
void print_mna_system(mna_system_t *mna) {
	printf("MNA A array:\n\n");
	print_array(mna->A);
	printf("MNA b vector:\n\n");
	print_vector(mna->b);
}

/* Print the array */
void print_array(gsl_matrix *A) {
	gsl_matrix_fprintf(stdout, A, "%lf");
	printf("\n");
}

/* Print the vector */
void print_vector(gsl_vector *b) {
	gsl_vector_fprintf(stdout, b, "%lf");
    printf("\n");
}

/* Free all the memory allocated for the MNA system */
void free_mna_system(mna_system_t **mna) {
	/* Free everything we allocated from GSL */
	gsl_matrix_free((*mna)->A);
	gsl_vector_free((*mna)->b);
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
