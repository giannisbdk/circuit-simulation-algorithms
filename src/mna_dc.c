#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "mna_dc.h"

/* Allocate memory for the MNA system */
mna_system_t *init_mna_system(int num_nodes, int num_g2_elem, options_t *options, int nz) {
	/* Allocate for the whole struct */
	mna_system_t *mna = (mna_system_t *)malloc(sizeof(mna_system_t));
	assert(mna != NULL);
	mna->dimension = num_nodes + num_g2_elem;
	if (options->SPARSE) {
		init_sparse_matrix(mna, options, nz);
	}
	else {
		init_dense_matrix(mna, options);
	}
	/* In case we will use iterative methods allocate memory for the prerequisites */
	if (options->ITER) {
		mna->M = init_vector(mna->dimension);
		if(options->TRAN) {
			mna->M_trans = init_vector(mna->dimension);
		}
		else {
			mna->M_trans = NULL;
		}
	}
	else {
		mna->M = NULL;
		mna->M_trans = NULL;
	}
	mna->b = init_vector(mna->dimension);
	mna->is_decomp = false;
	mna->num_nodes = num_nodes;
	mna->num_g2_elem = num_g2_elem;
	mna->g2_indx = (g2_indx_t *)malloc(num_g2_elem * sizeof(g2_indx_t));
	assert(mna->g2_indx != NULL);

	mna->tr_analysis_init = false;

	return mna;
}

/* Allocate memory for a dense representation of MNA */
void init_dense_matrix(mna_system_t *mna, options_t *options) {
	/* Allocate for the matrices struct */
	mna->matrix = (matrix_t *)malloc(sizeof(matrix_t));
	assert(mna->matrix != NULL);
	/* Allocate for every different matrix/vector */
	mna->matrix->A = init_array(mna->dimension, mna->dimension);
	mna->matrix->P = init_permutation(mna->dimension);
	mna->sp_matrix = NULL;
	if (options->TRAN) {
		mna->matrix->G    = init_array(mna->dimension, mna->dimension);
		mna->matrix->hC   = init_array(mna->dimension, mna->dimension);
		mna->matrix->aGhC = init_array(mna->dimension, mna->dimension);
		/* If method is not the default Trapezoidal, we dont need the sGhC matrix, we use hC insteed */
		if (!options->TR) {
			mna->matrix->sGhC = NULL;
		}
		else {
			mna->matrix->sGhC = init_array(mna->dimension, mna->dimension);
		}
		mna->matrix->resp = (resp_t *)malloc(sizeof(resp_t));
		mna->matrix->resp->value = init_vector(mna->dimension);
		mna->matrix->resp->nodes = (list1_t **)malloc(mna->dimension * sizeof(list1_t *));
		assert(mna->matrix->resp->nodes != NULL);
		mna->matrix->resp->nodes[0] = (list1_t *)malloc(mna->dimension * sizeof(list1_t));
		assert(mna->matrix->resp->nodes[0] != NULL);
		for (int i = 1; i < mna->dimension; i++) {
			mna->matrix->resp->nodes[i] = mna->matrix->resp->nodes[i-1] + 1;
		}
	}
}

/* Allocate memory for a sparse representation of MNA */
void init_sparse_matrix(mna_system_t *mna, options_t *options, int nz) {
	mna->sp_matrix = (sp_matrix_t *)malloc(sizeof(sp_matrix_t));
	assert(mna->sp_matrix != NULL);
	mna->sp_matrix->A = cs_spalloc(mna->dimension, mna->dimension, nz, 1 , 1);
	assert(mna->sp_matrix->A != NULL);
	mna->sp_matrix->A->nz = 0;
	mna->matrix = NULL;

	if (options->TRAN) {
		mna->sp_matrix->G    = cs_spalloc(mna->dimension, mna->dimension, nz, 1 , 1);
		mna->sp_matrix->hC   = cs_spalloc(mna->dimension, mna->dimension, nz, 1 , 1);
		mna->sp_matrix->aGhC = cs_spalloc(mna->dimension, mna->dimension, nz, 1 , 1);
		/* If method is not the default Trapezoidal, we dont need the sGhC matrix, we use hC insteed */
		if (!options->TR) {
			mna->sp_matrix->sGhC = NULL;
		}
		else {
			mna->sp_matrix->sGhC = cs_spalloc(mna->dimension, mna->dimension, nz, 1 , 1);
		}
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

/* Allocate memory for the GSL permutation matrix of the MNA system */
gsl_permutation *init_permutation(int dimension) {
	gsl_permutation *P;
	P = gsl_permutation_calloc(dimension);
	return P;
}

/* Constructs the MNA system either sparse or dense */
void create_mna_system(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, double tr_step, int offset) {
	if (options->SPARSE) {
		create_sparse_mna(mna, index, hash_table, options, offset);
		if(options->TRAN) {
			create_sparse_trans_mna(mna, index, hash_table, options, tr_step, offset);
		}
	}
	else {
		create_dense_mna(mna, index, hash_table, options, offset);
		if(options->TRAN) {
			create_dense_trans_mna(mna, index, hash_table, options, tr_step, offset);
		}
	}
	printf("Creation of MNA system... OK\n");
}

/* Constructs the dense MNA system */
void create_dense_mna(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, int offset) {
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
				mna->matrix->A[j][j] += value; 
			}
			else if (probe2_id == 0) {
				mna->matrix->A[i][i] += value; 
			}
			else {
				mna->matrix->A[i][i] += value;
				mna->matrix->A[j][j] += value;
				mna->matrix->A[i][j] -= value;
				mna->matrix->A[j][i] -= value;
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
				mna->matrix->A[j][offset + volt_sources_cnt] = -1.0;
				mna->matrix->A[offset + volt_sources_cnt][j] = -1.0;
				mna->b[offset + volt_sources_cnt] += value; 
			}
			else if (probe2_id == 0) {
				mna->matrix->A[i][offset + volt_sources_cnt] = 1.0;
				mna->matrix->A[offset + volt_sources_cnt][i] = 1.0;
				mna->b[offset + volt_sources_cnt] += value;
			}
			else {
				mna->matrix->A[i][offset + volt_sources_cnt] =  1.0;
				mna->matrix->A[j][offset + volt_sources_cnt] = -1.0;
				mna->matrix->A[offset + volt_sources_cnt][i] =  1.0;
				mna->matrix->A[offset + volt_sources_cnt][j] = -1.0;
				mna->b[offset + volt_sources_cnt] += value;
			}
			/* Keep track of how many voltage sources or inductors (which are treated like voltages with 0), we have already found */
			volt_sources_cnt++;
		}
	}
	if (options->ITER) {
		/* Compute the M preconditioner */
		jacobi_precond(mna->M, mna->matrix->A, NULL, mna->dimension, options->SPARSE);
	}
}

void create_dense_trans_mna(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, double tr_step, int offset) {
	list1_t *curr;
	double value, h;
	int volt_sources_cnt = 0;

	memcpy(*mna->matrix->G, *mna->matrix->A, mna->dimension * mna->dimension * sizeof(double));
	assert(mna->matrix->G != NULL);

	if (!options->BE) {
		h = 2 / tr_step;
	}
	else {
		h = 1 / tr_step;
	}

	for (curr = index->head1; curr != NULL; curr = curr->next) {
		int probe1_id = ht_get_id(hash_table, curr->probe1);
		int probe2_id = ht_get_id(hash_table, curr->probe2);
		int i = probe1_id - 1;
		int j = probe2_id - 1;
		if (curr->type == 'C' || curr->type == 'c') {
			value = h * curr->value;
			if (probe1_id == 0) {
				mna->matrix->hC[j][j] += value;
			}
			else if (probe2_id == 0) {
				mna->matrix->hC[i][i] += value;
			}
			else {
				mna->matrix->hC[i][i] += value;
				mna->matrix->hC[j][j] += value;
				mna->matrix->hC[i][j] -= value;
				mna->matrix->hC[j][i] -= value;
			}
		}
		else if(curr->type == 'L' || curr->type == 'l') {
			/* Set the L value in the diagonal of g2 area in matrices */
			mna->matrix->hC[offset + volt_sources_cnt][offset + volt_sources_cnt] = -curr->value * h;
			/* Keep track of how many voltage sources or inductors, we have already found */
			volt_sources_cnt++;
		}
		else if (curr->type == 'V' || curr->type == 'v') {
			value = get_response_value(curr);
			/* Keep track of the node that contributes */
			mna->matrix->resp->nodes[offset + volt_sources_cnt] = curr;
			/* Keep track of how many voltage sources or inductors, we have already found */
			if (probe1_id == 0) {
				mna->matrix->resp->value[offset + volt_sources_cnt] += value;
			}
			else if (probe2_id == 0) {
				mna->matrix->resp->value[offset + volt_sources_cnt] += value;
			}
			else {
				mna->matrix->resp->value[offset + volt_sources_cnt] += value;
			}
			/* Keep track of how many voltage sources or inductors, we have already found */
			volt_sources_cnt++;
		}
		else if (curr->type == 'I' || curr->type == 'i') {
			value = get_response_value(curr);
			if (probe1_id == 0) {
				mna->matrix->resp->value[j] += value;
				mna->matrix->resp->nodes[j]  = curr;
			}
			else if (probe2_id == 0) {
				mna->matrix->resp->value[i] -= value;
				mna->matrix->resp->nodes[i]  = curr;
			}
			else {
				mna->matrix->resp->value[i] -= value;
				mna->matrix->resp->value[j] += value;
				mna->matrix->resp->nodes[i]  = curr;
				mna->matrix->resp->nodes[j]  = curr;
			}
		}
	}
	// printf("Printing resp struct\n");
	// for(int i = 0; i < mna->dimension; i++) {
	// 	if(mna->matrix->resp->nodes[i]->trans_spec != NULL) {
	// 		printf("curr element %s, type is %d, init value: %lf\n", mna->matrix->resp->nodes[i]->element, 
	// 			mna->matrix->resp->nodes[i]->trans_spec->type, mna->matrix->resp->value[i]);
	// 	}
	// 	else {
	// 		printf("curr element %s, no type, init value: %lf\n", mna->matrix->resp->nodes[i]->element, mna->matrix->resp->value[i]);
	// 	}
	// }

	/* Compute the matrices G + 1/h*C and G - 2/h*C */
	for (int i = 0; i < mna->dimension; i++) {
		for (int j = 0; j < mna->dimension; j++) {
			mna->matrix->aGhC[i][j] = mna->matrix->G[i][j] + mna->matrix->hC[i][j];
			if (options->TR) {
				mna->matrix->sGhC[i][j] = mna->matrix->G[i][j] - mna->matrix->hC[i][j];
			}
		}
	}

	if (options->ITER) {
		/* Compute the M preconditioner */
		jacobi_precond(mna->M_trans, mna->matrix->aGhC, NULL, mna->dimension, options->SPARSE);
	}

	/* Free what is no longer needed */
	// free(mna->matrix->G);
	// if (options->TR) {
	// 	free(mna->matrix->hC);
	// }
}

/* Constructs the sparse MNA system */
void create_sparse_mna(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, int offset) {
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
				cs_entry(mna->sp_matrix->A, j, j, value);
			}
			else if (probe2_id == 0) {
				cs_entry(mna->sp_matrix->A, i, i, value);
			}
			else {
				cs_entry(mna->sp_matrix->A, i, i, value);
				cs_entry(mna->sp_matrix->A, j, j, value);
				cs_entry(mna->sp_matrix->A, i, j, -value);
				cs_entry(mna->sp_matrix->A, j, i, -value);
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
				cs_entry(mna->sp_matrix->A, j, offset + volt_sources_cnt, -1.0);
				cs_entry(mna->sp_matrix->A, offset + volt_sources_cnt, j, -1.0);
				mna->b[offset + volt_sources_cnt] += value; 
			}
			else if (probe2_id == 0) {
				cs_entry(mna->sp_matrix->A, i, offset + volt_sources_cnt,  1.0);
				cs_entry(mna->sp_matrix->A, offset + volt_sources_cnt, i,  1.0);
				mna->b[offset + volt_sources_cnt] += value;
			}
			else {
				cs_entry(mna->sp_matrix->A, i, offset + volt_sources_cnt,  1.0);
				cs_entry(mna->sp_matrix->A, j, offset + volt_sources_cnt, -1.0);
				cs_entry(mna->sp_matrix->A, offset + volt_sources_cnt, i,  1.0);
				cs_entry(mna->sp_matrix->A, offset + volt_sources_cnt, j, -1.0);
				mna->b[offset + volt_sources_cnt] += value;
			}
			/* Keep track of how many voltage sources or inductors (which are treated like voltages with 0), we have already found */
			volt_sources_cnt++;
		}
	}
	cs *C = cs_compress(mna->sp_matrix->A);
	cs_spfree(mna->sp_matrix->A);
	mna->sp_matrix->A = C;
	cs_dupl(mna->sp_matrix->A);
	if (options->ITER) {
		/* Compute the M Jacobi Preconditioner */
		jacobi_precond(mna->M, NULL, C, mna->dimension, options->SPARSE);
	}
}

void create_sparse_trans_mna(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, double tr_step, int offset) {


}

/* Gets the response value according to the transient spec */
double get_response_value(list1_t *curr) {
	if (curr->type == 'I' || curr-> type == 'i' || curr->type == 'V' || curr->type == 'v') {
		if(curr->trans_spec != NULL) {
			switch (curr->trans_spec->type) {
				case EXP:
					return curr->trans_spec->exp->i1;
				case SIN:
					return curr->trans_spec->sin->i1;
				case PULSE:
					return curr->trans_spec->pulse->i1;
				case PWL:
					return curr->trans_spec->pwl->i[0];
			}
		}
		else {
			return curr->value;
		}
	}
	return 0.0;
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
	int maxiter, iterations;
	if (options->SPARSE) {
		if (options->ITER) {
			if (options->SPD) {
				/* Set the maximum number of iterations CG worst case is O(n) */
				maxiter = mna->dimension;
			 	iterations = conj_grad(NULL, mna->sp_matrix->A, *x, mna->b, mna->M, mna->dimension,
			 						   options->ITOL, maxiter, options->SPARSE);
				//printf("Conjugate gradient method did %d iterations.\n", iterations);
			}
			else {
				/* Set the maximum number of iterations Bi-CG worst case is O(2n) */
				maxiter = 2 * mna->dimension;
				iterations = bi_conj_grad(NULL, mna->sp_matrix->A, *x, mna->b, mna->M, mna->dimension, 
										  options->ITOL, maxiter, options->SPARSE);
				if (iterations == FAILURE) {
					fprintf(stderr, "Bi-Conjugate gradient method failed.\n");
					exit(EXIT_FAILURE);
				}
				//printf("Bi-Conjugate gradient method did %d iterations.\n", iterations);
			}
		}
		else {
			if (options->SPD) {
				solve_sparse_cholesky(mna, x);
			}
			else {
				solve_sparse_lu(mna, x);
			}
		}
	}
	else {
		/* Pointer to set the appropriate matrix */
		double **matrix_ptr;
		double *M_precond;
		if (mna->tr_analysis_init) {
			matrix_ptr = mna->matrix->aGhC;
			M_precond  = mna->M_trans;
		}
		else {
			matrix_ptr = mna->matrix->A;
			M_precond  = mna->M;
		}
		gsl_vector_view view_x = gsl_vector_view_array(*x, mna->dimension);
		if (options->ITER) {
			if (options->SPD) {
				/* Set the maximum number of iterations CG worst case is O(n) */
				maxiter = mna->dimension;
			 	iterations = conj_grad(matrix_ptr, NULL, *x, mna->b, M_precond, mna->dimension,
			 						   options->ITOL, maxiter, options->SPARSE);
				//printf("Conjugate gradient method did %d iterations.\n", iterations);
			}
			else {
				/* Set the maximum number of iterations Bi-CG worst case is O(2n) */
				maxiter = 2 * mna->dimension;
				iterations = bi_conj_grad(matrix_ptr, NULL, *x, mna->b, M_precond, mna->dimension, 
										  options->ITOL, maxiter, options->SPARSE);
				if (iterations == FAILURE) {
					fprintf(stderr, "Bi-Conjugate gradient method failed.\n");
					exit(EXIT_FAILURE);
				}
				//printf("Bi-Conjugate gradient method did %d iterations.\n", iterations);
			}
		}
		else {
			if (options->SPD) {
				solve_cholesky(matrix_ptr, mna->b, view_x, mna->dimension, mna->is_decomp);
			}
			else {
				solve_lu(matrix_ptr, mna->b, view_x, mna->matrix->P, mna->dimension, mna->is_decomp);
			}
		}
	}
	if (!mna->is_decomp) {
		mna->is_decomp = true;
		if (!options->TRAN) {
    		printf("Solution of MNA system... OK\n");
    	}
	}
}

/* Solves the sparse mna system with LU factorization */
void solve_sparse_lu(mna_system_t *mna, double **x) {
	double *temp_b = (double *)malloc(mna->dimension * sizeof(double));
	memcpy(temp_b, mna->b, mna->dimension * sizeof(double));
	if (!mna->is_decomp) {
		mna->sp_matrix->A_symbolic = cs_sqr(2, mna->sp_matrix->A, 0);
		mna->sp_matrix->A_numeric = cs_lu(mna->sp_matrix->A, mna->sp_matrix->A_symbolic, 1);
		cs_spfree(mna->sp_matrix->A);
	}
	cs_ipvec(mna->sp_matrix->A_numeric->pinv, temp_b, *x, mna->dimension);
	cs_lsolve(mna->sp_matrix->A_numeric->L, *x);
	cs_usolve(mna->sp_matrix->A_numeric->U, *x);
	cs_ipvec(mna->sp_matrix->A_symbolic->q, *x, temp_b, mna->dimension);
	memcpy(*x, temp_b, mna->dimension * sizeof(double));
	free(temp_b);
}

/* Solves the sparse mna system with Cholesky factorization */
void solve_sparse_cholesky(mna_system_t *mna, double **x) {
	double *temp_b = (double *)malloc(mna->dimension * sizeof(double));
	memcpy(temp_b, mna->b, mna->dimension * sizeof(double));
	if (!mna->is_decomp) {
		mna->sp_matrix->A_symbolic = cs_schol(1, mna->sp_matrix->A);
		mna->sp_matrix->A_numeric = cs_chol(mna->sp_matrix->A, mna->sp_matrix->A_symbolic);
		cs_spfree(mna->sp_matrix->A);
	}
	cs_ipvec(mna->sp_matrix->A_symbolic->pinv, temp_b, *x, mna->dimension);
	cs_lsolve(mna->sp_matrix->A_numeric->L, *x);
	cs_ltsolve(mna->sp_matrix->A_numeric->L, *x);
	cs_pvec(mna->sp_matrix->A_symbolic->pinv, *x, temp_b, mna->dimension);
	memcpy(*x, temp_b, mna->dimension * sizeof(double));
	free(temp_b);
}

/* Solve the MNA system using LU decomposition */
void solve_lu(double **A, double *b, gsl_vector_view x, gsl_permutation *P, int dimension, bool is_decomp) {
	/* The sign of the permutation matrix */
	int signum;
	/* Allocate memory for the solution vector */
	gsl_matrix_view view_A = gsl_matrix_view_array(A[0], dimension, dimension);
	gsl_vector_view view_b = gsl_vector_view_array(b, dimension);
	if (!is_decomp) {
		/* LU decomposition on A, PA = LU */
		gsl_linalg_LU_decomp(&view_A.matrix, P, &signum);
		// printf("LU Matrix:\n\n");
		// print_array(mna->matrix->A, mna->dimension);
		// printf("Permutation Vector:\n\n");
		// print_permutation(mna->matrix->P);
	}
	/* Solve the LU system */
	gsl_linalg_LU_solve(&view_A.matrix, P, &view_b.vector, &x.vector);
}

/* Solve the MNA system using cholesky decomposition */
void solve_cholesky(double **A, double *b, gsl_vector_view x, int dimension, bool is_decomp) {
	/* Allocate memory for the solution vector */
	gsl_matrix_view view_A = gsl_matrix_view_array(A[0], dimension, dimension);
	gsl_vector_view view_b = gsl_vector_view_array(b, dimension);
	if (!is_decomp) {
		/* Cholesky decomposition A = LL^T*/
		gsl_linalg_cholesky_decomp(&view_A.matrix);
		// printf("Cholesky Matrix:\n\n");
		// print_array(mna->matrix->A, mna->dimension);
		// printf("\n\n");
	}
	/* Solve the cholesky system */
	gsl_linalg_cholesky_solve(&view_A.matrix, &view_b.vector, &x.vector);
}

/* Print the MNA system */
void print_mna_system(mna_system_t *mna, options_t *options) {
	if (!options->SPARSE) {
		printf("\nMNA A array:\n\n");
		print_array(mna->matrix->A, mna->dimension);
		printf("MNA b vector:\n\n");
		print_vector(mna->b, mna->dimension);
		if (options->TRAN) {
			printf("MNA init_response vector:\n\n");
			print_vector(mna->matrix->resp->value, mna->dimension);
			printf("\nMNA G array:\n\n");
			print_array(mna->matrix->G, mna->dimension);
			printf("\nMNA hC array:\n\n");
			print_array(mna->matrix->hC, mna->dimension);
			printf("\nMNA aGhC array:\n\n");
			print_array(mna->matrix->aGhC, mna->dimension);
			if (options->TR) {
				printf("\nMNA sGhC array:\n\n");
				print_array(mna->matrix->sGhC, mna->dimension);
			}
		}
	}
	else {
		// cs_print(mna->sp_matrix->A, "sparse.txt", 0);
		printf("\nMNA b vector:\n\n");
		print_vector(mna->b, mna->dimension);
	}
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
void free_mna_system(mna_system_t **mna, options_t *options) {
	/* Free everything we allocated for the MNA and GSL */
	if (!options->SPARSE) {
		free((*mna)->matrix->A[0]);
		free((*mna)->matrix->A);
		gsl_permutation_free((*mna)->matrix->P);
		/* Free the matrix struct */
		free((*mna)->matrix);
	}
	else {
		if (!options->ITER) {
			cs_sfree((*mna)->sp_matrix->A_symbolic);
			cs_nfree((*mna)->sp_matrix->A_numeric);
		}
		else {
			cs_spfree((*mna)->sp_matrix->A);
		}
		free((*mna)->sp_matrix);
	}
	if (options->ITER) {
		free((*mna)->M);
	}
	free((*mna)->b);
	/* Free every string allocated for the group2 elements */
	for (int i = 0; i < (*mna)->num_g2_elem; i++) {
		free((*mna)->g2_indx[i].element);
	}
	/* Free the g2 array */
	free((*mna)->g2_indx);
	free(*mna);
	/* Set mna to NULL to limit further acesses */
	*mna = NULL;
}