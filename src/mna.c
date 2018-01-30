#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "mna.h"

/* Allocate memory for the MNA system */
mna_system_t *init_mna_system(int num_nodes, int num_g2_elem, options_t *options, int nz) {
	/* Allocate for the whole struct */
	mna_system_t *mna = (mna_system_t *)malloc(sizeof(mna_system_t));
	assert(mna != NULL);

	/* Dimension will be (n-1) + m2 */
	mna->dimension = num_nodes + num_g2_elem;

	if (options->SPARSE) {
		init_sparse_matrix(mna, options, nz);
	}
	else {
		init_dense_matrix(mna, options);
	}

	mna->resp = NULL;
	/* Set the preconditioner pointers to NULL */
	mna->M 	  = mna->M_trans   = NULL; 
	mna->M_ac = mna->M_ac_conj = NULL;

	/* In case we will use iterative methods allocate memory for the prerequisites */
	if (options->ITER) {
		// TODO PERHAPS REMOVE THE BELOW and just allocate preconditioner sets by default to 1.0 in case it founds 0.0
		mna->M = init_val_vector(mna->dimension, 1.0);
		if (options->TRAN) {
			mna->M_trans = init_val_vector(mna->dimension, 1.0);
		}
		if (options->AC) {
			mna->M_ac      = init_gsl_complex_vector(mna->dimension);
			mna->M_ac_conj = init_gsl_complex_vector(mna->dimension);
		}
	}

	/* 
	 * In case netlist has .TRAN, allocate the resp_t struct that holds the transient response or DC value
	 * and the nodes that contribute to it. It's the same for dense and sparse matrices.
	 */
	if (options->TRAN) {
		mna->resp = (resp_t *)malloc(sizeof(resp_t));
		assert(mna->resp != NULL);
		mna->resp->value = init_vector(mna->dimension);
		/* Allocate enough memory to hold at max n pointers to elements */
		mna->resp->nodes = (list1_t **)malloc(mna->dimension * sizeof(list1_t *));
		assert(mna->resp->nodes != NULL);
		/* Initialize all the pointers to NULL, only those that contribute will have a non-NULL value at the end */
		for (int i = 0; i < mna->dimension; i++) {
			mna->resp->nodes[i] = NULL;
		}
	}
	/* Allocate the rhs vector */
	mna->b = init_vector(mna->dimension);

	/* Initialize the other fields of the mna system */
	mna->is_decomp        = false;
	mna->tr_analysis_init = false;
	mna->ac_analysis_init = false;
	mna->num_nodes        = num_nodes;
	mna->num_g2_elem      = num_g2_elem;

	/* Allocate a vector that holds the names of g2 elements */
	mna->g2_indx = (g2_indx_t *)malloc(num_g2_elem * sizeof(g2_indx_t));
	assert(mna->g2_indx != NULL);

	return mna;
}

/* Allocate memory for a dense representation of MNA */
void init_dense_matrix(mna_system_t *mna, options_t *options) {
	/* Allocate for the matrices struct */
	mna->matrix = (matrix_t *)malloc(sizeof(matrix_t));
	assert(mna->matrix != NULL);

	mna->matrix->A = init_array(mna->dimension, mna->dimension);
	mna->matrix->P = init_permutation(mna->dimension);

	/* Initialize to NULL everything */
	mna->matrix->hC     = NULL;
	mna->matrix->aGhC   = NULL;
	mna->matrix->sGhC   = NULL;
	mna->matrix->G_ac   = NULL;
	mna->matrix->e_ac   = NULL;
	mna->matrix->A_base = NULL;

	mna->sp_matrix = NULL;

	if (options->TRAN) {
		mna->matrix->hC   = init_array(mna->dimension, mna->dimension);
		mna->matrix->aGhC = init_array(mna->dimension, mna->dimension);
		/* If method is Backward Euler we don't need sGhC matrix, we use hC instead */
		if (options->TR) {
			mna->matrix->sGhC = init_array(mna->dimension, mna->dimension);
		}
	}

	if (options->AC) {
		mna->matrix->G_ac = init_gsl_complex_array(mna->dimension, mna->dimension);
		mna->matrix->e_ac = init_gsl_complex_vector(mna->dimension);
	}

	if (options->TRAN || options->AC) {
		mna->matrix->A_base = init_array(mna->dimension, mna->dimension);
	}
}

/* Allocate memory for a sparse representation of MNA */
void init_sparse_matrix(mna_system_t *mna, options_t *options, int nz) {
	/* Allocate for the matrices struct */
	mna->sp_matrix = (sp_matrix_t *)malloc(sizeof(sp_matrix_t));
	assert(mna->sp_matrix != NULL);

	mna->sp_matrix->A = cs_spalloc(mna->dimension, mna->dimension, nz, 1, 1);
	assert(mna->sp_matrix->A != NULL);

	/* Initialize to NULL everything */
	mna->sp_matrix->A->nz  = 0;
	mna->sp_matrix->hC     = NULL;
	mna->sp_matrix->aGhC   = NULL;
	mna->sp_matrix->sGhC   = NULL;
	mna->sp_matrix->G_ac   = NULL;
	mna->sp_matrix->e_ac   = NULL;
	mna->sp_matrix->A_base = NULL;

	mna->matrix = NULL;

	if (options->TRAN) {
		mna->sp_matrix->hC = cs_di_spalloc(mna->dimension, mna->dimension, nz, 1, 1);
	}

	if (options->AC) {
		mna->sp_matrix->e_ac = (cs_complex_t *)malloc(mna->dimension * sizeof(cs_complex_t));
	}
}

/* Constructs the MNA system either sparse or dense for DC or TRAN analysis */
void create_mna_system(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, double tr_step, int offset) {
	if (options->SPARSE) {
		create_sparse_mna(mna, index, hash_table, options, offset);
		if (options->TRAN) {
			create_sparse_trans_mna(mna, index, hash_table, options, tr_step, offset);
		}
	}
	else {
		create_dense_mna(mna, index, hash_table, options, offset);
		if (options->TRAN) {
			create_dense_trans_mna(mna, index, hash_table, options, tr_step, offset);
		}
	}
	printf("\nCreation of MNA system...OK\n");
}

/* Constructs the AC MNA system either sparse or dense for AC analysis */
void create_ac_mna_system(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, int offset, double omega) {
	if (options->SPARSE) {
		create_sparse_ac_mna(mna, index, hash_table, options, offset, omega);
	}
	else {
		create_dense_ac_mna(mna, index, hash_table, options, offset, omega);
	}
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
			mna->g2_indx[volt_sources_cnt].element = (char *)malloc((strlen(curr->element) + 1) * sizeof(char));
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

	/* Copy the data from A to A_base before LU factorization to build matrices in AC and TRAN analysis */
	if (options->TRAN || options->AC) {
		for (int i = 0; i < mna->dimension; i++) {
			memcpy(mna->matrix->A_base[i], mna->matrix->A[i], mna->dimension * sizeof(double));
		}
	}
}

/* Constructs the dense transient MNA system */
void create_dense_trans_mna(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, double tr_step, int offset) {
	list1_t *curr;
	double value, h;
	int volt_sources_cnt = 0;

	/* Set the sampling step h according to the specified method (TRAP or BE) */
	if (options->TR) {
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
			/* Set the L value in the diagonal of g2 area in matrix */
			mna->matrix->hC[offset + volt_sources_cnt][offset + volt_sources_cnt] = -curr->value * h;
			/* Keep track of how many voltage sources or inductors, we have already found */
			volt_sources_cnt++;
		}
		else if (curr->type == 'V' || curr->type == 'v') {
			value = get_response_value(curr);
			/* Keep track of the node that contributes */
			mna->resp->nodes[offset + volt_sources_cnt] = curr;
			/* Keep track of how many voltage sources or inductors, we have already found */
			if (probe1_id == 0) {
				mna->resp->value[offset + volt_sources_cnt] += value;
			}
			else if (probe2_id == 0) {
				mna->resp->value[offset + volt_sources_cnt] += value;
			}
			else {
				mna->resp->value[offset + volt_sources_cnt] += value;
			}
			/* Keep track of how many voltage sources or inductors, we have already found */
			volt_sources_cnt++;
		}
		else if (curr->type == 'I' || curr->type == 'i') {
			value = get_response_value(curr);
			if (probe1_id == 0) {
				mna->resp->value[j] += value;
				mna->resp->nodes[j]  = curr;
			}
			else if (probe2_id == 0) {
				mna->resp->value[i] -= value;
				mna->resp->nodes[i]  = curr;
			}
			else {
				mna->resp->value[i] -= value;
				mna->resp->value[j] += value;
				mna->resp->nodes[i]  = curr;
				mna->resp->nodes[j]  = curr;
			}
		}
	}

	/* Compute the matrices G + 1/h*C and G - 2/h*C */
	for (int i = 0; i < mna->dimension; i++) {
		for (int j = 0; j < mna->dimension; j++) {
			mna->matrix->aGhC[i][j] = mna->matrix->A_base[i][j] + mna->matrix->hC[i][j];
			if (options->TR) {
				mna->matrix->sGhC[i][j] = mna->matrix->A_base[i][j] - mna->matrix->hC[i][j];
			}
		}
	}

	if (options->ITER) {
		/* Compute the M preconditioner */
		jacobi_precond(mna->M_trans, mna->matrix->aGhC, NULL, mna->dimension, options->SPARSE);
	}

	/* Free what is no longer needed, in case it's trapezoidal method */
	if (options->TR) {
		free(mna->matrix->hC[0]);
		free(mna->matrix->hC);
	}
}

/* Constructs the dense AC MNA system */
void create_dense_ac_mna(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, int offset, double omega) {
	list1_t *curr;
	int volt_sources_cnt = 0;

	/* Copy the base for the AC matrix and zero-out the vector for the next step */
	for (int i = 0; i < mna->dimension; i++) {
		gsl_complex z;
		for (int j = 0; j < mna->dimension; j++) {
			// mna->matrix->G_ac[i][j] = mna->matrix->A_base[i][j];
			z = gsl_complex_rect(mna->matrix->A_base[i][j], 0);
			gsl_matrix_complex_set(mna->matrix->G_ac, i, j, z);
		}
		// mna->matrix->e_ac[i] = 0.0 + I * 0.0;
	}
	gsl_vector_complex_set_zero(mna->matrix->e_ac);

	for (curr = index->head1; curr != NULL; curr = curr->next) {
		int probe1_id = ht_get_id(hash_table, curr->probe1);
		int probe2_id = ht_get_id(hash_table, curr->probe2);
		int i = probe1_id - 1;
		int j = probe2_id - 1;
		gsl_complex z, z_prev;
		if (curr->type == 'C' || curr->type == 'c') {
			// value = I * omega * curr->value;
			z = gsl_complex_rect(0.0, omega * curr->value);
			if (probe1_id == 0) {
				// mna->matrix->G_ac[j][j] += value;
				z_prev = gsl_matrix_complex_get(mna->matrix->G_ac, j, j);
				gsl_matrix_complex_set(mna->matrix->G_ac, j, j, gsl_complex_add(z_prev, z));
			}
			else if (probe2_id == 0) {
				// mna->matrix->G_ac[i][i] += value;
				z_prev = gsl_matrix_complex_get(mna->matrix->G_ac, i, i);
				gsl_matrix_complex_set(mna->matrix->G_ac, i, i, gsl_complex_add(z_prev, z));
			}
			else {
				// mna->matrix->G_ac[i][i] += value;
				z_prev = gsl_matrix_complex_get(mna->matrix->G_ac, i, i);
				gsl_matrix_complex_set(mna->matrix->G_ac, i, i, gsl_complex_add(z_prev, z));
				// mna->matrix->G_ac[j][j] += value;
				z_prev = gsl_matrix_complex_get(mna->matrix->G_ac, j, j);
				gsl_matrix_complex_set(mna->matrix->G_ac, j, j, gsl_complex_add(z_prev, z));
				// mna->matrix->G_ac[i][j] -= value;
				z_prev = gsl_matrix_complex_get(mna->matrix->G_ac, i, j);
				gsl_matrix_complex_set(mna->matrix->G_ac, i, j, gsl_complex_sub(z_prev, z));
				// mna->matrix->G_ac[j][i] -= value;
				z_prev = gsl_matrix_complex_get(mna->matrix->G_ac, j, i);
				gsl_matrix_complex_set(mna->matrix->G_ac, j, i, gsl_complex_sub(z_prev, z));
			}
		}
		else if (curr->type == 'I' || curr->type == 'i') {
			if (curr->ac == NULL) {
				// value = 0;
				z = gsl_complex_rect(0.0, 0.0);
			}
			else {
				// value = pol_to_rect(curr->ac->magnitude, curr->ac->phase);
				z = gsl_complex_polar(curr->ac->magnitude, curr->ac->phase);
			}
			if (probe1_id == 0) {
				// mna->matrix->e_ac[j] += value;
				z_prev = gsl_vector_complex_get(mna->matrix->e_ac, j);
				gsl_vector_complex_set(mna->matrix->e_ac, j, gsl_complex_add(z_prev, z));
			}
			else if (probe2_id == 0) {
				// mna->matrix->e_ac[i] -= value;
				z_prev = gsl_vector_complex_get(mna->matrix->e_ac, i);
				gsl_vector_complex_set(mna->matrix->e_ac, i, gsl_complex_sub(z_prev, z));
			}
			else {
				// mna->matrix->e_ac[i] -= value;
				z_prev = gsl_vector_complex_get(mna->matrix->e_ac, i);
				gsl_vector_complex_set(mna->matrix->e_ac, i, gsl_complex_sub(z_prev, z));
				// mna->matrix->e_ac[j] += value;
				z_prev = gsl_vector_complex_get(mna->matrix->e_ac, j);
				gsl_vector_complex_set(mna->matrix->e_ac, j, gsl_complex_add(z_prev, z));
			}
		}
		else if (curr->type == 'L' || curr->type == 'l') {
			/* Set the L value in the diagonal of g2 area in matrices */
			// mna->matrix->G_ac[offset + volt_sources_cnt][offset + volt_sources_cnt] = - I * omega * curr->value;
			z = gsl_complex_rect(0.0, -omega * curr->value);
			gsl_matrix_complex_set(mna->matrix->G_ac, offset + volt_sources_cnt, offset + volt_sources_cnt, z);
			/* Keep track of how many voltage sources or inductors (which are treated like voltages with 0), we have already found */
			volt_sources_cnt++;
		}
		else if (curr->type == 'V' || curr->type == 'v') {
			if (curr->ac == NULL) {
				// value = 0;
				z = gsl_complex_rect(0.0, 0.0);
			}
			else {
				// value = pol_to_rect(curr->ac->magnitude, curr->ac->phase);
				z = gsl_complex_polar(curr->ac->magnitude, curr->ac->phase);
			}
			/* Keep track of how many voltage sources or inductors, we have already found */
			// mna->matrix->e_ac[offset + volt_sources_cnt] += value;
			z_prev = gsl_vector_complex_get(mna->matrix->e_ac, offset + volt_sources_cnt);
			gsl_vector_complex_set(mna->matrix->e_ac, offset + volt_sources_cnt, gsl_complex_add(z_prev, z));
			/* Keep track of how many voltage sources or inductors, we have already found */
			volt_sources_cnt++;
		}
	}

	if (options->ITER) {
		/* Compute the M_ac M_ac_conj preconditioners */
		complex_jacobi_precond(mna->M_ac, mna->matrix->G_ac, NULL, mna->dimension, options->SPARSE);
		vector_conjugate(mna->M_ac_conj, mna->M_ac, mna->dimension);
	}
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
				cs_di_entry(mna->sp_matrix->A, j, j, value);
			}
			else if (probe2_id == 0) {
				cs_di_entry(mna->sp_matrix->A, i, i, value);
			}
			else {
				cs_di_entry(mna->sp_matrix->A, i, i,  value);
				cs_di_entry(mna->sp_matrix->A, j, j,  value);
				cs_di_entry(mna->sp_matrix->A, i, j, -value);
				cs_di_entry(mna->sp_matrix->A, j, i, -value);
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
			mna->g2_indx[volt_sources_cnt].element = (char *)malloc((strlen(curr->element) + 1) * sizeof(char));
			assert(mna->g2_indx[volt_sources_cnt].element != NULL);
			strcpy(mna->g2_indx[volt_sources_cnt].element, curr->element);
			if (probe1_id == 0) {
				cs_di_entry(mna->sp_matrix->A, j, offset + volt_sources_cnt, -1.0);
				cs_di_entry(mna->sp_matrix->A, offset + volt_sources_cnt, j, -1.0);
				mna->b[offset + volt_sources_cnt] += value; 
			}
			else if (probe2_id == 0) {
				cs_di_entry(mna->sp_matrix->A, i, offset + volt_sources_cnt,  1.0);
				cs_di_entry(mna->sp_matrix->A, offset + volt_sources_cnt, i,  1.0);
				mna->b[offset + volt_sources_cnt] += value;
			}
			else {
				cs_di_entry(mna->sp_matrix->A, i, offset + volt_sources_cnt,  1.0);
				cs_di_entry(mna->sp_matrix->A, j, offset + volt_sources_cnt, -1.0);
				cs_di_entry(mna->sp_matrix->A, offset + volt_sources_cnt, i,  1.0);
				cs_di_entry(mna->sp_matrix->A, offset + volt_sources_cnt, j, -1.0);
				mna->b[offset + volt_sources_cnt] += value;
			}
			/* Keep track of how many voltage sources or inductors (which are treated like voltages with 0), we have already found */
			volt_sources_cnt++;
		}
	}

	/* Copy the data from A to A_base before LU factorization to build matrices in AC */
	if (options->AC) {
		mna->sp_matrix->A_base = _cs_di_copy(mna->sp_matrix->A);
	}

	/* Convert to compressed-column format (CCF) from triplet */
	cs_di *C = cs_di_compress(mna->sp_matrix->A);
	cs_di_spfree(mna->sp_matrix->A);
	mna->sp_matrix->A = C;
	cs_di_dupl(mna->sp_matrix->A);

	if (options->ITER) {
		/* Compute the M Jacobi Preconditioner */
		jacobi_precond(mna->M, NULL, C, mna->dimension, options->SPARSE);
	}
}

void create_sparse_trans_mna(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, double tr_step, int offset) {
	list1_t *curr;
	double value, h;
	int volt_sources_cnt = 0;

	// memcpy(*mna->sp_matrix->G, *mna->sp_matrix->A, mna->dimension * mna->dimension * sizeof(double));
	// assert(mna->sp_matrix->G != NULL);

	/* Set the sampling step h according to the specified method (TRAP or BE) */
	if (options->TR) {
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
				cs_di_entry(mna->sp_matrix->hC, j, j, value);
			}
			else if (probe2_id == 0) {
				cs_di_entry(mna->sp_matrix->hC, i, i, value);
			}
			else {
				cs_di_entry(mna->sp_matrix->hC, i, i,  value);
				cs_di_entry(mna->sp_matrix->hC, j, j,  value);
				cs_di_entry(mna->sp_matrix->hC, i, j, -value);
				cs_di_entry(mna->sp_matrix->hC, j, i, -value);
			}
		}
		else if(curr->type == 'L' || curr->type == 'l') {
			/* Set the L value in the diagonal of g2 area in matrices */
			cs_di_entry(mna->sp_matrix->hC, offset + volt_sources_cnt, offset + volt_sources_cnt, -curr->value * h);
			/* Keep track of how many voltage sources or inductors, we have already found */
			volt_sources_cnt++;
		}
		else if (curr->type == 'V' || curr->type == 'v') {
			value = get_response_value(curr);
			/* Keep track of the node that contributes */
			mna->resp->nodes[offset + volt_sources_cnt] = curr;
			/* Keep track of how many voltage sources or inductors, we have already found */
			if (probe1_id == 0) {
				mna->resp->value[offset + volt_sources_cnt] += value;
			}
			else if (probe2_id == 0) {
				mna->resp->value[offset + volt_sources_cnt] += value;
			}
			else {
				mna->resp->value[offset + volt_sources_cnt] += value;
			}
			/* Keep track of how many voltage sources or inductors, we have already found */
			volt_sources_cnt++;
		}
		else if (curr->type == 'I' || curr->type == 'i') {
			value = get_response_value(curr);
			if (probe1_id == 0) {
				mna->resp->value[j] += value;
				mna->resp->nodes[j]  = curr;
			}
			else if (probe2_id == 0) {
				mna->resp->value[i] -= value;
				mna->resp->nodes[i]  = curr;
			}
			else {
				mna->resp->value[i] -= value;
				mna->resp->value[j] += value;
				mna->resp->nodes[i]  = curr;
				mna->resp->nodes[j]  = curr;
			}
		}
	}

	/* Compress the hC matrix i.e. convert it to CCF from triplet */
	cs *hC_comp = cs_compress(mna->sp_matrix->hC);
	cs_spfree(mna->sp_matrix->hC);
	mna->sp_matrix->hC = hC_comp;
	cs_dupl(mna->sp_matrix->hC);

	/* Compute the matrices G + 1/h*C and G - 2/h*C */
	/* cs_add takes compressed column matrices, the result is also compressed column matrix */
	mna->sp_matrix->aGhC = cs_add(mna->sp_matrix->A, mna->sp_matrix->hC, 1, 1);
	assert(mna->sp_matrix->aGhC != NULL);

	/* In case Trapezoidal method is specified */
	if (options->TR) {
		/* cs_add takes compressed column matrices, the result is also compressed column matrix */
		mna->sp_matrix->sGhC = cs_add(mna->sp_matrix->A, mna->sp_matrix->hC, 1, -1);
		assert(mna->sp_matrix->aGhC != NULL);
		/* Free what is no longer needed */
		cs_spfree(mna->sp_matrix->hC);
	}
	
	if (options->ITER) {
		/* Compute the M preconditioner */
		jacobi_precond(mna->M_trans, NULL, mna->sp_matrix->aGhC, mna->dimension, options->SPARSE);
	}
}

/* Constructs the sparse AC MNA system */
void create_sparse_ac_mna(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, int offset, double omega) {
	list1_t *curr;
	int volt_sources_cnt = 0;
	
	/*
	 * We need to free everytime the G_ac matrix because in every sweep step we convert the A_base
	 * into a complex sparse matrix G_ac. Thus, every time the last pointer value would be lost.
	 */
	if (mna->sp_matrix->G_ac != NULL) {
		cs_ci_spfree(mna->sp_matrix->G_ac);
	}

	/* Copy the base for the AC matrix and zero-out the vector for the next step */
	mna->sp_matrix->G_ac = cs_i_complex(mna->sp_matrix->A_base, 1);

	//TODO create a function for the below
	for (int i = 0; i < mna->dimension; i++) {
		mna->sp_matrix->e_ac[i] = 0.0 + I * 0.0;
	}

	for (curr = index->head1; curr != NULL; curr = curr->next) {
		int probe1_id = ht_get_id(hash_table, curr->probe1);
		int probe2_id = ht_get_id(hash_table, curr->probe2);
		int i = probe1_id - 1;
		int j = probe2_id - 1;
		cs_complex_t z;
		if (curr->type == 'C' || curr->type == 'c') {
			// value = I * omega * curr->value;
			z = 0.0 + (omega * curr->value * I);
			if (probe1_id == 0) {
				// mna->matrix->G_ac[j][j] += value;
				cs_ci_entry(mna->sp_matrix->G_ac, j, j, z);
			}
			else if (probe2_id == 0) {
				// mna->matrix->G_ac[i][i] += value;
				cs_ci_entry(mna->sp_matrix->G_ac, i, i, z);
			}
			else {
				// mna->matrix->G_ac[i][i] += value;
				cs_ci_entry(mna->sp_matrix->G_ac, i, i, z);
				// mna->matrix->G_ac[j][j] += value;
				cs_ci_entry(mna->sp_matrix->G_ac, j, j, z);
				// mna->matrix->G_ac[i][j] -= value;
				cs_ci_entry(mna->sp_matrix->G_ac, i, j, CS_COMPLEX_NEG(z));
				// mna->matrix->G_ac[j][i] -= value;
				cs_ci_entry(mna->sp_matrix->G_ac, j, i, CS_COMPLEX_NEG(z));
			}
		}
		else if (curr->type == 'I' || curr->type == 'i') {
			/* Check if ac spec exists for the current element */
			if (curr->ac == NULL) {
				z = 0.0 + 0.0 * I;
			}
			else {
				z = pol_to_rect(curr->ac->magnitude, curr->ac->phase);
			}
			if (probe1_id == 0) {
				mna->sp_matrix->e_ac[j] += z;
			}
			else if (probe2_id == 0) {
				mna->sp_matrix->e_ac[j] -= z;
			}
			else {
				mna->sp_matrix->e_ac[i] -= z;
				mna->sp_matrix->e_ac[j] += z;
			}
		}
		else if (curr->type == 'L' || curr->type == 'l') {
			/* Set the L value in the diagonal of g2 area in matrices */
			// mna->matrix->G_ac[offset + volt_sources_cnt][offset + volt_sources_cnt] = - I * omega * curr->value;
			z = 0.0 - (omega * curr->value * I);
			cs_ci_entry(mna->sp_matrix->G_ac, offset + volt_sources_cnt, offset + volt_sources_cnt, z);
			/* Keep track of how many voltage sources or inductors (which are treated like voltages with 0), we have already found */
			volt_sources_cnt++;
		}
		else if (curr->type == 'V' || curr->type == 'v') {
			/* Check if ac spec exists for the current element */
			if (curr->ac == NULL) {
				z = 0.0 + 0.0 * I;
			}
			else {
				z = pol_to_rect(curr->ac->magnitude, curr->ac->phase);
			}
			/* Keep track of how many voltage sources or inductors, we have already found */
			mna->sp_matrix->e_ac[offset + volt_sources_cnt] += z;
			/* Keep track of how many voltage sources or inductors, we have already found */
			volt_sources_cnt++;
		}
	}

	/* Convert to compressed-column format (CCF) from triplet */
	cs_ci *C = cs_ci_compress(mna->sp_matrix->G_ac);
	cs_ci_spfree(mna->sp_matrix->G_ac);
	mna->sp_matrix->G_ac = C;
	cs_ci_dupl(mna->sp_matrix->G_ac);

	if (options->ITER) {
		/* Compute the M Jacobi Preconditioner */
		gsl_vector_complex_set_all(mna->M_ac, GSL_COMPLEX_ONE);
		complex_jacobi_precond(mna->M_ac, NULL, mna->sp_matrix->G_ac, mna->dimension, options->SPARSE);
		vector_conjugate(mna->M_ac_conj, mna->M_ac, mna->dimension);
	}
}

/* Gets the response value according to the transient spec */
double get_response_value(list1_t *curr) {
	if (curr->type == 'I' || curr-> type == 'i' || curr->type == 'V' || curr->type == 'v') {
		if (curr->trans_spec != NULL) {
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
void solve_mna_system(mna_system_t *mna, double **x, gsl_vector_complex *x_complex, options_t *options) {
	/* Set the maximum number of iterations CG/bi-CG */
	int iterations, maxiter = mna->dimension;

	if (options->SPARSE) {
		/* Pointer to set the appropriate matrix */
		cs *matrix_ptr = NULL;
		double *M_precond = NULL;;
		gsl_vector_complex *gsl_e_ac = NULL;
		cs_complex_t *cs_x_complex = NULL;

		/* In case we are not currently in an AC analysis */
		if (!mna->ac_analysis_init) {
			if (mna->tr_analysis_init) {
				matrix_ptr = mna->sp_matrix->aGhC;
				M_precond  = mna->M_trans;
			}
			/* We are in DC analysis */
			else {
				matrix_ptr = mna->sp_matrix->A;
				M_precond  = mna->M;
			}
		}
		if (options->ITER) {
			if (mna->ac_analysis_init) {
				/* Convert x vector (which is the DC operating point) to a complex one in case we're in an AC analysis */
				real_to_gsl_complex_vector(x_complex, *x, mna->dimension);
				/* Also conver cs_complex_t *e_ac into a gsl_vector for the iterative solvers */
				gsl_e_ac = init_gsl_complex_vector(mna->dimension);
				cs_complex_to_gsl(gsl_e_ac, mna->sp_matrix->e_ac, mna->dimension);
			}
			if (options->SPD) {
				/* Check if AC analysis has started (complex), otherwise call ordinary solvers */
				if (mna->ac_analysis_init) {
					iterations = complex_conj_grad(NULL, mna->sp_matrix->G_ac, x_complex, gsl_e_ac,
					 							   mna->M_ac, mna->dimension, options->ITOL, maxiter, options->SPARSE);
				}
				else {
					iterations = conj_grad(NULL, matrix_ptr, *x, mna->b, M_precond, mna->dimension,
			 							   options->ITOL, maxiter, options->SPARSE);
				}
				//printf("Conjugate gradient method did %d iterations.\n", iterations);
			}
			else {
				/* Check if AC analysis has started (complex), otherwise call ordinary solvers */
				if (mna->ac_analysis_init) {
					iterations = complex_bi_conj_grad(NULL, mna->sp_matrix->G_ac, x_complex, gsl_e_ac,
													  mna->M_ac, mna->M_ac_conj, mna->dimension, options->ITOL,
													  maxiter, options->SPARSE);
				}
				else {
					iterations = bi_conj_grad(NULL, matrix_ptr, *x, mna->b, M_precond, mna->dimension, 
										  	  options->ITOL, maxiter, options->SPARSE);
				}
				if (iterations == FAILURE) {
					fprintf(stderr, "Bi-Conjugate gradient method failed.\n");
					exit(EXIT_FAILURE);
				}
				// printf("Bi-Conjugate gradient method did %d iterations.\n", iterations);
			}
		}
		else { /* Non iterative solvers */
			if (mna->ac_analysis_init) {
				cs_x_complex = (cs_complex_t *)malloc(mna->dimension * sizeof(cs_complex_t));
			}
			if (options->SPD) {
				/* Check if AC analysis has started (complex), otherwise call ordinary solvers */
				if (mna->ac_analysis_init) {
					solve_complex_sparse_cholesky(mna, cs_x_complex);
				}
				else {
					solve_sparse_cholesky(mna, matrix_ptr, x);
				}
			}
			else {
				/* Check if AC analysis has started (complex), otherwise call ordinary solvers */
				if (mna->ac_analysis_init) {
					solve_complex_sparse_lu(mna, cs_x_complex);
				}
				else {
					solve_sparse_lu(mna, matrix_ptr, x);
				}
			}
			if (mna->ac_analysis_init) {
				cs_complex_to_gsl(x_complex, cs_x_complex, mna->dimension);
			}
		}
		/* Free the converted e_ac to gsl and cs_x_complex in case it was allocated for AC analysis */
		if (mna->ac_analysis_init) {
			gsl_vector_complex_free(gsl_e_ac);
			free(cs_x_complex);
		}
	}
	else { /* Dense */
		/* Pointer to set the appropriate matrix */
		double **matrix_ptr = NULL;
		double *M_precond   = NULL;
		gsl_vector_view view_x;
		/* Set general pointers for matrices, vectors to use for the solvers */
		if (!mna->ac_analysis_init) {
			if (mna->tr_analysis_init) {
				matrix_ptr = mna->matrix->aGhC;
				M_precond  = mna->M_trans;
			}
			else {
				matrix_ptr = mna->matrix->A;
				M_precond  = mna->M;
			}
			view_x = gsl_vector_view_array(*x, mna->dimension);
		}
		if (options->ITER) {
			/* Convert x vector (which is the DC operating point) to a complex one in case we're in an AC analysis */
			if (mna->ac_analysis_init) {
				real_to_gsl_complex_vector(x_complex, *x, mna->dimension);
			}
			if (options->SPD) {
				/* Check if AC analysis has started (complex), otherwise call ordinary solvers */
				if (mna->ac_analysis_init) {
					iterations = complex_conj_grad(mna->matrix->G_ac, NULL, x_complex, mna->matrix->e_ac, mna->M_ac, 
												   mna->dimension, options->ITOL, maxiter, options->SPARSE);
				}
				else {
			 		iterations = conj_grad(matrix_ptr, NULL, *x, mna->b, M_precond, mna->dimension,
			 							   options->ITOL, maxiter, options->SPARSE);
				}
				//printf("Conjugate gradient method did %d iterations.\n", iterations);
			}
			else {
				/* Check if AC analysis has started (complex), otherwise call ordinary solvers */
				if (mna->ac_analysis_init) {
					iterations = complex_bi_conj_grad(mna->matrix->G_ac, NULL, x_complex, mna->matrix->e_ac, mna->M_ac,
												   	  mna->M_ac_conj, mna->dimension, options->ITOL, maxiter, options->SPARSE);
				}
				else {
					iterations = bi_conj_grad(matrix_ptr, NULL, *x, mna->b, M_precond, mna->dimension, 
											  options->ITOL, maxiter, options->SPARSE);
				}
				if (iterations == FAILURE) {
					fprintf(stderr, "Bi-Conjugate gradient method failed.\n");
					exit(EXIT_FAILURE);
				}
				else if ((maxiter < MAX_ITER_THRESHOLD && iterations == MAX_ITER_THRESHOLD) || 
						 (maxiter > MAX_ITER_THRESHOLD && iterations == maxiter)) {
					printf("Bi-Conjugate gradient reached max iterations without convergence.\n");
				}
				//printf("Bi-Conjugate gradient method did %d iterations.\n", iterations);
			}
		}
		else {
			// TODO perhaps add INLINE MACRO for these flags i.e. mna->ac_analysis_init, mna->tr_analysis_init
			if (options->SPD) {
				/* Check if AC analysis has started (complex), otherwise call ordinary solvers */
				if (mna->ac_analysis_init) {
					solve_complex_cholesky(mna->matrix->G_ac, mna->matrix->e_ac, x_complex, mna->dimension);
				}
				else {
					solve_cholesky(matrix_ptr, mna->b, view_x, mna->dimension, mna->is_decomp);
				}
			}
			else {
				/* Check if AC analysis has started (complex), otherwise call ordinary solvers */
				if (mna->ac_analysis_init) {
					solve_complex_lu(mna->matrix->G_ac, mna->matrix->e_ac, x_complex, mna->matrix->P, mna->dimension);
				}
				else {
					solve_lu(matrix_ptr, mna->b, view_x, mna->matrix->P, mna->dimension, mna->is_decomp);
				}
			}
		}
	}
	if (!mna->is_decomp) {
		mna->is_decomp = true;
		if (!mna->ac_analysis_init && !mna->tr_analysis_init) {
    		printf("Solution of MNA system...OK\n");
    	}
	}
}

/* Solve the MNA system using LU decomposition and store the result in vector x */
void solve_lu(double **A, double *b, gsl_vector_view x, gsl_permutation *P, int dimension, bool is_decomp) {
	/* The sign of the permutation matrix */
	int signum;
	gsl_matrix_view view_A = gsl_matrix_view_array(A[0], dimension, dimension);
	gsl_vector_view view_b = gsl_vector_view_array(b, dimension);
	if (!is_decomp) {
		/* LU decomposition on A, PA = LU */
		gsl_linalg_LU_decomp(&view_A.matrix, P, &signum);
	}
	gsl_linalg_LU_solve(&view_A.matrix, P, &view_b.vector, &x.vector);
}

/* Solve the MNA system using complex LU decomposition and store the result in vector x */
void solve_complex_lu(gsl_matrix_complex *A, gsl_vector_complex *b, gsl_vector_complex *x, gsl_permutation *P, int dimension) {
	/* The sign of the permutation matrix */
	int signum;
	/* LU decomposition on A, PA = LU */
	gsl_linalg_complex_LU_decomp(A, P, &signum);
	gsl_linalg_complex_LU_solve(A, P, b, x);
}

/* Solves the sparse mna system with LU factorization and store the result in vector x */
void solve_sparse_lu(mna_system_t *mna, cs *A, double **x) {
	double *temp_b = (double *)malloc(mna->dimension * sizeof(double));
	memcpy(temp_b, mna->b, mna->dimension * sizeof(double));

	if (!mna->is_decomp) {
		mna->sp_matrix->A_symbolic = cs_sqr(2, A, 0);
		mna->sp_matrix->A_numeric  = cs_lu(A, mna->sp_matrix->A_symbolic, 1);
		cs_spfree(A);
	}

	cs_ipvec(mna->sp_matrix->A_numeric->pinv, temp_b, *x, mna->dimension);
	cs_lsolve(mna->sp_matrix->A_numeric->L, *x);
	cs_usolve(mna->sp_matrix->A_numeric->U, *x);
	cs_ipvec(mna->sp_matrix->A_symbolic->q, *x, temp_b, mna->dimension);

	memcpy(*x, temp_b, mna->dimension * sizeof(double));
	free(temp_b);
}

/* Solves the sparse mna system with LU factorization and store the result in vector x */
void solve_complex_sparse_lu(mna_system_t *mna, cs_complex_t *x) {
	cs_complex_t *temp_b = (cs_complex_t *)malloc(mna->dimension * sizeof(cs_complex_t));
	memcpy(temp_b, mna->sp_matrix->e_ac, mna->dimension * sizeof(cs_complex_t));

	mna->sp_matrix->G_ac_symbolic = cs_ci_sqr(2, mna->sp_matrix->G_ac, 0);
	mna->sp_matrix->G_ac_numeric  = cs_ci_lu(mna->sp_matrix->G_ac, mna->sp_matrix->G_ac_symbolic, 1);

	cs_ci_ipvec(mna->sp_matrix->G_ac_numeric->pinv, temp_b, x, mna->dimension);
	cs_ci_lsolve(mna->sp_matrix->G_ac_numeric->L, x);
	cs_ci_usolve(mna->sp_matrix->G_ac_numeric->U, x);
	cs_ci_ipvec(mna->sp_matrix->G_ac_symbolic->q, x, temp_b, mna->dimension);

	memcpy(x, temp_b, mna->dimension * sizeof(cs_complex_t));
	free(temp_b);

	/*
	 * Symbolic and Numeric factorizations have to be freed in every step
	 * because in the next step/factorization the pointers will be lost otherwise
	 */
	cs_ci_sfree(mna->sp_matrix->G_ac_symbolic);
	cs_ci_nfree(mna->sp_matrix->G_ac_numeric);
}

/* Solve the MNA system using cholesky decomposition and store the result in vector x */
void solve_cholesky(double **A, double *b, gsl_vector_view x, int dimension, bool is_decomp) {
	gsl_matrix_view view_A = gsl_matrix_view_array(A[0], dimension, dimension);
	gsl_vector_view view_b = gsl_vector_view_array(b, dimension);
	if (!is_decomp) {
		/* Cholesky decomposition A = LL^T*/
		gsl_linalg_cholesky_decomp(&view_A.matrix);
	}
	gsl_linalg_cholesky_solve(&view_A.matrix, &view_b.vector, &x.vector);
}

/* Solve the MNA system using complex cholesky decomposition and store the result in vector x */
void solve_complex_cholesky(gsl_matrix_complex *A, gsl_vector_complex *b, gsl_vector_complex *x, int dimension) {
	/* Cholesky decomposition A = LL^T*/
	gsl_linalg_complex_cholesky_decomp(A);
	gsl_linalg_complex_cholesky_solve(A, b, x);
}

/* Solves the sparse mna system with Cholesky factorization and store the result in vector x */
void solve_sparse_cholesky(mna_system_t *mna, cs *A, double **x) {
	double *temp_b = (double *)malloc(mna->dimension * sizeof(double));
	memcpy(temp_b, mna->b, mna->dimension * sizeof(double));

	if (!mna->is_decomp) {
		mna->sp_matrix->A_symbolic = cs_schol(1, A);
		mna->sp_matrix->A_numeric  = cs_chol(A, mna->sp_matrix->A_symbolic);
		if (mna->sp_matrix->A_numeric == NULL || mna->sp_matrix->A_symbolic == NULL) {
			fprintf(stderr, "\nCholesky method failed...non SPD matrix\n");
			exit(EXIT_FAILURE);
		}
		cs_spfree(A);
	}

	cs_ipvec(mna->sp_matrix->A_symbolic->pinv, temp_b, *x, mna->dimension);
	cs_lsolve(mna->sp_matrix->A_numeric->L, *x);
	cs_ltsolve(mna->sp_matrix->A_numeric->L, *x);
	cs_pvec(mna->sp_matrix->A_symbolic->pinv, *x, temp_b, mna->dimension);

	memcpy(*x, temp_b, mna->dimension * sizeof(double));
	free(temp_b);
}

/* Solves the sparse mna system with Cholesky factorization and store the result in vector x */
void solve_complex_sparse_cholesky(mna_system_t *mna, cs_complex_t *x) {
	cs_complex_t *temp_b = (cs_complex_t *)malloc(mna->dimension * sizeof(cs_complex_t));
	memcpy(temp_b, mna->sp_matrix->e_ac, mna->dimension * sizeof(cs_complex_t));

	mna->sp_matrix->G_ac_symbolic = cs_ci_schol(1, mna->sp_matrix->G_ac);
	mna->sp_matrix->G_ac_numeric  = cs_ci_chol(mna->sp_matrix->G_ac, mna->sp_matrix->G_ac_symbolic);

	if (mna->sp_matrix->G_ac_numeric == NULL || mna->sp_matrix->G_ac_symbolic == NULL) {
		fprintf(stderr, "\nCholesky method failed...non SPD matrix\n");
		exit(EXIT_FAILURE);
	}

	cs_ci_ipvec(mna->sp_matrix->G_ac_symbolic->pinv, temp_b, x, mna->dimension);
	cs_ci_lsolve(mna->sp_matrix->G_ac_numeric->L, x);
	cs_ci_ltsolve(mna->sp_matrix->G_ac_numeric->L, x);
	cs_ci_pvec(mna->sp_matrix->G_ac_symbolic->pinv, x, temp_b, mna->dimension);

	memcpy(x, temp_b, mna->dimension * sizeof(cs_complex_t));
	free(temp_b);

	/*
	 * Symbolic and Numeric factorizations have to be freed in every step
	 * because in the next step/factorization the pointers will be lost otherwise
	 */
	cs_ci_sfree(mna->sp_matrix->G_ac_symbolic);
	cs_ci_nfree(mna->sp_matrix->G_ac_numeric);
}

/* Print the MNA system */
void print_mna_system(mna_system_t *mna, options_t *options) {
	if (!options->SPARSE) {
		printf("\nMNA A array:\n\n");
		print_array(mna->matrix->A, mna->dimension);
		printf("MNA b vector:\n\n");
		print_vector(mna->b, mna->dimension);
		if (options->TRAN) {
			// printf("MNA init_response vector:\n\n");
			// print_vector(mna->resp->value, mna->dimension);
			printf("\nMNA G array:\n\n");
			print_array(mna->matrix->A, mna->dimension);
			if (!options->TR) {
				printf("\nMNA hC array:\n\n");
				print_array(mna->matrix->hC, mna->dimension);
			}
			printf("\nMNA aGhC array:\n\n");
			print_array(mna->matrix->aGhC, mna->dimension);
			if (options->TR) {
				printf("\nMNA sGhC array:\n\n");
				print_array(mna->matrix->sGhC, mna->dimension);
			}
		}
		if (options->AC) {
			printf("\nMNA G_ac array:\n\n");
			print_complex_array(mna->matrix->G_ac, mna->dimension);
			printf("\nMNA e_ac vector:\n\n");
			print_complex_vector(mna->matrix->e_ac, mna->dimension);
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

/* Print the complex array */
void print_complex_array(gsl_matrix_complex *A, int dimension) {
	int rows = dimension;
	int cols = dimension;
	gsl_complex z;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			z = gsl_matrix_complex_get(A, i, j);
			printf("% 8.4lf %s %-8.4lfi ", GSL_REAL(z), GSL_IMAG(z) >= 0 ? "+" : "-", ABS(GSL_IMAG(z)));
		}
		printf("\n");
	}
	printf("\n");
}

/* Print the vector */
void print_vector(double *b, int dimension) {
	for (int i = 0; i < dimension; i++) {
		printf("%lf\n", b[i]);
	}
    printf("\n");
}

/* Print the complex vector */
void print_complex_vector(gsl_vector_complex *b, int dimension) {
	gsl_complex z;
	for (int i = 0; i < dimension; i++) {
		z = gsl_vector_complex_get(b, i);
		printf("% lf %s %lfi\n", GSL_REAL(z), GSL_IMAG(z) >= 0.0 ? "+" : "-", ABS(GSL_IMAG(z)));
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

/* Copy cs_di *A into a newly allocated struct */
cs_di *_cs_di_copy(cs_di *A) {
    cs_di *C ;
    int n, triplet, nn, p, nz, *Ap, *Ai, *Cp, *Ci ;
    double *Ax, *Cx;
    if (!A || !A->x) return NULL;    /* return if A NULL or pattern-only */
    n = A->n; Ap = A->p; Ai = A->i; Ax = A->x;
    triplet = (A->nz >= 0);          /* true if A is a triplet matrix */
    nz = triplet ? A->nz : Ap[n];
    C = cs_di_spalloc(A->m, n, A->nzmax, 1, triplet);
    if (!C) return NULL;
    Cp = C->p; Ci = C->i; Cx = C->x;
    nn = triplet ? nz : (n + 1);
    for (p = 0; p < nz; p++) Ci[p] = Ai[p];
    for (p = 0; p < nn; p++) Cp[p] = Ap[p];
    for (p = 0; p < nz; p++) Cx[p] = Ax[p];
    if (triplet) C->nz = nz;
    return C;
}

/* Free all the memory allocated for the MNA system */
void free_mna_system(mna_system_t **mna, options_t *options) {
	if (options->SPARSE) {
		if (options->ITER) {
			/* This is only needed in case we have iterative solvers otherwise it is free'd inside direct methods */
			cs_di_spfree((*mna)->sp_matrix->A);
			if (options->AC) {
				cs_ci_spfree((*mna)->sp_matrix->G_ac);
			}
			if (options->TRAN) {
				cs_di_spfree((*mna)->sp_matrix->aGhC);
			}
		}
		else { /* Direct Methods */
			/*
			 * Notice that the equivalent symbolic,numeric for complex AC are being freed inside the solvers,
			 * because every time we're going to solve the AC MNA, we're gonna get new symbolic and numeric
			 * factorizations. Thus, it's required to free them inside the functions/solvers.
			 */
			cs_di_sfree((*mna)->sp_matrix->A_symbolic);
			cs_di_nfree((*mna)->sp_matrix->A_numeric);
			if (options->AC) {
				/* Free G_ac for the last step, this is freed by create_sparse_ac_mna function in the previous steps */
				cs_ci_spfree((*mna)->sp_matrix->G_ac);
			}
		}
		/* In case there is an TRAN analysis in the netlist */
		if (options->TRAN) {
			/* sGhC only exists in Trapezoidal method */
			if (options->TR) {
				cs_di_spfree((*mna)->sp_matrix->sGhC);
			}
			/* This is free'd in the creation of the MNA in case it is Trapezoidal */
			if (options->BE) {
				cs_di_spfree((*mna)->sp_matrix->hC);
			}
		}
		/* Free the A_base in case it exists */
		if (options->TRAN || options->AC) {
			cs_di_spfree((*mna)->sp_matrix->A_base);
		}
		/* Free the complex struct for the cs routines */
		if (options->AC) {
			free((*mna)->sp_matrix->e_ac);
		}
		/* Free the whole sp matrix struct */
		free((*mna)->sp_matrix);
	}
	else { /* Free everything from sparse */
		free((*mna)->matrix->A[0]);
		free((*mna)->matrix->A);
		gsl_permutation_free((*mna)->matrix->P);

		/* In case there is an TRAN analysis in the netlist */
		if (options->TRAN) {
			free((*mna)->matrix->aGhC[0]);
			free((*mna)->matrix->aGhC);
			/* sGhC only exists in Trapezoidal method */
			if (options->TR) {
				free((*mna)->matrix->sGhC[0]);
				free((*mna)->matrix->sGhC);
			}
			/* This is free'd in the creation of the MNA in case it is Trapezoidal */
			if (options->BE) {
				free((*mna)->matrix->hC[0]);
				free((*mna)->matrix->hC);
			}
		}

		/* In case there is an AC analysis in the netlist */
		if (options->AC) {
			/* Free all the allocated complex structures from GSL */
			gsl_matrix_complex_free((*mna)->matrix->G_ac);
			gsl_vector_complex_free((*mna)->matrix->e_ac);
		}
		
		/* Free the A_base in case it exists */
		if (options->TRAN || options->AC) {
			free((*mna)->matrix->A_base[0]);
			free((*mna)->matrix->A_base);
		}

		/* Free the matrix struct */
		free((*mna)->matrix);
	}

	/* Free the preconditioners for the iterative solvers */
	if (options->ITER) {
		free((*mna)->M);
		if (options->TRAN) {
			free((*mna)->M_trans);
		}
		if (options->AC) {
			gsl_vector_complex_free((*mna)->M_ac);
			gsl_vector_complex_free((*mna)->M_ac_conj);
		}
	}
	/* Free b*/
	free((*mna)->b);

	/* Free every string allocated for the group2 elements */
	for (int i = 0; i < (*mna)->num_g2_elem; i++) {
		free((*mna)->g2_indx[i].element);
	}
	/* Free the g2 array */
	free((*mna)->g2_indx);

	/* Free the resp struct */
	if (options->TRAN) {
		free((*mna)->resp->value);
		free((*mna)->resp->nodes);
		free((*mna)->resp);
	}

	/* Free the whole MNA struct */
	free(*mna);
	/* Set mna to NULL to limit further acesses */
	*mna = NULL;
}