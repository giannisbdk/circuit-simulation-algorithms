#include "mna_dc.h"

mna_system_t *init_mna_system(int dimension) {
	mna_system_t *mna = (mna_system_t *)malloc(sizeof(mna_system_t));
	assert(mna != NULL);
	mna->dimension = dimension;
	mna->A = init_array(dimension, dimension);
	mna->b = init_vector(dimension);
	return mna;
}

gsl_matrix *init_array(int row, int col) {
    gsl_matrix *A;
    A = gsl_matrix_calloc(row, col);
    return A;
}

gsl_vector *init_vector(int row) {
	gsl_vector *b;
	b = gsl_vector_calloc(row);
	return b;
}

void create_mna_system(mna_system_t *mna, index_t *index, hash_table_t *hash_table, int offset) {

	list1_t *curr;
	double value;
	int volt_sources_cnt = 0; 

	for (curr = index->head1; curr != NULL; curr = curr->next) {

		int probe1_id = ht_get_id(hash_table, curr->probe1);
		int probe2_id = ht_get_id(hash_table, curr->probe2);
		int i = probe1_id - 1;
		int j = probe2_id - 1;

		if(curr->type == 'C' || curr->type == 'c') {
			continue;
		}
		else if(curr->type == 'R' || curr->type == 'r') {
			value = 1 / curr->value;
			if(probe1_id == 0) {
				gsl_matrix_set(mna->A, j, j, gsl_matrix_get(mna->A, j, j) + value);
			}
			else if(probe2_id == 0) {
				gsl_matrix_set(mna->A, i, i, gsl_matrix_get(mna->A, i, i) + value);
			}
			else {
				gsl_matrix_set(mna->A, i, i, gsl_matrix_get(mna->A, i, i) + value);
				gsl_matrix_set(mna->A, j, j, gsl_matrix_get(mna->A, j, j) + value);
				gsl_matrix_set(mna->A, i, j, gsl_matrix_get(mna->A, i, j) - value);
				gsl_matrix_set(mna->A, j, i, gsl_matrix_get(mna->A, j, i) - value);
			}
		}
		else if(curr->type == 'I' || curr->type == 'i') {
			value = curr->value;
			if(probe1_id == 0) {
				gsl_vector_set(mna->b, j, gsl_vector_get(mna->b, j) + value);
			}
			else if(probe2_id == 0) {
				gsl_vector_set(mna->b, i, gsl_vector_get(mna->b, i) + value);
			}
			else {
				gsl_vector_set(mna->b, i, gsl_vector_get(mna->b, i) + value);
				gsl_vector_set(mna->b, j, gsl_vector_get(mna->b, j) + value);
			}
		}
		else {
			if(curr->type == 'L' || curr->type == 'l') {
				/* We treat the Inductor like a voltage source with value 0 */
				value = 0.0;
			}
			else if(curr->type == 'V' || curr->type == 'v') {
				value = curr->value;
			}
			if(probe1_id == 0) {
				gsl_matrix_set(mna->A, j, offset + volt_sources_cnt, 1.0);
				gsl_matrix_set(mna->A, offset + volt_sources_cnt, j, 1.0);
				gsl_vector_set(mna->b, offset + volt_sources_cnt, gsl_vector_get(mna->b, offset + volt_sources_cnt) + value);
			}
			else if(probe2_id == 0) {
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

void print_mna_system(mna_system_t *mna) {
	print_array(mna->A);
	print_vector(mna->b);
}

void print_array(gsl_matrix *A) {
	printf("MNA A array:\n\n");
	gsl_matrix_fprintf(stdout, A, "%lf");
	printf("\n");
}

void print_vector(gsl_vector *b) {
	printf("MNA b vector:\n\n");
	gsl_vector_fprintf(stdout, b, "%lf");
    printf("\n");
}