#ifndef MNA_DC_H
#define MNA_DC_H

#include <gsl/gsl_linalg.h>
#include <stdbool.h>

#include "list.h"

typedef struct mna_system {
	gsl_matrix *A;
	gsl_vector *b;
	gsl_permutation *P;
	int dimension;
	bool is_decomp;
} mna_system_t;

void create_mna_system(mna_system_t *mna, index_t *index, hash_table_t *hash_table, int offset);
mna_system_t *init_mna_system(int dimension);
gsl_vector *solve_mna_system(mna_system_t *mna, bool SPD);
gsl_vector *solve_lu(mna_system_t *mna);
gsl_vector *solve_cholesky(mna_system_t *mna);
gsl_matrix *init_array(int row, int col);
gsl_vector *init_vector(int row);
gsl_permutation *init_permutation(int dimension);
void print_mna_system(mna_system_t *mna);
void print_array(gsl_matrix *A);
void print_vector(gsl_vector *b);
void free_mna_system(mna_system_t *mna);

#endif