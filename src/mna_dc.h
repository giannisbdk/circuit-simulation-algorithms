#ifndef MNA_DC_H
#define MNA_DC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "list.h"

typedef struct mna_system {
	gsl_matrix *A;
	gsl_vector *b;
	int dimension;
} mna_system_t;

void create_mna_system(mna_system_t *mna, index_t *index, hash_table_t *hash_table, int offset);
mna_system_t *init_mna_system(int dimension);
gsl_matrix *init_array(int row, int col);
gsl_vector *init_vector(int row);
void print_mna_system(mna_system_t *mna);
void print_array(gsl_matrix *A);
void print_vector(gsl_vector *b);

#endif