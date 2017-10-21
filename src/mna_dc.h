#ifndef MNA_DC_H
#define MNA_DC_H

#include "list.h"

typedef struct mna_arrays {
	double **left;
	double **right;
} mna_arrays_t;

void create_mna_arrays(mna_arrays_t *mna, index_t *index, hash_table_t *hash_table, int offset);
double **init_array(int row, int col);
void print_mna_left(mna_arrays_t *mna, int size);
void print_mna_right(mna_arrays_t *mna, int size);

#endif