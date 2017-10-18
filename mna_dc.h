#ifndef MNA_DC_H
#define MNA_DC_H

#include "list.h"

typedef struct mna_arrays {
	double **left;
	double **right;
} mna_arrays_t;

void create_mna_arrays(mna_arrays_t *, index_t *, hash_table_t *, int);
double **init_array(int, int);
void print_mna_left(mna_arrays_t *, int);
void print_mna_right(mna_arrays_t *, int);

#endif