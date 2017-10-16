#ifndef MNA_DC_H
#define MNA_DC_H

#include "list.h"

typedef struct mna_arrays {
	double **left;
	double **right;
} mna_arrays_t;

void create_mna_arrays(mna_arrays_t *, index_t *, int);
double **init_array(int, int);

#endif