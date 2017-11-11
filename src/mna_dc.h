#ifndef MNA_DC_H
#define MNA_DC_H

#include <gsl/gsl_linalg.h>
#include <stdbool.h>

#include "list.h"

/* Keeps the indexing for the sources of group 2 */
typedef struct g2_indx {
	char *element;
} g2_indx_t;

typedef struct mna_system {
	double **A;
	double *b;
	gsl_permutation *P;
	bool is_decomp;
	int num_nodes;
	/* Keep some info about g2 elements */
	int num_g2_elem;
	g2_indx_t *g2_indx;
} mna_system_t;

mna_system_t *init_mna_system(int num_nodes, int num_g2_elem);
void create_mna_system(mna_system_t *mna, index_t *index, hash_table_t *hash_table, int offset);
int g2_elem_indx(g2_indx_t *g2_indx, int num_nodes, int num_g2_elem, char *element);
gsl_vector *solve_mna_system(mna_system_t *mna, bool SPD);
gsl_vector *solve_lu(mna_system_t *mna);
gsl_vector *solve_cholesky(mna_system_t *mna);
double **init_array(int row, int col);
double *init_vector(int row);
gsl_permutation *init_permutation(int dimension);
void print_mna_system(mna_system_t *mna);
void print_array(double **A, int dimension);
void print_vector(double *b, int dimension);
void print_permutation(gsl_permutation *P);
void free_mna_system(mna_system_t **mna);

#endif