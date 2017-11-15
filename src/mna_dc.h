#ifndef MNA_DC_H
#define MNA_DC_H

#include <gsl/gsl_linalg.h>
#include <stdbool.h>

#include "list.h"
#include "parser.h"
#include "iter.h"
#include "routines.h"

/* The Matrices,Vectors below should keep their values among different invocations of iter solvers */
typedef struct matrix {
	/* Jacobi Preconditioner M and M transpose (vectors), we only store the diagonal of A not zeros */
	double *M;
	double *M_trans;
	/* MNA A and A transpose matrices */
	double **A;
	double **A_trans;
	/* Right hand side vector b for the Ax=b and Permutation vector */
	double *b;
	gsl_permutation *P;
} matrix_t;

/* Keeps the indexing for the sources of group 2 */
typedef struct g2_indx {
	char *element;
} g2_indx_t;

typedef struct mna_system {
	/* Mat is a pointer to a struct that contains all the matrices of the MNA system */
	matrix_t *matrix;
	bool is_decomp;
	int dimension;
	int num_nodes;
	/* Keep some info about g2 elements */
	int num_g2_elem;
	g2_indx_t *g2_indx;
} mna_system_t;

mna_system_t *init_mna_system(int num_nodes, int num_g2_elem, options_t *options);
void create_mna_system(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, int offset);
void solve_mna_system(mna_system_t *mna, double **x, options_t *options);
void solve_lu(mna_system_t *mna, gsl_vector_view x);
void solve_cholesky(mna_system_t *mna, gsl_vector_view x);
int g2_elem_indx(g2_indx_t *g2_indx, int num_nodes, int num_g2_elem, char *element);
double **init_array(int row, int col);
double *init_vector(int row);
gsl_permutation *init_permutation(int dimension);
void print_mna_system(mna_system_t *mna);
void print_array(double **A, int dimension);
void print_vector(double *b, int dimension);
void print_permutation(gsl_permutation *P);
void free_mna_system(mna_system_t **mna, options_t *options);

#endif