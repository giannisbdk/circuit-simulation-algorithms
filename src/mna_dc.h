#ifndef MNA_DC_H
#define MNA_DC_H

#include <gsl/gsl_linalg.h>
#include <stdbool.h>

#include "list.h"
#include "parser.h"
#include "iter.h"
#include "routines.h"
#include "csparse.h"

/* Holds the transient response and the nodes that contribute to it */
typedef struct resp {
	double  *value;
	list1_t **nodes;
} resp_t;

/* Holds the dense representation of the MNA */
typedef struct matrix {
	double **A;
	gsl_permutation *P;
	double **G;
	double **hC;
	/* aGhC = G + hC */
	double **aGhC;
	/* aGhC = G - hC */
	double **sGhC;
	resp_t *resp;
} matrix_t;

/* Holds the sparse representation of the MNA */
typedef struct sp_matrix {
	cs 	*A;
	css *A_symbolic;
	csn *A_numeric;
} sp_matrix_t;

/* Keeps the indexing for the sources of group 2 */
typedef struct g2_indx {
	char *element;
} g2_indx_t;

typedef struct mna_system {
	/* Pointers to the dense and sparse matrix of the MNA */
	matrix_t 	*matrix;
	sp_matrix_t *sp_matrix;
	/* Jacobi Preconditioner M, we only store the diagonal of A not zeros */
	double *M;
	/* Right hand side vector b for the Ax=b */
	double *b;
	bool is_decomp;
	int dimension;
	int num_nodes;
	/* Keep some info about g2 elements */
	int num_g2_elem;
	g2_indx_t *g2_indx;
	bool tr_analysis_init;
} mna_system_t;

mna_system_t *init_mna_system(int num_nodes, int num_g2_elem, options_t *options, int nz);
void init_sparse_matrix(mna_system_t *mna, options_t *options, int nz);
void init_dense_matrix(mna_system_t *mna, options_t *options);
void create_mna_system(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, double tr_step, int offset);
void create_sparse_mna(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, int offset);
void create_dense_mna(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, int offset);
void create_dense_trans_mna(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, double tr_step, int offset);
void solve_mna_system(mna_system_t *mna, double **x, options_t *options);
// void solve_lu(mna_system_t *mna, gsl_vector_view x);
void solve_lu(double **A, double *b, gsl_vector_view x, gsl_permutation *P, int dimension, bool is_decomp);
// void solve_cholesky(mna_system_t *mna, gsl_vector_view x);
void solve_cholesky(double **A, double *b, gsl_vector_view x, int dimension, bool is_decomp);
void solve_sparse_lu(mna_system_t *mna, double **x);
void solve_sparse_cholesky(mna_system_t *mna, double **x);
int g2_elem_indx(g2_indx_t *g2_indx, int num_nodes, int num_g2_elem, char *element);
double **init_array(int row, int col);
double *init_vector(int row);
gsl_permutation *init_permutation(int dimension);
double get_response_value(list1_t *curr);
void print_mna_system(mna_system_t *mna, options_t *options);
void print_array(double **A, int dimension);
void print_vector(double *b, int dimension);
void print_permutation(gsl_permutation *P);
void free_mna_system(mna_system_t **mna, options_t *options);

#endif