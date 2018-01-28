#ifndef MNA_H
#define MNA_H

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <stdbool.h>

#include "list.h"
#include "parser.h"
#include "iter.h"
#include "routines.h"
#include "../cx_sparse/Include/cs.h"

/* Holds the transient response and the nodes that contribute to it */
typedef struct resp {
	double  *value;
	list1_t **nodes;
} resp_t;

/* Holds the dense representation of the MNA */
typedef struct matrix {
	/* Main Matrix of the MNA for DC */
	double **A;
	/* A_base is A before factorization */
	double **A_base;

	/* Permutation for LU factorizations either complex or normal */
	gsl_permutation *P;

	/* Matrices for the Transient Analysis */
	double **hC;
	/* aGhC = G + hC, G = A */
	double **aGhC;
	/* aGhC = G - hC, G = A */
	double **sGhC;

	/* Complex matrix G for the AC Analysis and RHS vector */
	gsl_matrix_complex *G_ac;
	gsl_vector_complex *e_ac;
} matrix_t;

/* Holds the sparse representation of the MNA */
typedef struct sp_matrix {
	/* Main Matrix of the MNA for DC */
	cs *A;
	/* The A matrix before any modifications to use in AC analysis */
	cs *A_base;

	/* Matrices for the Transient Analysis */
	cs *hC;
	/* aGhC = G + hC, G = A */
	cs *aGhC;
	/* aGhC = G - hC, G = A */
	cs *sGhC;

	/* Hold the symbolic and numeric representation of the LU factorization */
	css *A_symbolic;
	csn *A_numeric;

	/* Sparse data structures for AC analysis */
	cs_ci *G_ac;

	/* This is necessary for the sparse routines, instead of using gsl complex */
	cs_complex_t *e_ac;

	/* Hold the symbolic and numeric representation of the LU factorization */
	cs_cis *G_ac_symbolic;
	cs_cin *G_ac_numeric;
} sp_matrix_t;

/* Keeps the indexing for the sources of group 2 */
typedef struct g2_indx {
	char *element;
} g2_indx_t;

typedef struct mna_system {
	/* Pointers to the dense and sparse matrix of the MNA */
	matrix_t 	*matrix;
	sp_matrix_t *sp_matrix;

	/* Pointer to the transient response and the nodes that contribute to it */
	resp_t *resp;

	/* Jacobi Preconditioner M, we only store the diagonal of A not zeros */
	/* For DC and TRANSIENT analysis */
	double *M;
	double *M_trans;
	/* For AC analysis */
	gsl_vector_complex *M_ac;
	gsl_vector_complex *M_ac_conj;

	/* Right hand side vector b for the Ax=b */
	double *b;

	/* General info about MNA */
	bool is_decomp;
	int dimension;
	int num_nodes;

	/* Keep some info about g2 elements */
	int num_g2_elem;
	g2_indx_t *g2_indx;

	//TODO change the below flag names or create a new enum called current state? or something
	/* Flags that indicate in which analysis we're currently at */
	bool tr_analysis_init;
	bool ac_analysis_init;
} mna_system_t;

mna_system_t *init_mna_system(int num_nodes, int num_g2_elem, options_t *options, int nz);
void init_sparse_matrix(mna_system_t *mna, options_t *options, int nz);
void init_dense_matrix(mna_system_t *mna, options_t *options);
void create_mna_system(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, double tr_step, int offset);
void create_ac_mna_system(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, int offset, double omega);
void create_sparse_mna(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, int offset);
void create_dense_mna(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, int offset);
void create_dense_trans_mna(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, double tr_step, int offset);
void create_dense_ac_mna(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, int offset, double omega);
void create_sparse_ac_mna(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, int offset, double omega);
void create_sparse_trans_mna(mna_system_t *mna, index_t *index, hash_table_t *hash_table, options_t *options, double tr_step, int offset);
void solve_mna_system(mna_system_t *mna, double **x, gsl_vector_complex *x_complex, options_t *options);
void solve_lu(double **A, double *b, gsl_vector_view x, gsl_permutation *P, int dimension, bool is_decomp);
void solve_complex_lu(gsl_matrix_complex *A, gsl_vector_complex *b, gsl_vector_complex *x, gsl_permutation *P, int dimension);
void solve_sparse_lu(mna_system_t *mna, cs *A, double **x);
void solve_complex_sparse_lu(mna_system_t *mna, cs_complex_t *x);
void solve_cholesky(double **A, double *b, gsl_vector_view x, int dimension, bool is_decomp);
void solve_complex_cholesky(gsl_matrix_complex *A, gsl_vector_complex *b, gsl_vector_complex *x, int dimension);
void solve_sparse_cholesky(mna_system_t *mna, cs *A, double **x);
void solve_complex_sparse_cholesky(mna_system_t *mna, cs_complex_t *x);
int g2_elem_indx(g2_indx_t *g2_indx, int num_nodes, int num_g2_elem, char *element);
double get_response_value(list1_t *curr);
void print_mna_system(mna_system_t *mna, options_t *options);
void print_array(double **A, int dimension);
void print_complex_array(gsl_matrix_complex*A, int dimension);
void print_vector(double *b, int dimension);
void print_complex_vector(gsl_vector_complex *b, int dimension);
void print_permutation(gsl_permutation *P);
cs_di *_cs_di_copy (cs_di *A);
void free_mna_system(mna_system_t **mna, options_t *options);

#endif