#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "parser.h"
#include "list.h"
#include "hash_table.h"
#include "mna.h"
#include "routines.h"
#include "dc_analysis.h"
#include "transient_analysis.h"
#include "ac_analysis.h"
#include "time_tools.h"

/* This is defined here because we specify the length from main */
#define HASH_TABLE_SIZE 59999

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "You must specify the input file (netlist) from cmd arguments.\nExiting....\n");
        exit(EXIT_FAILURE);
    }
    /* Initialize our double linked lists */
    index_t *index = init_lists();

    /* Initialize the parser */
    parser_t *parser = init_parser();

    /* Create the hashtable */
    hash_table_t *hash_table = ht_create(HASH_TABLE_SIZE);

    /* Parse the netlist filename is in argv[1] */
    parse_netlist(parser, argv[1], index, hash_table);

    /* Initialize the MNA_system */
    mna_system_t *mna = init_mna_system(parser->netlist->num_nodes, parser->netlist->num_g2_elem, parser->options, parser->netlist->nz);
    
    /* Create the MNA system */
    create_mna_system(mna, index, hash_table, parser->options, parser->tr_analysis->time_step, parser->netlist->num_nodes);
    
    /* Print the MNA system */
    /* Dimension of MNA system */
    int dimension = parser->netlist->num_nodes + parser->netlist->num_g2_elem;

    /* Sol_x will hold the solution of the MNA system */
    double *sol_x = (double *)calloc(dimension, sizeof(double));
    double *dc_op = (double *)calloc(dimension, sizeof(double));

    /* Solve the MNA system */
    solve_mna_system(mna, &sol_x, NULL, parser->options);

    /* DC Operating Point to file */
    dc_operating_point(hash_table, sol_x);
    
    /* Hold DC operating point values for transient and ac analyses */
    memcpy(dc_op, sol_x, mna->dimension * sizeof(double));

    /* DC Sweep analysis to file */
    dc_analysis(index->head1, hash_table, mna, parser, sol_x);

    /* Transient analysis to file */
    tr_analysis(hash_table, mna, parser, dc_op, sol_x);

    /* AC analysis to file */
    gsl_vector_complex *x_complex = init_gsl_complex_vector(mna->dimension);
    ac_analysis(index, hash_table, mna, parser, dc_op, x_complex);

    /* Free all the dynamic allocated memory */
    free_index(&index);
    free_mna_system(&mna, parser->options);
    free_parser(&parser);
    ht_free(&hash_table);
    gsl_vector_complex_free(x_complex);
    free(sol_x);
    free(dc_op);

	return 0;
}