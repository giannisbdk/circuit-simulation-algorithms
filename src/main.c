#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include "parser.h"
#include "list.h"
#include "hash_table.h"
#include "mna_dc.h"
#include "routines.h"
#include "dc_analysis.h"

int errno;

int main(int argc, char *argv[]) {
	FILE *file_input;

    if (argc < 2) {
        printf("You must specify input file from cmd arguments.\nExiting....\n");
        exit(EXIT_FAILURE);
    }
    printf("\nInput file is: %s\n", argv[1]);
    file_input = fopen(argv[1], "rb");
    if (file_input == NULL) {
        fprintf(stderr, "Error: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    index_t *index = init_lists();
    hash_table_t *hash_table = ht_create(HASH_TABLE_SIZE);
    /* Parse the netlist */
    parser_t *parser = parse_netlist(file_input, index, hash_table);
    fclose(file_input);
    /* Initialize the MNA_system */
    mna_system_t *mna = init_mna_system(parser->netlist_elem->num_nodes, parser->netlist_elem->num_g2_elem, parser->options);
    create_mna_system(mna, index, hash_table, parser->options, parser->netlist_elem->num_nodes);
    print_mna_system(mna);

    int size = parser->netlist_elem->num_nodes + parser->netlist_elem->num_g2_elem;
    double *sol_x = (double *)calloc(size, sizeof(double));
    solve_mna_system(mna, &sol_x, parser->options);
    printf("Solution of the MNA system:\n\n");
    print_vector(sol_x, size);

    /* DC Operating Point to file */
    dc_operating_point(hash_table, sol_x);
    /* DC Sweep to file */
    dc_sweep(index->head1, hash_table, mna, parser->dc_analysis, parser->options, parser->netlist_elem, sol_x);

    //TODO Need to free the dc_analysis mallocs we did before
    // TODO Need to free the parser 
    /* Free all the dynamic allocated memory */
    free_index(&index);
    free_mna_system(&mna, parser->options);
    ht_free(&hash_table);
	return 0;
}