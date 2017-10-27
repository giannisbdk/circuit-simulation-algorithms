#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parser.h"
#include <errno.h>
#include "list.h"
#include "hash_table.h"
#include "mna_dc.h"

#define HASH_TABLE_SIZE 65536

int errno;

int main(int argc, char *argv[]) {

	FILE *file_input;
    char *line = NULL;
    int num_tokens = 0;
    size_t len = 0;
    ssize_t read;
    char **tokens;
    int num_branches = 0;

    index_t *index = init_lists();
    hash_table_t *hash_table = ht_create(HASH_TABLE_SIZE);

    if(argc < 2) {
        printf("You must specify input file from cmd arguments.\nExiting....\n");
        exit(EXIT_FAILURE);
    }

    printf("Input file is: %s\n", argv[1]);

    file_input = fopen(argv[1], "rb");
    if (file_input == NULL) {
        exit(EXIT_FAILURE);
    }

    while((read = getline(&line, &len, file_input)) != -1) {

    	tokens = tokenizer(line, &num_tokens);
    	if(tokens == NULL) {
    		continue;
    	}
        if(add_to_list(index, tokens, hash_table) == FAILURE) {
            exit(EXIT_FAILURE);
        }
        if(tokens[1][0] == 'V' || tokens[1][0] == 'v' || tokens[1][0] == 'L' || tokens[1][0] == 'l') {
            num_branches++;
        }
        /* Free all the memory we allocated */
        for (int i = 0; i < num_tokens; i++) {
            free(tokens[i]);
        }
        free(tokens);
    }
#ifdef DEBUGL
    printf("Printing the lists\n");
    print_lists(index, hash_table);
#endif
    printf("Finished parsing %d circuit elements.\n", index->size1 + index->size2);

    mna_system_t *mna;

    int num_nodes = hash_table->seq - 1;
    int size = num_nodes + num_branches;
    printf("\nsize: %d\nnum_nodes(w/o ground): %d\nnum_branches_g2: %d\n\n", size, num_nodes, num_branches);
    
    /* Initialize the MNA_system */
    mna = init_mna_system(size);
    int SPD = 0;
    create_mna_system(mna, index, hash_table, num_nodes);
    print_mna_system(mna);
    gsl_vector *sol_x = solve_mna_system(mna, SPD);
    printf("Solution of the MNA system:\n\n");
    print_vector(sol_x);

    FILE *file_out = fopen("DC_Opearting_Point.txt", "w");
    if (file_out == NULL) {
        fprintf(stderr, "Error opening file: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
    fprintf(file_out, "%-15s%-15s\n", "Node", "Value");
    entry_t *curr;
    double value;
    int id;
    for (int i = 0; i < hash_table->size; i++) {
        for (curr = hash_table->table[i]; curr != NULL; curr = curr->next) {
            /* Get node name */
            id = curr->id;
            if (id == 0) {
                continue;
            }
            else {
                id -= 1;
            }
            /* Get the corresponding cell of the solution vector */
            value = gsl_vector_get(sol_x, id);
            /* Output to the file */
            fprintf(file_out, "%-15s%-15lf\n", curr->key, value);
        }
    }
    fclose(file_out);
    fclose(file_input);
    if (line) {
        free(line);
    }
    //TODO free the lists
    /* Free all the dynamic allocated memory */
    free_mna_system(mna);
    ht_free(hash_table);
	return 0;
}