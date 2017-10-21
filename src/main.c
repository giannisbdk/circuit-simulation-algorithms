#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parser.h"

#include "list.h"
#include "hash_table.h"
#include "mna_dc.h"

#define HASH_TABLE_SIZE 65536

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

    mna_arrays_t mna;

    int num_nodes = hash_table->seq - 1;
    int size = num_nodes + num_branches;
    printf("\nsize: %d\nnum_nodes(w/o ground): %d\nnum_branches_g2: %d\n\n", size, num_nodes, num_branches);
    
    // Initialize the arrays
    mna.left = init_array(size, size);
    mna.right = init_array(size, 1);

    create_mna_arrays(&mna, index, hash_table, num_nodes);

    print_mna_left(&mna, size);
    print_mna_right(&mna, size);

    fclose(file_input);
    if (line) {
        free(line);
    }
    
    /* Free the hashtable memory */
    ht_free(hash_table);

    exit(EXIT_SUCCESS);

	return 0;
}