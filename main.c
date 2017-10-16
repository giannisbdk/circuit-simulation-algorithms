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

    printf("%s\n", argv[1]);

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
        if(tokens[1][0] == 'V' || tokens[1][0] == 'v') {
            num_branches++;
        }
        /* Free all the memory we allocated */
        for (int i = 0; i < num_tokens; i++) {
            free(tokens[i]);
        }
        free(tokens);
    }

    printf("Printing the lists\n");
    print_lists(index, hash_table);

    for(int i=0; i<=hash_table->size; i++) {
        if(ht_get_value(hash_table, i) != NULL)
            printf("Key: %d, Value: %s\n", i, ht_get_value(hash_table, i));
    }

    mna_arrays_t mna;

    int num_nodes = hash_table->size - 1;
    int size = num_nodes + num_branches;
    printf("size: %d, num_nodes: %d, num_branches: %d\n", size, num_nodes, num_branches);
    
    // Initialize the arrays
    mna.left = init_array(size, size);
    mna.right = init_array(size, 1);

    create_mna_arrays(&mna, index, num_nodes);

    printf("\n");
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            printf("%lf\t", mna.left[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < 1; ++j) {
            printf("%lf\t", mna.right[i][j]);
        }
        printf("\n");
    }

    fclose(file_input);
    if (line) {
        free(line);
    }

    exit(EXIT_SUCCESS);

	return 0;
}