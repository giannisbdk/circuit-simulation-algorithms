#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parser.h"
#include "list.h"
#include "hash_table.h"

#define HASH_TABLE_SIZE 65536

int main(int argc, char *argv[]) {

	FILE *file_input;
    char *line = NULL;
    int num_tokens = 0;
    size_t len = 0;
    ssize_t read;
    char **tokens;

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
        /* Free all the memory we allocated */
        for (int i = 0; i < num_tokens; i++) {
            free(tokens[i]);
        }
        free(tokens);
    }

    printf("Printing the lists\n");
    print_lists(index, hash_table);

    fclose(file_input);
    if (line) {
        free(line);
    }

    exit(EXIT_SUCCESS);

	return 0;
}