#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parser.h"
#include "list.h"

int main(int argc, char *argv[]) {

	FILE *file_input;
    char *line = NULL;
    int num_tokens = 0;
    size_t len = 0;
    ssize_t read;
    char **tokens;

    printf("%s\n", argv[1]);

    file_input = fopen(argv[1], "rb");
    if (file_input == NULL) {
        exit(EXIT_FAILURE);
    }

    index_t *index = init_lists();

    while((read = getline(&line, &len, file_input)) != -1) {

    	tokens = tokenizer(line, &num_tokens);
    	if(tokens == NULL) {
    		continue;
    	}
        if(add_to_list(index, tokens) == 0) {
            exit(EXIT_FAILURE);
        }
    }

    printf("Printing the lists\n");
    print_lists(index);

    fclose(file_input);

    if (line) {
        free(line);
    }

    exit(EXIT_SUCCESS);

	return 0;
}