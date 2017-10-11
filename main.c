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
    if (file_input == NULL)
        exit(EXIT_FAILURE);

    list_uno_t *head1 = init_list1();
    list_dos_t *head2 = init_list2();

    while((read = getline(&line, &len, file_input)) != -1) {

    	tokens = tokenizer(line, &num_tokens);
    	if(tokens == NULL) {
    		continue;
    	}
        if(add_to_list(head1, head2, tokens) == -1) {
            exit(EXIT_FAILURE);
        }
    }

    printf("Printing the lists\n");
    print_lists(head1, head2);

    fclose(file_input);

    if (line) {
        free(line);
    }

    exit(EXIT_SUCCESS);

	return 0;
}