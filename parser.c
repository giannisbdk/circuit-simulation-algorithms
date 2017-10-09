#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum {
	RESISTOR,
	CURRENT_SOURCE,
	VOLTAGE_SOURCE,
	CAPACITOR,
	INDUCTOR,
	DIODES,
	MOS,
	BJT
} TYPE;

struct list {
	int id;
	int *probes;
	float value;
	struct list *nxt;
	TYPE type;
};

typedef struct list list_t;

int get_num_tokens(char *line) {

	int tokens = 0;
	int len = strlen(line);
	line[0] = (line[0] == '\t') ? ' ' : line[0];
	char found_char = (line[0] == ' ' ? 0 : 1);

	if(found_char == 1) {
		tokens++;
	}

	for (int i = 1; i < len; ++i) {
		if(line[i] == '\t') {
			line[i] = ' ';
		}
		if(line[i] == ' ') {
			found_char = 0;
		}
		else if(line[i] > 32) {
			if(found_char == 0) {
				tokens++;
				found_char = 1;
			}
		}
	}
	return tokens;
}

char *remove_comments(char *line) {

	char *ret, *token;

	if (line[0] == '*') {
		return NULL;
	}

	ret = strchr(line, '*');

	if(ret != NULL) {
		token = strtok(line, "*");
		ret = token;
	}
	else {
		ret = line;
	}

	return ret;
}

char **tokenizer(char *line, int *num_tokens) {

	const char delim[2] = " ";
	char **tokens, *token;
	int i = 0;

	line = remove_comments(line);

	if(line == NULL) {
		return NULL;
	}

	*num_tokens = get_num_tokens(line);

	if(*num_tokens == 0) {
		return NULL;
	}

	tokens = (char **) malloc((*num_tokens) * sizeof(char *));

	token = strtok(line, delim);

	while(token != NULL && i < *num_tokens) {
		tokens[i] = (char *) malloc(strlen(token) * sizeof(char));
		strcpy(tokens[i], token);
		token = strtok(NULL, delim);
		i++;
	}

	/* Trim '\n' if necessary */
	if(tokens[i-1][strlen(tokens[i-1])-1] == '\n') {
		tokens[i-1][strlen(tokens[i-1])-1] = '\0';
	}
	return tokens;
}

void save_node(list_t *head, char **tokens) {

	list_t *curr;

	curr = (list_t *) malloc(sizeof(list_t));

	while(curr->nxt != NULL) {

	}

}

int main(int argc, char *argv[]) {

	FILE *file_input;
    char *line = NULL;
    int num_tokens = 0;
    size_t len = 0;
    ssize_t read;
    list_t *head = NULL;
    char **tokens;

    printf("%s\n", argv[1]);

    file_input = fopen(argv[1], "rb");
    if (file_input == NULL)
        exit(EXIT_FAILURE);

    while((read = getline(&line, &len, file_input)) != -1) {
    	// save_node(head, tokenizer(line, &num_tokens));
    	tokens = tokenizer(line, &num_tokens);
    	if(tokens == NULL) {
    		continue;
    	}
    	for (int i = 0; i < num_tokens; i++) {
    		printf("%s\t", tokens[i]);
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