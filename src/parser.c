#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parser.h"

int get_num_tokens(char *line) {

	int tokens = 0;
	int len = strlen(line);
	line[0] = (line[0] == '\t') ? ' ' : line[0];
	char found_char = (line[0] == ' ' ? 0 : 1);

	if(found_char == 1) {
		tokens++;
	}
	for (int i = 1; i < len; ++i) {
		if(line[i] == ' ' || line[i] == '\t') {
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

/* Tokenizes the current line and returns the number of tokens and the tokens */
char **tokenizer(char *line, int *num_tokens) {

	const char delim[] = " \t";
	const char CRLF[] = "\r\n";
	char **tokens, *token;
	short int i = 1;

	/* In case the line starts with a comment or it's an empty line ignore */
	if(line[0] == '*' || line[0] == '\n' || strcmp(line, CRLF) == 0) {
		return NULL;
	}

	*num_tokens = get_num_tokens(line);
	if(*num_tokens == 0) {
		return NULL;
	}

	tokens = (char **)malloc((*num_tokens + 1) * sizeof(char *));
	tokens[0] = (char *)malloc(sizeof(char));
	sprintf(tokens[0], "%d", *num_tokens);

	token = strtok(line, delim);

	while(token != NULL && i <= *num_tokens) {
		tokens[i] = (char *)malloc(strlen(token) * sizeof(char));
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

/* Initializes all options to the default values = false */
void init_options(options_t *options) {
	options->SPD 	= false;
	options->ITER 	= false;
	options->SPARSE = false;
}