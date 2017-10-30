#ifndef PARSER_H
#define PARSER_H

#include <stdbool.h>

/* Struct to hold the different options for analysis */
typedef struct options {
	bool SPD;
	bool SPARSE;
	bool ITER;
} options_t;

/* Struct to hold the different DC analysis options */
typedef struct dc_analysis {
	char *volt_source;
	double start;
	double end;
	double increment;
	/* Plot / Print */
	char **nodes;
	int num_nodes;
} dc_analysis_t;

int get_num_tokens(char *line);
char **tokenizer(char *line, int *num_tokens);
void init_options(options_t *options);

#endif