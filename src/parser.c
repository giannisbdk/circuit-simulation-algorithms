#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>

#include "parser.h"

/* Initializes the parser */
parser_t *init_parser() {
    /* Allocate memory for the parser struct */
    parser_t *parser = (parser_t *)malloc(sizeof(parser_t));
    assert(parser != NULL);

    /* Initializes all options to the default values */
    parser->options = (options_t *)malloc(sizeof(options_t));
    assert(parser->options != NULL);
    parser->options->SPD    = false;
    parser->options->ITER   = false;
    parser->options->SPARSE = false;
    parser->options->ITOL   = DEFAULT_ITOL;

    /* Initializes the netlist struct that holds info about the elements */
    parser->netlist = (netlist_t *)malloc(sizeof(netlist_t));
    assert(parser->netlist != NULL);
    parser->netlist->num_nodes   = 0;
    parser->netlist->num_g2_elem = 0;
    parser->netlist->dc_counter  = 0;

    /* Initializes the dc_analysis array with a DC_ANALYSIS_NUM value that holds .DC options */
    parser->dc_analysis = (dc_analysis_t *)malloc(DC_ANALYSIS_NUM * sizeof(dc_analysis_t));
    assert(parser->dc_analysis != NULL);
    return parser;
}

/* Returns the number of tokens that exist in the current line */
int get_num_tokens(char *line) {
	int tokens = 0;
	int len = strlen(line);
	line[0] = (line[0] == '\t') ? ' ' : line[0];
	char found_char = (line[0] == ' ' ? 0 : 1);

	if (found_char == 1) {
		tokens++;
	}
	for (int i = 1; i < len; i++) {
		if (line[i] == ' ' || line[i] == '\t') {
			found_char = 0;
		}
		else if (line[i] > 32) {
			if (found_char == 0) {
				tokens++;
				found_char = 1;
			}
		}
	}
	return tokens;
}

/* Tokenizes the current line and returns the number of tokens and the tokens */
char **tokenizer(char *line) {
	const char delim[] = " \t";
	const char CRLF[]  = "\r\n";
	char **tokens, *token;
	short int i = 1;

	/* In case the line starts with a comment or it's an empty line ignore */
	if (line[0] == '*' || line[0] == '\n' || strcmp(line, CRLF) == 0) {
		return NULL;
	}

	int num_tokens = get_num_tokens(line);
	if (num_tokens == 0) {
		return NULL;
	}

	tokens = (char **)malloc((num_tokens + 1) * sizeof(char *));
	tokens[0] = (char *)malloc(sizeof(char));
	assert(tokens != NULL);
	assert(tokens[0] != NULL);
	sprintf(tokens[0], "%d", num_tokens);

	token = strtok(line, delim);

	while (token != NULL && i <= num_tokens) {
		tokens[i] = (char *)malloc(strlen(token) * sizeof(char));
		assert(tokens[i] != NULL);
		strcpy(tokens[i], token);
		token = strtok(NULL, delim);
		i++;
	}
	/* Trim '\n' if necessary */
	if (tokens[i-1][strlen(tokens[i-1])-1] == '\n') {
		tokens[i-1][strlen(tokens[i-1])-1] = '\0';
	}
	return tokens;
}

/* Parses the given netlist with <file_name> and returns the parser struct that holds all the required information */
void parse_netlist(parser_t *parser, char *file_name, index_t *index, hash_table_t *hash_table) {
    FILE *file_input;
    ssize_t read;
    size_t len = 0;
    int num_tokens = 0, dc_counter = 0;
    char **tokens = NULL, *line = NULL;

    printf("\nInput file is: %s\n", file_name);
    file_input = fopen(file_name, "rb");
    if (file_input == NULL) {
        fprintf(stderr, "Error: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
    /* Read line by line the netlist */
    while ((read = getline(&line, &len, file_input)) != -1) {
        tokens = tokenizer(line);
        if (tokens == NULL) {
            continue;
        }
        sscanf(tokens[0], "%d", &num_tokens);
        if (tokens[1][0] == '.') {
            if (strcmp(".OPTIONS", &tokens[1][0]) == 0) {
                for (int i = 2; i <= num_tokens; i++) {
                    if (strcmp("SPD", &tokens[i][0]) == 0) {
                        parser->options->SPD = true;
                    }
                    if (strcmp("ITER", &tokens[i][0]) == 0) {
                        parser->options->ITER = true;
                    }
                    if (strcmp("SPARSE", &tokens[i][0]) == 0) {
                        parser->options->SPARSE = true;
                    }
                    if (strncmp("ITOL", &tokens[i][0], 4) == 0) {
                        sscanf((&tokens[i][0]) + 5, "%lf", &parser->options->ITOL);
                    }
                }
            }
            else if (strcmp(".DC", &tokens[1][0]) == 0) {
                parser->dc_analysis[dc_counter].volt_source = (char *)malloc(strlen(&tokens[2][0]) * sizeof(char));
                assert(parser->dc_analysis[dc_counter].volt_source != NULL);
                sscanf(tokens[2], "%s",   parser->dc_analysis[dc_counter].volt_source);
                sscanf(tokens[3], "%lf", &parser->dc_analysis[dc_counter].start);
                sscanf(tokens[4], "%lf", &parser->dc_analysis[dc_counter].end);
                sscanf(tokens[5], "%lf", &parser->dc_analysis[dc_counter].increment);
            }
            else if (strcmp(".PLOT", &tokens[1][0]) == 0 || strcmp(".PRINT", &tokens[1][0]) == 0) {
                /* Allocate a 2d array that each row contains the node name */
                parser->dc_analysis[dc_counter].nodes = (char **)malloc(num_tokens * sizeof(char *));
                assert(parser->dc_analysis[dc_counter].nodes != NULL);
                parser->dc_analysis[dc_counter].num_nodes = 0;
                for (int i = 2; i <= num_tokens; i++) {
                    /* Allocate memory for the node name ommiting the parentheses and the V */
                    parser->dc_analysis[dc_counter].nodes[i-2] = (char *)malloc((strlen(tokens[i]) - 3) * sizeof(char));
                    assert(parser->dc_analysis[dc_counter].nodes[i-2] != NULL);
                    /* Strip V and the parentheses around node name */
                    strncpy(parser->dc_analysis[dc_counter].nodes[i-2], tokens[i] + 2, (strlen(tokens[i]) - 3));
                    parser->dc_analysis[dc_counter].num_nodes++;
                }
                dc_counter++;
            }
        }
        else {
            if (add_to_list(index, tokens, hash_table) == FAILURE) {
                exit(EXIT_FAILURE);
            } 
            if (tokens[1][0] == 'V' || tokens[1][0] == 'v' || tokens[1][0] == 'L' || tokens[1][0] == 'l') {
                parser->netlist->num_g2_elem++;
            }
        }
        /* Free all the memory we allocated */
        for (int i = 0; i < num_tokens; i++) {
            free(tokens[i]);
        }
        free(tokens);
    }
    if (line) {
        free(line);
    }
    fclose(file_input);

    /* DC counter holds the number of .DC we found in the netlist */
    parser->netlist->nz = get_nz();
    parser->netlist->dc_counter = dc_counter;
    parser->netlist->num_nodes  = hash_table->seq - 1;

#ifdef DEBUGL
    printf("Printing the lists\n");
    print_lists(index, hash_table);
#endif

    printf("Finished parsing %d circuit elements.\n", index->size1 + index->size2);
    print_options(parser->options);
    print_netlist_info(parser->netlist);
    print_dc_analysis_options(parser->dc_analysis, dc_counter);
}

/* Print all the specified options from the netlist */
void print_options(options_t *options) {
	printf("\nNetlist Specified Options:\n");
	printf("SPD:\t%s\n",    options->SPD    ? "true" : "false");
	printf("ITER:\t%s\n",   options->ITER   ? "true" : "false");
	printf("SPARSE:\t%s\n", options->SPARSE ? "true" : "false");
	printf("ITOL:\t%lf\n",  options->ITOL);
}

/* Print the number of the different netlist elements info */
void print_netlist_info(netlist_t *netlist) {
    printf("\nNetlist Elements:\n");
    printf("Number of nodes without ground:\t%d\n", netlist->num_nodes);
    printf("Number of group 2 elements:\t%d\n",     netlist->num_g2_elem);
    printf("Number of dc analysis targets:\t%d\n",  netlist->dc_counter);
}

void print_dc_analysis_options(dc_analysis_t *dc_analysis, int dc_counter) {
    printf("\nDC Analysis Summary:\n");
    for (int i = 0; i < dc_counter; i++) {
        printf("Volt_source: %s\n", dc_analysis[i].volt_source);
        printf("Start: %lf\n", dc_analysis[i].start);
        printf("End: %lf\n", dc_analysis[i].end);
        printf("Increment: %lf\n", dc_analysis[i].increment);
        printf("Node: ");
        for (int j = 0; j < dc_analysis[i].num_nodes; j++) {
            printf("%s ", dc_analysis[i].nodes[j]);
        }
        printf("\n");
    }
}

/* Free all the memory we allocated for the parser */
void free_parser(parser_t **parser) {
    /* Free options struct */
    free((*parser)->options);
    /* Free everything malloc'd in dc_analysis struct array */
    for (int i = 0; i < (*parser)->netlist->dc_counter; i++) {
        free((*parser)->dc_analysis[i].volt_source);
        for (int j = 0; j < (*parser)->dc_analysis[i].num_nodes; j++) {
            free((*parser)->dc_analysis[i].nodes[j]);
        }
        free((*parser)->dc_analysis[i].nodes);
    }
    /* Free the dc_analysis struct */
    free((*parser)->dc_analysis);
    /* Free netlist struct */
    free((*parser)->netlist);
    /* Free the parser */
    free((*parser));
    /* Set to NULL to limit further accesses */
    *parser = NULL;
}
