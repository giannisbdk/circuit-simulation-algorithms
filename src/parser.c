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
    parser->options->TR     = false;
    parser->options->BE     = false;
    parser->options->TRAN   = false;
    parser->options->AC     = false;
    parser->options->ITOL   = DEFAULT_ITOL;

    /* Initializes the netlist struct that holds info about the elements */
    parser->netlist = (netlist_t *)malloc(sizeof(netlist_t));
    assert(parser->netlist != NULL);
    parser->netlist->num_nodes   = 0;
    parser->netlist->num_g2_elem = 0;
    parser->netlist->dc_counter  = 0;

    /* Initializes the dc_analysis array with a DC_ANALYSIS_NUM value that holds .DC options */
    parser->dc_analysis = (dc_analysis_t *)malloc(ANALYSIS_NUM * sizeof(dc_analysis_t));
    assert(parser->dc_analysis != NULL);
    parser->dc_analysis->start       = 0.0;
    parser->dc_analysis->end         = 0.0;
    parser->dc_analysis->increment   = 0.0;
    parser->dc_analysis->num_nodes   = 0;
    parser->dc_analysis->volt_source = NULL;
    parser->dc_analysis->nodes = NULL;

    /* Initializes the tr_analsysis array */
    parser->tr_analysis = (tr_analysis_t *)malloc(ANALYSIS_NUM * sizeof(tr_analysis_t));
    assert(parser->tr_analysis != NULL);
    parser->tr_analysis->time_step = 0.0;
    parser->tr_analysis->fin_time  = 0.0;
    parser->tr_analysis->num_nodes = 0;
    parser->tr_analysis->nodes = NULL;

    /* Initializes the ac_analsysis array */
    parser->ac_analysis = (ac_analysis_t *)malloc(ANALYSIS_NUM * sizeof(ac_analysis_t));
    assert(parser->ac_analysis != NULL);
    parser->ac_analysis->start_freq = 0.0;
    parser->ac_analysis->end_freq   = 0.0;
    parser->ac_analysis->num_nodes  = 0;
    parser->ac_analysis->nodes = NULL;
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
    /* Find the length required to store the num_tokens in a string */
    int length = snprintf(NULL, 0, "%d", num_tokens);
    /* Allocate +1 for NULL termination */
	tokens[0] = (char *)malloc((length + 1) * sizeof(char));
	assert(tokens != NULL);
	assert(tokens[0] != NULL);
    /* Copy the integer to string */
    snprintf(tokens[0], length + 1, "%d", num_tokens);

	token = strtok(line, delim);

	while (token != NULL && i <= num_tokens) {
		tokens[i] = (char *)malloc((strlen(token) + 1) * sizeof(char));
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
    int num_tokens = 0, dc_counter = 0, tr_counter = 0, ac_counter = 0;
    bool tran = false, ac = false;
    char **tokens = NULL, *line = NULL;

    printf("\nInput file: %s\n", file_name);
    file_input = fopen(file_name, "rb");
    if (file_input == NULL) {
        fprintf(stderr, "Error: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    /* Read line by line the netlist */
    while ((read = getline(&line, &len, file_input)) != -1) {
        tokens = tokenizer(line);
        if (tokens == NULL) continue;
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
                    if (strcmp("METHOD=BE", &tokens[i][0]) == 0) {
                        parser->options->BE = true;
                    }
                    else if(strcmp("METHOD=TR", &tokens[i][0]) == 0) {
                        parser->options->TR = true;
                    }
                }
            }
            else if (strcmp(".DC", &tokens[1][0]) == 0) {
                parser->dc_analysis[dc_counter].volt_source = (char *)malloc((strlen(&tokens[2][0]) + 1) * sizeof(char));
                assert(parser->dc_analysis[dc_counter].volt_source != NULL);
                sscanf(tokens[2], "%s",   parser->dc_analysis[dc_counter].volt_source);
                sscanf(tokens[3], "%lf", &parser->dc_analysis[dc_counter].start);
                sscanf(tokens[4], "%lf", &parser->dc_analysis[dc_counter].end);
                sscanf(tokens[5], "%lf", &parser->dc_analysis[dc_counter].increment);
            }
            else if (strcmp(".TRAN", &tokens[1][0]) == 0) {
                sscanf(tokens[2], "%lf", &parser->tr_analysis[tr_counter].time_step);
                sscanf(tokens[3], "%lf", &parser->tr_analysis[tr_counter].fin_time);
                tran = true;
            }
            else if (strcmp(".AC", &tokens[1][0]) == 0) {
                if (strcmp(tokens[2], "LIN") == 0) {
                    parser->ac_analysis[ac_counter].sweep = LIN;
                } 
                else {
                    parser->ac_analysis[ac_counter].sweep = LOG;
                }
                sscanf(tokens[3], "%d",  &parser->ac_analysis[ac_counter].points);
                sscanf(tokens[4], "%lf", &parser->ac_analysis[ac_counter].start_freq);
                sscanf(tokens[5], "%lf", &parser->ac_analysis[ac_counter].end_freq);
                ac = true;
            }
            else if (strcmp(".PLOT", &tokens[1][0]) == 0 || strcmp(".PRINT", &tokens[1][0]) == 0) {
                if (tran) {
                    /* Allocate a 2d array that each row contains the node name */
                    parser->tr_analysis[tr_counter].nodes = (char **)malloc(num_tokens * sizeof(char *));
                    assert(parser->tr_analysis[tr_counter].nodes != NULL);
                    parser->tr_analysis[tr_counter].num_nodes = 0;
                    for (int i = 2; i <= num_tokens; i++) {
                        /* Allocate memory for the node name ommiting the parentheses and the V */
                        /* Add +1 to include null termination */
                        int length = strlen(tokens[i]) - 3 + 1;
                        parser->tr_analysis[tr_counter].nodes[i-2] = (char *)malloc(length * sizeof(char));
                        assert(parser->tr_analysis[tr_counter].nodes[i-2] != NULL);
                        /* Strip V and the parentheses around node name */
                        snprintf(parser->tr_analysis[tr_counter].nodes[i-2], length * sizeof(char), "%s", tokens[i] + 2);
                        parser->tr_analysis[tr_counter].num_nodes++;
                    }
                    tr_counter++;
                    tran = false;
                }
                else if (ac) {
                    /* Allocate a 2d array that each row contains the node name */
                    parser->ac_analysis[ac_counter].nodes = (char **)malloc(num_tokens * sizeof(char *));
                    assert(parser->ac_analysis[ac_counter].nodes != NULL);
                    parser->ac_analysis[ac_counter].num_nodes = 0;
                    for (int i = 2; i <= num_tokens; i++) {
                        /* Allocate memory for the node name ommiting the parentheses and the V */
                        /* Add +1 to include null termination */
                        int length = strlen(tokens[i]) - 3 + 1;
                        parser->ac_analysis[ac_counter].nodes[i-2] = (char *)malloc(length * sizeof(char));
                        assert(parser->ac_analysis[ac_counter].nodes[i-2] != NULL);
                        /* Strip V and the parentheses around node name */
                        snprintf(parser->ac_analysis[ac_counter].nodes[i-2], length * sizeof(char), "%s", tokens[i] + 2);
                        parser->ac_analysis[ac_counter].num_nodes++;
                    }
                    ac_counter++;
                    ac = false;
                }
                else {
                    /* Allocate a 2d array that each row contains the node name */
                    parser->dc_analysis[dc_counter].nodes = (char **)malloc(num_tokens * sizeof(char *));
                    assert(parser->dc_analysis[dc_counter].nodes != NULL);
                    parser->dc_analysis[dc_counter].num_nodes = 0;
                    for (int i = 2; i <= num_tokens; i++) {
                        /* Allocate memory for the node name ommiting the parentheses and the V */
                        /* Add +1 to include null termination */
                        int length = strlen(tokens[i]) - 3 + 1;
                        parser->dc_analysis[dc_counter].nodes[i-2] = (char *)malloc(length * sizeof(char));
                        assert(parser->dc_analysis[dc_counter].nodes[i-2] != NULL);
                        /* Strip V and the parentheses around node name */
                        snprintf(parser->dc_analysis[dc_counter].nodes[i-2], length * sizeof(char), "%s", tokens[i] + 2);
                        parser->dc_analysis[dc_counter].num_nodes++;
                    }
                    dc_counter++;
                }
            }
        }
        else {
            if (add_to_list(index, tokens, hash_table) == FAILURE) {
                exit(EXIT_FAILURE);
            }
            /* Keep track of how many group 2 elements there are in the nelist */
            if (tokens[1][0] == 'V' || tokens[1][0] == 'v' || tokens[1][0] == 'L' || tokens[1][0] == 'l') {
                parser->netlist->num_g2_elem++;
            }
        }
        /* Free all the memory we allocated */
        for (int i = 0; i <= num_tokens; i++) {
            free(tokens[i]);
        }
        free(tokens);
    }
    if (line) {
        free(line);
    }
    fclose(file_input);

    /* dc,tr,ac counters holds the number of .DC/TRAN/AC we found in the netlist */
    parser->netlist->nz = get_nz();
    parser->netlist->dc_counter = dc_counter;
    parser->netlist->tr_counter = tr_counter;
    parser->netlist->ac_counter = ac_counter;
    parser->netlist->num_nodes  = hash_table->seq - 1;

    /* Set the appropriate flags according to the counters */
    if (parser->netlist->tr_counter) {
        parser->options->TRAN = true;
    }
    if (parser->netlist->ac_counter) {
        parser->options->AC = true;
    }

    /* In case there exists a .TRAN analysis and METHOD field is absent we set it to Trapezoidal as default */
    if (parser->options->TRAN && !parser->options->BE && !parser->options->TR) {
        parser->options->TR = true;
    }

#ifdef DEBUGL
    printf("Printing the lists\n");
    print_lists(index, hash_table);
#endif

    printf("\nFinished parsing %d circuit elements.\n", index->size1 + index->size2);
    print_options(parser->options);
    print_netlist_info(parser->netlist);
    print_dc_analysis_options(parser->dc_analysis, dc_counter);
    print_tr_analysis_options(parser->tr_analysis, tr_counter);
    print_ac_analysis_options(parser->ac_analysis, ac_counter);
}

/* Print all the specified options from the netlist */
void print_options(options_t *options) {
	printf("\n--- Netlist Specified Options ---\n");
	printf("SPD:\t%s\n",    options->SPD    ? "true" : "false");
	printf("ITER:\t%s\n",   options->ITER   ? "true" : "false");
	printf("SPARSE:\t%s\n", options->SPARSE ? "true" : "false");
    printf("TRAN:\t%s\n",   options->TRAN   ? "true" : "false");
    printf("TR:\t%s\n",     options->TR     ? "true" : "false");
    printf("BE:\t%s\n",     options->BE     ? "true" : "false");
    printf("AC:\t%s\n",     options->AC     ? "true" : "false");
    printf("ITOL:\t%g\n",   options->ITOL);
}

/* Print the number of the different netlist elements info */
void print_netlist_info(netlist_t *netlist) {
    printf("\n--- Netlist Elements and Analyses Summary ---\n");
    printf("Number of nodes without ground: %d\n", netlist->num_nodes);
    printf("Number of group 2 elements:     %d\n", netlist->num_g2_elem);
    printf("Number of DC analyses:          %d\n", netlist->dc_counter);
    printf("Number of Transient analyses:   %d\n", netlist->tr_counter);
    printf("Number of AC analyses:          %d\n", netlist->ac_counter);
}

/* Prints the dc analysis options */
void print_dc_analysis_options(dc_analysis_t *dc_analysis, int dc_counter) {
    if (dc_counter <= 0) {
        return;
    }
    else if (dc_counter == 1) {
        printf("\n--- DC Analysis Summary ---\n");
    }
    else {
        printf("\n--- DC Analyses Summary ---\n");
    }
    for (int i = 0; i < dc_counter; i++) {
        printf("Voltage source: %s\n",  dc_analysis[i].volt_source);
        printf("Start:          %lf\n", dc_analysis[i].start);
        printf("End:            %lf\n", dc_analysis[i].end);
        printf("Increment:      %lf\n", dc_analysis[i].increment);
        printf("Nodes:          ");
        for (int j = 0; j < dc_analysis[i].num_nodes; j++) {
            printf("%s ", dc_analysis[i].nodes[j]);
        }
        printf("\n");
        if (i < (dc_counter - 1)) printf("\n");
    }
}

/* Prints the transient analysis options */
void print_tr_analysis_options(tr_analysis_t *tr_analysis, int tr_counter) {
    if (tr_counter <= 0) {
        return;
    }
    else if (tr_counter == 1) {
        printf("\n--- Transient Analysis Summary ---\n");
    }
    else {
        printf("\n--- Transient Analyses Summary ---\n");
    }
    for (int i = 0; i < tr_counter; i++) {
        printf("Time Step:  %lf\n", tr_analysis[i].time_step);
        printf("Final Time: %lf\n", tr_analysis[i].fin_time);
        printf("Nodes:      ");
        for (int j = 0; j < tr_analysis[i].num_nodes; j++) {
            printf("%s ", tr_analysis[i].nodes[j]);
        }
        printf("\n");
        if (i < (tr_counter - 1)) printf("\n");
    }
}

/* Prints the dc analysis options */
void print_ac_analysis_options(ac_analysis_t *ac_analysis, int ac_counter) {
    if (ac_counter <= 0) {
        return;
    }
    else if (ac_counter == 1) {
        printf("\n--- AC Analysis Summary ---\n");
    }
    else {
        printf("\n--- AC Analyses Summary ---\n");
    }
    for (int i = 0; i < ac_counter; i++) {
        printf("Sweep:      %s\n",  ac_analysis[i].sweep ? "LOG" : "LIN");
        printf("Points:     %d\n",  ac_analysis[i].points);
        printf("Start freq: %g\n", ac_analysis[i].start_freq);
        printf("End freq:   %g\n", ac_analysis[i].end_freq);
        printf("Nodes:      ");
        for (int j = 0; j < ac_analysis[i].num_nodes; j++) {
            printf("%s ", ac_analysis[i].nodes[j]);
        }
        printf("\n");
        if (i < (ac_counter - 1)) printf("\n");
    }
}

/* Free all the memory we allocated for the parser */
void free_parser(parser_t **parser) {
    /* Free options struct */
    free((*parser)->options);

    /* Free everything allocated for DC analysis */
    if ((*parser)->netlist->dc_counter) {
        /* Free everything malloc'd in dc_analysis struct array */
        for (int i = 0; i < (*parser)->netlist->dc_counter; i++) {
            free((*parser)->dc_analysis[i].volt_source);
            for (int j = 0; j < (*parser)->dc_analysis[i].num_nodes; j++) {
                free((*parser)->dc_analysis[i].nodes[j]);
            }
            free((*parser)->dc_analysis[i].nodes);
        }
    }
    /* Free the dc_analysis struct */
    free((*parser)->dc_analysis);

    /* Free everything allocated for TRANSIENT analyiss */
    if ((*parser)->netlist->tr_counter) {
        /* Free everything malloc'd in tr_analysis struct array */
        for (int i = 0; i < (*parser)->netlist->tr_counter; i++) {
            for (int j = 0; j < (*parser)->tr_analysis[i].num_nodes; j++) {
                free((*parser)->tr_analysis[i].nodes[j]);
            }
            free((*parser)->tr_analysis[i].nodes);
        }
    }
    /* Free the tr_analysis struct */
    free((*parser)->tr_analysis);  

    /* Free everything allocated for AC analysis */
    if ((*parser)->netlist->ac_counter) {
        /* Free everything malloc'd in ac_analysis struct array */
        for (int i = 0; i < (*parser)->netlist->ac_counter; i++) {
            for (int j = 0; j < (*parser)->ac_analysis[i].num_nodes; j++) {
                free((*parser)->ac_analysis[i].nodes[j]);
            }
            free((*parser)->ac_analysis[i].nodes);
        }
    }
    /* Free the tr_analysis struct */
    free((*parser)->ac_analysis); 

    /* Free netlist struct */
    free((*parser)->netlist);

    /* Free the parser */
    free((*parser));

    /* Set to NULL to limit further accesses */
    *parser = NULL;
}
