#ifndef PARSER_H
#define PARSER_H

#include <stdbool.h>

#include "list.h"
#include "hash_table.h"

#define DEFAULT_ITOL    0.001
#define DC_ANALYSIS_NUM 25

extern int errno;

/* Struct to hold the different options for analysis */
typedef struct options {
	bool SPD;
	bool SPARSE;
	bool ITER;
	double itol;
} options_t;

/* Struct to hold the elements of netlists, such as  *
 * number of nodes, number of g2 elements and number *
 * of dc analysis                                    */
typedef struct netlist_elem {
	int num_nodes;
	int num_g2_elem;
	int dc_counter;
} netlist_elem_t;

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

/* Struct to hold all the previous values, so that we are able *
 * to return all the information we read in the netlist        */
typedef struct parser {
	options_t *options;
	dc_analysis_t *dc_analysis;
	netlist_elem_t *netlist_elem;
} parser_t;

parser_t *init_parser();
void parse_netlist(parser_t *parser, char *file_name, index_t *index, hash_table_t *hash_table);
void print_options(options_t *options);
void print_netlist_elem(netlist_elem_t *netlist_elem);
int get_num_tokens(char *line);
char **tokenizer(char *line);
void free_parser(parser_t **parser);

#endif