#ifndef PARSER_H
#define PARSER_H

#include <stdbool.h>

#include "list.h"
#include "hash_table.h"

#define DEFAULT_ITOL    0.001
#define ANALYSIS_NUM 	5

extern int errno;

/* Struct to hold the different options for the analyses */
typedef struct options {
	bool SPD;
	bool SPARSE;
	bool ITER;
	bool TRAN;
	bool TR;
	bool BE;
	bool AC;
	double ITOL;
} options_t;

/*
 * Struct to hold the information of netlists, such as
 * number of nodes, number of group2 elements and
 * number of dc analysis targets
 */
typedef struct netlist {
	int num_nodes;
	int num_g2_elem;
	int dc_counter;
	int tr_counter;
	int ac_counter;
	int nz;
} netlist_t;

/* Struct to hold the different DC analysis options */
typedef struct dc_analysis {
	char  *volt_source;
	double start;
	double end;
	double increment;
	/* Plot / Print */
	char **nodes;
	int num_nodes;
} dc_analysis_t;

/* Struct to hold the different TRAN analysis options */
typedef struct tr_analysis {
	double time_step;
	double fin_time;
	/* Plot / Print */
	char **nodes;
	int num_nodes;
} tr_analysis_t;

/* Sweep enumeration for AC analysis */
typedef enum ac_sweep {
	LIN,
	LOG
} ac_sweep_t;

/* Struct to hold the different AC analysis options */
typedef struct ac_analysis {
	double start_freq;
	double end_freq;
	int points;
	ac_sweep_t sweep;
	/* Plot / Print */
	char **nodes;
	int num_nodes;
} ac_analysis_t;

/* 
 * Struct to hold all the previous values, so that we are able
 * to return all the information we read in the netlist
 */
typedef struct parser {
	options_t 	  *options;
	dc_analysis_t *dc_analysis;
	tr_analysis_t *tr_analysis;
	ac_analysis_t *ac_analysis;
	netlist_t 	  *netlist;
} parser_t;

parser_t *init_parser();
int get_num_tokens(char *line);
char **tokenizer(char *line);
void parse_netlist(parser_t *parser, char *file_name, index_t *index, hash_table_t *hash_table);
void print_options(options_t *options);
void print_netlist_info(netlist_t *netlist);
void print_dc_sweep_analysis_options(dc_analysis_t *dc_analysis, int dc_counter);
void print_tr_analysis_options(tr_analysis_t *tr_analysis, int tr_counter);
void print_ac_analysis_options(ac_analysis_t *ac_analysis, int ac_counter);
void free_parser(parser_t **parser);

#endif