#include <stdio.h>
#include <string.h>

#include "ac_analysis.h"

/* Do AC analysis */
void ac_analysis(index_t *index, hash_table_t *hash_table, mna_system_t *mna, parser_t *parser, double *sol_x) {
	/* Set the prefix name for the files */
	char prefix[] = "ac_analysis_";
	char file_name[MAX_FILE_NAME];

	int ac_counter = parser->netlist->ac_counter;

	for (int i = 0; i < ac_counter; i++) {
		FILE *files[parser->ac_analysis[i].num_nodes];
		/* Open different files for each node in plot/print array */
		for (int j = 0; j < parser->ac_analysis[i].num_nodes; j++) {
			/* Temp buffer */
			char name[10];
			/* Construct the file name */
			strcpy(file_name, prefix);
			strcat(file_name, parser->ac_analysis[i].nodes[j]);
			strcat(file_name, "_");
			sprintf(name, "%g", parser->ac_analysis[i].start_freq);
			strcat(file_name, name);
			strcat(file_name, "_");
			sprintf(name, "%g", parser->ac_analysis[i].end_freq);
			strcat(file_name, name);
			strcat(file_name, ".txt");
			/* Open the output file */
			files[j] = fopen(file_name, "w");
			if (files[j] == NULL) {
				fprintf(stderr, "Error opening file: %s\n", strerror(errno));
				exit(EXIT_FAILURE);
			}
			fprintf(files[j], "%-15s%-15s%-15s\n", "Frequency", "Magnitude", "Phase");
		}
		double omega = parser->ac_analysis[ac_counter].start_freq * 2 * M_PI;
		create_dense_ac_mna(mna, index, hash_table, parser->options, parser->netlist->num_nodes, omega);
		/* Find how many steps are required */
		double freq_step = get_freq_step(parser->ac_analysis[ac_counter]);
		int n_steps = parser->ac_analysis[ac_counter].points;
		for (int step = 0; step < n_steps; step++) {
			/* Solve the system */
			solve_mna_system(mna, &sol_x, parser->options);
			/* Print current solution to file */
			for (int j = 0; j < parser->ac_analysis[i].num_nodes; j++) {
				int offset = ht_get_id(hash_table, parser->ac_analysis[i].nodes[j]) - 1;
				ac_t curr_ac = rect_to_polar(sol_x[offset]);
				fprintf(files[j], "%-15lf%-15lf%-15lf\n", omega / (2 * M_PI), curr_ac.magnitude, curr_ac.phase);
			}
			/* Prepare matrices, with new frequency, for the next step */
			omega += freq_step;
			create_dense_ac_mna(mna, index, hash_table, parser->options, parser->netlist->num_nodes, omega);
		}
		/* Close the file descriptors for the current transient analysis */
		for (int j = 0; j < parser->dc_analysis[i].num_nodes; j++) {
		    fclose(files[j]);
		}
	}
}

double get_freq_step(ac_analysis_t ac_analysis) {
	switch (ac_analysis.sweep) {
		case LIN:
			return (ac_analysis.end_freq - ac_analysis.start_freq) / (ac_analysis.points - 1);
		case LOG:
			return 0.0;
		default:
			return 0.0;
	}
}