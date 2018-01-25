#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "ac_analysis.h"

/* Do AC analysis */
void ac_analysis(index_t *index, hash_table_t *hash_table, mna_system_t *mna, parser_t *parser,
				 double *dc_op, gsl_vector_complex *sol_x) {
	/* Set the flag that we're currently on an AC analysis */
	int ac_counter = parser->netlist->ac_counter;
	mna->ac_analysis_init = true;

	/* Run all the AC analysis according to the ac_counter */
	for (int i = 0; i < ac_counter; i++) {
		/* Create an array with files for every node and create/open them */
		FILE *files[parser->ac_analysis[i].num_nodes];
		create_ac_out_files(files, parser->ac_analysis[i]);;

		/* Find how many steps are required */
		int n_steps = parser->ac_analysis[i].points;

		/* Allocate a temporary array to hold the linspace and logspace sweep points */
		double *sweep_points_freq = (double *)malloc(n_steps * sizeof(double));
		assert(sweep_points_freq != NULL);

		/* Fill sweep_points_freq buffer. This buffer contains all the sweep points (freq) for the AC analysis */
		get_sweep_points(sweep_points_freq, parser->ac_analysis[i]);

		/* For every sweep point/step solve the according MNA AC system */
		for (int step = 0; step < n_steps; step++) {
			/* Find the current ω = 2πf */
			double omega = 2 * M_PI * sweep_points_freq[step];
			/* Create the AC MNA matrix for the AC analysis */
			create_ac_mna_system(mna, index, hash_table, parser->options, parser->netlist->num_nodes, omega);
			/* Solve the system */
			solve_mna_system(mna, &dc_op, sol_x, parser->options);
			/* Print the output to files */
			write_ac_out_files(files, parser->ac_analysis[i], hash_table, sol_x, sweep_points_freq[step]);
		}
		/* Close the file descriptors for the current transient analysis */
		for (int j = 0; j < parser->dc_analysis[i].num_nodes; j++) {
			fclose(files[j]);
		}
		/* Free sweep points matrix because in next AC analysis number of points might differ, we've to allocate again */
		free(sweep_points_freq);
	}
	/* Indicate that we stopped the AC analysis */
	mna->ac_analysis_init = false;
	if (parser->netlist->ac_counter) {
		printf("AC Analysis..............OK\n");
	}
}

/* According to the specified sweep type it fills an array containing all the sweep points */
void get_sweep_points(double *array, ac_analysis_t ac_analysis) {
	switch (ac_analysis.sweep) {
		case LIN:
			lin_sweep(array, ac_analysis.start_freq, ac_analysis.end_freq, ac_analysis.points);
			return;
		case LOG:
			log_sweep(array, log10(ac_analysis.start_freq), log10(ac_analysis.end_freq), ac_analysis.points);
			return;
		default:
			fprintf(stderr, "Error: Wrong sweep type.\n");
			exit(EXIT_FAILURE);
	}
}

/* Fills the array with the linear steps, for a linear sweep between start-end freq with n points */
void lin_sweep(double *array, double start, double end, int points) {
	/* Check for valid points */
	if (points < 2) {
		fprintf(stderr, "Error: Creating linear step, points must be >= 2.\n");
	}
	double step = (end - start) / (points - 1);
	for (int i = 0; i < points; i++) {
		array[i] = start + i * step;
	}
}

/*
 * Fills the array with the log steps, for a log sweep between start-end freq with n points
 * Common usage is to supply the function with arguments in log scale as log10(start), log10(end).
 */
void log_sweep(double *array, double start, double end, int points) {
	/* Check for valid points */
	if (points < 2) {
		fprintf(stderr, "Error: Creating log step, points must be >= 2.\n");
	}
	/* Logarithmic base of 10 */
	double base = 10.0;
	double step = (end - start) / (points - 1);
	/* Set initial and last value to the start and end */
	for (int i = 0; i < points; i++) {
		array[i] = pow(base, start + i * step);
	}
}

/* Creates and opens output files for every node included in the current AC analysis */
void create_ac_out_files(FILE *files[], ac_analysis_t ac_analysis) {
	/* Set the prefix name for the files */
	char prefix[] = "ac_analysis_V(";
	char file_name[MAX_FILE_NAME];

	/* Open different files for each node in plot/print array */
	for (int j = 0; j < ac_analysis.num_nodes; j++) {
		/* Temp buffer */
		char name[10];
		/* Will contain a string either Magnitude (voltage) or Magnitude (dB), according to the sweep type */
		char magn_output[] = "Magnitude ";
		/* Construct the file name */
		strcpy(file_name, prefix);
		strcat(file_name, ac_analysis.nodes[j]);
		strcat(file_name, ")_");
		sprintf(name, "%g", ac_analysis.start_freq);
		strcat(file_name, name);
		strcat(file_name, "_");
		sprintf(name, "%g", ac_analysis.end_freq);
		strcat(file_name, name);
		/* Set the output string for the magnitude according to the sweep type */
		switch (ac_analysis.sweep) {
			case LIN:
				strcat(file_name, "_LIN");
				strcat(magn_output, "(voltage)");
				break;
			case LOG:
				strcat(file_name, "_LOG");
				strcat(magn_output, "(dB)");
				break;
			default:
				fprintf(stderr, "Wrong sweep type.\n");
				exit(EXIT_FAILURE);
		}
		strcat(file_name, ".txt");
		/* Open the output file */
		files[j] = fopen(file_name, "w");
		if (files[j] == NULL) {
			fprintf(stderr, "Error opening file: %s\n", strerror(errno));
			exit(EXIT_FAILURE);
		}
		fprintf(files[j], "%-30s%-30s%-30s\n", "Frequency (hertz)", magn_output, "Phase (degrees)");
	}
}

/* Writes the output of the current step of the AC analysis to the output files */
void write_ac_out_files(FILE *files[], ac_analysis_t ac_analysis, hash_table_t *hash_table, gsl_vector_complex *sol_x, double freq_step) {
	/* Print current solution to file */
	for (int j = 0; j < ac_analysis.num_nodes; j++) {
		int offset = ht_get_id(hash_table, ac_analysis.nodes[j]) - 1;
		/* Convert from complex to polar with magnitude and phase */
		ac_t curr_ac = rect_to_polar(gsl_vector_complex_get(sol_x, offset));
		/* In case it is a LOG sweep we need to output 20*log10(magnitude) */
		switch (ac_analysis.sweep) {
			case LIN:
				fprintf(files[j], "%-30lf%-30lf%-30lf\n", freq_step, curr_ac.magnitude, curr_ac.phase);
				break;
			case LOG:
				fprintf(files[j], "%-30lf%-30lf%-30lf\n", freq_step, 20 * log10(curr_ac.magnitude), curr_ac.phase);
				break;
			default:
				fprintf(stderr, "Wrong sweep type.\n");
				exit(EXIT_FAILURE);
		}
	}
}