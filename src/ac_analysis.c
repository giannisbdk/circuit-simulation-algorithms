#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "ac_analysis.h"

/* Do AC analysis */
void ac_analysis(index_t *index, hash_table_t *hash_table, mna_system_t *mna, parser_t *parser,
				 double *dc_op, gsl_vector_complex *sol_x) {
	/* Set the prefix name for the files */
	char prefix[] = "ac_analysis_";
	char file_name[MAX_FILE_NAME];

	int ac_counter = parser->netlist->ac_counter;
	mna->ac_analysis_init = true;

	/* Run all the AC analysis according to the ac_counter */
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
			create_dense_ac_mna(mna, index, hash_table, parser->options, parser->netlist->num_nodes, omega);
			/* Solve the system */
			solve_mna_system(mna, &dc_op, sol_x, parser->options);
			/* Print current solution to file */
			for (int j = 0; j < parser->ac_analysis[i].num_nodes; j++) {
				int offset = ht_get_id(hash_table, parser->ac_analysis[i].nodes[j]) - 1;
				/* Convert from complex to polar with magnitude and phase */
				ac_t curr_ac = rect_to_polar(gsl_vector_complex_get(sol_x, offset));
				fprintf(files[j], "%-15lf%-15lf%-15lf\n", sweep_points_freq[step], curr_ac.magnitude, curr_ac.phase);
			}
		}
		/* Close the file descriptors for the current transient analysis */
		for (int j = 0; j < parser->dc_analysis[i].num_nodes; j++) {
			fclose(files[j]);
		}
		/* Free sweep points matrix because in next AC analysis number of points might differ, we've to allocate again */
		free(sweep_points_freq);
	}
	mna->ac_analysis_init = false;
	if (parser->netlist->ac_counter) {
		printf("AC Analysis..............OK\n");
	}
}

/* According to the specified sweep type it fills an array containing all the sweep points */
void get_sweep_points(double *array, ac_analysis_t ac_analysis) {
	switch (ac_analysis.sweep) {
		case LIN:
			lin_step(array, ac_analysis.start_freq, ac_analysis.end_freq, ac_analysis.points);
			return;
		case LOG:
			log_step(array, ac_analysis.start_freq, ac_analysis.end_freq, ac_analysis.points);
			return;
		default:
			fprintf(stderr, "Error: Wrong sweep type.\n");
			exit(EXIT_FAILURE);
	}
}

/* Fills the array with the linear steps, for a linear sweep between start-end freq with n points */
void lin_step(double *array, double start, double end, int points) {
	/* Check for valid points */
	if (points < 2) {
		fprintf(stderr, "Error: Creating linear step, points must be >= 2.\n");
	}
	double step = (end - start) / (points - 1);
	for (int i = 0; i < points; i++) {
		array[i] = start + i * step;
	}
}

/* Fills the array with the log steps, for a log sweep between start-end freq with n points */
void log_step(double *array, double start, double end, int points) {
	/* Check for valid points */
	if (points < 2) {
		fprintf(stderr, "Error: Creating log step, points must be >=2.\n");
	}
	double step = (end - start) / (points - 1);
	/* Logarithmic base of 10 */
	double base = 10.0;
	for (int i = 0; i < points; i++) {
		array[i] = pow(base, start + i * step);
	}
}