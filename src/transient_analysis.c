#include <math.h>
#include <string.h>

#include "transient_analysis.h"

/* Do transient analysis */
void tr_analysis(hash_table_t *hash_table, mna_system_t *mna, parser_t *parser, double *init_sol, double *sol_x) {
	/* Clear the decomposition flag for the mna system */
	mna->is_decomp = false;
	mna->tr_analysis_init = true;

	/* b is the RHS of the trapezoidal/backward euler */
	double *prev_sol = init_vector(mna->dimension);
	double *curr_response = init_vector(mna->dimension);
	double *prev_response = NULL;

	/* In case it is trapezoidal method previous response e(tkn-1) is required */
	if (parser->options->TR) {
		prev_response = init_vector(mna->dimension);
	}

	for (int i = 0; i < parser->netlist->tr_counter; i++) {
		/* Create an array with files for every node and create/open them */
		FILE *files[parser->tr_analysis[i].num_nodes];
		create_tr_out_files(files, parser->tr_analysis[i]);

		/* Store the initial values to the prev_ vectors */
		memcpy(prev_sol, init_sol, mna->dimension * sizeof(double));
		if (parser->options->TR) {
			memcpy(prev_response, mna->resp->value, mna->dimension * sizeof(double));
		}

		/* Find how many steps are required */
		int n_steps = parser->tr_analysis->fin_time / parser->tr_analysis->time_step;

		for (int step = 0; step <= n_steps; step++) {
			if (parser->options->TR) {
				set_trapezoidal_rhs(mna, curr_response, prev_response, prev_sol, parser->tr_analysis->time_step, step, parser->options->SPARSE);
			}
			else {
				set_backward_euler_rhs(mna, curr_response, prev_sol, parser->tr_analysis->time_step, step, parser->options->SPARSE);
			}
			/* Solve the system */
			solve_mna_system(mna, &sol_x, NULL, parser->options);
			/* Print the output to files */
			write_tr_out_files(files, parser->tr_analysis[i], hash_table, sol_x, step);
			/* Copy current solution to prev to use for next iteration */
			memcpy(prev_sol, sol_x, mna->dimension * sizeof(double));
		}
		/* Close the file descriptors for the current transient analysis */
		for (int j = 0; j < parser->tr_analysis[i].num_nodes; j++) {
		    fclose(files[j]);
		}
	}

	/* Set flag to false to indicate that we're no longer inside TRANSIENT analysis */
	mna->tr_analysis_init = false;
	if (parser->netlist->tr_counter) {
		printf("Transient Analysis.......OK\n");
	}

	/* Free everything we allocated */
	free(prev_sol);
	free(curr_response);
	if (parser->options->TR) {
		free(prev_response);
	}
}

/*
 * Computes and stores the right hand side of the trapezoidal method at vector mna->b
 * h: is time_step and k: is the current iteration
 */ 
void set_trapezoidal_rhs(mna_system_t *mna, double *curr_response, double *prev_response, double *prev_sol, double h, int k, bool SPARSE) {
	/* curr_response is e(tk) and prev_response is e(tk-1) */
	/* Set the values of the e(tk) vector */
	set_response_vector(curr_response, mna->resp, h * k, mna->dimension);
	/* Compute: e(tk) + e(tk-1) - sGhC*x(tk-1) and save it to the provided b vector */
	double *response_add = init_vector(mna->dimension);
	double *sGhc_x = init_vector(mna->dimension);
	add_vector(response_add, curr_response, prev_response, mna->dimension);
	if (SPARSE) {
		// cs_mat_vec_mul(sGhc_x, mna->sp_matrix->sGhC, prev_sol);
		cs_gaxpy(mna->sp_matrix->sGhC, prev_sol, sGhc_x);
	}
	else {
		mat_vec_mul(sGhc_x, mna->matrix->sGhC, prev_sol, mna->dimension);
	}
	sub_vector(mna->b, response_add, sGhc_x, mna->dimension);

	/* Copy the curr_response to prev_response to use for the next call */
	memcpy(prev_response, curr_response, mna->dimension * sizeof(double));

	/* Free everything we allocated */
	free(response_add);
	free(sGhc_x);
}

/* 
 * Computes and stores the right hand side of the backward euler method at vector mna->b
 * h: is time_step and k: is the current iteration
 */
void set_backward_euler_rhs(mna_system_t *mna, double *curr_response, double *prev_sol, double h, int k, bool SPARSE) {
	/* Set the values of the curr_response e(tk) vector */
	set_response_vector(curr_response, mna->resp, h * k, mna->dimension);
	double *hC_x = init_vector(mna->dimension);
	/* Compute: e(tk) + (1/h)C*x(tk-1) and save it to the mna->b vector */
	if (SPARSE) {
		cs_gaxpy(mna->sp_matrix->hC, prev_sol, hC_x);
		// cs_mat_vec_mul(hC_x, mna->sp_matrix->hC, prev_sol);
	}
	else {
		mat_vec_mul(hC_x, mna->matrix->hC, prev_sol, mna->dimension);
	}
	add_vector(mna->b, curr_response, hC_x, mna->dimension);

	/* Free everything we allocated */
	free(hC_x);
}

/* Set the values of the e(tk) vector */
void set_response_vector(double *response, resp_t *resp, double time, int dimension) {
	list1_t *node;
	for (int i = 0; i < dimension; i++) {
		node = resp->nodes[i];
		if (node != NULL) {
			if (node->trans_spec != NULL) {
				switch (node->trans_spec->type) {
					case EXP:
						response[i] = eval_exp(node->trans_spec->exp, time);
						break;
					case SIN:
						response[i] = eval_sin(node->trans_spec->sin, time);
						break;
					case PULSE:
						response[i] = eval_pulse(node->trans_spec->pulse, time);
						break;
					case PWL:
						response[i] = eval_pwl(node->trans_spec->pwl, time);
						break;
				}
			}
			else {
				/* If transient spec is absent set the DC value */
				response[i] = resp->value[i];
			}
		}
	}
}

/* Evaluates the exponential transient at given time t */
double eval_exp(exp_t *expon, double t) {
	if (0 <= t && t < expon->td1) {
		return expon->i1;
	}
	else if (expon->td1 <= t && t < expon->td2) {
		return expon->i1 + (expon->i2 - expon->i1) * (1 - exp(-(t - expon->td1) / expon->tc1));
	}
	else {
		return expon->i1 + (expon->i2 - expon->i1) *
				(exp(-(t - expon->td2) / expon->tc2) - exp(-(t - expon->td1) / expon->tc1));
	}
}

/* Evaluates the sinusodial transient at given time t */
double eval_sin(sin_t *sinus, double t) {
	if (0 <= t && t < sinus->td) {
		return sinus->i1 + sinus->ia * sin(2 * M_PI * sinus->ph / 360);
	}
	else {
		return sinus->i1 + sinus->ia *
				sin(2 * M_PI * sinus->fr * (t - sinus->td) + 2 * M_PI * sinus->ph / 360) *
				exp(-(t - sinus->td) * sinus->df);
	}
}

/* Evaluates the pulse transient at given time t */
double eval_pulse(pulse_t *pulse, double t) {
	int curr_per = t / pulse->per;
	double period_off = curr_per * pulse->per;

	if (0 <= t && t < pulse->td) {
		return pulse->i1;
	}
	else if ((pulse->td + period_off) <= t && t < (pulse->td + pulse->tr + period_off)) {
		/* Find the line equation */
		return pulse->i1 + ((pulse->i2 - pulse->i1) / pulse->tr) * (t - (curr_per * pulse->per + pulse->td));
	}
	else if ((pulse->td + pulse->tr + period_off) <= t && t < (pulse->td + pulse->tr + pulse->pw + period_off)) {
		return pulse->i2;
	}
	else if ((pulse->td + pulse->tr + pulse->pw + period_off) <= t && t < (pulse->td + pulse->tr + pulse->pw + pulse->tf + period_off)) {
		/* Find the line equation */
		return pulse->i2 + ((pulse->i1 - pulse->i2) / pulse->tf) * (t - (pulse->td + curr_per * pulse->per) - pulse->tr - pulse->pw);
	}
	else {
		return pulse->i1;
	}
}

/* Evaluates the pwl transient at given time t */
double eval_pwl(pwl_t *pwl, double t) {
	/*  t1 > 0, 0 <= t <= t1 */
	if (pwl->t[0] > 0 && t < pwl->t[0]) {
		return pwl->i[0];
	}
	else {
		for (int i = 0; i < pwl->n; i++) {
			if (pwl->t[i] > t) {
				/* Find the line equation */
				return pwl->i[i - 1] + (pwl->i[i] - pwl->i[i - 1]) / (pwl->t[i] - pwl->t[i - 1]) * (t - pwl->t[i - 1]);
			}
		}
		/* tn < fin_time, tn <= t <= fin_time */
		return pwl->i[pwl->n - 1];
	}
}

/* Creates and opens output files for every node included in the current TRAN analysis */
void create_tr_out_files(FILE *files[], tr_analysis_t tr_analysis) {
	/* Set the prefix name for the files */
	char prefix[] = "tr_analysis_V(";
	char file_name[MAX_FILE_NAME];

	/* Open different files for each node in plot/print array */
	for (int j = 0; j < tr_analysis.num_nodes; j++) {
		/* Temp buffer */
		char name[10];
		/* Construct the file name */
		strcpy(file_name, prefix);
		strcat(file_name, tr_analysis.nodes[j]);
		strcat(file_name, ")_");
		sprintf(name, "%g", tr_analysis.time_step);
		strcat(file_name, name);
		strcat(file_name, "_");
		sprintf(name, "%g", tr_analysis.fin_time);
		strcat(file_name, name);
		strcat(file_name, ".txt");
		/* Open the output file */
		files[j] = fopen(file_name, "w");
		if (files[j] == NULL) {
			fprintf(stderr, "Error opening file: %s\n", strerror(errno));
			exit(EXIT_FAILURE);
		}
		fprintf(files[j], "%-30s%-30s\n", "Time (seconds)", "Value (voltage)");
	}
}

/* Writes the output of the current step of the TRAN analysis to the output files */
void write_tr_out_files(FILE *files[], tr_analysis_t tr_analysis, hash_table_t *hash_table, double *sol_x, int step) {
	for (int j = 0; j < tr_analysis.num_nodes; j++) {
		int offset = ht_get_id(hash_table, tr_analysis.nodes[j]) - 1;
		fprintf(files[j], "%-30lf%-30lf\n", step * tr_analysis.time_step, sol_x[offset]);
	}
}