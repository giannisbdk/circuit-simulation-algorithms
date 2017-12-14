#include <math.h>
#include <string.h>

#include "transient_analysis.h"
#include "mna_dc.h"

/* Do transient analysis */
void tr_analysis(hash_table_t *hash_table, mna_system_t *mna, parser_t *parser, double *init_sol, double *sol_x) {
	/* Set the prefix name for the files */
	char prefix[] = "tr_analysis_";
	char file_name[MAX_FILE_NAME];

	FILE *gnuplot = popen("gnuplot", "w");
    fprintf(gnuplot, "plot '-'\n");

	/* Clear the decomposition flag for the mna system */
	mna->is_decomp = false;
	mna->tr_analysis_init = true;
	/* b is the RHS of the trapezoidal/backward euler */
	double *prev_sol = init_vector(mna->dimension);
	double *curr_response = init_vector(mna->dimension);
	double *prev_response;
	if (parser->options->TR) {
		prev_response = init_vector(mna->dimension);
	}

	for (int i = 0; i < parser->netlist->tr_counter; i++) {
		memcpy(prev_sol, init_sol, mna->dimension * sizeof(double));
		if (parser->options->TR) {
			memcpy(prev_response, mna->matrix->resp->value, mna->dimension * sizeof(double));
		}
		FILE *files[parser->tr_analysis[i].num_nodes];
        /* Open different files for each node in plot/print array */
        for (int j = 0; j < parser->tr_analysis[i].num_nodes; j++) {
        	/* Temp buffer */
        	char name[10];
            /* Construct the file name */
            strcpy(file_name, prefix);
            strcat(file_name, parser->tr_analysis[i].nodes[j]);
            strcat(file_name, "_");
    		sprintf(name, "%g", parser->tr_analysis[i].time_step);
            strcat(file_name, name);
            strcat(file_name, "_");
    		sprintf(name, "%g", parser->tr_analysis[i].fin_time);
            strcat(file_name, name);
            strcat(file_name, ".txt");
            /* Open the output file */
            files[j] = fopen(file_name, "w");
            if (files[j] == NULL) {
                fprintf(stderr, "Error opening file: %s\n", strerror(errno));
                exit(EXIT_FAILURE);
            }
            //TODO perhaps add fin_time and time_step inside file
            fprintf(files[j], "%-15s%-15s\n", "Time", "Value");
        }
		/* Find how many steps are required */
		int n_steps = parser->tr_analysis->fin_time / parser->tr_analysis->time_step;
		for (int step = 0; step <= n_steps; step++) {
			if (parser->options->TR) {
				set_trapezoidal_rhs(mna, curr_response, prev_response, prev_sol, parser->tr_analysis->time_step, step);
			}
			else {
				set_backward_euler_rhs(mna, curr_response, prev_sol, parser->tr_analysis->time_step, step);
			}
			/* Solve the system */
			solve_mna_system(mna, &sol_x, parser->options);
            for (int j = 0; j < parser->tr_analysis[i].num_nodes; j++) {
                int offset = ht_get_id(hash_table, parser->tr_analysis[i].nodes[j]) - 1;
                fprintf(files[j], "%-15lf%-15lf\n", step * parser->tr_analysis->time_step, sol_x[offset]);
			    fprintf(gnuplot, "%g %g\n", step * parser->tr_analysis->time_step, sol_x[offset]);
            }
			/* Copy current solution to prev to use for next iteration */
			memcpy(prev_sol, sol_x, mna->dimension * sizeof(double));
		}
		/* Close the file descriptors for the current transient analysis */
        for (int j = 0; j < parser->dc_analysis[i].num_nodes; j++) {
            fclose(files[j]);
        }
	}
	mna->tr_analysis_init = false;
	printf("Transient Analysis....... OK\n");
	fprintf(gnuplot, "e\n");
	fflush(gnuplot);
}

/* Computes and returns the right hand side of the trapezoidal */
/* h is time_step and k is the iteration */
void set_trapezoidal_rhs(mna_system_t *mna, double *curr_response, double *prev_response, double *prev_sol, double h, int k) {
	/* curr_response is e(tk) and prev_response is e(tk-1) */
	/* Set the values of the e(tk) vector */
	for (int i = 0; i < mna->dimension; i++) {
		list1_t *curr = mna->matrix->resp->nodes[i];
		if (curr->trans_spec != NULL) {
			switch (curr->trans_spec->type) {
				case EXP:
					curr_response[i] = eval_exp(curr->trans_spec->exp, h * k);
					break;
				case SIN:
					curr_response[i] = eval_sin(curr->trans_spec->sin, h * k);
					break;
				case PULSE:
					curr_response[i] = eval_pulse(curr->trans_spec->pulse, h * k);
					// printf("PULSE IS %lf\n\n\n", curr_response[i]);
					break;
				case PWL:
					curr_response[i] = eval_pwl(curr->trans_spec->pwl, h * k);
					break;
			}
		}
		else {
			/* If transient spec is absent set the DC value */
			curr_response[i] = mna->matrix->resp->value[i];
		}
	}
	/* Compute: e(tk) + e(tk-1) - sGhC*x(tk-1) and save it to the provided b vector */
	double *response_add = init_vector(mna->dimension);
	double *sGhc_x = init_vector(mna->dimension);
	add_vector(response_add, curr_response, prev_response, mna->dimension);
	mat_vec_mul(sGhc_x, mna->matrix->sGhC, prev_sol, mna->dimension);
	sub_vector(mna->b, response_add, sGhc_x, mna->dimension);
	memcpy(prev_response, curr_response, mna->dimension * sizeof(double));
}

void set_backward_euler_rhs(mna_system_t *mna, double *curr_response, double *prev_sol, double h, int k) {
	/* curr_response is e(tk) and prev_response is e(tk-1) */
	/* Set the values of the e(tk) vector */
	for (int i = 0; i < mna->dimension; i++) {
		list1_t *curr = mna->matrix->resp->nodes[i];
		if (curr->trans_spec != NULL) {
			switch (curr->trans_spec->type) {
				case EXP:
					curr_response[i] = eval_exp(curr->trans_spec->exp, h * k);
					break;
				case SIN:
					curr_response[i] = eval_sin(curr->trans_spec->sin, h * k);
					break;
				case PULSE:
					curr_response[i] = eval_pulse(curr->trans_spec->pulse, h * k);
					break;
				case PWL:
					curr_response[i] = eval_pwl(curr->trans_spec->pwl, h * k);
					break;
			}
		}
		else {
			/* If transient spec is absent set the DC value */
			curr_response[i] = mna->matrix->resp->value[i];
		}
	}
	/* Compute: e(tk) + (1/h)C*x(tk-1) and save it to the provided b vector */
	double *hC_x = init_vector(mna->dimension);
	mat_vec_mul(hC_x, mna->matrix->hC, prev_sol, mna->dimension);
	add_vector(mna->b, curr_response, hC_x, mna->dimension);
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
	// printf("curr_time is %lf\n", t);
	// printf("curr_period is %d\n", curr_per);
	// printf("period_off is %lf\n", period_off);

	if (0 <= t && t < pulse->td) {
		// printf("in 1st if %lf - %lf\n", t, pulse->td);
		return pulse->i1;
	}
	else if ((pulse->td + period_off) <= t && t < (pulse->td + pulse->tr + period_off)) {
		// printf("in 2nd if %lf - %lf\n", (pulse->td + period_off), (pulse->td + pulse->tr + period_off));
		/* Find the line equation */
		// printf("i1: %lf, i2: %lf, tr: %lf, t: %lf, td: %lf\n", pulse->i1, pulse->i2, pulse->tr, t, pulse->td);
		return pulse->i1 + ((pulse->i2 - pulse->i1) / pulse->tr) * (t - (curr_per * pulse->per + pulse->td));
	}
	else if ((pulse->td + pulse->tr + period_off) <= t && t < (pulse->td + pulse->tr + pulse->pw + period_off)) {
		// printf("in 3rd if %lf - %lf\n", (pulse->td + pulse->tr + period_off), (pulse->td + pulse->tr + pulse->pw + period_off));
		return pulse->i2;
	}
	else if ((pulse->td + pulse->tr + pulse->pw + period_off) <= t && t < (pulse->td + pulse->tr + pulse->pw + pulse->tf + period_off)) {
		// printf("in 4th if %lf - %lf\n", (pulse->td + pulse->tr + pulse->pw + period_off), (pulse->td + pulse->tr + pulse->pw + pulse->tf + period_off));
		return pulse->i2 + ((pulse->i1 - pulse->i2) / pulse->tf) * (t - (pulse->td + curr_per * pulse->per) - pulse->tr - pulse->pw);
	}
	else {
		// printf("in 5th if %lf - %lf\n", pulse->td+pulse->tr+pulse->pw+pulse->tf+period_off, pulse->td+pulse->per+period_off);
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
		return pwl->i[pwl->n];
	}
}