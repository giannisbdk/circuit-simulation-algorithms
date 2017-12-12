#include <math.h>

#include "transient_analysis.h"

void tr_analysis(list1_t *head, hash_table_t *hash_table, mna_system_t *mna, parser_t *parser, double *sol_x) {
	/* Set the prefix name for the files */
    char prefix[] = "tr_analysis_";
    char file_name[MAX_FILE_NAME];

    if (!parser->options->BE) {
    	// trapezoidal(head, hash_table, mna, parser, sol_x);
    }
    else {
    	// backward_euler(head, hash_table, mna, parser, sol_x);
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
		return 	expon->i1 + (expon->i2 - expon->i1) *
				(exp(-(t - expon->td2) / expon->tc2) - exp(-(t - expon->td1) / expon->tc1));
	}
}

/* Evaluates the sinusodial transient at given time t */
double eval_sin(sin_t *sinus, double t) {
	if (0 <= t && t < sinus->td) {
		return sinus->i1 + sinus->ia * sin(2 * M_PI * sinus->ph / 360);
	}
	else {
		return	sinus->i1 + sinus->ia *
				sin(2 * M_PI * sinus->fr * (t - sinus->td) + 2 * M_PI * sinus->ph / 360) *
				exp(-(t - sinus->td) * sinus->df);
	}
}

/* Evaluates the pulse transient at given time t */
double eval_pulse(pulse_t *pulse, double t, int k) {
	double curr_per = k * pulse->per;
	if (0 <= t && t < pulse->td) {
		return pulse->i1;
	}
	else if ((pulse->td + curr_per) <= t && t < (pulse->td + pulse->tr + curr_per)) {
		/* Find the line equation */
		return pulse->i1 + ((pulse->i2 - pulse->i1) / pulse->tr) * (t - pulse->td);
	}
	else if ((pulse->td + pulse->tr + curr_per) <= t && t < (pulse->td + pulse->tr + pulse->pw + curr_per)) {
		return pulse->i2;
	}
	else if ((pulse->td + pulse->tr + pulse->pw + curr_per) <= t && t < (pulse->td + pulse->tr + pulse->pw + pulse->tf + curr_per)) {
		return pulse->i2 + ((pulse->i1 - pulse->i2) / pulse->tf) * (t - pulse->td + pulse->tr + pulse->pw);
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
		return pwl->i[pwl->n];
	}
}