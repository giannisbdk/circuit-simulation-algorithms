#ifndef TRANSIENT_ANALYSIS_H
#define TRANSIENT_ANALYSIS_H

#include "list.h"
#include "mna_dc.h"
#include "routines.h"

#define MAX_FILE_NAME 50

void tr_analysis(hash_table_t *hash_table, mna_system_t *mna, parser_t *parser, double *init_sol, double *sol_x);
void set_trapezoidal_rhs(mna_system_t *mna, double *curr_response, double *prev_response, double *prev_sol, double h, int k);
double eval_exp(exp_t *expon, double t);
double eval_sin(sin_t *sinus, double t);
double eval_pulse(pulse_t *pulse, double t);
double eval_pwl(pwl_t *pwl, double t);

#endif