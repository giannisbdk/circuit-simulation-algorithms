#ifndef TRANSIENT_ANALYSIS_H
#define TRANSIENT_ANALYSIS_H

#include "mna.h"
#include "list.h"
#include "routines.h"

#define MAX_FILE_NAME 50

void tr_analysis(index_t *index, hash_table_t *hash_table, mna_system_t *mna, parser_t *parser, double *dc_op, double *sol_x);
void set_trapezoidal_rhs(mna_system_t *mna, double *curr_response, double *prev_response, double *prev_sol, double h,
                         int k, bool SPARSE);
void set_backward_euler_rhs(mna_system_t *mna, double *curr_response, double *prev_sol, double h, int k, bool SPARSE);
void set_response_vector(double *response, resp_t *resp, double time, int dimension);
double eval_exp(exp_t *expon, double t);
double eval_sin(sin_t *sinus, double t);
double eval_pulse(pulse_t *pulse, double t);
double eval_pwl(pwl_t *pwl, double t);
void init_transient_state(mna_system_t *mna, parser_t *parser, int dimension);
void create_tr_out_files(FILE *files[], tr_analysis_t tr_analysis);
void write_tr_out_files(FILE *files[], tr_analysis_t tr_analysis, hash_table_t *hash_table, double *sol_x, int step);

#endif