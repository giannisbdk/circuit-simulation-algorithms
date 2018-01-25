#ifndef AC_ANALYSIS_H
#define AC_ANALYSIS_H

#include "mna.h"
#include "routines.h"
#include "ac_spec.h"

#define MAX_FILE_NAME 50

void ac_analysis(index_t *index, hash_table_t *hash_table, mna_system_t *mna, parser_t *parser,
                 double *dc_op, gsl_vector_complex *sol_x);
void get_sweep_points(double *array, ac_analysis_t ac_analysis);
void lin_sweep(double *array, double start, double end, int points);
void log_sweep(double *array, double start, double end, int points);

#endif