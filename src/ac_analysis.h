#ifndef AC_ANALYSIS_H
#define AC_ANALYSIS_H

#include "mna.h"
#include "routines.h"
#include "ac_spec.h"

#define MAX_FILE_NAME 50

void ac_analysis(index_t *index, hash_table_t *hash_table, mna_system_t *mna, parser_t *parser, double *sol_x);
double get_freq_step(ac_analysis_t ac_analysis);

#endif