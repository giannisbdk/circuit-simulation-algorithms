#ifndef DC_ANALYSIS_H
#define DC_ANALYSIS_H

#include "hash_table.h"
#include "mna.h"
#include "routines.h"

#define MAX_FILE_NAME 50

void dc_operating_point(hash_table_t *hash_table, double *sol_x);
void dc_sweep(list1_t *head, hash_table_t *hash_table, mna_system_t *mna, parser_t *parser, double *sol_x);

#endif
