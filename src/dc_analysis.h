#ifndef DC_ANALYSIS_H
#define DC_ANALYSIS_H

#include "hash_table.h"
#include "mna_dc.h"
#include "routines.h"

#define MAX_FILE_NAME 50

void dc_operating_point(hash_table_t *hash_table, double *sol_x);
void dc_sweep(list1_t *head, hash_table_t *hash_table, mna_system_t *mna, dc_analysis_t *dc_analysis,
			options_t *options, netlist_elem_t *netlist_elem, double *sol_x);

#endif
