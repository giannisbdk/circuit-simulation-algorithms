#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "dc_analysis.h"

void dc_operating_point(hash_table_t *hash_table, double *sol_x) {
    FILE *file_out;
    /* DC Operating Point to file */
    file_out = fopen("dc_opearting_point.txt", "w");
    if (file_out == NULL) {
        fprintf(stderr, "Error opening file: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
    fprintf(file_out, "%-15s%-15s\n", "Node", "Value");
    entry_t *curr;
    double value;
    int id;
    for (int i = 0; i < hash_table->size; i++) {
        for (curr = hash_table->table[i]; curr != NULL; curr = curr->next) {
            /* Get node name */
            id = curr->id;
            if (id == 0) {
                continue;
            }
            else {
                id -= 1;
            }
            /* Get the corresponding cell of the solution vector */
            value = sol_x[id];
            /* Output to the file */
            fprintf(file_out, "%-15s%-15lf\n", curr->key, value);
        }
    }
    fclose(file_out);
}

void dc_sweep(list1_t *head, hash_table_t *hash_table, mna_system_t *mna, dc_analysis_t *dc_analysis,
        options_t *options, netlist_elem_t *netlist_elem, double *sol_x) {
    char prefix[] = "dc_analysis_";
    char file_name[MAX_FILE_NAME];
    /* Cycle through dc analyisis targets */
    for (int i = 0; i < netlist_elem->dc_counter; i++) {
        list1_t *curr;
        for (curr = head; curr != NULL; curr = curr->next) {
            /* Find the voltage source for the analysis */
            if (strcmp(dc_analysis[i].volt_source, curr->element) == 0) {
                /* Create an array with file names for every node */
                FILE *files[dc_analysis[i].num_nodes];
                /* Open different files for each node in plot/print array */
                for (int j = 0; j < dc_analysis[i].num_nodes; j++) {
                    /* Construct the file name */
                    strcpy(file_name, prefix);
                    strcat(file_name, dc_analysis[i].volt_source);
                    strcat(file_name, "_");
                    strcat(file_name, dc_analysis[i].nodes[j]);
                    strcat(file_name, ".txt");
                    /* Open the output file */
                    files[j] = fopen(file_name, "w");
                    if (files[j] == NULL) {
                        fprintf(stderr, "Error opening file: %s\n", strerror(errno));
                        exit(EXIT_FAILURE);
                    }
                    fprintf(files[j], "%-15s%-15s\n", "Step", "Value");
                }
                /* Run the DC analysis with the step */
                int n_steps = (dc_analysis[i].end - dc_analysis[i].start) / dc_analysis[i].increment;
                double val  = dc_analysis[i].start;
                int volt_indx = g2_elem_indx(mna->g2_indx, mna->num_nodes, mna->num_g2_elem, dc_analysis[i].volt_source);
                int probe1_id = ht_get_id(hash_table, curr->probe1);
                int probe2_id = ht_get_id(hash_table, curr->probe2);
                /* We need to zero out the sol_x vector from the previous operating point analysis value */
                /* Important notice that we don't use memset because we are dealing with double */
                int size = netlist_elem->num_nodes + netlist_elem->num_g2_elem;
                zero_out_vec(sol_x, size);
                //TODO Add a method in mna_dc.c to set the vector, so that we don't copy-pate the below
                for (int step = 0; step <= n_steps; step++) {
                    if (dc_analysis[i].volt_source[0] == 'V' || dc_analysis[i].volt_source[0] == 'v') {
                        mna->matrix->b[volt_indx] = val;
                    }
                    else if (dc_analysis[i].volt_source[0] == 'I' || dc_analysis[i].volt_source[0] == 'i') {
                        if (probe1_id == 0) {
                            mna->matrix->b[probe2_id - 1] = val;
                        }
                        else if (probe2_id == 0) {
                            mna->matrix->b[probe1_id - 1] = -val;
                        }
                        else {
                            mna->matrix->b[probe1_id - 1] = -val;
                            mna->matrix->b[probe2_id - 1] =  val;
                        }
                    }
                    /* Solve the system */
                    solve_mna_system(mna, &sol_x, options);
                    /* DC analysis output to every file */
                    int offset;
                    for (int j = 0; j < dc_analysis[i].num_nodes; j++) {
                        offset = ht_get_id(hash_table, dc_analysis[i].nodes[j]) - 1;
                        fprintf(files[j], "%-15lf%-15lf\n", val, sol_x[offset]);
                    }
                    val += dc_analysis[i].increment;
                }
                /* Close the file descriptors for the current dc analysis */
                for (int j = 0; j < dc_analysis[i].num_nodes; j++) {
                    fclose(files[j]);
                }
                /* We found the voltage element and we did the DC sweep. That means we have to stop the iteration through the list */
                break;
            }
        }
    }
}