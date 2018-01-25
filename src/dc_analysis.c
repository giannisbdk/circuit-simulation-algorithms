#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "dc_analysis.h"

/* DC Operating Point Analysis and prints the output to a file */
void dc_operating_point(hash_table_t *hash_table, double *sol_x) {
    entry_t *curr;
    double value;
    int id;
    /* DC Operating Point file */
    FILE *file_out = fopen("dc_operating_point.txt", "w");
    if (file_out == NULL) {
        fprintf(stderr, "Error opening file: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
    fprintf(file_out, "%-30s%-30s\n", "Node", "Value (voltage)");
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
            fprintf(file_out, "%-30s%-30lf\n", curr->key, value);
        }
    }
    fclose(file_out);
    printf("DC Operating Point.......OK\n");
}

/* DC Sweep analysis and outputs the results to a file */
void dc_sweep(list1_t *head, hash_table_t *hash_table, mna_system_t *mna, parser_t *parser, double *sol_x) {
    /* Cycle through dc analyisis targets */
    for (int i = 0; i < parser->netlist->dc_counter; i++) {
        list1_t *curr;
        for (curr = head; curr != NULL; curr = curr->next) {
            /* Find the voltage source for the analysis */
            if (strcmp(parser->dc_analysis[i].volt_source, curr->element) == 0) {
                /* Create an array with files for every node and create/open them */
                FILE *files[parser->dc_analysis[i].num_nodes];
                create_dc_out_files(files, parser->dc_analysis[i]);

                /* Run the DC analysis with the step */
                double value  = parser->dc_analysis[i].start;
                int n_steps   = (parser->dc_analysis[i].end - parser->dc_analysis[i].start) / parser->dc_analysis[i].increment;
                int volt_indx = g2_elem_indx(mna->g2_indx, mna->num_nodes, mna->num_g2_elem, parser->dc_analysis[i].volt_source);
                int probe1_id = ht_get_id(hash_table, curr->probe1);
                int probe2_id = ht_get_id(hash_table, curr->probe2);

                /* We need to zero out the sol_x vector from the previous operating point analysis value */
                int size = parser->netlist->num_nodes + parser->netlist->num_g2_elem;
                zero_out_vec(sol_x, size);

                //TODO Add a method in mna_dc.c to set the vector, so that we don't copy-pate the below
                for (int step = 0; step <= n_steps; step++) {
                    if (parser->dc_analysis[i].volt_source[0] == 'V' || parser->dc_analysis[i].volt_source[0] == 'v') {
                        mna->b[volt_indx] = value;
                    }
                    else if (parser->dc_analysis[i].volt_source[0] == 'I' || parser->dc_analysis[i].volt_source[0] == 'i') {
                        if (probe1_id == 0) {
                            mna->b[probe2_id - 1] = value;
                        }
                        else if (probe2_id == 0) {
                            mna->b[probe1_id - 1] = -value;
                        }
                        else {
                            mna->b[probe1_id - 1] = -value;
                            mna->b[probe2_id - 1] =  value;
                        }
                    }
                    /* Solve the system */
                    solve_mna_system(mna, &sol_x, NULL, parser->options);
                    /* DC analysis output to every file */
                    write_dc_out_files(files, parser->dc_analysis[i], hash_table, sol_x, value);
                    /* Increment value of voltage according to the increment step */
                    value += parser->dc_analysis[i].increment;
                }
                /* Close the file descriptors for the current dc analysis */
                for (int j = 0; j < parser->dc_analysis[i].num_nodes; j++) {
                    fclose(files[j]);
                }
                /* We found the voltage element and we did the DC sweep. That means we have to stop the iteration through the list */
                break;
            }
        }
    }
    if (parser->netlist->dc_counter) {
        printf("DC Sweep.................OK\n");
    }
}

/* Creates and opens output files for every node included in the current DC sweep analysis */
void create_dc_out_files(FILE *files[], dc_analysis_t dc_analysis) {
    /* Set the prefix name for the files */
    char prefix[] = "dc_analysis_";
    char file_name[MAX_FILE_NAME];
    /* Open different files for each node in plot/print array */
    for (int j = 0; j < dc_analysis.num_nodes; j++) {
        /* Construct the file name */
        strcpy(file_name, prefix);
        strcat(file_name, dc_analysis.volt_source);
        strcat(file_name, "_V(");
        strcat(file_name, dc_analysis.nodes[j]);
        strcat(file_name, ").txt");
        /* Open the output file */
        files[j] = fopen(file_name, "w");
        if (files[j] == NULL) {
            fprintf(stderr, "Error opening file: %s\n", strerror(errno));
            exit(EXIT_FAILURE);
        }
        fprintf(files[j], "%-30s%-30s\n", "Step (voltage)", "Value (voltage)");
    }
}

/* Writes the output of the current step of the DC analysis to the output files */
void write_dc_out_files(FILE *files[], dc_analysis_t dc_analysis, hash_table_t *hash_table, double *sol_x, double value) {
    int offset;
    for (int j = 0; j < dc_analysis.num_nodes; j++) {
        offset = ht_get_id(hash_table, dc_analysis.nodes[j]) - 1;
        fprintf(files[j], "%-30lf%-30lf\n", value, sol_x[offset]);
    }
}