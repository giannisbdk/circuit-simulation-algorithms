#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include "parser.h"
#include "list.h"
#include "hash_table.h"
#include "mna_dc.h"
#include "parser.h"

#define HASH_TABLE_SIZE 65536
#define DC_ANALYSIS_NUM 25
#define MAX_FILE_NAME   50

int errno;

int main(int argc, char *argv[]) {
	FILE *file_input;
    char *line = NULL;
    int num_tokens = 0;
    size_t len = 0;
    ssize_t read;
    char **tokens;
    int num_g2_elem = 0;

    options_t options;
    
    /* Array to hold .DC options */
    dc_analysis_t dc_analysis[DC_ANALYSIS_NUM];
    int dc_cnt = 0;

    index_t *index = init_lists();
    hash_table_t *hash_table = ht_create(HASH_TABLE_SIZE);

    if (argc < 2) {
        printf("You must specify input file from cmd arguments.\nExiting....\n");
        exit(EXIT_FAILURE);
    }

    printf("Input file is: %s\n", argv[1]);

    file_input = fopen(argv[1], "rb");
    if (file_input == NULL) {
        exit(EXIT_FAILURE);
    }

    init_options(&options);

    //TODO define all the below to the parser
    //TODO parser should return num_g2_elem
    //TODO num_tokens is redundant we have it stored at &tokens[0][0]
    while ((read = getline(&line, &len, file_input)) != -1) {
    	tokens = tokenizer(line, &num_tokens);
    	if (tokens == NULL) {
    		continue;
    	}
        if (tokens[1][0] == '.') {
            if (strcmp(".OPTIONS", &tokens[1][0]) == 0) {
                for (int i = 2; i <= num_tokens; i++) {
                    if (strcmp("SPD", &tokens[i][0]) == 0) {
                        options.SPD = true;
                    }
                    if (strcmp("ITER", &tokens[i][0]) == 0) {
                        options.ITER = true;
                    }
                    if (strncmp("ITOL", &tokens[i][0], 4) == 0) {
                        sscanf((&tokens[i][0]) + 5, "%lf", &options.itol);
                    }
                }
            }
            else if (strcmp(".DC", &tokens[1][0]) == 0) {
                dc_analysis[dc_cnt].volt_source = (char *)malloc(strlen(&tokens[2][0]) * sizeof(char));
                assert(dc_analysis[dc_cnt].volt_source != NULL);
                sscanf(tokens[2], "%s", dc_analysis[dc_cnt].volt_source);
                sscanf(tokens[3], "%lf", &dc_analysis[dc_cnt].start);
                sscanf(tokens[4], "%lf", &dc_analysis[dc_cnt].end);
                sscanf(tokens[5], "%lf", &dc_analysis[dc_cnt].increment);
            }
            else if (strcmp(".PLOT", &tokens[1][0]) == 0 || strcmp(".PRINT", &tokens[1][0]) == 0) {
                dc_analysis[dc_cnt].nodes = (char **)malloc(num_tokens * sizeof(char *));
                assert(dc_analysis[dc_cnt].nodes != NULL);
                dc_analysis[dc_cnt].num_nodes = 0;
                for (int i = 2; i <= num_tokens; i++) {
                    /* Allocate memory for the node name ommiting the parentheses and the V */
                    dc_analysis[dc_cnt].nodes[i-2] = (char *)malloc((strlen(tokens[i]) - 3) * sizeof(char));
                    assert(dc_analysis[dc_cnt].nodes[i-2] != NULL);
                    /* Strip V and the parentheses around node name */
                    strncpy(dc_analysis[dc_cnt].nodes[i-2], tokens[i] + 2, (strlen(tokens[i]) - 3));
                    dc_analysis[dc_cnt].num_nodes++;
                }
                dc_cnt++;
            }
        }
        else {
            if (add_to_list(index, tokens, hash_table) == FAILURE) {
                exit(EXIT_FAILURE);
            } 
            if (tokens[1][0] == 'V' || tokens[1][0] == 'v' || tokens[1][0] == 'L' || tokens[1][0] == 'l') {
                num_g2_elem++;
            }
        }
       
        /* Free all the memory we allocated */
        for (int i = 0; i < num_tokens; i++) {
            free(tokens[i]);
        }
        free(tokens);
    }

#ifdef DEBUGL
    printf("Printing the lists\n");
    print_lists(index, hash_table);
#endif
    printf("Finished parsing %d circuit elements.\n", index->size1 + index->size2);

    mna_system_t *mna;

    int num_nodes = hash_table->seq - 1;
    int size = num_nodes + num_g2_elem;
    printf("\nsize: %d\nnum_nodes(w/o ground): %d\nnum_branches_g2: %d\n\n", size, num_nodes, num_g2_elem);
    
    /* Initialize the MNA_system */
    mna = init_mna_system(num_nodes, num_g2_elem);
    create_mna_system(mna, index, hash_table, num_nodes);
    print_mna_system(mna);
    double *sol_x = (double *)calloc(size, sizeof(double));
    solve_mna_system(mna, &sol_x, &options);
    printf("Solution of the MNA system:\n\n");
    print_vector(sol_x, size);
    printf("\n");

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
    fclose(file_input);
    if (line) {
        free(line);
    }

    char prefix[] = "dc_analysis_";
    char file_name[MAX_FILE_NAME];
    /* Cycle through dc analyisis targets */
    for (int i = 0; i < dc_cnt; i++) {
        list1_t *curr;
        for (curr = index->head1; curr != NULL; curr = curr->next) {
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
                //TODO Add a method in mna_dc.c to set the vector, so that we don't copy-pate the below
                for (int step = 0; step <= n_steps; step++) {
                    if (dc_analysis[i].volt_source[0] == 'V' || dc_analysis[i].volt_source[0] == 'v') {
                        mna->b[volt_indx] = val;
                    }
                    else if (dc_analysis[i].volt_source[0] == 'I' || dc_analysis[i].volt_source[0] == 'i') {
                        if (probe1_id == 0) {
                            mna->b[probe2_id - 1] = val;
                        }
                        else if (probe2_id == 0) {
                            mna->b[probe1_id - 1] = -val;
                        }
                        else {
                            mna->b[probe1_id - 1] = -val;
                            mna->b[probe2_id - 1] =  val;
                        }
                    }
                    /* Solve the system */
                    solve_mna_system(mna, &sol_x, &options);
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
                break;
            }
        }
    }
    /* Free all the dynamic allocated memory */
    free_index(&index);
    free_mna_system(&mna);
    ht_free(&hash_table);
	return 0;
}