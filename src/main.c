#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parser.h"
#include <errno.h>

#include "list.h"
#include "hash_table.h"
#include "mna_dc.h"
#include "parser.h"

#define HASH_TABLE_SIZE 65536
#define DC_ANALYSIS_NUM 25

int errno;

int main(int argc, char *argv[]) {

	FILE *file_input;
    char *line = NULL;
    int num_tokens = 0;
    size_t len = 0;
    ssize_t read;
    char **tokens;
    int num_branches = 0;

    options_t options;
    /* Array to hold .DC options */
    dc_analysis_t dc_analysis[DC_ANALYSIS_NUM];
    int dc_analysis_cnt = 0;

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
    //TODO parser should return num_branches
    //TODO num_tokens is redundant we have it stored at &tokens[0][0]
    while((read = getline(&line, &len, file_input)) != -1) {

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
                }
            }
            else if (strcmp(".DC", &tokens[1][0]) == 0) {
                dc_analysis[dc_analysis_cnt].volt_source = 
                    (char *)malloc(strlen(&tokens[2][0]) * sizeof(char));
                //TODO change strcpy to sscanf in list.h
                sscanf(tokens[2], "%s", dc_analysis[dc_analysis_cnt].volt_source);
                sscanf(tokens[3], "%lf", &dc_analysis[dc_analysis_cnt].start);
                sscanf(tokens[4], "%lf", &dc_analysis[dc_analysis_cnt].end);
                sscanf(tokens[5], "%lf", &dc_analysis[dc_analysis_cnt].increment);
            }
            else if (strcmp(".PLOT", &tokens[1][0]) == 0 ||
                     strcmp(".PRINT", &tokens[1][0]) == 0) {
                dc_analysis[dc_analysis_cnt].nodes = (char **)malloc(num_tokens * sizeof(char *));
                dc_analysis[dc_analysis_cnt].num_nodes = 0;
                for (int i = 2; i <= num_tokens; i++) {
                    /* Allocate memory for the node name ommiting the parentheses and the V */
                    dc_analysis[dc_analysis_cnt].nodes[i-2] = (char *)malloc((strlen(tokens[i]) - 3) * sizeof(char));
                    /* Strip V and the parentheses around node name */
                    strncpy(dc_analysis[dc_analysis_cnt].nodes[i-2], tokens[i] + 2, (strlen(tokens[i]) - 3));
                    dc_analysis[dc_analysis_cnt].num_nodes++;
                }
                dc_analysis_cnt++;
            }
        }
        else {
            if (add_to_list(index, tokens, hash_table) == FAILURE) {
                exit(EXIT_FAILURE);
            } 
            if (tokens[1][0] == 'V' || tokens[1][0] == 'v' || 
                tokens[1][0] == 'L' || tokens[1][0] == 'l') {
                num_branches++;
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
    int size = num_nodes + num_branches;
    printf("\nsize: %d\nnum_nodes(w/o ground): %d\nnum_branches_g2: %d\n\n", size, num_nodes, num_branches);
    
    /* Initialize the MNA_system */
    mna = init_mna_system(size);
    create_mna_system(mna, index, hash_table, num_nodes);
    print_mna_system(mna);
    gsl_vector *sol_x = solve_mna_system(mna, options.SPD);
    printf("Solution of the MNA system:\n\n");
    print_vector(sol_x);

    FILE *file_out;

    /* DC Operating Point to file */
    file_out = fopen("DC_Opearting_Point.txt", "w");
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
            value = gsl_vector_get(sol_x, id);
            /* Output to the file */
            fprintf(file_out, "%-15s%-15lf\n", curr->key, value);
        }
    }
    fclose(file_out);
    fclose(file_input);
    if (line) {
        free(line);
    }

    char file_name[50] = "dc_analysis_";
    /* Cycle through dc analyisis targets */
    for (int i = 0; i < dc_analysis_cnt; i++) {
        list1_t *curr;
        for (curr = index->head1; curr != NULL; curr = curr->next) {
            /* Find the voltage source for the analysis */
            if (strcmp(dc_analysis[i].volt_source, curr->element) == 0) {
                // int probe1 = ht_get_id(hash_table, curr->probe1) - 1;
                // int probe2 = ht_get_id(hash_table, curr->probe2) - 1;
                int probe1 = ht_get_id(hash_table, curr->probe1);
                int probe2 = ht_get_id(hash_table, curr->probe2);
                /* Create an array with file names for every node */
                char *file_names[dc_analysis[i].num_nodes];
                FILE *files[dc_analysis[i].num_nodes];
                /* Open different files for each node in plot/print array */
                for (int j = 0; j < dc_analysis[i].num_nodes; j++) {
                    strcpy(file_names[j], file_name);
                    strcat(file_names[j], dc_analysis[i].nodes[j]);
                    strcat(file_names[j], ".txt");
                    /* Open the output file */
                    files[j] = fopen(file_names[j], "w");
                    if (files[j] == NULL) {
                        fprintf(stderr, "Error opening file: %s\n", strerror(errno));
                        exit(EXIT_FAILURE);
                    }
                    fprintf(files[j], "%-15s%-15s\n", "Step", "Value");
                }
                /* Run the DC analysis with the step */
                for (double val = dc_analysis[i].start; val <= dc_analysis[i].end; val += dc_analysis[i].increment) {
                    if (probe1 == 0) {
                        gsl_vector_set(mna->b, probe2-1, val);
                    }
                    else if (probe2 == 0) {
                        gsl_vector_set(mna->b, probe1-1, val);
                    }
                    else {
                        gsl_vector_set(mna->b, probe1-1, val);
                        gsl_vector_set(mna->b, probe2-1, val);
                    }
                    /* Solve the system */
                    sol_x = solve_mna_system(mna, options.SPD);
                    /* DC analysis to every file */
                    int offset;
                    for (int j = 0; j < dc_analysis[i].num_nodes; j++) {
                        offset = (ht_get_id(hash_table, dc_analysis[i].nodes[j]) - 1);
                        fprintf(files[j], "%-15lf%-15lf\n", val, gsl_vector_get(sol_x, offset));
                    }
                }
                // for (int j = 0; j < dc_analysis[i].num_nodes; j++) {
                //     // fclose(files[j]);
                // }
            }
        }
    }

    //TODO free the lists
    /* Free all the dynamic allocated memory */
    free_mna_system(mna);
    ht_free(hash_table);
	return 0;
}