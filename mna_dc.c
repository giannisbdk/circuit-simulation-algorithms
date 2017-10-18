#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mna_dc.h"

double **init_array(int row, int col) {

	double **array;;

	array = (double **) calloc(row, sizeof(double *));
	for(int i=0; i<row; i++) {
		array[i] = (double *) calloc(col, sizeof(double));
	}
	return array;
}

void create_mna_arrays(mna_arrays_t *mna, index_t *index, hash_table_t *hash_table, int offset) {

	list1_t *curr;
	double value;
	int volt_sources_cnt = 0; 

	for (curr = index->head1; curr != NULL; curr = curr->next) {

		int probe1_id = ht_get_id(hash_table, curr->probe1);
		int probe2_id = ht_get_id(hash_table, curr->probe2);
		int i = probe1_id - 1;
		int j = probe2_id - 1;
		// int i = curr->probe1 - 1;
		// int j = curr->probe2 - 1;
		if(curr->type == 'C' || curr->type == 'c') {
			continue;
		}
		else if(curr->type == 'L' || curr->type == 'l') {
			value = 0;
		}
		else if(curr->type == 'R' || curr->type == 'r') {
			value = 1 / curr->value;
			if(probe1_id == 0) {
				mna->left[j][j] += value;
			}
			else if(probe2_id == 0) {
				mna->left[i][i] += value;
			}
			else {
				mna->left[i][i] +=  value;
				mna->left[j][j] +=  value;
				mna->left[i][j] += -value;
				mna->left[j][i] += -value;
			}
		}
		else if(curr->type == 'V' || curr->type == 'v') {
			value = curr->value;
			if(probe1_id == 0) {
				mna->left[j][offset + volt_sources_cnt] = 1;
				mna->left[offset + volt_sources_cnt][j] = 1;
				mna->right[offset + volt_sources_cnt][0] += value;
			}
			else if(probe2_id == 0) {
				mna->left[i][offset + volt_sources_cnt] = 1;
				mna->left[offset + volt_sources_cnt][i] = 1;
				mna->right[offset + volt_sources_cnt][0] += value;
			}
			else {
				mna->left[i][offset + volt_sources_cnt] =  1;
				mna->left[j][offset + volt_sources_cnt] = -1;
				mna->left[offset + volt_sources_cnt][i] =  1;
				mna->left[offset + volt_sources_cnt][j] = -1;
				mna->right[offset + volt_sources_cnt][0] += value;
			}
			// Keep track of how many voltage sources we have already found
			volt_sources_cnt++;
		}
		else if(curr->type == 'I' || curr->type == 'i') {
			value = curr->value;
			if(probe1_id == 0) {
				mna->right[j][0] += value;
			}
			else if(probe2_id == 0) {
				mna->right[i][0] += value;
			}
			else {
				mna->right[i][0] +=  value;
				mna->right[j][0] +=  value;
			}
		}
	}
}

void print_mna_left(mna_arrays_t *mna, int size) {
	printf("MNA left array:\n\n");
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            printf("%lf\t", mna->left[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_mna_right(mna_arrays_t *mna, int size) {
	printf("MNA right array:\n\n");
	for (int i = 0; i < size; ++i) {
        for (int j = 0; j < 1; ++j) {
            printf("%lf\t", mna->right[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}