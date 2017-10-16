#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "list.h"
#include "hash_table.h"

index_t* init_lists() {
	index_t *index = (index_t *)malloc(sizeof(index_t));
	index->head1 = NULL;
	index->head2 = NULL;
	index->size1 = 0;
	index->size2 = 0;
	return index;
}

/* General method to add new node to the appropriate list */
int add_to_list(index_t *index, char **tokens, hash_table_t *hash_table) {

	char type = tokens[1][0];
	if(type == 'R' || type == 'C' || type == 'I' || type == 'V' || type == 'L' ||
		type == 'r' || type == 'c' || type == 'i' || type == 'v' || type == 'l') {
		index->head1 = add_to_list1(index->head1, index->size1, tokens, hash_table);
		if(index->head1 != NULL) {
			index->size1++;
		}
		else {
			return FAILURE;
		}
	}
	else if(type == 'M' || type == 'Q' || type == 'D' ||
		type == 'm' || type == 'q' || type == 'd') {
		index->head2 = add_to_list2(index->head2, index->size2, tokens, hash_table);
		if(index->head2 != NULL) {
			index->size2++;
		}
		else {
			return FAILURE;
		}
	}
	else {
		printf("Something is totally wrong with the input.\n");
		return FAILURE;
	}
	return SUCCESS;
}

list1_t *add_to_list1(list1_t *head, int size, char **tokens, hash_table_t *hash_table) {

	list1_t *new_node = (list1_t *)malloc(sizeof(list1_t));
	if(new_node == NULL) {
		printf("Something went wrong with malloc\n");
		return NULL;
	}
	/* In case the list is empty */
	if(size == 0) {
		new_node->prev = NULL;
		new_node->next = NULL;
		head = new_node;
	}
	else {
		list1_t *curr;
		/* Find the last node */
		for(curr = head; curr->next != NULL; curr = curr->next);
		curr->next = new_node;
		new_node->prev = curr;
		new_node->next = NULL;
	}
	/* Set new node struct fields */
	strncpy(&new_node->type, &tokens[1][0], 1);
	new_node->id = (char *)malloc(strlen(&tokens[1][1]) * sizeof(char));
	strcpy(new_node->id, &tokens[1][1]);
	new_node->probe1 = ht_add(hash_table, &tokens[2][0]);
	new_node->probe2 = ht_add(hash_table, &tokens[3][0]);
	if(new_node->probe1 == FAILURE || new_node->probe2 == FAILURE) {
		return NULL;
	}
	sscanf(tokens[4], "%Lf", &new_node->value);
	return head;
}

list2_t *add_to_list2(list2_t *head, int size, char **tokens, hash_table_t *hash_table) {

	//TODO error handling
	int num_tokens = atoi(tokens[0]);

	list2_t *new_node = (list2_t *)malloc(sizeof(list2_t));
	if(new_node == NULL) {
		return NULL;
	}
	/* In case the list is empty */
	if(size == 0) {
		new_node->prev = NULL;
		new_node->next = NULL;
		head = new_node;
	}
	else {
		list2_t *curr;
		/* Find the last node */
		for(curr = head; curr->next != NULL; curr = curr->next);
		curr->next = new_node;
		new_node->prev = curr;
		new_node->next = NULL;
	}

	/* Set new node struct fields */
	strncpy(&new_node->type, &tokens[1][0], 1);
	new_node->id = (char *)malloc(strlen(&tokens[1][1]) * sizeof(char));
	strcpy(new_node->id, &tokens[1][1]);
	if(new_node->type == 'D' || new_node->type == 'd') {
		new_node->probe1 = ht_add(hash_table, &tokens[2][0]);
		new_node->probe2 = ht_add(hash_table, &tokens[3][0]);
		if(new_node->probe1 == FAILURE || new_node->probe2 == FAILURE) {
			return NULL;
		}
		new_node->probe3 = -1;
		new_node->probe4 = -1;
		new_node->model_name = -1;
		new_node->length = new_node->width = -1;
		/* For now we ignore mode_name, thus we use 3 */
		new_node->area = num_tokens > 3 ? atoi(&tokens[4][0]) : 1;
	}
	else if(new_node->type == 'M' || new_node->type == 'm') {
		new_node->probe1 = ht_add(hash_table, &tokens[2][0]);
		new_node->probe2 = ht_add(hash_table, &tokens[3][0]);
		new_node->probe3 = ht_add(hash_table, &tokens[4][0]);
		new_node->probe4 = ht_add(hash_table, &tokens[5][0]);
		if(new_node->probe1 == FAILURE || new_node->probe2 == FAILURE ||
			new_node->probe3 == FAILURE || new_node->probe4 == FAILURE) {
			return NULL;
		}
		new_node->model_name = -1;
		sscanf(tokens[6], "%Lf", &new_node->length);
		sscanf(tokens[7], "%Lf", &new_node->width);
	}
	else if(new_node->type == 'Q' || new_node->type == 'q') {
		new_node->probe1 = ht_add(hash_table, &tokens[2][0]);
		new_node->probe2 = ht_add(hash_table, &tokens[3][0]);
		new_node->probe3 = ht_add(hash_table, &tokens[4][0]);
		if(new_node->probe1 == FAILURE || new_node->probe2 == FAILURE ||
			new_node->probe3 == FAILURE) {
			return NULL;
		}
		new_node->probe4 = -1;
		new_node->model_name = -1;
		new_node->length = new_node->width = -1;
		/* For now we ignore mode_name, thus we use 4 */
		new_node->area = num_tokens > 4 ? atoi(&tokens[5][0]) : 1;
	}
	else {
		printf("ERROR: Unknown circuit element\n");
		return 0;
	}
	return head;
}

void print_lists(index_t *index, hash_table_t *hash_table) {

	print_list1(index->head1, hash_table);
	print_list2(index->head2, hash_table);
}

void print_list1(list1_t *head, hash_table_t *hash_table) {

	list1_t *curr = head;

	while(curr != NULL) {
		printf("Type: %c\n", curr->type);
		printf("Id: %s\n", curr->id);
		printf("Probe1: %s\n", ht_get_value(hash_table, curr->probe1));
		printf("Probe2: %s\n", ht_get_value(hash_table, curr->probe2));
		printf("Value: %.20Lf\n", curr->value);
		printf("\n");
		curr = curr->next;
	}
}

void print_list2(list2_t *head, hash_table_t *hash_table) {

	list2_t *curr = head;

	while(curr != NULL) {
		printf("Type: %c\n", curr->type);
		printf("Id: %s\n", curr->id);
		printf("Probe1: %s\n", ht_get_value(hash_table, curr->probe1));
		printf("Probe2: %s\n", ht_get_value(hash_table, curr->probe2));
		printf("Probe3: %s\n", ht_get_value(hash_table, curr->probe3));
		printf("Probe4: %s\n", ht_get_value(hash_table, curr->probe4));
		printf("Model_name: %d\n", curr->model_name);
		printf("Area: %d\n", curr->area);
		printf("Length %.20Lf\nWidth: %.20Lf\n", curr->length, curr->width);
		printf("\n");
		curr = curr->next;
	}
}