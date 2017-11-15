#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "list.h"

index_t *init_lists() {
	/* Allocate memory for the index that stores the lists */
	index_t *index = (index_t *)malloc(sizeof(index_t));
	assert(index != NULL);
	/* Initiliaze the fields */
	index->head1 = index->tail1 = NULL;
	index->head2 = index->tail2 = NULL;
	index->size1 = index->size2 = 0;
	return index;
}

int add_to_list(index_t *index, char **tokens, hash_table_t *hash_table) {
	int res;
	/* Get the type of the circuit element */
	char type = tokens[1][0];
	if (type == 'R' || type == 'C' || type == 'I' || type == 'V' || type == 'L' ||
		type == 'r' || type == 'c' || type == 'i' || type == 'v' || type == 'l') {
		res = add_to_list1(index, tokens, hash_table);
		if (res == SUCCESS) {
			index->size1++;
		}
		else {
			return FAILURE;
		}
	}
	else if (type == 'M' || type == 'Q' || type == 'D' ||
		type == 'm' || type == 'q' || type == 'd') {
		res = add_to_list2(index, tokens, hash_table);
		if (res == SUCCESS) {
			index->size2++;
		}
		else {
			return FAILURE;
		}
	}
	else {
		printf("Circuit element has wrong type %c\n", type);
		return FAILURE;
	}
	return SUCCESS;
}

/* Add new node to the end of the list */
int add_to_list1(index_t *index, char **tokens, hash_table_t *hash_table) {
	list1_t *new_node = (list1_t *)malloc(sizeof(list1_t));
	if (new_node == NULL) {
		printf("Could not allocate memory for new node!\n");
		return FAILURE;
	}
	/* In case list is empty */
	if (index->size1 == 0) {
		new_node->prev = NULL;
		new_node->next = NULL;
		index->head1 = new_node;
		index->tail1 = new_node;
#ifdef DEBUGL
		printf("Adding node to the head of list1, new_node: %p, head1: %p, tail1: %p\n", new_node, index->head1, index->tail1);
#endif
	}
	/* If it is not empty add element to the end of the list */
	else {
		new_node->prev = index->tail1;
		index->tail1->next = new_node;
		index->tail1 = new_node;
		new_node->next = NULL;
#ifdef DEBUGL
		printf("Adding node to the end of list1, new_node: %p, new_node->prev: %p\n", new_node, new_node->prev);
#endif
	}
	/* Set new node struct fields */
	strncpy(&new_node->type, &tokens[1][0], 1);
	new_node->element = (char *)malloc(strlen(&tokens[1][0]) * sizeof(char));
	new_node->probe1  = (char *)malloc(strlen(&tokens[2][0]) * sizeof(char));
	new_node->probe2  = (char *)malloc(strlen(&tokens[3][0]) * sizeof(char));
	assert(new_node->element != NULL);
	assert(new_node->probe1 != NULL);
	assert(new_node->probe2 != NULL);
	strcpy(new_node->element, &tokens[1][0]);
	strcpy(new_node->probe1,  &tokens[2][0]);
	strcpy(new_node->probe2,  &tokens[3][0]);
	/* Add probes to the hashtable */
	// TODO perhaps replace &tokens with new_node probes
	ht_set(hash_table, &tokens[2][0]);
	ht_set(hash_table, &tokens[3][0]);
	sscanf(tokens[4], "%Lf", &new_node->value);
	return SUCCESS;
}

/* Add new node to the end of the list */
int add_to_list2(index_t *index, char **tokens, hash_table_t *hash_table) {
	int num_tokens = atoi(tokens[0]);
	list2_t *new_node = (list2_t *)malloc(sizeof(list2_t));
	if (new_node == NULL) {
		printf("Could not allocate memory for new node!\n");
		return FAILURE;
	}
	/* In case list is empty */
	if (index->size2 == 0) {
		new_node->prev = NULL;
		new_node->next = NULL;
		index->head2 = new_node;
		index->tail2 = new_node;
#ifdef DEBUGL
		printf("Adding node to the head of list2, new_node: %p, head1: %p, tail1: %p\n", new_node, index->head1, index->tail1);
#endif
	}
	/* If it is not empty add element to the end of the list */
	else {
		new_node->prev = index->tail2;
		index->tail2->next = new_node;
		index->tail2 = new_node;
		new_node->next = NULL;
#ifdef DEBUGL
		printf("Adding node to the end of list2, new_node: %p, new_node->prev: %p\n", new_node, new_node->prev);
#endif
	}
	/* Set new node struct fields */
	strncpy(&new_node->type, &tokens[1][0], 1);
	new_node->element = (char *)malloc(strlen(&tokens[1][0]) * sizeof(char));
	assert(new_node->element != NULL);
	strcpy(new_node->element, &tokens[1][0]);
	/* Init probes to NULL */
	new_node->probe1 = NULL;
	new_node->probe2 = NULL;
	new_node->probe3 = NULL;
	new_node->probe4 = NULL;
	if (new_node->type == 'D' || new_node->type == 'd') {
		new_node->probe1 = (char *)malloc(strlen(&tokens[2][0]) * sizeof(char));
		new_node->probe2 = (char *)malloc(strlen(&tokens[2][0]) * sizeof(char));
		assert(new_node->probe1 != NULL);
		assert(new_node->probe2 != NULL);
		strcpy(new_node->probe1, &tokens[2][0]);
		strcpy(new_node->probe1, &tokens[3][0]);
		/* Add probes to the hashtable */
		ht_set(hash_table, &tokens[2][0]);
		ht_set(hash_table, &tokens[3][0]);
		new_node->model_name = -1;
		new_node->length = new_node->width = -1;
		/* For now we ignore model_name, thus we use 3 */
		new_node->area = num_tokens > 3 ? atoi(&tokens[4][0]) : 1;
	}
	else if (new_node->type == 'M' || new_node->type == 'm') {
		new_node->probe1 = (char *)malloc(strlen(&tokens[2][0]) * sizeof(char));
		new_node->probe2 = (char *)malloc(strlen(&tokens[3][0]) * sizeof(char));
		new_node->probe3 = (char *)malloc(strlen(&tokens[4][0]) * sizeof(char));
		new_node->probe4 = (char *)malloc(strlen(&tokens[5][0]) * sizeof(char));
		assert(new_node->probe1 != NULL);
		assert(new_node->probe2 != NULL);
		assert(new_node->probe3 != NULL);
		assert(new_node->probe4 != NULL);
		strcpy(new_node->probe1, &tokens[2][0]);
		strcpy(new_node->probe1, &tokens[3][0]);
		strcpy(new_node->probe1, &tokens[4][0]);
		strcpy(new_node->probe1, &tokens[5][0]);
		/* Add probes to the hashtable */
		ht_set(hash_table, &tokens[2][0]);
		ht_set(hash_table, &tokens[3][0]);
		ht_set(hash_table, &tokens[4][0]);
		ht_set(hash_table, &tokens[5][0]);
		new_node->model_name = -1;
		sscanf(tokens[6], "%Lf", &new_node->length);
		sscanf(tokens[7], "%Lf", &new_node->width);
	}
	else if (new_node->type == 'Q' || new_node->type == 'q') {
		new_node->probe1 = (char *)malloc(strlen(&tokens[2][0]) * sizeof(char));
		new_node->probe2 = (char *)malloc(strlen(&tokens[3][0]) * sizeof(char));
		new_node->probe3 = (char *)malloc(strlen(&tokens[4][0]) * sizeof(char));
		assert(new_node->probe1 != NULL);
		assert(new_node->probe2 != NULL);
		assert(new_node->probe3 != NULL);
		strcpy(new_node->probe1, &tokens[2][0]);
		strcpy(new_node->probe1, &tokens[3][0]);
		strcpy(new_node->probe1, &tokens[4][0]);
		/* Add probes to the hashtable */
		ht_set(hash_table, &tokens[2][0]);
		ht_set(hash_table, &tokens[3][0]);
		ht_set(hash_table, &tokens[4][0]);
		new_node->model_name = -1;
		new_node->length = new_node->width = -1;
		/* For now we ignore mode_name, thus we use 4 */
		new_node->area = num_tokens > 4 ? atoi(&tokens[5][0]) : 1;
	}
	return SUCCESS;
}

/* Free all the dynamic memory allocated for the lists index */
void free_index(index_t **index) {
	/* Free the the lists */
	free_list1(&((*index)->head1), &((*index)->tail1));
	free_list2(&((*index)->head2), &((*index)->tail2));
	free(*index);
	/* Set index to NULL to limit further acesses */
	*index = NULL;
}

/* Free the elements from list1 */
void free_list1(list1_t **head, list1_t **tail) {
	list1_t *curr = *head;
	list1_t *next;
	while (curr != NULL) {
		/* Keep the next node */
		next = curr->next;
		/* Free everything from the current node */
		free(curr->element);
		free(curr->probe1);
		free(curr->probe2);
		free(curr);
		curr = next;
	}
	/* Set head/tail to NULL to limit further acesses */
	*head = *tail = NULL;
}

/* Free the elements from list2 */
void free_list2(list2_t **head, list2_t **tail) {
	list2_t *curr = *head;
	list2_t *next;
	while (curr != NULL) {
		/* Keep the next node*/
		next = curr->next;
		/* Free everything from the current node */
		free(curr->element);
		free(curr->probe1);
		free(curr->probe2);
		free(curr->probe3);
		free(curr->probe4);
		free(curr);
		curr = next;
	}
	/* Set head/tail to NULL to limit further acesses */
	*head = *tail = NULL;
}

/* Print the lists */
void print_lists(index_t *index, hash_table_t *hash_table) {
	print_list1(index->head1, hash_table);
	print_list2(index->head2, hash_table);
}

/* Print list1 elements */
void print_list1(list1_t *head, hash_table_t *hash_table) {
	list1_t *curr = head;
	while (curr != NULL) {
		printf("Type: %c\n", curr->type);
		printf("Element: %s\n", curr->element);
		printf("Probe1: %s\n", curr->probe1);
		printf("Probe2: %s\n", curr->probe2);
		printf("Value: %.20Lf\n", curr->value);
		printf("Probe1 id: %d\n", ht_get_id(hash_table, curr->probe1));
		printf("Probe2 id: %d\n", ht_get_id(hash_table, curr->probe2));
		printf("\n");
		curr = curr->next;
	}
}

/* Print list2 elements */
void print_list2(list2_t *head, hash_table_t *hash_table) {
	list2_t *curr = head;
	while (curr != NULL) {
		printf("Type: %c\n", curr->type);
		printf("Element: %s\n", curr->element);
		printf("Probe1: %s\n", curr->probe1);
		printf("Probe2: %s\n", curr->probe2);
		printf("Probe1 id: %d\n", ht_get_id(hash_table, curr->probe1));
		printf("Probe2 id: %d\n", ht_get_id(hash_table, curr->probe2));
		if (curr->probe3 != NULL) {
			printf("Probe3: %s\n", curr->probe3);
			printf("Probe3 id: %d\n", ht_get_id(hash_table, curr->probe3));
		}
		else {
			printf("Probe3: NULL\n");
		}
		if (curr->probe4 != NULL) {
			printf("Probe4: %s\n", curr->probe4);
			printf("Probe4 id: %d\n", ht_get_id(hash_table, curr->probe4));
		}
		else {
			printf("Probe4: NULL\n");
		}
		printf("Model_name: %d\n", curr->model_name);
		printf("Area: %d\n", curr->area);
		printf("Length %.20Lf\nWidth: %.20Lf\n", curr->length, curr->width);
		printf("\n");
		curr = curr->next;
	}
}