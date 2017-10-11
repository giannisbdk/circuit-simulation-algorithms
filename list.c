#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "list.h"

list_uno_t *init_list1() {

	list_uno_t *head = (list_uno_t *)malloc(sizeof(list_uno_t));
	head->type = '\0';
	head->next = NULL;
	head->prev = NULL;
	return head;
}

list_dos_t *init_list2() {

	list_dos_t *head = (list_dos_t *)malloc(sizeof(list_dos_t));
	head->type = '\0';
	head->area = 1;
	/* Ignore for now the model_name */
	head->model_name = -1;
	head->prev = NULL;
	head->next = NULL;
	return head;
}

/* General method to add new node to the appropriate list */
int add_to_list(list_uno_t *head_uno, list_dos_t *head_dos, char **tokens) {

	char type = tokens[1][0];

	if(type == 'R' || type == 'C' || type == 'I' || type == 'V' || type == 'L') {
		add_to_list1(head_uno, tokens);
	}
	else if(type == 'M' || type == 'Q' || type == 'D') {
		add_to_list2(head_dos, tokens);
	}
	else {
		printf("Something is totally wrong with the input.\n");
		return -1;
	}
	return 1;
}

void add_to_list1(list_uno_t *head, char **tokens) {

	//TODO perhaps add error handling with token number
	// int num_tokens = atoi(tokens[0]);
	list_uno_t *new_node;

	/* In case the list is empty */
	if(head->type == '\0') {
		new_node = head;
	}
	else {
		new_node = (list_uno_t *)malloc(sizeof(list_uno_t));
		list_uno_t *curr;
		/* Find the last node */
		for(curr = head; curr->next != NULL; curr = curr->next);
		curr->next = new_node;
		new_node->prev = curr;
		new_node->next = NULL;
	}
	/* Set new node struct fields */
	strncpy(&new_node->type, &tokens[1][0], 1);
	new_node->id = atoi(&tokens[1][1]);
	new_node->probe1 = atoi(&tokens[2][0]);
	new_node->probe2 = atoi(&tokens[3][0]);
	sscanf(tokens[4], "%Lf", &new_node->value);
}

void add_to_list2(list_dos_t *head, char **tokens) {

	//TODO error handling
	int num_tokens = atoi(tokens[0]);
	list_dos_t *new_node;

	/* In case the list is empty */
	if(head->type == '\0') {
		new_node = head;
	}
	else {
		new_node = (list_dos_t *)malloc(sizeof(list_dos_t));
		list_dos_t *curr;
		/* Find the last node */
		for(curr = head; curr->next != NULL; curr = curr->next);
		curr->next = new_node;
		new_node->prev = curr;
		new_node->next = NULL;
	}

	/* Set new node struct fields */
	strncpy(&new_node->type, &tokens[1][0], 1);
	new_node->id = atoi(&tokens[1][1]);
	switch(new_node->type) {
		case 'D':
			new_node->probe1 = atoi(&tokens[2][0]);
			new_node->probe2 = atoi(&tokens[3][0]);
			new_node->probe3 = -1;
			new_node->probe4 = -1;
			new_node->model_name = -1;
			new_node->length = new_node->width = -1;
			/* For now we ignore mode_name, thus we use 3 */
			new_node->area = num_tokens > 3 ? atoi(&tokens[4][0]) : 1;
			break;
		case 'M':
			new_node->probe1 = atoi(&tokens[2][0]);
			new_node->probe2 = atoi(&tokens[3][0]);
			new_node->probe3 = atoi(&tokens[4][0]);
			new_node->probe4 = atoi(&tokens[5][0]);
			new_node->model_name = -1;
			sscanf(tokens[6], "%Lf", &new_node->length);
			sscanf(tokens[7], "%Lf", &new_node->width);
			break;
		case 'Q':
			new_node->probe1 = atoi(&tokens[2][0]);
			new_node->probe2 = atoi(&tokens[3][0]);
			new_node->probe3 = atoi(&tokens[4][0]);
			new_node->probe4 = -1;
			new_node->model_name = -1;
			new_node->length = new_node->width = -1;
			/* For now we ignore mode_name, thus we use 4 */
			new_node->area = num_tokens > 4 ? atoi(&tokens[5][0]) : 1;
			break;
		default:
			printf("Unknown circuit element\n");
	}
}

void remove_from_list1(list_uno_t *head, char **tokens);
void remove_from_list2(list_dos_t *head, char **tokens);

void print_lists(list_uno_t *head_uno, list_dos_t *head_dos) {

	print_list1(head_uno);
	print_list2(head_dos);
}

void print_list1(list_uno_t *head) {

	list_uno_t *curr = head;

	while(curr != NULL) {
		printf("Type: %c\n", curr->type);
		printf("Id: %d\n", curr->id);
		printf("Probe1: %d\n", curr->probe1);
		printf("Probe2: %d\n", curr->probe2);
		printf("Value: %.20Lf\n", curr->value);
		printf("\n");
		curr = curr->next;
	}
}

void print_list2(list_dos_t *head) {

	list_dos_t *curr = head;

	while(curr != NULL) {
		printf("Type: %c\n", curr->type);
		printf("Id: %d\n", curr->id);
		printf("Probe1: %d\n", curr->probe1);
		printf("Probe2: %d\n", curr->probe2);
		printf("Probe3: %d\n", curr->probe3);
		printf("Probe4: %d\n", curr->probe4);
		printf("Model_name: %d\n", curr->model_name);
		printf("Area: %d\n", curr->area);
		printf("Length %.20Lf\nWidth: %.20Lf\n", curr->length, curr->width);
		printf("\n");
		curr = curr->next;
	}
}
