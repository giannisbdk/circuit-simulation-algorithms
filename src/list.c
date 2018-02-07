#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "list.h"

/* Holds the value of the non-zero elements */
static int nz = 0;

/* Initializes the lists */
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

/* Adds an element to the lists */
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
		fprintf(stderr, "Circuit element has wrong type %c\n", type);
		return FAILURE;
	}
	return SUCCESS;
}

/* Computes the number of contributions (e.g nonzeros) to the MNA for the current element */
void set_nz(char type, char probe1, char probe2) {
	int mul = 1;
	if (probe1 == '0' || probe2 == '0') {
		mul = 0;		
	}
	if (type == 'r' || type == 'R') {
		nz += 3*mul + 1;
	}
	else if (type == 'v' || type == 'V' || type == 'l' || type == 'L') {
		nz += 2*mul + 2;
	}
}

/* Simple getter for the nonzeros */
int get_nz() {
	return nz;
}

/* Add new node to the end of the list */
int add_to_list1(index_t *index, char **tokens, hash_table_t *hash_table) {
	list1_t *new_node = (list1_t *)malloc(sizeof(list1_t));
	if (new_node == NULL) {
		fprintf(stderr, "Could not allocate memory for new node!\n");
		return FAILURE;
	}
	/* In case list is empty */
	if (index->size1 == 0) {
		new_node->prev = NULL;
		new_node->next = NULL;
		index->head1   = new_node;
		index->tail1   = new_node;
#ifdef DEBUGL
		printf("Adding node to the head of list1, new_node: %p, head1: %p, tail1: %p\n", new_node, index->head1, index->tail1);
#endif
	}
	/* If it is not empty add element to the end of the list */
	else {
		new_node->prev     = index->tail1;
		index->tail1->next = new_node;
		index->tail1       = new_node;
		new_node->next     = NULL;
#ifdef DEBUGL
		printf("Adding node to the end of list1, new_node: %p, new_node->prev: %p\n", new_node, new_node->prev);
#endif
	}
	/* Set new node struct fields */
	strncpy(&new_node->type, &tokens[1][0], 1);
	new_node->element = (char *)malloc((strlen(&tokens[1][0]) + 1) * sizeof(char));
	new_node->probe1  = (char *)malloc((strlen(&tokens[2][0]) + 1) * sizeof(char));
	new_node->probe2  = (char *)malloc((strlen(&tokens[3][0]) + 1) * sizeof(char));
	assert(new_node->element != NULL);
	assert(new_node->probe1  != NULL);
	assert(new_node->probe2  != NULL);
	strcpy(new_node->element, &tokens[1][0]);
	strcpy(new_node->probe1,  &tokens[2][0]);
	strcpy(new_node->probe2,  &tokens[3][0]);

	set_nz(new_node->type, new_node->probe1[0], new_node->probe2[0]);
	/* Add probes to the hashtable */
	// TODO perhaps replace &tokens with new_node probes
	ht_set(hash_table, &tokens[2][0]);
	ht_set(hash_table, &tokens[3][0]);
	sscanf(tokens[4], "%Lf", &new_node->value);

	/* Initialize to null both */
	new_node->trans_spec = NULL;
	new_node->ac_spec    = NULL;

	/* Get the total number of tokens */
	int num_tokens = atoi(tokens[0]);

	//TODO create function to set the transient and ac specs
	/* Check if transient or ac spec exists */
	if (num_tokens > 4) {
		trans_type type;
		int ac_index = 5;
		bool AC = check_ac(tokens, num_tokens);
		if (is_transient(tokens[5], &type)) {
			new_node->trans_spec = (trans_spec_t *)malloc(sizeof(trans_spec_t));
			new_node->trans_spec->exp 	= NULL;
			new_node->trans_spec->sin 	= NULL;
			new_node->trans_spec->pulse = NULL;
			new_node->trans_spec->pwl 	= NULL;
			/* Set the transient spec according to the type */
			switch (type) {
				case EXP:
					new_node->trans_spec->exp  = (exp_t *)malloc(sizeof(exp_t));
					new_node->trans_spec->type = type;
					/* Strip parentheses of first token and last token */
					sscanf(tokens[6],  "(%lf", &(new_node->trans_spec->exp->i1));
					sscanf(tokens[7],  "%lf",  &(new_node->trans_spec->exp->i2));
					sscanf(tokens[8],  "%lf",  &(new_node->trans_spec->exp->td1));
					sscanf(tokens[9],  "%lf",  &(new_node->trans_spec->exp->tc1));
					sscanf(tokens[10], "%lf",  &(new_node->trans_spec->exp->td2));
					sscanf(tokens[11], "%lf)", &(new_node->trans_spec->exp->tc2));
					ac_index = 12;
					break;
				case SIN:
					new_node->trans_spec->sin  = (sin_t *)malloc(sizeof(sin_t));
					new_node->trans_spec->type = type;
					/* Strip parentheses of first token and last token */
					sscanf(tokens[6],  "(%lf", &(new_node->trans_spec->sin->i1));
					sscanf(tokens[7],  "%lf",  &(new_node->trans_spec->sin->ia));
					sscanf(tokens[8],  "%lf",  &(new_node->trans_spec->sin->fr));
					sscanf(tokens[9],  "%lf",  &(new_node->trans_spec->sin->td));
					sscanf(tokens[10], "%lf",  &(new_node->trans_spec->sin->df));
					sscanf(tokens[11], "%lf)", &(new_node->trans_spec->sin->ph));
					ac_index = 12;
					break;
				case PULSE:
					new_node->trans_spec->pulse = (pulse_t *)malloc(sizeof(pulse_t));
					new_node->trans_spec->type  = type;
					/* Strip parentheses of first token and last token */
					sscanf(tokens[6],  "(%lf", &(new_node->trans_spec->pulse->i1));
					sscanf(tokens[7],  "%lf",  &(new_node->trans_spec->pulse->i2));
					sscanf(tokens[8],  "%lf",  &(new_node->trans_spec->pulse->td));
					sscanf(tokens[9],  "%lf",  &(new_node->trans_spec->pulse->tr));
					sscanf(tokens[10], "%lf",  &(new_node->trans_spec->pulse->tf));
					sscanf(tokens[11], "%lf",  &(new_node->trans_spec->pulse->pw));
					sscanf(tokens[12], "%lf)", &(new_node->trans_spec->pulse->per));
					ac_index = 13;
					break;
				case PWL:
					new_node->trans_spec->pwl  = (pwl_t *)malloc(sizeof(pwl_t));
					new_node->trans_spec->type = type;
					/* Start shows the index of trans_spec tokens */
					int start = 6;
					/* Get the number of states in the pwl, this differ whether AC exists or not, because of num_tokens */
					if (AC) {
						new_node->trans_spec->pwl->n = ((num_tokens - 3) - (start -1)) / 2;
					}
					else {
						new_node->trans_spec->pwl->n = (num_tokens - (start - 1)) / 2;
					}
					new_node->trans_spec->pwl->t = (double *)malloc(new_node->trans_spec->pwl->n * sizeof(double));
					new_node->trans_spec->pwl->i = (double *)malloc(new_node->trans_spec->pwl->n * sizeof(double));
					for (int i = 0; i < new_node->trans_spec->pwl->n; i++) {
						sscanf(tokens[start + (i * 2)],     "(%lf", &(new_node->trans_spec->pwl->t[i]));
						sscanf(tokens[start + (i * 2) + 1], "%lf)", &(new_node->trans_spec->pwl->i[i]));
						ac_index = start + (i * 2) + 2;
					}
					break;
				default:
					fprintf(stderr, "Wrong transient spec type: %s\n", tokens[5]);
					exit(EXIT_FAILURE);
			}
		}
		/* Set AC struct in case it exists in the tokens */
		if (AC) {
			new_node->ac_spec = (ac_spec_t *)malloc(sizeof(ac_spec_t));
			sscanf(tokens[ac_index + 1], "%lf",  &(new_node->ac_spec->magnitude));
			sscanf(tokens[ac_index + 2], "%lf",  &(new_node->ac_spec->phase));
		}
	}
	return SUCCESS;
}

/* Add new node to the end of the list */
int add_to_list2(index_t *index, char **tokens, hash_table_t *hash_table) {
	int num_tokens = atoi(tokens[0]);
	list2_t *new_node = (list2_t *)malloc(sizeof(list2_t));
	if (new_node == NULL) {
		fprintf(stderr, "Could not allocate memory for new node!\n");
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

/* Check if AC exists inside tokens */
bool check_ac(char **tokens, int num_tokens) {
	if (strcasecmp(tokens[num_tokens - 2], "AC") == 0) {
		return true;
	}
	return false;
}

/* Check if it is a transient or ac and in case it is transient it sets type to the appropriate enum */
bool is_transient(char *spec, trans_type *type) {
	if (strcasecmp(spec, "EXP") == 0) {
		*type = EXP;
		return true;
	}
	else if (strcasecmp(spec, "SIN") == 0) {
		*type = SIN;
		return true;
	}
	else if (strcasecmp(spec, "PULSE") == 0) {
		*type = PULSE;
		return true;
	}
	else if (strcasecmp(spec, "PWL") == 0) {
		*type = PWL;
		return true;
	}
	else {
		/* It means its AC */
		type = NULL;
		return false;
	}
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
		printf("\nType:      %c\n", curr->type);
		printf("Element:   %s\n",   curr->element);
		printf("Probe1:    %s\n",   curr->probe1);
		printf("Probe2:    %s\n",   curr->probe2);
		printf("Value:     %Lf\n",  curr->value);
		printf("Probe1 id: %d\n",   ht_get_id(hash_table, curr->probe1));
		printf("Probe2 id: %d\n",   ht_get_id(hash_table, curr->probe2));
		if (curr->trans_spec != NULL) {
			print_trans_spec(curr->trans_spec);
		}
		if (curr->ac_spec != NULL) {
			print_ac_spec(curr->ac_spec);
		}
		printf("\n");
		curr = curr->next;
	}
}

/* Prints the given transient spec */
void print_trans_spec(trans_spec_t *trans_spec) {
	switch (trans_spec->type) {
		case EXP:
			printf("Transient Spec: EXP\n");
			printf("i1:  %lf\n", trans_spec->exp->i1);
			printf("i2:  %lf\n", trans_spec->exp->i2);
			printf("td1: %lf\n", trans_spec->exp->td1);
			printf("tc1: %lf\n", trans_spec->exp->tc1);
			printf("td2: %lf\n", trans_spec->exp->td2);
			printf("tc2: %lf\n", trans_spec->exp->tc2);
			break;
		case SIN:
			printf("Transient Spec: SIN\n");
			printf("i1: %lf\n", trans_spec->sin->i1);
			printf("ia: %lf\n", trans_spec->sin->ia);
			printf("fr: %lf\n", trans_spec->sin->fr);
			printf("td: %lf\n", trans_spec->sin->td);
			printf("df: %lf\n", trans_spec->sin->df);
			printf("ph: %lf\n", trans_spec->sin->ph);
			break;
		case PULSE:
			printf("Transient Spec: PULSE\n");
			printf("i1:  %lf\n", trans_spec->pulse->i1);
			printf("i2:  %lf\n", trans_spec->pulse->i2);
			printf("td:  %lf\n", trans_spec->pulse->td);
			printf("tr:  %lf\n", trans_spec->pulse->tr);
			printf("tf:  %lf\n", trans_spec->pulse->tf);
			printf("pw:  %lf\n", trans_spec->pulse->pw);
			printf("per: %lf\n", trans_spec->pulse->per);
			break;
		case PWL:
			printf("Transient Spec: PWL\n");
			for (int i = 0; i < trans_spec->pwl->n; i++) {
				printf("t[%d]: %6g\ti[%d]: %6g\n", i, trans_spec->pwl->t[i], i, trans_spec->pwl->i[i]);
			}
			break;
		default:
			fprintf(stderr, "Wrong transient type %d\n", trans_spec->type);
	}
}

/* Prints the given ac_spec */
void print_ac_spec(ac_spec_t *ac_spec) {
	printf("AC Spec:\n");
	printf("Magnitude: %lf\n", ac_spec->magnitude);
	printf("Phase:     %lf\n", ac_spec->phase);
}

/* Print list2 elements */
void print_list2(list2_t *head, hash_table_t *hash_table) {
	list2_t *curr = head;
	while (curr != NULL) {
		printf("\nType:      %c\n", curr->type);
		printf("Element:   %s\n",   curr->element);
		printf("Probe1:    %s\n",   curr->probe1);
		printf("Probe2:    %s\n",   curr->probe2);
		printf("Probe1 id: %d\n",   ht_get_id(hash_table, curr->probe1));
		printf("Probe2 id: %d\n",   ht_get_id(hash_table, curr->probe2));
		if (curr->probe3 != NULL) {
			printf("Probe3:    %s\n", curr->probe3);
			printf("Probe3 id: %d\n", ht_get_id(hash_table, curr->probe3));
		}
		else {
			printf("Probe3:    NULL\n");
		}
		if (curr->probe4 != NULL) {
			printf("Probe4:    %s\n", curr->probe4);
			printf("Probe4 id: %d\n", ht_get_id(hash_table, curr->probe4));
		}
		else {
			printf("Probe4: NULL\n");
		}
		printf("Model_name: %d\n",   curr->model_name);
		printf("Area:       %d\n",   curr->area);
		printf("Length:     %.Lf\n", curr->length);
		printf("Width:      %.Lf\n", curr->width);
		printf("\n");
		//TODO Check if TRAN/AC spec is required to be printed
		curr = curr->next;
	}
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

		/* Remove the transient spec in case it exists */
		if (curr->trans_spec != NULL) {
			switch (curr->trans_spec->type) {
				case EXP:
					free(curr->trans_spec->exp);
					break;
				case SIN:
					free(curr->trans_spec->sin);
					break;
				case PULSE:
					free(curr->trans_spec->pulse);
					break;
				case PWL:
					free(curr->trans_spec->pwl->t);
					free(curr->trans_spec->pwl->i);
					free(curr->trans_spec->pwl);
					break;
			}
			free(curr->trans_spec);
		}
		/* Remove AC struct in case it exists */
		if (curr->ac_spec != NULL) {
			free(curr->ac_spec);
		}
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