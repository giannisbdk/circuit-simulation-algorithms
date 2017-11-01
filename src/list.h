#ifndef LIST_H
#define LIST_H

#include "hash_table.h"

#define SUCCESS	 1
#define FAILURE	-1

typedef struct list1 {
	char type;
	char *element;
	char *probe1;
	char *probe2;
	long double value;
	struct list1 *next;
	struct list1 *prev;
} list1_t;

typedef struct list2 {
	char type;
	char *element;
	char *probe1;
	char *probe2;
	char *probe3;
	char *probe4;
	int model_name;
	int area;	
	long double length;
	long double width;
	struct list2 *next;
	struct list2 *prev;
} list2_t;

typedef struct index {
	/* Pointers to the head and tail of the list1 */
	list1_t *head1;
	list1_t *tail1;
	/* Pointers to the head and tail of the list2 */
	list2_t *head2;
	list2_t *tail2;
	unsigned int size1;
	unsigned int size2;
} index_t;


index_t *init_lists();
int add_to_list(index_t *index, char **tokens, hash_table_t *hash_table);
int add_to_list1(index_t *index, char **tokens, hash_table_t *hash_table);
int add_to_list2(index_t *index, char **tokens, hash_table_t *hash_table);
void free_index(index_t **index);
void free_list1(list1_t **head, list1_t **tail);
void free_list2(list2_t **head, list2_t **tail);
void print_lists(index_t *index, hash_table_t *hash_table);
void print_list1(list1_t *head, hash_table_t *hash_table);
void print_list2(list2_t *head, hash_table_t *hash_table);

#endif