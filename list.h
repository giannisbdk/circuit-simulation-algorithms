#ifndef LIST_H
#define LIST_H

#include "hash_table.h"

#define SUCCESS 	 1
#define FAILURE -1

typedef struct list1 {
	char type;
	char *id;
	char *probe1;
	char *probe2;
	// int probe1;
	// int probe2;
	long double value;
	struct list1 *next;
	struct list1 *prev;
} list1_t;

typedef struct list2 {
	char type;
	char *id;
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
	list1_t *head1;
	list2_t *head2;
	unsigned int size1;
	unsigned int size2;
} index_t;


index_t *init_lists();
int add_to_list(index_t *, char **, hash_table_t *);
list1_t *add_to_list1(list1_t *, int, char **, hash_table_t *);
list2_t *add_to_list2(list2_t *, int, char **, hash_table_t *);
void print_lists(index_t *, hash_table_t *);
void print_list1(list1_t *, hash_table_t *);
void print_list2(list2_t *, hash_table_t *);

#endif