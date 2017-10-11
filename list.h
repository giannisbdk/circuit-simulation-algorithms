#ifndef LIST_H
#define LIST_H

typedef struct list_uno {
	char type;
	int id;
	int probe1;
	int probe2;
	long double value;
	struct list_uno *next;
	struct list_uno *prev;
} list_uno_t;

typedef struct list_dos {
	char type;
	int id;
	int probe1;
	int probe2;
	int probe3;
	int probe4;
	int model_name;
	int area;	
	long double length;
	long double width;
	struct list_dos *next;
	struct list_dos *prev;
} list_dos_t;

list_uno_t *init_list1();
list_dos_t *init_list2();

int add_to_list(list_uno_t *head_uno, list_dos_t *head_dos, char **tokens);
void add_to_list1(list_uno_t *head, char **tokens);
void add_to_list2(list_dos_t *head, char **tokens);
void remove_from_list1(list_uno_t *head, char **tokens);
void remove_from_list2(list_dos_t *head, char **tokens);
void print_lists(list_uno_t *head_uno, list_dos_t *head_dos);
void print_list1(list_uno_t *head);
void print_list2(list_dos_t *head);

#endif