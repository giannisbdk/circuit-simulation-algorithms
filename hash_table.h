#ifndef HASHTABLE_H
#define HASHTABLE_H

#define SUCCESS 1
#define FAILURE -1

typedef struct hash_table_s {
	unsigned int size;
  	unsigned int capacity;
  	double upper_bound_ratio;
  	char **table;
} hash_table_t;

hash_table_t *ht_create(int);
int ht_add(hash_table_t *, char *);
char *ht_get_value(hash_table_t *, int);

#endif