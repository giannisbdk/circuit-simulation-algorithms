#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <assert.h>

#include "hash_table.h"

/* Create a new hash_table. */
hash_table_t *ht_create(int size) {
	hash_table_t *hash_table = NULL;
	int i;
	if(size < 1) {
		return NULL;
	}

	/* Allocate the table itself. */
	if((hash_table = (hash_table_t *)malloc(sizeof(hash_table_t))) == NULL) {
		return NULL;
	}

	/* Allocate pointers to the head nodes. */
	if((hash_table->table = (entry_t **)malloc(sizeof(entry_t *) * size)) == NULL) {
		free(hash_table);
		return NULL;
	}
	for(i = 0; i < size; i++) {
		hash_table->table[i] = NULL;
	}
	hash_table->size = size;
	hash_table->seq = 1;
	return hash_table;	
}

/* Hash a string for a particular hash table. */
int ht_hash( hash_table_t *hash_table, char *key ) {
	unsigned long int hashval = 0;
	int i = 0;
	/* Convert our string to an integer */
	while(hashval < ULONG_MAX && i < strlen(key)) {
		hashval = hashval << 8;
		hashval += key[ i ];
		i++;
	}
	return hashval % hash_table->size;
}

// unsigned long ht_hash(char* str) {

// 	unsigned long hash = 5381;
// 	int c;

// 	while ((c = *str++))
// 		hash = ((hash << 5) + hash) + c;  /* hash * 33 + c */

// 	return hash;
// }

/* Create a key-value pair. */
entry_t *ht_new_node(hash_table_t *hash_table, char *key) {

	entry_t *new_node;

	if((new_node = (entry_t *)malloc(sizeof(entry_t))) == NULL) {
		return NULL;
	}

	if((new_node->key = strdup(key)) == NULL) {
		return NULL;
	}
	/* Give to ground the id 0 */
	if(strcmp(key, "0") == 0) {
		new_node->id = 0;
	}
	else {
		new_node->id = hash_table->seq;
		hash_table->seq++;
	}
	new_node->next = NULL;
	return new_node;
}

/* Insert a key-value pair into a hash table. */
void ht_set(hash_table_t *hash_table, char *key) {
	int bin = 0;
	bin = ht_hash(hash_table, key);
	entry_t *new_node = NULL;
	entry_t *curr = NULL;
	entry_t *last = NULL;

	curr = hash_table->table[bin];
	/* In case its empty create new node */
	if(curr == NULL) {
		new_node = ht_new_node(hash_table, key);
#ifdef DEBUGH
		printf("Inserted new node with key: %s and id: %d at bin: %d\n", key, new_node->id, bin);
#endif
		hash_table->table[bin] = new_node;
		return;
	}
	else {
		for(; curr != NULL; curr = curr->next) {
			if(strcmp(key, curr->key) == 0) {
				/* Means we have already stored that string */
#ifdef DEBUGH
				printf("Key: %s, has already been stored.\n", key);
#endif
				return;
			}
			last = curr;
		}
	}
	/* Last contains the last node in the list */
	new_node = ht_new_node(hash_table, key);
	if(new_node == NULL) {
		printf("Something went wrong\n");
		return;
	}
	last->next = new_node;
#ifdef DEBUGH
	printf("Inserted new node at the end of list at bin: %d with key: %s and id: %d\n", bin, key, new_node->id);
#endif
}

/* Retrieve a key-value pair from a hash table. */
int ht_get_id(hash_table_t *hash_table, char *key) {

	int bin = 0;
	entry_t *curr = NULL;
	bin = ht_hash(hash_table, key);
	
	/* Step through the bin, looking for our value. */
	curr = hash_table->table[bin];
	while(curr != NULL && strcmp(key, curr->key) != 0) {
		curr = curr->next;
	}
	if(curr == NULL) {
		return -1;
	}
	else {
		return curr->id;
	}
}

/* Free the memory allocated */
void ht_free(hash_table_t **hash_table) {
	entry_t *curr, *prev;
	/* Cycle through the hash_table */
	for (int i = 0; i < (*hash_table)->size; i++) {
		curr = (*hash_table)->table[i];
		/* Free every node in the list in a cell of the hash_table */
		while(curr != NULL) {
			prev = curr;
			curr = curr->next;
			free(prev);
		}
	}
	free(*hash_table);
	*hash_table = NULL;
}