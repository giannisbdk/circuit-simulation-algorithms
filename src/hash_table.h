#ifndef HASH_TABLE_H
#define HASH_TABLE_H

typedef struct entry {
	char *key;
	unsigned int id;
	struct entry *next;
} entry_t;


typedef struct hash_table {
	unsigned int size;
	unsigned int seq;
	entry_t **table;
} hash_table_t;

hash_table_t *ht_create(int size);
entry_t *ht_new_node(hash_table_t *hash_table, char *key);
int ht_hash(hash_table_t *hash_table, char *key );
void ht_set(hash_table_t *hash_table, char *key);
int ht_get_id(hash_table_t *hash_table, char *key);
void ht_free(hash_table_t **hash_table);

#endif