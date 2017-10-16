//This is a small implementation of a extensible hashtable
//Initialy the user must specify the capacity of the hashtable in order to perform any additions to it
//If the size becomes bigger than the upper_bound_ratio then the hashtable capacity becomes twice as big

//Functions to use in main are:
//	-hash_table_t *ht_create(int capacity) *Creation of hashtable
//	-int ht_add(hash_table_t *hashtable, char *str) *This function enters the input string into the hashtable and returns the a unique key whitch is a positive number. If the element already exists in the hashtable the function returns the key of it.
// 	-char *get_ht_value(hash_table_t *hashtable,int hashcode) *This function simply returns the inner of the specified hashtable cell


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hash_table.h"

hash_table_t *ht_create(int capacity) {
	hash_table_t *hashtable;

	if((hashtable = (hash_table_t *)malloc(sizeof(hash_table_t))) == NULL) {
		return NULL;
	}
	if((hashtable->table = (char **)malloc(capacity * sizeof(char *))) == NULL) {
		return NULL;
	}
	for(int i = 0; i < capacity; i++) {
		hashtable->table[i] = NULL;
	}
	hashtable->size = 0;
	hashtable->capacity = capacity;
	hashtable->upper_bound_ratio=0.9;
	return hashtable;
}

int get_hashcode(hash_table_t *hashtable, char *str) {
	int sum=0;
	for(int i=0; i<strlen(str); i++) {
		sum += str[i];
	}
	return sum % hashtable->capacity;
}

int ht_contains(hash_table_t *hashtable, char *str) {
	for (int n = 0; n < hashtable->capacity-1; n++){
		if(hashtable->table[n] != NULL){
			if(strcmp(hashtable->table[n], str) == 0) {
				return n;
			}
		}
	}	
	return FAILURE;
}

int ht_inner_add(hash_table_t *hashtable, char *str) {
	int hashcode = get_hashcode(hashtable, str);
	int contains = ht_contains(hashtable, str);
	if (contains != FAILURE) {
		return contains;
	} else{
		for (int n = 0; n < hashtable->capacity-1; n++){
			if(hashtable->table[hashcode] == NULL) {
				hashtable->table[hashcode] = (char *)malloc(sizeof(char) * strlen(str));
				strcpy(hashtable->table[hashcode], str);
				hashtable->size++;

				return hashcode;
			}
			else if(n == hashtable->capacity-1) {
				return FAILURE;
			}
			if(hashcode == hashtable->capacity-1) {
				hashcode = 0;
			}
			else {
				hashcode++;
			}
		}
	}
	return FAILURE;
}

int extend_ht(hash_table_t *hashtable){
	if(hashtable->size/(double)(hashtable->capacity) > hashtable->upper_bound_ratio) {
		hashtable->capacity = hashtable->capacity*2;
		hashtable->table = (char **)realloc(hashtable->table,hashtable->capacity * sizeof(char *));
		if(hashtable->table == NULL) {
			return FAILURE;
		}else{
			for(int i = hashtable->capacity/2; i < hashtable->capacity; i++) {
				hashtable->table[i] = NULL;
			}
		}
		return SUCCESS;

	}else{
		return 0;
	}
}

int ht_add(hash_table_t *hashtable, char *str){
	int rehash = 0;
	int add = ht_inner_add(hashtable,str);
	if (add != FAILURE){
		rehash = extend_ht(hashtable);
		if(rehash == FAILURE) {
			return FAILURE;
		}
	}
	return add;
}

char *ht_get_value(hash_table_t *hashtable, int hashcode) {
	// Check if value is NULL or NOT?

	if(hashcode == -1) {
		return "empty";
	}

	if(hashtable->table[hashcode] != NULL){
		return hashtable->table[hashcode];
	} else {
		printf("The hashcode provided leads to a not existing cell\n");
		return NULL;
	}
}
