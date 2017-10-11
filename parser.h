#ifndef PARSER_H
#define PARSER_H

int get_num_tokens(char *line);
char **tokenizer(char *line, int *num_tokens);

#endif