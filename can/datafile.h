#ifndef DATAFILE_H
#define DATAFILE_H

typedef struct {
    char *name;
    double value;
} variable;

char** read_datafile(char*);
int delete_buffer(char**);
variable read_line(char*);
char* remove_comments(char*);
int store_data(variable);
int get_vars(char*);

#endif
