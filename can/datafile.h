#ifndef DATAFILE_H
#define DATAFILE_H

#ifdef __cplusplus
extern "C" {
#endif

struct var {
    char *name;
    double value;
    struct var *next;
} variable;

typedef struct var varible;

char** read_datafile(char*);
int delete_buffer(char**);
struct var* read_line(char*);
char* remove_comments(char*);
int store_data(struct var*);
int get_vars(char*);
struct var* new_var();
void destroy_var(struct var*);
int print_global_vars(void);
int get_vars(char*);

#ifdef __cplusplus
}
#endif

#endif
