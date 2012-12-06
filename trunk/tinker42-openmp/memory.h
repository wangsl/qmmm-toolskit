 
/* $Id: memory.h 2 2008-02-28 22:52:52Z wangsl $ */

#ifndef MEMORY_H
#define MEMORY_H

#define alloc_double(x) allocate_double_loc(x, __FILE__, __LINE__)
#define alloc_int(x) allocate_int_loc(x, __FILE__, __LINE__)
#define alloc_char(x) allocate_char_loc(x, __FILE__, __LINE__)

#define cfree(x) free_loc(x, __FILE__, __LINE__)

#define memory_information() show_memory_allocation_information(__FILE__, __LINE__)

void *allocate_void_loc(const int &n, const char *file_name, const int &line);
double *allocate_double_loc(const int &n, const char *file_name, const int &line);
int *allocate_int_loc(const int &n, const char *file_name, const int &line);
char *allocate_char_loc(const int &n, const char *file_name, const int &line);

void *allocate_void_loc(const long int &n, const char *file_name, const int &line);
double *allocate_double_loc(const long int &n, const char *file_name, const int &line);
int *allocate_int_loc(const long int &n, const char *file_name, const int &line);
char *allocate_char_loc(const long int &n, const char *file_name, const int &line);

void free_loc(void * &p, const char *file_name, const int &line);
void free_loc(double * &p, const char *file_name, const int &line);
void free_loc(int * &p, const char *file_name, const int &line);
void free_loc(char * &p, const char *file_name, const int &line);

void show_memory_allocation_information(const char *file_name, const int &line);

#endif /* MEMORY_H */
