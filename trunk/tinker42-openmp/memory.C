 
/* $Id$ */

#include <iostream>
using namespace std;
#include <malloc.h>

#include "fort.h"
#include "fint.h"
#include "memptr.h"
#include "die.h"

#define MAXMEMNUMBER 1000

static MemoryPointer mp[MAXMEMNUMBER];
static long int memory_in_use = 0;

// C++: memory_information()
void show_memory_allocation_information(const char *file, const int &line)
{
  cout << "\n"
       << " Memory allocation information file: " <<  file << " line: " << line << "\n"
       << " Total memory in use: " << memory_in_use << " byte, " 
       << memory_in_use/(1024.0*1024.0) << " megabyte" << endl;
  int ii = 0;
  for(int i = 0; i < MAXMEMNUMBER; i++)
    if(mp[i].pointer)
      cout << " " << ii++ << "  " << i << "  " << mp[i] << "\n";
  cout << endl;
}

// Fortran version: FortranMemInfWrapper
extern "C" {
  void FORT(fortranmeminf)(const char *file, const FInt &len_file,
			   const int &line)
  {
    show_memory_allocation_information(Str(file, len_file), line);
  }
}

void free_loc(void * &p, const char *file_name, const int &line)
{
  if(!p) {
    cout << " *** Empty pointer to be freed: file: " << file_name 
	 << " line: " << line << endl;
    exit(1);
  }

  int have_found = 0;
  for(int i = 0; i < MAXMEMNUMBER; i++)
    if(p == mp[i].pointer) {
      have_found = 1;
      memory_in_use -= mp[i].memory_size();
      mp[i].deallocate();
      p = 0;
      break;
    }
  
  if(!have_found) {
    cout << " *** Failed to find pointer to free memory, file: " << file_name 
	 << " line: " << line << " pointer: " << p << endl;
    exit(1);
  }
}

void free_loc(char * &p, const char *file_name, const int &line)
{
  void *p_ = (void *) p;
  free_loc(p_, file_name, line);
  p = 0;
}

void free_loc(double * &p, const char *file_name, const int &line)
{
  void *p_ = (void *) p;
  free_loc(p_, file_name, line);
  p = 0;
}

void free_loc(int * &p, const char *file_name, const int &line)
{
  void *p_ = (void *) p;
  free_loc(p_, file_name, line);
  p = 0;
}

// Fortran version: FORTFreeLoc
extern "C" {
  void FORT(fortfreeloc)(void **p, 
			 const char *file_name, const FInt &len_file_name,
			 const FInt &line)
  {
    const int line_ = line;
    free_loc(*p, Str(file_name, len_file_name), line_);
  }
}

double *allocate_double_loc(const long int &n, const char *file_name, const int &line)
{
  int i;
  for(i = 0; i < MAXMEMNUMBER; i++)
    if(!mp[i].pointer) 
      break;  

  if(i == MAXMEMNUMBER) {
    cout << "Too many memory number: " << file_name << "  " << line << endl;
    exit(1);
  }
  
  mp[i].file = file_name;
  mp[i].line = line;
  mp[i].type = "DOUBLE";
  mp[i].size = n;
  
  mp[i].allocate();
  memory_in_use += mp[i].memory_size();
  return (double *) mp[i].pointer;
}

double *allocate_double_loc(const int &n, const char *file_name, const int &line)
{
  const long int n_ = n;
  return allocate_double_loc(n_, file_name, line);
}

// Fortran version: FORTAllocReal8Loc
extern "C" {
  void FORT(fortallocreal8loc)(double **p, const FInt &n,
			       const char *file_name_, const FInt &len_file_name_,
			       const FInt &line)
  {
    Str file_name(file_name_, len_file_name_);
    
    if(*p) {
      for(int i = 0; i < MAXMEMNUMBER; i++)
	if(*p == mp[i].pointer) {
	  cout << "*** Error memory allocation: memory has been allocated\n"
	       << "file: " << file_name << " line: " << line << endl
	       << mp[i] << endl;
	  exit(1);
	}
      
      cout << " *** Error memory allocation: non-zero pointer and no allocates memory pointer match\n"
	   << " file: " << file_name << " line: " << line << endl;
    }
    
    const long int n_ = n;
    const int line_ = line;
    
    *p = allocate_double_loc(n_, file_name, line_);
  }
}

int *allocate_int_loc(const long int &n, const char *file_name, const int &line)
{
  int i;
  for(i = 0; i < MAXMEMNUMBER; i++)
    if(!mp[i].pointer) 
      break;  

  if(i == MAXMEMNUMBER) {
    cout << " *** Too many memory number: " << file_name << "  " << line << endl;
    exit(1);
  }
  
  mp[i].file = file_name;
  mp[i].line = line;
  mp[i].type = "INT";
  mp[i].size = n;
  
  mp[i].allocate();
  memory_in_use += mp[i].memory_size();
  return (int *) mp[i].pointer;
}

int *allocate_int_loc(const int &n, const char *file_name, const int &line)
{
  const long int n_ = n;
  return allocate_int_loc(n_, file_name, line);
}

// Fortran version: FORTAllocIntLoc
extern "C" {
  void FORT(fortallocintloc)(FInt **p, const FInt &n,
			     const char *file_name_, const FInt &len_file_name_,
			     const FInt &line)
  {
    Str file_name(file_name_, len_file_name_);
    
    if(*p) {
      for(int i = 0; i < MAXMEMNUMBER; i++)
	if(*p == mp[i].pointer) {
	  cout << " *** Error memory allocation: memory has been allocated\n"
	       << "file: " << file_name << " line: " << line << endl
	       << mp[i] << endl;
	  exit(1);
	}
      
      cout << " *** Error memory allocation: non-zero pointer and no allocates memory pointer match\n"
	   << " file: " << file_name << " line: " << line << endl;

      exit(1);
    }
    
    int i;
    for(i = 0; i < MAXMEMNUMBER; i++)
      if(!mp[i].pointer) 
	break;
    
    if(i == MAXMEMNUMBER) {
      cout << " *** Too many memory number: " << file_name << "  " << line << endl;
      exit(1);
    }
    
    mp[i].file = file_name;
    mp[i].line = line;
    mp[i].type = "FINT";
    mp[i].size = n;
    
    mp[i].allocate();
    *p = (FInt *) mp[i].pointer;

    memory_in_use += mp[i].memory_size();
  }
}

char *allocate_char_loc(const long int &n, const char *file_name, const int &line)
{
  int i;
  for(i = 0; i < MAXMEMNUMBER; i++)
    if(!mp[i].pointer) 
      break;  
  
  if(i == MAXMEMNUMBER) {
    cout << " *** Too many memory number: " << file_name << "  " << line << endl;
    exit(1);
  }
  
  mp[i].file = file_name;
  mp[i].line = line;
  mp[i].type = "CHAR";
  mp[i].size = n;
  
  mp[i].allocate();
  memory_in_use += mp[i].memory_size();

  char *p = (char *) mp[i].pointer;
  p[0] = '\0';

  return p;
}

char *allocate_char_loc(const int &n, const char *file_name, const int &line)
{
  const long int n_ = n;
  return allocate_char_loc(n_, file_name, line);
}

// Fortran version: FORTAllocCharLoc
extern "C" {
  void FORT(fortalloccharloc)(char **p, const FInt &n,
			      const char *file_name_, const FInt &len_file_name_,
			      const FInt &line)
  {
    Str file_name(file_name_, len_file_name_);

    if(*p) {
      for(int i = 0; i < MAXMEMNUMBER; i++)
	if(*p == mp[i].pointer) {
	  cout << " *** Error memory allocation: memory has been allocated\n"
	       << "file: " << file_name << " line: " << line << endl
	       << mp[i] << endl;
	  exit(1);
	}
      
      cout << " *** Error memory allocation: non-zero pointer"
	   << " and no allocates memory pointer match\n"
	   << " file: " << file_name << " line: " << line << endl;
    }
    
    const long int n_ = n;
    const int line_ = line;
    *p = allocate_char_loc(n_, file_name, line_);
    memset(*p, 0, n*sizeof(char));
  }
}

