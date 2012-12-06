 
/* $Id: memptr.h 2 2008-02-28 22:52:52Z wangsl $ */

#ifndef MEMPTR_H
#define MEMPTR_H

#include <iostream>
using namespace std;

#include "str.h"

class MemoryPointer
{
public:
  void *pointer; // out
  Str file; // out
  int line; // out
  Str type; // out
  long int size; // out
  
  MemoryPointer();
  MemoryPointer(void *pointer_, Str file_, int line_,
		Str type_, long int size_);
  
  ~MemoryPointer();
  
  void allocate();
  void deallocate();
  
  long int memory_size() const;

private:

  /* IO */
  friend ostream & operator <<(ostream &s, const MemoryPointer &c);
  void write_fields(ostream &s) const;
};

#endif /* MEMPTR_H */
