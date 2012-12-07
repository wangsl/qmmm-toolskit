 
/* $Id$ */

#include <cassert>

#include "fint.h"
#include "memptr.h"

MemoryPointer::MemoryPointer() : pointer(0), line(0), size(0)
{ }

MemoryPointer::MemoryPointer(void *pointer_, Str file_, int line_,
			     Str type_, long int size_) :
  pointer(pointer_), file(file_), line(line_), type(type_),
  size(size_)
{ }

MemoryPointer::~MemoryPointer()
{
  deallocate();
}

void MemoryPointer::deallocate()
{
  if(pointer) {
    free(pointer);
    pointer = 0;
  }
  file = "";
  type = "";
  line = 0;
  size = 0;
}

void MemoryPointer::allocate()
{
  if(pointer) {
    cout << *this << endl;
    die("pointer is not empty");
  }
  
  pointer = malloc(memory_size());
  assert(pointer);
}

long int MemoryPointer::memory_size() const
{
  long int n = 0;
  if(!strcmp(type, "INT"))
    n = size*sizeof(int);
  else if(!strcmp(type, "FINT")) 
    n = size*sizeof(FInt);
  else if(!strcmp(type, "DOUBLE"))
    n = size*sizeof(double);
  else if(!strcmp(type, "CHAR"))
    n = size*sizeof(char);
  else {
    cout << *this << endl;
    die("Unknow type");
  }
  return n;
}

ostream & operator <<(ostream &s, const MemoryPointer &c)
{
  c.write_fields(s);
  return s;
}

void MemoryPointer::write_fields(ostream &s) const
{
  if (pointer) {
    s << "pointer: " << pointer << ", ";
    if (strlen(file) > 0)
      s << "file: " << file << ", ";
    s << "line: " << line << ", ";
    if (strlen(type) > 0)
      s << "type: " << type << ", ";
    s << "size: " << size << ", ";
    s << "memory size: " << memory_size();
  }
}
  
