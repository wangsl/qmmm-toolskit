 
/* $Id: die.C 2 2008-02-28 22:52:52Z wangsl $ */

#include <unistd.h>
#include <cstdlib>
#include <iostream>
using namespace std;
#include "die.h"

void die_at(const char *s, const char *file, int line)
{
  cout << s << " at " << file << ":" << line << "\n" << flush;
  exit(13);
}

/*
void die(const char *s)
{
  cout << s << "\n" << flush;
  exit(13);
}
*/



