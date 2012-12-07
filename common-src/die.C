
/* $Id: die.C 147 2009-03-10 19:15:31Z wangsl $ */

#include <unistd.h>
#include <cstdlib>
#include <iostream>
using namespace std;
#include "die.h"

void die_at(const char *s, const char *file, int line)
{
  cout << " " << s << " at file: '" << file << "' line: " << line << "\n" << flush;
  FORT(fatal)();
}

void die(const char *s)
{
  cout << s << "\n" << flush;
  FORT(fatal)();
}



