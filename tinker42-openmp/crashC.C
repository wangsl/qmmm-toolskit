 
/* $Id$ */

#include "tinkercpp.h"

void CrashLoc(const char *message, const char *module, const int line)
{
  cout << " Tinker error occurred in module " << module << ", line " << line << "\n"
       << "  " << message << endl;
  FORT(fatal)();
}
