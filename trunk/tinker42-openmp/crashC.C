 
/* $Id: crashC.C 2 2008-02-28 22:52:52Z wangsl $ */

#include "tinkercpp.h"

void CrashLoc(const char *message, const char *module, const int line)
{
  cout << " Tinker error occurred in module " << module << ", line " << line << "\n"
       << "  " << message << endl;
  FORT(fatal)();
}
