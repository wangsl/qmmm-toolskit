 
/* $Id$ */

#ifndef TINKERCPP_H
#define TINKERCPP_H

#include <iostream>
using namespace std;

#include "fort.h"

#define QCrash(x) Crash(x)
#define Crash(x) CrashLoc(x, __FILE__, __LINE__)

extern "C" {
  void FORT(fatal)();
}

void CrashLoc(const char *message, const char *module, const int line);

#endif /* TINKERCPP_H */
