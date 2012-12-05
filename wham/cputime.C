
#include <time.h>
#include "fort.h"

double cpu_time()
{
  clock_t start;
  double cpu;
  start = clock();
  /* cpu = ((double) start) / CLOCKS_PER_SEC; */
  cpu = ((double) start) * 1.0e-6;
  return cpu;
}

// Fortran version: CPUTime
extern "C" {
  double FORT(cputime)()
  {
    return cpu_time();
  }
}
