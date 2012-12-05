
#ifndef WHAM_H
#define WHAM_H

#include "fort.h"
#include "fint.h"

void FORT(wham_iteration_with_fortran)(double *ebw, FInt *Nwind, FInt *Niter, 
				       double *tolerance, FInt *nt, FInt *maxnt, 
				       double *kBT, double *ebf, double *ebf2,
				       char *outtput, FInt *lenoutput);

#endif /* WHAM_H */
