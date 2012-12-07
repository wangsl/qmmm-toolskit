
/* $Id$ */

#ifndef DIE_H
#define DIE_H

#include <cstdlib>

void die_at(const char *s, const char *file, int line);
void die(const char *s);

#ifdef QCHEM
#include "tinker.h"

#else // QCHEM

#include "fort.h"
inline void FORT(fatal)() { exit(13); }
#define insist(x) if (!(x)) die_at("insist failed: " #x, __FILE__, __LINE__)
#define QCrash(x) die_at(x, __FILE__, __LINE__)

#endif /* QCHEM */

#endif /* DIE_H */
