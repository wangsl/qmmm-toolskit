

/* $Id: fort.h 175 2010-09-06 15:19:30Z wangsl $ */

#ifndef FMACRO_H
#define FMACRO_H

#ifdef NO_UNDERSCORE
#define FORT(x) x
#else
#define FORT(x) x##_
#endif

#endif /* FORT_H */
