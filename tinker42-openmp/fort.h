 
/* $Id: fort.h 2 2008-02-28 22:52:52Z wangsl $ */

#ifndef FMACRO_H
#define FMACRO_H

#ifdef NO_UNDERSCORE
#define FORT(x) x
#else
#define FORT(x) x##_
#endif

#endif /* FORT_H */
