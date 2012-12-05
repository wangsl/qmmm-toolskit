
/**********************************************
  Note by Shenglong Wang, 2-13-2006
  the Integer variables for 64-bits compiler 
  -I64 or -qintsize=8

  Fortran sizef(integer) = 8
  Fortran: sizeof(Logical) = 8
  Fortran: sizeof(Real*8) = 8
  Fortran: sizeof(Double Precision) = 16

  C++: sizeof(long) = 8
  C++: sizeof(int) = 4
  C++: sizeof((void *) 0) = 8, pointer

**********************************************/

#ifndef FINT_H
#define FINT_H

#ifdef I64
typedef long int f77_int;
typedef long int f77int;
typedef long int FtranInt;
typedef long int F77Int;
typedef long int FInt8;
typedef long int FInt;

#define PyArray_F77INT PyArray_LONG
#define PyArray_FINT PyArray_LONG

#else

typedef int f77_int;
typedef int f77int;
typedef int FtranInt;
typedef int F77Int;
typedef int FInt8;
typedef int FInt;

#define PyArray_F77INT PyArray_INT
#define PyArray_FINT PyArray_INT

#endif

#endif /* FINT_H */
