/* $Id: out.h 2 2007-10-20 14:59:10Z wangsl $ */

#ifndef OUT_H
#define OUT_H

#include <iostream>
using namespace std;

void OutInit();
void OutSetToFile(const char *);
void SetPrecision(int);
void OutSetToFixed();
void OutUnsetToFixed();
ostream &Out();
ostream *FileStream(const char *);
istream *InFileStream(const char *);

#endif /* OUT_H */
