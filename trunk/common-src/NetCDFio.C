#include <iostream>
using namespace std;
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "indent.h"
#include "NetCDF.h"
#include "die.h"

ostream & operator <<(ostream &s, const NetCDFFile &c)
{
  s << "{\n";
  IndentPush();
  c.write_fields(s);
  IndentPop();
  return s << Indent() << "} ";
}

void NetCDFFile::write_fields(ostream &s) const
{
  if (file_name)
    s << Indent() << "file_name " << file_name << "\n";
  s << Indent() << "file_mode " << file_mode << "\n";
  s << Indent() << "dimensions " << dimensions << "\n";
  s << Indent() << "variables " << variables << "\n";
}

