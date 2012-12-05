#include <iostream>
using namespace std;
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "indent.h"
#include "ReadNetCDFtraj.h"
#include "die.h"
#include "out.h"

ostream & operator<<(ostream &s, const NetCDFTrajectoryAtom &c)
{
  s << "{\n";
  IndentPush();
  c.write_fields(s);
  IndentPop();
  return s << Indent() << "}\n";
}

void NetCDFTrajectoryAtom::write_fields(ostream &s) const
{
  s << Indent() << "atom_index " << atom_index << "\n";
  s << Indent() << "atom_type_index " << atom_type_index << "\n";
  if (strlen(atom_name) > 0)
    s << Indent() << "atom_name " << atom_name << "\n";
  s << Indent() << "coordinate " << coordinate << "\n";
  if (gradient)
    s << Indent() << "gradient " << *gradient << "\n";
  if (esp_charge)
    s << Indent() << "esp_charge " << *esp_charge << "\n";
}

