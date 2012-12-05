

#include <iostream>
using namespace std;
#include <netcdfcpp.h>

#include "NetCDF.h"
#include "ReadNetCDFtraj.h"
#include "vec.h"
#include "cvec.h"

void add_char_end(char *word, int len)
{
  for(int i = len-1; i >= 0; i--) {
    if(!isspace(word[i])) {
      word[i+1] = '\0';
      return;
    }
  }
}

int main()
{
  cout << "\n Test by Shenglong Wang: " << __FILE__<< "  " << __LINE__ 
       << "\n" << endl;

  ReadNetCDFTrajectoryFile ReadNC("CH3ClCl-MD.tnc");
  
  cout << " n_snapshot = " << ReadNC.n_snapshots << endl;
  cout << " n_atoms = " << ReadNC.n_atoms << endl;

  for(int i = 0; i < ReadNC.n_snapshots; i++) {
    ReadNC.read_snapshot(i);
    cout << i << endl
	 << CVec(ReadNC.n_atoms, ReadNC.atom_coordinates) 
	 << endl;
  }
  
  return 0;
}

