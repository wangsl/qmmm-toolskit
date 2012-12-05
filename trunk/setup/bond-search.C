
/* $Id: bond-search.C 4 2012-12-05 04:19:17Z wangsl2001@gmail.com $ */


#include <iostream>
#include <fstream>
using namespace std;
#include <cassert>

// The bond serach method is same as in VMD source code

#ifndef MAXNBONDS 
#define MAXNBONDS 10
#endif

inline double sq(const double x) { return x*x; }

inline double dist2(const double *r1, const double *r2)
{
  return sq(r1[0]-r2[0]) + sq(r1[1]-r2[1]) + sq(r1[2]-r2[2]);
}

int main(int argc, char *argv[])
{
  const char *input_file = argv[1];
  const char *output_file = argv[2];

  cout << " Input file: " << input_file << "\n"
       << " Output file: " << output_file << endl;

  char buf[1024];
  FILE *fin = fopen(input_file, "r");
  assert(fin);

  assert(fgets(buf, 1024, fin) && buf[strlen(buf)-1] == '\n');
  int n_atoms = 0;
  assert(sscanf(buf,"%d", &n_atoms) == 1);
  cout << " Total atom number: " << n_atoms << endl;
  
  int *atomic_number = new int [n_atoms];
  assert(atomic_number);
  double *xyz = new double [3*n_atoms];
  assert(xyz);
  double *radii = new double [n_atoms];
  assert(radii);
  
  int atom_index = 0;
  for(int i = 0; i < n_atoms; i++) {
    assert(fgets(buf, 1024, fin) && buf[strlen(buf)-1] == '\n');
    assert(sscanf(buf,"%d%d%lf%lf%lf%lf", &atom_index, &atomic_number[i],
		  &xyz[3*i], &xyz[3*i+1], &xyz[3*i+2], &radii[i]) == 6);
    assert(atom_index == i+1);
  }
  fclose(fin);
  
  int *bond_number = new int [n_atoms];
  memset(bond_number, 0, n_atoms*sizeof(int));
  int *bond_list = new int [n_atoms*MAXNBONDS];
  memset(bond_list, 0, MAXNBONDS*n_atoms*sizeof(int));

  int n_bonds = 0;
  for(int i = 0; i < n_atoms; i++) {
    const int atomic_number_i = atomic_number[i];
    const double rad_i = radii[i];
    const double *coord_i = xyz+3*i;
    int &n_bond_i = bond_number[i];
    int *bond_list_i = bond_list + i*MAXNBONDS;
    for(int j = i+1; j < n_atoms; j++) {
      const int atomic_number_j = atomic_number[j];
      if(atomic_number_i == 1 && atomic_number_j == 1)
	continue;
      const double rad_j = radii[j];
      const double *coord_j = xyz+3*j;
      int &n_bond_j = bond_number[j];
      int *bond_list_j = bond_list + j*MAXNBONDS;
      const double cutoff = 0.60*(rad_i + rad_j);
      const double rij2 = dist2(coord_i, coord_j);
      if(rij2 < cutoff*cutoff) {
	bond_list_i[n_bond_i++] = j+1;
	bond_list_j[n_bond_j++] = i+1;
	n_bonds++;
      }
    }
  }

  int n_bonds_ = 0;
  for(int i = 0; i < n_atoms; i++) {
    if(bond_number[i] < MAXNBONDS)
      n_bonds_ += bond_number[i];
    else {
      cout << " Bond number error for atom index: " << i+1 << "\n"
	   << " MAXNBONDS = " << MAXNBONDS << "\n"
	   << " Please use -DMAXNBONDS=?? to compile source code: " << __FILE__ << endl; 
      exit(13);
    }
  }
  assert(2*n_bonds == n_bonds_);
  cout << " There are " << n_bonds << " bonds have been constructed" << endl;
  
  if(atomic_number) { delete [] atomic_number; atomic_number = 0; }
  if(xyz) { delete [] xyz; xyz = 0; }
  if(radii) { delete [] radii; radii = 0; }

  ofstream fout(output_file, ios::out);
  fout << n_atoms << "\n";
  for(int i = 0; i < n_atoms; i++) {
    const int n_bond = bond_number[i];
    const int *bond = bond_list + i*MAXNBONDS;
    fout << i+1 << "  ";
    for(int ib = 0; ib < n_bond; ib++) 
      fout << bond[ib] << " ";
    fout << "\n";
  }
  fout.close();

  if(bond_number) { delete [] bond_number; bond_number = 0; }
  if(bond_list) { delete [] bond_list; bond_list = 0; }

  cout << " Bond list contruction finished" << endl;

  return 0;
}
