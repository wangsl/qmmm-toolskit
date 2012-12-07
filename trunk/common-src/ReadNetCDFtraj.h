
/* $Id$ */

#ifndef READNETCDFTRAJ_H
#define READNETCDFTRAJ_H

#include "coord.h"
#include "str.h"
#include "NetCDF.h"
#include "NCtrajdata.h"

class NetCDFTrajectoryAtom
{
public:
  NetCDFTrajectoryAtom();

  NetCDFTrajectoryAtom(int atom_index_, int atom_type_index_,
		       const Str &atom_name_, const Cartesian &coordinate_);

  ~NetCDFTrajectoryAtom();

  void set_esp_charge(const double esp_charge_) 
  {
    if(!esp_charge) esp_charge = new double;
    *esp_charge = esp_charge_;
  }
  
  void set_gradient(const Cartesian &gradient_) 
  {
    if(!gradient) gradient = new Cartesian;
    *gradient = gradient_;
  }

  void write(ostream &fout) const;
  
private:

  int atom_index; // out
  int atom_type_index; // out
  Str atom_name; // out
  Cartesian coordinate; // out

  Cartesian *gradient; // out
  double *esp_charge; // out

private:
  const static char *atom_type [];

  // copy constructor to prevent pass-by-value
  NetCDFTrajectoryAtom(const NetCDFTrajectoryAtom &);

  friend ostream & operator <<(ostream &s, const NetCDFTrajectoryAtom &c);
  void write_fields(ostream &s) const;  
};

class ReadNetCDFTrajectoryFile
{
public:
  ReadNetCDFTrajectoryFile(const char *MD_trajectory_file_name);
  
  ~ReadNetCDFTrajectoryFile();

  void read_and_write_all_snapshots();
  void read_snapshot(const int irecord);

  int n_atoms;
  int n_snapshots;
  double *atom_coordinates;
  char *atom_names;
  int *atom_types;

  int n_use;
  int *active_atom_list;
  double *qm_esp_charges;
  double *gradients; 
  double *velocities;

  int have_esp_charges_saved() const 
  { return save_esp_charges ? 1 : 0; }

  int have_gradients_saved() const
  { return save_gradients ? 1 : 0; }

  int have_velocities_saved() const
  { return save_velocities ? 1 : 0; }

  const char *file_name() const
  { return MD_trajectory_file ? MD_trajectory_file->NetCDF_file_name() : 0; }

  int n_total_qm_atoms() const
  { return  n_qm_atoms + n_pseudo_atoms; }

private:
  
  NetCDFFile *MD_trajectory_file;

  int save_esp_charges; 
  int save_gradients;
  int save_velocities;
  int snapshot_size;  

  //int n_use;
  int n_qm_atoms;
  int n_pseudo_atoms;
  int n_zero_charge_atoms;
  int n_mm_atoms;

  //int *active_atom_list;
  int *use;

  void initialize();

  void check_dimensions_and_variables();

  void calculate_snapshot_size();
  void check_snapshot_information();

  void read_qmmm_groups();
  void read_active_atom_list();
  void read_fixed_atom_coordinates();
  void read_atom_names();
  void read_atom_types();

  void write_full_snapshot(ostream &fout, const double *snapshot);

private:
  // copy constructor to prevent pass-by-value
  ReadNetCDFTrajectoryFile(const ReadNetCDFTrajectoryFile &);

  void create_full_atom_coordinates(const double *snapshot);

  void check_version() const;
};

#endif /* READNETCDFTRAJ_H */
