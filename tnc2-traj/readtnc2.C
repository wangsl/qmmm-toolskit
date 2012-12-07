
/* $Id$ */


#include "readtnc2.h"
#include "ReadNetCDFtraj.h"

static ReadNetCDFTrajectoryFile *read_tnc2 = 0;

void create_read_tnc2(const char *file_name)
{
  if(!read_tnc2)
    read_tnc2 = new ReadNetCDFTrajectoryFile(file_name);
  else {
    cout << " NetCDF MD trajectory has been created already with file name: " 
	 << read_tnc2->file_name() << endl;
    QCrash("Error");
  }
}

void delete_read_tnc2() 
{
  if(read_tnc2) {
    delete read_tnc2;
    read_tnc2 = 0; 
  }
}

int snapshot_size() 
{
  return read_tnc2 ? read_tnc2->n_snapshots : 0; 
}

int n_atoms()
{
  assert(read_tnc2);
  return read_tnc2->n_atoms;
}

const int *atom_types()
{
  assert(read_tnc2 && read_tnc2->atom_types);
  return read_tnc2->atom_types;
}

const char *atom_names()
{
  assert(read_tnc2 && read_tnc2->atom_names);
  return read_tnc2->atom_names;
}

const double *snapshot(const int i)
{
  assert(read_tnc2);
  read_tnc2->read_snapshot(i);
  return read_tnc2->atom_coordinates;
}

int n_active_atoms() 
{
  assert(read_tnc2);
  return read_tnc2->n_use;
}

const int *active_atom_list()
{
  assert(read_tnc2 && read_tnc2->active_atom_list);
  return read_tnc2->active_atom_list;
}

const double *qm_esp_charges()
{
  assert(read_tnc2);
  return read_tnc2->qm_esp_charges;
}

const double *gradients()
{
  assert(read_tnc2);
  return read_tnc2->gradients;
}

const double *velocities()
{
  assert(read_tnc2);
  return read_tnc2->velocities;
}

int qm_esp_charges_saved()
{
  assert(read_tnc2);
  return read_tnc2->have_esp_charges_saved();
}

int gradients_saved()
{
  assert(read_tnc2);
  return read_tnc2->have_gradients_saved();
}

int velocities_saved()
{
  assert(read_tnc2);
  return read_tnc2->have_velocities_saved();
}

int n_total_qm_atoms()
{
  assert(read_tnc2);
  return read_tnc2->n_total_qm_atoms();
}
