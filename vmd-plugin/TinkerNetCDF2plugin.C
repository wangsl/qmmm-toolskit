
/*
  QChem-Tinker NetCDF file .tnc2, which means Tinker NetCDF file
  by Shenglong Wang
  Department of Chemistry, New York University
*/

#include <iostream>
using namespace std;

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctype.h>
#include "molfile_plugin.h"
#include "ReadNetCDFtraj.h"

static inline void add_char_end(char *word, int len)
{
  for(int i = len-1; i >= 0; i--) {
    if(!isspace(word[i])) {
      word[i+1] = '\0';
      return;
    }
  }
}

static void *open_Tinker_NetCDF(const char *file_name, const char *filetype, int *n_atoms)
{
  cout << "***********************************************\n"
       << " QChem-Tinker NetCDF MD trajectory plugin by \n"
       << "   Shenglong Wang and Yingkai Zhang\n"
       << "     Version 2\n" 
       << " Department of Chemistry, New York University\n"
       << "***********************************************"
       << endl;

  ReadNetCDFTrajectoryFile *NC_traj_file = new ReadNetCDFTrajectoryFile(file_name);
  assert(NC_traj_file);
  
  *n_atoms = NC_traj_file->n_atoms;
  
  return NC_traj_file;
}

static void close_Tinker_NetCDF(void *my_data) 
{
  ReadNetCDFTrajectoryFile *NC_traj_file = (ReadNetCDFTrajectoryFile *) my_data;
  if(NC_traj_file) {
    delete NC_traj_file;
    NC_traj_file = 0;
  }
}

static int read_Tinker_NetCDF_structure(void *my_data, int *opt_flags, molfile_atom_t *atoms)
{
  ReadNetCDFTrajectoryFile *NC_traj_file = (ReadNetCDFTrajectoryFile *) my_data;

  *opt_flags = MOLFILE_NOOPTIONS;

  const int n_atoms = NC_traj_file->n_atoms;
  assert(n_atoms > 0);

  const char *atom_names = NC_traj_file->atom_names;
  assert(atom_names);

  const int *atom_types = NC_traj_file->atom_types;
  assert(atom_types);

  for(int i = 0; i < n_atoms; i++) {
    molfile_atom_t *atom = atoms + i;
    memcpy(atom->name, atom_names+3*i, 3*sizeof(char));
    add_char_end(atom->name, 3);
    sprintf(atom->type, "%d", *(atom_types+i));
    atom->resname[0] = '\0';
    atom->resid = 1;
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';
  }

  return MOLFILE_SUCCESS;
}

static int read_Tinker_NetCDF_timestep(void *my_data, int n_atoms, molfile_timestep_t *ts)
{
  static int snapshot_index = 0;
  snapshot_index++;
  
  ReadNetCDFTrajectoryFile *NC_traj_file = (ReadNetCDFTrajectoryFile *) my_data;
  assert(n_atoms == NC_traj_file->n_atoms);
  
  if(snapshot_index > NC_traj_file->n_snapshots)
    return MOLFILE_ERROR;
  
  NC_traj_file->read_snapshot(snapshot_index-1);
  const double *xyz = (double *) NC_traj_file->atom_coordinates;
  
  if(ts) {
    for(int i = 0; i < 3*n_atoms; i++) 
      ts->coords[i] = xyz[i];
  }
  
  return MOLFILE_SUCCESS;
}

/* registration stuff */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "tnc2";
  plugin.prettyname = "Tinker-NetCDF-2";
  plugin.author = "Shenglong Wang";
  plugin.majorv = 0;
  plugin.minorv = 0;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "tnc2";
  plugin.open_file_read = open_Tinker_NetCDF; // open_tinker_read;
  plugin.read_structure = read_Tinker_NetCDF_structure; // read_tinker_structure;
  plugin.read_next_timestep = read_Tinker_NetCDF_timestep; //read_tinker_timestep;
  plugin.close_file_read = close_Tinker_NetCDF; //close_tinker_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

#ifdef TEST_PLUGIN

int main(int argc, char *argv[]) {
  molfile_timestep_t timestep;
  void *v;
  int natoms;
  int i, nsets, set;

  while (--argc) {
    ++argv;
    v = open_tinker_read(*argv, "tinker", &natoms);
    if (!v) {
      fprintf(stderr, "open_tinker_read failed for file %s\n", *argv);
      return 1;
    }
    fprintf(stderr, "open_tinker_read succeeded for file %s\n", *argv);
    fprintf(stderr, "number of atoms: %d\n", natoms);

    i = 0;
    timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
    while (!read_tinker_timestep(v, natoms, &timestep)) {
      i++;
    }
    fprintf(stderr, "ended read_next_timestep on frame %d\n", i);

    close_tinker_read(v);
  }
  return 0;
}

#endif

