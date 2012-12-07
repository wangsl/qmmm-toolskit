
/* $Id$ */

/*
  QChem-Tinker NetCDF file .tnc, which means Tinker NetCDF file
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
#include "NetCDF.h"

typedef struct {
  NetCDFFile *MD_traj;
  NcVar *xyz;
  NcVar *xyz_esp;
  int natoms;
  int have_esp_charges;
  int nrecords;
  int record_size;
} Tinker_NetCDF_data;

static void add_char_end(char *word, int len)
{
  for(int i = len-1; i >= 0; i--) {
    if(!isspace(word[i])) {
      word[i+1] = '\0';
      return;
    }
  }
}

static void *open_Tinker_NetCDF(const char *file_name, const char *filetype, int *natoms)
{
  cout << "**********************************************\n"
       << " QChem-Tinker NetCDF trajectory plugin by \n"
       << "   Shenglong Wang and Yingkai Zhang\n"
       << " Department of Chemistry, New York University\n"
       << "**********************************************"
       << endl;
  
  Tinker_NetCDF_data *NetCDF_data = (Tinker_NetCDF_data *) malloc(sizeof(Tinker_NetCDF_data));
  assert(NetCDF_data);

  NetCDF_data->MD_traj = new NetCDFFile(file_name, NCREADONLY);
  assert(NetCDF_data->MD_traj);
  
  NetCDF_data->xyz = 0;
  NetCDF_data->xyz_esp = 0;
  NetCDF_data->natoms = 0;
  NetCDF_data->have_esp_charges = 0;
  NetCDF_data->nrecords = 0;
  NetCDF_data->record_size = 0;
  
  cout << " QChem-Tinker NetCDF trajectory file name: " 
       << NetCDF_data->MD_traj->NetCDF_file_name() << endl;

  int *natoms_ = new int [7];
  assert(natoms_);
  NetCDF_data->MD_traj->get_int_var("natoms", natoms_);
  
  cout << " There are " << natoms_[0] << " atoms in total\n"
       << " " << natoms_[1] << " QM atoms\n"
       << " " << natoms_[2] << " Pseudo atoms\n"
       << " " << natoms_[3] << " Zero charge atoms\n"
       << " " << natoms_[4] << " MM atoms\n"
       << " QM total charge " << natoms_[5] << "\n"
       << " QM spin multiplicity " << natoms_[6] << endl;
  
  *natoms = natoms_[0];
  NetCDF_data->natoms = *natoms;
  
  const int save_esp_charges = NetCDF_data->MD_traj->int_var("save_esp_charges");
  if(save_esp_charges) {
    cout << " ESP charges have been saved in NetCDF trajectory file" << endl;
    NetCDF_data->have_esp_charges = 1;
    NetCDF_data->xyz_esp = NetCDF_data->MD_traj->var("coordinates_esp_charges");
    NetCDF_data->nrecords = NetCDF_data->xyz_esp->get_dim(0)->size();
    NetCDF_data->record_size = NetCDF_data->xyz_esp->get_dim(1)->size();
    assert(NetCDF_data->record_size == 3*natoms_[0] + natoms_[1] + natoms_[2]);
  } else {
    cout << " ESP charges have not been saved in NetCDF trajectory file" << endl;
    NetCDF_data->have_esp_charges = 0;
    NetCDF_data->xyz = NetCDF_data->MD_traj->var("coordinates");
    NetCDF_data->nrecords = NetCDF_data->xyz->get_dim(0)->size();
    NetCDF_data->record_size = NetCDF_data->xyz->get_dim(1)->size();
    assert(NetCDF_data->record_size == 3*natoms_[0]);
  }

  cout << " There are " << NetCDF_data->nrecords << " snapshots have been saved,"
       << " each record size is " << NetCDF_data->record_size << endl;

  if(natoms_) { delete natoms_; natoms_ = 0; }

  return NetCDF_data;
}

static void close_Tinker_NetCDF(void *my_data) 
{ 
  Tinker_NetCDF_data *data = (Tinker_NetCDF_data *) my_data;

  if(data->MD_traj) {
    delete data->MD_traj;
    data->MD_traj = 0;
  }

  if(data->xyz) { data->xyz = 0; }
  if(data->xyz_esp) { data->xyz_esp = 0; }
    
  if(data) { free(data); data = 0; }
}

static int read_Tinker_NetCDF_structure(void *my_data, int *opt_flags, molfile_atom_t *atoms)
{
  Tinker_NetCDF_data *data = (Tinker_NetCDF_data *) my_data;

  *opt_flags = MOLFILE_NOOPTIONS;

  const int natoms = data->natoms;

  char *atom_names = new char [3*natoms+1];
  assert(atom_names);
  assert(data->MD_traj->var("atom_names")->get(atom_names, 3*natoms));
  atom_names[3*natoms] = '\0';

  int *atom_types = new int [3*natoms];
  assert(atom_types);
  data->MD_traj->get_int_var("atom_types", atom_types);

  for(int i = 0; i < natoms; i++) {
    molfile_atom_t *atom = atoms + i;
    memcpy(atom->name, atom_names+3*i, 3*sizeof(char));
    add_char_end(atom->name, 3);
    sprintf(atom->type, "%d", *(atom_types+i));
    atom->resname[0] = '\0';
    atom->resid = 1;
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';
  }

  if(atom_names) { delete [] atom_names; atom_names = 0; }
  if(atom_types) { delete [] atom_types; atom_types = 0; }

  return MOLFILE_SUCCESS;
}

static int read_Tinker_NetCDF_timestep(void *my_data, int natoms, molfile_timestep_t *ts)
{
  static int snapshot_index = 0;
  snapshot_index++;
  
  Tinker_NetCDF_data *data = (Tinker_NetCDF_data *) my_data;
  
  if(snapshot_index > data->nrecords) 
    return MOLFILE_ERROR;
  
  assert(natoms == data->natoms);

  NcVar *snaptshot_var = data->xyz;
  if(data->xyz_esp) snaptshot_var = data->xyz_esp;
  
  NcValues *record = snaptshot_var->get_rec(long(snapshot_index-1));
  double *xyz_ = (double *) record->base();
  
  if(ts) {
    const int &natoms = data->natoms;
    for(int i = 0; i < 3*natoms; i++) 
      ts->coords[i] = xyz_[i];
  }

  if(record) { delete record; record = 0; }
  
  return MOLFILE_SUCCESS;
}
 
/* registration stuff */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "tnc";
  plugin.prettyname = "Tinker-NetCDF";
  plugin.author = "Shenglong Wang";
  plugin.majorv = 0;
  plugin.minorv = 0;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "tnc";
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

