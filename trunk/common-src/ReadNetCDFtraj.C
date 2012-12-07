
/* $Id: ReadNetCDFtraj.C 158 2009-11-28 19:33:18Z wangsl $ */

#include <fstream>
#include "ReadNetCDFtraj.h"
#include "NCrestartdata.h"
#include "utility01.h"
#include "die.h"

// Macros for memory allocation 

#ifdef QCHEM /* memory allocation */

#define AllocateInt(n) QAllocInt(n)
#define AllocateDouble(n) QAllocDouble(n)
#define AllocateChar(n) QAllocChar(n)

#define DeAllocate(x) if(x) { QFree(x); x = 0; }

#else

#define AllocateInt(n) new int [n]
#define AllocateDouble(n) new double [n]
#define AllocateChar(n) new char [n]

#define DeAllocate(x) if(x) { delete [] x; x = 0; }

#endif /* QCHEM, memory allocation */

#define _QMAtom_  0
#define _PseudoAtom_  1
#define _ZeroChargeAtom_ 2
#define _FreeMMAtom_ 3
#define _FixedMMAtom_ 4

const char *NetCDFTrajectoryAtom::atom_type [] = {"QM", "PS", "ZO", "FR", "FX"};

static inline void add_char_end(char *word, const int len)
{
  for(int i = len-1; i >= 0; i--) {
    if(!isspace(word[i])) {
      word[i+1] = '\0';
      return;
    }
  }
}

// function definitions for NetCDFTrajectoryAtom

NetCDFTrajectoryAtom::NetCDFTrajectoryAtom() :
  atom_index(0), atom_type_index(0), gradient(0), esp_charge(0)
{ }

NetCDFTrajectoryAtom::NetCDFTrajectoryAtom(int atom_index_, int atom_type_index_,
					   const Str &atom_name_, 
					   const Cartesian &coordinate_) :
  atom_index(atom_index_), atom_type_index(atom_type_index_),
  atom_name(atom_name_), coordinate(coordinate_),
  gradient(0), esp_charge(0)
{ }

NetCDFTrajectoryAtom::~NetCDFTrajectoryAtom() 
{
  if(gradient) { delete gradient; gradient = 0; }
  if(esp_charge) { delete esp_charge; esp_charge = 0; }
}

void NetCDFTrajectoryAtom::write(ostream &fout) const
{
  char *line = AllocateChar(1024);
  assert(line);
  memset(line, 0, 1024*sizeof(char));
  //                    1234567890123456789
  const char NA12 [] = "    --NA--  ";
  const char NA14 [] = "    --NA--    ";

  sprintf(line, "%6d  %-3s%12.6f%12.6f%12.6f  %s  ", 
	  atom_index, (const char *) atom_name,
	  coordinate.x, coordinate.y, coordinate.z, atom_type[atom_type_index]);

  if(esp_charge) {
    sprintf(&line[strlen(line)], "%12.4f", *esp_charge);
  } else {
    sprintf(&line[strlen(line)], "%s", NA12);
  }
  
  if(gradient) {
    sprintf(&line[strlen(line)], "%14.6e%14.6e%14.6e", 
	    gradient->x, gradient->y, gradient->z);
  } else {
    sprintf(&line[strlen(line)], "%s%s%s", NA14, NA14, NA14);
  }

  fout << line << endl;

  DeAllocate(line);
}

// function definitions for ReadNetCDFTrajectoryFile

ReadNetCDFTrajectoryFile::ReadNetCDFTrajectoryFile(const char *MD_trajectory_file_name) :
  n_atoms(0), n_snapshots(0),
  atom_coordinates(0), atom_names(0), atom_types(0),
  n_use(0), active_atom_list(0), qm_esp_charges(0), gradients(0), velocities(0),
  MD_trajectory_file(0), 
  save_esp_charges(0), save_gradients(0), save_velocities(0),
  snapshot_size(0),
  n_qm_atoms(0), n_pseudo_atoms(0), n_zero_charge_atoms(0), n_mm_atoms(0),
  use(0)
{
  const int exist = file_exist(MD_trajectory_file_name);
  if(!exist) {
    cout << " NetCDF MD trajectory file " << MD_trajectory_file_name 
	 << " does not exist" << endl;
    QCrash(" NetCDF MD trajectory file doe not exist");
  }
  
  MD_trajectory_file = new NetCDFFile(MD_trajectory_file_name, NCREADONLY);
  assert(MD_trajectory_file);

  cout << " NetCDF MD trajectory file: " << file_name() << endl;

  initialize();
}

ReadNetCDFTrajectoryFile::~ReadNetCDFTrajectoryFile()
{
  if(MD_trajectory_file) {
    delete MD_trajectory_file;
    MD_trajectory_file = 0;
  }

  DeAllocate(atom_coordinates);
  DeAllocate(atom_names);
  DeAllocate(atom_types);
  DeAllocate(active_atom_list);
  DeAllocate(use);

  qm_esp_charges = 0;
  gradients = 0;
}
  
void ReadNetCDFTrajectoryFile::check_dimensions_and_variables()
{
  cout << " QChem-Tinker will check NetCDF dimensions and variables" << endl;

  assert(MD_trajectory_file->dimension_exist(_One_) &&
	 MD_trajectory_file->dim(_One_)->size() == 1);

  assert(MD_trajectory_file->dimension_exist(_Seven_) &&
	 MD_trajectory_file->dim(_Seven_)->size() == 7);

  assert(MD_trajectory_file->dimension_exist(_NAtoms_));
  n_atoms = MD_trajectory_file->dim(_NAtoms_)->size();
  cout << " There are " << n_atoms << " atoms in total" << endl; 

  assert(MD_trajectory_file->dimension_exist(_NUse_));
  n_use = MD_trajectory_file->dim(_NUse_)->size();
  cout << " There are " << n_use << " active atoms" << endl;

  assert(MD_trajectory_file->dimension_exist(_NAtoms3_) &&
	 MD_trajectory_file->dim(_NAtoms3_)->size() == 3*n_atoms);

  assert(MD_trajectory_file->dimension_exist(_NFixedAtoms3_) &&
	 MD_trajectory_file->dim(_NFixedAtoms3_)->size() == 3*(n_atoms - n_use));

  assert(MD_trajectory_file->variable_exist(_SaveESPCharges_));
  save_esp_charges = MD_trajectory_file->int_var(_SaveESPCharges_);
  if(save_esp_charges) 
    cout << " ESP charges have been saved" << endl;

  assert(MD_trajectory_file->variable_exist(_SaveGradients_));
  save_gradients = MD_trajectory_file->int_var(_SaveGradients_);
  if(save_gradients)
    cout << " Gradients have been saved" << endl;

  save_velocities = 0;
  if(MD_trajectory_file->variable_exist(_SaveVelocities_))
    save_velocities = MD_trajectory_file->int_var(_SaveVelocities_);
  if(save_velocities)
    cout << " Velocities have been saved" << endl;
}

void ReadNetCDFTrajectoryFile::read_qmmm_groups()
{
  const NcVar *qmmm_group_numbers_var = MD_trajectory_file->var(_QMMMGroup_);
  assert(qmmm_group_numbers_var->type() == ncInt);
  assert(qmmm_group_numbers_var->num_dims() == 1);
  assert(qmmm_group_numbers_var->get_dim(0)->size() == 7);
  
  int *qmmm_group_numbers = AllocateInt(7);
  assert(qmmm_group_numbers);
  
  assert(qmmm_group_numbers_var->get(qmmm_group_numbers, 7));

  assert(n_atoms == qmmm_group_numbers[0]);
  
  n_qm_atoms = qmmm_group_numbers[1];
  n_pseudo_atoms = qmmm_group_numbers[2];
  n_zero_charge_atoms =  qmmm_group_numbers[3];
  n_mm_atoms =  qmmm_group_numbers[4];

  DeAllocate(qmmm_group_numbers);

  cout << " There are " << n_qm_atoms << " QM atoms\n"
       << "           " << n_pseudo_atoms << " pseudo atoms\n" 
       << "           " << n_zero_charge_atoms << " zero charge atoms\n"
       << "           " << n_mm_atoms << " MM atoms" 
       << endl;
}

void ReadNetCDFTrajectoryFile::calculate_snapshot_size()
{
  snapshot_size = 3*n_use;
  
  if(save_esp_charges)
    snapshot_size += n_qm_atoms + n_pseudo_atoms;
  
  if(save_gradients)
    snapshot_size += 3*n_use;

  if(save_velocities)
    snapshot_size += 3*n_use;

  assert(MD_trajectory_file->dimension_exist(_SnapshotSize_) &&
	 MD_trajectory_file->dim(_SnapshotSize_)->size() == snapshot_size);
 
  cout << " NetCDF trajectory file snapshot size: " << snapshot_size << endl;
}

void ReadNetCDFTrajectoryFile::read_active_atom_list()
{
  const NcVar *active_atom_list_var = MD_trajectory_file->var(_ActiveAtomList_);
  assert(active_atom_list_var->type() == ncInt);
  assert(active_atom_list_var->num_dims() == 1);
  assert(active_atom_list_var->get_dim(0)->size() == n_use);
  
  DeAllocate(active_atom_list);
  active_atom_list = AllocateInt(n_use);
  assert(active_atom_list);
  
  MD_trajectory_file->get_int_var(_ActiveAtomList_, active_atom_list);
  
  cout << " Active atom list have been read" << endl;
}

void ReadNetCDFTrajectoryFile::read_fixed_atom_coordinates()
{
  const int n_fixed_atoms = n_atoms - n_use;

  const NcVar *xyz_fixed_var = MD_trajectory_file->var(_FixedAtomCoordinates_);
  assert(xyz_fixed_var->type() == ncDouble);
  assert(xyz_fixed_var->num_dims() == 1);
  assert(xyz_fixed_var->get_dim(0)->size() == 3*n_fixed_atoms);
  
  double *fixed_atom_coordinates = AllocateDouble(3*n_fixed_atoms);
  assert(fixed_atom_coordinates);
  assert(xyz_fixed_var->get(fixed_atom_coordinates, 3*n_fixed_atoms));
  
  cout << " Fixed atom coordinates have been read" << endl;
  
  DeAllocate(use);
  use = AllocateInt(n_atoms);
  assert(use);
  memset(use, 0, n_atoms*sizeof(int));
  for(int i = 0; i < n_use; i++)
    use[active_atom_list[i] - 1] = 1;

  DeAllocate(atom_coordinates);
  atom_coordinates = AllocateDouble(3*n_atoms);
  assert(atom_coordinates);
  memset(atom_coordinates, 0, 3*n_atoms*sizeof(double));
  int i3 = 0;
  for(int i = 0; i < n_atoms; i++) {
    if(!use[i]) {
      atom_coordinates[3*i] = fixed_atom_coordinates[i3++];
      atom_coordinates[3*i+1] = fixed_atom_coordinates[i3++];
      atom_coordinates[3*i+2] = fixed_atom_coordinates[i3++];
    }
  }
  assert(i3 == 3*n_fixed_atoms);

  DeAllocate(fixed_atom_coordinates);
}

void ReadNetCDFTrajectoryFile::read_atom_names()
{
  const NcVar *atom_names_var = MD_trajectory_file->var(_AtomNames_);
  assert(atom_names_var->type() == ncChar);
  assert(atom_names_var->num_dims() == 1);
  assert(atom_names_var->get_dim(0)->size() == 3*n_atoms);

  DeAllocate(atom_names);
  atom_names = AllocateChar(3*n_atoms);
  assert(atom_names);

  assert(atom_names_var->get(atom_names, 3*n_atoms));
  
  cout << " Atom names have been read" << endl;
}

void ReadNetCDFTrajectoryFile::read_atom_types()
{
  const NcVar *atom_types_var = MD_trajectory_file->var(_AtomTypes_);
  assert(atom_types_var->type() == ncInt);
  assert(atom_types_var->num_dims() == 1);
  assert(atom_types_var->get_dim(0)->size() == n_atoms);
  
  DeAllocate(atom_types);
  atom_types = AllocateInt(n_atoms);
  assert(atom_types);
  
  assert(atom_types_var->get(atom_types, n_atoms));

  cout << " Atom types have been read" << endl;
}

void ReadNetCDFTrajectoryFile::check_snapshot_information()
{
  const NcVar *snapshot_var = MD_trajectory_file->var(_Snapshot_);
  assert(snapshot_var->type() == ncDouble);
  assert(snapshot_var->num_dims() == 2);
  assert(snapshot_var->get_dim(0)->is_unlimited());
  assert(snapshot_var->get_dim(1)->size() == snapshot_size);

  n_snapshots = MD_trajectory_file->dim(_SnapshotUnlimit_)->size();
  
  cout << " There are " << n_snapshots << " snapshots have been saved" << endl;
}

void ReadNetCDFTrajectoryFile::read_and_write_all_snapshots()
{
  const int n_record = MD_trajectory_file->dim(_SnapshotUnlimit_)->size();
  NcVar *snaptshot_var = MD_trajectory_file->var(_Snapshot_);
  
  ofstream fout("md.xyz", ios::out);
  assert(fout.is_open());
  
  for(int i = 0; i < n_record; i++) {
    const NcValues *record = snaptshot_var->get_rec(i);
    const double *s = (double *) record->base();
    create_full_atom_coordinates(s);
    write_full_snapshot(fout, s);
  }
  
  fout.close();
}

void ReadNetCDFTrajectoryFile::read_snapshot(const int irecord)
{
  const int n_record = MD_trajectory_file->dim(_SnapshotUnlimit_)->size();
  assert(irecord < n_record);

  NcVar *snaptshot_var = MD_trajectory_file->var(_Snapshot_);
  const NcValues *record = snaptshot_var->get_rec(irecord);
  const double *s = (double *) record->base();
  create_full_atom_coordinates(s);

  s += 3*n_use;
  
  qm_esp_charges = 0;
  if(save_esp_charges) {
    qm_esp_charges = const_cast<double *> (s);
    s += n_qm_atoms + n_pseudo_atoms;
  }
  
  gradients = 0;
  if(save_gradients) {
    gradients = const_cast<double *> (s);
    s += 3*n_use;
  }

  velocities = 0;
  if(save_velocities) {
    velocities = const_cast<double *> (s);
    s += 3*n_use;
  }
}

void ReadNetCDFTrajectoryFile::create_full_atom_coordinates(const double *snapshot)
{
  int i3 = 0;
  for(int i = 0; i < n_use; i++) {
    const int i_atom = active_atom_list[i] - 1;
    atom_coordinates[3*i_atom] = snapshot[i3++];
    atom_coordinates[3*i_atom+1] = snapshot[i3++];
    atom_coordinates[3*i_atom+2] = snapshot[i3++];
  }
  assert(i3 == 3*n_use);
}

void ReadNetCDFTrajectoryFile::write_full_snapshot(ostream &fout, const double *snapshot)
{
  fout << "  " << n_atoms 
       << "  " << n_qm_atoms 
       << "  " << n_pseudo_atoms 
       << "  " << n_zero_charge_atoms 
       << "  " << n_mm_atoms << endl;

  char *atom_name = AllocateChar(3);
  assert(atom_name);
  
  double *s = const_cast<double *>(snapshot + 3*n_use);

  //double *qm_esp_charges = 0;
  qm_esp_charges = 0;
  if(save_esp_charges) {
    qm_esp_charges = s;
    s += n_qm_atoms + n_pseudo_atoms;
  }
  
  //double *gradients = 0;
  gradients = 0;
  if(save_gradients) gradients = s;

  int iuse = 0;
  
  for(int iatom = 0; iatom < n_atoms; iatom++) {
    if(iatom < n_qm_atoms + n_pseudo_atoms + n_zero_charge_atoms)
      assert(use[iatom]);
    
    memcpy(atom_name, atom_names+3*iatom, 3*sizeof(char));
    add_char_end(atom_name, 3);
    Cartesian r(atom_coordinates[3*iatom], atom_coordinates[3*iatom+1], 
		atom_coordinates[3*iatom+2]);

    int itype = -10000;
    if(iatom < n_qm_atoms)
      itype = _QMAtom_;
    else if(n_qm_atoms <= iatom && iatom < n_qm_atoms + n_pseudo_atoms)
      itype = _PseudoAtom_;
    else if(n_qm_atoms + n_pseudo_atoms <= iatom &&
	    iatom < n_qm_atoms + n_pseudo_atoms + n_zero_charge_atoms)
      itype = _ZeroChargeAtom_;
    else {
      if(use[iatom]) 
	itype = _FreeMMAtom_;
      else
	itype = _FixedMMAtom_;
    }
    
    NetCDFTrajectoryAtom a(iatom+1, itype, Str(atom_name), r);
    
    if(itype == _QMAtom_ || itype == _PseudoAtom_)
      if(qm_esp_charges)
	a.set_esp_charge(qm_esp_charges[iuse]);
    
    if(gradients && use[iatom]) a.set_gradient(Cartesian(gradients[3*iuse],
							 gradients[3*iuse+1],
							 gradients[3*iuse+2]));
    
    if(use[iatom] && (qm_esp_charges || gradients)) iuse++;
    
    a.write(fout);
  }
  
  DeAllocate(atom_name);
}

void ReadNetCDFTrajectoryFile::initialize()
{
  check_version();

  check_dimensions_and_variables();
  read_qmmm_groups();  
  calculate_snapshot_size();
  check_snapshot_information();

  read_atom_names();
  read_atom_types();
  read_active_atom_list();
  read_fixed_atom_coordinates();
}

void ReadNetCDFTrajectoryFile::check_version() const
{
  if(MD_trajectory_file->variable_exist(_NCMDTrajectoryVersion_)) {
    int version_number = 0;
    MD_trajectory_file->get_int_var(_NCMDTrajectoryVersion_, &version_number);
    assert(version_number >= _NCMDTrajectoryVersionNumber_1_);
  } else {
    cout << " There is no _NCMDTrajectoryVersionNumber_, which should be " 
	 << _NCMDTrajectoryVersionNumber_ << "\n"
	 << " It seems the MD trajectory file: " << MD_trajectory_file->NetCDF_file_name() 
	 << " is not a valid trajectory file" << endl;
    QCrash(" NetCDFTrajectoryFile::_check_version_: error");
  }
}
