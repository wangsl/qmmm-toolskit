
#include <stdio.h>
#include <assert.h>
#include <Python.h>
#include <numpy/ndarrayobject.h>

#include "readtnc2.h"

/* 
   From "Guide to NumPy", Travis E. Oliphant
   12.1.1 PyArray_Type

typedef struct PyArrayObject {
  PyObject_HEAD
  char *data;
  int nd;
  npy_intp *dimensions;
  npy_intp *strides;
  PyObject *base;
  PyArray_Descr *descr;
  int flags;
  PyObject *weakreflist;
} PyArrayObject;
*/

#define DOUBLEDATA(p) ((double *)(((PyArrayObject*) p)->data))
#define FINTDATA(p) ((F77Int *)(((PyArrayObject *)p)->data))
#define INTDATA(p) ((int *)(((PyArrayObject *)p)->data))

#define NDIMS(p) ((int)(((PyArrayObject *)p)->nd))
#define DIMS(p) ((npy_intp *)(((PyArrayObject *)p)->dimensions))

void die(char *message) {
  fprintf(stderr, "%s\n", message);
  exit(13);
}

static PyObject *hello(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  printf("Hello world from file: %s, line: %d\n", __FILE__, __LINE__);
  
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject *_create_read_tnc2_(PyObject *self, PyObject *args)
{
  char *file_name = 0;
  if(!PyArg_ParseTuple(args, "s", &file_name))
    return 0;
  assert(file_name);

  create_read_tnc2(file_name);

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject *_snapshot_size_(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  const int n = snapshot_size();
  
  return Py_BuildValue("i", n);
}

static PyObject *_n_atoms_(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  const int n = n_atoms();
  
  return Py_BuildValue("i", n);
}


static PyObject *_snapshot_(PyObject *self, PyObject *args)
{
  int index = -1;
  if(!PyArg_ParseTuple(args, "i", &index))
    return 0;
  assert(index >= 0);
 
  const double *xyz = snapshot(index);
  assert(xyz);
 
  const int natoms = n_atoms();
  int dim2 [] = {natoms, 3};
  PyArrayObject *coords = (PyArrayObject *) PyArray_FromDimsAndData(2, dim2, PyArray_DOUBLE, 
								    (char *) (xyz));
  return PyArray_Return(coords);
}

#if 0

static PyObject *get_xyz(PyObject *self, PyObject *args)
{
  char *xyz_file_name = 0;
  if(!PyArg_ParseTuple(args, "s", &xyz_file_name))
    return 0;
  assert(xyz_file_name);

  const INTEGER len = strlen(xyz_file_name);
  
  FORT(getxyzwrapper)(xyz_file_name, &len);
    
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject *show_qm_atoms_as_xyz(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  FORT(showqmxyz)();
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject *coordinates_x(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  const int natoms = FORT(atoms).n_atoms;
  int dim1 [] = {natoms};
  PyArrayObject *x = (PyArrayObject *) PyArray_FromDimsAndData(1, dim1, PyArray_DOUBLE, 
							       (char *) (FORT(atoms).x));
  return PyArray_Return(x);
}

static PyObject *coordinates_y(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  const int natoms = FORT(atoms).n_atoms;
  int dim1 [] = {natoms};
  PyArrayObject *y = (PyArrayObject *) PyArray_FromDimsAndData(1, dim1, PyArray_DOUBLE, 
							       (char *) (FORT(atoms).y));
  return PyArray_Return(y);
}

static PyObject *coordinates_z(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  const int natoms = FORT(atoms).n_atoms;
  int dim1 [] = {natoms};
  PyArrayObject *z = (PyArrayObject *) PyArray_FromDimsAndData(1, dim1, PyArray_DOUBLE, 
							       (char *) (FORT(atoms).z));
  return PyArray_Return(z);
}

static PyObject *mechanic(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  FORT(mechanic)();
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject *create_NetCDF_trajectory_file_wrapper(PyObject *self, PyObject *args)
{
  char *nc_file_name = 0;
  int save_esp_charges = 0;
  if(!PyArg_ParseTuple(args, "si", &nc_file_name, &save_esp_charges))
    return 0; 
  
  FORT(netcdftrajectorysetup)();
  
  FORT(netcdftrajectorydata).NetCDF_trajectory_file_len = strlen(nc_file_name);
  FORT(netcdftrajectorydata).use_NetCDF_trajectory_file = 1;
  FORT(netcdftrajectorydata).save_esp_charges = save_esp_charges;
  memcpy(FORT(netcdftrajectorydata).NetCDF_trajectory_file, nc_file_name, 
	 strlen(nc_file_name)*sizeof(char));
    
  create_NetCDF_trajectory_file();

  if(nc_file_name) { nc_file_name = 0; }

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject *write_current_snapshot_to_NetCDF_trajectory_file_wrapper(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0;

  write_current_snapshot_to_NetCDF_trajectory_file();
  
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject *close_NetCDF_trajectory_file_wrapper(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0;

  close_NetCDF_trajectory_file();

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject *Tinker_test(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0;
  
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject *atom_numbers(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 

  int dim [] = {7};
  PyArrayObject *num_py = (PyArrayObject *) PyArray_FromDims(1, dim, PyArray_INT);
  assert(num_py);

  int *num = INTDATA(num_py);
  assert(num);

  num[0] = FORT(atoms).n_atoms;
  num[1] = FORT(qmmmgroup).n_qm_atoms;
  num[2] = FORT(qmmmgroup).n_pseudo_atoms;
  num[3] = FORT(qmmmgroup).n_zero_charge_atoms;
  num[4] = FORT(qmmmgroup).n_mm_atoms;
  num[5] = FORT(qmmmgroup).qm_total_charge;
  num[6] = FORT(qmmmgroup).qm_spin_multiplicity;
  
  return PyArray_Return(num_py);
}

static PyObject *open_NetCDF_trajectory_file_to_read_wrapper(PyObject *self, PyObject *args)
{
  char *nc_file_name = 0;
  if(!PyArg_ParseTuple(args, "s", &nc_file_name))
    return 0; 

  open_Tinker_NetCDF_trajector_file_to_read(nc_file_name);
  
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject *close_NetCDF_trajectory_file_2(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject *read_Tinker_NetCDF_timestep_wrapper(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  int OK = read_Tinker_NetCDF_timestep();

  return Py_BuildValue("i", OK);
}

static PyObject *NetCDF_trajectory_file_setup_wrapper(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  Py_INCREF(Py_None); 
  return Py_None;
}


static PyMethodDef ExtensionFunctions [] = {
  {"get_xyz", get_xyz, METH_VARARGS, "read Tinker xyz"},
  {"show_qm_atoms_as_xyz", show_qm_atoms_as_xyz, METH_VARARGS, "show qm atoms"},
  {"coordinates_x", coordinates_x, METH_VARARGS, "x coordinates"},
  {"coordinates_y", coordinates_y, METH_VARARGS, "y coordinates"},
  {"coordinates_z", coordinates_z, METH_VARARGS, "z coordinates"},
  {"mechanic", mechanic, METH_VARARGS, "mechanic in Tinker"},
  {"atom_numbers", atom_numbers, METH_VARARGS, "atom numbers"},  
  {"Tinker_test", Tinker_test, METH_VARARGS, "Tinker test"},
  {"create_NetCDF_trajectory_file", create_NetCDF_trajectory_file_wrapper, 
   METH_VARARGS, "create NetCDF trajectory file"},
  {"write_current_snapshot_to_NetCDF_trajectory_file", 
   write_current_snapshot_to_NetCDF_trajectory_file_wrapper,
   METH_VARARGS, "write current snapshots to NetCDF trajectory file"},
  {"close_NetCDF_trajectory_file", close_NetCDF_trajectory_file_wrapper,
   METH_VARARGS, "close NetCDF trajectory file"},
  {"open_NetCDF_trajectory_file_to_read", open_NetCDF_trajectory_file_to_read_wrapper,
   METH_VARARGS, "open a NetCDF trajectory file to read"},
  {"close_NetCDF_trajectory_file_2", close_NetCDF_trajectory_file_2,
   METH_VARARGS, "close NetCDF trajectory file"},
  {"read_Tinker_NetCDF_timestep", read_Tinker_NetCDF_timestep_wrapper,
   METH_VARARGS, "read NetCDF trajectory snapshot"},

  {0, 0}
};

#endif

static PyMethodDef ExtensionFunctions [] = {
  {"hello", hello, METH_VARARGS, "Hello world test"},
  {"create_read_tnc2_wrapper", _create_read_tnc2_, METH_VARARGS, "create tnc2 object"},
  {"snapshot_size_wrapper", _snapshot_size_, METH_VARARGS, " "},
  {"n_atoms_wrapper", _n_atoms_, METH_VARARGS, " "},
  {"snapshot_wrapper", _snapshot_,  METH_VARARGS, " "},
  { 0, 0 }
};


void initQMMMMDtnc2()
{
  char cmnt [] = "Python extension for QMMM MD trajectory tnc2";
  Py_InitModule3("QMMMMDtnc2", ExtensionFunctions, cmnt);
  import_array();
}
  
