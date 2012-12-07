
/* $Id$ */

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

static PyObject *_delete_read_tnc2_(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 

  delete_read_tnc2();

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

static PyObject *_atom_types_(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  const int *_atom_types_ = atom_types();
  assert(_atom_types_);
  const int natoms = n_atoms();
  long int dim1 [] = {natoms};
  PyArrayObject *types = (PyArrayObject *) PyArray_SimpleNewFromData(1, dim1, PyArray_INT, (void *) _atom_types_);
  return PyArray_Return(types);
}

static PyObject *_atom_names_(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  const char *_atom_names_ = atom_names();
  assert(_atom_names_);
  const int natoms = n_atoms();
  long int dim2 [] = {natoms, 3};
  PyArrayObject *names = (PyArrayObject *) PyArray_SimpleNewFromData(2, dim2, NPY_CHAR, (void *) _atom_names_);
  return PyArray_Return(names);
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
  long int dim2 [] = {natoms, 3};
  PyArrayObject *coords = (PyArrayObject *) PyArray_SimpleNewFromData(2, dim2, NPY_DOUBLE, (void *) xyz);
  return PyArray_Return(coords);
}

static PyObject *_n_active_atoms_(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  const int n = n_active_atoms();
  
  return Py_BuildValue("i", n);
}

static PyObject *_qm_esp_charges_saved_(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  return Py_BuildValue("i", qm_esp_charges_saved());
}

static PyObject *_gradients_saved_(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  return Py_BuildValue("i", gradients_saved());
}

static PyObject *_velocities_saved_(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  return Py_BuildValue("i", velocities_saved());
}

static PyObject *_active_atom_list_(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  const int *use = active_atom_list();
  assert(use);
  
  const int n_use = n_active_atoms();
  long int dim1 [] = {n_use};
  PyArrayObject *_use_ = (PyArrayObject *) PyArray_SimpleNewFromData(1, dim1, NPY_INT, (void *) use);
  return PyArray_Return(_use_);
}

static PyObject *_qm_esp_charges_(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  const double *esp = qm_esp_charges();
  
  if(!esp) {
    Py_INCREF(Py_None); 
    return Py_None;
  }

  const int n_qm_atoms = n_total_qm_atoms();
  long int dim1 [] = {n_qm_atoms};
  PyArrayObject *_esp_ = (PyArrayObject *) PyArray_SimpleNewFromData(1, dim1, NPY_DOUBLE, (void *) esp);
  return PyArray_Return(_esp_);
}

static PyObject *_gradients_(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  const double *g = gradients();
  
  if(!g) {
    Py_INCREF(Py_None); 
    return Py_None;
  }

  const int n_use = n_active_atoms();
  long int dim2 [] = {n_use, 3};
  PyArrayObject *_g_ = (PyArrayObject *) PyArray_SimpleNewFromData(2, dim2, PyArray_DOUBLE, (void *) g);
  return PyArray_Return(_g_);
}

static PyObject *_velocities_(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0; 
  
  const double *v = velocities();
  
  if(!v) {
    Py_INCREF(Py_None); 
    return Py_None;
  }

  const int n_use = n_active_atoms();
  long int dim2 [] = {n_use, 3};
  PyArrayObject *_v_ = (PyArrayObject *) PyArray_SimpleNewFromData(2, dim2, PyArray_DOUBLE, (void *) v);
  return PyArray_Return(_v_);
}

static PyMethodDef ExtensionFunctions [] = {
  {"create_read_tnc2_wrapper", _create_read_tnc2_, METH_VARARGS, "create tnc2 object"},
  {"snapshot_size_wrapper", _snapshot_size_, METH_VARARGS, " "},
  {"n_atoms_wrapper", _n_atoms_, METH_VARARGS, " "},
  {"atom_types_wrapper", _atom_types_, METH_VARARGS, " "},
  {"atom_names_wrapper", _atom_names_, METH_VARARGS, " "},
  {"snapshot_wrapper", _snapshot_, METH_VARARGS, " "},
  {"n_active_atoms_wrapper", _n_active_atoms_, METH_VARARGS, " "},
  {"active_atom_list_wrapper", _active_atom_list_, METH_VARARGS, " "},
  {"qm_esp_charges_saved_wrapper", _qm_esp_charges_saved_, METH_VARARGS, " "},
  {"gradients_saved_wrapper", _gradients_saved_, METH_VARARGS, " "},
  {"velocities_saved_wrapper", _velocities_saved_, METH_VARARGS, " "},
  {"qm_esp_charges_wrapper", _qm_esp_charges_, METH_VARARGS, " "},
  {"gradients_wrapper", _gradients_, METH_VARARGS, " "},
  {"velocities_wrapper", _velocities_, METH_VARARGS, " "},
  {"delete_read_tnc2_wrapper", _delete_read_tnc2_, METH_VARARGS, " "},
  { 0, 0 }
};

void initQMMMMDtnc2()
{
  char cmnt [] = "Python extension for QMMM MD trajectory tnc2";
  Py_InitModule3("QMMMMDtnc2", ExtensionFunctions, cmnt);
  import_array();
}
  
