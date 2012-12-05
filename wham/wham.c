
#include <stdio.h>
#include <assert.h>
#include <Python.h>
#include <numpy/ndarrayobject.h>
#include "fint.h"
#include "fort.h"
#include "wham.h"
				      
/*****
 * For this new version NumPy is used
 * Guide to NumPy
 *****/

/*
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
#define FINTDATA(p) ((FInt *)(((PyArrayObject *)p)->data))
#define INTDATA(p) ((int *)(((PyArrayObject *)p)->data))

#define NDIMS(p) ((int)(((PyArrayObject *)p)->nd))
#define DIMS(p) ((npy_intp *)(((PyArrayObject *)p)->dimensions))

static PyObject *wham_iteration(PyObject *self, PyObject *args)
{
  PyObject *ebw_py = 0;
  PyObject *nt_py = 0;
  PyObject *kBT_py = 0;
  PyObject *ebf_py = 0;
  int maximum_iterations = -100;
  double tolerance = -1.0e+23;
  int Nwind = -100;
  int max_nt = -100;
  char *iteration_output = 0;
  if(!PyArg_ParseTuple(args, "OidOOOiis",
		       &ebw_py, &maximum_iterations, &tolerance,
		       &nt_py, &kBT_py, &ebf_py, &Nwind, &max_nt, 
		       &iteration_output))
    return 0;
  
  assert(iteration_output);
  
  double *ebw = DOUBLEDATA(ebw_py);
  assert(ebw);
  assert(NDIMS(ebw_py) == 3);
  npy_intp *ebw_dim = DIMS(ebw_py);
  assert(ebw_dim);
  assert(ebw_dim[0] == Nwind && ebw_dim[1] == max_nt && ebw_dim[2] == Nwind);
  
  // Be careful of integer arrays in Python
  FInt *nt = FINTDATA(nt_py);
  assert(nt);
  assert(NDIMS(nt_py) == 1);
  assert(DIMS(nt_py)[0] == Nwind);
  
  double *kBT = DOUBLEDATA(kBT_py);
  assert(kBT);
  assert(NDIMS(kBT_py) == 1);
  assert(DIMS(nt_py)[0] == Nwind);
  
  double *ebf = DOUBLEDATA(ebf_py);
  assert(ebf);
  assert(NDIMS(nt_py) == 1);
  assert(DIMS(nt_py)[0] == Nwind);
  
  double *ebf2 = (double *) malloc(Nwind*sizeof(double));
  assert(ebf2);

  int len_iteration_output = strlen(iteration_output);

  // For 64-byte computer, Integer 8 bytes
  FInt NWindL = Nwind;
  FInt MaxIterL = maximum_iterations;
  FInt MaxNTL = max_nt;
  FInt LenIterOutputL = len_iteration_output;
  
  Py_BEGIN_ALLOW_THREADS;
  FORT(wham_iteration_with_fortran)(ebw, &NWindL, &MaxIterL, &tolerance, 
				    nt, &MaxNTL, kBT, ebf, ebf2, 
				    iteration_output, &LenIterOutputL);
  Py_END_ALLOW_THREADS;
  
  PyObject *ebf2_lst = PyList_New(Nwind);
  assert(ebf2_lst);
  
  int i;
  for(i = 0; i < Nwind; i++)
    PyList_SET_ITEM(ebf2_lst, i, PyFloat_FromDouble(ebf2[i]));
  
  free(ebf2);
  ebf2 = 0;

  return ebf2_lst;
}

static PyObject *cpu_numbers(PyObject *self, PyObject *args)
{
  if(!PyArg_ParseTuple(args, ""))
    return 0;
  // number of active/running CPUs 
  int ncpu = sysconf(_SC_NPROCESSORS_ONLN);
  printf("There are %4i cpus detected.", ncpu);
  return Py_BuildValue("i", ncpu);
}

static PyMethodDef WHAM_Functions [] = {
  {"wham_iteration", wham_iteration, METH_VARARGS, 
   "WHAM iteration from c and Fortran with OpenMP"},
  {"cpu_numbers", cpu_numbers, METH_VARARGS, "total cpu numbers available"},
  {0}
};

void initwhamit()
{
  char cmmt [] = "WHAM iterations";
  Py_InitModule3("whamit", WHAM_Functions, cmmt);
  import_array();
}



