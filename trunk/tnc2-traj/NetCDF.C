

/* $Id: NetCDF.C 137 2008-12-22 02:35:06Z wangsl $ */

#include <cassert>
#include "NetCDF.h"

#ifdef QCHEM
#include "tinker.h"
#endif

#include "str.h"

NetCDFFile::NetCDFFile(const char *file_name_, int file_mode_) : 
  nc_file(0), file_mode(file_mode_)
{ 
  int len = strlen(file_name_);
#ifdef QCHEM
  file_name = QAllocChar(len+1);
#else
  file_name = new char [len+1];
#endif
  insist(file_name);
  strcpy(file_name, file_name_);
  file_name[len] = '\0';
  initialize();
}

void NetCDFFile::initialize()
{
  insist(file_name);
  
  insist(file_mode == NCREPLACE || file_mode == NCWRITE ||
	 file_mode == NCREADONLY);

  insist(!nc_file);
  
  if(file_mode == NCREPLACE)
    nc_file = new NcFile(file_name, NcFile::Replace);
  else if(file_mode == NCREADONLY)
    nc_file = new NcFile(file_name, NcFile::ReadOnly);
  else if(file_mode == NCWRITE)
    nc_file = new NcFile(file_name, NcFile::Write);
  else 
    QCrash("unknown FileMode for NetCDFFile");
  
  insist(nc_file->is_valid());
  
  dimensions.default_value = 0;
  variables.default_value = 0;
  
  if(file_mode == NCREADONLY || file_mode == NCWRITE) {
    update_dimensions_table_from_nc_file_dimensions();
    update_variables_table_from_nc_file_variables();
  }
}

NetCDFFile::~NetCDFFile()
{
  dimensions.remove_all_entries();
  variables.remove_all_entries();
  delete nc_file;
  nc_file = 0;
  if(file_name) {
#ifdef QCHEM
    QFree(file_name);
#else    
    delete [] file_name;
#endif
    file_name = 0;
  }
}

void NetCDFFile::update_dimensions_table_from_nc_file_dimensions() 
{
  dimensions.remove_all_entries();
  for(int i = 0; i < nc_file->num_dims(); i++) {
    NcDim *nc_dim = nc_file->get_dim(i);
    insist(nc_dim->is_valid());
    dimensions[nc_dim->name()] = nc_dim;
  }
}

void NetCDFFile::update_variables_table_from_nc_file_variables()
{
  variables.remove_all_entries();
  for(int i = 0; i < nc_file->num_vars(); i++) {
    NcVar *nc_var = nc_file->get_var(i);
    insist(nc_var->is_valid());
    variables[nc_var->name()] = nc_var;
  }
}

int NetCDFFile::dimension_exist(const char *dim_name) const
{
  return dimensions.exists(dim_name);
}

int NetCDFFile::variable_exist(const char *var_name) const
{
  return variables.exists(var_name);
}

void NetCDFFile::add_dim(const char *dim_name, long dim_size)
{
  if(dimension_exist(dim_name)) {
    NcDim *dim = dimensions[dim_name];
    insist(dim->is_valid());
    if(dim->size() != dim_size) {
      cout << "NetCDFFile::add_dim dimension name '" << dim_name 
	   << "' exist with size " << dim->size() 
	   << " size added is " << dim_size << endl;
      QCrash("NetCDFFile::add_dim error");
    }
  } else {
    NcDim *new_dim = nc_file->add_dim(dim_name, dim_size);
    insist(new_dim->is_valid());
    dimensions[dim_name] = new_dim;
  }
}

void NetCDFFile::flush() const 
{ 
  insist(nc_file->sync());
}

const NcDim *NetCDFFile::dim(const char *dim_name) const
{
  if(!dimension_exist(dim_name)) {
    cout << "NetCDFFile::dim: dimension name '" << dim_name << "' does not exist" << endl;
    QCrash("NetCDFFile::dim error");
  }
  return dimensions[dim_name];
}

const NcVar *NetCDFFile::var(const char *var_name) const
{
  if(!variable_exist(var_name)) {
    cout << "NetCDFFile::var: variable name '" << var_name << "' does not exist" << endl;
    QCrash("NetCDFFile::var error");
  }
  return variables[var_name];
}

NcVar *NetCDFFile::var(const char *var_name)
{
  if(!variable_exist(var_name)) {
    cout << "NetCDFFile::var: variable name '" << var_name << "' does not exist" << endl;
    QCrash(" ");
  }
  return variables[var_name];
}

NcVar *NetCDFFile::variable_with_name(const char *var_name) const
{
  if(!variable_exist(var_name)) {
    cout << "NetCDFFile::var: variable name '" << var_name << "' does not exist" << endl;
    QCrash("NetCDFFile::var error");
  }
  return variables[var_name];
}

void NetCDFFile::add_double_var(const char *var_name, const char *dim_name, 
				const double *var_data)
{
  insist(dimension_exist(dim_name));
  const NcDim *dim = dimensions[dim_name];
  insist(dim->is_valid());
  
  insist(!variable_exist(var_name));
  NcVar *var = nc_file->add_var(var_name, ncDouble, dim);
  insist(var->is_valid());

  variables[var_name] = var;
  
  if(var_data) {
    const long &size = dim->size();
    insist(var->put(var_data, &size));
  }
}

void NetCDFFile::add_double_var(const char *var_name) 
{
  add_dim("one", 1);
  add_double_var(var_name, "one");
}

void NetCDFFile::add_double_var(const char *var_name, double var_data)
{
  add_dim("one", 1);
  add_double_var(var_name, "one", &var_data);
}

void NetCDFFile::add_int_var(const char *var_name, const char *dim_name, 
			     const int *var_data)
{
  insist(dimension_exist(dim_name));
  const NcDim *dim = dimensions[dim_name];
  insist(dim->is_valid());

  insist(!variable_exist(var_name));
  NcVar *var = nc_file->add_var(var_name, ncInt, dim);
  insist(var->is_valid());

  variables[var_name] = var;
  
  if(var_data) {
    const long &size = dim->size();
    insist(var->put(var_data, &size));
  }
}

void NetCDFFile::add_int_var(const char *var_name)
{
  add_dim("one", 1);
  add_int_var(var_name, "one");
}

void NetCDFFile::add_int_var(const char *var_name, int var_data)
{
  add_dim("one", 1);
  add_int_var(var_name, "one", &var_data);
}

double NetCDFFile::double_var(const char *var_name) const
{
  insist(variable_exist(var_name));
  const NcVar *var = nc_file->get_var(var_name);
  insist(var->is_valid());
  insist(var->num_dims() == 1);
  insist(var->get_dim(0)->size() == 1);
  double val = 0.0;
  insist(var->get(&val, 1));
  return val;
}

void NetCDFFile::get_double_var(const char *var_name, double *val) const
{
  insist(variable_exist(var_name));
  const NcVar *var = nc_file->get_var(var_name);
  insist(var->is_valid());
  insist(var->num_dims() == 1);
  long size = var->get_dim(0)->size();
  insist(size > 0);
  insist(var->get(val, size));
}

int NetCDFFile::int_var(const char *var_name) const
{
  insist(variable_exist(var_name));
  const NcVar *var = nc_file->get_var(var_name);
  insist(var->is_valid());
  insist(var->num_dims() == 1);
  insist(var->get_dim(0)->size() == 1);
  int val = 0;
  insist(var->get(&val, 1));
  return val;
}

void NetCDFFile::get_int_var(const char *var_name, int *val) const
{
  insist(variable_exist(var_name));
  const NcVar *var = nc_file->get_var(var_name);
  insist(var->is_valid());
  insist(var->num_dims() == 1);
  long size = var->get_dim(0)->size();
  insist(size > 0);
  insist(var->get(val, size));
}

void NetCDFFile::put_double_var(const char *var_name, const double val)
{
  insist(variable_exist(var_name));
  NcVar *var = nc_file->get_var(var_name);
  insist(var->is_valid());
  insist(var->num_dims() == 1);
  insist(var->get_dim(0)->size() == 1);
  const long size = 1;
  insist(var->put(&val, &size));
}

void NetCDFFile::put_double_var(const char *var_name, const double *val)
{
  insist(variable_exist(var_name));
  NcVar *var = nc_file->get_var(var_name);
  insist(var->is_valid());
  insist(var->num_dims() == 1);
  const long size = var->get_dim(0)->size();
  insist(var->put(val, &size));
}

void NetCDFFile::put_int_var(const char *var_name, const int val)
{
  insist(variable_exist(var_name));
  NcVar *var = nc_file->get_var(var_name);
  insist(var->is_valid());
  insist(var->num_dims() == 1);
  insist(var->get_dim(0)->size() == 1);
  const long size = 1;
  insist(var->put(&val, &size));
}

void NetCDFFile::put_int_var(const char *var_name, const int *val)
{
  insist(variable_exist(var_name));
  NcVar *var = nc_file->get_var(var_name);
  insist(var->is_valid());
  insist(var->num_dims() == 1);
  const long size = var->get_dim(0)->size();
  insist(var->put(val, &size));
}

void NetCDFFile::close() const
{
  insist(nc_file->close());
}

int NetCDFFile::is_valid() const
{ 
  return nc_file->is_valid() ? 1 : 0;
}

const char *NetCDFFile::NetCDF_file_name() const
{
  return file_name;
}

void NetCDFFile::add_unlimit_double_var(const char *var_name, const char *unlimit_dim_name,
					const char *dim_name)
{
  if(variable_exist(var_name)) 
    QCrash("NetCDF variable " + Str(var_name) + " has already existed");
  
  NcDim *unlimit_dim = 0;
  if(!dimension_exist(unlimit_dim_name)) {
    unlimit_dim = nc_file->add_dim(unlimit_dim_name);
    assert(unlimit_dim);
    dimensions[unlimit_dim_name] = unlimit_dim;
  } else 
    unlimit_dim = dimensions[unlimit_dim_name];
  
  assert(unlimit_dim->is_valid());
  assert(dimension_exist(dim_name) && dimensions[dim_name]->is_valid());
  
  NcVar *var = nc_file->add_var(var_name, ncDouble, unlimit_dim, dimensions[dim_name]);
  assert(var);
  
  variables[var_name] = var;
}

void NetCDFFile::add_char_var(const char *var_name, const char *dim_name, 
			      const char *var_data)
{
  assert(dimension_exist(dim_name));
  const NcDim *dim = dimensions[dim_name];
  assert(dim->is_valid());

  if(variable_exist(var_name))
    QCrash("NetCDF variable" + Str(var_name) + " has already existed");
  NcVar *var = nc_file->add_var(var_name, ncChar, dim);
  assert(var->is_valid());
  
  variables[var_name] = var;

  if(var_data) {
    const long &size = dim->size();
    insist(var->put(var_data, &size));
  }
}

void NetCDFFile::put_char_var(const char *var_name, const char *val)
{
  assert(variable_exist(var_name));
  NcVar *var = nc_file->get_var(var_name);
  assert(var == variables[var_name]);
  assert(var->is_valid());
  assert(var->num_dims() == 1);
  const long size = var->get_dim(0)->size();
  assert(var->put(val, &size));
}
