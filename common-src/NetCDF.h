
/* $Id$ */

#ifndef NET_CDF_H
#define NET_CDF_H

#include <netcdfcpp.h>
#include "table.h"

#define NCREPLACE 1
#define NCWRITE 2
#define NCREADONLY 3 

class NetCDFFile
{
  friend class NetCDFRestartFile;
public:
  NetCDFFile(const char *file_name_, int file_mode_ = NCREPLACE);
  ~NetCDFFile();

  void add_dim(const char *dim_name, long dim_size);
  
  void add_double_var(const char *var_name, const char *dim_name, 
		      const double *var_data = 0);
  void add_double_var(const char *var_name, double var_data);
  void add_double_var(const char *var_name);
  
  void add_int_var(const char *var_name, const char *dim_name, 
		   const int *var_data = 0);
  void add_int_var(const char *var_name, int var_data);
  void add_int_var(const char *var_name);

  void put_double_var(const char *var_name, const double val);
  void put_double_var(const char *var_name, const double *val);

  void put_int_var(const char *var_name, const int val);
  void put_int_var(const char *var_name, const int *val);
  
  double double_var(const char *var_name) const;
  void get_double_var(const char *var_name, double *val) const;
  
  int int_var(const char *var_name) const;
  void get_int_var(const char *var_name, int *val) const;
  
  int dimension_exist(const char *) const;
  int variable_exist(const char *) const;

  const NcDim *dim(const char *dim_name) const;
  const NcVar *var(const char *dim_name) const;

  NcVar *var(const char *dim_name);

  void add_unlimit_double_var(const char *var_name, const char *unlimit_dim_name,
			      const char *dim_name);
  
  void add_char_var(const char *var_name, const char *dim_name, 
		    const char *var_data = 0);
  void put_char_var(const char *var_name, const char *val);
  
  void flush() const;
  void close() const;
  int is_valid() const;
  const char *NetCDF_file_name() const;

  NcVar *variable_with_name(const char *var_name) const;

  int nc_file_mode() const { return file_mode; }

private :
  NcFile *nc_file; 
  char *file_name; // out
  int file_mode; // out
  
  Table<NcDim *> dimensions; // out
  Table<NcVar *> variables; // out

  void initialize();

  void update_dimensions_table_from_nc_file_dimensions();
  void update_variables_table_from_nc_file_variables();
  
  friend istream & operator>>(istream &s, NetCDFFile &c);
  friend ostream & operator<<(ostream &s, const NetCDFFile &c);
  void write_fields(ostream &s) const;
  int read_field(istream &s, const char *buf);
  void write_field_names() const;
};

#endif /* NET_CDF_H */
