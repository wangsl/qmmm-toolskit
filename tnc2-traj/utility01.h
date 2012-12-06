
/* $Id: utility01.h 2 2007-10-20 14:59:10Z wangsl $ */

#ifndef UTILITY01_H
#define UTILITY01_H

#include "fort.h"

int file_exist(const char *file_name);
void unlink_file(const char *file_name);
void link_file(const char *fname_old, const char *fname_new);
void rename_file(const char *fname_old, const char *fname_new);
void copy_file(const char *in_file, const char *out_file);

void show_file_content(const char *in_file);

#ifdef __cplusplus
extern "C" {
#endif
  
  void FORT(unlink_file_wrapper)(const char *file_name, const int &fname_len);
  void FORT(link_file_wrapper)(const char *fname_old, const int &len_old,
			       const char *fname_new, const int &len_new);
  void FORT(file_exist_wrapper)(const char *file_name, const int &fname_len,
				int &status);
  void FORT(rename_file_wrapper)(const char *fname_old, const int &len_old,
				 const char *fname_new, const int &len_new);
  void FORT(copy_file_wrapper)(const char *fname_old, const int &len_old,
			       const char *fname_new, const int &len_new);
  
#ifdef __cplusplus
}
#endif

#endif /* UTILITY01_H */
