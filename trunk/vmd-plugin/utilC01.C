
/* $Id: utilC01.C 161 2010-03-08 17:05:12Z wangsl $ */

#include <iostream>
#include <fstream>
using namespace std;
#include <unistd.h>
#include "fort.h"
#include "die.h"
#include "str.h"

#ifdef QCHEM
#include "tinker.h"
#endif

int file_exist(const char *file_name)
{
  return access(file_name, F_OK) ? 0 : 1;
}

void unlink_file(const char *file_name)
{
  if(file_exist(file_name))
    if(unlink(file_name))
      cerr << "Failed to remove file " << file_name << endl;
}

void link_file(const char *fname_old, const char *fname_new)
{
  if(!file_exist(fname_old))
    cerr << "Source file " << fname_old << " does not exist" << endl;
  
  if(file_exist(fname_new))
    cerr << "New file " << fname_new << " exist already" << endl;
  
  if(link(fname_old, fname_new)) 
    cerr << "Failed to link file " << fname_old << " to " << fname_new << endl;
}

void rename_file(const char *fname_old, const char *fname_new)
{
  if(!file_exist(fname_old))
    cerr << "Source file " << fname_old << " does not exist" << endl;
  
  if(file_exist(fname_new))
    cerr << "New file " << fname_new << " exist already" << endl;
  
  if(rename(fname_old, fname_new)) 
    cerr << "Failed to rename file " << fname_old << " to " << fname_new << endl;
}

void copy_file(const char *in_file, const char *out_file)
{
  if(!file_exist(in_file))
    cerr << "Input file " << in_file << " does not exist" << endl;

  if(file_exist(out_file))
    cerr << "Output file " << out_file << " exist already" << endl;
  
  char c;
  ifstream in(in_file);
  if(!in)
    cerr << "Failed to open " << in_file << endl;
  
  ofstream out(out_file);
  if(!out)
    cerr << "Failed to create " << out_file << endl;
  
  while(out && in.get(c))
    out.put(c);
  
  out.close();
  in.close();
}

void FORT(unlink_file_wrapper)(const char *file_name, const int &fname_len)
{  
  char *fname = new char [fname_len+2];
  insist(fname);
  memcpy(fname, file_name, fname_len*sizeof(char));
  fname[fname_len] = '\0';

  unlink_file(fname);
  
  if(fname) {
    delete [] fname;
    fname = 0;
  }
}

void FORT(link_file_wrapper)(const char *fname_old, const int &len_old,
			     const char *fname_new, const int &len_new)
{
  char *fold = new char [len_old+2];
  insist(fold);
  memcpy(fold, fname_old, len_old*sizeof(char));
  fold[len_old] = '\0';

  char *fnew = new char [len_new+2];
  insist(fnew);
  memcpy(fnew, fname_new, len_new*sizeof(char));
  fnew[len_new] = '\0';
  
  link_file(fold, fnew);

  if(fold) {
    delete [] fold;
    fold = 0;
  }

  if(fnew) {
    delete [] fnew;
    fnew = 0;
  }
}

void FORT(file_exist_wrapper)(const char *file_name, const int &fname_len,
			      int &status)
{  
  char *fname = new char [fname_len+2];
  insist(fname);
  memcpy(fname, file_name, fname_len*sizeof(char));
  fname[fname_len] = '\0';

  status = file_exist(fname);

  if(fname) {
    delete [] fname;
    fname = 0;
  }
}

void FORT(rename_file_wrapper)(const char *fname_old, const int &len_old,
			       const char *fname_new, const int &len_new)
{
  char *fold = new char [len_old+2];
  insist(fold);
  memcpy(fold, fname_old, len_old*sizeof(char));
  fold[len_old] = '\0';

  char *fnew = new char [len_new+2];
  insist(fnew);
  memcpy(fnew, fname_new, len_new*sizeof(char));
  fnew[len_new] = '\0';
  
  rename_file(fold, fnew);

  if(fold) {
    delete [] fold;
    fold = 0;
  }

  if(fnew) {
    delete [] fnew;
    fnew = 0;
  }
}

void FORT(copy_file_wrapper)(const char *fname_old, const int &len_old,
			     const char *fname_new, const int &len_new)
{
  char *fold = new char [len_old+2];
  insist(fold);
  memcpy(fold, fname_old, len_old*sizeof(char));
  fold[len_old] = '\0';

  char *fnew = new char [len_new+2];
  insist(fnew);
  memcpy(fnew, fname_new, len_new*sizeof(char));
  fnew[len_new] = '\0';
  
  copy_file(fold, fnew);

  if(fold) {
    delete [] fold;
    fold = 0;
  }

  if(fnew) {
    delete [] fnew;
    fnew = 0;
  }
}

void show_file_content(const char *in_file)
{
  if(!file_exist(in_file)) {
    cout << "* Input file " << in_file << " does not exist *" << endl;
    return;
  }
  
  char c;
  ifstream fin(in_file);
  if(!fin)
    cerr << "Failed to open " << in_file << endl;
  
  while(fin.get(c))
    cout.put(c);
  
  fin.close();
}

// Fortran version: ShowFileContent

extern "C" void FORT(showfilecontent)(const char *file_name, const int &len_file_name)
{
  show_file_content(Str(file_name, len_file_name));
}
