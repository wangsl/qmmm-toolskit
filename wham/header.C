
#include <iostream>
using namespace std;
#include <ctime>
#include <cassert>
#include <unistd.h>

#define FORT(x) x##_

extern "C" {
  void FORT(write_header)();
}

void FORT(write_header)()
{
  /* Get the date/time and hostname */
  time_t t = time(0);
  char *hostname = new char [256];
  assert(hostname);
  assert(!gethostname(hostname, 256));
  char *pwd = new char [1024];
  assert(pwd);
  assert(getcwd(pwd, 1024));
  char *user_name = new char [1024];
  assert(user_name);
  char *user_name_tmp = getenv("USER");
  if(user_name_tmp)
    strcpy(user_name, user_name_tmp);
  else
    strcpy(user_name, "Unknown");
  cout << "****** WHAM Runing Information from Fotran *****************\n"
       << "* Hostname: " << hostname << "\n"
       << "* User: " << user_name << "\n"
       << "* Current directory: " << pwd << "\n"
       << "* Job begin time: " << ctime(&t)
       << "* Process id: " << getpid() << "\n"
       << "****** End of WHAM Running Information from Fortran ********\n"
       << flush;
  delete [] hostname;
  hostname = 0;
  delete [] pwd;
  pwd = 0;
  delete [] user_name;
  user_name = 0;
  user_name_tmp = 0;
}
