 
/* $Id$ */

#include <iostream>
using namespace std;
#include <cassert>
#include "fint.h"
#include "tinkercpp.h"

void atomic_numbers_scan(const char *buf, INTEGER &n, INTEGER &qgrp, INTEGER &ygrp,
			 INTEGER &m0grp, INTEGER &mgrp, INTEGER &qmchg, INTEGER &spinmlt)
{
  assert(buf && strlen(buf) > 0);
  
  n = 0;
  qgrp = 0;
  ygrp = 0;
  m0grp = 0;
  mgrp = 0;
  
  qmchg = 0;
  spinmlt = 1;
  
#ifdef I64
  int m = sscanf(buf, "%ld%ld%ld%ld%ld%ld%ld", &n, &qgrp, &ygrp, &m0grp, &mgrp, &qmchg, &spinmlt);
#else
  int m = sscanf(buf, "%d%d%d%d%d%d%d", &n, &qgrp, &ygrp, &m0grp, &mgrp, &qmchg, &spinmlt);
#endif
  
  if(m != 1 && m != 5 && m != 7)
    QCrash("Error in reading atomic number information");
  
  if(m == 1) {
    qgrp = 0;
    ygrp = 0;
    m0grp = 0;
    mgrp = n;
  } 
  
  if(n != qgrp + ygrp + m0grp + mgrp)
    QCrash("atomic number information is not correct");
  
  assert(n > 0 && qgrp >= 0 && ygrp >= 0 && m0grp >= 0 && mgrp >= 0);
  
  if(!qgrp) 
    assert(!ygrp && !m0grp);
  
  return;
}

// Fortran version: AtomicNumbersScan
extern "C" {
  void FORT(atomicnumbersscan)(const char *buf, INTEGER &n, INTEGER &qgrp, INTEGER &ygrp,
                               INTEGER &m0grp, INTEGER &mgrp, INTEGER &qmchg, INTEGER &spinmlt)
  {
    atomic_numbers_scan(buf, n, qgrp, ygrp, m0grp, mgrp, qmchg, spinmlt);
  }
}
