 
C     $Id: kpolr.i 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  kpolr.i  --  forcefield parameters for polarizability  ##
c     ##                                                         ##
c     #############################################################
c
c
c     polr   dipole polarizability parameters for each atom type
c     pgrp   connected types in polarization group of each atom type
c
c
      integer pgrp
      real*8 polr
      common /kpolr/ polr(maxtyp),pgrp(maxval,maxtyp)
