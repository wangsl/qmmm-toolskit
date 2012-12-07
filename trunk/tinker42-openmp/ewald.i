 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  ewald.i  --  parameters for regular or PM Ewald summation  ##
c     ##                                                             ##
c     #################################################################
c
c
c     aewald    Ewald convergence coefficient value (Ang-1)
c     frecip    fractional cutoff value for reciprocal sphere
c     tinfoil   flag governing use of tinfoil boundary conditions
c
c
      real*8 aewald,frecip
      logical tinfoil
      common /ewald/ aewald,frecip,tinfoil
