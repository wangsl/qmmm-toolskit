 
C     $Id: ksttor.i 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  ksttor.i  --  forcefield parameters for stretch-torsions  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxnbt   maximum number of stretch-torsion parameter entries
c
c     btcon    force constant parameters for stretch-torsion
c     kbt      string of atom classes for bonds in stretch-torsion
c
c
      integer maxnbt
      parameter (maxnbt=200)
      real*8 btcon
      character*8 kbt
      common /ksttor/ btcon(3,maxnbt),kbt(maxnbt)
