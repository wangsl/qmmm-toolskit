 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  angle.i  --  bond angles within the current structure  ##
c     ##                                                         ##
c     #############################################################
c
c
c     ak       harmonic angle force constant (kcal/mole/rad**2)
c     anat     ideal bond angle or phase shift angle (degrees)
c     afld     periodicity for Fourier bond angle term
c     nangle   total number of bond angles in the system
c     iang     numbers of the atoms in each bond angle
c     angtyp   potential energy function type for each bond angle
c
c
      integer nangle,iang
      real*8 ak,anat,afld
      character*8 angtyp
      common /angle/ ak(maxang),anat(maxang),afld(maxang),nangle,
     &               iang(4,maxang),angtyp(maxang)
