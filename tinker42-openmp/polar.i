 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  polar.i  --  polarizabilities and induced dipole moments  ##
c     ##                                                            ##
c     ################################################################
c
c
c     polarity  dipole polarizability for each multipole site (Ang**3)
c     pdamp     value of polarizability damping factor for each site
c     uind      induced dipole components at each multipole site
c     uinp      induced dipoles in field used for energy interactions
c     npolar    total number of polarizable sites in the system
c
c
      integer npolar
      real*8 polarity,pdamp
      real*8 uind,uinp
      common /polar/ polarity(maxatm),pdamp(maxatm),uind(3,maxatm),
     &               uinp(3,maxatm),npolar
