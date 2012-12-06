 
C     $Id: korbs.i 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  korbs.i  --  forcefield parameters for pisystem orbitals  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxnpi     maximum number of pisystem bond parameter entries
c
c     electron   number of pi-electrons for each atom class
c     ionize     ionization potential for each atom class
c     repulse    repulsion integral value for each atom class
c     sslope     slope for bond stretch vs. pi-bond order
c     tslope     slope for 2-fold torsion vs. pi-bond order
c     kpi        string of atom classes for pisystem bonds
c
c
      integer maxnpi
      parameter (maxnpi=100)
      real*8 electron,ionize,repulse
      real*8 sslope,tslope
      character*8 kpi
      common /korbs/ electron(maxclass),ionize(maxclass),
     &               repulse(maxclass),sslope(maxnpi),tslope(maxnpi),
     &               kpi(maxnpi)
