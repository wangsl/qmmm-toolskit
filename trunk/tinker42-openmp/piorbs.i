 
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
c     ##  piorbs.i  --  conjugated system in the current structure  ##
c     ##                                                            ##
c     ################################################################
c
c
c     norbit    total number of pisystem orbitals in the system
c     iorbit    numbers of the atoms containing pisystem orbitals
c     reorbit   number of evaluations between orbital updates
c     piperp    atoms defining a normal plane to each orbital
c     nbpi      total number of bonds affected by the pisystem
c     bpi       bond and piatom numbers for each pisystem bond
c     ntpi      total number of torsions affected by the pisystem
c     tpi       torsion and pibond numbers for each pisystem torsion
c     listpi    atom list indicating whether each atom has an orbital
c
c
      integer norbit,iorbit
      integer reorbit,piperp
      integer nbpi,bpi
      integer ntpi,tpi
      logical listpi
      common /piorbs/ norbit,iorbit(maxpi),reorbit,piperp(3,maxpi),
     &                nbpi,bpi(3,maxpib),ntpi,tpi(2,maxpit),
     &                listpi(maxatm)
