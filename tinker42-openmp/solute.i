 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  solute.i  --  parameters for continuum solvation models  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     rsolv    atomic radius of each atom for continuum solvation
c     vsolv    atomic volume of each atom for continuum solvation
c     asolv    atomic solvation parameters (kcal/mole/Ang**2)
c     rborn    Born radius of each atom for GB/SA solvation
c     drb      solvation derivatives with respect to Born radii
c     doffset  dielectric offset to continuum solvation atomic radii
c     p1       single-atom scale factor for analytical Still GB/SA
c     p2       1-2 interaction scale factor for analytical Still GB/SA
c     p3       1-3 interaction scale factor for analytical Still GB/SA
c     p4       nonbonded scale factor for analytical Still GB/SA
c     p5       soft cutoff parameter for analytical Still GB/SA
c     gpol     polarization self-energy values for each atom
c     shct     overlap scaling factors for Hawkins-Cramer-Truhlar GB/SA
c     wace     "omega" values for atom class pairs for use with ACE
c     s2ace    "sigma^2" values for atom class pairs for use with ACE
c     uace     "mu" values for atom class pairs for use with ACE
c     solvtyp  solvation model (ASP, SASA, ONION, STILL, HCT, ACE)
c
c
      real*8 rsolv,vsolv,asolv
      real*8 rborn,drb,doffset
      real*8 p1,p2,p3,p4,p5,gpol
      real*8 shct,wace,s2ace,uace
      character*5 solvtyp
      common /solute/ rsolv(maxatm),vsolv(maxatm),asolv(maxatm),
     &                rborn(maxatm),drb(maxatm),doffset,p1,p2,p3,p4,p5,
     &                gpol(maxatm),shct(maxatm),wace(maxclass,maxclass),
     &                s2ace(maxclass,maxclass),uace(maxclass,maxclass),
     &                solvtyp
