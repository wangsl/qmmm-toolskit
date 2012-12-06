 
C     $Id: chgpot.i 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  chgpot.i  --  specifics of charge-charge functional form  ##
c     ##                                                            ##
c     ################################################################
c
c
c     dielec    dielectric constant for electrostatic interactions
c     c2scale   factor by which 1-2 charge interactions are scaled
c     c3scale   factor by which 1-3 charge interactions are scaled
c     c4scale   factor by which 1-4 charge interactions are scaled
c     c5scale   factor by which 1-5 charge interactions are scaled
c     neutnbr   logical flag governing use of neutral group neighbors
c     neutcut   logical flag governing use of neutral group cutoffs
c
c
      real*8 dielec
      real*8 c2scale,c3scale
      real*8 c4scale,c5scale
      logical neutnbr,neutcut
      common /chgpot/ dielec,c2scale,c3scale,c4scale,c5scale,neutnbr,
     &                neutcut
