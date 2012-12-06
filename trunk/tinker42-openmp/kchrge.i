 
C     $Id: kchrge.i 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  kchrge.i  --  forcefield parameters for partial charges  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     chg   partial charge parameters for each atom type
c
c
      real*8 chg
      common /kchrge/ chg(maxtyp)
