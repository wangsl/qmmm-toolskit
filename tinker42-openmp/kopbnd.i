 
C     $Id: kopbnd.i 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  kopbnd.i  --  forcefield parameters for out-of-plane bend  ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxnopb   maximum number of out-of-plane bending entries
c
c     copb      force constant parameters for out-of-plane bending
c     kaopb     string of atom classes for out-of-plane bending
c
c
      integer maxnopb
      parameter (maxnopb=200)
      real*8 copb
      character*8 kaopb
      common /kopbnd/ copb(maxnopb),kaopb(maxnopb)
