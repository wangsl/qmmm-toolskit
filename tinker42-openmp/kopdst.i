 
C     $Id: kopdst.i 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  kopdst.i  --  forcefield parameters for out-plane distance  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     maxnopb   maximum number of out-of-plane distance entries
c
c     copb      force constant parameters for out-of-plane distance
c     kaopb     string of atom classes for out-of-plane distance
c
c
      integer maxnopd
      parameter (maxnopd=200)
      real*8 copd
      character*16 kaopd
      common /kopdst/ copd(maxnopd),kaopd(maxnopd)
