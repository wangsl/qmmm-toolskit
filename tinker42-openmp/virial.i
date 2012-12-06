 
C     $Id: virial.i 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  virial.i  --  components of internal virial tensor  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     vir    total internal virial Cartesian tensor components
c
c
      real*8 vir
      common /virial/ vir(3,3)
