 
C     $Id: border.i 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  border.i  --  bond orders for a conjugated pisystem  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     pbpl   pi-bond orders for bonds in "planar" pisystem
c     pnpl   pi-bond orders for bonds in "nonplanar" pisystem
c
c
      real*8 pbpl,pnpl
      common /border/ pbpl(maxpib),pnpl(maxpib)
