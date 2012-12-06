 
C     $Id: gyrate.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine gyrate  --  compute the radius of gyration  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "gyrate" computes the radius of gyration of a molecular
c     system from its atomic coordinates
c
c
      subroutine gyrate (rg)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer i
      real*8 rg,xc,yc,zc
c
c
c     find the centroid of the atomic coordinates
c
      xc = 0.0d0
      yc = 0.0d0
      zc = 0.0d0
      do i = 1, n
         xc = xc + x(i)
         yc = yc + y(i)
         zc = zc + z(i)
      end do
      xc = xc / dble(n)
      yc = yc / dble(n)
      zc = zc / dble(n)
c
c     compute and print out the radius of gyration
c
      rg = 0.0d0
      do i = 1, n
         rg = rg + (x(i)-xc)**2 + (y(i)-yc)**2 + (z(i)-zc)**2
      end do
      rg = sqrt(rg/dble(n))
      return
      end
