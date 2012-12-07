 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine pressure  --  constant pressure via coupling  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "pressure" uses the internal virial to find the pressure
c     in a periodic box and maintains a constant desired pressure
c     by scaling the coordinates via coupling to an external
c     constant pressure bath
c
c     literature references:
c
c     H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren,
c     A. DiNola and J. R. Hauk, "Molecular Dynamics with Coupling
c     to an External Bath", Journal of Chemical Physics, 81,
c     3684-3690 (1984)
c
c     S. E. Feller, Y. Zhang, R. W. Pastor, B. R. Brooks, "Constant
c     Pressure Molecular Dynamics Simulation: The Langevin Piston
c     Method", Journal of Chemical Physics, 103, 4613-4621 (1995)
c
c
      subroutine pressure (dt,ekin,pres,stress)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'boxes.i'
      include 'bound.i'
      include 'group.i'
      include 'mdstuf.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,k
      integer start,stop
      real*8 dt,pres
      real*8 scale,third
      real*8 term,weigh
      real*8 xcm,xmove
      real*8 ycm,ymove
      real*8 zcm,zmove
      real*8 ekin(3,3)
      real*8 stress(3,3)
c
c
c     only necessary if periodic boundaries are in use
c
      if (.not. use_bounds)  return
c
c     calculate the pressor tensor for anisotropic systems
c
      term = prescon / volbox
      do i = 1, 3
         do j = 1, 3
            stress(j,i) = term * (2.0d0*ekin(j,i)-vir(j,i))
         end do
      end do
c
c     set isotropic pressure to the average of tensor diagonal
c
      pres = (stress(1,1)+stress(2,2)+stress(3,3)) / 3.0d0
c
c     find the scale factor to maintain constant pressure
c
      if (isobaric) then
         scale = 1.0d0
         third = 1.0d0 / 3.0d0
         if (barostat .eq. 'BERENDSEN') then
            scale = (1.0d0 + (dt*compress/taupres)*(pres-atmsph))**third
         end if
c
c     alter the size of the periodic box dimensions
c
         xbox = scale * xbox
         ybox = scale * ybox
         zbox = scale * zbox
c
c     update the lattice parameters for the periodic box
c
         call lattice
c
c     couple to pressure bath via atom scaling in Cartesian space
c
         if (integrate .ne. 'RIGIDBODY') then
            do i = 1, n
               if (use(i)) then
                  x(i) = scale * x(i)
                  y(i) = scale * y(i)
                  z(i) = scale * z(i)
               end if
            end do
c
c     couple to pressure bath via center of mass of rigid bodies
c
         else
            scale = scale - 1.0d0
            do i = 1, ngrp
               start = igrp(1,i)
               stop = igrp(2,i)
               xcm = 0.0d0
               ycm = 0.0d0
               zcm = 0.0d0
               do j = start, stop
                  k = kgrp(j)
                  weigh = mass(k)
                  xcm = xcm + x(k)*weigh
                  ycm = ycm + y(k)*weigh
                  zcm = zcm + z(k)*weigh
               end do
               xmove = scale * xcm/grpmass(i)
               ymove = scale * ycm/grpmass(i)
               zmove = scale * zcm/grpmass(i)
               do j = start, stop
                  k = kgrp(j)
                  x(k) = x(k) + xmove
                  y(k) = y(k) + ymove
                  z(k) = z(k) + zmove
               end do
            end do
         end if
      end if
      return
      end
