 
C     $Id$
c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine rotpole  --  rotate multipoles to global frame  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "rotpole" constructs the set of atomic multipoles in the global
c     frame by applying the correct rotation matrix for each site
c
c
      subroutine rotpole
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      integer i
      real*8 a(3,3)
c
c
c     rotate the atomic multipoles at each site in turn
c
      do i = 1, npole
         call rotmat (i,a)
         call rotsite (i,a)
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rotmat  --  find global frame rotation matrix  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rotmat" finds the rotation matrix that converts from the local
c     coordinate system to the global frame at a multipole site
c
c
      subroutine rotmat (i,a)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      integer i
      real*8 r,dotxz
      real*8 dx,dy,dz
      real*8 dx1,dy1,dz1
      real*8 dx2,dy2,dz2
      real*8 a(3,3)
c
c
c     rotation matrix elements for z- and x-axes, first z then x
c
      if (polaxe(i) .eq. 'Z-then-X') then
         dx = x(zaxis(i)) - x(ipole(i))
         dy = y(zaxis(i)) - y(ipole(i))
         dz = z(zaxis(i)) - z(ipole(i))
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dx = x(xaxis(i)) - x(ipole(i))
         dy = y(xaxis(i)) - y(ipole(i))
         dz = z(xaxis(i)) - z(ipole(i))
         dotxz = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
         dx = dx - dotxz*a(1,3)
         dy = dy - dotxz*a(2,3)
         dz = dz - dotxz*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     rotation matrix elements for z- and x-axes, bisector method
c
      else if (polaxe(i) .eq. 'Bisector') then
         dx = x(zaxis(i)) - x(ipole(i))
         dy = y(zaxis(i)) - y(ipole(i))
         dz = z(zaxis(i)) - z(ipole(i))
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = x(xaxis(i)) - x(ipole(i))
         dy = y(xaxis(i)) - y(ipole(i))
         dz = z(xaxis(i)) - z(ipole(i))
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = dx1 + dx2
         dy = dy1 + dy2
         dz = dz1 + dz2
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dotxz = dx2*a(1,3) + dy2*a(2,3) + dz2*a(3,3)
         dx = dx2 - dotxz*a(1,3)
         dy = dy2 - dotxz*a(2,3)
         dz = dz2 - dotxz*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
      end if
c
c     finally, find rotation matrix elements for the y-axis
c
      a(1,2) = a(3,1)*a(2,3) - a(2,1)*a(3,3)
      a(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
      a(3,2) = a(2,1)*a(1,3) - a(1,1)*a(2,3)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rotsite  --  rotate multipoles at single site  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rotsite" computes the atomic multipoles at a specified site
c     in the global coordinate frame by applying a rotation matrix
c
c
      subroutine rotsite (imdq,a)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      integer i,j,k,m
      integer imdq
      real*8 a(3,3)
      real*8 m2(3,3)
      real*8 r2(3,3)
c
c
c     monopoles have the same value in any coordinate frame
c
      rpole(1,imdq) = pole(1,imdq)
c
c     rotate the dipoles to the global coordinate frame
c
      do i = 2, 4
         rpole(i,imdq) = 0.0d0
         do j = 2, 4
            rpole(i,imdq) = rpole(i,imdq) + pole(j,imdq)*a(i-1,j-1)
         end do
      end do
c
c     rotate the quadrupoles to the global coordinate frame
c
      k = 5
      do i = 1, 3
         do j = 1, 3
            m2(i,j) = pole(k,imdq)
            r2(i,j) = 0.0d0
            k = k + 1
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            if (j .lt. i) then
               r2(i,j) = r2(j,i)
            else
               do k = 1, 3
                  do m = 1, 3
                     r2(i,j) = r2(i,j) + a(i,k)*a(j,m)*m2(k,m)
                  end do
               end do
            end if
         end do
      end do
      k = 5
      do i = 1, 3
         do j = 1, 3
            rpole(k,imdq) = r2(i,j)
            k = k + 1
         end do
      end do
      return
      end
