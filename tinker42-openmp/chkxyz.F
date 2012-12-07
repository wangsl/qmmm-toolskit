 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine chkxyz  --  check for coincident coordinates  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "chkxyz" finds any pairs of atoms with identical Cartesian
c     coordinates, and prints a warning message
c
c
      subroutine chkxyz
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'iounit.i'
      integer i,j
      real*8 xi,yi,zi,r2
      logical header
c
c
c     loop over atom pairs testing for identical coordinates
c
      header = .true.
C$OMP Parallel Do
C$OMP& Default(Shared)
C$OMP& Private(I, Xi, Yi, Zi, J, R2)
      do i = 1, n-1
         xi = x(i)
         yi = y(i)
         zi = z(i)
         do j = i+1, n
            r2 = (x(j)-xi)**2 + (y(j)-yi)**2 + (z(j)-zi)**2
            if (r2 .le. 0.0d0) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format ()
               end if
               write (iout,20)  i,j
   20          format (' CHKXYZ  --  Warning, Atoms',i6,' and',
     &                    i6,' have Identical Coordinates')
            end if
         end do
      end do
C$OMP End Parallel Do
      return
      end
