 
C     $Id: vibrot.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program vibrot  --  vibrational analysis over torsions  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "vibrot" computes the eigenvalues and eigenvectors of the
c     torsional Hessian matrix
c
c
      program vibrot
      implicit none
      include 'sizes.i'
      include 'iounit.i'
      include 'omega.i'
      integer i,j,ihess
      real*8 hrot(maxrot,maxrot)
      real*8 matrix((maxrot+1)*maxrot/2)
      real*8 eigen(maxrot)
      real*8 vects(maxrot,maxrot)
      real*8 a(maxrot+1),b(maxrot+1)
      real*8 p(maxrot+1),w(maxrot+1)
      real*8 ta(maxrot+1),tb(maxrot+1)
      real*8 ty(maxrot+1)
c
c
c     compute torsional Hessian matrix elements
c
      call initial
      call getint
      call mechanic
      call initrot
      call hessrot ('FULL',hrot)
c
c     write out the torsional Hessian diagonal
c
      write (iout,10)
   10 format (/,' Diagonal of the Torsional Hessian :',/)
      write (iout,20)  (hrot(i,i),i=1,nomega)
   20 format (6f13.4)
c
c     write out the torsional Hessian elements
c
      write (iout,30)
   30 format (/,' Torsional Hessian Matrix Elements :',/)
      do i = 1, nomega
         write (iout,40)  (hrot(j,i),j=1,nomega)
   40    format (6f13.4)
      end do
c
c     place Hessian elements into triangular form
c
      ihess = 0
      do i = 1, nomega
         do j = i, nomega
            ihess = ihess + 1
            matrix(ihess) = hrot(i,j)
         end do
      end do
c
c     perform diagonalization to get Hessian eigenvalues
c
      call diagq (nomega,maxrot,nomega,matrix,eigen,vects,
     &                      a,b,p,w,ta,tb,ty)
      write (iout,50)
   50 format (/,' Eigenvalues of the Hessian Matrix :',/)
      write (iout,60)  (eigen(i),i=1,nomega)
   60 format (6d13.4)
c
c     perform any final tasks before program exit
c
      call final
      end
