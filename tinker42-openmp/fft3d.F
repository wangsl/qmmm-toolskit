 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine fftsetup  --  3-D FFT initialization  ##
c     ##                                                   ##
c     #######################################################
c
c
      subroutine fftsetup (n1,n2,n3,ntable,table)
      implicit none
      integer n1,n2,n3,ntable
      real*8 table(ntable,3)
c
c
      call cffti (n1,table(1,1))
      call cffti (n2,table(1,2))
      call cffti (n3,table(1,3))
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine fftfront  --  3-D FFT forward transform  ##
c     ##                                                      ##
c     ##########################################################
c
c
      subroutine fftfront (n1,n2,n3,ndim1,ndim2,ndim3,
     &                     ntable,nwork,table,work,w)
      implicit none
      integer i,j,k
      integer n1,n2,n3
      integer ndim1,ndim2,ndim3
      integer ntable,nwork
      real*8 table(ntable,3)
      double complex work(nwork)
      double complex w(ndim1,ndim2,ndim3)
c
c
c     forward transform along X, Y and then Z
c
      do k = 1, n3
         do j = 1, n2
            do i = 1, n1
               work(i) = w(i,j,k)
            end do
            call cfftf (n1,work,table(1,1))
            do i = 1, n1
               w(i,j,k) = work(i)
            end do
         end do
      end do
      do k = 1, n3
         do i = 1, n1
            do j = 1, n2
               work(j) = w(i,j,k)
            end do
            call cfftf (n2,work,table(1,2))
            do j = 1, n2
               w(i,j,k) = work(j)
            end do
         end do
      end do
      do i = 1, n1
         do j = 1, n2
            do k = 1, n3
               work(k) = w(i,j,k)
            end do
            call cfftf (n3,work,table(1,3))
            do k = 1, n3
               w(i,j,k) = work(k)
            end do
         end do
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine fftback  --  3-D FFT backward transform  ##
c     ##                                                      ##
c     ##########################################################
c
c
      subroutine fftback (n1,n2,n3,ndim1,ndim2,ndim3,
     &                    ntable,nwork,table,work,w)
      implicit none
      integer i,j,k
      integer n1,n2,n3
      integer ndim1,ndim2,ndim3
      integer ntable,nwork
      real*8 table(ntable,3)
      double complex work(nwork)
      double complex w(ndim1,ndim2,ndim3)
c
c
c     backward transform along X, Y and then Z
c
      do k = 1, n3
         do j = 1, n2
            do i = 1, n1
               work(i) = w(i,j,k)
            end do
            call cfftb (n1,work,table(1,1))
            do i = 1, n1
               w(i,j,k) = work(i)
            end do
         end do
      end do
      do k = 1, n3
         do i = 1, n1
            do j = 1, n2
               work(j) = w(i,j,k)
            end do
            call cfftb (n2,work,table(1,2))
            do j = 1, n2
               w(i,j,k) = work(j)
            end do
         end do
      end do
      do i = 1, n1
         do j = 1, n2
            do k = 1, n3
               work(k) = w(i,j,k)
            end do
            call cfftb (n3,work,table(1,3))
            do k = 1, n3
               w(i,j,k) = work(k)
            end do
         end do
      end do
      return
      end
