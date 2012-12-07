 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program testlight  --  time different neighbor schemes  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "testlight" performs a set of timing tests to compare the
c     evaluation of potential energy and energy/gradient using
c     the method of lights with a double loop over all atom pairs
c
c
      program testlight
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'cutoff.i'
      include 'deriv.i'
      include 'energi.i'
      include 'iounit.i'
      include 'potent.i'
      include 'vdwpot.i'
      integer i,j,k
      integer ncalls,npair
      integer map(maxlight)
      real*8 xi,yi,zi,r2
      real*8 elapsed,eps
      real*8 cev(3,maxatm)
      real*8 cec(3,maxatm)
      real*8 xsort(maxlight)
      real*8 ysort(maxlight)
      real*8 zsort(maxlight)
      logical exist,query
      logical header,match
      character*1 axis(3)
      character*120 string
      data axis  / 'X','Y','Z' /
c
c
c     read the molecular system and setup molecular mechanics
c
      call initial
      call getxyz
      call mechanic
c
c     get the number calculation cycles to be performed
c
      ncalls = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  ncalls
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Enter Desired Number of Repetitions [1] :  ',$)
         read (input,30)  ncalls
   30    format (i10)
      end if
      if (ncalls .eq. 0)  ncalls = 1
c
c     get the timing for setup of method of lights
c
      call setime
      do k = 1, ncalls
         do i = 1, n
            xsort(i) = x(i)
            ysort(i) = y(i)
            zsort(i) = z(i)
         end do
         call lights (n,map,xsort,ysort,zsort)
      end do
      call getime (elapsed)
      write (iout,40)  ncalls
   40 format (/,' Computation Overhead :',11x,'Time',
     &           7x,i5,' Evaluations')
      write (iout,50)  elapsed
   50 format (/,' Method of Lights',11x,f10.3)
c
c     get the timing for setup of double nested loop
c
      call setime
      npair = 0
      do k = 1, ncalls
         do i = 1, n-1
            xi = x(i)
            yi = y(i)
            zi = z(i)
            do j = i+1, n
               r2 = (x(j)-xi)**2 + (y(j)-yi)**2 + (z(j)-zi)**2
               if (r2 .lt. 100.0d0)  npair = npair + 1
            end do
         end do
      end do
      call getime (elapsed)
      write (iout,60)  elapsed
   60 format (' Double Nest Loop',11x,f10.3)
c
c     get the timing for energy terms via method of lights
c
      use_lights = .true.
      call setime
      do k = 1, ncalls
         if (use_vdw) then
            if (vdwtyp .eq. 'LENNARD-JONES')  call elj
            if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck
            if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb
            if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal
            if (vdwtyp .eq. 'GAUSSIAN')  call egauss
         end if
         if (use_charge)  call echarge
      end do
      call getime (elapsed)
      write (iout,70)
   70 format (/,' Potential Energy Only :',10x,'Time',14x,'EV',10x,'EC')
      write (iout,80)  elapsed,ev,ec
   80 format (/,' Method of Lights',11x,f10.3,6x,2f12.4)
c
c     get the timing for energy terms via double nested loop
c
      use_lights = .false.
      call setime
      do k = 1, ncalls
         if (use_vdw) then
            if (vdwtyp .eq. 'LENNARD-JONES')  call elj
            if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck
            if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb
            if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal
            if (vdwtyp .eq. 'GAUSSIAN')  call egauss
         end if
         if (use_charge)  call echarge
      end do
      call getime (elapsed)
      write (iout,90)  elapsed,ev,ec
   90 format (' Double Nest Loop',11x,f10.3,6x,2f12.4)
c
c     get the timing for gradient via method of lights
c
      use_lights = .true.
      call setime
      do k = 1, ncalls
         if (use_vdw) then
            if (vdwtyp .eq. 'LENNARD-JONES')  call elj1
            if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck1
            if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb1
            if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal1
            if (vdwtyp .eq. 'GAUSSIAN')  call egauss1
         end if
         if (use_charge)  call echarge1
      end do
      call getime (elapsed)
      write (iout,100)
  100 format (/,' Energy and Gradient :',12x,'Time')
      write (iout,110)  elapsed
  110 format (/,' Method of Lights',11x,f10.3)
c
c     store the method of lights gradient
c
      do i = 1, n
         do j = 1, 3
            cev(j,i) = dev(j,i)
            cec(j,i) = dec(j,i)
         end do
      end do
c
c     get the timing for gradient via double nested loop
c
      use_lights = .false.
      call setime
      do k = 1, ncalls
         if (use_vdw) then
            if (vdwtyp .eq. 'LENNARD-JONES')  call elj1
            if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck1
            if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb1
            if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal1
            if (vdwtyp .eq. 'GAUSSIAN')  call egauss1
         end if
         if (use_charge)  call echarge1
      end do
      call getime (elapsed)
      write (iout,120)  elapsed
  120 format (' Double Nest Loop',11x,f10.3)
c
c     compare method of lights and double loop vdw gradient
c
      eps = 1.0d-4
      match = .true.
      header = .true.
      if (use_vdw) then
         do i = 1, n
            do j = 1, 3
               if (abs(cev(j,i)-dev(j,i)) .gt. eps) then
                  if (header) then
                     match = .false.
                     header = .false.
                     write (iout,130)
  130                format (/,' Comparison of Lights and Double Loop',
     &                          ' vdw Gradients :',/)
                  end if
                  write (iout,140)  i,axis(j),cev(j,i),dev(j,i)
  140             format (10x,i6,' (',a1,')',4x,2f16.4)
               end if
            end do
         end do
      end if
c
c     compare method of lights and double loop charge gradient
c
      header = .true.
      if (use_charge) then
         do i = 1, n
            do j = 1, 3
               if (abs(cec(j,i)-dec(j,i)) .gt. eps) then
                  if (header) then
                     match = .false.
                     header = .false.
                     write (iout,150)
  150                format (/,' Comparison of Lights and Double Loop',
     &                          ' Charge Gradients :',/)
                  end if
                  write (iout,160)  i,axis(j),cec(j,i),dec(j,i)
  160             format (10x,i6,' (',a1,')',4x,2f16.4)
               end if
            end do
         end do
      end if
c
c     success if method of lights and double loop give same gradient
c
      if (match) then
         write (iout,170)
  170    format (/,' Method of Lights and Double Loop Gradients',
     &              ' are Identical')
      end if
c
c     perform any final tasks before program exit
c
      call final
      end
