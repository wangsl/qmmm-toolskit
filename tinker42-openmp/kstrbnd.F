 
C     $Id: kstrbnd.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kstrbnd  --  assign stretch-bend parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kstrbnd" assigns the parameters for the stretch-bend
c     interactions and processes new or changed parameter values
c
c
      subroutine kstrbnd
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'iounit.i'
      include 'keys.i'
      include 'kstbnd.i'
      include 'potent.i'
      include 'strbnd.i'
      integer i,j,k,next
      integer ia,ib,ic
      integer it,itb,ina,inc
      integer nh,nb1,nb2
      real*8 sbk(3)
      logical header
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     process keywords containing stretch-bend parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'STRBND ') then
            it = 0
            do j = 1, 3
               sbk(j) = 0.0d0
            end do
            string = record(next:120)
            read (string,*,err=10,end=10)  it,(sbk(j),j=1,3)
   10       continue
            if (header) then
               header = .false.
               write (iout,20)
   20          format (/,' Additional Stretch-Bend Parameters :',
     &                 //,5x,'Atom Class',8x,'K(SB)-1',5x,'K(SB)-2',
     &                    5x,'K(SB)-3',/)
            end if
            write (iout,30)  it,(sbk(j),j=1,3)
   30       format (9x,i3,6x,3f12.3)
            do j = 1, 3
               stbn(j,it) = sbk(j)
            end do
         end if
      end do
c
c     assign stretch-bend parameters for each angle
c
      nstrbnd = 0
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         itb = class(ib)
         ina = atomic(ia)
         inc = atomic(ic)
         do k = 1, n12(ib)
            if (i12(k,ib) .eq. ia)  nb1 = bndlist(k,ib)
            if (i12(k,ib) .eq. ic)  nb2 = bndlist(k,ib)
         end do
         nh = 1
         if (ina .le. 1)  nh = nh + 1
         if (inc .le. 1)  nh = nh + 1
         if (stbn(nh,itb) .ne. 0.0d0) then
            nstrbnd = nstrbnd + 1
            isb(1,nstrbnd) = i
            isb(2,nstrbnd) = nb1
            isb(3,nstrbnd) = nb2
            ksb(nstrbnd) = stbn(nh,itb)
            if (mm2stbn) then
               if (ina .le. 1)  isb(2,nstrbnd) = 0
               if (inc .le. 1)  isb(3,nstrbnd) = 0
            end if
         end if
      end do
c
c     turn off the stretch-bend potential if it is not used
c
      if (nstrbnd .eq. 0)  use_strbnd = .false.
      return
      end
