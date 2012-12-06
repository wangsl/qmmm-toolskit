 
C     $Id: kopdist.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine kopdist  --  out-of-plane distance parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "kopdist" assigns the force constants for out-of-plane
c     distance at trigonal centers via the central atom height;
c     also processes any new or changed parameter values
c
c
      subroutine kopdist
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'kopdst.i'
      include 'opdist.i'
      include 'potent.i'
      integer i,j,k,nopd
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer size,next
      real*8 fopd
      logical header
      character*4 pa,pb,pc,pd
      character*16 blank,pti
      character*16 pt(6)
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     process keywords containing out-of-plane distance parameters
c
      blank = '                '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'OPDIST ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            fopd = 0.0d0
            string = record(next:120)
            read (string,*,err=10,end=10)  ia,ib,ic,id,fopd
   10       continue
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            pti = pa//pb//pc//pd
            if (header) then
               header = .false.
               write (iout,20)
   20          format (/,' Additional Out-of-Plane Distance',
     &                    ' Parameters :',
     &                 //,5x,'Atom Classes',7x,'K(OPD)',/)
            end if
            write (iout,30)  ia,ib,ic,id,fopd
   30       format (1x,4i4,1x,2f12.3)
            do j = 1, maxnopd
               if (kaopd(j).eq.blank .or. kaopd(j).eq.pti) then
                  kaopd(j) = pti
                  copd(j) = fopd
                  goto 50
               end if
            end do
            write (iout,40)
   40       format (/,' KOPDIST  --  Too many Out-of-Plane Distance',
     &                 ' Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      nopd = maxnopd
      do i = maxnopd, 1, -1
         if (kaopd(i) .eq. blank)  nopd = i - 1
      end do
c
c     assign out-of-plane distance parameters for trigonal sites
c
      nopdist = 0
      if (nopd .ne. 0) then
         do i = 1, n
            if (n12(i) .eq. 3) then
               ia = i
               ib = i12(1,i)
               ic = i12(2,i)
               id = i12(3,i)
               ita = class(ia)
               itb = class(ib)
               itc = class(ic)
               itd = class(id)
               size = 4
               call numeral (ita,pa,size)
               call numeral (itb,pb,size)
               call numeral (itc,pc,size)
               call numeral (itd,pd,size)
               pt(1) = pa//pb//pc//pd
               pt(2) = pa//pb//pd//pc
               pt(3) = pa//pc//pb//pd
               pt(4) = pa//pc//pd//pb
               pt(5) = pa//pd//pb//pc
               pt(6) = pa//pd//pc//pb
               do j = 1, nopd
                  if (kaopd(j)(1:4) .eq. pa) then
                     do k = 1, 6
                        if (kaopd(j) .eq. pt(k)) then
                           nopdist = nopdist + 1
                           iopd(1,nopdist) = ia
                           if (k .eq. 1) then
                              iopd(2,nopdist) = ib
                              iopd(3,nopdist) = ic
                              iopd(4,nopdist) = id
                           else if (k .eq. 2) then
                              iopd(2,nopdist) = ib
                              iopd(3,nopdist) = id
                              iopd(4,nopdist) = ic
                           else if (k .eq. 3) then
                              iopd(2,nopdist) = ic
                              iopd(3,nopdist) = ib
                              iopd(4,nopdist) = id
                           else if (k .eq. 4) then
                              iopd(2,nopdist) = ic
                              iopd(3,nopdist) = id
                              iopd(4,nopdist) = ib
                           else if (k .eq. 5) then
                              iopd(2,nopdist) = id
                              iopd(3,nopdist) = ib
                              iopd(4,nopdist) = ic
                           else if (k .eq. 6) then
                              iopd(2,nopdist) = id
                              iopd(3,nopdist) = ic
                              iopd(4,nopdist) = ib
                           end if
                           kopd(nopdist) = copd(j)
                           goto 60
                        end if
                     end do
                  end if
               end do
   60          continue
            end if
         end do
      end if
c
c     mark angles at trigonal sites to use projected in-plane values
c
      do i = 1, nopdist
         ia = iopd(1,i)
         do j = 1, 3
            k = anglist(j,ia)
            if (angtyp(k) .eq. 'HARMONIC')  angtyp(k) = 'IN-PLANE'
         end do
      end do
c
c     turn off out-of-plane distance potential if it is not used
c
      if (nopdist .eq. 0)  use_opdist = .false.
      return
      end
