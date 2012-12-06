 
C     $Id: kopbend.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kopbend  --  out-of-plane bending parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kopbend" assigns the force constants for out-of-plane
c     bending at trigonal centers via Wilson-Decius-Cross angle
c     bends; also processes any new or changed parameter values
c
c
      subroutine kopbend
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'kopbnd.i'
      include 'opbend.i'
      include 'potent.i'
      integer i,j,nopb
      integer ia,ib,id,it
      integer itb,itd,number
      integer next,size
      real*8 fopb
      logical header
      logical jopb(maxclass)
      character*4 pa,pb,pd
      character*8 blank,ps
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     process keywords containing out-of-plane bend parameters
c
      blank = '        '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'OPBEND ') then
            ia = 0
            ib = 0
            fopb = 0.0d0
            string = record(next:120)
            read (string,*,err=10,end=10)  ia,ib,fopb
   10       continue
            if (header) then
               header = .false.
               write (iout,20)
   20          format (/,' Additional Out-of-Plane Bend Parameters :',
     &                 //,5x,'Atom Classes',7x,'K(OPB)',/)
            end if
            write (iout,30)  ia,ib,fopb
   30       format (6x,2i4,4x,f12.3)
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            ps = pa//pb
            do j = 1, maxnopb
               if (kaopb(j).eq.blank .or. kaopb(j).eq.ps) then
                  kaopb(j) = ps
                  copb(j) = fopb
                  goto 50
               end if
            end do
            write (iout,40)
   40       format (/,' KOPBEND --  Too many Out-of-Plane',
     &                 ' Angle Bending Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      nopb = maxnopb
      do i = maxnopb, 1, -1
         if (kaopb(i) .eq. blank)  nopb = i - 1
      end do
c
c     make list of atom classes using out-of-plane bending
c
      do i = 1, maxclass
         jopb(i) = .false.
      end do
      do i = 1, maxnopb
         if (kaopb(i) .eq. blank)  goto 60
         it = number (kaopb(i)(1:4))
         jopb(it) = .true.
      end do
   60 continue
c
c     assign out-of-plane bending parameters for each angle
c
      nopbend = 0
      if (nopb .ne. 0) then
         header = .true.
         do i = 1, nangle
            ib = iang(2,i)
            itb = class(ib)
            if (jopb(itb) .and. n12(ib).eq.3) then
               id = iang(4,i)
               itd = class(id)
               size = 4
               call numeral (itb,pb,size)
               call numeral (itd,pd,size)
               ps = pb//pd
               do j = 1, nopb
                  if (kaopb(j) .eq. ps) then
                     nopbend = nopbend + 1
                     iopb(nopbend) = i
                     kopb(nopbend) = copb(j)
                     goto 90
                  end if
               end do
               abort = .true.
               if (header) then
                  header = .false.
                  write (iout,70)
   70             format (/,' Undefined Out-of-Plane Bending',
     &                       ' Parameters :',
     &                    //,' Type',13x,'Atom Names',11x,
     &                       'Atom Classes',/)
               end if
               write (iout,80)  ib,name(ib),id,name(id),itb,itd
   80          format (' Angle-OP',3x,i6,'-',a3,i6,'-',a3,7x,2i5)
   90          continue
            else
               iang(4,i) = ib
            end if
         end do
      end if
c
c     mark angles at trigonal sites to use projected in-plane values
c
      do i = 1, nopbend
         j = iopb(i)
         if (angtyp(j) .eq. 'HARMONIC')  angtyp(j) = 'IN-PLANE'
      end do
c
c     turn off the out-of-plane bending term if it is not used
c
      if (nopbend .eq. 0)  use_opbend = .false.
      return
      end
