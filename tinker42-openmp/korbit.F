 
C     $Id: korbit.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine korbit  --  pisystem parameter assignment  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "korbit" assigns pi-orbital parameters to conjugated
c     systems and processes any new or changed parameters
c
c
      subroutine korbit
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'korbs.i'
      include 'orbits.i'
      include 'piorbs.i'
      include 'pistuf.i'
      include 'tors.i'
      include 'units.i'
      integer i,j,k,npi
      integer ia,ib,it
      integer ita,itb
      integer size,next
      real*8 elect,ioniz
      real*8 repuls
      real*8 sslop,tslop
      logical header
      character*4 pa,pb
      character*8 blank,pt
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     process keywords containing pisystem atom parameters
c
      blank = '        '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'PIATOM ') then
            if (header) then
               header = .false.
               write (iout,10)
   10          format (/,' Additional Pisystem Atom Parameters :',
     &                 //,6x,'Atom Type',3x,'Electron',3x,'Ionization',
     &                    3x,'Repulsion',/)
            end if
            ia = 0
            elect = 0.0d0
            ioniz = 0.0d0
            repuls = 0.0d0
            string = record(next:120)
            read (string,*,err=20)  ia,elect,ioniz,repuls
   20       continue
            write (iout,30)  ia,elect,ioniz,repuls
   30       format (8x,i4,3x,f10.3,3x,f10.3,2x,f10.3)
            if (ia.gt.0 .and. ia.le.maxclass) then
               electron(ia) = elect
               ionize(ia) = ioniz
               repulse(ia) = repuls
            else
   40          format (/,' KORBIT  --  Too many Atom Classes;',
     &                    ' Increase MAXCLASS')
               abort = .true.
            end if
         end if
      end do
c
c     process keywords containing pisystem bond parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'PIBOND ') then
            if (header) then
               header = .false.
               write (iout,50)
   50          format (/,' Additional Pisystem Bond Parameters :',
     &                 //,6x,'Atom Types',7x,'d Force',4x,'d Length',/)
            end if
            ia = 0
            ib = 0
            sslop = 0.0d0
            tslop = 0.0d0
            string = record(next:120)
            read (string,*,err=60)  ia,ib,sslop,tslop
   60       continue
            write (iout,70)  ia,ib,sslop,tslop
   70       format (6x,2i4,5x,2f11.3)
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            do j = 1, maxnpi
               if (kpi(j).eq.blank .or. kpi(j).eq.pt) then
                  kpi(j) = pt
                  sslope(j) = sslop
                  tslope(j) = tslop
                  goto 90
               end if
            end do
            write (iout,80)
   80       format (/,' KORBIT  --  Too many Pi-System Bond',
     &                 ' Type Parameters')
            abort = .true.
   90       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      npi = maxnpi
      do i = maxnpi, 1, -1
         if (kpi(i) .eq. blank)  npi = i - 1
      end do
c
c     assign the values characteristic of the piatom types;
c     count the number of filled pi molecular orbitals
c
      nfill = 0
      do i = 1, norbit
         it = type(iorbit(i))
         q(i) = electron(it)
         w(i) = ionize(it) / evolt
         em(i) = repulse(it) / evolt
         nfill = nfill + nint(q(i))
      end do
      nfill = nfill / 2
c
c     assign parameters for all bonds between piatoms;
c     store the original bond lengths and force constants
c
      do i = 1, nbpi
         j = bpi(1,i)
         ia = bpi(2,i)
         ib = bpi(3,i)
         ita = class(iorbit(ia))
         itb = class(iorbit(ib))
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         if (ita .le. itb) then
            pt = pa//pb
         else
            pt = pb//pa
         end if
         do k = 1, npi
            if (kpi(k) .eq. pt) then
               bkpi(i) = bk(j)
               blpi(i) = bl(j)
               kslope(i) = sslope(k)
               lslope(i) = tslope(k)
               goto 110
            end if
         end do
         write (iout,100)  ita,itb
  100    format (/,' KORBIT  --  No Parameters for Pi-Bond',
     &              ' between Atom Types',2i4)
  110    continue
      end do
c
c     store the original torsional constants across pibonds
c
      do i = 1, ntpi
         j = tpi(1,i)
         torsp2(i) = tors2(1,j)
      end do
      return
      end
