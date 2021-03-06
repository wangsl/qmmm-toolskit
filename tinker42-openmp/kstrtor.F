 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kstrtor  --  find stretch-torsion parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kstrtor" assigns stretch-torsion parameters to torsions
c     needing them, and processes any new or changed values
c
c
      subroutine kstrtor
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'ksttor.i'
      include 'potent.i'
      include 'strtor.i'
      include 'tors.i'
      integer i,j,k,nbt
      integer ia,ib
      integer ita,itb
      integer size,next
      real*8 bt1,bt2,bt3
      logical header
      character*4 pa,pb
      character*8 blank,pt
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     process keywords containing stretch-torsion parameters
c
      blank = '        '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'STRTORS ') then
            ia = 0
            ib = 0
            bt1 = 0.0d0
            bt2 = 0.0d0
            bt3 = 0.0d0
            string = record(next:120)
            read (string,*,err=10,end=10)  ia,ib,bt1,bt2,bt3
   10       continue
            if (header) then
               header = .false.
               write (iout,20)
   20          format (/,' Additional Stretch-Torsion Parameters :',
     &                 //,5x,'Atom Classes',7x,'1-Fold',6x,'2-Fold',
     &                    6x,'3-Fold',/)
            end if
            write (iout,30)  ia,ib,bt1,bt2,bt3
   30       format (6x,2i4,4x,3f12.3)
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            do j = 1, maxnbt
               if (kbt(j).eq.blank .or. kbt(j).eq.pt) then
                  kbt(j) = pt
                  btcon(1,j) = bt1
                  btcon(2,j) = bt2
                  btcon(3,j) = bt3
                  goto 50
               end if
            end do
            write (iout,40)
   40       format (/,' KSTRTOR  --  Too many Stretch-Torsion',
     &                 ' Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      nbt = maxnbt
      do i = maxnbt, 1, -1
         if (kbt(i) .eq. blank)  nbt = i - 1
      end do
c
c     assign stretch-torsion parameters as required
c
      nstrtor = 0
      if (nbt .ne. 0) then
         do i = 1, ntors
            ia = itors(2,i)
            ib = itors(3,i)
            ita = class(ia)
            itb = class(ib)
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            if (ita .le. itb) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            do j = 1, nbt
               if (kbt(j) .eq. pt) then
                  nstrtor = nstrtor + 1
                  kst(1,nstrtor) = btcon(1,j)
                  kst(2,nstrtor) = btcon(2,j)
                  kst(3,nstrtor) = btcon(3,j)
                  ist(1,nstrtor) = i
                  do k = 1, n12(ia)
                     if (i12(k,ia) .eq. ib) then
                        ist(2,nstrtor) = bndlist(k,ia)
                        goto 60
                     end if
                  end do
               end if
            end do
   60       continue
         end do
      end if
c
c     turn off the stretch-torsion potential if it is not used
c
      if (nstrtor .eq. 0)  use_strtor = .false.
      return
      end
