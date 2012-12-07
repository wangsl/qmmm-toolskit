 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine kmpole  --  multipole parameter assignment  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "kmpole" assigns atomic multipole moments to the atoms of
c     the structure and processes any new or changed values
c
c
      subroutine kmpole
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'kmulti.i'
      include 'mpole.i'
      include 'polgrp.i'
      include 'potent.i'
      include 'units.i'
      integer i,j,k,m
      integer jt,kk,mm
      integer it,mt,imp,nmp
      integer size,next
      integer big,number
      integer kz,kx,ky
      integer ztyp,xtyp,ytyp
      integer mpt(maxnmp)
      integer mpkey(maxnmp)
      integer mpz(maxnmp)
      integer mpx(maxnmp)
      integer mpy(maxnmp)
      integer start(maxtyp)
      integer stop(maxtyp)
      real*8 random,mpl(13)
      logical header
      character*4 pa,pb,pc,pd
      character*8 axt
      character*16 blank,pt
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     process keywords containing atomic multipole parameters
c
      blank = '                '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:10) .eq. 'MULTIPOLE ') then
            k = 0
            kz = 0
            kx = 0
            ky = 0
            axt = 'Z-then-X'
            do j = 1, 13
               mpl(j) = 0.0d0
            end do
            string = record(next:120)
            read (string,*,err=10,end=10)  k,kz,kx,ky,mpl(1)
            goto 20
   10       continue
            read (string,*,err=70,end=70)  k,kz,kx,mpl(1)
   20       continue
            if (k .gt. 0) then
               if (kz.lt.0 .or. kx.lt.0)  axt = 'Bisector'
               kz = abs(kz)
               kx = abs(kx)
               ky = abs(ky)
               record = keyline(i+1)
               read (record,*,err=70,end=70)  (mpl(j),j=2,4)
               record = keyline(i+2)
               read (record,*,err=70,end=70)  mpl(5)
               record = keyline(i+3)
               read (record,*,err=70,end=70)  (mpl(j),j=8,9)
               record = keyline(i+4)
               read (record,*,err=70,end=70)  (mpl(j),j=11,13)
               mpl(6) = mpl(8)
               mpl(7) = mpl(11)
               mpl(10) = mpl(12)
               if (header) then
                  header = .false.
                  write (iout,30)
   30             format (/,' Additional Atomic Multipole Parameters :',
     &                    //,5x,'Atom Type',5x,'Coordinate Frame',
     &                       ' Definition',8x,'Multipole Moments')
               end if
               write (iout,40)  k,kz,kx,ky,axt,(mpl(j),j=1,5),
     &                          mpl(8),mpl(9),(mpl(j),j=11,13)
   40          format (/,4x,i6,5x,i6,1x,i6,1x,i6,3x,a8,2x,f9.5,
     &                    /,48x,3f9.5,/,48x,f9.5,
     &                    /,48x,2f9.5,/,48x,3f9.5)
               size = 4
               call numeral (k,pa,size)
               call numeral (kz,pb,size)
               call numeral (kx,pc,size)
               call numeral (ky,pd,size)
               pt = pa//pb//pc//pd
               do j = 1, maxnmp
                  if (kmp(j).eq.blank .or. kmp(j).eq.pt) then
                     kmp(j) = pt
                     mpaxis(j) = axt
                     do m = 1, 13
                        multip(m,j) = mpl(m)
                     end do
                     goto 60
                  end if
               end do
               write (iout,50)
   50          format (/,' KMPOLE  --  Too many Atomic Multipole',
     &                    ' Parameters')
               abort = .true.
   60          continue
            end if
   70       continue
         end if
      end do
c
c     zero out local axes, multipoles and polarization attachments
c
      npole = n
      do i = 1, n
         zaxis(i) = 0
         xaxis(i) = 0
         yaxis(i) = 0
         polaxe(i) = '        '
         do j = 1, 13
            pole(j,i) = 0.0d0
         end do
         np11(i) = 0
         np12(i) = 0
         np13(i) = 0
         np14(i) = 0
      end do
c
c     determine the total number of forcefield parameters
c
      nmp = maxnmp
      do i = maxnmp, 1, -1
         if (kmp(i) .eq. blank)  nmp = i - 1
      end do
c
c     set indices into the list of atomic multipole parameters
c
      do i = 1, nmp
         mpt(i) = number(kmp(i)(1:4))
         mpz(i) = number(kmp(i)(5:8))
         mpx(i) = number(kmp(i)(9:12))
         mpy(i) = number(kmp(i)(13:16))
      end do
      call sort3 (nmp,mpt,mpkey)
      do i = 1, maxtyp
         start(i) = 1
         stop(i) = 0
      end do
      do i = 1, nmp
         k = mpt(i)
         if (start(k) .eq. 1)  start(k) = i
      end do
      do i = nmp, 1, -1
         k = mpt(i)
         if (stop(k) .eq. 0)  stop(k) = i
      end do
c
c     assign atomic multipole parameters to each atom
c
      do i = 1, n
         it = type(i)
         do jt = start(it), stop(it)
            imp = mpkey(jt)
            ztyp = mpz(imp)
            xtyp = mpx(imp)
            ytyp = mpy(imp)
            if (ztyp .eq. 0) then
               zaxis(i) = n + 1
               xaxis(i) = n + 2
               polaxe(i) = mpaxis(imp)
               do j = 1, 13
                  pole(j,i) = multip(j,imp)
               end do
               goto 80
            end if
            do k = 1, n12(i)
               if (type(i12(k,i)) .eq. ztyp) then
                  if (xtyp .eq. 0) then
                     zaxis(i) = i12(k,i)
                     xaxis(i) = n + 1
                     polaxe(i) = mpaxis(imp)
                     do j = 1, 13
                        pole(j,i) = multip(j,imp)
                     end do
                     goto 80
                  end if
                  do m = 1, n12(i)
                     mt = type(i12(m,i))
                     if (mt.eq.xtyp .and. m.ne.k) then
                        zaxis(i) = i12(k,i)
                        xaxis(i) = i12(m,i)
                        if (ytyp .eq. 0) then
                           polaxe(i) = mpaxis(imp)
                           do j = 1, 13
                              pole(j,i) = multip(j,imp)
                           end do
                           goto 80
                        end if
                        do kk = 1, n12(i)
                           mm = i12(kk,i)
                           mt = type(mm)
                           if (mm.ne.zaxis(i) .and. mm.ne.xaxis(i)) then
                              if (mt .eq. ytyp) then
                                 yaxis(i) = mm
                                 polaxe(i) = mpaxis(imp)
                                 do j = 1, 13
                                    pole(j,i) = multip(j,imp)
                                 end do
                                 goto 80
                              end if
                           end if
                        end do
                        do kk = 1, n13(i)
                           mm = i13(kk,i)
                           mt = type(mm)
                           if (mm.ne.zaxis(i) .and. mm.ne.xaxis(i)) then
                              if (mt .eq. ytyp) then
                                 yaxis(i) = mm
                                 polaxe(i) = mpaxis(imp)
                                 do j = 1, 13
                                    pole(j,i) = multip(j,imp)
                                 end do
                                 goto 80
                              end if
                           end if
                        end do
                     end if
                  end do
               end if
            end do
         end do
         do jt = start(it), stop(it)
            imp = mpkey(jt)
            ztyp = mpz(imp)
            xtyp = mpx(imp)
            do k = 1, n12(i)
               if (type(i12(k,i)) .eq. ztyp) then
                  do m = 1, n13(i)
                     mt = type(i13(m,i))
                     if (mt .eq. xtyp) then
                        zaxis(i) = i12(k,i)
                        xaxis(i) = i13(m,i)
                        polaxe(i) = mpaxis(imp)
                        do j = 1, 13
                           pole(j,i) = multip(j,imp)
                        end do
                        goto 80
                     end if
                  end do
               end if
            end do
         end do
   80    continue
      end do
c
c     process keywords with multipole parameters for specific atoms
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:10) .eq. 'MULTIPOLE ') then
            k = 0
            kz = 0
            kx = 0
            ky = 0
            axt = 'Z-then-X'
            do j = 1, 13
               mpl(j) = 0.0d0
            end do
            string = record(next:120)
            read (string,*,err=90,end=90)  k,kz,kx,ky,mpl(1)
            goto 100
   90       continue
            read (string,*,err=130,end=130)  k,kz,kx,mpl(1)
  100       continue
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               if (kz.lt.0 .or. kx.lt.0)  axt = 'Bisector'
               kz = abs(kz)
               kx = abs(kx)
               ky = abs(ky)
               record = keyline(i+1)
               read (record,*,err=130,end=130)  (mpl(j),j=2,4)
               record = keyline(i+2)
               read (record,*,err=130,end=130)  mpl(5)
               record = keyline(i+3)
               read (record,*,err=130,end=130)  (mpl(j),j=8,9)
               record = keyline(i+4)
               read (record,*,err=130,end=130)  (mpl(j),j=11,13)
               mpl(6) = mpl(8)
               mpl(7) = mpl(11)
               mpl(10) = mpl(12)
               if (header) then
                  header = .false.
                  write (iout,110)
  110             format (/,' Additional Atomic Multipoles',
     &                       ' for Specific Atoms :',
     &                    //,6x,'Atom',9x,'Coordinate Frame',
     &                       ' Definition',8x,'Multipole Moments')
               end if
               write (iout,120)  k,kz,kx,ky,axt,(mpl(j),j=1,5),
     &                           mpl(8),mpl(9),(mpl(j),j=11,13)
  120          format (/,4x,i6,5x,i6,1x,i6,1x,i6,3x,a8,2x,f9.5,
     &                    /,48x,3f9.5,/,48x,f9.5,
     &                    /,48x,2f9.5,/,48x,3f9.5)
               if (kz .eq. 0)  kz = n + 1
               if (kx .eq. 0)  kx = n + 2
               if (ky .eq. 0)  ky = n + 3
               zaxis(k) = kz
               xaxis(k) = kx
               yaxis(k) = ky
               polaxe(k) = axt
               do j = 1, 13
                  pole(j,k) = mpl(j)
               end do
            end if
  130       continue
         end if
      end do
c
c     convert the dipole and quadrupole moments to Angstroms,
c     quadrupole divided by 3 for use as traceless values
c
      do i = 1, n
         do k = 2, 4
            pole(k,i) = pole(k,i) * bohr
         end do
         do k = 5, 13
            pole(k,i) = pole(k,i) * bohr**2 / 3.0d0
         end do
      end do
c
c     get the order of the multipole expansion at each site
c
      do i = 1, n
         size = 0
         do k = 1, maxpole
            if (pole(k,i) .ne. 0.0d0)  size = max(k,size)
         end do
         if (size .gt. 4) then
            size = 13
         else if (size .gt. 1) then
            size = 4
         end if
         polsiz(i) = size
      end do
c
c     if needed, get random coordinates for dummy axis defining atoms
c
      big = 0
      do i = 1, n
         big = max(big,zaxis(i),xaxis(i),yaxis(i))
      end do
      if (big .gt. n) then
         do i = n+1, n+3
            x(i) = random ()
            y(i) = random ()
            z(i) = random ()
         end do
      end if
c
c     if polarization not used, remove zero and undefined multipoles
c
      if (.not. use_polar) then
         npole = 0
         do i = 1, n
            if (polsiz(i) .ne. 0) then
               if (zaxis(i).ne.0 .and. xaxis(i).ne.0) then
                  npole = npole + 1
                  ipole(npole) = i
                  zaxis(npole) = zaxis(i)
                  xaxis(npole) = xaxis(i)
                  yaxis(npole) = yaxis(i)
                  polaxe(npole) = polaxe(i)
                  polsiz(npole) = polsiz(i)
                  do j = 1, maxpole
                     pole(j,npole) = pole(j,i)
                  end do
               end if
            end if
         end do
c
c     test multipoles at chiral sites and invert if necessary
c
         call chkpole
c
c     turn off the atomic multipole potential if it is not used
c
         if (npole .eq. 0)  use_mpole = .false.
      end if
      return
      end
