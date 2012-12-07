 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kpolar  --  assign polarizability parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kpolar" assigns atomic dipole polarizabilities to the atoms
c     within the structure and processes any new or changed values
c
c
      subroutine kpolar
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'kpolr.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'potent.i'
      integer i,j,k
      integer npg,next
      integer pg(maxval)
      real*8 pol,sixth
      logical header
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     process keywords containing polarizability parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'POLARIZE ') then
            k = 0
            pol = 0.0d0
            do j = 1, maxval
               pg(j) = 0
            end do
            string = record(next:120)
            read (string,*,err=10,end=10)  k,pol,(pg(j),j=1,maxval)
   10       continue
            if (k .gt. 0) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Atomic Dipole',
     &                       ' Polarizability Parameters :',
     &                    //,5x,'Atom Type',11x,'Alpha',9x,
     &                       'Group Atom Types'/)
               end if
               if (k .le. maxtyp) then
                  polr(k) = pol
                  do j = 1, maxval
                     pgrp(j,k) = pg(j)
                     if (pg(j) .eq. 0) then
                        npg = j - 1
                        goto 30
                     end if
                  end do
   30             continue
                  write (iout,40)  k,pol,(pg(j),j=1,npg)
   40             format (4x,i6,8x,f12.3,9x,20i5)
               else
                  write (iout,50)
   50             format (/,' KPOLAR  --  Too many Dipole',
     &                       ' Polarizability Parameters')
                  abort = .true.
               end if
            end if
         end if
      end do
c
c     find and store all the atomic dipole polarizabilities
c
      do i = 1, n
         polarity(i) = polr(type(i))
      end do
c
c     process keywords containing atom specific polarizabilities
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'POLARIZE ') then
            k = 0
            pol = 0.0d0
            string = record(next:120)
            read (string,*,err=80,end=80)  k,pol
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               if (header) then
                  header = .false.
                  write (iout,60)
   60             format (/,' Additional Dipole Polarizabilities',
     &                       ' for Specific Atoms :',
     &                    //,6x,'Atom',15x,'Alpha',/)
               end if
               write (iout,70)  k,pol
   70          format (4x,i6,8x,f12.3)
               polarity(k) = pol
            end if
   80       continue
         end if
      end do
c
c     remove zero and undefined polarizable multipoles from the list
c
      npole = 0
      npolar = 0
      do i = 1, n
         if (polsiz(i).ne.0 .or. polarity(i).ne.0.0d0) then
            if (zaxis(i).ne.0 .and. xaxis(i).ne.0) then
               npole = npole + 1
               ipole(npole) = i
               zaxis(npole) = zaxis(i)
               xaxis(npole) = xaxis(i)
               yaxis(npole) = yaxis(i)
               polaxe(npole) = polaxe(i)
               polsiz(npole) = polsiz(i)
               do k = 1, maxpole
                  pole(k,npole) = pole(k,i)
               end do
               if (polarity(i) .ne. 0.0d0)  npolar = npolar + 1
               polarity(npole) = polarity(i)
            end if
         end if
      end do
c
c     set the values used in the damping of the polarizability
c
      if (pgamma .eq. 0.0d0) then
         do i = 1, npole
            pdamp(i) = 0.0d0
         end do
      else
         sixth = 1.0d0 / 6.0d0
         do i = 1, npole
            pdamp(i) = polarity(i)**sixth
         end do
      end if
c
c     assign polarization group connectivity of each atom
c
      call polargrp
c
c     test multipoles at chiral sites and invert if necessary
c
      call chkpole
c
c     turn off polarizable multipole potential if it is not used
c
      if (npole .eq. 0)  use_mpole = .false.
      if (npolar .eq. 0)  use_polar = .false.
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine polargrp  --  polarization group connectivity  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "polargrp" generates members of the polarization group of
c     each atom and separate lists of the 1-2, 1-3 and 1-4 group
c     connectivities
c
c
      subroutine polargrp
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'kpolr.i'
      include 'mpole.i'
      include 'polgrp.i'
      integer i,j,k
      integer it,jt
      integer jj,kk
      integer start,stop
      integer nlist,nkeep
      integer list(maxatm)
      integer keep(maxatm)
      integer mask(maxatm)
      logical done
c
c
c     find the directly connected group members for each atom
c
      do i = 1, n
         np11(i) = 1
         ip11(1,i) = i
         it = type(i)
         do j = 1, n12(i)
            jj = i12(j,i)
            jt = type(jj)
            do k = 1, maxval
               kk = pgrp(k,it)
               if (kk .eq. 0)  goto 20
               if (pgrp(k,it) .eq. jt) then
                  np11(i) = np11(i) + 1
                  if (np11(i) .le. maxp11) then
                     ip11(np11(i),i) = jj
                  else
                     write (iout,10)
   10                format (/,' POLARGRP  --  Too many Atoms',
     &                          ' in Polarization Group')
                     abort = .true.
                  end if
               end if
            end do
   20       continue
         end do
      end do
c
c     find any other group members for each atom in turn
c
      do i = 1, n
         list(i) = 0
      end do
      do i = 1, n
         done = .false.
         start = 1
         stop = np11(i)
         do j = start, stop
            jj = ip11(j,i)
            if (jj .lt. i) then
               done = .true.
               np11(i) = np11(jj)
               do k = 1, np11(i)
                  ip11(k,i) = ip11(k,jj)
               end do
            else
               list(jj) = i
            end if
         end do
         dowhile (.not. done)
            done = .true.
            do j = start, stop
               jj = ip11(j,i)
               do k = 1, np11(jj)
                  kk = ip11(k,jj)
                  if (list(kk) .ne. i) then
                     np11(i) = np11(i) + 1
                     if (np11(i) .le. maxp11) then
                        ip11(np11(i),i) = kk
                     else
                        write (iout,30)
   30                   format (/,' POLARGRP  --  Too many Atoms',
     &                             ' in Polarization Group')
                        abort = .true.
                     end if
                     list(kk) = i
                  end if
               end do
            end do
            if (np11(i) .ne. stop) then
               done = .false.
               start = stop + 1
               stop = np11(i)
            end if
         end do
         call sort (np11(i),ip11(1,i))
      end do
c
c     loop over atoms finding all the 1-2 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         nlist = 0
         do j = 1, np11(i)
            jj = ip11(j,i)
            nlist = nlist + 1
            list(nlist) = jj
            mask(jj) = i
         end do
         nkeep = 0
         do j = 1, nlist
            jj = list(j)
            do k = 1, n12(jj)
               kk = i12(k,jj)
               if (mask(kk) .ne. i) then
                  nkeep = nkeep + 1
                  keep(nkeep) = kk
               end if
            end do
         end do
         nlist = 0
         do j = 1, nkeep
            jj = keep(j)
            do k = 1, np11(jj)
               kk = ip11(k,jj)
               nlist = nlist + 1
               list(nlist) = kk
            end do
         end do
         call sort8 (nlist,list)
         np12(i) = nlist
         do j = 1, nlist
            ip12(j,i) = list(j)
         end do
      end do
c
c     loop over atoms finding all the 1-3 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         do j = 1, np12(i)
            jj = ip12(j,i)
            mask(jj) = i
         end do
         nlist = 0
         do j = 1, np12(i)
            jj = ip12(j,i)
            do k = 1, np12(jj)
               kk = ip12(k,jj)
               if (mask(kk) .ne. i) then
                  nlist = nlist + 1
                  list(nlist) = kk
               end if
            end do
         end do
         call sort8 (nlist,list)
         np13(i) = nlist
         do j = 1, nlist
            ip13(j,i) = list(j)
         end do
      end do
c
c     loop over atoms finding all the 1-4 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         do j = 1, np12(i)
            jj = ip12(j,i)
            mask(jj) = i
         end do
         do j = 1, np13(i)
            jj = ip13(j,i)
            mask(jj) = i
         end do
         nlist = 0
         do j = 1, np13(i)
            jj = ip13(j,i)
            do k = 1, np12(jj)
               kk = ip12(k,jj)
               if (mask(kk) .ne. i) then
                  nlist = nlist + 1
                  list(nlist) = kk
               end if
            end do
         end do
         call sort8 (nlist,list)
         np14(i) = nlist
         do j = 1, nlist
            ip14(j,i) = list(j)
         end do
      end do
      return
      end
