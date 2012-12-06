 
C     $Id: rings.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine rings  --  locate and store small rings  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "rings" searches the structure for small rings and stores
c     their constituent atoms
c
c     note code to remove reducible rings consisting of smaller
c     rings is commented in this version since reducible rings
c     are needed for parameter assignment
c
c
      subroutine rings
      implicit none
      include 'sizes.i'
      include 'angle.i'
c     include 'atoms.i'
      include 'bitor.i'
      include 'bond.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'ring.i'
      include 'tors.i'
      integer i,j,k,imax
      integer ia,ib,ic
      integer id,ie,ig
c     integer list1,list2
c     integer list3,list4
c     integer m,list(maxatm)
c
c
c     zero out the number of small rings in the structure
c
      nring3 = 0
      nring4 = 0
      nring5 = 0
      nring6 = 0
c
c     parse to find bonds, angles, torsions and bitorsions
c
      if (nbond .eq. 0)  call bonds
      if (nangle .eq. 0)  call angles
      if (ntors .eq. 0)  call torsions
      if (nbitor .eq. 0)  call bitors
c
c     search for and store all of the 3-membered rings
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         if (ib.lt.ia .and. ib.lt.ic) then
            do j = 1, n12(ia)
               if (i12(j,ia) .eq. ic) then
                  nring3 = nring3 + 1
                  if (nring3 .gt. maxring) then
                     write (iout,10)
   10                format (/,' RINGS  --  Too many 3-Membered Rings')
                     call fatal
                  end if
                  iring3(1,nring3) = ia
                  iring3(2,nring3) = ib
                  iring3(3,nring3) = ic
                  goto 20
               end if
            end do
   20       continue
         end if
      end do
c
c     search for and store all of the 4-membered rings
c
c     do i = 1, n
c        list(i) = 0
c     end do
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         if (ia.lt.ic .and. id.lt.ib) then
            do j = 1, n12(ia)
               if (i12(j,ia) .eq. id) then
                  nring4 = nring4 + 1
                  if (nring4 .gt. maxring) then
                     write (iout,30)
   30                format (/,' RINGS  --  Too many 4-Membered Rings')
                     call fatal
                  end if
                  iring4(1,nring4) = ia
                  iring4(2,nring4) = ib
                  iring4(3,nring4) = ic
                  iring4(4,nring4) = id
c
c     remove the ring if it is reducible into smaller rings
c
c                 list(ia) = nring4
c                 list(ib) = nring4
c                 list(ic) = nring4
c                 list(id) = nring4
c                 do m = 1, nring3
c                    list1 = list(iring3(1,m))
c                    list2 = list(iring3(2,m))
c                    list3 = list(iring3(3,m))
c                    if (list1.eq.nring4 .and. list2.eq.nring4
c    &                      .and. list3.eq.nring4) then
c                       nring4 = nring4 - 1
c                       list(ia) = 0
c                       list(ib) = 0
c                       list(ic) = 0
c                       list(id) = 0
c                       goto 40
c                    end if
c                 end do
                  goto 40
               end if
            end do
   40       continue
         end if
      end do
c
c     search for and store all of the 5-membered rings
c
c     do i = 1, n
c        list(i) = 0
c     end do
      do i = 1, nbitor
         ia = ibitor(1,i)
         ib = ibitor(2,i)
         ic = ibitor(3,i)
         id = ibitor(4,i)
         ie = ibitor(5,i)
         if (ia.lt.id .and. ie.lt.ib .and. min(ia,ie).lt.ic) then
            do j = 1, n12(ia)
               if (i12(j,ia) .eq. ie) then
                  nring5 = nring5 + 1
                  if (nring5 .gt. maxring) then
                     write (iout,50)
   50                format (/,' RINGS  --  Too many 5-Membered Rings')
                     call fatal
                  end if
                  iring5(1,nring5) = ia
                  iring5(2,nring5) = ib
                  iring5(3,nring5) = ic
                  iring5(4,nring5) = id
                  iring5(5,nring5) = ie
c
c     remove the ring if it is reducible into smaller rings
c
c                 list(ia) = nring5
c                 list(ib) = nring5
c                 list(ic) = nring5
c                 list(id) = nring5
c                 list(ie) = nring5
c                 do m = 1, nring3
c                    list1 = list(iring3(1,m))
c                    list2 = list(iring3(2,m))
c                    list3 = list(iring3(3,m))
c                    if (list1.eq.nring5 .and. list2.eq.nring5
c    &                      .and. list3.eq.nring5) then
c                       nring5 = nring5 - 1
c                       list(ia) = 0
c                       list(ib) = 0
c                       list(ic) = 0
c                       list(id) = 0
c                       list(ie) = 0
c                       goto 60
c                    end if
c                 end do
                  goto 60
               end if
            end do
   60       continue
         end if
      end do
c
c     search for and store all of the 6-membered rings
c
c     do i = 1, n
c        list(i) = 0
c     end do
      do i = 1, nbitor
         ia = ibitor(1,i)
         ib = ibitor(2,i)
         ic = ibitor(3,i)
         id = ibitor(4,i)
         ie = ibitor(5,i)
         imax = max(ia,ib,ic,id,ie)
         do j = 1, n12(ia)
            ig = i12(j,ia)
            if (ig .gt. imax) then
               do k = 1, n12(ie)
                  if (i12(k,ie) .eq. ig) then
                     nring6 = nring6 + 1
                     if (nring6 .gt. maxring) then
                        write (iout,70)
   70                   format (/,' RINGS  --  Too many',
     &                             ' 6-Membered Rings')
                        call fatal
                     end if
                     iring6(1,nring6) = ia
                     iring6(2,nring6) = ib
                     iring6(3,nring6) = ic
                     iring6(4,nring6) = id
                     iring6(5,nring6) = ie
                     iring6(6,nring6) = ig
c
c     remove the ring if it is reducible into smaller rings
c
c                    list(ia) = nring6
c                    list(ib) = nring6
c                    list(ic) = nring6
c                    list(id) = nring6
c                    list(ie) = nring6
c                    list(ig) = nring6
c                    do m = 1, nring3
c                       list1 = list(iring3(1,m))
c                       list2 = list(iring3(2,m))
c                       list3 = list(iring3(3,m))
c                       if (list1.eq.nring6 .and. list2.eq.nring6
c    &                         .and. list3.eq.nring6) then
c                          nring6 = nring6 - 1
c                          list(ia) = 0
c                          list(ib) = 0
c                          list(ic) = 0
c                          list(id) = 0
c                          list(ie) = 0
c                          list(ig) = 0
c                          goto 80
c                       end if
c                    end do
c                    do m = 1, nring4
c                       list1 = list(iring4(1,m))
c                       list2 = list(iring4(2,m))
c                       list3 = list(iring4(3,m))
c                       list4 = list(iring4(4,m))
c                       if (list1.eq.nring6 .and. list2.eq.nring6 .and.
c    &                      list3.eq.nring6 .and. list4.eq.nring6) then
c                          nring6 = nring6 - 1
c                          list(ia) = 0
c                          list(ib) = 0
c                          list(ic) = 0
c                          list(id) = 0
c                          list(ie) = 0
c                          list(ig) = 0
c                          goto 80
c                       end if
c                    end do
   80                continue
                  end if
               end do
            end if
         end do
      end do
c
c     print out lists of the small rings in the structure
c
      if (debug) then
         if (nring3 .gt. 0) then
            write (iout,90)
   90       format (/,' Three-Membered Rings Contained',
     &                 ' in the Structure :',
     &              //,11x,'Ring',14x,'Atoms in Ring',/)
            do i = 1, nring3
               write (iout,100)  i,(iring3(j,i),j=1,3)
  100          format (9x,i5,10x,3i6)
            end do
         end if
         if (nring4 .gt. 0) then
            write (iout,110)
  110       format (/,' Four-Membered Rings Contained',
     &                 ' in the Structure :',
     &              //,11x,'Ring',17x,'Atoms in Ring',/)
            do i = 1, nring4
               write (iout,120)  i,(iring4(j,i),j=1,4)
  120          format (9x,i5,10x,4i6)
            end do
         end if
         if (nring5 .gt. 0) then
            write (iout,130)
  130       format (/,' Five-Membered Rings Contained',
     &                 ' in the Structure :',
     &              //,11x,'Ring',20x,'Atoms in Ring',/)
            do i = 1, nring5
               write (iout,140)  i,(iring5(j,i),j=1,5)
  140          format (9x,i5,10x,5i6)
            end do
         end if
         if (nring6 .gt. 0) then
            write (iout,150)
  150       format (/,' Six-Membered Rings Contained',
     &                 ' in the Structure :',
     &              //,11x,'Ring',23x,'Atoms in Ring',/)
            do i = 1, nring6
               write (iout,160)  i,(iring6(j,i),j=1,6)
  160          format (9x,i5,10x,6i6)
            end do
         end if
      end if
      return
      end
