 
C     $Id: active.F 3 2008-03-04 17:47:30Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine active  --  set the list of active atoms  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "active" sets the list of atoms that are used during
c     each potential energy function calculation
c
c
      subroutine active
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'usage.i'
      integer i,j,next
      integer nmobile,nfixed
      integer center,nsphere
      integer nlist,list(maxatm)
      integer mobile(maxatm)
      integer fixed(maxatm)
      real*8 xcenter,ycenter,zcenter
      real*8 radius,radius2,dist2
      character*20 keyword
      character*120 record
      character*120 string

#if 0
      Include 'couple.i'
      Logical Done
      Integer IFixed, II
#endif

c
c
c     set defaults for the numbers and lists of active atoms
c
      nuse = n
      do i = 1, n
         use(i) = .true.
      end do
      nmobile = 0
      nfixed = 0
      do i = 1, n
         mobile(i) = 0
         fixed(i) = 0
      end do
      nsphere = 0
c
c     get any keywords containing active atom parameters
c
      do j = 1, nkey
         next = 1
         record = keyline(j)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
c
c     get any lists of atoms whose coordinates are active
c
         if (keyword(1:7) .eq. 'ACTIVE ') then
            read (string,*,err=10,end=10)  (mobile(i),i=nmobile+1,n)
   10       continue
            dowhile (mobile(nmobile+1) .ne. 0)
               nmobile = nmobile + 1
            end do
c
c     get any lists of atoms whose coordinates are inactive
c
         else if (keyword(1:9) .eq. 'INACTIVE ') then
            read (string,*,err=20,end=20)  (fixed(i),i=nfixed+1,n)
   20       continue
            dowhile (fixed(nfixed+1) .ne. 0)
               nfixed = nfixed + 1
            end do
c
c     get the center and radius of the sphere of active atoms
c
         else if (keyword(1:7) .eq. 'SPHERE ') then
            center = 0
            xcenter = 0.0d0
            ycenter = 0.0d0
            zcenter = 0.0d0
            radius = 0.0d0
            read (string,*,err=30,end=30)  xcenter,ycenter,
     &                                     zcenter,radius
   30       continue
            if (radius .eq. 0.0d0) then
               read (string,*,err=60,end=60)  center,radius
               xcenter = x(center)
               ycenter = y(center)
               zcenter = z(center)
            end if
            nsphere = nsphere + 1
            if (nsphere .eq. 1) then
               nuse = 0
               do i = 1, n
                  use(i) = .false.
               end do
               if (verbose) then
                  write (iout,40)
   40             format (/,' Active Site Spheres used to',
     &                        ' Select Active Atoms :',
     &                     //,3x,'Atom Center',11x,'Coordinates',
     &                        12x,'Radius',6x,'# Active Atoms')
               end if
            end if
            radius2 = radius * radius
            do i = 1, n
               if (.not. use(i)) then
                  dist2 = (x(i)-xcenter)**2 + (y(i)-ycenter)**2
     &                            + (z(i)-zcenter)**2
                  if (dist2 .le. radius2) then
                     nuse = nuse + 1
                     use(i) = .true.
                  end if
               end if
            end do

#if 0
C     Modified by Shenglong Wang accroding to Yingkai
C     an atom attached to two inactive atoms will be inactive
C     It seems this part need to be fixed, Shenglong Wang
C     How about terminal atoms, such H atoms in H-O-H molecule,
C     if O is fixed and H is free
            Done = .False.
            Do while(.not. Done) 
               Done = .True.
               Do I = 1, N
                  If(Use(I)) then
                     IFixed = 0
                     Do II = 1, N12(I)
                        If(.not. Use(I12(II,I))) IFixed = IFixed + 1
                     End Do
                     If(IFixed .ge. 2) then
                        NUse = NUse - 1
                        Use(I) = .False.
                        Done = .False.
                     End If
                  End If
               End Do
            End Do
C     End of modification
#endif

            if (verbose) then
               write (iout,50)  center,xcenter,ycenter,
     &                          zcenter,radius,nuse
   50          format (2x,i8,6x,3f9.2,2x,f9.2,7x,i8)
            end if
   60       continue
         end if
      end do
c
c     set active atoms to those not on the inactive atom list
c
      i = 1
      dowhile (fixed(i) .ne. 0)
         if (fixed(i) .gt. 0) then
            use(fixed(i)) = .false.
            nuse = nuse - 1
            i = i + 1
         else
            do j = abs(fixed(i)), abs(fixed(i+1))
               use(j) = .false.
               nuse = nuse - 1
            end do
            i = i + 2
         end if
      end do
c
c     set active atoms to only those on the active atom list
c
      i = 1
      dowhile (mobile(i) .ne. 0)
         if (i .eq. 1) then
            nuse = 0
            do j = 1, n
               use(j) = .false.
            end do
         end if
         if (mobile(i) .gt. 0) then
            use(mobile(i)) = .true.
            nuse = nuse + 1
            i = i + 1
         else
            do j = abs(mobile(i)), abs(mobile(i+1))
               use(j) = .true.
               nuse = nuse + 1
            end do
            i = i + 2
         end if
      end do
c
c     output the final list of the active atoms
c
      if (debug .and. nuse.gt.0 .and. nuse.lt.n) then
         nlist = 0
         do i = 1, n
            if (use(i)) then
               nlist = nlist + 1
               list(nlist) = i
            end if
         end do
         write (iout,70)
   70    format (/,' List of Active Atoms for Energy',
     &              ' Calculations :',/)
         write (iout,80)  (list(i),i=1,nlist)
   80    format (3x,10i7)
      end if
      return
      end
