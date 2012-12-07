 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program archive  --  create or extract from an archive  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "archive" is a utility program for coordinate files which
c     concatenates multiple coordinate sets into a single archive
c     file, or extracts individual coordinate sets from an archive
c
c
      program archive
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'iounit.i'
      integer i,k,iarc,ixyz
      integer start,stop,step
      integer lext,next,now
      integer last,nold
      integer leng1,leng2
      integer freeunit,trimtext
      logical exist,query
      character*1 letter,answer
      character*4 output
      character*7 ext
      character*8 mode
      character*120 record
      character*120 string
      character*120 arcfile
      character*120 coordfile
c
c
c     get the name to use for the coordinate archive file
c
      call initial
      call nextarg (filename,exist)
      if (.not. exist) then
         write (iout,10)
   10    format (/,' Enter Name of the Coordinate Archive File :  ',$)
         read (input,20)  filename
   20    format (a120)
      end if
      call basefile (filename)
      iarc = freeunit ()
c
c     decide whether to compress or extract files
c
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,30)
   30    format (/,' Create (C) or Extract (E) from an Archive',
     &              ' [C] :  ',$)
         read (input,40)  record
   40    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'E') then
         mode = 'EXTRACT'
      else
         mode = 'COMPRESS'
      end if
c
c     make separate files into a single concatenated file
c
      if (mode .eq. 'COMPRESS') then
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,50)
   50       format (/,' Use TINKER (T), XMol XYZ (X) or',
     &                 ' InsightII (I) File Format [T] :  ',$)
            read (input,60)  record
   60       format (a120)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         last = leng
         do i = 1, leng
            letter = filename(i:i)
            if (letter .eq. '/')  last = leng
c           if (letter .eq. '\')  last = leng
            if (ichar(letter) .eq. 92)  last = leng
            if (letter .eq. ']')  last = leng
            if (letter .eq. ':')  last = leng
            if (letter .eq. '~')  last = leng
            if (letter .eq. '.')  last = i - 1
         end do
         leng = min(leng,last)
         if (answer .eq. 'X') then
            output = 'XMOL'
            arcfile = filename(1:leng)
            call suffix (arcfile,'xmol')
         else if (answer .eq. 'I') then
            output = 'CAR'
            arcfile = filename(1:leng)
            call suffix (arcfile,'car')
         else
            output = 'XYZ'
            arcfile = filename(1:leng)
            call suffix (arcfile,'arc')
         end if
         start = 0
         stop = 0
         step = 0
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=70,end=70)  start
            query = .false.
         end if
         call nextarg (string,exist)
         if (exist)  read (string,*,err=70,end=70)  stop
         call nextarg (string,exist)
         if (exist)  read (string,*,err=70,end=70)  step
   70    continue
         if (query) then
            write (iout,80)
   80       format (/,' Numbers of First & Last File and Step',
     &                 ' Increment :  ',$)
            read (input,90)  record
   90       format (a120)
            read (record,*,err=100,end=100)  start,stop,step
  100       continue
         end if
         if (stop .eq. 0)  stop = start
         if (step .eq. 0)  step = 1
         call version (arcfile,'new')
         open (unit=iarc,file=arcfile,status='new')
c
c     cycle over the user specified coordinate files
c
         nold = 0
         i = start
         dowhile (i.ge.start .and. i.le.stop)
            ixyz = freeunit ()
            lext = 3
            call numeral (i,ext,lext)
            coordfile = filename(1:leng)//'.'//ext(1:lext)
            call version (coordfile,'old')
            inquire (file=coordfile,exist=exist)
            if (.not.exist .and. i.lt.100) then
               lext = 2
               call numeral (i,ext,lext)
               coordfile = filename(1:leng)//'.'//ext(1:lext)
               call version (coordfile,'old')
               inquire (file=coordfile,exist=exist)
            end if
            if (.not.exist .and. i.lt.10) then
               lext = 1
               call numeral (i,ext,lext)
               coordfile = filename(1:leng)//'.'//ext(1:lext)
               call version (coordfile,'old')
               inquire (file=coordfile,exist=exist)
            end if
            if (exist) then
               open (unit=ixyz,file=coordfile,status='old')
               rewind (unit=ixyz)
               call readxyz (ixyz)
               close (unit=ixyz)
               if (n .ne. nold)  call active
               if (output .eq. 'XMOL') then
                  call prtxmol (iarc)
               else if (output .eq. 'CAR') then
                  call prtcar (iarc)
               else
                  call prtarc (iarc)
               end if
            end if
            i = i + step
            nold = n
         end do
c
c     extract coordinates from a concatenated file
c
      else if (mode .eq. 'EXTRACT') then
  110    continue
         arcfile = filename(1:leng)
         call suffix (arcfile,'arc')
         call version (arcfile,'old')
         inquire (file=arcfile,exist=exist)
         if (exist) then
            open (unit=iarc,file=arcfile,status='old')
         else
            write (iout,120)
  120       format (/,' Enter Name of the Coordinate Archive',
     &                 ' File :  ',$)
            read (input,130)  filename
  130       format (a120)
            leng = trimtext (filename)
            goto 110
         end if
         rewind (unit=iarc)
         call readxyz (iarc)
         rewind (unit=iarc)
         now = 1
         leng1 = 1
         leng2 = leng
         do i = 1, leng
            if (filename(i:i) .eq. '/')  leng1 = i+1
            if (filename(i:i) .eq. ']')  leng1 = i+1
            if (filename(i:i) .eq. ':')  leng1 = i+1
         end do
         do i = leng, leng1, -1
            if (filename(i:i) .eq. '.')  leng2 = i-1
         end do
         leng = leng2 - leng1 + 1
         filename(1:leng) = filename(leng1:leng2)
         start = 0
         stop = 0
         step = 0
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=140,end=140)  start
            query = .false.
         end if
         call nextarg (string,exist)
         if (exist)  read (string,*,err=140,end=140)  stop
         call nextarg (string,exist)
         if (exist)  read (string,*,err=140,end=140)  step
  140    continue
         if (query) then
            write (iout,150)
  150       format (/,' Numbers of First & Last File and Step',
     &                 ' [<CR>=Exit] :  ',$)
            read (input,160)  record
  160       format (a120)
            read (record,*,err=170,end=170)  start,stop,step
  170       continue
         end if
         if (stop .eq. 0)  stop = start
         if (step .eq. 0)  step = 1
         dowhile (start .ne. 0)
            if (start .le. now) then
               now = 1
               rewind (unit=iarc)
            end if
            do k = 1, start-now
               call readxyz (iarc)
            end do
            i = start
            dowhile (i.ge.start .and. i.le.stop)
               lext = 3
               call numeral (i,ext,lext)
               call readxyz (iarc)
               if (n .eq. 0)  goto 180
               ixyz = freeunit ()
               coordfile = filename(1:leng)//'.'//ext(1:lext)
               call version (coordfile,'new')
               open (unit=ixyz,file=coordfile,status='new')
               call prtxyz (ixyz)
               close (unit=ixyz)
               i = i + step
               do k = 1, step-1
                  call readxyz (iarc)
               end do
            end do
  180       continue
            now = stop
            start = 0
            stop = 0
            step = 0
            query = .true.
            call nextarg (string,exist)
            if (exist) then
               read (string,*,err=190,end=190)  start
               query = .false.
            end if
            call nextarg (string,exist)
            if (exist)  read (string,*,err=190,end=190)  stop
            call nextarg (string,exist)
            if (exist)  read (string,*,err=190,end=190)  step
  190       continue
            if (query) then
               write (iout,200)
  200          format (/,' Numbers of First & Last File and Step',
     &                    ' [<CR>=Exit] :  ',$)
               read (input,210)  record
  210          format (a120)
               read (record,*,err=220,end=220)  start,stop,step
  220          continue
            end if
            if (stop .eq. 0)  stop = start
            if (step .eq. 0)  step = 1
         end do
      end if
      close (unit=iarc)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine prtarc  --  output of a TINKER archive file  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "prtarc" writes out a set of Cartesian coordinates for
c     all active atoms in the TINKER XYZ archive format
c
c
      subroutine prtarc (iarc)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'files.i'
      include 'inform.i'
      include 'titles.i'
      include 'usage.i'
      integer i,k,iarc
      logical opened
      character*120 coordfile
c
c
c     open output unit if not already done
c
      inquire (unit=iarc,opened=opened)
      if (.not. opened) then
         coordfile = filename(1:leng)//'.arc'
         call version (coordfile,'new')
         open (unit=iarc,file=coordfile,status='new')
      end if
c
c     write out the number of atoms and the title
c
      if (ltitle .eq. 0) then
         write (iarc,10)  nuse
   10    format (i6)
      else
         write (iarc,20)  nuse,title(1:ltitle)
   20    format (i6,2x,a)
      end if
c
c     finally, write the coordinates for each atom
c
      if (digits .le. 6) then
         do i = 1, n
            if (use(i)) then
               write (iarc,30)  i,name(i),x(i),y(i),z(i),type(i),
     &                          (i12(k,i),k=1,n12(i))
   30          format (i6,2x,a3,3f12.6,9i6)
            end if
         end do
      else if (digits .le. 8) then
         do i = 1, n
            if (use(i)) then
               write (iarc,40)  i,name(i),x(i),y(i),z(i),type(i),
     &                          (i12(k,i),k=1,n12(i))
   40          format (i6,2x,a3,3f14.8,9i6)
            end if
         end do
      else
         do i = 1, n
            if (use(i)) then
               write (iarc,50)  i,name(i),x(i),y(i),z(i),type(i),
     &                          (i12(k,i),k=1,n12(i))
   50          format (i6,2x,a3,3f16.10,9i6)
            end if
         end do
      end if
      if (.not. opened)  close (unit=iarc)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtxmol  --  output of a generic XYZ archive  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtxmol" writes out a set of Cartesian coordinates for
c     all active atoms in a simple, generic XYZ format originally
c     used by the XMOL program
c
c
      subroutine prtxmol (iarc)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'files.i'
      include 'titles.i'
      include 'usage.i'
      integer i,k,iarc
      logical opened
      character*3 atmnam
      character*120 coordfile
c
c
c     open output unit if not already done
c
      inquire (unit=iarc,opened=opened)
      if (.not. opened) then
         coordfile = filename(1:leng)//'.xmol'
         call version (coordfile,'new')
         open (unit=iarc,file=coordfile,status='new')
      end if
c
c     write out the number of active atoms and the title
c
      if (ltitle .eq. 0) then
         write (iarc,10)  nuse
   10    format (i6,/)
      else
         write (iarc,20)  nuse,title(1:ltitle)
   20    format (i6,2x,a,/)
      end if
c
c     write the atom name and coordinates for each active atom
c
      do i = 1, n
         if (use(i)) then
            atmnam = name(i)
            do k = 1, 3
               if ((atmnam(k:k).lt.'a'.or.atmnam(k:k).gt.'z') .and.
     &             (atmnam(k:k).lt.'A'.or.atmnam(k:k).gt.'Z'))
     &            atmnam(k:k) = ' '
            end do
            if (atmnam(1:1) .eq. 'H')  atmnam(2:3) = '  '
            if (atmnam(1:1) .eq. 'C')  atmnam(2:3) = '  '
            if (atmnam(1:1) .eq. 'N')  atmnam(2:3) = '  '
            if (atmnam(1:1) .eq. 'O')  atmnam(2:3) = '  '
            if (atmnam(1:1) .eq. 'S')  atmnam(2:3) = '  '
            write (iarc,30)  atmnam,x(i),y(i),z(i)
   30       format (5x,a3,2x,3f12.6)
         end if
      end do
      if (.not. opened)  close (unit=iarc)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtcar  --  output of InsightII archive file  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtcar" writes out a set of Cartesian coordinates for
c     all active atoms in the Accelerys InsightII .car format
c
c
      subroutine prtcar (iarc)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'files.i'
      include 'titles.i'
      include 'usage.i'
      integer i,k,iarc
      logical opened
      character*120 coordfile
c
c
c     open output unit if not already done
c
      inquire (unit=iarc,opened=opened)
      if (.not. opened) then
         coordfile = filename(1:leng)//'.car'
         call version (coordfile,'new')
         open (unit=iarc,file=coordfile,status='new')
      end if
c
c     write out the title and the date/time line
c
      if (ltitle .eq. 0) then
         write (iarc,10)
   10    format (/,'date:')
      else
         write (iarc,20)  title(1:ltitle)
   20    format (a,/,'date:')
      end if
c
c     write the coordinates for each active atom
c
      k = 0
      do i = 1, n
         if (use(i)) then
            k = k + 1
            write (iarc,30)  name(i),x(i),y(i),z(i),k
   30       format (a3,2x,3f15.6,4x,i6)
         end if
      end do
c
c     finally, write the end lines for this frame
c
      write (iarc,40)
   40 format ('end',/,'end')
      if (.not. opened)  close (unit=iarc)
      return
      end
