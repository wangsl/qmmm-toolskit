 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program document  --  make documentation lists from source  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "document" generates a formatted description of all the code
c     modules or common blocks, an index of routines called by each
c     source code module, a listing of all valid keywords, a list of
c     include file dependencies as needed by a Unix-style Makefile,
c     or a formatted force field parameter set summary
c
c     note the logical variable "guide" should be set true to make
c     output suitable for inclusion in the TINKER User's Guide
c
c
      program document
      implicit none
      integer maxline
      integer maxunit
      integer maxword
      integer maxfunc
      parameter (maxline=100)
      parameter (maxunit=1000)
      parameter (maxword=1000)
      parameter (maxfunc=65)
      include 'iounit.i'
      integer i,j,k,mode
      integer idoc,isrc
      integer nkey,nunit
      integer next,leng
      integer start,last
      integer freeunit
      integer trimtext
      integer nexttext
      integer nline(maxunit)
      integer link(maxunit)
      logical exist,done,guide
      character*20 module
      character*20 keyword
      character*20 keylast
      character*20 fname1,fname2
      character*20 key(maxword)
      character*20 fname(maxfunc)
      character*120 docfile
      character*120 srcfile
      character*120 record
      character*120 string
      character*120 routine(maxunit)
      character*120 info(maxline,maxunit)
      character*2048 field
c
c     list of the Fortran functions in the TINKER package
c
      data fname / 'ADJACENT',  'ANORM',     'BETACF',    'BETAI',
     &             'BMAX',      'BNDERR',    'CHIRER',    'CJKM',
     &             'D1D2',      'DEPTH',     'DIST2',     'DOT',
     &             'ENERGY',    'ERF',       'ERFC',      'ERFINV',
     &             'FREEUNIT',  'GAMMLN',    'GDA1',      'GDA2',
     &             'GEOMETRY',  'INITERR',   'INVBETA',   'LOCERR',
     &             'MAXWELL',   'MCM1',      'MCMSTEP',   'MIDERR',
     &             'MINIMIZ1',  'MINIROT1',  'MINRIGID1', 'NEWTON1',
     &             'NEWTROT1',  'NEXTTEXT',  'NORMAL',    'NUMBER',
     &             'OPTIMIZ1',  'OPTIROT1',  'OPTRIGID1', 'PATH1',
     &             'POTNRG',    'PRECISE',   'PROPERTY',  'PSS1',
     &             'PSSRGD1',   'PSSROT1',   'PTINCY',    'RANDOM',
     &             'RMSFIT',    'ROTANG',    'ROTCHECK',  'SADDLE1',
     &             'SCAN1',     'SIGMOID',   'SNIFFER1',  'TORSER',
     &             'TOTERR',    'TRANSIT',   'TRIMTEXT',  'TRIPLE',
     &             'VDWERR',    'VECANG',    'WATSON1',   'XTALLAT1',
     &             'XTALMOL1' /
c
c
c     set flag to create output formatted for TINKER User's Guide
c
      guide = .true.
c
c     find out what documentation the user wants to generate
c
      call initial
      write (iout,10)
   10 format (/,' The TINKER Document Facility can Provide :',
     &        //,4x,'(1) List of Routines from a Source File',
     &        /,4x,'(2) List of Calls made by each Routine',
     &        /,4x,'(3) List of Common Blocks from Source',
     &        /,4x,'(4) List of the TINKER Option Keywords',
     &        /,4x,'(5) List of Include File Dependencies',
     &        /,4x,'(6) Documentation from a Parameter File')
      mode = 0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=20,end=20)  mode
   20 continue
      dowhile (mode.lt.1 .or. mode.gt.6)
         write (iout,30)
   30    format (/,' Enter the Number of the Desired Choice :  ',$)
         read (input,40,err=20)  mode
   40    format (i10)
      end do
c
c     get the filename and open the input source code file
c
      if (mode .ne. 6) then
         isrc = freeunit ()
         call nextarg (srcfile,exist)
         if (exist) then
            call suffix (srcfile,'txt')
            call version (srcfile,'old')
            inquire (file=srcfile,exist=exist)
         end if
         dowhile (.not. exist)
            write (iout,50)
   50       format (/,' Enter Name of Source Code Listing File :  ',$)
            read (input,60)  srcfile
   60       format (a120)
            call suffix (srcfile,'txt')
            call version (srcfile,'old')
            inquire (file=srcfile,exist=exist)
         end do
         open (unit=isrc,file=srcfile,status='old')
      end if
c
c     get a list of routines and descriptions from source code
c
      if (mode .eq. 1) then
         nunit = 0
         dowhile (.true.)
            read (isrc,70,err=100,end=100)  record
   70       format (a120)
            if (record(1:9) .eq. 'c     ## ') then
               next = 10
               call getword (record,module,next)
               call lowcase (module)
               call upcase (module(1:1))
               if (module.eq.'Subroutine' .or. module.eq.'Function'
     &                      .or. module.eq.'Program') then
                  nunit = nunit + 1
                  call getword (record,routine(nunit),next)
                  call upcase (routine(nunit))
                  leng = trimtext (routine(nunit))
                  routine(nunit) = routine(nunit)(1:leng)//' '//module
                  read (isrc,80,err=100,end=100)
   80             format (///)
                  k = 0
                  done = .false.
                  dowhile (.not. done)
                     read (isrc,90,err=100,end=100)  record
   90                format (a120)
                     leng = trimtext (record)
                     if (leng .lt. 7) then
                        done = .true.
                     else if (record(1:1) .eq. ' ') then
                        done = .true.
                     else
                        k = k + 1
                        info(k,nunit) = record(7:leng)
                     end if
                  end do
                  nline(nunit) = k
               end if
            end if
         end do
  100    continue
         close (unit=isrc)
         call sort7 (nunit,routine,link)
         idoc = freeunit ()
         docfile = 'routines.doc'
         call version (docfile,'new')
         open (unit=idoc,file=docfile,status='new')
         do i = 1, nunit
            string = routine(i)
            leng = trimtext (string)
            write (idoc,110)  string(1:leng)
  110       format (a,/)
            last = 0
            j = link(i)
            do k = 1, nline(j)
               string = info(k,j)
               if (guide) then
                  if (k .eq. 1) then
                     leng = trimtext (string)
                     field(1:leng) = string(1:leng)
                     last = leng
                  else
                     last = last + 1
                     field(last:last) = ' '
                     leng = trimtext (string)
                     field(last+1:last+leng) = string(1:leng)
                     last = last + leng
                  end if
               else
                  string = info(k,j)
                  leng = trimtext (string)
                  write (idoc,120)  string(1:leng)
  120             format (a)
               end if
            end do
            if (guide .and. last.ne.0) then
               write (idoc,130)  field(1:last)
  130          format (a)
            end if
            if (nline(j) .ne. 0) then
               write (idoc,140)
  140          format ()
            end if
         end do
         close (unit=idoc)
         write (iout,150)  docfile(1:trimtext(docfile))
  150    format (/,' Source Documentation Written To:  ',a)
      end if
c
c     get a list of the calls made by each source code routine
c
      if (mode .eq. 2) then
         nunit = 0
         dowhile (.true.)
            read (isrc,160,err=170,end=170)  record
  160       format (a120)
            call upcase (record)
            if (record(1:1) .ne. 'C') then
               next = 1
               call getword (record,module,next)
               if (module.eq.'SUBROUTINE' .or. module.eq.'FUNCTION'
     &                      .or. module.eq.'PROGRAM') then
                  nunit = nunit + 1
                  call getword (record,routine(nunit),next)
                  nline(nunit) = 0
               else
                  next = index (record,' CALL ')
                  if (next .ne. 0) then
                     next = next + 6
                     call getword (record,keyword,next)
                     nline(nunit) = nline(nunit) + 1
                     info(nline(nunit),nunit) = keyword
                  else
                     do i = 1, maxfunc
                        leng = trimtext (fname(i))
                        fname1 = fname(i)(1:leng)//'('
                        fname2 = fname(i)(1:leng)//' ('
                        if (index(record,fname1(1:leng+1)).ne.0 .or.
     &                      index(record,fname2(1:leng+2)).ne.0) then
                           nline(nunit) = nline(nunit) + 1
                           info(nline(nunit),nunit) = fname(i)
                        end if
                     end do
                  end if
               end if
            end if
         end do
  170    continue
         close (unit=isrc)
         call sort7 (nunit,routine,link)
         idoc = freeunit ()
         docfile = 'calls.doc'
         call version (docfile,'new')
         open (unit=idoc,file=docfile,status='new')
         do i = 1, nunit
            string = routine(i)
            leng = trimtext (string)
            j = link(i)
            call sort10 (nline(j),info(1,j))
            if (guide) then
               field = string(1:leng)
               do k = 1, nline(j)
                  leng = trimtext (info(k,j))
                  last = trimtext (field)
                  field = field(1:last)//'    '//info(k,j)(1:leng)
               end do
               leng = trimtext (field)
               write (idoc,180)  field(1:leng)
  180          format (a)
            else
               write (idoc,190)  string(1:leng)
  190          format (a)
               do k = 1, nline(j)
                  string = info(k,j)
                  leng = trimtext (string)
                  write (idoc,200)  string(1:leng)
  200             format (5x,a)
               end do
            end if
         end do
         close (unit=idoc)
         write (iout,210)  docfile(1:trimtext(docfile))
  210    format (/,' Source Documentation Written To:  ',a)
      end if
c
c     get a list of common block contents from source code
c
      if (mode .eq. 3) then
         nunit = 0
         dowhile (.true.)
            read (isrc,220,err=250,end=250)  record
  220       format (a120)
            if (record(1:9) .eq. 'c     ## ') then
               next = index (record,'.i  --')
               if (next .ne. 0) then
                  nunit = nunit + 1
                  leng = trimtext (record)
                  call upcase (record(11:next-1))
                  string = record(11:next-1)
                  start = 20
                  if (guide)  start = trimtext(string) + 5
                  string(start:120) = record(next+8:leng-4)
                  routine(nunit) = string
                  read (isrc,230,err=250,end=250)
  230             format (///)
                  k = 0
                  done = .false.
                  dowhile (.not. done)
                     read (isrc,240,err=250,end=250)  record
  240                format (a120)
                     leng = trimtext (record)
                     if (record(1:1) .eq. ' ') then
                        done = .true.
                     else if (leng .ge. 7) then
                        k = k + 1
                        next = 7
                        call getword (record,string,next)
                        record = record(next:leng)
                        next = nexttext (record)
                        leng = trimtext (record)
                        start = 20
                        if (guide)  start = trimtext(string) + 5
                        string(start:120) = record(next:leng)
                        info(k,nunit) = string
                     end if
                  end do
                  nline(nunit) = k
               end if
            end if
         end do
  250    continue
         close (unit=isrc)
         call sort7 (nunit,routine,link)
         idoc = freeunit ()
         docfile = 'common.doc'
         call version (docfile,'new')
         open (unit=idoc,file=docfile,status='new')
         do i = 1, nunit
            string = routine(i)
            leng = trimtext (string)
            write (idoc,260)  string(1:leng)
  260       format (a,/)
            j = link(i)
            do k = 1, nline(j)
               string = info(k,j)
               leng = trimtext (string)
               write (idoc,270)  string(1:leng)
  270          format (a)
            end do
            if (nline(j) .ne. 0) then
               write (idoc,280)
  280          format ()
            end if
         end do
         close (unit=idoc)
         write (iout,290)  docfile(1:trimtext(docfile))
  290    format (/,' Source Documentation Written To:  ',a)
      end if
c
c     get the keyword values from the source code listing
c
      if (mode .eq. 4) then
         nkey = 0
         dowhile (.true.)
            read (isrc,300,err=310,end=310)  record
  300       format (a120)
            next = index (record,'if (keyword(')
            if (next .ne. 0) then
               next = index (record,'.eq.')
               if (next .ne. 0) then
                  next = index (record,'''')
                  if (next .ne. 0) then
                     next = next + 1
                     call getword (record,keyword,next)
                     call upcase (keyword)
                     nkey = nkey + 1
                     key(nkey) = keyword
                  end if
               end if
            end if
         end do
  310    continue
         close (unit=isrc)
         call sort6 (nkey,key)
         keylast = '                    '
         idoc = freeunit ()
         docfile = 'keyword.doc'
         call version (docfile,'new')
         open (unit=idoc,file=docfile,status='new')
         do i = 1, nkey
            keyword = key(i)
            leng = trimtext (keyword)
            if (keyword .ne. keylast) then
               write (idoc,320)  keyword(1:leng)
  320          format (a)
               keylast = keyword
            end if
         end do
         close (unit=idoc)
         write (iout,330)  docfile(1:trimtext(docfile))
  330    format (/,' Keyword Listing Written To:  ',a)
      end if
c
c     get the included files from the source code listing
c
      if (mode .eq. 5) then
         nkey = 0
         dowhile (.true.)
            read (isrc,340,err=350,end=350)  record
  340       format (a120)
            next = 1
            call getword (record,keyword,next)
            if (keyword .eq. 'include') then
               call gettext (record,keyword,next)
               keyword = keyword(2:trimtext(keyword)-1)
               nkey = nkey + 1
               key(nkey) = keyword
            end if
         end do
  350    continue
         close (unit=isrc)
         call sort6 (nkey,key)
         keylast = '                    '
         leng = index (srcfile,'.')
         field = srcfile(1:leng-1)//'.o:'
         do i = 1, nkey
            keyword = key(i)
            leng = trimtext (keyword)
            if (keyword .ne. keylast) then
               last = trimtext (field)
               field = field(1:last)//' '//keyword(1:leng)
               keylast = keyword
            end if
         end do
         write (iout,360)
  360    format (/,' File Dependencies in Makefile Format :',/)
         leng = trimtext (field)
         write (iout,370)  field (1:leng)
  370    format (a)
      end if
c
c     get a force field parameter file and write a listing
c
      if (mode .eq. 6) then
         call getprm
         idoc = freeunit ()
         docfile = 'parameter.doc'
         call version (docfile,'new')
         open (unit=idoc,file=docfile,status='new')
         call prtprm (idoc)
         close (unit=idoc)
         write (iout,380)  docfile(1:trimtext(docfile))
  380    format (/,' Parameter Listing Written To:  ',a)
      end if
c
c     perform any final tasks before program exit
c
      call final
      end
