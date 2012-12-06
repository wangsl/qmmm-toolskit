 
C     $Id: prmedit.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2004  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program prmedit  --  edit and renumber parameter files  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "prmedit" reformats an existing parameter file, and revises
c     type and class numbers based on the "atom" parameter ordering
c
c
      program prmedit
      implicit none
      include 'iounit.i'
      integer iprm
      integer mode,nmode
      integer freeunit
      integer trimtext
      logical dotype,doclass
      character*120 prmfile
c
c
c     read and store the original force field parameter file
c
      call initial
      call getprm
c
c     get the desired type of coordinate file modification
c
      write (iout,10)
   10 format (/,' The Parameter File Edit Facility can Provide :',
     &        //,4x,'(1) Renumber Only the Atom Types',
     &        /,4x,'(2) Renumber Only the Atom Classes',
     &        /,4x,'(3) Renumber Both Atom Types and Classes'
     &        /,4x,'(4) Sort and Format Multipole Parameters')
   20 continue
      nmode = 4
      mode = 0
      write (iout,30)
   30 format (/,' Enter the Number of the Desired Choice :  ',$)
      read (input,40,err=20)  mode
   40 format (i10)
      if (mode.lt.0 .or. mode.gt.nmode)  goto 20
c
c     set the renumbering operations to be performed
c
      dotype = .false.
      doclass = .false.
      if (mode .eq. 1)  dotype = .true.
      if (mode .eq. 2)  doclass = .true.
      if (mode .eq. 3) then
         dotype = .true.
         doclass = .true.
      end if
c
c     format and renumber the atom types and atom classes
c
      if (mode.ge.1 .and. mode.le.3) then
         iprm = freeunit ()
         prmfile = 'parameter.prm'
         call version (prmfile,'new')
         open (unit=iprm,file=prmfile,status='new')
         call prmorder (iprm,dotype,doclass)
         write (iout,50)  prmfile(1:trimtext(prmfile))
   50    format (/,' Renumbered Parameter File Written To:  ',a)
         close (unit=iprm)
      end if
c
c     sort the atomic multipole parameters by atom type
c
      if (mode .eq. 4) then
         iprm = freeunit ()
         prmfile = 'multipole.prm'
         call version (prmfile,'new')
         open (unit=iprm,file=prmfile,status='new')
         call polesort (iprm)
         write (iout,60)  prmfile(1:trimtext(prmfile))
   60    format (/,' Sorted Multipole Values Written To:  ',a)
         close (unit=iprm)
      end if
      call final
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prmorder  --  reorder atom types and classes  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prmorder" places a list of atom type or class numbers into
c     canonical order for potential energy parameter definitions
c
c
      subroutine prmorder (iprm,dotype,doclass)
      implicit none
      include 'sizes.i'
      include 'iounit.i'
      include 'params.i'
      integer i,j,iprm
      integer it,ic,kt,kc
      integer ia,ib,id,ie
      integer offset,next
      integer length
      integer trimtext
      integer kg,ig(maxval)
      integer itype(0:maxtyp)
      integer iclass(0:maxclass)
      real*8 value
      logical dotype,doclass
      character*20 keyword
      character*30 blank
      character*120 record
      character*120 string
c
c
c     zero out the storage for atom types and classes
c
      ia = 0
      ib = 0
      ic = 0
      id = 0
      ie = 0
      kt = 0
      kc = 0
      do i = 0, maxtyp
         itype(i) = 0
      end do
      do i = 0, maxclass
         iclass(i) = 0
      end do
      blank = '                              '
c
c     get the starting numbers for atom types and classes
c
      if (dotype) then
         write (iout,10)
   10    format (/,' Enter Starting Number for Atom Types [1] :  ',$)
         read (input,20)  offset
   20    format (i10)
         if (offset .gt. 0)  kt = offset - 1
      end if
      if (doclass) then
         write (iout,30)
   30    format (/,' Enter Starting Number for Atom Classes [1] :  ',$)
         read (input,40)  offset
   40    format (i10)
         if (offset .gt. 0)  kc = offset - 1
      end if
c
c     count and order the atom types and atom classes
c
      do i = 1, nprm
         record = prmline(i)
         next = 1
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'ATOM ') then
            it = 0
            ic = 0
            call getnumb (record,it,next)
            call getnumb (record,ic,next)
            if (ic .eq. 0)  ic = it
            if (itype(it) .eq. 0) then
               kt = kt + 1
               if (dotype) then
                  itype(it) = kt
               else
                  itype(it) = it
               end if
            end if
            if (iclass(ic) .eq. 0) then
               kc = kc + 1
               if (doclass) then
                  iclass(ic) = kc
               else
                  iclass(ic) = ic
               end if
            end if
         end if
      end do
c
c     renumber and print the various parameter types
c
      do i = 1, nprm
         record = prmline(i)
         length = trimtext (record)
         next = 1
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'ATOM ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            if (ib .eq. 0)  ib = ia
            ia = itype(ia)
            ib = iclass(ib)
            write (iprm,50)  ia,ib,record(next:length)
   50       format ('atom',6x,2i5,a)
         else if (keyword(1:4) .eq. 'VDW ') then
            ia = 0
            call getnumb (record,ia,next)
            ia = iclass(ia)
            write (iprm,60)  ia,record(next:length)
   60       format ('vdw',7x,i5,a)
         else if (keyword(1:6) .eq. 'VDWPR ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,70)  ia,ib,record(next:length)
   70       format ('vdwpr',5x,2i5,a)
         else if (keyword(1:5) .eq. 'BOND ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,80)  ia,ib,record(next:length)
   80       format ('bond',6x,2i5,a)
         else if (keyword(1:6) .eq. 'BOND5 ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,90)  ia,ib,record(next:length)
   90       format ('bond5',5x,2i5,a)
         else if (keyword(1:6) .eq. 'BOND4 ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,100)  ia,ib,record(next:length)
  100       format ('bond4',5x,2i5,a)
         else if (keyword(1:6) .eq. 'BOND3 ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,110)  ia,ib,record(next:length)
  110       format ('bond3',5x,2i5,a)
         else if (keyword(1:9) .eq. 'ELECTNEG ') then
            ia = 0
            ib = 0
            ic = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            write (iprm,120)  ia,ib,ic,record(next:length)
  120       format ('electneg',2x,3i5,a)
         else if (keyword(1:6) .eq. 'ANGLE ') then
            ia = 0
            ib = 0
            ic = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            call prmsort (3,ia,ib,ic,0,0)
            write (iprm,130)  ia,ib,ic,record(next:length)
  130       format ('angle',5x,3i5,a)
         else if (keyword(1:7) .eq. 'ANGLE5 ') then
            ia = 0
            ib = 0
            ic = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            call prmsort (3,ia,ib,ic,0,0)
            write (iprm,140)  ia,ib,ic,record(next:length)
  140       format ('angle5',4x,3i5,a)
         else if (keyword(1:7) .eq. 'ANGLE4 ') then
            ia = 0
            ib = 0
            ic = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            call prmsort (3,ia,ib,ic,0,0)
            write (iprm,150)  ia,ib,ic,record(next:length)
  150       format ('angle4',4x,3i5,a)
         else if (keyword(1:7) .eq. 'ANGLE3 ') then
            ia = 0
            ib = 0
            ic = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            call prmsort (3,ia,ib,ic,0,0)
            write (iprm,160)  ia,ib,ic,record(next:length)
  160       format ('angle3',4x,3i5,a)
         else if (keyword(1:7) .eq. 'ANGLEF ') then
            ia = 0
            ib = 0
            ic = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            call prmsort (3,ia,ib,ic,0,0)
            write (iprm,170)  ia,ib,ic,record(next:length)
  170       format ('anglef',4x,3i5,a)
         else if (keyword(1:7) .eq. 'STRBND ') then
            ia = 0
            call getnumb (record,ia,next)
            ia = iclass(ia)
            write (iprm,180)  ia,record(next:length)
  180       format ('strbnd',4x,i5,a)
         else if (keyword(1:9) .eq. 'UREYBRAD ') then
            ia = 0
            ib = 0
            ic = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            call prmsort (3,ia,ib,ic,0,0)
            write (iprm,190)  ia,ib,ic,record(next:length)
  190       format ('ureybrad',2x,3i5,a)
         else if (keyword(1:7) .eq. 'ANGANG ') then
            ia = 0
            call getnumb (record,ia,next)
            ia = iclass(ia)
            write (iprm,200)  ia,record(next:length)
  200       format ('angang',4x,i5,a)
         else if (keyword(1:7) .eq. 'OPBEND ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            write (iprm,210)  ia,ib,record(next:length)
  210       format ('opbend',4x,2i5,a)
         else if (keyword(1:7) .eq. 'OPDIST ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            id = iclass(id)
            call prmsort (2,ib,ic,0,0,0)
            call prmsort (2,ib,id,0,0,0)
            call prmsort (2,ic,id,0,0,0)
            write (iprm,220)  ia,ib,ic,id,record(next:length)
  220       format ('opdist',4x,4i5,a)
         else if (keyword(1:9) .eq. 'IMPROPER ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            id = iclass(id)
            write (iprm,230)  ia,ib,ic,id,record(next:length)
  230       format ('improper',2x,4i5,a)
         else if (keyword(1:8) .eq. 'IMPTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            id = iclass(id)
            write (iprm,240)  ia,ib,ic,id,record(next:length)
  240       format ('imptors',3x,4i5,a)
         else if (keyword(1:8) .eq. 'TORSION ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            id = iclass(id)
            call prmsort (4,ia,ib,ic,id,0)
            write (iprm,250)  ia,ib,ic,id,record(next:length)
  250       format ('torsion',3x,4i5,a)
         else if (keyword(1:9) .eq. 'TORSION5 ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            id = iclass(id)
            call prmsort (4,ia,ib,ic,id,0)
            write (iprm,260)  ia,ib,ic,id,record(next:length)
  260       format ('torsion5',2x,4i5,a)
         else if (keyword(1:9) .eq. 'TORSION4 ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            id = iclass(id)
            call prmsort (4,ia,ib,ic,id,0)
            write (iprm,270)  ia,ib,ic,id,record(next:length)
  270       format ('torsion4',2x,4i5,a)
         else if (keyword(1:7) .eq. 'PITORS ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,280)  ia,ib,record(next:length)
  280       format ('pitors',4x,2i5,a)
         else if (keyword(1:8) .eq. 'STRTORS ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,290)  ia,ib,record(next:length)
  290       format ('strtors',3x,2i5,a)
         else if (keyword(1:8) .eq. 'TORTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            ie = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            call getnumb (record,ie,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            id = iclass(id)
            ie = iclass(ie)
            write (iprm,300)  ia,ib,ic,id,ie,record(next:length)
  300       format ('tortors',3x,5i5,a)
         else if (keyword(1:7) .eq. 'CHARGE ') then
            ia = 0
            call getnumb (record,ia,next)
            ia = itype(ia)
            write (iprm,310)  ia,record(next:length)
  310       format ('charge',4x,i5,a)
         else if (keyword(1:7) .eq. 'DIPOLE ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = itype(ia)
            ib = itype(ib)
            write (iprm,320)  ia,ib,record(next:length)
  320       format ('dipole',4x,2i5,a)
         else if (keyword(1:8) .eq. 'DIPOLE5 ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = itype(ia)
            ib = itype(ib)
            write (iprm,330)  ia,ib,record(next:length)
  330       format ('dipole5',3x,2i5,a)
         else if (keyword(1:8) .eq. 'DIPOLE4 ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = itype(ia)
            ib = itype(ib)
            write (iprm,340)  ia,ib,record(next:length)
  340       format ('dipole4',3x,2i5,a)
         else if (keyword(1:8) .eq. 'DIPOLE3 ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = itype(ia)
            ib = itype(ib)
            write (iprm,350)  ia,ib,record(next:length)
  350       format ('dipole3',3x,2i5,a)
         else if (keyword(1:10) .eq. 'MULTIPOLE ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = itype(ia)
            ib = isign(1,ib) * itype(abs(ib))
            ic = isign(1,ic) * itype(abs(ic))
            id = isign(1,id) * itype(abs(id))
            if (id .eq. 0) then
               write (iprm,360)  ia,ib,ic,record(next:length)
  360          format ('multipole',1x,3i5,a)
            else
               write (iprm,370)  ia,ib,ic,id,record(next:length)
  370          format ('multipole',1x,4i5,a)
            end if
         else if (keyword(1:9) .eq. 'POLARIZE ') then
            ia = 0
            value = 0.0d0
            kg = 0
            do j = 1, maxval
               ig(j) = 0
            end do
            string = record(next:120)
            read (string,*,err=380,end=380)  ia,value,(ig(j),j=1,maxval)
  380       continue
            ia = itype(ia)
            do j = 1, maxval
               if (ig(j) .ne. 0) then
                  kg = j
                  ig(j) = itype(ig(j))
               end if
            end do
            call sort (kg,ig)
            write (iprm,390)  ia,value,(ig(j),j=1,kg)
  390       format ('polarize',2x,i5,f20.3,5x,8i5)
         else if (keyword(1:7) .eq. 'PIATOM ') then
            ia = 0
            call getnumb (record,ia,next)
            ia = iclass(ia)
            write (iprm,400)  ia,record(next:length)
  400       format ('piatom',4x,i5,a)
         else if (keyword(1:7) .eq. 'PIBOND ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,410)  ia,ib,record(next:length)
  410       format ('pibond',4x,2i5,a)
         else if (keyword(1:8) .eq. 'PIBOND5 ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,420)  ia,ib,record(next:length)
  420       format ('pibond5',3x,2i5,a)
         else if (keyword(1:8) .eq. 'PIBOND4 ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,430)  ia,ib,record(next:length)
  430       format ('pibond4',3x,2i5,a)
         else if (keyword(1:6) .eq. 'METAL ') then
            ia = 0
            call getnumb (record,ia,next)
            ia = iclass(ia)
            write (iprm,440)  ia,record(next:length)
  440       format ('metal',5x,i5,a)
         else if (keyword(1:8) .eq. 'BIOTYPE ') then
            ia = 0
            ib = 0
            string = record(next:120)
            read (string,*,err=450,end=450)  ia
            call getword (record,string,next)
            call getstring (record,string,next)
            string = record(next:120)
            read (string,*,err=450,end=450)  ib
  450       continue
            if (ib .ne. 0) then
               ib = itype(ib)
               if (ib.eq.0 .and. kt.lt.999)  ib = 999
            end if
            length = min(30,max(1,51-next))
            write (iprm,460)  record(8:next)//blank(1:length),ib
  460       format ('biotype',a,i5)
         else if (length .eq. 0) then
            write (iprm,470)
  470       format ()
         else
            write (iprm,480)  record(1:length)
  480       format (a)
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine prmsort  --  reorder atom types and classes  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "prmsort" places a list of atom type or class numbers into
c     canonical order for potential energy parameter definitions
c
c
      subroutine prmsort (index,ia,ib,ic,id,ie)
      implicit none
      integer ia,ib,ic,id,ie
      integer index,temp
c
c
c     put the atom type or class numbers into canonical order
c
      if (index .eq. 2) then
         if (ia .gt. ib) then
            temp = ia
            ia = ib
            ib = temp
         end if
      else if (index .eq. 3) then
         if (ia .gt. ic) then
            temp = ia
            ia = ic
            ic = temp
         end if
      else if (index .eq. 4) then
         if (ib.gt.ic .or. (ib.eq.ic.and.ia.gt.id)) then
            temp = ib
            ib = ic
            ic = temp
            temp = ia
            ia = id
            id = temp
         end if
      else if (index .eq. 5) then
         if (ib.gt.id .or. (ib.eq.id.and.ia.gt.ie)) then
            temp = ib
            ib = id
            id = temp
            temp = ia
            ia = ie
            ie = temp
         end if
      end if
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine polesort  --  sort multipoles by atom type  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "polesort" sorts a set of atomic multipole parameters based
c     on the atom types of centers involved
c
c
      subroutine polesort (iprm)
      implicit none
      include 'sizes.i'
      include 'params.i'
      integer i,j,n,iprm
      integer size,next
      integer length
      integer trimtext
      integer ia,ib,ic,id
      integer key(maxprm)
      integer line(maxprm)
      real*8 v1,v2,v3
      character*4 pa,pb,pc,pd
      character*16 list(maxprm)
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     find and store atom types for the multipole parameters
c
      n = 0
      do i = 1, nprm
         record = prmline(i)
         length = trimtext (record)
         next = 1
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:10) .eq. 'MULTIPOLE ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = abs(ia)
            ib = abs(ib)
            ic = abs(ic)
            id = abs(id)
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            n = n + 1
            line(n) = i
            list(n) = pa//pb//pc//pd
         end if
      end do
c
c     sort the parameters based on the atom type numbers
c
      call sort7 (n,list,key)
c
c     format and output the sorted multipole parameters
c
      do i = 1, n
         j = line(key(i))
         record = prmline(j)
         next = 1
         call gettext (record,keyword,next)
         ia = 0
         ib = 0
         ic = 0
         id = 0
         string = record(next:120)
         read (string,*,err=20,end=20)  ia,ib,ic,id,v1
         write (iprm,10)  ia,ib,ic,id,v1
   10    format ('multipole ',4i5,5x,f11.5)
         goto 40
   20    continue
         read (string,*,err=90,end=90)  ia,ib,ic,v1
         write (iprm,30)  ia,ib,ic,v1
   30    format ('multipole ',3i5,10x,f11.5)
   40    continue
         j = j + 1
         record = prmline(j)
         read (record,*,err=90,end=90)  v1,v2,v3
         write (iprm,50)  v1,v2,v3
   50    format (35x,3f11.5)
         j = j + 1
         record = prmline(j)
         read (record,*,err=90,end=90)  v1
         write (iprm,60)  v1
   60    format (35x,f11.5)
         j = j + 1
         record = prmline(j)
         read (record,*,err=90,end=90)  v1,v2
         write (iprm,70)  v1,v2
   70    format (35x,2f11.5)
         j = j + 1
         record = prmline(j)
         read (record,*,err=90,end=90)  v1,v2,v3
         write (iprm,80)  v1,v2,v3
   80    format (35x,3f11.5)
   90    continue
      end do
      return
      end
