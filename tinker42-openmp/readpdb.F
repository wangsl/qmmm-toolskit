 
C     $Id: readpdb.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine readpdb  --  input of Protein Data Bank file  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "readpdb" gets a set of Protein Data Bank coordinates
c     from an external disk file
c
c
      subroutine readpdb (ipdb)
      implicit none
      include 'sizes.i'
      include 'files.i'
      include 'iounit.i'
      include 'pdb.i'
      include 'sequen.i'
      include 'titles.i'
      integer i,ipdb
      integer next,nres
      integer nalt,nins
      integer index,trimtext
      integer length,serial
      integer residue,reslast
      real*8 xx,yy,zz
      logical exist,opened
      logical model
      character*1 altloc,chain
      character*1 chnsym,altsym
      character*1 altlast,chnlast
      character*1 insert,inslast
      character*1 letter,chnatm(maxatm)
      character*3 resname
      character*3 namelast
      character*4 atmname
      character*6 remark
      character*20 blank,text
      character*20 alttyp
      character*20 chntyp,chntemp
      character*20 instyp,instemp
      character*120 pdbfile
      character*120 record
      character*120 string
c
c
c     open the input file if it has not already been done
c
      inquire (unit=ipdb,opened=opened)
      if (.not. opened) then
         pdbfile = filename(1:leng)//'.pdb'
         call version (pdbfile,'old')
         inquire (file=pdbfile,exist=exist)
         if (exist) then
            open (unit=ipdb,file=pdbfile,status='old')
            rewind (unit=ipdb)
         else
            write (iout,10)
   10       format (/,' READPDB  --  Unable to Find the Protein',
     &                 ' Data Bank File')
            call fatal
         end if
      end if
c
c     initialize the residue counter and residue name
c
      nres = 0
      reslast = maxres
      namelast = '   '
      blank = '                    '
c
c     initialize alternate site, chain and insertion lists
c
      nalt = 0
      nchain = 0
      nins = 0
      altlast = '#'
      chnlast = '#'
      inslast = '#'
      alttyp = blank
      chntyp = '####################'
      instyp = blank
c
c     extract header information from the PDB file
c
      dowhile (.true.)
         read (ipdb,20,err=30,end=30) record
   20    format (a120)
         call upcase (record(1:6))
         remark = record(1:6)
         if (remark .eq. 'HEADER') then
            title = record(11:70)
            ltitle = trimtext (title)
            goto 30
         end if
      end do
   30 continue
c
c     scan for alternate locations, multiple chains and inserts
c
      rewind (unit=ipdb)
      dowhile (.true.)
         read (ipdb,40,err=60,end=60)  record
   40    format (a120)
         call upcase (record)
         remark = record(1:6)
         string = record(7:120)
         if (remark.eq.'ATOM  ' .or. remark.eq.'HETATM') then
            read (string,50)  altloc,chain,insert
   50       format (10x,a1,4x,a1,4x,a1)
            if (altloc .ne. altlast) then
               if (index(alttyp,altloc) .eq. 0) then
                  nalt = nalt + 1
                  alttyp(nalt:nalt) = altloc
                  altlast = altloc
               end if
            end if
            if (chain .ne. chnlast) then
               if (index(chntyp,chain) .eq. 0) then
                  nchain = nchain + 1
                  chntyp(nchain:nchain) = chain
                  chnlast = chain
               end if
            end if
            if (insert .ne. inslast) then
               if (index(instyp,insert) .eq. 0) then
                  nins = nins + 1
                  instyp(nins:nins) = insert
                  inslast = insert
               end if
            end if
         end if
      end do
   60 continue
c
c     find out which of the alternate locations will be used
c
      altsym = ' '
      if (nalt .gt. 0) then
         call nextarg (altsym,exist)
         if (.not. exist) then
            string(1:3) = '['//alttyp(1:1)//']'
            length = 3
            do i = 2, nalt
               string = string(1:length)//' '//alttyp(i:i)
               length = length + 2
            end do
            write (iout,70)  string(1:length)
   70       format (/,' Enter a Set of Alternate Atom Locations',
     &                 ' from (',a,') :  ',$)
            read (input,80)  record
   80       format (a120)
            next = 1
            call gettext (record,altsym,next)
         end if
         if (altsym .eq. ' ')  altsym = alttyp(1:1)
         call upcase (altsym)
      end if
c
c     find out which of the multiple chains will be used
c
      if (nchain .gt. 1) then
         call nextarg (chntemp,exist)
         if (.not. exist) then
            chntemp = blank
            if (chntyp(1:1) .eq. ' ') then
               string = 'BLANK'
               length = 5
            else
               string(1:1) = chntyp(1:1)
               length = 1
            end if
            do i = 2, nchain
               if (chntyp(i:i) .eq. ' ') then
                  string = string(1:length)//' BLANK'
                  length = length + 6
               else
                  string = string(1:length)//' '//chntyp(i:i)
                  length = length + 2
               end if
            end do
            string = string(1:length)//' [ALL]'
            length = length + 6
            write (iout,90)  string(1:length)
   90       format (/,' Enter the Chain Names to Include',
     &                 ' (',a,') :  ',$)
            read (input,100)  chntemp
  100       format (a20)
         end if
         call upcase (chntemp)
         next = 1
         call gettext (chntemp,text,next)
         if (text.eq.blank .or. text.eq.'ALL ') then
            chntyp = chntyp(1:nchain)
         else
            nchain = 1
            chntyp = chntemp(1:1)
         end if
      end if
c
c     find out which of the insert records will be used
c
      if (nins .gt. 0) then
         call nextarg (instemp,exist)
         if (.not. exist) then
            instemp = blank
            string(1:1) = instyp(1:1)
            length = 1
            do i = 2, nins
               string = string(1:length)//' '//instyp(i:i)
               length = length + 2
            end do
            string = string(1:length)//' [ALL] NONE'
            length = length + 11
            write (iout,110)  string(1:length)
  110       format (/,' Enter the Insert Records to Include',
     &                 ' (',a,') :  ',$)
            read (input,120)  instemp
  120       format (a20)
         end if
         call upcase (instemp)
         next = 1
         call gettext (instemp,text,next)
         if (text.eq.blank .or. text.eq.'ALL ') then
            instyp = instyp(1:nins)
         else if (text .eq. 'NONE ') then
            instyp = blank
         else
            instyp = instemp
         end if
      end if
c
c     process individual atoms from the Protein Data Bank file
c
      do i = 1, nchain
         rewind (unit=ipdb)
         model = .false.
         chnsym = chntyp(i:i)
         dowhile (.true.)
            read (ipdb,130,err=210,end=210)  record
  130       format (a120)
            call upcase (record)
            remark = record(1:6)
            if (remark .eq. 'ATOM  ') then
               string = record(7:120)
               read (string,140)  serial,atmname,altloc,resname,
     &                            chain,residue,insert,xx,yy,zz
  140          format (i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3)
               if (chain .ne. chnsym)  goto 160
               if (altloc.ne.' ' .and. altloc.ne.altsym)  goto 160
               if (insert.ne.' ' .and. index(instyp,insert).eq.0)
     &            goto 160
               call fixpdb (resname,atmname)
               if (residue.ne.reslast .or. resname.ne.namelast
     &                     .or. insert.ne.inslast) then
                  nres = nres + 1
                  reslast = residue
                  namelast = resname
                  inslast = insert
                  if (nres .gt. maxres) then
                     write (iout,150)  maxres
  150                format (/,' READPDB  --  The Maximum of',i6,
     &                          ' Residues has been Exceeded')
                     call fatal
                  end if
               end if
               npdb = npdb + 1
               xpdb(npdb) = xx
               ypdb(npdb) = yy
               zpdb(npdb) = zz
               pdbtyp(npdb) = remark
               atmnam(npdb) = atmname
               resnam(npdb) = resname
               resnum(npdb) = nres
               chnatm(npdb) = chain
  160          continue
            else if (remark .eq. 'HETATM') then
               string = record(7:120)
               read (string,170)  serial,atmname,altloc,resname,
     &                            chain,residue,insert,xx,yy,zz
  170          format (i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3)
               if (chain .ne. chnsym)  goto 180
               if (altloc.ne.' ' .and. altloc.ne.altsym)  goto 180
               if (insert.ne.' ' .and. index(instyp,insert).eq.0)
     &            goto 180
               call fixpdb (resname,atmname)
               npdb = npdb + 1
               xpdb(npdb) = xx
               ypdb(npdb) = yy
               zpdb(npdb) = zz
               pdbtyp(npdb) = remark
               atmnam(npdb) = atmname
               resnam(npdb) = resname
               resnum(npdb) = 0
               chnatm(npdb) = chain
  180          continue
            else if (remark .eq. 'MODEL ') then
               if (model) then
                  write (iout,190)
  190             format (/,' READPDB  --  File contains Multiple',
     &                       ' Models; First one Used')
                  goto 210
               else
                  model = .true.
               end if
            end if
            if (npdb .gt. maxatm) then
               write (iout,200)  maxatm
  200          format (/,' READPDB  --  The Maximum of',i6,
     &                    ' Atoms has been Exceeded')
               call fatal
            end if
         end do
  210    continue
      end do
c
c     set the total sequence length and chain termini information
c
      nseq = npdb
      nchain = 0
      chnlast = '#'
      do i = 1, npdb
         if (pdbtyp(i) .eq. 'ATOM  ') then
            letter = chnatm(i)
            if (letter .ne. chnlast) then
               nchain = nchain + 1
               ichain(1,nchain) = resnum(i)
               chnnam(nchain) = letter
               chnlast = letter
            else
               ichain(2,nchain) = resnum(i)
            end if
         end if
      end do
      if (.not. opened)  close (unit=ipdb)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine fixpdb  --  standardize atom and residue names  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "fixpdb" corrects problems with PDB files by converting residue
c     and atom names to the forms used by TINKER
c
c
      subroutine fixpdb (resname,atmname)
      implicit none
      include 'sizes.i'
      include 'resdue.i'
      integer i
      character*3 resname
      character*4 atmname
      character*7 restype
c
c
c     convert 1-letter base names to 3-letter base names
c
      if (resname .eq. '  A')  resname = 'ADE'
      if (resname .eq. '  G')  resname = 'GUA'
      if (resname .eq. '  C')  resname = 'CYT'
      if (resname .eq. '  U')  resname = 'URA'
      if (resname .eq. '  T')  resname = 'THY'
c
c     convert any unusual histidine residue names
c
      if (resname .eq. 'HSD')  resname = 'HID'
      if (resname .eq. 'HSE')  resname = 'HIE'
      if (resname .eq. 'HSP')  resname = 'HIS'
      if (resname .eq. 'HIP')  resname = 'HIS'
c
c     decide whether residue is protein or nucleic acid
c
      restype = 'UNKNOWN'
      do i = 1, maxamino
         if (resname .eq. amino(i))  restype = 'PROTEIN'
      end do
      do i = 1, maxnuc
         if (resname .eq. nuclz(i))  restype = 'NUCLEIC'
      end do
c
c     convert any generically used unusual atom names
c
      if (atmname .eq. ' HN ')  atmname = ' H  '
c
c     convert any unusual atom names in terminal residues
c
      if (restype .eq. 'PROTEIN') then
         if (atmname .eq. ' H1 ')  atmname = '1H  '
         if (atmname .eq. ' HN1')  atmname = '1H  '
         if (atmname .eq. ' HT1')  atmname = '1H  '
         if (atmname .eq. ' H2 ')  atmname = '2H  '
         if (atmname .eq. ' HN2')  atmname = '2H  '
         if (atmname .eq. ' HT2')  atmname = '2H  '
         if (atmname .eq. ' H3 ')  atmname = '3H  '
         if (atmname .eq. ' HN3')  atmname = '3H  '
         if (atmname .eq. ' HT3')  atmname = '3H  '
         if (atmname .eq. ' O1 ')  atmname = ' O  '
         if (atmname .eq. ' OT1')  atmname = ' O  '
         if (atmname .eq. 'OCT1')  atmname = ' O  '
         if (atmname .eq. ' O2 ')  atmname = ' OXT'
         if (atmname .eq. ' OT2')  atmname = ' OXT'
         if (atmname .eq. 'OCT2')  atmname = ' OXT'
         if (atmname .eq. ' OT ')  atmname = ' OXT'
      end if
c
c     glycine residue  (GLY)
c
      if (resname .eq. 'GLY') then
         if (atmname .eq. ' HA1')  atmname = '1HA '
         if (atmname .eq. ' HA2')  atmname = '2HA '
c
c     alanine residue  (ALA)
c
      else if (resname .eq. 'ALA') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HB3')  atmname = '3HB '
c
c     valine residue  (VAL)
c
      else if (resname .eq. 'VAL') then
         if (atmname .eq. 'HG11')  atmname = '1HG1'
         if (atmname .eq. 'HG12')  atmname = '2HG1'
         if (atmname .eq. 'HG13')  atmname = '3HG1'
         if (atmname .eq. 'HG21')  atmname = '1HG2'
         if (atmname .eq. 'HG22')  atmname = '2HG2'
         if (atmname .eq. 'HG23')  atmname = '3HG2'
c
c     leucine residue  (LEU)
c
      else if (resname .eq. 'LEU') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. 'HD11')  atmname = '1HD1'
         if (atmname .eq. 'HD12')  atmname = '2HD1'
         if (atmname .eq. 'HD13')  atmname = '3HD1'
         if (atmname .eq. 'HD21')  atmname = '1HD2'
         if (atmname .eq. 'HD22')  atmname = '2HD2'
         if (atmname .eq. 'HD23')  atmname = '3HD2'
c
c     isoleucine residue  (ILE)
c
      else if (resname .eq. 'ILE') then
         if (atmname .eq. ' CD ')  atmname = ' CD1'
         if (atmname .eq. 'HG11')  atmname = '1HG1'
         if (atmname .eq. 'HG13')  atmname = '1HG1'
         if (atmname .eq. 'HG12')  atmname = '2HG1'
         if (atmname .eq. 'HG21')  atmname = '1HG2'
         if (atmname .eq. 'HG22')  atmname = '2HG2'
         if (atmname .eq. 'HG23')  atmname = '3HG2'
         if (atmname .eq. 'HD11')  atmname = '1HD1'
         if (atmname .eq. ' HD1')  atmname = '1HD1'
         if (atmname .eq. 'HD12')  atmname = '2HD1'
         if (atmname .eq. ' HD2')  atmname = '2HD1'
         if (atmname .eq. 'HD13')  atmname = '3HD1'
         if (atmname .eq. ' HD3')  atmname = '3HD1'
c
c     serine residue  (SER)
c
      else if (resname .eq. 'SER') then
         if (atmname .eq. ' OG1')  atmname = ' OG '
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = ' HG '
         if (atmname .eq. ' HOG')  atmname = ' HG '
c
c     threonine residue  (THR)
c
      else if (resname .eq. 'THR') then
         if (atmname .eq. ' OG ')  atmname = ' OG1'
         if (atmname .eq. ' CG ')  atmname = ' CG2'
         if (atmname .eq. ' HOG')  atmname = ' HG1'
         if (atmname .eq. 'HOG1')  atmname = ' HG1'
         if (atmname .eq. 'HG21')  atmname = '1HG2'
         if (atmname .eq. 'HG22')  atmname = '2HG2'
         if (atmname .eq. 'HG23')  atmname = '3HG2'
c
c     cysteine residue  (CYS)
c
      else if (resname .eq. 'CYS') then
         if (atmname .eq. ' SG1')  atmname = ' SG '
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = ' HG '
         if (atmname .eq. ' HSG')  atmname = ' HG '
c
c     proline residue  (PRO)
c
      else if (resname .eq. 'PRO') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = '1HG '
         if (atmname .eq. ' HG3')  atmname = '1HG '
         if (atmname .eq. ' HG2')  atmname = '2HG '
         if (atmname .eq. ' HD1')  atmname = '1HD '
         if (atmname .eq. ' HD3')  atmname = '1HD '
         if (atmname .eq. ' HD2')  atmname = '2HD '
c
c     phenylalanine residue  (PHE)
c
      else if (resname .eq. 'PHE') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
c
c     tyrosine residue  (TYR)
c
      else if (resname .eq. 'TYR') then
         if (atmname .eq. ' HOH')  atmname = ' HH '
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
c
c     tryptophan residue  (TRP)
c
      else if (resname .eq. 'TRP') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HNE')  atmname = ' HE1'
c
c     histidine (HD and HE) residue  (HIS)
c
      else if (resname .eq. 'HIS') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HD ')  atmname = ' HD2'
         if (atmname .eq. ' HE ')  atmname = ' HE1'
         if (atmname .eq. ' HND')  atmname = ' HD1'
         if (atmname .eq. 'HND1')  atmname = ' HD1'
         if (atmname .eq. ' HNE')  atmname = ' HE2'
         if (atmname .eq. 'HNE2')  atmname = ' HE2'
c
c     histidine (HD only) residue  (HID)
c
      else if (resname .eq. 'HID') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HD ')  atmname = ' HD2'
         if (atmname .eq. ' HE ')  atmname = ' HE1'
         if (atmname .eq. ' HND')  atmname = ' HD1'
         if (atmname .eq. 'HND1')  atmname = ' HD1'
c
c     histidine (HE only) residue  (HIE)
c
      else if (resname .eq. 'HIE') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HD ')  atmname = ' HD2'
         if (atmname .eq. ' HE ')  atmname = ' HE1'
         if (atmname .eq. ' HNE')  atmname = ' HE2'
         if (atmname .eq. 'HNE2')  atmname = ' HE2'
c
c     aspartic acid residue  (ASP)
c
      else if (resname .eq. 'ASP') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
c
c     asparagine residue  (ASN)
c
      else if (resname .eq. 'ASN') then
         if (atmname .eq. ' OD ')  atmname = ' OD1'
         if (atmname .eq. ' ND ')  atmname = ' ND2'
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. 'HD21')  atmname = '1HD2'
         if (atmname .eq. 'HND1')  atmname = '1HD2'
         if (atmname .eq. 'HD22')  atmname = '2HD2'
         if (atmname .eq. 'HND2')  atmname = '2HD2'
c
c     glutamic acid residue  (GLU)
c
      else if (resname .eq. 'GLU') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = '1HG '
         if (atmname .eq. ' HG3')  atmname = '1HG '
         if (atmname .eq. ' HG2')  atmname = '2HG '
c
c     glutamine residue  (GLN)
c
      else if (resname .eq. 'GLN') then
         if (atmname .eq. ' OE ')  atmname = ' OE1'
         if (atmname .eq. ' NE ')  atmname = ' NE2'
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = '1HG '
         if (atmname .eq. ' HG3')  atmname = '1HG '
         if (atmname .eq. ' HG2')  atmname = '2HG '
         if (atmname .eq. 'HE21')  atmname = '1HE2'
         if (atmname .eq. 'HNE1')  atmname = '1HE2'
         if (atmname .eq. 'HE22')  atmname = '2HE2'
         if (atmname .eq. 'HNE2')  atmname = '2HE2'
c
c     methionine residue  (MET)
c
      else if (resname .eq. 'MET') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = '1HG '
         if (atmname .eq. ' HG3')  atmname = '1HG '
         if (atmname .eq. ' HG2')  atmname = '2HG '
         if (atmname .eq. ' HE1')  atmname = '1HE '
         if (atmname .eq. ' HE2')  atmname = '2HE '
         if (atmname .eq. ' HE3')  atmname = '3HE '
c
c     lysine residue  (LYS)
c
      else if (resname .eq. 'LYS') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = '1HG '
         if (atmname .eq. ' HG3')  atmname = '1HG '
         if (atmname .eq. ' HG2')  atmname = '2HG '
         if (atmname .eq. ' HD1')  atmname = '1HD '
         if (atmname .eq. ' HD3')  atmname = '1HD '
         if (atmname .eq. ' HD2')  atmname = '2HD '
         if (atmname .eq. ' HE1')  atmname = '1HE '
         if (atmname .eq. ' HE3')  atmname = '1HE '
         if (atmname .eq. ' HE2')  atmname = '2HE '
         if (atmname .eq. ' HZ1')  atmname = '1HZ '
         if (atmname .eq. 'HNZ1')  atmname = '1HZ '
         if (atmname .eq. ' HZ2')  atmname = '2HZ '
         if (atmname .eq. 'HNZ2')  atmname = '2HZ '
         if (atmname .eq. ' HZ3')  atmname = '3HZ '
         if (atmname .eq. 'HNZ3')  atmname = '3HZ '
c
c     arginine residue  (ARG)
c
      else if (resname .eq. 'ARG') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = '1HG '
         if (atmname .eq. ' HG3')  atmname = '1HG '
         if (atmname .eq. ' HG2')  atmname = '2HG '
         if (atmname .eq. ' HD1')  atmname = '1HD '
         if (atmname .eq. ' HD3')  atmname = '1HD '
         if (atmname .eq. ' HD2')  atmname = '2HD '
         if (atmname .eq. 'HH11')  atmname = '1HH1'
         if (atmname .eq. 'HN11')  atmname = '1HH1'
         if (atmname .eq. 'HH12')  atmname = '2HH1'
         if (atmname .eq. 'HN12')  atmname = '2HH1'
         if (atmname .eq. 'HH21')  atmname = '1HH2'
         if (atmname .eq. 'HN21')  atmname = '1HH2'
         if (atmname .eq. 'HH22')  atmname = '2HH2'
         if (atmname .eq. 'HN22')  atmname = '2HH2'
c
c     ornithine residue  (ORN)
c
      else if (resname .eq. 'ORN') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = '1HG '
         if (atmname .eq. ' HG3')  atmname = '1HG '
         if (atmname .eq. ' HG2')  atmname = '2HG '
         if (atmname .eq. ' HD1')  atmname = '1HD '
         if (atmname .eq. ' HD3')  atmname = '1HD '
         if (atmname .eq. ' HD2')  atmname = '2HD '
         if (atmname .eq. ' HE1')  atmname = '1HE '
         if (atmname .eq. 'HNE1')  atmname = '1HE '
         if (atmname .eq. ' HE2')  atmname = '2HE '
         if (atmname .eq. 'HNE2')  atmname = '2HE '
         if (atmname .eq. ' HE3')  atmname = '3HE '
         if (atmname .eq. 'HNE3')  atmname = '3HE '
c
c     methylalanine residue  (AIB)
c
      else if (resname .eq. 'AIB') then
         if (atmname .eq. 'HB11')  atmname = '1HB1'
         if (atmname .eq. 'HB12')  atmname = '2HB1'
         if (atmname .eq. 'HB13')  atmname = '3HB1'
         if (atmname .eq. 'HB21')  atmname = '1HB2'
         if (atmname .eq. 'HB22')  atmname = '2HB2'
         if (atmname .eq. 'HB23')  atmname = '3HB2'
c
c     pyroglutamic acid residue  (PCA)
c
      else if (resname .eq. 'PCA') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB3')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = '1HG '
         if (atmname .eq. ' HG3')  atmname = '1HG '
         if (atmname .eq. ' HG2')  atmname = '2HG '
c
c     N-terminal acetyl residue  (ACE)
c
      else if (resname .eq. 'ACE') then
         if (atmname .eq. ' CY ')  atmname = ' C  '
         if (atmname .eq. ' CAY')  atmname = ' CH3'
         if (atmname .eq. ' CA ')  atmname = ' CH3'
         if (atmname .eq. ' OY ')  atmname = ' O  '
         if (atmname .eq. ' H1 ')  atmname = '1H  '
         if (atmname .eq. ' H2 ')  atmname = '2H  '
         if (atmname .eq. ' H3 ')  atmname = '3H  '
         if (atmname .eq. ' HY1')  atmname = '1H  '
         if (atmname .eq. ' HY2')  atmname = '2H  '
         if (atmname .eq. ' HY3')  atmname = '3H  '
         if (atmname .eq. 'HH31')  atmname = '1H  '
         if (atmname .eq. 'HH32')  atmname = '2H  '
         if (atmname .eq. 'HH33')  atmname = '3H  '
c
c     N-terminal formyl residue  (FOR)
c
      else if (resname .eq. 'FOR') then
         if (atmname .eq. ' CY ')  atmname = ' C  '
         if (atmname .eq. ' OY ')  atmname = ' O  '
         if (atmname .eq. ' HY ')  atmname = ' H  '
c
c     C-terminal N-methylamide residue  (NME)
c
      else if (resname .eq. 'NME') then
         if (atmname .eq. ' NT ')  atmname = ' N  '
         if (atmname .eq. ' CT ')  atmname = ' CH3'
         if (atmname .eq. ' CAT')  atmname = ' CH3'
         if (atmname .eq. ' CA ')  atmname = ' CH3'
         if (atmname .eq. ' HNT')  atmname = ' H  '
         if (atmname .eq. ' H1 ')  atmname = '1H  '
         if (atmname .eq. ' H2 ')  atmname = '2H  '
         if (atmname .eq. ' H3 ')  atmname = '3H  '
         if (atmname .eq. ' HT1')  atmname = '1H  '
         if (atmname .eq. ' HT2')  atmname = '2H  '
         if (atmname .eq. ' HT3')  atmname = '3H  '
         if (atmname .eq. 'HH31')  atmname = '1H  '
         if (atmname .eq. 'HH32')  atmname = '2H  '
         if (atmname .eq. 'HH33')  atmname = '3H  '
c
c     C-terminal amide residue  (NH2)
c
      else if (resname .eq. 'NH2') then
         if (atmname .eq. ' NT ')  atmname = ' N  '
         if (atmname .eq. ' H1 ')  atmname = '1H  '
         if (atmname .eq. ' H2 ')  atmname = '2H  '
         if (atmname .eq. ' HT1')  atmname = '1H  '
         if (atmname .eq. ' HT2')  atmname = '2H  '
c
c     deoxythymidine residue  (THY)
c
      else if (resname .eq. 'THY') then
         if (atmname .eq. 'H5M1')  atmname = '1H5M'
         if (atmname .eq. 'H5M2')  atmname = '2H5M'
         if (atmname .eq. 'H5M3')  atmname = '3H5M'
      end if
      return
      end
