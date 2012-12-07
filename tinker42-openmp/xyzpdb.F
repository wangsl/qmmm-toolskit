 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  program xyzpdb  --  Cartesian to Protein Data Bank file  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "xyzpdb" takes as input a Cartesian coordinates file,
c     then converts to and writes out a Protein Data Bank file
c
c
      program xyzpdb
      implicit none
      include 'files.i'
      integer ipdb,freeunit
      character*120 pdbfile
c
c
c     get the Cartesian coordinates file for the system
c
      call initial
      call getxyz
c
c     get atomic number of each atom and count the molecules
c
      call field
      call katom
      call molecule
c
c     convert from Cartesian coordinates to PDB format
c
      call makepdb
c
c     write out the Protein Data Bank file
c
      pdbfile = filename(1:leng)//'.pdb'
      call version (pdbfile,'new')
      ipdb = freeunit ()
      open (unit=ipdb,file=pdbfile,status='new')
      call prtpdb (ipdb)
      close (unit=ipdb)
c
c     perform any final tasks before program exit
c
      call final
      end
