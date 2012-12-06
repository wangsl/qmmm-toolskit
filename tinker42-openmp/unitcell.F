 
C     $Id: unitcell.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine unitcell  --  get periodic boundary conditions  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "unitcell" gets the periodic boundary box size and related
c     values from an external keyword file
c
c
      subroutine unitcell
      implicit none
      include 'sizes.i'
      include 'bound.i'
      include 'boxes.i'
      include 'iounit.i'
      include 'keys.i'
      integer i,next
      real*8 boxmax
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     set the default values for periodic boundary conditions
c
      use_bounds = .false.
      use_image = .false.
      use_replica = .false.
c
c     set the default values for the unitcell variables
c
      xbox = 0.0d0
      ybox = 0.0d0
      zbox = 0.0d0
      alpha = 0.0d0
      beta = 0.0d0
      gamma = 0.0d0
      orthogonal = .false.
      monoclinic = .false.
      triclinic = .false.
      octahedron = .false.
      spacegrp = '          '
c
c     get keywords containing crystal lattice dimensions
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:7) .eq. 'A-AXIS ') then
            read (string,*,err=10)  xbox
         else if (keyword(1:7) .eq. 'B-AXIS ') then
            read (string,*,err=10)  ybox
         else if (keyword(1:7) .eq. 'C-AXIS ') then
            read (string,*,err=10)  zbox
         else if (keyword(1:6) .eq. 'ALPHA ') then
            read (string,*,err=10)  alpha
         else if (keyword(1:5) .eq. 'BETA ') then
            read (string,*,err=10)  beta
         else if (keyword(1:6) .eq. 'GAMMA ') then
            read (string,*,err=10)  gamma
         else if (keyword(1:11) .eq. 'OCTAHEDRON ') then
            octahedron = .true.
         else if (keyword(1:11) .eq. 'SPACEGROUP ') then
            call getword (record,spacegrp,next)
         end if
   10    continue
      end do
c
c     use periodic boundary conditions if a cell was defined
c
      boxmax = max(xbox,ybox,zbox)
      if (boxmax .ne. 0.0d0) then
         use_bounds = .true.
         use_image = .true.
      end if
c
c     set unspecified periodic boundary box lengths and angles
c
      if (use_bounds) then
         if (xbox .eq. 0.0d0)  xbox = boxmax
         if (ybox .eq. 0.0d0)  ybox = xbox
         if (zbox .eq. 0.0d0)  zbox = xbox
         if (alpha .eq. 0.0d0)  alpha = 90.0d0
         if (beta .eq. 0.0d0)  beta = 90.0d0
         if (gamma .eq. 0.0d0)  gamma = 90.0d0
c
c     determine the general periodic boundary lattice type
c
         if (alpha.eq.90.0d0 .and. beta.eq.90.0d0
     &          .and. gamma.eq.90.0d0) then
            orthogonal = .true.
         else if (alpha.eq.90.0d0 .and. gamma.eq.90.0d0) then
            monoclinic = .true.
         else
            triclinic = .true.
         end if
      end if
c
c     check for proper use of truncated octahedron boundary
c
      if (octahedron) then
         if (xbox.eq.ybox .and. xbox.eq.zbox .and. orthogonal) then
            orthogonal = .false.
         else
            write (iout,20)
   20       format (/,' UNITCELL  --  Truncated Octahedron',
     &                 ' Incompatible with Defined Cell')
            call fatal
         end if
      end if
      return
      end
