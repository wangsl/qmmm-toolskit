 
C     $Id: cutoffs.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine cutoffs  --  set distance and Hessian cutoffs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "cutoffs" initializes and stores spherical energy cutoff
c     distance windows, Hessian element and Ewald sum cutoffs,
c     and the pairwise neighbor generation method
c
c
      subroutine cutoffs
      implicit none
      include 'sizes.i'
      include 'bound.i'
      include 'cutoff.i'
      include 'hescut.i'
      include 'keys.i'
      integer i,next
      real*8 big,value
      logical truncate
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     set defaults for spherical energy cutoff distances
c
      big = 1.0d12
      if (use_bounds) then
         vdwcut = 9.0d0
         chgcut = 9.0d0
         dplcut = 9.0d0
         mpolecut = 9.0d0
         ewaldcut = 9.0d0
      else
         vdwcut = big
         chgcut = big
         dplcut = big
         mpolecut = big
         ewaldcut = big
      end if
c
c     set defaults for switching function starting distance
c
      vdwtaper = 0.90d0
      chgtaper = 0.65d0
      dpltaper = 0.75d0
      mpoletaper = 0.65d0
c
c     set default for Hessian element magnitude cutoff
c
      hesscut = 0.0d0
c
c     set default tapering style, method of lights and Ewald
c
      truncate = .false.
      use_lights = .false.
      use_ewald = .false.
c
c     search the keywords for various cutoff parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
c
c     get values related to use of Ewald summation
c
         if (keyword(1:6) .eq. 'EWALD ') then
            if (use_bounds)  use_ewald = .true.
         else if (keyword(1:13) .eq. 'EWALD-CUTOFF ') then
            read (string,*,err=10,end=10)  ewaldcut
c
c     get values for to tapering style and method of lights
c
         else if (keyword(1:9) .eq. 'TRUNCATE ') then
            truncate = .true.
         else if (keyword(1:7) .eq. 'LIGHTS ') then
            use_lights = .true.
c
c     get cutoff for the magnitude of Hessian elements
c
         else if (keyword(1:12) .eq. 'HESS-CUTOFF ') then
            read (string,*,err=10,end=10)  hesscut
c
c     get the cutoff radii for potential energy functions
c
         else if (keyword(1:7) .eq. 'CUTOFF ') then
            read (string,*,err=10,end=10)  value
            vdwcut = value
            chgcut = value
            dplcut = value
            mpolecut = value
            ewaldcut = value
         else if (keyword(1:11) .eq. 'VDW-CUTOFF ') then
            read (string,*,err=10,end=10)  vdwcut
         else if (keyword(1:11) .eq. 'CHG-CUTOFF ') then
            read (string,*,err=10,end=10)  chgcut
         else if (keyword(1:11) .eq. 'DPL-CUTOFF ') then
            read (string,*,err=10,end=10)  dplcut
         else if (keyword(1:13) .eq. 'MPOLE-CUTOFF ') then
            read (string,*,err=10,end=10)  mpolecut
c
c     get distance for initialization of energy switching
c
         else if (keyword(1:6) .eq. 'TAPER ') then
            read (string,*,err=10,end=10)  value
            vdwtaper = value
            chgtaper = value
            dpltaper = value
            mpoletaper = value
         else if (keyword(1:10) .eq. 'VDW-TAPER ') then
            read (string,*,err=10,end=10)  vdwtaper
         else if (keyword(1:10) .eq. 'CHG-TAPER ') then
            read (string,*,err=10,end=10)  chgtaper
         else if (keyword(1:10) .eq. 'DPL-TAPER ') then
            read (string,*,err=10,end=10)  dpltaper
         else if (keyword(1:12) .eq. 'MPOLE-TAPER ') then
            read (string,*,err=10,end=10)  mpoletaper
         end if
   10    continue
      end do
c
c     convert any distance percentages to absolute distances
c
      if (vdwtaper .lt. 1.0d0)  vdwtaper = vdwtaper * vdwcut
      if (chgtaper .lt. 1.0d0)  chgtaper = chgtaper * chgcut
      if (dpltaper .lt. 1.0d0)  dpltaper = dpltaper * dplcut
      if (mpoletaper .lt. 1.0d0)  mpoletaper = mpoletaper * mpolecut
c
c     apply truncation cutoffs if they were requested
c
      if (truncate) then
         vdwtaper = big
         chgtaper = big
         dpltaper = big
         mpoletaper = big
      end if
      return
      end
