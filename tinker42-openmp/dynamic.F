 
C     $Id: dynamic.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program dynamic  --  run molecular or stochastic dynamics  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dynamic" computes a molecular dynamics trajectory in any of
c     several statistical mechanical ensembles with optional periodic
c     boundaries and optional coupling to temperature and pressure baths
c     alternatively a stochastic dynamics trajectory can be generated
c
c
      program dynamic
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bond.i'
      include 'bound.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mdstuf.i'
      include 'potent.i'
      include 'solute.i'
      include 'stodyn.i'
      include 'usage.i'
      integer istep,nstep
      integer mode,modstep
      real*8 dt,dtdump
      logical exist,query
      character*120 string
c
c
c     set up the structure and molecular mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     initialize the temperature, pressure and coupling baths
c
      kelvin = 0.0d0
      atmsph = 0.0d0
      isothermal = .false.
      isobaric = .false.
c
c     initialize the simulation length as number of time steps
c
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  nstep
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Enter the Number of Dynamics Steps to be',
     &              ' Taken :  ',$)
         read (input,30)  nstep
   30    format (i10)
      end if
c
c     get the length of the dynamics time step in picoseconds
c
      dt = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=40,end=40)  dt
   40 continue
      if (dt .le. 0.0d0) then
         write (iout,50)
   50    format (/,' Enter the Time Step Length in Femtoseconds',
     &              ' [1.0] :  ',$)
         read (input,60)  dt
   60    format (f20.0)
      end if
      if (dt .le. 0.0d0)  dt = 1.0d0
      dt = 1.0d-3 * dt
c
c     set bounds on the Berendsen thermostat coupling parameters
c
      tautemp = max(tautemp,dt)
      taupres = max(taupres,dt)
c
c     set the time between trajectory snapshot coordinate dumps
c
      dtdump = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=70,end=70)  dtdump
   70 continue
      if (dtdump .le. 0.0d0) then
         write (iout,80)
   80    format (/,' Enter Time between Dumps in Picoseconds',
     &              ' [0.1] :  ',$)
         read (input,90)  dtdump
   90    format (f20.0)
      end if
      if (dtdump .le. 0.0d0)  dtdump = 0.1d0
      iwrite = nint(dtdump/dt)
c
c     get choice of statistical ensemble for periodic system
c
      if (use_bounds) then
         mode = 0
         call nextarg (string,exist)
         if (exist)  read (string,*,err=100,end=100)  mode
  100    continue
         if (mode.lt.1 .or. mode.gt.4) then
            write (iout,110)
  110       format (/,' Possible Statistical Mechanical Ensembles :',
     &              //,4x,'(1) Microcanonical (NVE)',
     &              /,4x,'(2) Canonical (NVT)',
     &              /,4x,'(3) Isoenthalpic-Isobaric (NPH)',
     &              /,4x,'(4) Isothermal-Isobaric (NPT)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,120)  mode
  120       format (i10)
         end if
         if (mode.eq.2 .or. mode.eq.4) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=130,end=130)  kelvin
  130       continue
            if (kelvin .le. 0.0d0) then
               write (iout,140)
  140          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,150)  kelvin
  150          format (f20.0)
            end if
            if (kelvin .le. 0.0d0)  kelvin = 298.0d0
            kelvin0 = kelvin
         end if
         if (mode.eq.3 .or. mode.eq.4) then
            isobaric = .true.
            atmsph = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=160,end=160)  atmsph
  160       continue
            if (atmsph .le. 0.0d0) then
               write (iout,170)
  170          format (/,' Enter the Desired Pressure in Atm',
     &                    ' [1.0] :  ',$)
               read (input,180)  atmsph
  180          format (f20.0)
            end if
            if (atmsph .le. 0.0d0)  atmsph = 1.0d0
         end if
      end if
c
c     get desired temperature for nonperiodic system
c
      if (.not. use_bounds) then
         isothermal = .true.
         kelvin = -1.0d0
         call nextarg (string,exist)
         if (exist)  read (string,*,err=190,end=190)  kelvin
  190    continue
         if (kelvin .le. 0.0d0) then
            write (iout,200)
  200       format (/,' Enter the Desired Temperature in Degrees',
     &                 ' K [298] :  ',$)
            read (input,210)  kelvin
  210       format (f20.0)
         end if
         if (kelvin .le. 0.0d0)  kelvin = 298.0d0
      end if
c
c     initialize any rattle constraints and setup dynamics
c
      call shakeup
      call mdinit
c
c     print out a header line for the dynamics computation
c
      if (integrate .eq. 'VERLET') then
         write (iout,220)
  220    format (/,' Molecular Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'STOCHASTIC') then
         write (iout,230)
  230    format (/,' Stochastic Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'RIGIDBODY') then
         write (iout,240)
  240    format (/,' Molecular Dynamics Trajectory via',
     &              ' Rigid Body Algorithm')
      else
         write (iout,250)
  250    format (/,' Molecular Dynamics Trajectory via',
     &              ' Modified Beeman Algorithm')
      end if
c
c     integrate equations of motion to take a time step
c
      do istep = 1, nstep
         if (integrate .eq. 'VERLET') then
            call verlet (istep,dt)
         else if (integrate .eq. 'STOCHASTIC') then
            call sdstep (istep,dt)
         else if (integrate .eq. 'RIGIDBODY') then
            call rgdstep (istep,dt)
         else
            call beeman (istep,dt)
         end if
c
c     remove center of mass translation and rotation if needed
c
         modstep = mod(istep,iprint)
         if (modstep.eq.0 .and. nuse.eq.n) then
            if (integrate.ne.'STOCHASTIC' .and.
     &          thermostat.ne.'ANDERSEN')  call mdrest
         end if
      end do
c
c     perform any final tasks before program exit
c
      call final
      end
