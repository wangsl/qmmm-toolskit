 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program anneal  --  molecular dynamics simulated annealing  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "anneal" performs a simulated annealing protocol by means of
c     variable temperature molecular dynamics using either linear,
c     exponential or sigmoidal cooling schedules
c
c
      program anneal
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bond.i'
      include 'bound.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mdstuf.i'
      include 'potent.i'
      include 'solute.i'
      include 'usage.i'
      include 'warp.i'
      integer i,nequil
      integer nstep,istep
      integer next,modstep
      real*8 logmass,factor
      real*8 ratio,sigmoid
      real*8 dt,dtdump
      real*8 hot,cold
      real*8 fuzzy,sharp
      real*8 loose,tight
      logical exist
      character*1 answer
      character*8 cooltyp
      character*120 record
      character*120 string
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     get choice of statistical ensemble for periodic system
c
      hot = -1.0d0
      cold = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  hot
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  cold
   10 continue
      if (hot.lt.0.0d0 .or. cold.lt.0.0d0) then
         hot = -1.0d0
         cold = -1.0d0
         write (iout,20)
   20    format (/,' Enter the Initial and Final Temperatures in',
     &              ' Degrees K [1000,0] :  ',$)
         read (input,30)  record
   30    format (a120)
         read (record,*,err=40,end=40)  hot,cold
   40    continue
      end if
      if (hot .le. 0.0d0)  hot = 1000.0d0
      if (cold .le. 0.0d0)  cold = 0.0d0
c
c     set the number of steps of initial equilibration
c
      nequil = -1
      call nextarg (string,exist)
      if (exist)  read (string,*,err=50,end=50)  nequil
   50 continue
      if (nequil .lt. 0) then
         write (iout,60)
   60    format (/,' Enter the Number of Equilibration Steps [0] :  ',$)
         read (input,70)  nequil
   70    format (i10)
      end if
      if (nequil .le. 0)  nequil = 0
c
c     set the number of dynamics steps for the cooling protocol
c
      nstep = -1
      call nextarg (string,exist)
      if (exist)  read (string,*,err=80,end=80)  nstep
   80 continue
      if (nstep .lt. 0) then
         write (iout,90)
   90    format (/,' Enter the Number of Cooling Protocol Steps',
     &              ' [2000] :  ',$)
         read (input,100)  nstep
  100    format (i10)
      end if
      if (nstep .le. 0)  nstep = 2000
c
c     decide which annealing cooling protocol to use
c
      cooltyp = 'LINEAR'
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,110)
  110    format (/,' Use Linear, Sigmoidal or Exponential Cooling',
     &              ' Protocol ([L], S or E) :  ',$)
         read (input,120)  record
  120    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'S')  cooltyp = 'SIGMOID'
      if (answer .eq. 'E')  cooltyp = 'EXPONENT'
c
c     get the length of the dynamics time step in picoseconds
c
      dt = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=130,end=130)  dt
  130 continue
      if (dt .le. 0.0d0) then
         write (iout,140)
  140    format (/,' Enter the Time Step Length in Femtoseconds',
     &              ' [1.0] :  ',$)
         read (input,150)  dt
  150    format (f20.0)
      end if
      if (dt .le. 0.0d0)  dt = 1.0d0
      dt = 1.0d-3 * dt
c
c     set the time between trajectory snapshot coordinate dumps
c
      dtdump = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=160,end=160)  dtdump
  160 continue
      if (dtdump .lt. 0.0d0) then
         write (iout,170)
  170    format (/,' Enter Time between Dumps in Picoseconds',
     &              ' [0.1] :  ',$)
         read (input,180)  dtdump
  180    format (f20.0)
      end if
      if (dtdump .le. 0.0d0)  dtdump = 0.1d0
      iwrite = nint(dtdump/dt)
c
c     get factor by which atomic weights are to be increased
c
      logmass = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=190,end=190)  logmass
  190 continue
      if (logmass .lt. 0.0d0) then
         write (iout,200)
  200    format (/,' Increase Atomic Weights by a Factor of',
     &              ' 10^x [x=0.0] :  ',$)
         read (input,210)  logmass
  210    format (f20.0)
      end if
      if (logmass .le. 0.0d0) then
         logmass = 0.0d0
      else
         factor = 10.0d0**(logmass)
         do i = 1, n
            mass(i) = mass(i) * factor
         end do
      end if
c
c     rate of deformation change for potential surface smoothing
c
      if (use_smooth) then
         sharp = -1.0d0
         call nextarg (string,exist)
         if (exist)  read (string,*,err=220,end=220)  sharp
  220    continue
         if (sharp .le. 0.0d0) then
            write (iout,230)
  230       format (/,' Enter Final Desired Deformation Parameter',
     &                 ' [0.0] :  ',$)
            read (input,240)  sharp
  240       format (f20.0)
         end if
         if (sharp .le. 0.0d0)  sharp = 0.0d0
         fuzzy = deform - sharp
         if (fuzzy .le. 0.0d0)  fuzzy = 0.0d0
      end if
c
c     set values for temperature, pressure and coupling baths
c
      isothermal = .true.
      isobaric = .false.
      loose = 100.0d0 * dt
      tight = 10.0d0 * dt
      kelvin = hot
      kelvin0 = kelvin
      tautemp = loose
c
c     initialize any rattle constraints and setup dynamics
c
      call shakeup
      call mdinit
c
c     print out a header lines for the equilibration phase
c
      if (nequil .ne. 0) then
         write (iout,250)
  250    format (/,' Simulated Annealing Equilibration Phase')
         write (iout,260)  nequil,dt,logmass,hot,hot
  260    format (/,' Steps:',i6,3x,'Time/Step:',f6.3,' ps',3x,
     &              'LogMass:',f5.2,3x,'Temp:',f7.1,' to',f7.1)
      end if
c
c     take the dynamics steps for the equilibration phase
c
      do istep = 1, nequil
         if (integrate .eq. 'VERLET') then
            call verlet (istep,dt)
         else if (integrate .eq. 'STOCHASTIC') then
            call sdstep (istep,dt)
         else if (integrate .eq. 'RIGIDBODY') then
            call rgdstep (istep,dt)
         else
            call beeman (istep,dt)
         end if
         modstep = mod(istep,iprint)
         if (nuse.eq.n .and. modstep.eq.0)  call mdrest
      end do
c
c     start the cooling phase from the end of equilibration phase
c
      if (nequil .ne. 0)  call mdinit
c
c     print out a header lines for the cooling protocol
c
      write (iout,270)
  270 format (/,' Simulated Annealing Cooling Protocol')
      write (iout,280)  nstep,dt,logmass,hot,cold
  280 format (/,' Steps:',i6,3x,'Time/Step:',f6.3,' ps',3x,
     &           'LogMass:',f5.2,3x,'Temp:',f7.1,' to',f7.1)
c
c     set target temperature using the desired cooling protocol
c
      do istep = 1, nstep
         ratio = dble(istep) / dble(nstep)
         if (cooltyp .eq. 'SIGMOID') then
            ratio = sigmoid (3.5d0,ratio)
         else if (cooltyp .eq. 'EXPONENT') then
            ratio = 1.0d0 - exp(-5.0d0*ratio)
         end if
         kelvin = hot*(1.0d0-ratio) + cold*ratio
         kelvin0 = kelvin
         tautemp = loose*(1.0d0-ratio) + tight*ratio
c
c     set the deformation value if potential smoothing is used
c
         if (use_smooth) then
            ratio = (1.0d0-dble(istep)/dble(nstep))**3
            deform = sharp + ratio*fuzzy
         end if
c
c     integrate equations of motion to take a time step
c
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
