 
C     $Id: monte.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ################################################################
c     ##                   COPYRIGHT (C) 2001 by                    ##
c     ##  Michael Schnieders, Alan Grossfield & Jay William Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program monte  --  Monte Carlo/MCM conformational search  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "monte" performs a Monte Carlo/MCM conformational search
c     using either Cartesian single atom or torsional move sets
c
c
      program monte
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'iounit.i'
      include 'omega.i'
      include 'units.i'
      include 'zcoord.i'
      integer i,k,next,keep
      integer istep,nstep
      integer imin,freeunit
      real*8 global,ratio
      real*8 minimum,pminimum
      real*8 oldtors,tsize
      real*8 beta,eps
      real*8 size,grdmin
      real*8 boltz,random
      real*8 vector(3)
      real*8 xi(maxatm)
      real*8 yi(maxatm)
      real*8 zi(maxatm)
      character*1 answer
      character*120 minfile
      character*120 record
      character*120 string
      logical exist,collide
      logical torsmove
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     choose either the torsional or single atom move set
c
      torsmove = .false.
      call nextarg (answer, exist)
      if (.not. exist) then
         write (iout, 10)
   10    format (/,' Use [S]ingle Atom or [T]orsional Moves [S] : ',$)
         read (input, 20)  record
   20    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'T')  torsmove = .true.
c
c     if necessary, generate the internal coordinates
c
      if (torsmove) then
         call makeint (0)
         call initrot
      end if
c
c     get the desired Cartesian or torsional stepsize
c
      size = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=30,end=30)  size
   30 continue
      if (size .lt. 0.0d0) then
         if (torsmove) then
            write (iout,40)
   40       format (/,' Enter Maximum Step in Degrees [90] : ', $)
         else
            write (iout,50)
   50       format (/,' Enter Maximum Step in Angstroms [3.0] : ', $)
         end if
         read (input, 60)  string
   60    format (a120)
         read (string,*,err=70,end=70)  size
   70    continue
         if (size .lt. 0.0d0) then
            if (torsmove) then
               size = 90.0d0
            else
               size = 3.0d0
            end if
         end if
      end if
c
c     get the desired number of Monte Carlo steps
c
      nstep = -1
      call nextarg (string,exist)
      if (exist)  read (string,*,err=80,end=80)  nstep
   80 continue
      if (nstep .le. 0) then
         write (iout,90)
   90    format (/,' Number of Monte Carlo Steps [1000] : ', $)
         read (input,100)  nstep
  100    format (i15)
         if (nstep .le. 0)  nstep = 1000
      end if
c
c     get the simulation temperature for Metropolis criterion
c
      beta = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=110,end=110)  beta
  110 continue
      if (beta .lt. 0.0d0) then
         write (iout,120)
  120    format (/,' Enter the Temperature in Degrees K [300] : ', $)
         read (input, 130)  string
  130    format (a120)
         read (string,*,err=140,end=140)  beta
  140    continue
         if (beta .lt. 0.0d0)  beta = 300.0d0
      end if
      beta = 1.0d0 / (gasconst*beta)
c
c     get the gradient convergence for MCM protocol
c
      grdmin = -1.0d0
      call nextarg (string,exist)
      if (exist)  read(string,*,err=110,end=110)  grdmin
  150 continue
      if (grdmin .lt. 0.0d0) then
         write (iout,160)
  160    format (/,' Enter RMS Gradient Criterion [0.01] : ', $)
         read (input, 170)  string
  170    format (a120)
         read (string,*,err=180,end=180)  grdmin
  180    continue
         if (grdmin .lt. 0.0d0)  grdmin = 0.01
      end if
c
c     save the coordinates, then perform a minimization
c
      imin = freeunit ()
      minfile = filename(1:leng)//'.xyz'
      call version (minfile,'new')
      open (unit=imin,file=minfile,status='new')
      call prtxyz (imin)
      close (unit=imin)
      call mcmstep (minimum,grdmin)
      global = minimum
      write (iout,190)  minimum
  190 format (/,' Energy of Initial Structure :',11x,f12.4,
     &           ' Kcal/mole',/)
c
c     beginning of the Monte Carlo minimization protocol
c
c     the procedure is to store the current coordinates as the
c     reference, then try a move and store the new coordinates in
c     (xi,yi,zi), then minimize and apply the Metropolis criterion;
c     if accept, go to the coordinates before the minimization,
c     if reject, go revert back to the reference coordinates
c
      keep = 0
      do istep = 1, nstep
         call makeref
c
c     generate a random torsional angle move
c
         if (torsmove) then
            collide = .false.
            dowhile (.not. collide)
               k = int(nomega * random()) + 1
               k = zline(k)
               oldtors = ztors(k)
               tsize = 2.0d0 * size * (random()-0.5d0)
               ztors(k) = ztors(k) + tsize
               if (ztors(k) .gt. 180.0d0) then
                  ztors(k) = ztors(k) - 360.0d0
               else if (ztors(k) .lt. -180.0d0) then
                  ztors(k) = ztors(k) + 360.0d0
               end if
               call makexyz
               call chkclash (collide)
            end do
c
c     generate a random single atom Cartesian move
c
         else
            call ranvec (vector)
            k = int(n * random()) + 1
            x(k) = x(k) + size*vector(1)
            y(k) = y(k) + size*vector(2)
            z(k) = z(k) + size*vector(3)
         end if
c
c     store the coordinates, then perform a minimization
c
         do i = 1, n
            xi(i) = x(i)
            yi(i) = y(i)
            zi(i) = z(i)
         end do
         pminimum = minimum
         call mcmstep (minimum,grdmin)
         collide = .false.
         call chkclash (collide)
c
c     if the energy is lower, accept the current move
c
         if (minimum.le.pminimum .and. .not.collide) then
            keep = keep + 1
            pminimum = minimum
c
c     save the coordinates if this is the best minimum found
c
            if (minimum .lt. global) then
               global = minimum
               imin = freeunit ()
               minfile = filename(1:leng)//'.xyz'
               call version (minfile,'old')
               open (unit=imin,file=minfile,status='old')
               call prtxyz (imin)
               close (unit=imin)
c
c     restore the coordinates to those before the minimization
c
               do i = 1, n
                  x(i) = xi(i)
                  y(i) = yi(i)
                  z(i) = zi(i)
               end do
               if (torsmove) then
                  call makeint (2)
               end if
            end if
c
c     if the energy is higher, apply the Metropolis criterion
c
         else
            boltz = exp(beta*(pminimum-minimum))
            eps = random ()
c
c     reject the step and reset to the prior coordinates
c
            if (boltz.ge.eps .or. collide) then
               if (torsmove) then
                  ztors(k) = oldtors
                  call makexyz
               else
                  call getref
               end if
c
c     accept the step and set to coordinates before minimization
c
            else
               keep = keep + 1
               pminimum = minimum
               do i = 1, n
                  x(i) = xi(i)
                  y(i) = yi(i)
                  z(i) = zi(i)
               end do
               if (torsmove) then
                  call makeint (2)
               end if
            end if
         end if
c
c     print out some information about the current step
c
         ratio = dble(keep) / dble(istep)
         write (iout, 200)  istep,ratio,minimum,global
  200    format (' Step:',i8,'  Accept:',f8.3,
     &              '  Current:',f12.4,'  Global:',f12.4)
      end do
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function mcmstep  --  minimization phase of an MCM step  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mcmstep" implements the minimization phase of an MCM step
c     via Cartesian minimization following a Monte Carlo step
c
c
      subroutine mcmstep (minimum,grdmin)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'inform.i'
      include 'output.i'
      include 'usage.i'
      integer i,nvar
      real*8 xx(maxvar)
      real*8 mcm1,minimum,grdmin
      character*6 mode,method
      external mcm1,mcm2,optsave
c
c
c     prepare for the truncated Newton minimization
c
      mode = 'AUTO'
      method = 'AUTO'
      verbose = .false.
      iprint = 0
      iwrite = 0
      coordtype = 'NONE'
c
c     translate the coordinates of each active atom
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            xx(nvar) = x(i)
            nvar = nvar + 1
            xx(nvar) = y(i)
            nvar = nvar + 1
            xx(nvar) = z(i)
         end if
      end do
c
c     make the call to the optimization routine
c
      call tncg (mode,method,nvar,xx,minimum,grdmin,
     &                  mcm1,mcm2,optsave)
c
c     untranslate the final coordinates for active atoms
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar)
         end if
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  function mcm1  --  energy and gradient for MCM search  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "mcm1" is a service routine that computes the energy
c     and gradient for truncated Newton optimization in Cartesian
c     coordinate space
c
c
      function mcm1 (xx,g)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      integer i,nvar
      real*8 mcm1,e
      real*8 xx(maxvar)
      real*8 g(maxvar)
      real*8 derivs(3,maxatm)
c
c
c     translate optimization parameters to atomic coordinates
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar)
         end if
      end do
c
c     compute and store the energy and gradient
c
      call gradient (e,derivs)
      mcm1 = e
c
c     store Cartesian gradient as optimization gradient
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            g(nvar) = derivs(1,i)
            nvar = nvar + 1
            g(nvar) = derivs(2,i)
            nvar = nvar + 1
            g(nvar) = derivs(3,i)
         end if
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine mcm2  --  Hessian values for MCM search  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "mcm2" is a service routine that computes the sparse
c     matrix Hessian elements for truncated Newton optimization
c     in Cartesian coordinate space
c
c
      subroutine mcm2 (mode,xx,h,hinit,hstop,hindex,hdiag)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      integer i,j,k,nvar
      integer hinit(maxvar)
      integer hstop(maxvar)
      integer hvar(maxvar)
      integer huse(maxvar)
      integer hindex(maxhess)
      real*8 xx(maxvar)
      real*8 hdiag(maxvar)
      real*8 h(maxhess)
      character*4 mode
c
c
c     translate optimization parameters to atomic coordinates
c
      if (mode .eq. 'NONE')  return
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar)
         end if
      end do
c
c     compute and store the Hessian elements
c
      call hessian (h,hinit,hstop,hindex,hdiag)
c
c     transform the sparse Hessian to use only active atoms
c
      nvar = 0
      if (nuse .ne. n) then
         do i = 1, n
            k = 3 * (i-1)
            if (use(i)) then
               do j = 1, 3
                  nvar = nvar + 1
                  hvar(nvar) = j + k
                  huse(j+k) = nvar
               end do
            else
               do j = 1, 3
                  huse(j+k) = 0
               end do
            end if
         end do
         do i = 1, nvar
            k = hvar(i)
            hinit(i) = hinit(k)
            hstop(i) = hstop(k)
            hdiag(i) = hdiag(k)
            do j = hinit(i), hstop(i)
               hindex(j) = huse(hindex(j))
            end do
         end do
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine chkclash  --  check for near atom collisions  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "chkclash" determines if there are any atom clashes which
c     might cause trouble on subsequent energy evaluation
c
c
      subroutine chkclash (collide)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer i,j
      real*8 xi,yi,zi
      real*8 r2,eps
      logical collide
c
c
c     loop over atom pairs testing for near collisions
c
      eps = 0.01d0
      collide = .false.
      do i = 1, n-1
         xi = x(i)
         yi = y(i)
         zi = z(i)
         do j = i+1, n
            r2 = (x(j)-xi)**2 + (y(j)-yi)**2 + (z(j)-zi)**2
            if (r2 .le. eps) then
               collide = .true.
               goto 10
            end if
         end do
      end do
   10 continue
      return
      end
