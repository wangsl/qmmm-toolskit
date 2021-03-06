 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine diffeq  --  differential equation integration  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "diffeq" performs the numerical integration of an ordinary
c     differential equation using an adaptive stepsize method to
c     solve the corresponding coupled first-order equations of the
c     general form dyi/dx = f(x,y1,...,yn) for yi = y1,...,yn
c
c     variables and parameters :
c
c     nvar      number of coupled first-order differential equations
c     y         contains the values of the dependent variables
c     x1        value of the beginning integration limit
c     x2        value of the ending integration limit
c     eps       relative accuracy required of the integration steps
c     h1        initial guess for the first integration stepsize
c     hmin      minimum allowed integration stepsize
c     nok       number of initially successful integration steps
c     nbad      number of integration steps that required retry
c
c     required external routines :
c
c     gvalue    subroutine to find the right-hand side of the
c                  first-order differential equations
c
c
      subroutine diffeq (nvar,y,x1,x2,eps,h1,hmin,nok,nbad,gvalue)
      include 'sizes.i'
      include 'iounit.i'
      integer maxgda,maxstep
      real*8 tiny
      parameter (maxgda=4*maxatm)
      parameter (maxstep=1000)
      parameter (tiny=1.0d-30)
      integer i,nvar,nstep
      integer nok,nbad
      real*8 x1,x2,eps,h1,hmin
      real*8 x,h,hdid,hnext
      real*8 y(maxgda)
      real*8 dydx(maxgda)
      real*8 yscal(maxgda)
      logical terminate
      character*7 status
      external gvalue
c
c
c     initialize starting limit, step size and status counters
c
      terminate = .false.
      x = x1
      h = sign(h1,x2-x1)
      nstep = 0
      nok = 0
      nbad = 0
c
c     perform a series of individual integration steps
c
      dowhile (.not. terminate)
         call gvalue (x,y,dydx)
         do i = 1, nvar
            yscal(i) = abs(y(i)) + abs(h*dydx(i)) + tiny
         end do
c
c     set the final step to stop at the integration limit
c
         if ((x+h-x2)*(x+h-x1) .gt. 0.0d0)  h = x2 - x
c
c     take a Bulirsch-Stoer integration step
c
         call bsstep (nvar,x,dydx,y,h,eps,yscal,hdid,hnext,gvalue)
c
c     mark the current step as either good or bad
c
         if (hdid .eq. h) then
            nok = nok + 1
            status = 'Success'
         else
            nbad = nbad + 1
            status = ' Retry '
         end if
c
c     update stepsize and get information about the current step
c
         h = hnext
         nstep = nstep + 1
         call gdastat (nstep,x,y,status)
c
c     test for convergence to the final integration limit
c
         if ((x-x2)*(x2-x1) .ge. 0.0d0) then
            write (iout,10)
   10       format (/,' DIFFEQ  --  Normal Termination',
     &                 ' at Integration Limit')
            terminate = .true.
         end if
c
c     test for a trial stepsize that is too small
c
         if (abs(hnext) .lt. hmin) then
            write (iout,20)
   20       format (/,' DIFFEQ  --  Incomplete Integration',
     &                 ' due to SmallStep')
            terminate = .true.
         end if
c
c     test for too many total integration steps
c
         if (nstep .ge. maxstep) then
            write (iout,30)
   30       format (/,' DIFFEQ  --  Incomplete Integration',
     &                 ' due to IterLimit')
            terminate = .true.
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine bsstep  --  Bulirsch-Stoer integration step  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "bsstep" takes a single Bulirsch-Stoer step with monitoring
c     of local truncation error to ensure accuracy
c
c     literature reference:
c
c     W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P.
c     Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge
c     University Press, 1992, Section 16.4
c
c
      subroutine bsstep (nv,x,dydx,y,htry,eps,yscal,hdid,hnext,gvalue)
      include 'sizes.i'
      include 'iounit.i'
      integer maxgda,kmaxx,imax
      real*8 safe1,safe2,redmax,redmin,tiny,scalmx
      parameter (maxgda=4*maxatm)
      parameter (kmaxx=8)
      parameter (imax=kmaxx+1)
      parameter (safe1=0.25d0)
      parameter (safe2=0.7d0)
      parameter (redmax=1.0d-5)
      parameter (redmin=0.7d0)
      parameter (tiny=1.0d-30)
      parameter (scalmx=0.1d0)
      integer i,iq,k,kk,km,kmax,kopt
      integer nv,nseq(imax)
      real*8 eps,hdid,hnext,htry,x
      real*8 eps1,epsold,errmax,fact,h,red
      real*8 scale,work,wrkmin,xest,xnew
      real*8 dydx(maxgda),y(maxgda),yscal(maxgda)
      real*8 a(imax),alf(kmaxx,kmaxx),err(kmaxx)
      real*8 yerr(maxgda),ysav(maxgda),yseq(maxgda)
      logical first,reduct
      save a,alf,epsold,first,kmax,kopt,nseq,xnew
      external gvalue
      data first  / .true. /
      data epsold / -1.0d0 /
      data nseq   / 2,4,6,8,10,12,14,16,18 /
c
c
      if (eps .ne. epsold) then
         hnext = -1.0d29
         xnew = -1.0d29
         eps1 = safe1 * eps
         a(1) = 1.0d0 + dble(nseq(1))
         do k = 1, kmaxx
            a(k+1) = a(k) + dble(nseq(k+1))
         end do
         do iq = 2, kmaxx
            do k = 1, iq-1
               alf(k,iq) = eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.0d0)
     &                                                 *(2*k+1)))
            end do
         end do
         epsold = eps
         do kopt = 2, kmaxx-1
            if (a(kopt+1) .gt. a(kopt)*alf(kopt-1,kopt))  goto 10
         end do
   10    continue
         kmax = kopt
      end if
      h = htry
      do i = 1, nv
         ysav(i) = y(i)
      end do
      if (h.ne.hnext .or. x.ne.xnew) then
         first = .true.
         kopt = kmax
      end if
      reduct = .false.
   20 continue
      do k = 1, kmax
         xnew = x + h
         if (xnew .eq. x) then
            write (iout,30)
   30       format (' BSSTEP  --  Underflow of Step Size')
            call fatal
         end if
         call mmid (nseq(k),h,nv,x,dydx,ysav,yseq,gvalue)
         xest = (h/dble(nseq(k)))**2
         call pzextr (k,nv,xest,yseq,y,yerr)
         if (k .ne. 1) then
            errmax = tiny
            do i = 1, nv
               errmax = max(errmax,abs(yerr(i)/yscal(i)))
            end do
            errmax = errmax / eps
            km = k - 1
            err(km) = (errmax/safe1)**(1.0d0/(2*km+1))
         end if
         if (k.ne.1 .and. (k.ge.kopt-1 .or. first)) then
            if (errmax .lt. 1.0d0)  goto 50
            if (k.eq.kmax .or. k.eq.kopt+1) then
               red = safe2 / err(km)
               goto 40
            else if (k .eq. kopt) then
               if (alf(kopt-1,kopt) .lt. err(km)) then
                  red = 1.0d0 / err(km)
                  goto 40
               end if
            else if (kopt .eq. kmax)then
               if (alf(km,kmax-1) .lt. err(km)) then
                  red = alf(km,kmax-1) * safe2 / err(km)
                  goto 40
               end if
            else if (alf(km,kopt) .lt. err(km)) then
               red = alf(km,kopt-1) / err(km)
               goto 40
            end if
         end if
      end do
   40 continue
      red = min(red,redmin)
      red = max(red,redmax)
      h = h * red
      reduct = .true.
      goto 20
   50 continue
      x = xnew
      hdid = h
      first = .false.
      wrkmin = 1.0d35
      do kk = 1, km
         fact = max(err(kk),scalmx)
         work = fact * a(kk+1)
         if (work .lt. wrkmin) then
            scale = fact
            wrkmin = work
            kopt = kk + 1
         end if
      end do
      hnext = h / scale
      if (kopt.ge.k .and. kopt.ne.kmax .and. .not.reduct) then
         fact = max(scale/alf(kopt-1,kopt),scalmx)
         if (a(kopt+1)*fact .le. wrkmin) then
            hnext = h / fact
            kopt = kopt + 1
         end if
      end if
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine mmid  --  takes a modified midpoint step  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "mmid" implements a modified midpoint method to advance the
c     integration of a set of first order differential equations
c
c
      subroutine mmid (nstep,htot,nvar,xs,dydx,y,yout,gvalue)
      include 'sizes.i'
      integer maxgda
      parameter (maxgda=4*maxatm)
      integer i,k
      integer nstep,nvar
      real*8 htot,h,h2
      real*8 xs,x,temp
      real*8 y(maxgda)
      real*8 yout(maxgda)
      real*8 dydx(maxgda)
      real*8 ym(maxgda)
      real*8 yn(maxgda)
      external gvalue
c
c
c     set substep size based on number of steps to be taken
c
      h = htot / dble(nstep)
      h2 = 2.0d0 * h
c
c     take the first substep and get values at ends of step
c
      do i = 1, nvar
         ym(i) = y(i)
         yn(i) = y(i) + h*dydx(i)
      end do
      x = xs + h
      call gvalue (x,yn,yout)
c
c     take the second and subsequent substeps
c
      do k = 2, nstep
         do i = 1, nvar
            temp = ym(i) + h2*yout(i)
            ym(i) = yn(i)
            yn(i) = temp
         end do
         x = x + h
         call gvalue (x,yn,yout)
      end do
c
c     complete the update of values for the last substep
c
      do i = 1, nvar
         yout(i) = 0.5d0 * (ym(i)+yn(i)+h*yout(i))
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine pzextr  --  polynomial extrapolation method  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "pzextr" is a polynomial extrapolation routine used during
c     Bulirsch-Stoer integration of ordinary differential equations
c
c
      subroutine pzextr (iest,nvar,xest,yest,yz,dy)
      include 'sizes.i'
      integer maxgda,imax
      parameter (maxgda=4*maxatm)
      parameter (imax=13)
      integer i,j,iest,nvar
      real*8 xest,delta,f1,f2,q
      real*8 yz(maxgda),dy(maxgda)
      real*8 yest(maxgda),d(maxgda)
      real*8 x(imax),qcol(maxgda,imax)
      save x,qcol
c
c
      x(iest) = xest
      do j = 1, nvar
         dy(j) = yest(j)
         yz(j) = yest(j)
      end do
      if (iest .eq. 1) then
         do j = 1, nvar
            qcol(j,1) = yest(j)
         end do
      else
         do j = 1, nvar
            d(j) = yest(j)
         end do
         do i = 1, iest-1
            delta = 1.0d0 / (x(iest-i)-xest)
            f1 = xest * delta
            f2 = x(iest-i) * delta
            do j = 1, nvar
               q = qcol(j,i)
               qcol(j,i) = dy(j)
               delta = d(j) - q
               dy(j) = f1 * delta
               d(j) = f2 * delta
               yz(j) = yz(j) + dy(j)
            end do
         end do
         do j = 1, nvar
            qcol(j,iest) = dy(j)
         end do
      end if
      return
      end
