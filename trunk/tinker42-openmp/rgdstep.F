 
C     $Id: rgdstep.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###########################################################
c     ##                 COPYRIGHT (C) 2001 by                 ##
c     ##  Marina A. Vorobieva, Andrey Kutepov & Jay W. Ponder  ##
c     ##                  All Rights Reserved                  ##
c     ###########################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine rgdstep  --  rigid body molecular dynamics step  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "rgdstep" performs a single molecular dynamics time step
c     for a rigid body calculation
c
c
      subroutine rgdstep (istep,dt)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'group.i'
      include 'rgddyn.i'
      include 'units.i'
      include 'virial.i'
      integer i,j,k
      integer istep,size
      integer start,stop
      integer iter,maxiter
      real*8 dt,etot,eps
      real*8 eksum,epot
      real*8 temp,pres
      real*8 delta,weigh
      real*8 fx,fy,fz
      real*8 xr,yr,zr
      real*8 x2,y2,z2
      real*8 xcm,ycm,zcm
      real*8 xcp,ycp,zcp
      real*8 fc(3),tc(3)
      real*8 dfi(3),dfip(3)
      real*8 wcp(3),inert(6)
      real*8 arot(3,3)
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 xm(maxatm),xp(maxatm)
      real*8 ym(maxatm),yp(maxatm)
      real*8 zm(maxatm),zp(maxatm)
      real*8 derivs(3,maxatm)
c
c
c     set iteration limit and convergence for angular velocity 
c
      maxiter = 15
      eps = 1.0d-12
c
c     get the energy and atomic forces prior to the step
c
      call gradient (epot,derivs)
c
c     get the extent and size of the current rigid body
c
      do i = 1, ngrp
         start = igrp(1,i)
         stop = igrp(2,i)
         size = stop - start + 1
c
c     compute the center of mass coordinates prior to step
c
         xcm = 0.0d0
         ycm = 0.0d0
         zcm = 0.0d0
         do j = start, stop
            k = kgrp(j)
            weigh = mass(k)
            xcm = xcm + x(k)*weigh
            ycm = ycm + y(k)*weigh
            zcm = zcm + z(k)*weigh
         end do
         xcm = xcm / grpmass(i)
         ycm = ycm / grpmass(i)
         zcm = zcm / grpmass(i)
         do j = start, stop
            k = kgrp(j)
            xm(k) = x(k) - xcm
            ym(k) = y(k) - ycm
            zm(k) = z(k) - zcm
         end do
c
c     make center of mass correction to virial for rigid body
c
         do j = start, stop
            k = kgrp(j)
            vir(1,1) = vir(1,1) - xm(k)*derivs(1,k)
            vir(2,1) = vir(2,1) - ym(k)*derivs(1,k)
            vir(3,1) = vir(3,1) - zm(k)*derivs(1,k)
            vir(1,2) = vir(1,2) - xm(k)*derivs(2,k)
            vir(2,2) = vir(2,2) - ym(k)*derivs(2,k)
            vir(3,2) = vir(3,2) - zm(k)*derivs(2,k)
            vir(1,3) = vir(1,3) - xm(k)*derivs(3,k)
            vir(2,3) = vir(2,3) - ym(k)*derivs(3,k)
            vir(3,3) = vir(3,3) - zm(k)*derivs(3,k)
         end do
c
c     compute the force and torque components for rigid body
c
         do j = 1, 3
            fc(j) = 0.0d0
            tc(j) = 0.0d0
         end do
         do j = start, stop
            k = kgrp(j)
            xr = xm(k)
            yr = ym(k)
            zr = zm(k)
            fx = -convert * derivs(1,k)
            fy = -convert * derivs(2,k)
            fz = -convert * derivs(3,k)
            fc(1) = fc(1) + fx
            fc(2) = fc(2) + fy
            fc(3) = fc(3) + fz
            tc(1) = tc(1) + yr*fz - zr*fy
            tc(2) = tc(2) + zr*fx - xr*fz
            tc(3) = tc(3) + xr*fy - yr*fx
         end do
c
c     update the translational velocities and center of mass
c
         do j = 1, 3
            vcm(j,i) = vcm(j,i) + dt*fc(j)/grpmass(i)
         end do
         xcp = xcm + dt*vcm(1,i)
         ycp = ycm + dt*vcm(2,i)
         zcp = zcm + dt*vcm(3,i)
c
c     single atom groups are treated as a separate case
c
         if (size .eq. 1) then
            k = kgrp(igrp(1,i))
            x(k) = xcm
            y(k) = ycm
            z(k) = zcm
            do j = 1, 3
               wcm(j,i) = 0.0d0
               lm(j,i) = 0.0d0
            end do
c
c     make the update to the rigid body angular momentum
c
         else
            do j = 1, 3
               lm(j,i) = lm(j,i) + dt*tc(j)
            end do
c
c     use iterative scheme to converge the angular velocity
c
            iter = 0
            delta = 1.0d0
            do j = 1, 3
               dfi(j) = dt * wcm(j,i)
            end do
c
c     first, make a prediction for the atomic coordinates
c
            dowhile (delta.gt.eps .and. iter.le.maxiter)
               iter = iter + 1
               call rotrgd (dfi,arot)
               do j = start, stop
                  k = kgrp(j)
                  xp(k) = arot(1,1)*xm(k) + arot(1,2)*ym(k)
     &                       + arot(1,3)*zm(k)
                  yp(k) = arot(2,1)*xm(k) + arot(2,2)*ym(k)
     &                       + arot(2,3)*zm(k)
                  zp(k) = arot(3,1)*xm(k) + arot(3,2)*ym(k)
     &                       + arot(3,3)*zm(k)
               end do
c
c     calculate the upper triangle of the inertia tensor
c
               do j = 1, 6
                  inert(j) = 0.0d0
               end do
               do j = start, stop
                  k = kgrp(j)
                  xr = xp(k)
                  yr = yp(k)
                  zr = zp(k)
                  x2 = xr * xr
                  y2 = yr * yr
                  z2 = zr * zr
                  weigh = mass(k)
                  inert(1) = inert(1) + weigh*(y2+z2)
                  inert(2) = inert(2) - weigh*xr*yr
                  inert(3) = inert(3) - weigh*xr*zr
                  inert(4) = inert(4) + weigh*(x2+z2)
                  inert(5) = inert(5) - weigh*yr*zr
                  inert(6) = inert(6) + weigh*(x2+y2)
               end do
c
c     compute the angular velocity from the relation L=Iw
c
               do j = 1, 3
                  wcp(j) = lm(j,i)
               end do
               if (linear(i)) then
                  call linbody (i,inert,wcp)
               else
                  call cholesky (3,inert,wcp)
               end if
c
c     check to see if angular distance change has converged
c
               delta = 0.0d0
               do j = 1, 3
                  dfip(j) = dt * wcp(j)
                  delta = delta + abs(dfip(j)-dfi(j))
                  dfi(j) = dfip(j)
               end do
            end do
c
c     set the final angular velocities and atomic coordinates
c
            do j = 1, 3
               wcm(j,i) = wcp(j)
            end do
            call rotrgd (dfi,arot)
            do j = start, stop
               k = kgrp(j)
               x(k) = arot(1,1)*xm(k) + arot(1,2)*ym(k)
     &                   + arot(1,3)*zm(k) + xcp
               y(k) = arot(2,1)*xm(k) + arot(2,2)*ym(k)
     &                   + arot(2,3)*zm(k) + ycp
               z(k) = arot(3,1)*xm(k) + arot(3,2)*ym(k)
     &                   + arot(3,3)*zm(k) + zcp
            end do
         end if
      end do
c
c     accumulate the kinetic energy and its outer product
c
      call kinetic (eksum,ekin)
c
c     compute and control the temperature and pressure
c
      call temper2 (dt,eksum,temp)
      call pressure (dt,ekin,pres,stress)
c
c     system energy is sum of kinetic and potential energies
c
      etot = eksum + epot
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine rotrgd  --  rigid dynamics rotation matrix  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "rotrgd" finds the rotation matrix for a rigid body due
c     to a single step of dynamics
c
c
      subroutine rotrgd (dfi,arot)
      implicit none
      real*8 x,xc,xs
      real*8 y,yc,ys
      real*8 z,zc,zs
      real*8 cosine,sine
      real*8 anorm,coterm
      real*8 dfi(3)
      real*8 arot(3,3)
c
c
c     construct rotation matrix from angular distance
c
      anorm = sqrt(dfi(1)**2 + dfi(2)**2 + dfi(3)**2)
      cosine = cos(anorm)
      sine = sin(anorm)
      coterm = 1.0d0 - cosine
      if (anorm .le. 0.0d0)  anorm = 1.0d0
      x = dfi(1) / anorm
      y = dfi(2) / anorm
      z = dfi(3) / anorm
      xc = x * coterm
      yc = y * coterm
      zc = z * coterm
      xs = x * sine
      ys = y * sine
      zs = z * sine
      arot(1,1) = xc*x + cosine
      arot(2,1) = xc*y + zs
      arot(3,1) = xc*z - ys
      arot(1,2) = yc*x - zs
      arot(2,2) = yc*y + cosine
      arot(3,2) = yc*z + xs
      arot(1,3) = zc*x + ys
      arot(2,3) = zc*y - xs
      arot(3,3) = zc*z + cosine
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine linbody  --  angular velocity of linear body  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "linbody" finds the angular velocity of a linear rigid body
c     given the inertia tensor and angular momentum
c
c
      subroutine linbody (i,inert,wcp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'group.i'
      integer i,j,k
      real*8 rinv,rmin
      real*8 a11,a12,a22
      real*8 b1,b2,w1,w2
      real*8 wcp(3),rmol(3)
      real*8 r1(3),r2(3),r3(3)
      real*8 inert(6)
c
c
c     construct a normalized vector along the molecular axis
c
      j = kgrp(igrp(1,i))
      k = kgrp(igrp(2,i))
      rmol(1) = x(k) - x(j)
      rmol(2) = y(k) - y(j)
      rmol(3) = z(k) - z(j)
      rinv = 1.0d0 / sqrt(rmol(1)**2+rmol(2)**2+rmol(3)**2)
      do j = 1, 3
         rmol(j) = rmol(j) * rinv
      end do
c
c     find two orthogonal vectors to complete coordinate frame
c
      k = 1
      rmin = abs(rmol(1))
      do j = 2, 3
         if (abs(rmol(j)) .lt. rmin) then
            k = j
            rmin = abs(rmol(j))
         end if
      end do
      do j = 1, 3
         r1(j) = -rmol(k) * rmol(j)
      end do
      r1(k) = 1.0d0 + r1(k)
      rinv = 1.0d0 / sqrt(r1(1)**2+r1(2)**2+r1(3)**2)
      do j = 1, 3
         r1(j) = r1(j) * rinv
      end do
      r2(1) = r1(2)*rmol(3) - r1(3)*rmol(2)
      r2(2) = r1(3)*rmol(1) - r1(1)*rmol(3)
      r2(3) = r1(1)*rmol(2) - r1(2)*rmol(1)
c
c     solve the 2-by-2 linear system for angular velocity
c
      r3(1) = inert(1)*r1(1) + inert(2)*r1(2) + inert(3)*r1(3)
      r3(2) = inert(2)*r1(1) + inert(4)*r1(2) + inert(5)*r1(3)
      r3(3) = inert(3)*r1(1) + inert(5)*r1(2) + inert(6)*r1(3)
      a11 = r1(1)*r3(1) + r1(2)*r3(2) + r1(3)*r3(3)
      r3(1) = inert(1)*r2(1) + inert(2)*r2(2) + inert(3)*r2(3)
      r3(2) = inert(2)*r2(1) + inert(4)*r2(2) + inert(5)*r2(3)
      r3(3) = inert(3)*r2(1) + inert(5)*r2(2) + inert(6)*r2(3)
      a12 = r1(1)*r3(1) + r1(2)*r3(2) + r1(3)*r3(3)
      a22 = r2(1)*r3(1) + r2(2)*r3(2) + r2(3)*r3(3)
      b1 = r1(1)*wcp(1) + r1(2)*wcp(2) + r1(3)*wcp(3)
      b2 = r2(1)*wcp(1) + r2(2)*wcp(2) + r2(3)*wcp(3)
      w1 = (a12*b2-a22*b1) / (a12*a12-a11*a22)
      w2 = (b2-a12*w1) / a22
      do j = 1, 3
         wcp(j) = w1*r1(j) + w2*r2(j)
      end do
      return
      end
