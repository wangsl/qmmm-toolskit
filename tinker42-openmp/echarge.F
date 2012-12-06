 
C     $Id: echarge.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge  --  charge-charge potential energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge" calculates the charge-charge interaction energy
c
c
      subroutine echarge
      implicit none
      include 'sizes.i'
      include 'cutoff.i'
      include 'warp.i'
c
c
c     choose pairwise double loop, method of lights, potential
c     smoothing or Ewald summation version
c
      if (use_ewald .and. use_lights) then
         call echarge0e
      else if (use_ewald) then
         call echarge0d
      else if (use_smooth) then
         call echarge0c
      else if (use_lights) then
         call echarge0b
      else
         call echarge0a
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine echarge0a  --  charge energy via double loop  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "echarge0a" calculates the charge-charge interaction energy
c     using a pairwise double loop
c
c
      subroutine echarge0a
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k
      integer ii,kk
      integer in,kn
      integer ic,kc
      real*8 e,r,r2,fgrp
      real*8 f,fi,fik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xc,yc,zc
      real*8 xic,yic,zic
      real*8 shift,taper,trans
      real*8 rc,rc2,rc3,rc4
      real*8 rc5,rc6,rc7
      real*8 cscale(maxatm)
      logical proceed,iuse
c
c
c     zero out the charge interaction energy
c
      ec = 0.0d0
      if (nion .eq. 0)  return
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('CHARGE')
c
c     calculate the charge interaction energy term
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(ii)
         ic = kion(ii)
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xi = x(i) - xic
         yi = y(i) - yic
         zi = z(i) - zic
         fi = f * pchg(ii)
         iuse = (use(i) .or. use(ic))
         do j = 1, nion
            cscale(iion(j)) = 1.0d0
         end do
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(kk)
            kc = kion(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kc))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xc = xic - x(kc)
               yc = yic - y(kc)
               zc = zic - z(kc)
               if (use_image)  call image (xc,yc,zc,0)
               rc2 = xc*xc + yc*yc + zc*zc
               if (rc2 .le. off2) then
                  xr = xc + xi - x(k) + x(kc)
                  yr = yc + yi - y(k) + y(kc)
                  zr = zc + zi - z(k) + z(kc)
                  r2 = xr*xr + yr*yr + zr*zr
                  r = sqrt(r2)
                  fik = fi * pchg(kk) * cscale(kn)
                  e = fik / r
c
c     use shifted energy switching if near the cutoff distance
c
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  if (rc2 .gt. cut2) then
                     rc = sqrt(rc2)
                     rc3 = rc2 * rc
                     rc4 = rc2 * rc2
                     rc5 = rc2 * rc3
                     rc6 = rc3 * rc3
                     rc7 = rc3 * rc4
                     taper = c5*rc5 + c4*rc4 + c3*rc3
     &                          + c2*rc2 + c1*rc + c0
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     e = e * taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall charge-charge energy component
c
                  ec = ec + e
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(ii)
         ic = kion(ii)
         iuse = (use(i) .or. use(ic))
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xi = x(i) - xic
         yi = y(i) - yic
         zi = z(i) - zic
         fi = f * pchg(ii)
         do j = 1, nion
            cscale(iion(j)) = 1.0d0
         end do
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii, nion
            k = iion(kk)
            kn = jion(kk)
            kc = kion(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kc))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  xc = xic - x(kc)
                  yc = yic - y(kc)
                  zc = zic - z(kc)
                  call image (xc,yc,zc,j)
                  rc2 = xc*xc + yc*yc + zc*zc
                  if (rc2 .le. off2) then
                     xr = xc + xi - x(k) + x(kc)
                     yr = yc + yi - y(k) + y(kc)
                     zr = zc + zi - z(k) + z(kc)
                     r2 = xr*xr + yr*yr + zr*zr
                     r = sqrt(r2)
                     fik = fi * pchg(kk)
                     if (use_polymer) then
                        if (r2 .le. polycut2)  fik = fik * cscale(kn)
                     end if
                     e = fik / r
c
c     use shifted energy switching if near the cutoff distance
c
                     shift = fik / (0.5d0*(off+cut))
                     e = e - shift
                     if (rc2 .gt. cut2) then
                        rc = sqrt(rc2)
                        rc3 = rc2 * rc
                        rc4 = rc2 * rc2
                        rc5 = rc2 * rc3
                        rc6 = rc3 * rc3
                        rc7 = rc3 * rc4
                        taper = c5*rc5 + c4*rc4 + c3*rc3
     &                             + c2*rc2 + c1*rc + c0
                        trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                                  + f3*rc3 + f2*rc2 + f1*rc + f0)
                        e = e * taper + trans
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group)  e = e * fgrp
c
c     increment the overall charge-charge energy component
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ec = ec + e
                  end if
               end do
            end if
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge0b  --  method of lights charge energy  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge0b" calculates the charge-charge interaction energy
c     using the method of lights to locate neighboring atoms
c
c
      subroutine echarge0b
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'iounit.i'
      include 'light.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,ii,kk
      integer in,ic,kn,kc
      integer kgy,kgz
      integer start,stop
      integer map(maxlight)
      real*8 e,r,r2,fgrp
      real*8 f,fi,fik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xc,yc,zc
      real*8 xic,yic,zic
      real*8 shift,taper,trans
      real*8 rc,rc2,rc3,rc4
      real*8 rc5,rc6,rc7
      real*8 cscale(maxatm)
      real*8 xsort(maxlight)
      real*8 ysort(maxlight)
      real*8 zsort(maxlight)
      logical proceed,iuse,repeat
c
c
c     zero out the charge interaction energy
c
      ec = 0.0d0
      if (nion .eq. 0)  return
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('CHARGE')
c
c     transfer the interaction site coordinates to sorting arrays
c
      do i = 1, nion
         k = kion(i)
         xsort(i) = x(k)
         ysort(i) = y(k)
         zsort(i) = z(k)
      end do
c
c     use the method of lights to generate neighbors
c
      call lights (nion,map,xsort,ysort,zsort)
c
c     now, loop over all atoms computing the interactions
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(ii)
         ic = kion(ii)
         xic = xsort(rgx(ii))
         yic = ysort(rgy(ii))
         zic = zsort(rgz(ii))
         xi = x(i) - x(ic)
         yi = y(i) - y(ic)
         zi = z(i) - z(ic)
         fi = f * pchg(ii)
         iuse = (use(i) .or. use(ic))
         do j = 1, nion
            cscale(iion(j)) = 1.0d0
         end do
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            k = iion(map(kk))
            kn = jion(map(kk))
            kc = kion(map(kk))
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kc))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xc = xic - xsort(j)
               yc = yic - ysort(kgy)
               zc = zic - zsort(kgz)
               if (use_image) then
                  if (abs(xc) .gt. xcell2)  xc = xc - sign(xcell,xc)
                  if (abs(yc) .gt. ycell2)  yc = yc - sign(ycell,yc)
                  if (abs(zc) .gt. zcell2)  zc = zc - sign(zcell,zc)
                  if (monoclinic) then
                     xc = xc + zc*beta_cos
                     zc = zc * beta_sin
                  else if (triclinic) then
                     xc = xc + yc*gamma_cos + zc*beta_cos
                     yc = yc*gamma_sin + zc*beta_term
                     zc = zc * gamma_term
                  end if
               end if
               rc2 = xc*xc + yc*yc + zc*zc
               if (rc2 .le. off2) then
                  xr = xc + xi - x(k) + x(kc)
                  yr = yc + yi - y(k) + y(kc)
                  zr = zc + zi - z(k) + z(kc)
                  r2 = xr*xr + yr*yr + zr*zr
                  r = sqrt(r2)
                  fik = fi * pchg(map(kk)) * cscale(kn)
                  if (use_polymer) then
                     if (r2 .gt. polycut2)  fik = fi * pchg(map(kk))
                  end if
                  e = fik / r
c
c     use shifted energy switching if near the cutoff distance
c
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  if (rc2 .gt. cut2) then
                     rc = sqrt(rc2)
                     rc3 = rc2 * rc
                     rc4 = rc2 * rc2
                     rc5 = rc2 * rc3
                     rc6 = rc3 * rc3
                     rc7 = rc3 * rc4
                     taper = c5*rc5 + c4*rc4 + c3*rc3
     &                          + c2*rc2 + c1*rc + c0
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     e = e * taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall charge-charge energy component
c
                  ec = ec + e
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = nlight
            goto 10
         end if
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine echarge0c  --  charge energy for smoothing  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "echarge0c" calculates the charge-charge interaction energy
c     for use with potential smoothing methods
c
c
      subroutine echarge0c
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'units.i'
      include 'usage.i'
      include 'warp.i'
      integer i,j,k
      integer ii,kk
      integer in,kn
      real*8 e,r,r2,fgrp
      real*8 f,fi,fik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 erf,wterm,width
      real*8 width2,width3
      real*8 cscale(maxatm)
      logical proceed,iuse
      external erf
c
c
c     zero out the charge interaction energy
c
      ec = 0.0d0
      if (nion .eq. 0)  return
c
c     set the energy units conversion factor
c
      f = electric / dielec
c
c     set the extent of smoothing to be performed
c
      width = deform * diffc
      if (use_dem) then
         if (width .gt. 0.0d0)  width = 0.5d0 / sqrt(width)
      else if (use_gda) then
         wterm = sqrt(3.0d0/(2.0d0*diffc))
      end if
      width2 = width * width
      width3 = width * width2
c
c     calculate the charge interaction energy term
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(ii)
         iuse = (use(i))
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         do j = ii+1, nion
            cscale(iion(j)) = 1.0d0
         end do
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               r = sqrt(r2)
               fik = fi * pchg(kk) * cscale(kn)
               e = fik / r
c
c     transform the potential function via smoothing
c
               if (use_dem) then
                  if (width .gt. 0.0d0)  e = e * erf(width*r)
               else if (use_gda) then
                  width = m2(i) + m2(k)
                  if (width .gt. 0.0d0) then
                     width = wterm / sqrt(width)
                     e = e * erf(width*r)
                  end if
               else if (use_tophat .and. r.lt.width) then
                  if (width .gt. r)
     &               e = fik * (3.0d0*width2-r2) / (2.0d0*width3)
               else if (use_stophat) then
                  e = fik / (r+width)
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the overall charge-charge energy component
c
               ec = ec + e
            end if
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine echarge0d  --  Ewald summation charge energy  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "echarge0d" calculates the charge-charge interaction energy
c     using a particle mesh Ewald summation
c
c
      subroutine echarge0d
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'ewald.i'
      include 'group.i'
      include 'math.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k
      integer ii,kk
      integer in,kn
      real*8 e,fs,fgrp
      real*8 f,fi,fik
      real*8 r,r2,rew
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xd,yd,zd
      real*8 erfc,erfterm
      real*8 scale,scaleterm
      real*8 cscale(maxatm)
      logical proceed,iuse
      external erfc
c
c
c     zero out the Ewald charge interaction energy
c
      ec = 0.0d0
      if (nion .eq. 0)  return
c
c     set conversion factor, cutoff and scaling coefficients
c
      f = electric / dielec
      call switch ('EWALD')
c
c     compute the Ewald self-energy term over all the atoms
c
      fs = -f * aewald / sqrtpi
      do ii = 1, nion
         e = fs * pchg(ii)**2
         ec = ec + e
      end do
c
c     compute the reciprocal space part of the Ewald summation
c
      call epme
c
c     compute the cell dipole boundary correction term
c
      if (.not. tinfoil) then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do ii = 1, nion
            i = iion(ii)
            xd = xd + pchg(ii)*x(i)
            yd = yd + pchg(ii)*y(i)
            zd = zd + pchg(ii)*z(i)
         end do
         e = (2.0d0/3.0d0) * f * (pi/volbox) * (xd*xd+yd*yd+zd*zd)
         ec = ec + e
      end if
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(ii)
         iuse = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         do j = 1, nion
            cscale(iion(j)) = 1.0d0
         end do
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(kk)
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (iuse .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  fik = fi * pchg(kk)
                  rew = aewald * r
                  erfterm = erfc (rew)
                  scale = cscale(kn)
                  if (use_group)  scale = scale * fgrp
                  scaleterm = scale - 1.0d0
                  e = (fik/r) * (erfterm+scaleterm)
c
c     increment the overall charge-charge energy component
c
                  ec = ec + e
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate real space portion involving other unit cells
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(ii)
         iuse = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         do j = 1, nion
            cscale(iion(j)) = 1.0d0
         end do
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii, nion
            k = iion(kk)
            kn = jion(kk)
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (iuse .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  xr = xi - x(k)
                  yr = yi - y(k)
                  zr = zi - z(k)
                  call image (xr,yr,zr,j)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     fik = fi * pchg(kk)
                     rew = aewald * r
                     erfterm = erfc (rew)
                     scale = 1.0d0
                     if (use_group)  scale = scale * fgrp
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           scale = scale * cscale(kn)
                        end if
                     end if
                     scaleterm = scale - 1.0d0
                     e = (fik/r) * (erfterm+scaleterm)
c
c     increment the overall charge-charge energy component;
c     interaction of an atom with its own image counts half
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ec = ec + e
                  end if
               end do
            end if
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge0e  --  Ewald charge energy via lights  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge0e" calculates the charge-charge interaction energy
c     using a particle mesh Ewald summation and the method of lights
c     to locate neighboring atoms
c
c
      subroutine echarge0e
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'ewald.i'
      include 'group.i'
      include 'light.i'
      include 'math.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k
      integer ii,kk
      integer in,kn
      integer kgy,kgz
      integer start,stop
      integer map(maxlight)
      real*8 e,fs,fgrp
      real*8 f,fi,fik
      real*8 r,r2,rew
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xd,yd,zd
      real*8 erfc,erfterm
      real*8 scale,scaleterm
      real*8 cscale(maxatm)
      real*8 xsort(maxlight)
      real*8 ysort(maxlight)
      real*8 zsort(maxlight)
      logical proceed,iuse,repeat
      external erfc
c
c
c     zero out the Ewald charge interaction energy
c
      ec = 0.0d0
      if (nion .eq. 0)  return
c
c     set conversion factor, cutoff and scaling coefficients
c
      f = electric / dielec
      call switch ('EWALD')
c
c     compute the Ewald self-energy term over all the atoms
c
      fs = -f * aewald / sqrtpi
      do ii = 1, nion
         e = fs * pchg(ii)**2
         ec = ec + e
      end do
c
c     compute the reciprocal space part of the Ewald summation
c
      call epme
c
c     compute the cell dipole boundary correction term
c
      if (.not. tinfoil) then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do ii = 1, nion
            i = iion(ii)
            xd = xd + pchg(ii)*x(i)
            yd = yd + pchg(ii)*y(i)
            zd = zd + pchg(ii)*z(i)
         end do
         e = (2.0d0/3.0d0) * f * (pi/volbox) * (xd*xd+yd*yd+zd*zd)
         ec = ec + e
      end if
c
c     compute the real space portion of the Ewald summation;
c     transfer the interaction site coordinates to sorting arrays
c
      do i = 1, nion
         k = iion(i)
         xsort(i) = x(k)
         ysort(i) = y(k)
         zsort(i) = z(k)
      end do
c
c     use the method of lights to generate neighbors
c
      call lights (nion,map,xsort,ysort,zsort)
c
c     now, loop over all atoms computing the interactions
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         iuse = (use(i))
         do j = 1, nion
            cscale(iion(j)) = 1.0d0
         end do
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            k = iion(map(kk))
            kn = jion(map(kk))
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)
               if (use_image) then
                  if (abs(xr) .gt. xcell2)  xr = xr - sign(xcell,xr)
                  if (abs(yr) .gt. ycell2)  yr = yr - sign(ycell,yr)
                  if (abs(zr) .gt. zcell2)  zr = zr - sign(zcell,zr)
                  if (monoclinic) then
                     xr = xr + zr*beta_cos
                     zr = zr * beta_sin
                  else if (triclinic) then
                     xr = xr + yr*gamma_cos + zr*beta_cos
                     yr = yr*gamma_sin + zr*beta_term
                     zr = zr * gamma_term
                  end if
               end if
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rew = aewald * r
                  erfterm = erfc (rew)
                  scale = cscale(kn)
                  if (use_group)  scale = scale * fgrp
                  fik = fi * pchg(map(kk))
                  if (use_polymer) then
                     if (r2 .gt. polycut2)  fik = fi * pchg(map(kk))
                  end if
                  scaleterm = scale - 1.0d0
                  e = (fik/r) * (erfterm+scaleterm)
c
c     increment the overall charge-charge energy component
c
                  ec = ec + e
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = nlight
            goto 10
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine epme  --  PME reciprocal space charge energy  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "epme" computes the reciprocal space energy for a particle
c     mesh Ewald summation over partial charges
c
c     literature reference:
c
c     U. Essmann, L. Perera, M. L Berkowitz, T. Darden, H. Lee and
c     L. G. Pedersen, "A Smooth Particle Mesh Ewald Method", Journal
c     of Chemical Physics, 103, 8577-8593 (1995)
c
c
      subroutine epme
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'energi.i'
      include 'ewald.i'
      include 'math.i'
      include 'pme.i'
      include 'units.i'
      integer i,j,k,m,ii
      integer i0,j0,k0
      integer it1,it2,it3
      integer k1,k2,k3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      integer nff,npoint
      integer ndim1,ndim2,ndim3
      integer ntable,nwork,ngrid
      real*8 xi,yi,zi
      real*8 w,chg,denom
      real*8 e,term,expterm
      real*8 pterm,volterm
      real*8 f,hsq,struc2
      real*8 t1,t2,t3
      real*8 h1,h2,h3
      real*8 r1,r2,r3
      real*8 fr1(maxatm)
      real*8 fr2(maxatm)
      real*8 fr3(maxatm)
      real*8 theta1(maxorder,maxatm)
      real*8 theta2(maxorder,maxatm)
      real*8 theta3(maxorder,maxatm)
      real*8 work(2*maxfft)
c
c
c     initialization of some array dimensions
c
      ndim1 = 2*(nfft1/2) + 1
      ndim2 = 2*(nfft2/2) + 1
      ndim3 = 2*(nfft3/2) + 1
      ntable = maxtable
      nwork = maxfft
      ngrid = maxgrid
      f = sqrt(electric/dielec)
c
c     zero out the 3-dimensional charge grid array
c
      do k = 1, ndim3
         do j = 1, ndim2
            do i = 1, ndim1
               qgrid(1,i,j,k) = 0.0d0
               qgrid(2,i,j,k) = 0.0d0
            end do
         end do
      end do
c
c     get the scaled and shifted fractional coordinates
c
      do ii = 1, nion
         i = iion(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         w = xi*recip(1,1) + yi*recip(2,1) + zi*recip(3,1)
         fr1(i) = nfft1 * (w-dble(nint(w))+0.5d0)
         w = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
         fr2(i) = nfft2 * (w-dble(nint(w))+0.5d0)
         w = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
         fr3(i) = nfft3 * (w-dble(nint(w))+0.5d0)
c
c     get the B-spline coefficients and their first derivatives
c
         w = fr1(i) - dble(int(fr1(i)))
         call bspline (w,bsorder,theta1(1,i))
         w = fr2(i) - dble(int(fr2(i)))
         call bspline (w,bsorder,theta2(1,i))
         w = fr3(i) - dble(int(fr3(i)))
         call bspline (w,bsorder,theta3(1,i))
      end do
c
c     fill the 3-dimensional charge grid array
c
      do ii = 1, nion
         m = iion(ii)
         chg = f * pchg(ii)
         k0 = int(fr3(m)) - bsorder
         do it3 = 1, bsorder
            t3 = theta3(it3,m) * chg
            k0 = k0 + 1
            k = k0 + 1 + (nfft3-sign(nfft3,k0))/2
            j0 = int(fr2(m)) - bsorder
            do it2 = 1, bsorder
               t2 = theta2(it2,m)
               term = t2 * t3
               j0 = j0 + 1
               j = j0 + 1 + (nfft2-sign(nfft2,j0))/2
               i0 = int(fr1(m)) - bsorder
               do it1 = 1, bsorder
                  t1 = theta1(it1,m)
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-sign(nfft1,i0))/2
                  qgrid(1,i,j,k) = qgrid(1,i,j,k) + term*t1
               end do
            end do
         end do
      end do
c
c     perform the 3-D FFT forward transformation
c
      call fftfront (nfft1,nfft2,nfft3,ngrid,ngrid,ngrid,
     &                  ntable,nwork,table,work,qgrid)
c
c     sum the reciprocal space interaction energy components
c
      npoint = nfft1 * nfft2 * nfft3
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      do i = 1, npoint-1
         k3 = i/nff + 1
         j = i - (k3-1)*nff
         k2 = j/nfft1 + 1
         k1 = j - (k2-1)*nfft1 + 1
         m1 = k1 - 1
         m2 = k2 - 1
         m3 = k3 - 1
         if (k1 .gt. nf1)  m1 = m1 - nfft1
         if (k2 .gt. nf2)  m2 = m2 - nfft2
         if (k3 .gt. nf3)  m3 = m3 - nfft3
         r1 = dble(m1)
         r2 = dble(m2)
         r3 = dble(m3)
         h1 = recip(1,1)*r1
         h2 = recip(2,1)*r1 + recip(2,2)*r2
         h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
         hsq = h1*h1 + h2*h2 + h3*h3
         term = -pterm * hsq
         expterm = 0.0d0
         if (term .gt. -50.0d0) then
            denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
            expterm = exp(term) / denom
            if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0 
            end if
            struc2 = qgrid(1,k1,k2,k3)**2 + qgrid(2,k1,k2,k3)**2
            e = 0.5d0 * expterm * struc2
            ec = ec + e
         end if
         qgrid(1,k1,k2,k3) = expterm * qgrid(1,k1,k2,k3)
         qgrid(2,k1,k2,k3) = expterm * qgrid(2,k1,k2,k3)
      end do
      return
      end
