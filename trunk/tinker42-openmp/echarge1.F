 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge1  --  charge-charge energy & derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge1" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine echarge1
      implicit none
      include 'sizes.i'
      include 'cutoff.i'
      include 'warp.i'
c
c
c     choose pairwise double loop, method of lights, potential
c     smoothing or Ewald summation version
c
      if (use_ewald) then
         call echarge1d
      else if (use_smooth) then
         call echarge1c
      else if (use_lights) then
         call echarge1b
      else
         call echarge1a
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine echarge1a  --  charge derivs via double loop  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "echarge1a" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c     using a pairwise double loop
c
c
      subroutine echarge1a
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'molcul.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,k
      integer ii,kk
      integer in,kn
      integer ic,kc
      real*8 e,de,dc
      real*8 f,fi,fik
      real*8 r,r2,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xc,yc,zc
      real*8 xic,yic,zic
      real*8 dedx,dedy,dedz
      real*8 dedxc,dedyc,dedzc
      real*8 shift,taper,trans
      real*8 dtaper,dtrans
      real*8 rc,rc2,rc3,rc4
      real*8 rc5,rc6,rc7
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 cscale(maxatm)
      logical proceed,iuse

      Real*8 EInterChg
c
c
c     zero out the charge interaction energy and derivatives
c
      ec = 0.0d0
C$OMP Parallel Do
C$OMP& Default(Shared), Private(I)
      do i = 1, n
         dec(1,i) = 0.0d0
         dec(2,i) = 0.0d0
         dec(3,i) = 0.0d0
      end do
C$OMP End Parallel Do
      if (nion .eq. 0)  return
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('CHARGE')
c
c     compute the charge interaction energy and first derivatives
c

C     add for OpenMP

      EInterChg = 0.0D0
      Vxx = 0.0D0
      Vyx = 0.0D0
      Vzx = 0.0D0
      Vyy = 0.0D0
      Vzy = 0.0D0
      Vzz = 0.0D0

C$OMP Parallel Do Schedule(Static, 1)
C$OMP& Default(Shared)
C$OMP& Private(II, I, IN, IC, XIC, YIC, ZIC, XI, YI, ZI, FI, IUse)
C$OMP& Private(J, CScale, KK, K, KN, KC, Proceed, FGrp)
C$OMP& Private(XC, YC, ZC, RC2, XR, YR, ZR, R2, R, FIK, E, DE, DC)
C$OMP& Private(Shift, RC, RC3, RC4, RC5, RC6, RC7)
C$OMP& Private(Taper, DTaper, Trans, DTrans)
C$OMP& Private(DEDX, DEDY, DEDZ, DEDXC, DEDYC, DEDZC)
C$OMP& Reduction(+:EC, DEC)
C$OMP& Reduction(+:VXX, VYX, VZX, VYY, VZY, VZZ, EInterChg)
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
                  de = -fik / r2
                  dc = 0.0d0
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
                     dtaper = 5.0d0*c5*rc4 + 4.0d0*c4*rc3
     &                           + 3.0d0*c3*rc2 + 2.0d0*c2*rc + c1
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     dtrans = fik * (7.0d0*f7*rc6 + 6.0d0*f6*rc5
     &                               + 5.0d0*f5*rc4 + 4.0d0*f4*rc3
     &                             + 3.0d0*f3*rc2 + 2.0d0*f2*rc + f1)
                     dc = (e * dtaper + dtrans) / rc
                     de = de * taper
                     e = e * taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     dc = dc * fgrp
                  end if
c
c     form the chain rule terms for derivative expressions
c
                  de = de / r
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  dedxc = dc * xc
                  dedyc = dc * yc
                  dedzc = dc * zc
c
c     increment the overall energy and derivative expressions
c
                  ec = ec + e
                  dec(1,i) = dec(1,i) + dedx
                  dec(2,i) = dec(2,i) + dedy
                  dec(3,i) = dec(3,i) + dedz
                  dec(1,ic) = dec(1,ic) + dedxc
                  dec(2,ic) = dec(2,ic) + dedyc
                  dec(3,ic) = dec(3,ic) + dedzc
                  dec(1,k) = dec(1,k) - dedx
                  dec(2,k) = dec(2,k) - dedy
                  dec(3,k) = dec(3,k) - dedz
                  dec(1,kc) = dec(1,kc) - dedxc
                  dec(2,kc) = dec(2,kc) - dedyc
                  dec(3,kc) = dec(3,kc) - dedzc
c
c     increment the internal virial tensor components
c
                  vxx = Vxx + xr*dedx + xc*dedxc
                  vyx = Vyx + yr*dedx + yc*dedxc
                  vzx = Vzx + zr*dedx + zc*dedxc
                  vyy = Vyy + yr*dedy + yc*dedyc
                  vzy = Vzy + zr*dedy + zc*dedyc
                  vzz = Vzz + zr*dedz + zc*dedzc
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
C                     einter = einter + e
                     EInterChg = EInterChg + E
                  end if
               end if
            end if
         end do
      end do
C$OMP End Parallel Do
      
      EInter = EInter + EInterChg
      vir(1,1) = vir(1,1) + vxx
      vir(2,1) = vir(2,1) + vyx
      vir(3,1) = vir(3,1) + vzx
      vir(1,2) = vir(1,2) + vyx
      vir(2,2) = vir(2,2) + vyy
      vir(3,2) = vir(3,2) + vzy
      vir(1,3) = vir(1,3) + vzx
      vir(2,3) = vir(2,3) + vzy
      vir(3,3) = vir(3,3) + vzz

C     Call WriteEGradients(EC, EInter, DEC, Vir)

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
                     de = -fik / r2
                     dc = 0.0d0
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
                        dtaper = 5.0d0*c5*rc4 + 4.0d0*c4*rc3
     &                              + 3.0d0*c3*rc2 + 2.0d0*c2*rc + c1
                        trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                                  + f3*rc3 + f2*rc2 + f1*rc + f0)
                        dtrans = fik * (7.0d0*f7*rc6 + 6.0d0*f6*rc5
     &                                  + 5.0d0*f5*rc4 + 4.0d0*f4*rc3
     &                                + 3.0d0*f3*rc2 + 2.0d0*f2*rc + f1)
                        dc = (e * dtaper + dtrans) / rc
                        de = de * taper
                        e = e * taper + trans
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        e = e * fgrp
                        de = de * fgrp
                        dc = dc * fgrp
                     end if
c
c     form the chain rule terms for derivative expressions
c
                     de = de / r
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
                     dedxc = dc * xc
                     dedyc = dc * yc
                     dedzc = dc * zc
c
c     increment the energy and gradient values
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ec = ec + e
                     dec(1,i) = dec(1,i) + dedx
                     dec(2,i) = dec(2,i) + dedy
                     dec(3,i) = dec(3,i) + dedz
                     dec(1,ic) = dec(1,ic) + dedxc
                     dec(2,ic) = dec(2,ic) + dedyc
                     dec(3,ic) = dec(3,ic) + dedzc
                     if (i .ne. k) then
                        dec(1,k) = dec(1,k) - dedx
                        dec(2,k) = dec(2,k) - dedy
                        dec(3,k) = dec(3,k) - dedz
                        dec(1,kc) = dec(1,kc) - dedxc
                        dec(2,kc) = dec(2,kc) - dedyc
                        dec(3,kc) = dec(3,kc) - dedzc
                     end if
c
c     increment the internal virial tensor components
c
                     vxx = xr*dedx + xc*dedxc
                     vyx = yr*dedx + yc*dedxc
                     vzx = zr*dedx + zc*dedxc
                     vyy = yr*dedy + yc*dedyc
                     vzy = zr*dedy + zc*dedyc
                     vzz = zr*dedz + zc*dedzc
                     vir(1,1) = vir(1,1) + vxx
                     vir(2,1) = vir(2,1) + vyx
                     vir(3,1) = vir(3,1) + vzx
                     vir(1,2) = vir(1,2) + vyx
                     vir(2,2) = vir(2,2) + vyy
                     vir(3,2) = vir(3,2) + vzy
                     vir(1,3) = vir(1,3) + vzx
                     vir(2,3) = vir(2,3) + vzy
                     vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
                     einter = einter + e
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
c     ##  subroutine echarge1b  --  method of lights charge derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge1b" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c     using the method of lights to locate neighboring atoms
c
c
      subroutine echarge1b
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'iounit.i'
      include 'light.i'
      include 'molcul.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,k,ii,kk
      integer in,ic,kn,kc
      integer kgy,kgz
      integer start,stop
      integer map(maxlight)
      real*8 e,de,dc
      real*8 f,fi,fik
      real*8 r,r2,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xc,yc,zc
      real*8 xic,yic,zic
      real*8 dedx,dedy,dedz
      real*8 dedxc,dedyc,dedzc
      real*8 shift,taper,trans
      real*8 dtaper,dtrans
      real*8 rc,rc2,rc3,rc4
      real*8 rc5,rc6,rc7
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 cscale(maxatm)
      real*8 xsort(maxlight)
      real*8 ysort(maxlight)
      real*8 zsort(maxlight)
      logical proceed,iuse,repeat
c
c
c     zero out the charge interaction energy and derivatives
c
      ec = 0.0d0
      do i = 1, n
         dec(1,i) = 0.0d0
         dec(2,i) = 0.0d0
         dec(3,i) = 0.0d0
      end do
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
                  de = -fik / r2
                  dc = 0.0d0
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
                     dtaper = 5.0d0*c5*rc4 + 4.0d0*c4*rc3
     &                           + 3.0d0*c3*rc2 + 2.0d0*c2*rc + c1
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     dtrans = fik * (7.0d0*f7*rc6 + 6.0d0*f6*rc5
     &                               + 5.0d0*f5*rc4 + 4.0d0*f4*rc3
     &                             + 3.0d0*f3*rc2 + 2.0d0*f2*rc + f1)
                     dc = (e * dtaper + dtrans) / rc
                     de = de * taper
                     e = e * taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     dc = dc * fgrp
                  end if
c
c     form the chain rule terms for derivative expressions
c
                  de = de / r
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  dedxc = dc * xc
                  dedyc = dc * yc
                  dedzc = dc * zc
c
c     increment the overall energy and derivative expressions
c
                  ec = ec + e
                  dec(1,i) = dec(1,i) + dedx
                  dec(2,i) = dec(2,i) + dedy
                  dec(3,i) = dec(3,i) + dedz
                  dec(1,ic) = dec(1,ic) + dedxc
                  dec(2,ic) = dec(2,ic) + dedyc
                  dec(3,ic) = dec(3,ic) + dedzc
                  dec(1,k) = dec(1,k) - dedx
                  dec(2,k) = dec(2,k) - dedy
                  dec(3,k) = dec(3,k) - dedz
                  dec(1,kc) = dec(1,kc) - dedxc
                  dec(2,kc) = dec(2,kc) - dedyc
                  dec(3,kc) = dec(3,kc) - dedzc
c
c     increment the internal virial tensor components
c
                  vxx = xr*dedx + xc*dedxc
                  vyx = yr*dedx + yc*dedxc
                  vzx = zr*dedx + zc*dedxc
                  vyy = yr*dedy + yc*dedyc
                  vzy = zr*dedy + zc*dedyc
                  vzz = zr*dedz + zc*dedzc
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
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
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine echarge1c  --  charge derivatives for smoothing  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "echarge1c" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c     for use with potential smoothing methods
c
c
      subroutine echarge1c
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'math.i'
      include 'molcul.i'
      include 'units.i'
      include 'usage.i'
      include 'warp.i'
      integer i,j,k
      integer ii,kk
      integer in,kn
      real*8 e,de,r,r2,fgrp
      real*8 f,fi,fik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 dedx,dedy,dedz
      real*8 erf,erfterm
      real*8 expcut,expterm
      real*8 wterm,width
      real*8 width2,width3
      real*8 cscale(maxatm)
      logical proceed,iuse
      external erf
c
c
c     zero out the charge interaction energy and derivatives
c
      ec = 0.0d0
C$OMP Parallel Do
C$OMP& Default(Shared), Private(I)
      do i = 1, n
         dec(1,i) = 0.0d0
         dec(2,i) = 0.0d0
         dec(3,i) = 0.0d0
      end do
C$OMP End Parallel Do
      if (nion .eq. 0)  return
c
c     set the energy units conversion factor
c
      f = electric / dielec
      expcut = -50.0d0
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
c     compute the charge interaction energy and first derivatives
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
               de = -fik / r2
c
c     transform the potential function via smoothing
c
               if (use_dem) then
                  if (width .gt. 0.0d0) then
                     erfterm = erf(width*r)
                     expterm = -r2 * width2
                     if (expterm .gt. expcut) then
                        expterm = 2.0d0*fik*width*exp(expterm)
     &                               / (sqrtpi*r)
                     else
                        expterm = 0.0d0
                     end if
                     e = e * erfterm
                     de = de*erfterm + expterm
                  end if
               else if (use_gda) then
                  width = m2(i) + m2(k)
                  if (width .gt. 0.0d0) then
                     width = wterm / sqrt(width)
                     erfterm = erf(width*r)
                     expterm = -r2 * width * width
                     if (expterm .gt. expcut) then
                        expterm = 2.0d0*fik*width*exp(expterm)
     &                               / (sqrtpi*r)
                     else
                        expterm = 0.0d0
                     end if
                     e = e * erfterm
                     de = de*erfterm + expterm
                  end if
               else if (use_tophat) then
                  if (width .gt. r) then
                     e = fik * (3.0d0*width2-r2) / (2.0d0*width3)
                     de = -fik * r / width3
                  end if
               else if (use_stophat) then
                  wterm = r + width
                  e = fik / wterm
                  de = -fik / (wterm*wterm)
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  de = de * fgrp
               end if
c
c     form the chain rule terms for derivative expressions
c
               de = de / r
               dedx = de * xr
               dedy = de * yr
               dedz = de * zr
c
c     increment the overall energy and derivative expressions
c
               ec = ec + e
               dec(1,i) = dec(1,i) + dedx
               dec(2,i) = dec(2,i) + dedy
               dec(3,i) = dec(3,i) + dedz
               dec(1,k) = dec(1,k) - dedx
               dec(2,k) = dec(2,k) - dedy
               dec(3,k) = dec(3,k) - dedz
c
c     increment the total intermolecular energy
c
               if (molcule(i) .ne. molcule(k)) then
                  einter = einter + e
               end if
            end if
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine echarge1d  --  Ewald summation charge derivs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "echarge1d" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c     using a particle mesh Ewald summation
c
c
      subroutine echarge1d
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'ewald.i'
      include 'group.i'
      include 'inter.i'
      include 'math.i'
      include 'molcul.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,k
      integer ii,kk
      integer in,kn
      real*8 e,de,eintra
      real*8 f,fi,fik,fs
      real*8 r,r2,rew
      real*8 fgrp,term
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xd,yd,zd
      real*8 erfc,erfterm
      real*8 scale,scaleterm
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 cscale(maxatm)
      logical proceed,iuse
      external erfc

      Real*8 ECTmp
c
c
c     zero out the Ewald summation energy and derivatives
c
      ec = 0.0d0
C$OMP Parallel Do
C$OMP& Default(Shared), Private(I)
      do i = 1, n
         dec(1,i) = 0.0d0
         dec(2,i) = 0.0d0
         dec(3,i) = 0.0d0
      end do
C$OMP End Parallel Do
      if (nion .eq. 0)  return
c
c     zero out the intramolecular portion of the Ewald energy
c
      eintra = 0.0d0
c
c     set conversion factor, cutoff and scaling coefficients
c
      f = electric / dielec
      call switch ('EWALD')
c
c     compute the Ewald self-energy term over all the atoms
c
      fs = -f * aewald / sqrtpi
C$OMP Parallel Do
C$OMP& Default(Shared), Private(II, E)
C$OMP& Reduction(+:EC)
      do ii = 1, nion
         e = fs * pchg(ii)**2
         ec = ec + e
      end do
C$OMP End Parallel Do
c
c     compute reciprocal space Ewald energy and first derivatives
c
      call epme1
c
c     compute the cell dipole boundary correction term
c
      if (.not. tinfoil) then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
C$OMP Parallel Do
C$OMP& Default(Shared), Private(II, I)
C$OMP& Reduction(+:XD, YD, ZD)
         do ii = 1, nion
            i = iion(ii)
            xd = xd + pchg(ii)*x(i)
            yd = yd + pchg(ii)*y(i)
            zd = zd + pchg(ii)*z(i)
         end do
C$OMP End Parallel Do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         e = term * (xd*xd+yd*yd+zd*zd)
         ec = ec + e
C$OMP Parallel Do
C$OMP& Default(Shared), Private(II, I, DE, DEDX, DEDY, DEDZ)
C$OMP& Reduction(+:DEC)
         do ii = 1, nion
            i = iion(ii)
            de = 2.0d0 * term * pchg(ii)
            dedx = de * xd
            dedy = de * yd
            dedz = de * zd
            dec(1,i) = dec(1,i) + dedx
            dec(2,i) = dec(2,i) + dedy
            dec(3,i) = dec(3,i) + dedz
         end do
C$OMP End Parallel Do
      end if
c
c     compute the real space Ewald energy and first derivatives
c
      Vxx = 0.0D0
      Vyx = 0.0D0
      Vzx = 0.0D0
      Vyy = 0.0D0
      Vzy = 0.0D0
      Vzz = 0.0D0
      ECTmp = EC
      EC = 0.0D0

C$OMP Parallel Do Schedule(Static, 1)
C$OMP& Default(Shared)
C$OMP& Private(II, I, IN, IUse, XI, YI, ZI, FI)
C$OMP& Private(J, CScale, KK, K, KN, FGrp, Proceed)
C$OMP& Private(XR, YR, ZR)
C$OMP& Private(R2, R, FIK, E)
C$OMP& Private(REW, ErfTerm, Scale, ScaleTerm, DE)
C$OMP& Private(DEDX, DEDY, DEDZ)
C$OMP& Reduction(+:EIntra, EC, DEC)
C$OMP& Reduction(+:Vxx, Vyx, Vzx, Vyy, Vzy, Vzz)
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
c
c     increment the total intramolecular energy
c
               if (molcule(i) .eq. molcule(k)) then
                  r2 = xr*xr + yr*yr + zr*zr
                  r = sqrt(r2)
                  fik = fi * pchg(kk) * cscale(kn)
                  e = fik / r
                  eintra = eintra + e
               end if
c
c     find minimum image and check the real space cutoff
c
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
                  de = -fik * ((erfterm+scaleterm)/r2
     &                    + (2.0d0*aewald/sqrtpi)*exp(-rew**2)/r)
c
c     form the chain rule terms for derivative expressions
c
                  de = de / r
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the overall energy and derivative expressions
c
                  ec = ec + e
                  dec(1,i) = dec(1,i) + dedx
                  dec(2,i) = dec(2,i) + dedy
                  dec(3,i) = dec(3,i) + dedz
                  dec(1,k) = dec(1,k) - dedx
                  dec(2,k) = dec(2,k) - dedy
                  dec(3,k) = dec(3,k) - dedz
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx + Vxx
                  vyx = yr * dedx + Vyx
                  vzx = zr * dedx + Vzx
                  vyy = yr * dedy + Vyy
                  vzy = zr * dedy + Vzy
                  vzz = zr * dedz + Vzz
                  
               end if
            end if
         end do
      end do
C$OMP End Parallel Do

      EC = ECTmp + EC
      vir(1,1) = vir(1,1) + vxx
      vir(2,1) = vir(2,1) + vyx
      vir(3,1) = vir(3,1) + vzx
      vir(1,2) = vir(1,2) + vyx
      vir(2,2) = vir(2,2) + vyy
      vir(3,2) = vir(3,2) + vzy
      vir(1,3) = vir(1,3) + vzx
      vir(2,3) = vir(2,3) + vzy
      vir(3,3) = vir(3,3) + vzz
c     
c     intermolecular energy is total minus intramolecular part
c
      einter = einter + ec - eintra
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
c
c     find appropriate image and check the real space cutoff
c
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
                     de = -fik * ((erfterm+scaleterm)/r2
     &                       + (2.0d0*aewald/sqrtpi)*exp(-rew**2)/r)
c
c     form the chain rule terms for derivative expressions
c
                     de = de / r
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
c
c     increment the energy and gradient values
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ec = ec + e
                     dec(1,i) = dec(1,i) + dedx
                     dec(2,i) = dec(2,i) + dedy
                     dec(3,i) = dec(3,i) + dedz
                     if (i .ne. k) then
                        dec(1,k) = dec(1,k) - dedx
                        dec(2,k) = dec(2,k) - dedy
                        dec(3,k) = dec(3,k) - dedz
                     end if
c
c     increment the internal virial tensor components
c
                     vxx = xr * dedx
                     vyx = yr * dedx
                     vzx = zr * dedx
                     vyy = yr * dedy
                     vzy = zr * dedy
                     vzz = zr * dedz
                     vir(1,1) = vir(1,1) + vxx
                     vir(2,1) = vir(2,1) + vyx
                     vir(3,1) = vir(3,1) + vzx
                     vir(1,2) = vir(1,2) + vyx
                     vir(2,2) = vir(2,2) + vyy
                     vir(3,2) = vir(3,2) + vzy
                     vir(1,3) = vir(1,3) + vzx
                     vir(2,3) = vir(2,3) + vzy
                     vir(3,3) = vir(3,3) + vzz
                  end if
               end do
            end if
         end do
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine epme1  --  PME reciprocal space energy & derivs  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epme1" computes the reciprocal space energy and first
c     derivatives for a particle mesh Ewald summation
c
c     literature reference:
c
c     U. Essmann, L. Perera, M. L Berkowitz, T. Darden, H. Lee and
c     L. G. Pedersen, "A Smooth Particle Mesh Ewald Method", Journal
c     of Chemical Physics, 103, 8577-8593 (1995)
c
c
      subroutine epme1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'deriv.i'
      include 'energi.i'
      include 'ewald.i'
      include 'math.i'
      include 'pme.i'
      include 'units.i'
      include 'virial.i'
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
      real*8 vterm,pterm
      real*8 volterm
      real*8 f,hsq,struc2
      real*8 de1,de2,de3
      real*8 dn1,dn2,dn3
      real*8 t1,t2,t3
      real*8 dt1,dt2,dt3
      real*8 h1,h2,h3
      real*8 r1,r2,r3
      real*8 fr1(maxatm)
      real*8 fr2(maxatm)
      real*8 fr3(maxatm)
      real*8 theta1(maxorder,maxatm)
      real*8 theta2(maxorder,maxatm)
      real*8 theta3(maxorder,maxatm)
      real*8 dtheta1(maxorder,maxatm)
      real*8 dtheta2(maxorder,maxatm)
      real*8 dtheta3(maxorder,maxatm)
      real*8 work(2*maxfft)

      Real*8 ECTmp, Vxx, Vyx, Vzx, Vyy, Vzy, Vzz

c
c
c     initialization of some array dimensions and parameters
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
C$OMP Parallel Do
C$OMP& Default(Shared), Private(K, J, I)
      do k = 1, ndim3
         do j = 1, ndim2
            do i = 1, ndim1
               qgrid(1,i,j,k) = 0.0d0
               qgrid(2,i,j,k) = 0.0d0
            end do
         end do
      end do
C$OMP End Parallel Do
c
c     get the scaled and shifted fractional coordinates
c
C$OMP Parallel Do Schedule(Static, 1)
C$OMP& Default(Shared)
C$OMP& Private(II, I, XI, YI, ZI, W)
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
         call bspline1 (w,bsorder,theta1(1,i),dtheta1(1,i))
         w = fr2(i) - dble(int(fr2(i)))
         call bspline1 (w,bsorder,theta2(1,i),dtheta2(1,i))
         w = fr3(i) - dble(int(fr3(i)))
         call bspline1 (w,bsorder,theta3(1,i),dtheta3(1,i))
      end do
C$OMP End Parallel Do
c
c     fill the 3-dimensional charge grid array
c
C$OMP Parallel Do Schedule(Static, 1)
C$OMP& Default(Shared)
C$OMP& Private(II, M, Chg, K0)
C$OMP& Private(It3, T3, K, J0)
C$OMP& Private(It2, T2, Term)
C$OMP& Private(J, I0, It1, T1, I)
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
C$OMP Atomic
                  qgrid(1,i,j,k) = qgrid(1,i,j,k) + term*t1
               end do
            end do
         end do
      end do
C$OMP End Parallel Do

c     perform the 3-D FFT forward transformation
c
      call fftfront (nfft1,nfft2,nfft3,ngrid,ngrid,ngrid,
     &                  ntable,nwork,table,work,qgrid)
c
c     sum the reciprocal space energy and virial components
c
      npoint = nfft1 * nfft2 * nfft3
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      
      ECTmp = 0.0D0
      Vxx = 0.0D0
      Vyx = 0.0D0
      Vzx = 0.0D0
      Vyy = 0.0D0
      Vzy = 0.0D0
      Vzz = 0.0D0
C$OMP Parallel Do Schedule(Static, 1)
C$OMP& Default(Shared)
C$OMP& Private(I, K3, J, K2, K1, M1, M2, M3)
C$OMP& Private(R1, R2, R3, H1, H2, H3, HSq, Term, ExpTerm)
C$OMP& Private(Denom, Struc2, E)
C$OMP& Reduction(+:ECTmp, Vxx, Vyx, Vzx, Vyy, Vzy, Vzz)
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
#if 0
            ec = ec + e
#endif
            ECTmp = ECTmp + E
            vterm = (2.0d0/hsq) * (1.0d0-term) * e
            Vxx = Vxx + h1*h1*vterm - e
            Vyx = Vyx + h1*h2*vterm
            Vzx = Vzx + h1*h3*vterm
            Vyy = Vyy + h2*h2*vterm - e
            Vzy = Vzy + h2*h3*vterm
            Vzz = Vzz + h3*h3*vterm - e
#if 0
            vir(1,1) = vir(1,1) + h1*h1*vterm - e
            vir(2,1) = vir(2,1) + h1*h2*vterm
            vir(3,1) = vir(3,1) + h1*h3*vterm
            vir(1,2) = vir(1,2) + h2*h1*vterm
            vir(2,2) = vir(2,2) + h2*h2*vterm - e
            vir(3,2) = vir(3,2) + h2*h3*vterm
            vir(1,3) = vir(1,3) + h3*h1*vterm
            vir(2,3) = vir(2,3) + h3*h2*vterm
            vir(3,3) = vir(3,3) + h3*h3*vterm - e
#endif
         end if
         qgrid(1,k1,k2,k3) = expterm * qgrid(1,k1,k2,k3)
         qgrid(2,k1,k2,k3) = expterm * qgrid(2,k1,k2,k3)
      end do
C$OMP End Parallel Do

      EC = EC + ECTmp
      vir(1,1) = vir(1,1) + Vxx
      vir(2,1) = vir(2,1) + Vyx
      vir(3,1) = vir(3,1) + Vzx
      vir(1,2) = vir(1,2) + Vyx
      vir(2,2) = vir(2,2) + Vyy
      vir(3,2) = vir(3,2) + Vzy
      vir(1,3) = vir(1,3) + Vzx
      vir(2,3) = vir(2,3) + Vzy
      vir(3,3) = vir(3,3) + Vzz
c
c     perform the 3-D FFT backward transformation
c
      call fftback (nfft1,nfft2,nfft3,ngrid,ngrid,ngrid,
     &                 ntable,nwork,table,work,qgrid)
c
c     get first derivatives of the reciprocal space energy
c
      dn1 = dble(nfft1)
      dn2 = dble(nfft2)
      dn3 = dble(nfft3)

C$OMP Parallel Do Schedule(Static, 1)
C$OMP& Default(Shared)
C$OMP& Private(II, M, Chg, DE1, DE2, DE3, K0)
C$OMP& Private(IT3, T3, DT3, K, J0)
C$OMP& Private(IT2, T2, DT2, J, I0)
C$OMP& Private(IT1, T1, DT1, I, Term)
      do ii = 1, nion
         m = iion(ii)
         chg = f * pchg(ii)
         de1 = 0.0d0
         de2 = 0.0d0
         de3 = 0.0d0
         k0 = int(fr3(m)) - bsorder
         do it3 = 1, bsorder
            t3 = theta3(it3,m)
            dt3 = dn3 * dtheta3(it3,m)
            k0 = k0 + 1
            k = k0 + 1 + (nfft3-sign(nfft3,k0))/2
            j0 = int(fr2(m)) - bsorder
            do it2 = 1, bsorder
               t2 = theta2(it2,m)
               dt2 = dn2 * dtheta2(it2,m)
               j0 = j0 + 1
               j = j0 + 1 + (nfft2-sign(nfft2,j0))/2
               i0 = int(fr1(m)) - bsorder
               do it1 = 1, bsorder
                  t1 = theta1(it1,m)
                  dt1 = dn1 * dtheta1(it1,m)
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-sign(nfft1,i0))/2
                  term = qgrid(1,i,j,k)
                  de1 = de1 + term*dt1*t2*t3
                  de2 = de2 + term*dt2*t1*t3
                  de3 = de3 + term*dt3*t1*t2
               end do
            end do
         end do
         dec(1,m) = dec(1,m) + chg*(recip(1,1)*de1)
         dec(2,m) = dec(2,m) + chg*(recip(2,1)*de1+recip(2,2)*de2)
         dec(3,m) = dec(3,m) + chg*(recip(3,1)*de1+recip(3,2)*de2
     &                                   +recip(3,3)*de3)
      end do
C$OMP End Parallel Do
      return
      end
