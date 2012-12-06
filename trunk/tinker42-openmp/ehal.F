 
C     $Id: ehal.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ehal  --  Buffered 14-7 van der Waals energy  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ehal" calculates the buffered 14-7 van der Waals energy
c
c
      subroutine ehal
      implicit none
      include 'cutoff.i'
c
c
c     choose pairwise double loop or method of lights version
c
      if (use_lights) then
         call ehal0b
      else
         call ehal0a
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal0a  --  Buffered 14-7 vdw via double loop  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal0a" calculates the buffered 14-7 van der Waals energy
c     using a pairwise double loop
c
c
      subroutine ehal0a
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k,ii,kk
      integer iv,kv,it,kt
      integer iv14(maxatm)
      real*8 e,rv,rv7
      real*8 eps,rdn,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rho,tau,taper
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,rik7
      real*8 xred(maxatm)
      real*8 yred(maxatm)
      real*8 zred(maxatm)
      real*8 vscale(maxatm)
      logical proceed,iuse
c
c
c     zero out the van der Waals energy contribution
c
      ev = 0.0d0
c
c     zero out the array marking 1-4 vdw interactions
c
      do i = 1, n
         iv14(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
c
c     apply any reduction factor to the atomic coordinates
c
      do k = 1, nvdw
         i = ivdw(k)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     find the van der Waals energy via double loop search
c
      do ii = 1, nvdw-1
         i = ivdw(ii)
         iv = ired(i)
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         iuse = (use(i) .or. use(iv))
         do j = ii+1, nvdw
            vscale(ivdw(j)) = 1.0d0
         end do
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nvdw
            k = ivdw(kk)
            kv = ired(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = class(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               if (use_image)  call image (xr,yr,zr,0)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  rik = sqrt(rik2)
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(k)
                  rv7 = rv**7
                  rik7 = rik**7
                  rho = rik7 + ghal*rv7
                  tau = (dhal+1.0d0) / (rik + dhal*rv)
                  e = eps * rv7 * tau**7
     &                   * ((ghal+1.0d0)*rv7/rho-2.0d0)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy components
c
                  ev = ev + e
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
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         iuse = (use(i) .or. use(iv))
         do j = ii, nvdw
            vscale(ivdw(j)) = 1.0d0
         end do
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii, nvdw
            k = ivdw(kk)
            kv = ired(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = class(k)
               do j = 1, ncell
                  xr = xi - xred(k)
                  yr = yi - yred(k)
                  zr = zi - zred(k)
                  call image (xr,yr,zr,j)
                  rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
                  if (rik2 .le. off2) then
                     rik = sqrt(rik2)
                     rv = radmin(kt,it)
                     eps = epsilon(kt,it)
                     if (use_polymer) then
                        if (rik2 .le. polycut2) then
                           if (iv14(k) .eq. i) then
                              rv = radmin4(kt,it)
                              eps = epsilon4(kt,it)
                           end if
                           eps = eps * vscale(k)
                        end if
                     end if
                     rv7 = rv**7
                     rik7 = rik**7
                     rho = rik7 + ghal*rv7
                     tau = (dhal+1.0d0) / (rik + dhal*rv)
                     e = eps * rv7 * tau**7
     &                      * ((ghal+1.0d0)*rv7/rho-2.0d0)
c
c     use energy switching if near the cutoff distance
c
                     if (rik2 .gt. cut2) then
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
                        rik5 = rik2 * rik3
                        taper = c5*rik5 + c4*rik4 + c3*rik3
     &                             + c2*rik2 + c1*rik + c0
                        e = e * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy component;
c     interaction of an atom with its own image counts half
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ev = ev + e
                  end if
               end do
            end if
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ehal0b  --  method of lights Buffered 14-7 vdw  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ehal0a" calculates the buffered 14-7 van der Waals energy
c     using the method of lights to locate neighboring atoms
c
c
      subroutine ehal0b
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'iounit.i'
      include 'light.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k,ii,kk
      integer iv,kv,it,kt
      integer kgy,kgz
      integer start,stop
      integer iv14(maxatm)
      integer map(maxlight)
      real*8 e,rv,rv7
      real*8 eps,rdn,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rho,tau,taper
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,rik7
      real*8 xred(maxatm)
      real*8 yred(maxatm)
      real*8 zred(maxatm)
      real*8 vscale(maxatm)
      real*8 xsort(maxlight)
      real*8 ysort(maxlight)
      real*8 zsort(maxlight)
      logical proceed,iuse,repeat
c
c
c     zero out the van der Waals energy contribution
c
      ev = 0.0d0
c
c     zero out the array marking 1-4 vdw interactions
c
      do i = 1, n
         iv14(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
c
c     apply any reduction factor to the atomic coordinates
c
      do j = 1, nvdw
         i = ivdw(j)
         iv = ired(i)
         rdn = kred(i)
         xred(j) = rdn*(x(i)-x(iv)) + x(iv)
         yred(j) = rdn*(y(i)-y(iv)) + y(iv)
         zred(j) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     transfer the interaction site coordinates to sorting arrays
c
      do i = 1, nvdw
         xsort(i) = xred(i)
         ysort(i) = yred(i)
         zsort(i) = zred(i)
      end do
c
c     use the method of lights to generate neighbors
c
      call lights (nvdw,map,xsort,ysort,zsort)
c
c     now, loop over all atoms computing the interactions
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         it = class(i)
         xi = xsort(rgx(ii))
         yi = ysort(rgy(ii))
         zi = zsort(rgz(ii))
         iuse = (use(i) .or. use(iv))
         do j = 1, nvdw
            vscale(ivdw(j)) = 1.0d0
         end do
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
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
            k = ivdw(map(kk))
            kv = ired(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = class(k)
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
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  rik = sqrt(rik2)
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(k)
                  rv7 = rv**7
                  rik7 = rik**7
                  rho = rik7 + ghal*rv7
                  tau = dhal / (rik + dhal*rv)
                  e = eps * rv7 * tau**7
     &                   * ((ghal+1.0d0)*rv7/rho-2.0d0)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                             + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy components
c
                  ev = ev + e
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
