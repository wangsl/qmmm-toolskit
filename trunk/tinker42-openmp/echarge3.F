 
C     $Id: echarge3.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge3  --  charge-charge energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge3" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms
c
c
      subroutine echarge3
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
         call echarge3e
      else if (use_ewald) then
         call echarge3d
      else if (use_smooth) then
         call echarge3c
      else if (use_lights) then
         call echarge3b
      else
         call echarge3a
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine echarge3a  --  charge analysis via double loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "echarge3a" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms using a pairwise
c     double loop
c
c
      subroutine echarge3a
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'molcul.i'
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
      logical header,huge
c
c
c     zero out the charge interaction energy and partitioning
c
      nec = 0
      ec = 0.0d0
      do i = 1, n
         aec(i) = 0.0d0
      end do
      if (nion .eq. 0)  return
      header = .true.
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('CHARGE')
c
c     compute and partition the charge interaction energy
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
            if (proceed)  proceed = (cscale(kn) .ne. 0.0d0)
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
                  nec = nec + 1
                  ec = ec + e
                  aec(i) = aec(i) + 0.5d0*e
                  aec(k) = aec(k) + 0.5d0*e
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 100.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Individual Charge-Charge',
     &                             ' Interactions :',
     &                          //,' Type',13x,'Atom Names',
     &                             16x,'Charges',5x,'Distance',
     &                             5x,'Energy',/)
                     end if
                     write (iout,20)  i,name(i),k,name(k),
     &                                pchg(ii),pchg(kk),r,e
   20                format (' Charge',5x,i5,'-',a3,1x,i5,'-',
     &                          a3,8x,2f7.2,f10.4,f12.4)
                  end if
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
                     if (e .ne. 0.0d0)  nec = nec + 1
                     ec = ec + e
                     aec(i) = aec(i) + 0.5d0*e
                     aec(k) = aec(k) + 0.5d0*e
c
c     increment the total intermolecular energy
c
                     einter = einter + e
c
c     print a message if the energy of this interaction is large
c
                     huge = (abs(e) .gt. 100.0d0)
                     if ((debug.and.e.ne.0.0d0)
     &                     .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,30)
   30                      format (/,' Individual Charge-Charge',
     &                                ' Interactions :',
     &                             //,' Type',13x,'Atom Names',
     &                                16x,'Charges',5x,'Distance',
     &                                5x,'Energy',/)
                        end if
                        write (iout,40)  i,name(i),k,name(k),
     &                                   pchg(ii),pchg(kk),r,e
   40                   format (' Charge',5x,i5,'-',a3,1x,i5,'-',a3,
     &                             1x,'(XTAL)',1x,2f7.2,f10.4,f12.4)
                     end if
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
c     ##  subroutine echarge3b  --  method of lights charge analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "echarge3b" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms using the method
c     of lights to locate neighboring atoms
c
c
      subroutine echarge3b
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'light.i'
      include 'molcul.i'
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
      logical header,huge
c
c
c     zero out the charge interaction energy and partitioning
c
      nec = 0
      ec = 0.0d0
      do i = 1, n
         aec(i) = 0.0d0
      end do
      if (nion .eq. 0)  return
      header = .true.
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
         iuse = (use(i) .or. use(ic))
         xic = xsort(rgx(ii))
         yic = ysort(rgy(ii))
         zic = zsort(rgz(ii))
         xi = x(i) - x(ic)
         yi = y(i) - y(ic)
         zi = z(i) - z(ic)
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
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 40
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 40
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 40
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 40
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
            if (proceed)  proceed = (cscale(kc) .ne. 0.0d0)
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
                  if (e .ne. 0.0d0)  nec = nec + 1
                  ec = ec + e
                  aec(i) = aec(i) + 0.5d0*e
                  aec(k) = aec(k) + 0.5d0*e
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 100.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Charge-Charge',
     &                             ' Interactions :',
     &                          //,' Type',13x,'Atom Names',
     &                             16x,'Charges',5x,'Distance',
     &                             5x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),pchg(ii),
     &                                pchg(map(kk)),r,e
   30                format (' Charge',5x,i5,'-',a3,1x,i5,'-',
     &                          a3,8x,2f7.2,f10.4,f12.4)
                  end if
               end if
            end if
   40       continue
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
c     ##  subroutine echarge3c  --  charge analysis for smoothing  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "echarge3c" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms for use with
c     potential smoothing methods
c
c
      subroutine echarge3c
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'molcul.i'
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
      logical header,huge
      external erf
c
c
c     zero out the charge interaction energy and partitioning
c
      nec = 0
      ec = 0.0d0
      do i = 1, n
         aec(i) = 0.0d0
      end do
      if (nion .eq. 0)  return
      header = .true.
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
c     compute and partition the charge interaction energy
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
            if (proceed)  proceed = (cscale(kn) .ne. 0.0d0)
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
               nec = nec + 1
               ec = ec + e
               aec(i) = aec(i) + 0.5d0*e
               aec(k) = aec(k) + 0.5d0*e
c
c     increment the total intermolecular energy
c
               if (molcule(i) .ne. molcule(k)) then
                  einter = einter + e
               end if
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(e) .gt. 100.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Individual Charge-Charge',
     &                          ' Interactions :',
     &                       //,' Type',11x,'Atom Names',
     &                          16x,'Charges',5x,'Distance',
     &                          5x,'Energy',/)
                  end if
                  write (iout,20)  i,name(i),k,name(k),pchg(ii),
     &                             pchg(kk),r,e
   20             format (' Charge',5x,i5,'-',a3,1x,i5,'-',
     &                       a3,8x,2f7.2,f10.4,f12.4)
               end if
            end if
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine echarge3d  --  Ewald summation charge analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "echarge3d" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms using a particle
c     mesh Ewald summation
c
c
      subroutine echarge3d
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
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
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'math.i'
      include 'molcul.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k
      integer ii,kk
      integer in,kn
      real*8 e,eintra
      real*8 f,fi,fik,fs
      real*8 r,r2,rew,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xd,yd,zd
      real*8 erfc,erfterm
      real*8 scale,scaleterm
      real*8 cscale(maxatm)
      logical proceed,iuse
      logical header,huge
      external erfc
c
c
c     zero out the Ewald summation energy and partitioning
c
      nec = 0
      ec = 0.0d0
      do i = 1, n
         aec(i) = 0.0d0
      end do
      if (nion .eq. 0)  return
      header = .true.
c
c     zero out the intramolecular portion of the Ewald energy
c
      eintra = 0.0d0
c
c     set conversion factor and cutoff coefficients
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
      call epme3
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
c
c     increment the overall charge-charge energy component
c
                  nec = nec + 1
                  ec = ec + e
                  aec(i) = aec(i) + 0.5d0*e
                  aec(k) = aec(k) + 0.5d0*e
c
c     print a message if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 100.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Individual Real Space Ewald',
     &                             ' Charge-Charge Interactions :',
     &                          //,' Type',13x,'Atom Names',
     &                             16x,'Charges',5x,'Distance',
     &                             5x,'Energy',/)
                     end if
                     write (iout,20)  i,name(i),k,name(k),
     &                                pchg(ii),pchg(kk),r,e
   20                format (' Charge',5x,i5,'-',a3,1x,i5,'-',
     &                          a3,8x,2f7.2,f10.4,f12.4)
                  end if
               end if
            end if
         end do
      end do
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
c     increment the overall charge-charge energy component
c
                     if (i .eq. k)  e = 0.5d0 * e
                     if (e .ne. 0.0d0)  nec = nec + 1
                     ec = ec + e
                     aec(i) = aec(i) + 0.5d0*e
                     aec(k) = aec(k) + 0.5d0*e
c
c     print a message if the energy of this interaction is large
c
                     huge = (abs(e) .gt. 100.0d0)
                     if ((debug.and.e.ne.0.0d0)
     &                     .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,30)
   30                      format (/,' Individual Real Space Ewald',
     &                                ' Charge-Charge Interactions :',
     &                             //,' Type',13x,'Atom Names',
     &                                16x,'Charges',5x,'Distance',
     &                                5x,'Energy',/)
                        end if
                        write (iout,40)  i,name(i),k,name(k),
     &                                   pchg(ii),pchg(kk),r,e
   40                   format (' Charge',5x,i5,'-',a3,1x,i5,'-',a3,
     &                             1x,'(XTAL)',1x,2f7.2,f10.4,f12.4)
                     end if
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
c     ##  subroutine echarge3e  --  Ewald charge analysis via lights  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "echarge3e" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms using a particle
c     mesh Ewald summation and the method of lights to locate
c     neighboring atoms
c
c
      subroutine echarge3e
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
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
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'light.i'
      include 'math.i'
      include 'molcul.i'
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
      logical header,huge
      external erfc
c
c
c     zero out the Ewald summation energy and partitioning
c
      nec = 0
      ec = 0.0d0
      do i = 1, n
         aec(i) = 0.0d0
      end do
      if (nion .eq. 0)  return
      header = .true.
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
      call epme3
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
         iuse = use(i)
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
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 40
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 40
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 40
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 40
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
                  if (e .ne. 0.0d0)  nec = nec + 1
                  ec = ec + e
                  aec(i) = aec(i) + 0.5d0*e
                  aec(k) = aec(k) + 0.5d0*e
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 100.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Real Space Ewald',
     &                             ' Charge-Charge Interactions :',
     &                          //,' Type',13x,'Atom Names',
     &                             16x,'Charges',5x,'Distance',
     &                             5x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),pchg(ii),
     &                                pchg(map(kk)),r,e
   30                format (' Charge',5x,i5,'-',a3,1x,i5,'-',
     &                          a3,8x,2f7.2,f10.4,f12.4)
                  end if
               end if
            end if
   40       continue
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
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine epme3  --  PME reciprocal energy & analysis  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "epme3" computes the reciprocal space energy for a particle
c     mesh Ewald summation over partial charges and prints information
c     about the energy over the charge grid points
c
c     literature reference:
c
c     U. Essmann, L. Perera, M. L Berkowitz, T. Darden, H. Lee and
c     L. G. Pedersen, "A Smooth Particle Mesh Ewald Method", Journal
c     of Chemical Physics, 103, 8577-8593 (1995)
c
c
      subroutine epme3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'energi.i'
      include 'ewald.i'
      include 'inform.i'
      include 'iounit.i'
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
      logical header,huge
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
      header = .true.
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
            nec = nec + 1
            ec = ec + e
         end if
         qgrid(1,k1,k2,k3) = expterm * qgrid(1,k1,k2,k3)
         qgrid(2,k1,k2,k3) = expterm * qgrid(2,k1,k2,k3)
c
c     print a message if the energy of this interaction is large
c
         huge = (abs(e) .gt. 10.0d0)
         if (debug .or. (verbose.and.huge)) then
            if (header) then
               header = .false.
               write (iout,10)
   10          format (/,' Reciprocal Space Ewald',
     &                    ' Charge-Charge Interactions :',
     &                 //,' Type',15x,'K-Vector',18x,'Fraction',
     &                    15x,'Energy',/)
            end if
            write (iout,20)  m1,m2,m3,hsq,e
   20       format (' Charge',7x,3i5,16x,f8.4,10x,f12.4)
         end if
      end do
      return
      end
