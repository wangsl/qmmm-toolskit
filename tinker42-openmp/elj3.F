 
C     $Id: elj3.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine elj3  --  Lennard-Jones vdw energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "elj3" calculates the Lennard-Jones 6-12 van der Waals energy
c     and also partitions the energy among the atoms
c
c
      subroutine elj3
      implicit none
      include 'sizes.i'
      include 'cutoff.i'
      include 'warp.i'
c
c
c     choose double loop, method of lights or smoothing version
c
      if (use_stophat) then
         call elj3d
      else if (use_smooth) then
         call elj3c
      else if (use_lights) then
         call elj3b
      else
         call elj3a
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine elj3a  --  double loop Lennard-Jones analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "elj3a" calculates the Lennard-Jones 6-12 van der Waals
c     energy and also partitions the energy among the atoms using
c     a pairwise double loop
c
c
      subroutine elj3a
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'molcul.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k,ii,kk
      integer iv,kv,it,kt
      integer iv14(maxatm)
      real*8 e,p6,p12,eps
      real*8 rv,rdn,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8 xred(maxatm)
      real*8 yred(maxatm)
      real*8 zred(maxatm)
      real*8 vscale(maxatm)
      logical proceed,iuse
      logical header,huge
c
c
c     zero out the van der Waals energy and partitioning terms
c
      nev = 0
      ev = 0.0d0
      do i = 1, n
         aev(i) = 0.0d0
      end do
      header = .true.
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
            if (proceed)  proceed = (vscale(k) .ne. 0.0d0)
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
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(k)
                  p6 = rv**6 / rik2**3
                  p12 = p6 * p6
                  e = eps * (p12 - 2.0d0*p6)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
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
                  nev = nev + 1
                  ev = ev + e
                  aev(i) = aev(i) + 0.5d0*e
                  aev(k) = aev(k) + 0.5d0*e
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 10.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Individual van der Waals',
     &                             ' Interactions :',
     &                          //,' Type',13x,'Atom Names',
     &                             18x,'Minimum',4x,'Actual',
     &                             6x,'Energy',/)
                     end if
                     write (iout,20)  i,name(i),k,name(k),
     &                                rv,sqrt(rik2),e
   20                format (' VDW-LJ',5x,i5,'-',a3,1x,i5,
     &                          '-',a3,12x,2f10.4,f12.4)
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
                     p6 = rv**6 / rik2**3
                     p12 = p6 * p6
                     e = eps * (p12 - 2.0d0*p6)
c
c     use energy switching if near the cutoff distance
c
                     if (rik2 .gt. cut2) then
                        rik = sqrt(rik2)
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
                        rik5 = rik2 * rik3
                        taper = c5*rik5 + c4*rik4 + c3*rik3
     &                                + c2*rik2 + c1*rik + c0
                        e = e * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy components
c
                     if (e .ne. 0.0d0)  nev = nev + 1
                     if (i .eq. k) then
                        ev = ev + 0.5d0*e
                        aev(i) = aev(i) + 0.5d0*e
                     else
                        ev = ev + e
                        aev(i) = aev(i) + 0.5d0*e
                        aev(k) = aev(k) + 0.5d0*e
                     end if
c
c     increment the total intermolecular energy
c
                     einter = einter + e
c
c     print a message if the energy of this interaction is large
c
                     huge = (e .gt. 10.0d0)
                     if ((debug.and.e.ne.0.0d0)
     &                     .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,30)
   30                      format (/,' Individual van der Waals',
     &                                ' Interactions :',
     &                             //,' Type',13x,'Atom Names',
     &                                18x,'Minimum',4x,'Actual',
     &                                6x,'Energy',/)
                        end if
                        write (iout,40)  i,name(i),k,name(k),
     &                                   rv,sqrt(rik2),e
   40                   format (' VDW-LJ',5x,i5,'-',a3,1x,i5,'-',a3,
     &                             3x,'(XTAL)',3x,2f10.4,f12.4)
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine elj3b  --  Lennard-Jones analysis via lights  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "elj3b" calculates the Lennard-Jones 6-12 van der Waals energy
c     and also partitions the energy among the atoms using the method
c     of lights to locate neighboring atoms
c
c
      subroutine elj3b
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'light.i'
      include 'molcul.i'
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
      real*8 e,p6,p12,eps
      real*8 rv,rdn,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8 xred(maxatm)
      real*8 yred(maxatm)
      real*8 zred(maxatm)
      real*8 vscale(maxatm)
      real*8 xsort(maxlight)
      real*8 ysort(maxlight)
      real*8 zsort(maxlight)
      logical proceed,iuse,repeat
      logical header,huge
c
c
c     zero out the van der Waals energy and partitioning terms
c
      nev = 0
      ev = 0.0d0
      do i = 1, n
         aev(i) = 0.0d0
      end do
      header = .true.
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
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(k)
                  p6 = rv**6 / rik2**3
                  p12 = p6 * p6
                  e = eps * (p12 - 2.0d0*p6)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
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
                  if (e .ne. 0.0d0)  nev = nev + 1
                  ev = ev + e
                  aev(i) = aev(i) + 0.5d0*e
                  aev(k) = aev(k) + 0.5d0*e
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 10.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual van der Waals',
     &                             ' Interactions :',
     &                          //,' Type',13x,'Atom Names',
     &                             18x,'Minimum',4x,'Actual',
     &                             6x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),
     &                                rv,sqrt(rik2),e
   30                format (' VDW-LJ',5x,i5,'-',a3,1x,i5,'-',
     &                          a3,12x,2f10.4,f12.4)
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
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine elj3c  --  Lennard-Jones analysis for smoothing  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "elj3c" calculates the Lennard-Jones 6-12 van der Waals energy
c     and also partitions the energy among the atoms via a Gaussian
c     approximation for potential energy smoothing
c
c
      subroutine elj3c
      implicit none
      include 'math.i'
      include 'vdwpot.i'
c
c
c     set coefficients for a two-Gaussian fit to Lennard-Jones
c
      ngauss = 2
      igauss(1,1) = 14487.1d0
      igauss(2,1) = 9.05148d0 * twosix**2
      igauss(1,2) = -5.55338d0
      igauss(2,2) = 1.22536d0 * twosix**2
c
c     compute Gaussian approximation to Lennard-Jones potential
c
      call egauss3
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine elj3d  --  Lennard-Jones analysis for stophat  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "elj3d" calculates the Lennard-Jones 6-12 van der Waals energy
c     and also partitions the energy among the atoms for use with
c     stophat potential energy smoothing
c
c
      subroutine elj3d
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'molcul.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'warp.i'
      integer i,j,k,ii,kk
      integer iv,kv,it,kt
      integer iv14(maxatm)
      real*8 e,rik2,rdn,p6
      real*8 eps,rv,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rik,rik3,rik4
      real*8 rik5,rik6
      real*8 width,width2
      real*8 width3,width4
      real*8 width5,width6
      real*8 xred(maxatm)
      real*8 yred(maxatm)
      real*8 zred(maxatm)
      real*8 vscale(maxatm)
      logical proceed,iuse
      logical header,huge
c
c
c     zero out the van der Waals energy and partitioning terms
c
      nev = 0
      ev = 0.0d0
      do i = 1, n
         aev(i) = 0.0d0
      end do
      header = .true.
c
c     zero out the array marking 1-4 vdw interactions
c
      do i = 1, n
         iv14(i) = 0
      end do
c
c     set the extent of smoothing to be performed
c
      width = deform * diffv
      width2 = width * width
      width3 = width2 * width
      width4 = width2 * width2
      width5 = width2 * width3
      width6 = width3 * width3
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
               rik2 = xr*xr + yr*yr + zr*zr
               eps = epsilon(kt,it)
               rv = radmin(kt,it)
               if (iv14(k) .eq. i) then
                  eps = epsilon4(kt,it)
                  rv = radmin4(kt,it)
               end if
               eps = eps * vscale(k)
               p6 = rv**6
               rik = sqrt(rik2)
               rik3 = rik2 * rik
               rik4 = rik2 * rik2
               rik5 = rik2 * rik3
               rik6 = rik3 * rik3
c
c     transform the potential function via smoothing
c
               e = rik6 * (30.0d0*rik6 + 360.0d0*rik5*width
     &                + 1800.0d0*rik4*width2 + 4800.0d0*rik3*width3
     &                + 7200.0d0*rik2*width4 + 5760.0d0*rik*width5
     &                + 1920.0d0*width6)
               e = -e + p6 * (15.0d0*rik6 + 90.0d0*rik5*width
     &                + 288.0d0*rik4*width2 + 552.0d0*rik3*width3
     &                + 648.0d0*rik2*width4 + 432.0d0*rik*width5
     &                + 128.0d0*width6)
               e = e*eps*p6 / (15.0d0*(rik*(rik+2.0d0*width))**9)
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy components
c
               if (e .ne. 0.0d0)  nev = nev + 1
               ev = ev + e
               aev(i) = aev(i) + 0.5d0*e
               aev(k) = aev(k) + 0.5d0*e
c
c     increment the total intermolecular energy
c
               if (molcule(i) .ne. molcule(k)) then
                  einter = einter + e
               end if
c
c     print a message if the energy of this interaction is large
c
               huge = (e .gt. 10.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Individual van der Waals',
     &                          ' Interactions :',
     &                       //,' Type',13x,'Atom Names',
     &                          18x,'Minimum',4x,'Actual',
     &                          6x,'Energy',/)
                  end if
                  write (iout,20)  i,name(i),k,name(k),
     &                             rv,sqrt(rik2),e
   20             format (' VDW-LJ',5x,i5,'-',a3,1x,i5,
     &                       '-',a3,12x,2f10.4,f12.4)
               end if
            end if
         end do
      end do
      return
      end
