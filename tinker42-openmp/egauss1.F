 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine egauss1  --  Gaussian vdw energy & derivatives  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "egauss1" calculates the Gaussian expansion van der Waals
c     interaction energy and its first derivatives with respect
c     to Cartesian coordinates
c
c
      subroutine egauss1
      implicit none
      include 'sizes.i'
      include 'warp.i'
c
c
c     choose standard or potential energy smoothing version
c
      if (use_smooth) then
         call egauss1b
      else
         call egauss1a
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine egauss1a  --  double-loop Gaussian vdw derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "egauss1a" calculates the Gaussian expansion van der Waals
c     interaction energy and its first derivatives using a pairwise
c     double loop
c
c
      subroutine egauss1a
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'molcul.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k,ii,kk
      integer iv,kv,it,kt
      integer iv14(maxatm)
      real*8 e,de,rdn
      real*8 rik,rik2
      real*8 eps,rad2,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 dedx,dedy,dedz
      real*8 expcut,expterm
      real*8 a(4),b(4)
      real*8 xred(maxatm)
      real*8 yred(maxatm)
      real*8 zred(maxatm)
      real*8 vscale(maxatm)
      logical proceed,iuse
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      do i = 1, n
         dev(1,i) = 0.0d0
         dev(2,i) = 0.0d0
         dev(3,i) = 0.0d0
      end do
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
      expcut = -50.0d0
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
c     find van der Waals energy and derivatives via double loop
c
      do ii = 1, nvdw-1
         i = ivdw(ii)
         iv = ired(i)
         redi = kred(i)
         rediv = 1.0d0 - redi
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
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  eps = epsilon(kt,it)
                  rad2 = radmin(kt,it)**2
                  if (iv14(k) .eq. i) then
                     eps = epsilon4(kt,it)
                     rad2 = radmin4(kt,it)**2
                  end if
                  eps = eps * vscale(k)
                  do j = 1, ngauss
                     a(j) = igauss(1,j) * eps
                     b(j) = igauss(2,j) / rad2
                  end do
                  e = 0.0d0
                  de = 0.0d0
                  rik = sqrt(rik2)
                  do j = 1, ngauss
                     expterm = -b(j) * rik2
                     if (expterm .gt. expcut) then
                        expterm = a(j)*exp(expterm)
                        e = e + expterm
                        de = de - 2.0d0*b(j)*rik*expterm
                     end if
                  end do
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if
c
c     find the chain rule terms for derivative components
c
                  de = de / rik
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                  ev = ev + e
                  if (i .eq. iv) then
                     dev(1,i) = dev(1,i) + dedx
                     dev(2,i) = dev(2,i) + dedy
                     dev(3,i) = dev(3,i) + dedz
                  else
                     dev(1,i) = dev(1,i) + dedx*redi
                     dev(2,i) = dev(2,i) + dedy*redi
                     dev(3,i) = dev(3,i) + dedz*redi
                     dev(1,iv) = dev(1,iv) + dedx*rediv
                     dev(2,iv) = dev(2,iv) + dedy*rediv
                     dev(3,iv) = dev(3,iv) + dedz*rediv
                  end if
                  if (k .eq. kv) then
                     dev(1,k) = dev(1,k) - dedx
                     dev(2,k) = dev(2,k) - dedy
                     dev(3,k) = dev(3,k) - dedz
                  else
                     redk = kred(k)
                     redkv = 1.0d0 - redk
                     dev(1,k) = dev(1,k) - dedx*redk
                     dev(2,k) = dev(2,k) - dedy*redk
                     dev(3,k) = dev(3,k) - dedz*redk
                     dev(1,kv) = dev(1,kv) - dedx*redkv
                     dev(2,kv) = dev(2,kv) - dedy*redkv
                     dev(3,kv) = dev(3,kv) - dedz*redkv
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
               end if
            end if
         end do
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine egauss1b  --  Gaussian vdw derivs for smoothing  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "egauss1b" calculates the Gaussian expansion van der Waals
c     interaction energy and its first derivatives for use with
c     stophat potential energy smoothing
c
c
      subroutine egauss1b
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'math.i'
      include 'molcul.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'warp.i'
      integer i,j,k,ii,kk
      integer iv,kv,it,kt
      integer iv14(maxatm)
      real*8 e,de,rdn
      real*8 rik,rik2
      real*8 eps,rad2,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 dedx,dedy,dedz
      real*8 expcut,broot
      real*8 erf,term,term2
      real*8 expterm,expterm2
      real*8 width,wterm
      real*8 t1,t2,a(4),b(4)
      real*8 xred(maxatm)
      real*8 yred(maxatm)
      real*8 zred(maxatm)
      real*8 vscale(maxatm)
      logical proceed,iuse
      external erf
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      do i = 1, n
         dev(1,i) = 0.0d0
         dev(2,i) = 0.0d0
         dev(3,i) = 0.0d0
      end do
c
c     zero out the array marking 1-4 vdw interactions
c
      do i = 1, n
         iv14(i) = 0
      end do
c
c     set the extent of smoothing to be performed
c
      expcut = -50.0d0
      width = 0.0d0
      if (use_dem) then
         width = 4.0d0 * diffv * deform
      else if (use_gda) then
         wterm = (2.0d0/3.0d0) * diffv
      else if (use_tophat) then
         width = max(diffv*deform,0.0001d0)
      end if
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
c     find van der Waals energy and derivatives via double loop
c
      do ii = 1, nvdw-1
         i = ivdw(ii)
         iv = ired(i)
         redi = kred(i)
         rediv = 1.0d0 - redi
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
c
c     check for an interaction distance less than the cutoff
c
               eps = epsilon(kt,it)
               rad2 = radmin(kt,it)**2
               if (iv14(k) .eq. i) then
                  eps = epsilon4(kt,it)
                  rad2 = radmin4(kt,it)**2
               end if
               eps = eps * vscale(k)
               do j = 1, ngauss
                  a(j) = igauss(1,j) * eps
                  b(j) = igauss(2,j) / rad2
               end do
               e = 0.0d0
               de = 0.0d0
               rik = sqrt(rik2)
c
c     transform the potential function via smoothing
c
               if (use_tophat) then
                  rik = sqrt(rik2)
                  do j = 1, ngauss
                     broot = sqrt(b(j))
                     expterm = -b(j) * (rik+width)**2
                     if (expterm .gt. expcut) then
                        expterm = exp(expterm)
                     else
                        expterm = 0.0d0
                     end if
                     expterm2 = -b(j) * (width-rik)**2
                     if (expterm2 .gt. expcut) then
                        expterm2 = exp(expterm2)
                     else
                        expterm2 = 0.0d0
                     end if
                     term = broot * (expterm-expterm2)
                     term = term + sqrtpi*b(j)*rik
     &                         * (erf(broot*(rik+width))
     &                           +erf(broot*(width-rik)))
                     e = e + term*a(j)/(b(j)*b(j)*broot)
                     term = expterm * (2.0d0*rik*b(j)*width+1.0d0)
                     term2 = expterm2 * (2.0d0*rik*b(j)*width-1.0d0)
                     de = de + a(j)*(term+term2)/(b(j)*b(j))
                  end do
                  term = 3.0d0 / (8.0d0*rik*width**3)
                  e = e * term
                  de = -de * term / rik
               else
                  if (use_gda)  width = wterm * (m2(i)+m2(k))
                  do j = 1, ngauss
                     t1 = 1.0d0 + b(j)*width
                     t2 = sqrt(t1**3)
                     expterm = -b(j) * rik2 / t1
                     if (expterm .gt. expcut) then
                        expterm = (a(j)/t2)*exp(expterm)
                        e = e + expterm
                        de = de - (2.0d0*b(j)*rik/t1)*expterm
                     end if
                  end do
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  de = de * fgrp
               end if
c
c     find the chain rule terms for derivative components
c
               de = de / rik
               dedx = de * xr
               dedy = de * yr
               dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
               ev = ev + e
               if (i .eq. iv) then
                  dev(1,i) = dev(1,i) + dedx
                  dev(2,i) = dev(2,i) + dedy
                  dev(3,i) = dev(3,i) + dedz
               else
                  dev(1,i) = dev(1,i) + dedx*redi
                  dev(2,i) = dev(2,i) + dedy*redi
                  dev(3,i) = dev(3,i) + dedz*redi
                  dev(1,iv) = dev(1,iv) + dedx*rediv
                  dev(2,iv) = dev(2,iv) + dedy*rediv
                  dev(3,iv) = dev(3,iv) + dedz*rediv
               end if
               if (k .eq. kv) then
                  dev(1,k) = dev(1,k) - dedx
                  dev(2,k) = dev(2,k) - dedy
                  dev(3,k) = dev(3,k) - dedz
               else
                  redk = kred(k)
                  redkv = 1.0d0 - redk
                  dev(1,k) = dev(1,k) - dedx*redk
                  dev(2,k) = dev(2,k) - dedy*redk
                  dev(3,k) = dev(3,k) - dedz*redk
                  dev(1,kv) = dev(1,kv) - dedx*redkv
                  dev(2,kv) = dev(2,kv) - dedy*redkv
                  dev(3,kv) = dev(3,kv) - dedz*redkv
               end if
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
