 
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
c     ##  subroutine egauss2  --  atom-by-atom Gaussian vdw Hessian  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "egauss2" calculates the Gaussian expansion van der Waals
c     second derivatives for a single atom at a time
c
c
      subroutine egauss2 (i,xred,yred,zred)
      implicit none
      include 'sizes.i'
      include 'warp.i'
      integer i
      real*8 xred(maxatm)
      real*8 yred(maxatm)
      real*8 zred(maxatm)
c
c
c     choose standard or potential energy smoothing version
c
      if (use_smooth) then
         call egauss2b (i,xred,yred,zred)
      else
         call egauss2a (i,xred,yred,zred)
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine egauss2a  --  double loop Gaussian vdw Hessian  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "egauss2a" calculates the Gaussian expansion van der Waals
c     second derivatives using a pairwise double loop
c
c
      subroutine egauss2a (iatom,xred,yred,zred)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'group.i'
      include 'hessn.i'
      include 'shunt.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k,iatom
      integer ii,iv,it
      integer kk,kv,kt
      integer nuse,use(5)
      integer iv14(maxatm)
      real*8 de,d2e,rik,rik2
      real*8 eps,rad2,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 redi2,rediv2
      real*8 rediiv
      real*8 redik,redivk
      real*8 redikv,redivkv
      real*8 d2edx,d2edy,d2edz
      real*8 expcut,expterm
      real*8 a(4),b(4)
      real*8 term(3,3)
      real*8 xred(maxatm)
      real*8 yred(maxatm)
      real*8 zred(maxatm)
      real*8 vscale(maxatm)
      logical proceed
c
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
c     check to see if the atom of interest is a vdw site
c
      nuse = 0
      do k = 1, nvdw
         if (ivdw(k) .eq. iatom) then
            nuse = nuse + 1
            use(nuse) = iatom
            goto 10
         end if
      end do
      return
   10 continue
c
c     determine the atoms involved via reduction factors
c
      nuse = 1
      use(nuse) = iatom
      do k = 1, n12(iatom)
         i = i12(k,iatom)
         if (ired(i) .eq. iatom) then
            nuse = nuse + 1
            use(nuse) = i
         end if
      end do
c
c     find van der Waals Hessian elements for involved atoms
c
      do ii = 1, nuse
         i = use(ii)
         iv = ired(i)
         redi = kred(i)
         if (i .ne. iv) then
            rediv = 1.0d0 - redi
            redi2 = redi * redi
            rediv2 = rediv * rediv
            rediiv = redi * rediv
         end if
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
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
c
c     decide whether to compute the current interaction
c
         do kk = 1, nvdw
            k = ivdw(kk)
            kv = ired(k)
            proceed = (k .ne. i)
            if (proceed .and. use_group)
     &         call groups (proceed,fgrp,i,k,0,0,0,0)
c
c     compute the Hessian elements for this interaction
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
                  de = 0.0d0
                  d2e = 0.0d0
                  rik = sqrt(rik2)
                  do j = 1, ngauss
                     expterm = -b(j) * rik2
                     if (expterm .gt. expcut) then
                        expterm = a(j)*b(j)*exp(expterm)
                        de = de - 2.0d0*rik*expterm
                        d2e = d2e + (4.0d0*b(j)*rik2-2.0d0)*expterm
                     end if
                  end do
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     de = de * fgrp
                     d2e = d2e * fgrp
                  end if
c
c     get chain rule terms for van der Waals Hessian elements
c
                  de = de / rik
                  d2e = (d2e-de) / rik2
                  d2edx = d2e * xr
                  d2edy = d2e * yr
                  d2edz = d2e * zr
                  term(1,1) = d2edx*xr + de
                  term(1,2) = d2edx*yr
                  term(1,3) = d2edx*zr
                  term(2,1) = term(1,2)
                  term(2,2) = d2edy*yr + de
                  term(2,3) = d2edy*zr
                  term(3,1) = term(1,3)
                  term(3,2) = term(2,3)
                  term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
                  if (i .eq. iatom) then
                     if (i.eq.iv .and. k.eq.kv) then
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)
                           hessy(j,i) = hessy(j,i) + term(2,j)
                           hessz(j,i) = hessz(j,i) + term(3,j)
                           hessx(j,k) = hessx(j,k) - term(1,j)
                           hessy(j,k) = hessy(j,k) - term(2,j)
                           hessz(j,k) = hessz(j,k) - term(3,j)
                        end do
                     else if (k .eq. kv) then
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)*redi2
                           hessy(j,i) = hessy(j,i) + term(2,j)*redi2
                           hessz(j,i) = hessz(j,i) + term(3,j)*redi2
                           hessx(j,k) = hessx(j,k) - term(1,j)*redi
                           hessy(j,k) = hessy(j,k) - term(2,j)*redi
                           hessz(j,k) = hessz(j,k) - term(3,j)*redi
                           hessx(j,iv) = hessx(j,iv) + term(1,j)*rediiv
                           hessy(j,iv) = hessy(j,iv) + term(2,j)*rediiv
                           hessz(j,iv) = hessz(j,iv) + term(3,j)*rediiv
                        end do
                     else if (i .eq. iv) then
                        redk = kred(k)
                        redkv = 1.0d0 - redk
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)
                           hessy(j,i) = hessy(j,i) + term(2,j)
                           hessz(j,i) = hessz(j,i) + term(3,j)
                           hessx(j,k) = hessx(j,k) - term(1,j)*redk
                           hessy(j,k) = hessy(j,k) - term(2,j)*redk
                           hessz(j,k) = hessz(j,k) - term(3,j)*redk
                           hessx(j,kv) = hessx(j,kv) - term(1,j)*redkv
                           hessy(j,kv) = hessy(j,kv) - term(2,j)*redkv
                           hessz(j,kv) = hessz(j,kv) - term(3,j)*redkv
                        end do
                     else
                        redk = kred(k)
                        redkv = 1.0d0 - redk
                        redik = redi * redk
                        redikv = redi * redkv
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)*redi2
                           hessy(j,i) = hessy(j,i) + term(2,j)*redi2
                           hessz(j,i) = hessz(j,i) + term(3,j)*redi2
                           hessx(j,k) = hessx(j,k) - term(1,j)*redik
                           hessy(j,k) = hessy(j,k) - term(2,j)*redik
                           hessz(j,k) = hessz(j,k) - term(3,j)*redik
                           hessx(j,iv) = hessx(j,iv) + term(1,j)*rediiv
                           hessy(j,iv) = hessy(j,iv) + term(2,j)*rediiv
                           hessz(j,iv) = hessz(j,iv) + term(3,j)*rediiv
                           hessx(j,kv) = hessx(j,kv) - term(1,j)*redikv
                           hessy(j,kv) = hessy(j,kv) - term(2,j)*redikv
                           hessz(j,kv) = hessz(j,kv) - term(3,j)*redikv
                        end do
                     end if
                  else if (iv .eq. iatom) then
                     if (k .eq. kv) then
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)*rediiv
                           hessy(j,i) = hessy(j,i) + term(2,j)*rediiv
                           hessz(j,i) = hessz(j,i) + term(3,j)*rediiv
                           hessx(j,k) = hessx(j,k) - term(1,j)*rediv
                           hessy(j,k) = hessy(j,k) - term(2,j)*rediv
                           hessz(j,k) = hessz(j,k) - term(3,j)*rediv
                           hessx(j,iv) = hessx(j,iv) + term(1,j)*rediv2
                           hessy(j,iv) = hessy(j,iv) + term(2,j)*rediv2
                           hessz(j,iv) = hessz(j,iv) + term(3,j)*rediv2
                        end do
                     else
                        redk = kred(k)
                        redkv = 1.0d0 - redk
                        redivk = rediv * redk
                        redivkv = rediv * redkv
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)*rediiv
                           hessy(j,i) = hessy(j,i) + term(2,j)*rediiv
                           hessz(j,i) = hessz(j,i) + term(3,j)*rediiv
                           hessx(j,k) = hessx(j,k) - term(1,j)*redivk
                           hessy(j,k) = hessy(j,k) - term(2,j)*redivk
                           hessz(j,k) = hessz(j,k) - term(3,j)*redivk
                           hessx(j,iv) = hessx(j,iv) + term(1,j)*rediv2
                           hessy(j,iv) = hessy(j,iv) + term(2,j)*rediv2
                           hessz(j,iv) = hessz(j,iv) + term(3,j)*rediv2
                           hessx(j,kv) = hessx(j,kv) - term(1,j)*redivkv
                           hessy(j,kv) = hessy(j,kv) - term(2,j)*redivkv
                           hessz(j,kv) = hessz(j,kv) - term(3,j)*redivkv
                        end do
                     end if
                  end if
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
c     ##  subroutine egauss2b  --  Gaussian Hessian for smoothing  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "egauss2b" calculates the Gaussian expansion van der Waals
c     second derivatives for stophat potential energy smoothing
c
c
      subroutine egauss2b (iatom,xred,yred,zred)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'group.i'
      include 'hessn.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'warp.i'
      integer i,j,k,iatom
      integer ii,iv,it
      integer kk,kv,kt
      integer nuse,use(5)
      integer iv14(maxatm)
      real*8 de,d2e,rik,rik2
      real*8 eps,rad2,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 redi2,rediv2
      real*8 rediiv
      real*8 redik,redivk
      real*8 redikv,redivkv
      real*8 d2edx,d2edy,d2edz
      real*8 expcut,b2
      real*8 term1,term2
      real*8 width,wterm
      real*8 expterm,expterm2
      real*8 t1,t2,a(4),b(4)
      real*8 term(3,3)
      real*8 xred(maxatm)
      real*8 yred(maxatm)
      real*8 zred(maxatm)
      real*8 vscale(maxatm)
      logical proceed
c
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
c     check to see if the atom of interest is a vdw site
c
      nuse = 0
      do k = 1, nvdw
         if (ivdw(k) .eq. iatom) then
            nuse = nuse + 1
            use(nuse) = iatom
            goto 10
         end if
      end do
      return
   10 continue
c
c     determine the atoms involved via reduction factors
c
      nuse = 1
      use(nuse) = iatom
      do k = 1, n12(iatom)
         i = i12(k,iatom)
         if (ired(i) .eq. iatom) then
            nuse = nuse + 1
            use(nuse) = i
         end if
      end do
c
c     find van der Waals Hessian elements for involved atoms
c
      do ii = 1, nuse
         i = use(ii)
         iv = ired(i)
         redi = kred(i)
         if (i .ne. iv) then
            rediv = 1.0d0 - redi
            redi2 = redi * redi
            rediv2 = rediv * rediv
            rediiv = redi * rediv
         end if
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
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
c
c     decide whether to compute the current interaction
c
         do kk = 1, nvdw
            k = ivdw(kk)
            kv = ired(k)
            proceed = (k .ne. i)
            if (proceed .and. use_group)
     &         call groups (proceed,fgrp,i,k,0,0,0,0)
c
c     compute the Hessian elements for this interaction
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
               de = 0.0d0
               d2e = 0.0d0
               rik = sqrt(rik2)
c
c     transform the potential function via smoothing
c
               if (use_tophat) then
                  rik = sqrt(rik2)
                  do j = 1, ngauss
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
                     b2 = b(j)*b(j)
                     term1 = expterm * (2.0d0*rik*b(j)*width+1.0d0)
                     term2 = expterm2 * (2.0d0*rik*b(j)*width-1.0d0)
                     de = de + a(j)*(term1+term2)/b2
                     term1 = 2.0d0*b(j)*width*rik * (b(j)*rik2+1.0d0)
     &                          * (expterm+expterm2)
                     term2 = (2.0d0*rik2*(b(j)*width)**2 + 1.0d0
     &                          + b(j)*rik2) * (expterm-expterm2)
                     d2e = d2e + a(j)*(term1+term2)/b2
                  end do
                  term1 = 3.0d0 / (8.0d0*rik*width**3)
                  de = -de * term1
                  d2e = 2.0d0 * d2e * term1 / rik
               else
                  if (use_gda)  width = wterm * (m2(i)+m2(k))
                  do j = 1, ngauss
                     t1 = 1.0d0 + b(j)*width
                     t2 = sqrt(t1**3)
                     expterm = -b(j) * rik2 / t1
                     if (expterm .gt. expcut) then
                        expterm = (a(j)*b(j)/(t2*t1))*exp(expterm)
                        de = de - 2.0d0*rik*expterm
                        d2e = d2e + (4.0d0*b(j)*rik2/t1-2.0d0)*expterm
                     end if
                  end do
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  de = de * fgrp
                  d2e = d2e * fgrp
               end if
c
c     get chain rule terms for van der Waals Hessian elements
c
               de = de / rik
               d2e = (d2e-de) / rik2
               d2edx = d2e * xr
               d2edy = d2e * yr
               d2edz = d2e * zr
               term(1,1) = d2edx*xr + de
               term(1,2) = d2edx*yr
               term(1,3) = d2edx*zr
               term(2,1) = term(1,2)
               term(2,2) = d2edy*yr + de
               term(2,3) = d2edy*zr
               term(3,1) = term(1,3)
               term(3,2) = term(2,3)
               term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
               if (i .eq. iatom) then
                  if (i.eq.iv .and. k.eq.kv) then
                     do j = 1, 3
                        hessx(j,i) = hessx(j,i) + term(1,j)
                        hessy(j,i) = hessy(j,i) + term(2,j)
                        hessz(j,i) = hessz(j,i) + term(3,j)
                        hessx(j,k) = hessx(j,k) - term(1,j)
                        hessy(j,k) = hessy(j,k) - term(2,j)
                        hessz(j,k) = hessz(j,k) - term(3,j)
                     end do
                  else if (k .eq. kv) then
                     do j = 1, 3
                        hessx(j,i) = hessx(j,i) + term(1,j)*redi2
                        hessy(j,i) = hessy(j,i) + term(2,j)*redi2
                        hessz(j,i) = hessz(j,i) + term(3,j)*redi2
                        hessx(j,k) = hessx(j,k) - term(1,j)*redi
                        hessy(j,k) = hessy(j,k) - term(2,j)*redi
                        hessz(j,k) = hessz(j,k) - term(3,j)*redi
                        hessx(j,iv) = hessx(j,iv) + term(1,j)*rediiv
                        hessy(j,iv) = hessy(j,iv) + term(2,j)*rediiv
                        hessz(j,iv) = hessz(j,iv) + term(3,j)*rediiv
                     end do
                  else if (i .eq. iv) then
                     redk = kred(k)
                     redkv = 1.0d0 - redk
                     do j = 1, 3
                        hessx(j,i) = hessx(j,i) + term(1,j)
                        hessy(j,i) = hessy(j,i) + term(2,j)
                        hessz(j,i) = hessz(j,i) + term(3,j)
                        hessx(j,k) = hessx(j,k) - term(1,j)*redk
                        hessy(j,k) = hessy(j,k) - term(2,j)*redk
                        hessz(j,k) = hessz(j,k) - term(3,j)*redk
                        hessx(j,kv) = hessx(j,kv) - term(1,j)*redkv
                        hessy(j,kv) = hessy(j,kv) - term(2,j)*redkv
                        hessz(j,kv) = hessz(j,kv) - term(3,j)*redkv
                     end do
                  else
                     redk = kred(k)
                     redkv = 1.0d0 - redk
                     redik = redi * redk
                     redikv = redi * redkv
                     do j = 1, 3
                        hessx(j,i) = hessx(j,i) + term(1,j)*redi2
                        hessy(j,i) = hessy(j,i) + term(2,j)*redi2
                        hessz(j,i) = hessz(j,i) + term(3,j)*redi2
                        hessx(j,k) = hessx(j,k) - term(1,j)*redik
                        hessy(j,k) = hessy(j,k) - term(2,j)*redik
                        hessz(j,k) = hessz(j,k) - term(3,j)*redik
                        hessx(j,iv) = hessx(j,iv) + term(1,j)*rediiv
                        hessy(j,iv) = hessy(j,iv) + term(2,j)*rediiv
                        hessz(j,iv) = hessz(j,iv) + term(3,j)*rediiv
                        hessx(j,kv) = hessx(j,kv) - term(1,j)*redikv
                        hessy(j,kv) = hessy(j,kv) - term(2,j)*redikv
                        hessz(j,kv) = hessz(j,kv) - term(3,j)*redikv
                     end do
                  end if
               else if (iv .eq. iatom) then
                  if (k .eq. kv) then
                     do j = 1, 3
                        hessx(j,i) = hessx(j,i) + term(1,j)*rediiv
                        hessy(j,i) = hessy(j,i) + term(2,j)*rediiv
                        hessz(j,i) = hessz(j,i) + term(3,j)*rediiv
                        hessx(j,k) = hessx(j,k) - term(1,j)*rediv
                        hessy(j,k) = hessy(j,k) - term(2,j)*rediv
                        hessz(j,k) = hessz(j,k) - term(3,j)*rediv
                        hessx(j,iv) = hessx(j,iv) + term(1,j)*rediv2
                        hessy(j,iv) = hessy(j,iv) + term(2,j)*rediv2
                        hessz(j,iv) = hessz(j,iv) + term(3,j)*rediv2
                     end do
                  else
                     redk = kred(k)
                     redkv = 1.0d0 - redk
                     redivk = rediv * redk
                     redivkv = rediv * redkv
                     do j = 1, 3
                        hessx(j,i) = hessx(j,i) + term(1,j)*rediiv
                        hessy(j,i) = hessy(j,i) + term(2,j)*rediiv
                        hessz(j,i) = hessz(j,i) + term(3,j)*rediiv
                        hessx(j,k) = hessx(j,k) - term(1,j)*redivk
                        hessy(j,k) = hessy(j,k) - term(2,j)*redivk
                        hessz(j,k) = hessz(j,k) - term(3,j)*redivk
                        hessx(j,iv) = hessx(j,iv) + term(1,j)*rediv2
                        hessy(j,iv) = hessy(j,iv) + term(2,j)*rediv2
                        hessz(j,iv) = hessz(j,iv) + term(3,j)*rediv2
                        hessx(j,kv) = hessx(j,kv) - term(1,j)*redivkv
                        hessy(j,kv) = hessy(j,kv) - term(2,j)*redivkv
                        hessz(j,kv) = hessz(j,kv) - term(3,j)*redivkv
                     end do
                  end if
               end if
            end if
         end do
      end do
      return
      end
