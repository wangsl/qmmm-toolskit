 
C     $Id$
c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder   ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole2  --  multipole & polarization Hessian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole2" calculates second derivatives of the multipole
c     and dipole polarization energy for a single atom at a time
c
c     note that since polarization effects are many body, it is not
c     really correct to neglect interactions with atoms that are not
c     directly involved as the multipole site or in axis definitions;
c     however, other sites are neglected in this version via the
c     passed atom list to quickly get approximate Hessian values;
c     to get exact values, "list" should include all multipole sites
c     and axis-defining atoms
c
c     also, the "reinduce" flag controls whether the induced dipoles
c     are recomputed every time an atom is moved during computation
c     of the numerical Hessian resulting in a further approximation;
c     setting the flag to "true" produces a much slower calculation,
c     but can greatly aid convergence of minimizations, etc.
c
c
      subroutine empole2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'hessn.i'
      include 'mpole.i'
      integer i,j,k
      integer nlist
      integer list(maxatm)
      real*8 eps,old
      real*8 d0(3,maxatm)
      logical reinduce
c
c
c     set the default stepsize and flag for induced dipoles
c
      eps = 1.0d-7
      reinduce = .false.
      if (n .le. 50)  reinduce = .true.
c
c     find the multipole definitions involving the current atom;
c     results in a faster but approximate Hessian calculation
c
      nlist = 0
      do k = 1, npole
         if (ipole(k).eq.i .or. zaxis(k).eq.i .or. xaxis(k).eq.i) then
            nlist = nlist + 1
            list(nlist) = k
         end if
      end do
c
c     get multipole first derivatives for the base structure
c
      call empole2a (nlist,list,reinduce)
      do k = 1, n
         do j = 1, 3
            d0(j,k) = dem(j,k) + dep(j,k)
         end do
      end do
c
c     find numerical x-components via perturbed structures
c
      old = x(i)
      x(i) = x(i) + eps
      call empole2a (nlist,list,reinduce)
      x(i) = old
      do k = 1, n
         do j = 1, 3
            hessx(j,k) = hessx(j,k) + (dem(j,k)+dep(j,k)-d0(j,k))/eps
         end do
      end do
c
c     find numerical y-components via perturbed structures
c
      old = y(i)
      y(i) = y(i) + eps
      call empole2a (nlist,list,reinduce)
      y(i) = old
      do k = 1, n
         do j = 1, 3
            hessy(j,k) = hessy(j,k) + (dem(j,k)+dep(j,k)-d0(j,k))/eps
         end do
      end do
c
c     find numerical z-components via perturbed structures
c
      old = z(i)
      z(i) = z(i) + eps
      call empole2a (nlist,list,reinduce)
      z(i) = old
      do k = 1, n
         do j = 1, 3
            hessz(j,k) = hessz(j,k) + (dem(j,k)+dep(j,k)-d0(j,k))/eps
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine empole2a  --  mpole & polar Hessian; numerical  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "empole2a" computes multipole and dipole polarization first
c     derivatives for a single atom with respect to Cartesian
c     coordinates; used to get finite difference second derivatives
c
c
      subroutine empole2a (nlist,list,reinduce)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'cutoff.i'
      include 'deriv.i'
      include 'group.i'
      include 'molcul.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k
      integer ii,jj
      integer ix,iz
      integer jx,jz
      integer iii,nlist
      integer list(maxatm)
      real*8 f,fgrp
      real*8 pterm,damp
      real*8 scale3,scale5
      real*8 scale7,scale9
      real*8 scale3i,scale5i
      real*8 scale7i
      real*8 xji,yji,zji
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 ci,di(3),qi(9)
      real*8 cj,dj(3),qj(9)
      real*8 fridmp(3),findmp(3)
      real*8 ftm2(3),ftm2i(3)
      real*8 ttm2(3),ttm3(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 dixdj(3),dixuj(3)
      real*8 djxui(3),fdir(3)
      real*8 uixr(3),ujxr(3)
      real*8 uixqjr(3),ujxqir(3)
      real*8 qiuj(3),qjui(3)
      real*8 rxqiuj(3),rxqjui(3)
      real*8 qidj(3),qjdi(3)
      real*8 qir(3),qjr(3)
      real*8 qiqjr(3),qjqir(3)
      real*8 qixqj(3),rxqir(3)
      real*8 dixr(3),djxr(3)
      real*8 dixqjr(3),djxqir(3)
      real*8 rxqjr(3),qjrxqir(3)
      real*8 rxqijr(3),rxqjir(3)
      real*8 rxqidj(3),rxqjdi(3)
      real*8 dscale3(3),dscale5(3)
      real*8 dscale7(3)
      real*8 sc(10),sci(10)
      real*8 gl(0:8),gli(0:8)
      real*8 gf(0:8),gfi(0:8)
      real*8 gfd(0:8),gti(0:8)
      real*8 mscale(maxatm)
      real*8 pscale(maxatm)
      real*8 uscale(maxatm)
      real*8 ftm1(3,maxatm)
      real*8 ftm1i(3,maxatm)
      real*8 ttm1(3,maxatm)
      real*8 ttm1i(3,maxatm)
      real*8 trq(3,maxatm)
      real*8 trqi(3,maxatm)
      logical proceed
      logical iuse,juse
      logical reinduce
c
c
c     zero out the multipole and polarization first derivatives
c
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
            dep(j,i) = 0.0d0
         end do
      end do
c
c     set conversion factor, cutoff and scaling coefficients
c
      if (nlist .eq. 0)  return
      f = electric / dielec
      call switch ('MPOLE')
c
c     check the sign of multipole components at chiral sites
c
      if (reinduce) then
         call chkpole
c
c     rotate the multipole components into the global frame
c
         call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
         call induce
      end if
c
c     zero out local accumulation arrays for derivatives
c
      do i = 1, n
         trq(1,i) = 0.0d0
         trq(2,i) = 0.0d0
         trq(3,i) = 0.0d0
         trqi(1,i) = 0.0d0
         trqi(2,i) = 0.0d0
         trqi(3,i) = 0.0d0
      end do
c
c     initialise temporary force and torque accumulators
c
      do i = 1, npole
         do k = 1, 3
            ftm1(k,i) = 0.0d0
            ftm1i(k,i) = 0.0d0
            ttm1(k,i) = 0.0d0
            ttm1i(k,i) = 0.0d0
         end do
      end do
c
c     set scale factors for permanent multipole and induced terms
c
      do iii = 1, nlist
         i = list(iii)
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         ci = rpole(1,i)
         di(1) = rpole(2,i)
         di(2) = rpole(3,i)
         di(3) = rpole(4,i)
         qi(1) = rpole(5,i)
         qi(2) = rpole(6,i)
         qi(3) = rpole(7,i)
         qi(4) = rpole(8,i)
         qi(5) = rpole(9,i)
         qi(6) = rpole(10,i)
         qi(7) = rpole(11,i)
         qi(8) = rpole(12,i)
         qi(9) = rpole(13,i)
         iuse = (use(ii) .or. use(iz) .or. use(ix))
         do j = 1, n
            mscale(j) = 1.0d0
            pscale(j) = 1.0d0
            uscale(j) = 1.0d0
         end do
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii)) 
     &            pscale(i14(j,ii)) = 0.5d0 * pscale(i14(j,ii))
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            uscale(ip14(j,ii)) = u4scale
         end do
         do j = i+1, npole
            jj = ipole(j)
            jz = zaxis(j)
            jx = xaxis(j)
            juse = (use(jj) .or. use(jz) .or. use(jx))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,jj,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (iuse .or. juse)
            if (.not. proceed)  goto 10
            cj = rpole(1,j)
            dj(1) = rpole(2,j)
            dj(2) = rpole(3,j)
            dj(3) = rpole(4,j)
            qj(1) = rpole(5,j)
            qj(2) = rpole(6,j)
            qj(3) = rpole(7,j)
            qj(4) = rpole(8,j)
            qj(5) = rpole(9,j)
            qj(6) = rpole(10,j)
            qj(7) = rpole(11,j)
            qj(8) = rpole(12,j)
            qj(9) = rpole(13,j)
            xji = x(jj) - x(ii)
            yji = y(jj) - y(ii)
            zji = z(jj) - z(ii)
            if (use_image)  call image (xji,yji,zji,0)
            r2 = xji*xji + yji*yji + zji*zji
            if (r2 .le. off2) then
               r = sqrt(r2)
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               scale9 = 1.0d0
               do k = 1, 3
                  dscale3(k) = 0.0d0
                  dscale5(k) = 0.0d0
                  dscale7(k) = 0.0d0
               end do
c
c     apply Thole polarization damping to scale factors
c
               pterm = pdamp(i) * pdamp(j)
               if (pterm .ne. 0.0d0) then
                  damp = -pgamma * (r/pterm)**3
                  if (damp .gt. -50.0d0) then
                     scale3 = 1.0d0 - exp(damp)
                     scale5 = 1.0d0 - (1.0d0-damp)*exp(damp)
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *exp(damp)
                     scale9 = 1.0d0 - (1.0d0-damp+(18.0d0*damp**2
     &                                 -9.0d0*damp**3)/35.0d0)*exp(damp)
                     dscale3(1) = -3.0d0*damp*exp(damp) * xji/r2
                     dscale3(2) = -3.0d0*damp*exp(damp) * yji/r2
                     dscale3(3) = -3.0d0*damp*exp(damp) * zji/r2
                     dscale5(1) = -damp * dscale3(1)
                     dscale5(2) = -damp * dscale3(2)
                     dscale5(3) = -damp * dscale3(3)
                     dscale7(1) = (-0.2d0-0.6d0*damp) * dscale5(1)
                     dscale7(2) = (-0.2d0-0.6d0*damp) * dscale5(2)
                     dscale7(3) = (-0.2d0-0.6d0*damp) * dscale5(3)
                  end if
               end if
               scale3i = scale3 * uscale(jj)
               scale5i = scale5 * uscale(jj)
               scale7i = scale7 * uscale(jj)
               scale3 = scale3 * pscale(jj)
               scale5 = scale5 * pscale(jj)
               scale7 = scale7 * pscale(jj)
               scale9 = scale9 * pscale(jj)
c
c     construct auxiliary vectors
c
               dixdj(1) = di(2)*dj(3) - di(3)*dj(2)
               dixdj(2) = di(3)*dj(1) - di(1)*dj(3)
               dixdj(3) = di(1)*dj(2) - di(2)*dj(1)
               dixuj(1) = di(2)*uind(3,j) - di(3)*uind(2,j)
               dixuj(2) = di(3)*uind(1,j) - di(1)*uind(3,j)
               dixuj(3) = di(1)*uind(2,j) - di(2)*uind(1,j)
               djxui(1) = dj(2)*uind(3,i) - dj(3)*uind(2,i)
               djxui(2) = dj(3)*uind(1,i) - dj(1)*uind(3,i)
               djxui(3) = dj(1)*uind(2,i) - dj(2)*uind(1,i)
               dixr(1) = di(2)*zji - di(3)*yji
               dixr(2) = di(3)*xji - di(1)*zji
               dixr(3) = di(1)*yji - di(2)*xji
               djxr(1) = dj(2)*zji - dj(3)*yji
               djxr(2) = dj(3)*xji - dj(1)*zji
               djxr(3) = dj(1)*yji - dj(2)*xji
               uixr(1) = uind(2,i)*zji - uind(3,i)*yji
               uixr(2) = uind(3,i)*xji - uind(1,i)*zji
               uixr(3) = uind(1,i)*yji - uind(2,i)*xji
               ujxr(1) = uind(2,j)*zji - uind(3,j)*yji
               ujxr(2) = uind(3,j)*xji - uind(1,j)*zji
               ujxr(3) = uind(1,j)*yji - uind(2,j)*xji
               qir(1) = qi(1)*xji + qi(4)*yji + qi(7)*zji
               qir(2) = qi(2)*xji + qi(5)*yji + qi(8)*zji
               qir(3) = qi(3)*xji + qi(6)*yji + qi(9)*zji
               qjr(1) = qj(1)*xji + qj(4)*yji + qj(7)*zji
               qjr(2) = qj(2)*xji + qj(5)*yji + qj(8)*zji
               qjr(3) = qj(3)*xji + qj(6)*yji + qj(9)*zji
               qiqjr(1) = qi(1)*qjr(1) + qi(4)*qjr(2) + qi(7)*qjr(3)
               qiqjr(2) = qi(2)*qjr(1) + qi(5)*qjr(2) + qi(8)*qjr(3)
               qiqjr(3) = qi(3)*qjr(1) + qi(6)*qjr(2) + qi(9)*qjr(3)
               qjqir(1) = qj(1)*qir(1) + qj(4)*qir(2) + qj(7)*qir(3)
               qjqir(2) = qj(2)*qir(1) + qj(5)*qir(2) + qj(8)*qir(3)
               qjqir(3) = qj(3)*qir(1) + qj(6)*qir(2) + qj(9)*qir(3)
               qixqj(1) = qi(2)*qj(3) + qi(5)*qj(6) + qi(8)*qj(9)
     &                       - qi(3)*qj(2) - qi(6)*qj(5) - qi(9)*qj(8)
               qixqj(2) = qi(3)*qj(1) + qi(6)*qj(4) + qi(9)*qj(7)
     &                       - qi(1)*qj(3) - qi(4)*qj(6) - qi(7)*qj(9)
               qixqj(3) = qi(1)*qj(2) + qi(4)*qj(5) + qi(7)*qj(8)
     &                       - qi(2)*qj(1) - qi(5)*qj(4) - qi(8)*qj(7)
               rxqir(1) = yji*qir(3) - zji*qir(2)
               rxqir(2) = zji*qir(1) - xji*qir(3)
               rxqir(3) = xji*qir(2) - yji*qir(1)
               rxqjr(1) = yji*qjr(3) - zji*qjr(2)
               rxqjr(2) = zji*qjr(1) - xji*qjr(3)
               rxqjr(3) = xji*qjr(2) - yji*qjr(1)
               rxqijr(1) = yji*qiqjr(3) - zji*qiqjr(2)
               rxqijr(2) = zji*qiqjr(1) - xji*qiqjr(3)
               rxqijr(3) = xji*qiqjr(2) - yji*qiqjr(1)
               rxqjir(1) = yji*qjqir(3) - zji*qjqir(2)
               rxqjir(2) = zji*qjqir(1) - xji*qjqir(3)
               rxqjir(3) = xji*qjqir(2) - yji*qjqir(1)
               qjrxqir(1) = qjr(2)*qir(3) - qjr(3)*qir(2)
               qjrxqir(2) = qjr(3)*qir(1) - qjr(1)*qir(3)
               qjrxqir(3) = qjr(1)*qir(2) - qjr(2)*qir(1)
               qidj(1) = qi(1)*dj(1) + qi(4)*dj(2) + qi(7)*dj(3)
               qidj(2) = qi(2)*dj(1) + qi(5)*dj(2) + qi(8)*dj(3)
               qidj(3) = qi(3)*dj(1) + qi(6)*dj(2) + qi(9)*dj(3)
               qjdi(1) = qj(1)*di(1) + qj(4)*di(2) + qj(7)*di(3)
               qjdi(2) = qj(2)*di(1) + qj(5)*di(2) + qj(8)*di(3)
               qjdi(3) = qj(3)*di(1) + qj(6)*di(2) + qj(9)*di(3)
               qiuj(1) = qi(1)*uind(1,j)+qi(4)*uind(2,j)+qi(7)*uind(3,j)
               qiuj(2) = qi(2)*uind(1,j)+qi(5)*uind(2,j)+qi(8)*uind(3,j)
               qiuj(3) = qi(3)*uind(1,j)+qi(6)*uind(2,j)+qi(9)*uind(3,j)
               qjui(1) = qj(1)*uind(1,i)+qj(4)*uind(2,i)+qj(7)*uind(3,i)
               qjui(2) = qj(2)*uind(1,i)+qj(5)*uind(2,i)+qj(8)*uind(3,i)
               qjui(3) = qj(3)*uind(1,i)+qj(6)*uind(2,i)+qj(9)*uind(3,i)
               dixqjr(1) = di(2)*qjr(3) - di(3)*qjr(2)
               dixqjr(2) = di(3)*qjr(1) - di(1)*qjr(3)
               dixqjr(3) = di(1)*qjr(2) - di(2)*qjr(1)
               djxqir(1) = dj(2)*qir(3) - dj(3)*qir(2)
               djxqir(2) = dj(3)*qir(1) - dj(1)*qir(3)
               djxqir(3) = dj(1)*qir(2) - dj(2)*qir(1)
               uixqjr(1) = uind(2,i)*qjr(3) - uind(3,i)*qjr(2)
               uixqjr(2) = uind(3,i)*qjr(1) - uind(1,i)*qjr(3)
               uixqjr(3) = uind(1,i)*qjr(2) - uind(2,i)*qjr(1)
               ujxqir(1) = uind(2,j)*qir(3) - uind(3,j)*qir(2)
               ujxqir(2) = uind(3,j)*qir(1) - uind(1,j)*qir(3)
               ujxqir(3) = uind(1,j)*qir(2) - uind(2,j)*qir(1)
               rxqidj(1) = yji*qidj(3) - zji*qidj(2)
               rxqidj(2) = zji*qidj(1) - xji*qidj(3)
               rxqidj(3) = xji*qidj(2) - yji*qidj(1)
               rxqjdi(1) = yji*qjdi(3) - zji*qjdi(2)
               rxqjdi(2) = zji*qjdi(1) - xji*qjdi(3)
               rxqjdi(3) = xji*qjdi(2) - yji*qjdi(1)
               rxqiuj(1) = yji*qiuj(3) - zji*qiuj(2)
               rxqiuj(2) = zji*qiuj(1) - xji*qiuj(3)
               rxqiuj(3) = xji*qiuj(2) - yji*qiuj(1)
               rxqjui(1) = yji*qjui(3) - zji*qjui(2)
               rxqjui(2) = zji*qjui(1) - xji*qjui(3)
               rxqjui(3) = xji*qjui(2) - yji*qjui(1)
c
c     get intermediate variables for permanent energy terms
c
               sc(2) = di(1)*dj(1) + di(2)*dj(2) + di(3)*dj(3)
               sc(3) = di(1)*xji + di(2)*yji + di(3)*zji
               sc(4) = dj(1)*xji + dj(2)*yji + dj(3)*zji
               sc(5) = qir(1)*xji + qir(2)*yji + qir(3)*zji
               sc(6) = qjr(1)*xji + qjr(2)*yji + qjr(3)*zji
               sc(7) = qir(1)*dj(1) + qir(2)*dj(2) + qir(3)*dj(3)
               sc(8) = qjr(1)*di(1) + qjr(2)*di(2) + qjr(3)*di(3)
               sc(9) = qir(1)*qjr(1) + qir(2)*qjr(2) + qir(3)*qjr(3)
               sc(10) = qi(1)*qj(1) + qi(2)*qj(2) + qi(3)*qj(3)
     &                     + qi(4)*qj(4) + qi(5)*qj(5) + qi(6)*qj(6)
     &                     + qi(7)*qj(7) + qi(8)*qj(8) + qi(9)*qj(9)
c
c     get intermediate variables for induction energy terms
c
               sci(1) = uind(1,i)*dj(1) + uind(2,i)*dj(2)
     &                     + uind(3,i)*dj(3) + di(1)*uind(1,j)
     &                     + di(2)*uind(2,j) + di(3)*uind(3,j)
               sci(2) = uind(1,i)*uind(1,j) + uind(2,i)*uind(2,j)
     &                     + uind(3,i)*uind(3,j)
               sci(3) = uind(1,i)*xji + uind(2,i)*yji + uind(3,i)*zji
               sci(4) = uind(1,j)*xji + uind(2,j)*yji + uind(3,j)*zji
               sci(7) = qir(1)*uind(1,j) + qir(2)*uind(2,j)
     &                     + qir(3)*uind(3,j)
               sci(8) = qjr(1)*uind(1,i) + qjr(2)*uind(2,i)
     &                     + qjr(3)*uind(3,i)
c
c     get the induced-induced derivative terms
c
               findmp(1) = uscale(jj) * (sci(2)*rr3*dscale3(1)
     &                            - sci(3)*sci(4)*rr5*dscale5(1))
               findmp(2) = uscale(jj) * (sci(2)*rr3*dscale3(2)
     &                            - sci(3)*sci(4)*rr5*dscale5(2))
               findmp(3) = uscale(jj) * (sci(2)*rr3*dscale3(3)
     &                            - sci(3)*sci(4)*rr5*dscale5(3))
c
c     calculate the gl functions for potential energy
c
               gl(0) = ci*cj
               gl(1) = cj*sc(3) - ci*sc(4)
               gl(2) = ci*sc(6) + cj*sc(5) - sc(3)*sc(4)
               gl(3) = sc(3)*sc(6) - sc(4)*sc(5)
               gl(4) = sc(5)*sc(6)
               gl(6) = sc(2)
               gl(7) = 2.0d0 * (sc(7)-sc(8))
               gl(8) = 2.0d0 * sc(10)
               gl(5) = -4.0d0 * sc(9)
               gli(1) = cj*sci(3) - ci*sci(4)
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4) - sci(3)*sci(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1) + sci(2)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
c
c    get the permanent and induced energy for this interaction
c
               fridmp(1) = (rr3*(gli(1)+gli(6)-sci(2))*dscale3(1)
     &                   + rr5*(gli(2)+gli(7)+sci(3)*sci(4))*dscale5(1)
     &                   + rr7*gli(3)*dscale7(1))*pscale(jj)
               fridmp(2) = (rr3*(gli(1)+gli(6)-sci(2))*dscale3(2)
     &                   + rr5*(gli(2)+gli(7)+sci(3)*sci(4))*dscale5(2)
     &                   + rr7*gli(3)*dscale7(2))*pscale(jj)
               fridmp(3) = (rr3*(gli(1)+gli(6)-sci(2))*dscale3(3)
     &                   + rr5*(gli(2)+gli(7)+sci(3)*sci(4))*dscale5(3)
     &                   + rr7*gli(3)*dscale7(3))*pscale(jj)
c
c     intermediate variables for the permanent multipole terms
c
               gf(1) = rr3*gl(0) + rr5*(gl(1)+gl(6))
     &                    + rr7*(gl(2)+gl(7)+gl(8))
     &                    + rr9*(gl(3)+gl(5)) + rr11*gl(4)
               gf(2) = -cj*rr3 + sc(4)*rr5 - sc(6)*rr7
               gf(3) = ci*rr3 + sc(3)*rr5 + sc(5)*rr7
               gf(4) = 2.0d0 * rr5
               gf(5) = 2.0d0 * (-cj*rr5+sc(4)*rr7-sc(6)*rr9)
               gf(6) = 2.0d0 * (-ci*rr5-sc(3)*rr7-sc(5)*rr9)
               gf(7) = 4.0d0 * rr7
c
c     intermediate variables for the induced-permanent terms
c
               gfi(1) = rr5*((gli(1)+sci(1))*scale3+sci(2)*scale3i)
     &                     + rr7*(gli(7)-sc(3)*sci(4)
     &                               -sc(4)*sci(3))*scale5
     &                     - rr7*sci(3)*sci(4)*scale5i
     &                     + rr9*gli(3)*scale7
               gfi(2) = -rr3*cj*scale3
     &                     + rr5*(sc(4)*scale5+sci(4)*scale5i)
     &                     - rr7*sc(6)*scale7
               gfi(3) = rr3*ci*scale3
     &                     + rr5*(sc(3)*scale5+sci(3)*scale5i)
     &                     + rr7*sc(5)*scale7
               gfi(4) = 2.0d0 * rr5 * scale5
               gfi(5) = 2.0d0 * rr7 * sci(4) * scale7
               gfi(6) = -2.0d0 * rr7 * sci(3) * scale7
c
c     get the permanent force
c
               ftm2(1) = gf(1)*xji + gf(2)*di(1) + gf(3)*dj(1)
     &                      + gf(4)*(qjdi(1)-qidj(1)) + gf(5)*qir(1)
     &                      + gf(6)*qjr(1) + gf(7)*(qiqjr(1)+qjqir(1))
               ftm2(2) = gf(1)*yji + gf(2)*di(2) + gf(3)*dj(2)
     &                      + gf(4)*(qjdi(2)-qidj(2)) + gf(5)*qir(2)
     &                      + gf(6)*qjr(2) + gf(7)*(qiqjr(2)+qjqir(2))
               ftm2(3) = gf(1)*zji + gf(2)*di(3) + gf(3)*dj(3)
     &                      + gf(4)*(qjdi(3)-qidj(3)) + gf(5)*qir(3)
     &                      + gf(6)*qjr(3) + gf(7)*(qiqjr(3)+qjqir(3))
c
c     get the induced force
c
               ftm2i(1) = gfi(1)*xji + gfi(2)*uind(1,i)
     &           + gfi(3)*uind(1,j) + sci(4)*rr5*di(1)*scale5
     &           + sci(3)*rr5*dj(1)*scale5 + gfi(4)*(qjui(1)-qiuj(1))
     &           + gfi(5)*qir(1) + gfi(6)*qjr(1)
               ftm2i(2) = gfi(1)*yji + gfi(2)*uind(2,i)
     &           + gfi(3)*uind(2,j) + sci(4)*rr5*di(2)*scale5
     &           + sci(3)*rr5*dj(2)*scale5 + gfi(4)*(qjui(2)-qiuj(2))
     &           + gfi(5)*qir(2) + gfi(6)*qjr(2)
               ftm2i(3) = gfi(1)*zji + gfi(2)*uind(3,i)
     &           + gfi(3)*uind(3,j) + sci(4)*rr5*di(3)*scale5
     &           + sci(3)*rr5*dj(3)*scale5 + gfi(4)*(qjui(3)-qiuj(3))
     &           + gfi(5)*qir(3) + gfi(6)*qjr(3)
c
c     handle of scaling for partially excluded interactions
c
               ftm2(1) = mscale(jj) * ftm2(1)
               ftm2(2) = mscale(jj) * ftm2(2)
               ftm2(3) = mscale(jj) * ftm2(3)
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)
c
c     correction to convert mutual to direct polarization force
c
               if (poltyp .eq. 'DIRECT') then
                  gfd(1) = rr5*sci(2)*scale3i
     &                        - rr7*sci(3)*sci(4)*scale5i
                  gfd(2) = rr5 * sci(4) * scale5i
                  gfd(3) = rr5 * sci(3) * scale5i
                  fdir(1) = gfd(1)*xji + gfd(2)*uind(1,i)
     &                         + gfd(3)*uind(1,j)
                  fdir(2) = gfd(1)*yji + gfd(2)*uind(2,i)
     &                         + gfd(3)*uind(2,j)
                  fdir(3) = gfd(1)*zji + gfd(2)*uind(3,i)
     &                         + gfd(3)*uind(3,j)
                  ftm2i(1) = ftm2i(1) - fdir(1) + findmp(1)
                  ftm2i(2) = ftm2i(2) - fdir(2) + findmp(2)
                  ftm2i(3) = ftm2i(3) - fdir(3) + findmp(3)
               end if
c
c     now perform the torque calculation
c     intermediate terms for torque between multipoles i and j
c
               gti(2) = sci(4) * rr5 * scale5
               gti(3) = sci(3) * rr5 * scale5
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
c
c     get the permanent interaction torques
c
               ttm2(1) = -rr3*dixdj(1) + gf(2)*dixr(1)-gf(5)*rxqir(1)
     &           + gf(4)*(dixqjr(1)+djxqir(1)+rxqidj(1)-2.0d0*qixqj(1))
     &           - gf(7)*(rxqijr(1)+qjrxqir(1))
               ttm2(2) = -rr3*dixdj(2) + gf(2)*dixr(2)-gf(5)*rxqir(2)
     &           + gf(4)*(dixqjr(2)+djxqir(2)+rxqidj(2)-2.0d0*qixqj(2))
     &           - gf(7)*(rxqijr(2)+qjrxqir(2))
               ttm2(3) = -rr3*dixdj(3) + gf(2)*dixr(3)-gf(5)*rxqir(3)
     &           + gf(4)*(dixqjr(3)+djxqir(3)+rxqidj(3)-2.0d0*qixqj(3))
     &           - gf(7)*(rxqijr(3)+qjrxqir(3))
               ttm3(1) = rr3*dixdj(1) + gf(3)*djxr(1) -gf(6)*rxqjr(1)
     &           - gf(4)*(dixqjr(1)+djxqir(1)+rxqjdi(1)-2.0d0*qixqj(1))
     &           - gf(7)*(rxqjir(1)-qjrxqir(1))
               ttm3(2) = rr3*dixdj(2) + gf(3)*djxr(2) -gf(6)*rxqjr(2)
     &           - gf(4)*(dixqjr(2)+djxqir(2)+rxqjdi(2)-2.0d0*qixqj(2))
     &           - gf(7)*(rxqjir(2)-qjrxqir(2))
               ttm3(3) = rr3*dixdj(3) + gf(3)*djxr(3) -gf(6)*rxqjr(3)
     &           - gf(4)*(dixqjr(3)+djxqir(3)+rxqjdi(3)-2.0d0*qixqj(3))
     &           - gf(7)*(rxqjir(3)-qjrxqir(3))
c
c     induced-torque on site i is due to the field of induced
c     dipole on site j interacting with permanent D and Q on site i;
c     di x fuindj + qi x fuindj
c
               ttm2i(1) = -rr3*dixuj(1)*scale3 + gti(2)*dixr(1)
     &              + gti(4)*(ujxqir(1)+rxqiuj(1)) - gti(5)*rxqir(1)
               ttm2i(2) = -rr3*dixuj(2)*scale3 + gti(2)*dixr(2)
     &              + gti(4)*(ujxqir(2)+rxqiuj(2)) - gti(5)*rxqir(2)
               ttm2i(3) = -rr3*dixuj(3)*scale3 + gti(2)*dixr(3)
     &              + gti(4)*(ujxqir(3)+rxqiuj(3)) - gti(5)*rxqir(3)
               ttm3i(1) = -rr3*djxui(1)*scale3 + gti(3)*djxr(1)
     &              - gti(4)*(uixqjr(1)+rxqjui(1)) - gti(6)*rxqjr(1)
               ttm3i(2) = -rr3*djxui(2)*scale3 + gti(3)*djxr(2)
     &              - gti(4)*(uixqjr(2)+rxqjui(2)) - gti(6)*rxqjr(2)
               ttm3i(3) = -rr3*djxui(3)*scale3 + gti(3)*djxr(3)
     &              - gti(4)*(uixqjr(3)+rxqjui(3)) - gti(6)*rxqjr(3)
c
c     handle the case where scaling is used
c
               ttm2(1) = ttm2(1) * mscale(jj)
               ttm2(2) = ttm2(2) * mscale(jj)
               ttm2(3) = ttm2(3) * mscale(jj)
               ttm3(1) = ttm3(1) * mscale(jj)
               ttm3(2) = ttm3(2) * mscale(jj)
               ttm3(3) = ttm3(3) * mscale(jj)
c
c     update force and torque on site j
c
               ftm1(1,j) = ftm1(1,j) + ftm2(1)
               ftm1(2,j) = ftm1(2,j) + ftm2(2)
               ftm1(3,j) = ftm1(3,j) + ftm2(3)
               ttm1(1,j) = ttm1(1,j) + ttm3(1)
               ttm1(2,j) = ttm1(2,j) + ttm3(2)
               ttm1(3,j) = ttm1(3,j) + ttm3(3)
               ftm1i(1,j) = ftm1i(1,j) + ftm2i(1)
               ftm1i(2,j) = ftm1i(2,j) + ftm2i(2)
               ftm1i(3,j) = ftm1i(3,j) + ftm2i(3)
               ttm1i(1,j) = ttm1i(1,j) + ttm3i(1)
               ttm1i(2,j) = ttm1i(2,j) + ttm3i(2)
               ttm1i(3,j) = ttm1i(3,j) + ttm3i(3)
c
c     update force and torque on site i
c
               ftm1(1,i) = ftm1(1,i) - ftm2(1)
               ftm1(2,i) = ftm1(2,i) - ftm2(2)
               ftm1(3,i) = ftm1(3,i) - ftm2(3)
               ttm1(1,i) = ttm1(1,i) + ttm2(1)
               ttm1(2,i) = ttm1(2,i) + ttm2(2)
               ttm1(3,i) = ttm1(3,i) + ttm2(3)
               ftm1i(1,i) = ftm1i(1,i) - ftm2i(1)
               ftm1i(2,i) = ftm1i(2,i) - ftm2i(2)
               ftm1i(3,i) = ftm1i(3,i) - ftm2i(3)
               ttm1i(1,i) = ttm1i(1,i) + ttm2i(1)
               ttm1i(2,i) = ttm1i(2,i) + ttm2i(2)
               ttm1i(3,i) = ttm1i(3,i) + ttm2i(3)
            end if
   10       continue
         end do
      end do
c
c     increment the total forces and torques
c
      do i = 1, npole
         ii = ipole(i)
         dem(1,ii) = dem(1,ii) - f*ftm1(1,i)
         dem(2,ii) = dem(2,ii) - f*ftm1(2,i)
         dem(3,ii) = dem(3,ii) - f*ftm1(3,i)
         dep(1,ii) = dep(1,ii) - f*ftm1i(1,i)
         dep(2,ii) = dep(2,ii) - f*ftm1i(2,i)
         dep(3,ii) = dep(3,ii) - f*ftm1i(3,i)
         trq(1,ii) = trq(1,ii) + f*ttm1(1,i)
         trq(2,ii) = trq(2,ii) + f*ttm1(2,i)
         trq(3,ii) = trq(3,ii) + f*ttm1(3,i)
         trqi(1,ii) = trqi(1,ii) + f*ttm1i(1,i)
         trqi(2,ii) = trqi(2,ii) + f*ttm1i(2,i)
         trqi(3,ii) = trqi(3,ii) + f*ttm1i(3,i)
      end do
      call torque (trq,dem)
      call torque (trqi,dep)
      return
      end
