 
C     $Id: empole1.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder   ##
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole1  --  mpole/polar energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole1" calculates the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c
c
      subroutine empole1
      implicit none
      include 'cutoff.i'
c
c
c     choose pairwise double loop or Ewald summation version
c
      if (use_ewald) then
         call empole1b
      else
         call empole1a
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empole1a  --  double loop multipole derivatives  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empole1a" calculates the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c     using a pairwise double loop
c
c
      subroutine empole1a
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'cutoff.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'molcul.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,k
      integer ii,jj
      integer ix,iz
      integer jx,jz
      real*8 e,ei,f,fgrp
      real*8 pterm,damp,gfd
      real*8 scale3,scale5
      real*8 scale7,scale9
      real*8 psc3,psc5,psc7,psc9
      real*8 dsc3,dsc5,dsc7,dsc9
      real*8 scale3i,scale5i
      real*8 scale7i
      real*8 xji,yji,zji
      real*8 xix,yix,zix
      real*8 xiz,yiz,ziz
      real*8 xjx,yjx,zjx
      real*8 xjz,yjz,zjz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 frcxi(3),frczi(3)
      real*8 frcxj(3),frczj(3)
      real*8 ci,di(3),qi(9)
      real*8 cj,dj(3),qj(9)
      real*8 fridmp(3),findmp(3)
      real*8 ftm2(3),ftm2i(3)
      real*8 ttm2(3),ttm3(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 dixdj(3),fdir(3)
      real*8 dixuj(3),djxui(3)
      real*8 dixujp(3),djxuip(3)
      real*8 uixqjr(3),ujxqir(3)
      real*8 uixqjrp(3),ujxqirp(3)
      real*8 qiuj(3),qjui(3)
      real*8 qiujp(3),qjuip(3)
      real*8 rxqiuj(3),rxqjui(3)
      real*8 rxqiujp(3),rxqjuip(3)
      real*8 qidj(3),qjdi(3)
      real*8 qir(3),qjr(3)
      real*8 qiqjr(3),qjqir(3)
      real*8 qixqj(3),rxqir(3)
      real*8 dixr(3),djxr(3)
      real*8 dixqjr(3),djxqir(3)
      real*8 rxqjr(3),qjrxqir(3)
      real*8 rxqijr(3),rxqjir(3)
      real*8 rxqidj(3),rxqjdi(3)
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3),tq(3)
      real*8 gl(0:8),gli(7),glip(7)
      real*8 sc(10),sci(8),scip(8)
      real*8 gf(7),gfi(6),gti(6)
      real*8 mscale(maxatm)
      real*8 pscale(maxatm)
      real*8 dscale(maxatm)
      real*8 uscale(maxatm)
      real*8 ftm1(3,maxatm)
      real*8 ftm1i(3,maxatm)
      real*8 ttm1(3,maxatm)
      real*8 ttm1i(3,maxatm)
      real*8 trq(3,maxatm)
      real*8 trqi(3,maxatm)
      logical proceed,iuse,juse
c
c
c     set conversion factor, cutoff and scaling coefficients
c
      if (npole .eq. 0)  return
      f = electric / dielec
      call switch ('MPOLE')
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
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
      do i = 1, npole-1
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
            dscale(j) = 1.0d0
            uscale(j) = 1.0d0
            pscale(j) = 1.0d0
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
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
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
                  ddsc3(k) = 0.0d0
                  ddsc5(k) = 0.0d0
                  ddsc7(k) = 0.0d0
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
                     ddsc3(1) = -3.0d0*damp*exp(damp) * xji/r2
                     ddsc3(2) = -3.0d0*damp*exp(damp) * yji/r2
                     ddsc3(3) = -3.0d0*damp*exp(damp) * zji/r2
                     ddsc5(1) = -damp * ddsc3(1)
                     ddsc5(2) = -damp * ddsc3(2)
                     ddsc5(3) = -damp * ddsc3(3)
                     ddsc7(1) = (-0.2d0-0.6d0*damp) * ddsc5(1)
                     ddsc7(2) = (-0.2d0-0.6d0*damp) * ddsc5(2)
                     ddsc7(3) = (-0.2d0-0.6d0*damp) * ddsc5(3)
                  end if
               end if
               scale3i = scale3 * uscale(jj)
               scale5i = scale5 * uscale(jj)
               scale7i = scale7 * uscale(jj)
               dsc3 = scale3 * dscale(jj)
               dsc5 = scale5 * dscale(jj)
               dsc7 = scale7 * dscale(jj)
               dsc9 = scale9 * dscale(jj)
               psc3 = scale3 * pscale(jj)
               psc5 = scale5 * pscale(jj)
               psc7 = scale7 * pscale(jj)
               psc9 = scale9 * pscale(jj)
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
               dixujp(1) = di(2)*uinp(3,j) - di(3)*uinp(2,j)
               dixujp(2) = di(3)*uinp(1,j) - di(1)*uinp(3,j)
               dixujp(3) = di(1)*uinp(2,j) - di(2)*uinp(1,j)
               djxuip(1) = dj(2)*uinp(3,i) - dj(3)*uinp(2,i)
               djxuip(2) = dj(3)*uinp(1,i) - dj(1)*uinp(3,i)
               djxuip(3) = dj(1)*uinp(2,i) - dj(2)*uinp(1,i)
               dixr(1) = di(2)*zji - di(3)*yji
               dixr(2) = di(3)*xji - di(1)*zji
               dixr(3) = di(1)*yji - di(2)*xji
               djxr(1) = dj(2)*zji - dj(3)*yji
               djxr(2) = dj(3)*xji - dj(1)*zji
               djxr(3) = dj(1)*yji - dj(2)*xji
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
               qiujp(1)= qi(1)*uinp(1,j)+qi(4)*uinp(2,j)+qi(7)*uinp(3,j)
               qiujp(2)= qi(2)*uinp(1,j)+qi(5)*uinp(2,j)+qi(8)*uinp(3,j)
               qiujp(3)= qi(3)*uinp(1,j)+qi(6)*uinp(2,j)+qi(9)*uinp(3,j)
               qjuip(1)= qj(1)*uinp(1,i)+qj(4)*uinp(2,i)+qj(7)*uinp(3,i)
               qjuip(2)= qj(2)*uinp(1,i)+qj(5)*uinp(2,i)+qj(8)*uinp(3,i)
               qjuip(3)= qj(3)*uinp(1,i)+qj(6)*uinp(2,i)+qj(9)*uinp(3,i)
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
               uixqjrp(1) = uinp(2,i)*qjr(3) - uinp(3,i)*qjr(2)
               uixqjrp(2) = uinp(3,i)*qjr(1) - uinp(1,i)*qjr(3)
               uixqjrp(3) = uinp(1,i)*qjr(2) - uinp(2,i)*qjr(1)
               ujxqirp(1) = uinp(2,j)*qir(3) - uinp(3,j)*qir(2)
               ujxqirp(2) = uinp(3,j)*qir(1) - uinp(1,j)*qir(3)
               ujxqirp(3) = uinp(1,j)*qir(2) - uinp(2,j)*qir(1)
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
               rxqiujp(1) = yji*qiujp(3) - zji*qiujp(2)
               rxqiujp(2) = zji*qiujp(1) - xji*qiujp(3)
               rxqiujp(3) = xji*qiujp(2) - yji*qiujp(1)
               rxqjuip(1) = yji*qjuip(3) - zji*qjuip(2)
               rxqjuip(2) = zji*qjuip(1) - xji*qjuip(3)
               rxqjuip(3) = xji*qjuip(2) - yji*qjuip(1)
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
               scip(1) = uinp(1,i)*dj(1) + uinp(2,i)*dj(2)
     &                      + uinp(3,i)*dj(3) + di(1)*uinp(1,j)
     &                      + di(2)*uinp(2,j) + di(3)*uinp(3,j)
               scip(2) = uind(1,i)*uinp(1,j)+uind(2,i)*uinp(2,j)
     &                      + uind(3,i)*uinp(3,j)+uinp(1,i)*uind(1,j)
     &                      + uinp(2,i)*uind(2,j)+uinp(3,i)*uind(3,j)
               scip(3) = uinp(1,i)*xji + uinp(2,i)*yji + uinp(3,i)*zji
               scip(4) = uinp(1,j)*xji + uinp(2,j)*yji + uinp(3,j)*zji
               scip(7) = qir(1)*uinp(1,j) + qir(2)*uinp(2,j)
     &                      + qir(3)*uinp(3,j)
               scip(8) = qjr(1)*uinp(1,i) + qjr(2)*uinp(2,i)
     &                      + qjr(3)*uinp(3,i)
c
c     get the induced-induced derivative terms
c
               findmp(1) = uscale(jj) * (scip(2)*rr3*ddsc3(1)
     &           - rr5*ddsc5(1)*(sci(3)*scip(4)+scip(3)*sci(4)))*0.5d0
               findmp(2) = uscale(jj) * (scip(2)*rr3*ddsc3(2)
     &           - rr5*ddsc5(2)*(sci(3)*scip(4)+scip(3)*sci(4)))*0.5d0
               findmp(3) = uscale(jj) * (scip(2)*rr3*ddsc3(3)
     &           - rr5*ddsc5(3)*(sci(3)*scip(4)+scip(3)*sci(4)))*0.5d0
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
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = cj*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c    get the permanent multipole and induced energies
c
               e = rr1*gl(0) + rr3*(gl(1)+gl(6))
     &                + rr5*(gl(2)+gl(7)+gl(8))
     &                + rr7*(gl(3)+gl(5)) + rr9*gl(4)
               ei = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &                   + rr5*(gli(2)+gli(7))*psc5
     &                   + rr7*gli(3)*psc7)
               e = f * mscale(jj) * e
               ei = f * ei
               em = em + e
               ep = ep + ei
c
c
c
               fridmp(1) = 0.5d0*(rr3*((gli(1)+gli(6))*pscale(jj)
     &           + (glip(1)+glip(6))*dscale(jj))*ddsc3(1)
     &           + rr5*((gli(2)+gli(7))*pscale(jj)
     &           + (glip(2)+glip(7))*dscale(jj))*ddsc5(1)
     &           + rr7*(gli(3)*pscale(jj)+glip(3)*dscale(jj))*ddsc7(1))
               fridmp(2) = 0.5d0*(rr3*((gli(1)+gli(6))*pscale(jj)
     &           + (glip(1)+glip(6))*dscale(jj))*ddsc3(2)
     &           + rr5*((gli(2)+gli(7))*pscale(jj)
     &           + (glip(2)+glip(7))*dscale(jj))*ddsc5(2)
     &           + rr7*(gli(3)*pscale(jj)+glip(3)*dscale(jj))*ddsc7(2))
               fridmp(3) = 0.5d0*(rr3*((gli(1)+gli(6))*pscale(jj)
     &           + (glip(1)+glip(6))*dscale(jj))*ddsc3(3)
     &           + rr5*((gli(2)+gli(7))*pscale(jj)
     &           + (glip(2)+glip(7))*dscale(jj))*ddsc5(3)
     &           + rr7*(gli(3)*pscale(jj)+glip(3)*dscale(jj))*ddsc7(3))
c
c     increment the total intermolecular energy
c
               if (molcule(ii) .ne. molcule(jj)) then
                  einter = einter + e + ei
               end if
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
               gfi(1) = rr5*0.5d0*((gli(1)+gli(6))*psc3
     &             + (glip(1)+glip(6))*dsc3 + scip(2)*scale3i)
     &             + rr7*((gli(7)+gli(2))*psc5 + (glip(7)+glip(2))*dsc5
     &             - (sci(3)*scip(4)+scip(3)*sci(4))*scale5i)*0.5d0
     &             + rr9*(gli(3)*psc7+glip(3)*dsc7)*0.5d0
               gfi(2) = -rr3*cj + rr5*sc(4) - rr7*sc(6)
               gfi(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfi(4) = 2.0d0 * rr5
               gfi(5) = rr7 * (sci(4)*psc7 + scip(4)*dsc7)
               gfi(6) = -rr7 * (sci(3)*psc7 + scip(3)*dsc7)
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
               ftm2i(1) = gfi(1)*xji + 0.5d0*
     &          (- rr3*cj*(uind(1,i)*psc3+uinp(1,i)*dsc3)
     &           + rr5*sc(4)*(uind(1,i)*psc5+uinp(1,i)*dsc5)
     &           - rr7*sc(6)*(uind(1,i)*psc7+uinp(1,i)*dsc7))
     &           +(rr3*ci*(uind(1,j)*psc3+uinp(1,j)*dsc3)
     &           + rr5*sc(3)*(uind(1,j)*psc5+uinp(1,j)*dsc5)
     &           + rr7*sc(5)*(uind(1,j)*psc7+uinp(1,j)*dsc7))*0.5d0
     &           + rr5*scale5i*(sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &           + sci(3)*uinp(1,j)+scip(3)*uind(1,j))*0.5d0
     &           + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &           + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dj(1)
     &           + 0.5d0*gfi(4)*((qjui(1)-qiuj(1))*psc5
     &           + (qjuip(1)-qiujp(1))*dsc5)
     &           + gfi(5)*qir(1) + gfi(6)*qjr(1)
               ftm2i(2) = gfi(1)*yji + 0.5d0*
     &          (- rr3*cj*(uind(2,i)*psc3+uinp(2,i)*dsc3)
     &           + rr5*sc(4)*(uind(2,i)*psc5+uinp(2,i)*dsc5)
     &           - rr7*sc(6)*(uind(2,i)*psc7+uinp(2,i)*dsc7))
     &           +(rr3*ci*(uind(2,j)*psc3+uinp(2,j)*dsc3)
     &           + rr5*sc(3)*(uind(2,j)*psc5+uinp(2,j)*dsc5)
     &           + rr7*sc(5)*(uind(2,j)*psc7+uinp(2,j)*dsc7))*0.5d0
     &           + rr5*scale5i*(sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &           + sci(3)*uinp(2,j)+scip(3)*uind(2,j))*0.5d0
     &           + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &           + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dj(2)
     &           + 0.5d0*gfi(4)*((qjui(2)-qiuj(2))*psc5
     &           + (qjuip(2)-qiujp(2))*dsc5)
     &           + gfi(5)*qir(2) + gfi(6)*qjr(2)
               ftm2i(3) = gfi(1)*zji  + 0.5d0*
     &          (- rr3*cj*(uind(3,i)*psc3+uinp(3,i)*dsc3)
     &           + rr5*sc(4)*(uind(3,i)*psc5+uinp(3,i)*dsc5)
     &           - rr7*sc(6)*(uind(3,i)*psc7+uinp(3,i)*dsc7))
     &           +(rr3*ci*(uind(3,j)*psc3+uinp(3,j)*dsc3)
     &           + rr5*sc(3)*(uind(3,j)*psc5+uinp(3,j)*dsc5)
     &           + rr7*sc(5)*(uind(3,j)*psc7+uinp(3,j)*dsc7))*0.5d0
     &           + rr5*scale5i*(sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &           + sci(3)*uinp(3,j)+scip(3)*uind(3,j))*0.5d0
     &           + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &           + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dj(3)
     &           + 0.5d0*gfi(4)*((qjui(3)-qiuj(3))*psc5
     &           + (qjuip(3)-qiujp(3))*dsc5)
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
                  gfd = (rr5*scip(2)*scale3i - rr7*(scip(3)*sci(4)
     &                 + sci(3)*scip(4))*scale5i)*0.5d0
                  fdir(1) = gfd*xji + 0.5d0*rr5*scale5i
     &                         * (sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &                           +sci(3)*uinp(1,j)+scip(3)*uind(1,j))
                  fdir(2) = gfd*yji + 0.5d0*rr5*scale5i
     &                         * (sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &                           +sci(3)*uinp(2,j)+scip(3)*uind(2,j))
                  fdir(3) = gfd*zji + 0.5d0*rr5*scale5i
     &                         * (sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &                           +sci(3)*uinp(3,j)+scip(3)*uind(3,j))
                  ftm2i(1) = ftm2i(1) - fdir(1) + findmp(1)
                  ftm2i(2) = ftm2i(2) - fdir(2) + findmp(2)
                  ftm2i(3) = ftm2i(3) - fdir(3) + findmp(3)
               end if
c
c     now perform the torque calculation
c     intermediate terms for torque between multipoles i and j
c
               gti(2) = 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5
               gti(3) = 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5
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
c     calculate the induced torque components
c
               ttm2i(1) = -rr3*(dixuj(1)*psc3+dixujp(1)*dsc3)*0.5d0
     &           + gti(2)*dixr(1) + gti(4)*((ujxqir(1)+rxqiuj(1))*psc5
     &           +(ujxqirp(1)+rxqiujp(1))*dsc5)*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -rr3*(dixuj(2)*psc3+dixujp(2)*dsc3)*0.5d0
     &           + gti(2)*dixr(2) + gti(4)*((ujxqir(2)+rxqiuj(2))*psc5
     &           +(ujxqirp(2)+rxqiujp(2))*dsc5)*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -rr3*(dixuj(3)*psc3+dixujp(3)*dsc3)*0.5d0
     &           + gti(2)*dixr(3) + gti(4)*((ujxqir(3)+rxqiuj(3))*psc5
     &           +(ujxqirp(3)+rxqiujp(3))*dsc5)*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -rr3*(djxui(1)*psc3+djxuip(1)*dsc3)*0.5d0
     &           + gti(3)*djxr(1) - gti(4)*((uixqjr(1)+rxqjui(1))*psc5
     &           +(uixqjrp(1)+rxqjuip(1))*dsc5)*0.5d0 - gti(6)*rxqjr(1)
               ttm3i(2) = -rr3*(djxui(2)*psc3+djxuip(2)*dsc3)*0.5d0
     &           + gti(3)*djxr(2) - gti(4)*((uixqjr(2)+rxqjui(2))*psc5
     &           +(uixqjrp(2)+rxqjuip(2))*dsc5)*0.5d0 - gti(6)*rxqjr(2)
               ttm3i(3) = -rr3*(djxui(3)*psc3+djxuip(3)*dsc3)*0.5d0
     &           + gti(3)*djxr(3) - gti(4)*((uixqjr(3)+rxqjui(3))*psc5
     &           +(uixqjrp(3)+rxqjuip(3))*dsc5)*0.5d0 - gti(6)*rxqjr(3)
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
c     now compute the virial components for this interaction
c
               xiz = x(zaxis(i)) - x(ii)
               yiz = y(zaxis(i)) - y(ii)
               ziz = z(zaxis(i)) - z(ii)
               xix = x(xaxis(i)) - x(ii)
               yix = y(xaxis(i)) - y(ii)
               zix = z(xaxis(i)) - z(ii)
               xjz = x(zaxis(j)) - x(jj)
               yjz = y(zaxis(j)) - y(jj)
               zjz = z(zaxis(j)) - z(jj)
               xjx = x(xaxis(j)) - x(jj)
               yjx = y(xaxis(j)) - y(jj)
               zjx = z(xaxis(j)) - z(jj)
               do k = 1, 3
                  tq(k) = ttm2(k) + ttm2i(k)
               end do
               call torque1 (i,tq,frczi,frcxi)
               do k = 1, 3
                  tq(k) = ttm3(k) + ttm3i(k)
               end do
               call torque1 (j,tq,frczj,frcxj)
               vxx = -xji*(ftm2(1)+ftm2i(1)) + xix*frcxi(1)
     &                  + xjz*frczj(1) + xjx*frcxj(1) + xiz*frczi(1)
               vyx = -yji*(ftm2(1)+ftm2i(1)) + yix*frcxi(1)
     &                  + yjz*frczj(1) + yjx*frcxj(1) + yiz*frczi(1)
               vzx = -zji*(ftm2(1)+ftm2i(1)) + zix*frcxi(1)
     &                  + zjz*frczj(1) + zjx*frcxj(1) + ziz*frczi(1)
               vyy = -yji*(ftm2(2)+ftm2i(2)) + yix*frcxi(2)
     &                  + yjz*frczj(2) + yjx*frcxj(2) + yiz*frczi(2)
               vzy = -zji*(ftm2(2)+ftm2i(2)) + zix*frcxi(2)
     &                  + zjz*frczj(2) + zjx*frcxj(2) + ziz*frczi(2)
               vzz = -zji*(ftm2(3)+ftm2i(3)) + zix*frcxi(3)
     &                  + zjz*frczj(3) + zjx*frcxj(3) + ziz*frczi(3)
               vir(1,1) = vir(1,1) + f*vxx
               vir(2,1) = vir(2,1) + f*vyx
               vir(3,1) = vir(3,1) + f*vzx
               vir(1,2) = vir(1,2) + f*vyx
               vir(2,2) = vir(2,2) + f*vyy
               vir(3,2) = vir(3,2) + f*vzy
               vir(1,3) = vir(1,3) + f*vzx
               vir(2,3) = vir(2,3) + f*vzy
               vir(3,3) = vir(3,3) + f*vzz
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
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine empole1b  --  Ewald summation multipole derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "empole1b" calculates the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c     using a regular Ewald summation
c
c
      subroutine empole1b
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'deriv.i'
      include 'energi.i'
      include 'ewald.i'
      include 'inter.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'units.i'
      include 'virial.i'
      integer i,j,ii
      real*8 e,ei,eintra
      real*8 f,term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 xup,yup,zup
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,xufield
      real*8 ydfield,yufield
      real*8 zdfield,zufield
      real*8 trqi(3,maxatm)
      real*8 trqd(3,maxatm)
c
c
c     zero out multipole and polarization energy and derivatives
c
      em = 0.0d0
      ep = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
            dep(j,i) = 0.0d0
         end do
      end do
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the induced dipole moment at each atom
c
      call induce
c
c     compute the reciprocal space part of the Ewald summation
c
      call erecip1
c
c     compute the real space part of the Ewald summation
c
      call ereal1 (eintra)
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         uii = dix*uix + diy*uiy + diz*uiz
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         ei = fterm * term * uii / 3.0d0
         em = em + e
         ep = ep + ei
      end do
c
c     compute the self-energy torque term due to induced dipole
c
      do i = 1, n
         trqi(1,i) = 0.0d0
         trqi(2,i) = 0.0d0
         trqi(3,i) = 0.0d0
         trqd(1,i) = 0.0d0
         trqd(2,i) = 0.0d0
         trqd(3,i) = 0.0d0
      end do
      term = (4.0d0/3.0d0) * f * aewald**3 / sqrtpi
      do i = 1, npole
         ii = ipole(i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         uix = 0.5d0 * (uind(1,i)+uinp(1,i))
         uiy = 0.5d0 * (uind(2,i)+uinp(2,i))
         uiz = 0.5d0 * (uind(3,i)+uinp(3,i))
         trqi(1,ii) = trqi(1,ii) + term*(diy*uiz-diz*uiy)
         trqi(2,ii) = trqi(2,ii) + term*(diz*uix-dix*uiz)
         trqi(3,ii) = trqi(3,ii) + term*(dix*uiy-diy*uix)
      end do
c
c     compute the cell dipole boundary correction term
c
      if (.not. tinfoil) then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xu = 0.0d0
         yu = 0.0d0
         zu = 0.0d0
         xup = 0.0d0
         yup = 0.0d0
         zup = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
            yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
            zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
            xu = xu + uind(1,i)
            yu = yu + uind(2,i)
            zu = zu + uind(3,i)
            xup = xup + uinp(1,i)
            yup = yup + uinp(2,i)
            zup = zup + uinp(3,i)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         em = em + term*(xd*xd+yd*yd+zd*zd)
         ep = ep + term*(xd*xu+yd*yu+zd*zu)
         do i = 1, npole
            ii = ipole(i)
            dem(1,ii) = dem(1,ii) + 2.0d0*term*rpole(1,i)*xd
            dem(2,ii) = dem(2,ii) + 2.0d0*term*rpole(1,i)*yd
            dem(3,ii) = dem(3,ii) + 2.0d0*term*rpole(1,i)*zd
            dep(1,ii) = dep(1,ii) + term*rpole(1,i)*(xu+xup)
            dep(2,ii) = dep(2,ii) + term*rpole(1,i)*(yu+yup)
            dep(3,ii) = dep(3,ii) + term*rpole(1,i)*(zu+zup)
         end do
         xdfield = -2.0d0 * term * xd
         ydfield = -2.0d0 * term * yd
         zdfield = -2.0d0 * term * zd
         xufield = -term * (xu+xup)
         yufield = -term * (yu+yup)
         zufield = -term * (zu+zup)
         do i = 1, npole
            ii = ipole(i)
            trqd(1,ii) = trqd(1,ii) + rpole(3,i)*zdfield
     &                      - rpole(4,i)*ydfield
            trqd(2,ii) = trqd(2,ii) + rpole(4,i)*xdfield
     &                      - rpole(2,i)*zdfield
            trqd(3,ii) = trqd(3,ii) + rpole(2,i)*ydfield
     &                      - rpole(3,i)*xdfield
            trqi(1,ii) = trqi(1,ii) + rpole(3,i)*zufield
     &                      - rpole(4,i)*yufield
            trqi(2,ii) = trqi(2,ii) + rpole(4,i)*xufield
     &                      - rpole(2,i)*zufield
            trqi(3,ii) = trqi(3,ii) + rpole(2,i)*yufield
     &                      - rpole(3,i)*xufield
         end do
c
c     boundary correction to virial due to overall cell dipole
c
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xq = 0.0d0
         yq = 0.0d0
         zq = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            xd = xd + rpole(2,i)
            yd = yd + rpole(3,i)
            zd = zd + rpole(4,i)
            xq = xq + rpole(1,i)*x(ii)
            yq = yq + rpole(1,i)*y(ii)
            zq = zq + rpole(1,i)*z(ii)
         end do
         xv = xq * (xd+0.5d0*(xu+xup))
         yv = yq * (yd+0.5d0*(yu+yup))
         zv = zq * (zd+0.5d0*(zu+zup))
         vterm = term * (xq*xq + yq*yq + zq*zq + 2.0d0*(xv+yv+zv)
     &                      + xu*xup + yu*yup + zu*zup
     &                      + xd*(xd+xu+xup) + yd*(yd+yu+yup)
     &                      + zd*(zd+zu+zup))
         vir(1,1) = vir(1,1) + 2.0d0*term*(xq*xq+xv) + vterm
         vir(2,1) = vir(2,1) + 2.0d0*term*(xq*yq+xv)
         vir(3,1) = vir(3,1) + 2.0d0*term*(xq*zq+xv)
         vir(1,2) = vir(1,2) + 2.0d0*term*(yq*xq+yv)
         vir(2,2) = vir(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
         vir(3,2) = vir(3,2) + 2.0d0*term*(yq*zq+yv)
         vir(1,3) = vir(1,3) + 2.0d0*term*(zq*xq+zv)
         vir(2,3) = vir(2,3) + 2.0d0*term*(zq*yq+zv)
         vir(3,3) = vir(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
         if (poltyp .eq. 'DIRECT') then
            vterm = term * (xu*xup+yu*yup+zu*zup)
            vir(1,1) = vir(1,1) + vterm
            vir(2,2) = vir(2,2) + vterm
            vir(3,3) = vir(3,3) + vterm
         end if
      end if
c
c     convert torque from self-interaction and boundary correction
c
      call torque (trqi,dep)
      call torque (trqd,dem)
c
c     intermolecular energy is total minus intramolecular part
c
      einter = einter + em + ep - eintra
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine erecip1  --  mpole Ewald recip energy & derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "erecip1" evaluates the reciprocal space portion of the regular
c     Ewald summation energy and gradient due to atomic multipole
c     interactions and dipole polarizability
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine erecip1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'deriv.i'
      include 'energi.i'
      include 'ewald.i'
      include 'ewreg.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'units.i'
      include 'virial.i'
      integer i,j,k,l
      integer ii,m1,m2
      integer jmin,jmax
      integer kmin,kmax
      integer lmin,lmax
      real*8 e,ei,etot,f,cut
      real*8 expterm,term,eterm
      real*8 uterm,vterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3,hsq
      real*8 ck,dk,qk,uk
      real*8 q1,q2,q3
      real*8 t1,t2,t3,t4
      real*8 ukp,t3p,t4p
      real*8 de,det1,det2
      real*8 dei,det1i,det2i
      real*8 wterm(3,3)
      real*8 t5(3,3),t6(3,3)
      real*8 t5u(3,3),t5p(3,3)
      real*8 t6u(3,3),t6p(3,3)
      real*8 qt(3,3),dt(3,3)
      real*8 dtu(3,3),dtp(3,3)
      real*8 ckr(maxatm),skr(maxatm)
      real*8 cjk(maxatm),sjk(maxatm)
      real*8 s1(maxatm),s2(maxatm)
      real*8 s3(maxatm),s4(maxatm)
      real*8 s3p(maxatm),s4p(maxatm)
      real*8 cm(maxatm),dm(3,maxatm)
      real*8 qm(9,maxatm),um(3,maxatm)
      real*8 trq(3,maxatm),trqi(3,maxatm)
      real*8 dkx(maxatm),qkx(maxatm)
      real*8 dky(maxatm),qky(maxatm)
      real*8 dkz(maxatm),qkz(maxatm)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      f = electric / dielec
      term = -0.25d0 / aewald**2
      eterm = 4.0d0 * pi * f / volbox
c
c     set the number of vectors based on box dimensions
c
      cut = 4.0d0 * pi * pi * frecip
      jmin = 0
      kmin = 0
      lmin = 1
      jmax = min(maxvec,int(frecip/recip(1,1)))
      kmax = min(maxvec,int(frecip/recip(2,2)))
      lmax = min(maxvec,int(frecip/recip(3,3)))
c
c     copy the multipole moments into local storage areas
c
      do i = 1, npole
         cm(i) = rpole(1,i)
         dm(1,i) = rpole(2,i)
         dm(2,i) = rpole(3,i)
         dm(3,i) = rpole(4,i)
         qm(1,i) = rpole(5,i)
         qm(2,i) = rpole(6,i)
         qm(3,i) = rpole(7,i)
         qm(4,i) = rpole(8,i)
         qm(5,i) = rpole(9,i)
         qm(6,i) = rpole(10,i)
         qm(7,i) = rpole(11,i)
         qm(8,i) = rpole(12,i)
         qm(9,i) = rpole(13,i)
         um(1,i) = uind(1,i)
         um(2,i) = uind(2,i)
         um(3,i) = uind(3,i)
      end do
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
c     calculate and store the exponential factors
c
      do i = 1, npole
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
         ejc(i,0) = 1.0d0
         ejs(i,0) = 0.0d0
         ekc(i,0) = 1.0d0
         eks(i,0) = 0.0d0
         elc(i,0) = 1.0d0
         els(i,0) = 0.0d0
         ejc(i,1) = cos(xfr)
         ejs(i,1) = sin(xfr)
         ekc(i,1) = cos(yfr)
         eks(i,1) = sin(yfr)
         elc(i,1) = cos(zfr)
         els(i,1) = sin(zfr)
         ekc(i,-1) = ekc(i,1)
         eks(i,-1) = -eks(i,1)
         elc(i,-1) = elc(i,1)
         els(i,-1) = -els(i,1)
         do j = 2, jmax
            ejc(i,j) = ejc(i,j-1)*ejc(i,1) - ejs(i,j-1)*ejs(i,1)
            ejs(i,j) = ejs(i,j-1)*ejc(i,1) + ejc(i,j-1)*ejs(i,1)
         end do
         do j = 2, kmax
            ekc(i,j) = ekc(i,j-1)*ekc(i,1) - eks(i,j-1)*eks(i,1)
            eks(i,j) = eks(i,j-1)*ekc(i,1) + ekc(i,j-1)*eks(i,1)
            ekc(i,-j) = ekc(i,j)
            eks(i,-j) = -eks(i,j)
         end do
         do j = 2, lmax
            elc(i,j) = elc(i,j-1)*elc(i,1) - els(i,j-1)*els(i,1)
            els(i,j) = els(i,j-1)*elc(i,1) + elc(i,j-1)*els(i,1)
            elc(i,-j) = elc(i,j)
            els(i,-j) = -els(i,j)
         end do
      end do
c
c     loop over all k vectors from the reciprocal lattice
c
      do j = jmin, jmax
         rj = 2.0d0 * pi * dble(j)
         do k = kmin, kmax
            rk = 2.0d0 * pi * dble(k)
            do i = 1, npole
               cjk(i) = ejc(i,j)*ekc(i,k) - ejs(i,j)*eks(i,k)
               sjk(i) = ejs(i,j)*ekc(i,k) + ejc(i,j)*eks(i,k)
            end do
            do l = lmin, lmax
               rl = 2.0d0 * pi * dble(l)
               h1 = recip(1,1)*rj
               h2 = recip(2,1)*rj + recip(2,2)*rk
               h3 = recip(3,1)*rj + recip(3,2)*rk + recip(3,3)*rl
               hsq = h1*h1 + h2*h2 + h3*h3
               if (hsq .le. cut) then
                  t1 = 0.0d0
                  t2 = 0.0d0
                  t3 = 0.0d0
                  t4 = 0.0d0
                  t3p = 0.0d0
                  t4p = 0.0d0
                  do m2 = 1, 3
                     do m1 = 1, 3
                        t5(m1,m2) = 0.0d0
                        t5u(m1,m2) = 0.0d0
                        t5p(m1,m2) = 0.0d0
                        t6(m1,m2) = 0.0d0
                        t6u(m1,m2) = 0.0d0
                        t6p(m1,m2) = 0.0d0
                     end do
                  end do
                  do i = 1, npole
                     ckr(i) = cjk(i)*elc(i,l) - sjk(i)*els(i,l)
                     skr(i) = sjk(i)*elc(i,l) + cjk(i)*els(i,l)
                     ck = cm(i)
                     dk = h1*dm(1,i) + h2*dm(2,i) + h3*dm(3,i)
                     dkx(i) = h3*dm(2,i) - h2*dm(3,i)
                     dky(i) = h1*dm(3,i) - h3*dm(1,i)
                     dkz(i) = h2*dm(1,i) - h1*dm(2,i)
                     q1 = h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i)
                     q2 = h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i)
                     q3 = h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i)
                     qk = h1*q1 + h2*q2 + h3*q3
                     qkx(i) = h3*q2 - h2*q3
                     qky(i) = h1*q3 - h3*q1
                     qkz(i) = h2*q1 - h1*q2
                     uk = h1*uind(1,i) + h2*uind(2,i) + h3*uind(3,i)
                     ukp = h1*uinp(1,i) + h2*uinp(2,i) + h3*uinp(3,i)
                     s1(i) = (ck-qk)*skr(i) + dk*ckr(i)
                     s2(i) = (ck-qk)*ckr(i) - dk*skr(i)
                     s3(i) = uk * ckr(i)
                     s4(i) = -uk * skr(i)
                     s3p(i) = ukp * ckr(i)
                     s4p(i) = -ukp * skr(i)
                     t1 = t1 + s1(i)
                     t2 = t2 + s2(i)
                     t3 = t3 + s3(i)
                     t4 = t4 + s4(i)
                     t3p = t3p + s3p(i)
                     t4p = t4p + s4p(i)
c
c     terms needed for subsequent virial tensor calculation
c
                     qt(1,1) = h1*(h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i))
                     qt(2,1) = h1*(h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i))
                     qt(3,1) = h1*(h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i))
                     qt(1,2) = h2*(h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i))
                     qt(2,2) = h2*(h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i))
                     qt(3,2) = h2*(h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i))
                     qt(1,3) = h3*(h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i))
                     qt(2,3) = h3*(h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i))
                     qt(3,3) = h3*(h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i))
                     dt(1,1) = h1 * dm(1,i)
                     dt(2,1) = h1 * dm(2,i)
                     dt(3,1) = h1 * dm(3,i)
                     dt(1,2) = h2 * dm(1,i)
                     dt(2,2) = h2 * dm(2,i)
                     dt(3,2) = h2 * dm(3,i)
                     dt(1,3) = h3 * dm(1,i)
                     dt(2,3) = h3 * dm(2,i)
                     dt(3,3) = h3 * dm(3,i)
                     dtu(1,1) = h1 * uind(1,i)
                     dtu(2,1) = h1 * uind(2,i)
                     dtu(3,1) = h1 * uind(3,i)
                     dtu(1,2) = h2 * uind(1,i)
                     dtu(2,2) = h2 * uind(2,i)
                     dtu(3,2) = h2 * uind(3,i)
                     dtu(1,3) = h3 * uind(1,i)
                     dtu(2,3) = h3 * uind(2,i)
                     dtu(3,3) = h3 * uind(3,i)
                     dtp(1,1) = h1 * uinp(1,i)
                     dtp(2,1) = h1 * uinp(2,i)
                     dtp(3,1) = h1 * uinp(3,i)
                     dtp(1,2) = h2 * uinp(1,i)
                     dtp(2,2) = h2 * uinp(2,i)
                     dtp(3,2) = h2 * uinp(3,i)
                     dtp(1,3) = h3 * uinp(1,i)
                     dtp(2,3) = h3 * uinp(2,i)
                     dtp(3,3) = h3 * uinp(3,i)
                     do m2 = 1, 3
                        do m1 = 1, 3
                           t5(m1,m2) = t5(m1,m2) - dt(m1,m2)*ckr(i)
     &                                    + 2.0d0*qt(m1,m2)*skr(i)
                           t5u(m1,m2) = t5u(m1,m2) - dtu(m1,m2)*ckr(i)
                           t5p(m1,m2) = t5p(m1,m2) - dtp(m1,m2)*ckr(i)
                           t6(m1,m2) = t6(m1,m2) + dt(m1,m2)*skr(i)
     &                                    + 2.0d0*qt(m1,m2)*ckr(i)
                           t6u(m1,m2) = t6u(m1,m2) + dtu(m1,m2)*skr(i)
                           t6p(m1,m2) = t6p(m1,m2) + dtp(m1,m2)*skr(i)
                        end do
                     end do
                  end do
c
c     get the energy contributions for current reciprocal vector
c
                  expterm = eterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0 
                  end if
                  e = expterm * (t1*t1+t2*t2)
                  ei = expterm * (t1*t3+t2*t4)
                  etot = e + ei
                  em = em + e
                  ep = ep + ei
c
c     get the virial contributions for current reciprocal vector
c
                  uterm = expterm * (t1*(t1+t3+t3p) + t3*t3p
     &                                 + t2*(t2+t4+t4p) + t4*t4p)
                  do m2 = 1, 3
                     do m1 = 1, 3
                        wterm(m1,m2) = 2.0d0 * expterm
     &                     * (t1*t5(m1,m2) + t2*t6(m1,m2)
     &                        + 0.5d0*(t1*(t5u(m1,m2)+t5p(m1,m2))
     &                        + t2*(t6u(m1,m2)+t6p(m1,m2))
     &                        + (t3+t3p)*t5(m1,m2)
     &                        + t3*t5p(m1,m2) + t3p*t5u(m1,m2)
     &                        + (t4+t4p)*t6(m1,m2)
     &                        + t4*t6p(m1,m2) + t4p*t6u(m1,m2)))
                     end do
                  end do
                  if (poltyp .eq. 'DIRECT') then
                     uterm = uterm - expterm*(t3*t3p+t4*t4p)
                     do m2 = 1, 3
                        do m1 = 1, 3
                           wterm(m1,m2) = wterm(m1,m2)
     &                        - expterm*(t3*t5p(m1,m2)+t3p*t5u(m1,m2)
     &                                  +t4*t6p(m1,m2)+t4p*t6u(m1,m2))
                        end do
                     end do
                  end if
                  wterm(2,1) = 0.5d0 * (wterm(2,1)+wterm(1,2))
                  wterm(3,1) = 0.5d0 * (wterm(3,1)+wterm(1,3))
                  wterm(3,2) = 0.5d0 * (wterm(3,2)+wterm(2,3))
                  wterm(1,2) = wterm(2,1)
                  wterm(1,3) = wterm(3,1)
                  wterm(2,3) = wterm(3,2)
                  vterm = 2.0d0 * uterm * (1.0d0-term*hsq) / hsq
                  vir(1,1) = vir(1,1) + h1*h1*vterm + wterm(1,1) - uterm
                  vir(2,1) = vir(2,1) + h2*h1*vterm + wterm(2,1)
                  vir(3,1) = vir(3,1) + h3*h1*vterm + wterm(3,1)
                  vir(1,2) = vir(1,2) + h1*h2*vterm + wterm(1,2)
                  vir(2,2) = vir(2,2) + h2*h2*vterm + wterm(2,2) - uterm
                  vir(3,2) = vir(3,2) + h3*h2*vterm + wterm(3,2)
                  vir(1,3) = vir(1,3) + h1*h3*vterm + wterm(1,3)
                  vir(2,3) = vir(2,3) + h2*h3*vterm + wterm(2,3)
                  vir(3,3) = vir(3,3) + h3*h3*vterm + wterm(3,3) - uterm
c
c     get the force contributions for current reciprocal vector
c
                  expterm = 2.0d0 * expterm
                  do i = 1, npole
                     ii = ipole(i)
                     de = expterm * (s2(i)*t1-s1(i)*t2)
                     dei = 0.5d0 * expterm * ((s4(i)+s4p(i))*t1
     &                                       -(s3(i)+s3p(i))*t2
     &                                +s2(i)*(t3+t3p)-s1(i)*(t4+t4p))
                     if (poltyp .eq. 'MUTUAL') then
                         dei = dei + 0.5d0 * expterm
     &                            * (s4p(i)*t3+s4(i)*t3p
     &                              -s3p(i)*t4-s3(i)*t4p)
                     end if
                     det1 = expterm * (skr(i)*t2-ckr(i)*t1)
                     det2 = 2.0d0 * expterm * (ckr(i)*t2+skr(i)*t1)
                     det1i = 0.5d0 * expterm * (skr(i)*(t4+t4p)
     &                                         -ckr(i)*(t3+t3p))
                     det2i = expterm * (ckr(i)*(t4+t4p)+skr(i)*(t3+t3p))
                     dem(1,ii) = dem(1,ii) + h1*de
                     dem(2,ii) = dem(2,ii) + h2*de
                     dem(3,ii) = dem(3,ii) + h3*de
                     dep(1,ii) = dep(1,ii) + h1*dei
                     dep(2,ii) = dep(2,ii) + h2*dei
                     dep(3,ii) = dep(3,ii) + h3*dei
                     trq(1,ii) = trq(1,ii) + dkx(i)*det1 + qkx(i)*det2
                     trq(2,ii) = trq(2,ii) + dky(i)*det1 + qky(i)*det2
                     trq(3,ii) = trq(3,ii) + dkz(i)*det1 + qkz(i)*det2
                     trqi(1,ii) = trqi(1,ii) + dkx(i)*det1i
     &                               + qkx(i)*det2i
                     trqi(2,ii) = trqi(2,ii) + dky(i)*det1i
     &                               + qky(i)*det2i
                     trqi(3,ii) = trqi(3,ii) + dkz(i)*det1i
     &                               + qkz(i)*det2i
                  end do
               end if
            end do
            lmin = -lmax
         end do
         kmin = -kmax
      end do
c
c     convert the torques to forces and increment the totals
c
      call torque (trq,dem)
      call torque (trqi,dep)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ereal1  --  mpole ewald real energy and derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ereal1" evaluates the real space portion of the regular Ewald
c     summation energy and gradient due to atomic multipole interactions
c     and dipole polarizability
c
c
      subroutine ereal1 (eintra)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'cell.i'
      include 'couple.i'
      include 'cutoff.i'
      include 'deriv.i'
      include 'energi.i'
      include 'ewald.i'
      include 'math.i'
      include 'molcul.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'units.i'
      include 'virial.i'
      integer i,j,k,m
      integer ii,jj,kn
      real*8 e,ei,f,bfac
      real*8 eintra,erfc
      real*8 pterm,damp
      real*8 scale3,scale5
      real*8 scale7,scale9
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 ruscale3,ruscale5
      real*8 ruscale7
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 gfdr,gfd
      real*8 xji,yji,zji
      real*8 xix,yix,zix
      real*8 xiz,yiz,ziz
      real*8 xjx,yjx,zjx
      real*8 xjz,yjz,zjz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 einin,erl,erli
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 vxy,vxz,vyz
      real*8 frcxi(3),frczi(3)
      real*8 frcxj(3),frczj(3)
      real*8 ci,di(3),qi(9)
      real*8 cj,dj(3),qj(9)
      real*8 fridmp(3),findmp(3)
      real*8 ftm2(3),ftm2i(3)
      real*8 ftm2r(3),ftm2ri(3)
      real*8 ttm2(3),ttm3(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 ttm2r(3),ttm3r(3)
      real*8 ttm2ri(3),ttm3ri(3)
      real*8 fdir(3),dixdj(3)
      real*8 djxui(3),dixuj(3)
      real*8 dixujp(3),djxuip(3)
      real*8 uixqjr(3),ujxqir(3)
      real*8 uixqjrp(3),ujxqirp(3)
      real*8 qiuj(3),qjui(3)
      real*8 qiujp(3),qjuip(3)
      real*8 rxqiuj(3),rxqjui(3)
      real*8 rxqiujp(3),rxqjuip(3)
      real*8 qidj(3),qjdi(3)
      real*8 qir(3),qjr(3)
      real*8 qiqjr(3),qjqir(3)
      real*8 qixqj(3),rxqir(3)
      real*8 dixr(3),djxr(3)
      real*8 dixqjr(3),djxqir(3)
      real*8 rxqjr(3),qjrxqir(3)
      real*8 rxqijr(3),rxqjir(3)
      real*8 rxqidj(3),rxqjdi(3)
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3),tq(3)
      real*8 bn(0:5)
      real*8 sc(10),gl(0:8)
      real*8 sci(8),scip(8)
      real*8 gli(7),glip(7)
      real*8 gf(7),gfi(6)
      real*8 gfr(7),gfri(6)
      real*8 gti(6),gtri(6)
      real*8 mscale(maxatm)
      real*8 pscale(maxatm)
      real*8 dscale(maxatm)
      real*8 uscale(maxatm)
      real*8 ftm1(3,maxatm)
      real*8 ftm1i(3,maxatm)
      real*8 ttm1(3,maxatm)
      real*8 ttm1i(3,maxatm)
      real*8 trq(3,maxatm)
      real*8 trqi(3,maxatm)
      external erfc
c
c
c     zero out the intramolecular portion of the Ewald energy
c
      eintra = 0.0d0
c
c     set conversion factor, cutoff and scaling coefficients
c
      if (npole .eq. 0)  return
      f = electric / dielec
      call switch ('EWALD')
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
      do i = 1, npole-1
         ii = ipole(i)
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
         do j = 1, n
            mscale(j) = 1.0d0
            dscale(j) = 1.0d0
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
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
         end do
         do j = i+1, npole
            jj = ipole(j)
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
c
c     calculate the real space error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do m = 1, 5
                  bfac = dble(m+m-1)
                  alsq2n = alsq2 * alsq2n
                  bn(m) = (bfac*bn(m-1)+alsq2n*exp2a) / r2
               end do
c
c     zero out the temporary force and torque
c
               do k = 1, 3
                  ftm2(k) = 0.0d0
                  ftm2i(k) = 0.0d0
                  ttm2(k) = 0.0d0
                  ttm2i(k) = 0.0d0
                  ttm3(k) = 0.0d0
                  ttm3i(k) = 0.0d0
               end do
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               scale9 = 1.0d0
               do k = 1, 3
                  ddsc3(k) = 0.0d0
                  ddsc5(k) = 0.0d0
                  ddsc7(k) = 0.0d0
               end do
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
                     ddsc3(1) = -3.0d0*damp*exp(damp) * xji/r2
                     ddsc3(2) = -3.0d0*damp*exp(damp) * yji/r2
                     ddsc3(3) = -3.0d0*damp*exp(damp) * zji/r2
                     ddsc5(1) = -damp * ddsc3(1)
                     ddsc5(2) = -damp * ddsc3(2)
                     ddsc5(3) = -damp * ddsc3(3)
                     ddsc7(1) = (-0.2d0-0.6d0*damp) * ddsc5(1)
                     ddsc7(2) = (-0.2d0-0.6d0*damp) * ddsc5(2)
                     ddsc7(3) = (-0.2d0-0.6d0*damp) * ddsc5(3)
                  end if
               end if
               dsc3 = 1.0d0 - scale3*dscale(jj)
               dsc5 = 1.0d0 - scale5*dscale(jj)
               dsc7 = 1.0d0 - scale7*dscale(jj)
               psc3 = 1.0d0 - scale3*pscale(jj)
               psc5 = 1.0d0 - scale5*pscale(jj)
               psc7 = 1.0d0 - scale7*pscale(jj)
               ruscale3 = 1.0d0 - scale3*uscale(jj)
               ruscale5 = 1.0d0 - scale5*uscale(jj)
               ruscale7 = 1.0d0 - scale7*uscale(jj)
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
               dixujp(1) = di(2)*uinp(3,j) - di(3)*uinp(2,j)
               dixujp(2) = di(3)*uinp(1,j) - di(1)*uinp(3,j)
               dixujp(3) = di(1)*uinp(2,j) - di(2)*uinp(1,j)
               djxuip(1) = dj(2)*uinp(3,i) - dj(3)*uinp(2,i)
               djxuip(2) = dj(3)*uinp(1,i) - dj(1)*uinp(3,i)
               djxuip(3) = dj(1)*uinp(2,i) - dj(2)*uinp(1,i)
               dixr(1) = di(2)*zji - di(3)*yji
               dixr(2) = di(3)*xji - di(1)*zji
               dixr(3) = di(1)*yji - di(2)*xji
               djxr(1) = dj(2)*zji - dj(3)*yji
               djxr(2) = dj(3)*xji - dj(1)*zji
               djxr(3) = dj(1)*yji - dj(2)*xji
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
               qiujp(1)= qi(1)*uinp(1,j)+qi(4)*uinp(2,j)+qi(7)*uinp(3,j)
               qiujp(2)= qi(2)*uinp(1,j)+qi(5)*uinp(2,j)+qi(8)*uinp(3,j)
               qiujp(3)= qi(3)*uinp(1,j)+qi(6)*uinp(2,j)+qi(9)*uinp(3,j)
               qjuip(1)= qj(1)*uinp(1,i)+qj(4)*uinp(2,i)+qj(7)*uinp(3,i)
               qjuip(2)= qj(2)*uinp(1,i)+qj(5)*uinp(2,i)+qj(8)*uinp(3,i)
               qjuip(3)= qj(3)*uinp(1,i)+qj(6)*uinp(2,i)+qj(9)*uinp(3,i)
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
               uixqjrp(1) = uinp(2,i)*qjr(3) - uinp(3,i)*qjr(2)
               uixqjrp(2) = uinp(3,i)*qjr(1) - uinp(1,i)*qjr(3)
               uixqjrp(3) = uinp(1,i)*qjr(2) - uinp(2,i)*qjr(1)
               ujxqirp(1) = uinp(2,j)*qir(3) - uinp(3,j)*qir(2)
               ujxqirp(2) = uinp(3,j)*qir(1) - uinp(1,j)*qir(3)
               ujxqirp(3) = uinp(1,j)*qir(2) - uinp(2,j)*qir(1)
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
               rxqiujp(1) = yji*qiujp(3) - zji*qiujp(2)
               rxqiujp(2) = zji*qiujp(1) - xji*qiujp(3)
               rxqiujp(3) = xji*qiujp(2) - yji*qiujp(1)
               rxqjuip(1) = yji*qjuip(3) - zji*qjuip(2)
               rxqjuip(2) = zji*qjuip(1) - xji*qjuip(3)
               rxqjuip(3) = xji*qjuip(2) - yji*qjuip(1)
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
c     get intermediate variables for induction energy;
c     interaction between induce and permanent multipoles
c     final induction energy exclude induce-induce term
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
               scip(1) = uinp(1,i)*dj(1) + uinp(2,i)*dj(2)
     &                      + uinp(3,i)*dj(3) + di(1)*uinp(1,j)
     &                      + di(2)*uinp(2,j) + di(3)*uinp(3,j)
               scip(2) = uind(1,i)*uinp(1,j)+uind(2,i)*uinp(2,j)
     &                      + uind(3,i)*uinp(3,j)+uinp(1,i)*uind(1,j)
     &                      + uinp(2,i)*uind(2,j)+uinp(3,i)*uind(3,j)
               scip(3) = uinp(1,i)*xji + uinp(2,i)*yji + uinp(3,i)*zji
               scip(4) = uinp(1,j)*xji + uinp(2,j)*yji + uinp(3,j)*zji
               scip(7) = qir(1)*uinp(1,j) + qir(2)*uinp(2,j)
     &                      + qir(3)*uinp(3,j)
               scip(8) = qjr(1)*uinp(1,i) + qjr(2)*uinp(2,i)
     &                      + qjr(3)*uinp(3,i)
c
c     get the induced-induced term
c
               einin = sci(2)*rr3 - sci(3)*sci(4)*rr5
               findmp(1) = uscale(jj) * (scip(2)*rr3*ddsc3(1)
     &          -(sci(3)*scip(4)+scip(3)*sci(4))*rr5*ddsc5(1))*0.5d0
               findmp(2) = uscale(jj) * (scip(2)*rr3*ddsc3(2)
     &          -(sci(3)*scip(4)+scip(3)*sci(4))*rr5*ddsc5(2))*0.5d0
               findmp(3) = uscale(jj) * (scip(2)*rr3*ddsc3(3)
     &          -(sci(3)*scip(4)+scip(3)*sci(4))*rr5*ddsc5(3))*0.5d0
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
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = cj*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     calculate energy for this interaction and increment total
c
               e = bn(0)*gl(0) + bn(1)*(gl(1)+gl(6))
     &                + bn(2)*(gl(2)+gl(7)+gl(8))
     &                + bn(3)*(gl(3)+gl(5)) + bn(4)*gl(4)
               ei = 0.5d0 * (bn(1)*(gli(1)+gli(6))
     &                      + bn(2)*(gli(2)+gli(7)) + bn(3)*gli(3))
c
c    real energy without screening function
c
               erl = rr1*gl(0) + rr3*(gl(1)+gl(6))
     &                  + rr5*(gl(2)+gl(7)+gl(8))
     &                  + rr7*(gl(3)+gl(5)) + rr9*gl(4)
               erli = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &                   + rr5*(gli(2)+gli(7))*psc5
     &                   + rr7*gli(3)*psc7)
               e = e - (1.0d0-mscale(jj))*erl
               ei = ei - erli
               e = f * e
               ei = f * ei
               em = em + e
               ep = ep + ei
c
c
c
               fridmp(1) = 0.5d0*(rr3*((gli(1)+gli(6))*pscale(jj)
     &           + (glip(1)+glip(6))*dscale(jj))*ddsc3(1)
     &           + rr5*((gli(2)+gli(7))*pscale(jj)
     &           + (glip(2)+glip(7))*dscale(jj))*ddsc5(1)
     &           + rr7*(gli(3)*pscale(jj)+glip(3)*dscale(jj))*ddsc7(1))
               fridmp(2) = 0.5d0*(rr3*((gli(1)+gli(6))*pscale(jj)
     &           + (glip(1)+glip(6))*dscale(jj))*ddsc3(2)
     &           + rr5*((gli(2)+gli(7))*pscale(jj)
     &           + (glip(2)+glip(7))*dscale(jj))*ddsc5(2)
     &           + rr7*(gli(3)*pscale(jj)+glip(3)*dscale(jj))*ddsc7(2))
               fridmp(3) = 0.5d0*(rr3*((gli(1)+gli(6))*pscale(jj)
     &           + (glip(1)+glip(6))*dscale(jj))*ddsc3(3)
     &           + rr5*((gli(2)+gli(7))*pscale(jj)
     &           + (glip(2)+glip(7))*dscale(jj))*ddsc5(3)
     &           + rr7*(gli(3)*pscale(jj)+glip(3)*dscale(jj))*ddsc7(3))
c
c     increment the intramolecular energy if appropriate;
c     assumes intramolecular distances are less than half
c     of cell length and less than the ewald cutoff
c
               if (molcule(ii) .eq. molcule(jj)) then
                  eintra = eintra + mscale(jj)*erl*f
                  eintra = eintra + 0.5d0*pscale(jj)
     &                        * (rr3*(gli(1)+gli(6))*scale3
     &                              + rr5*(gli(2)+gli(7))*scale5
     &                              + rr7*gli(3)*scale7)
               end if
c
c     intermediate variables for permanent multipole force terms
c
               gf(1) = bn(1)*gl(0) + bn(2)*(gl(1)+gl(6))
     &                    + bn(3)*(gl(2)+gl(7)+gl(8))
     &                    + bn(4)*(gl(3)+gl(5)) + bn(5)*gl(4)
               gf(2) = -cj*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gf(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gf(4) = 2.0d0 * bn(2)
               gf(5) = 2.0d0 * (-cj*bn(2)+sc(4)*bn(3)-sc(6)*bn(4))
               gf(6) = 2.0d0 * (-ci*bn(2)-sc(3)*bn(3)-sc(5)*bn(4))
               gf(7) = 4.0d0 * bn(3)
               gfr(1) = rr3*gl(0) + rr5*(gl(1)+gl(6))
     &                     + rr7*(gl(2)+gl(7)+gl(8))
     &                     + rr9*(gl(3)+gl(5)) + rr11*gl(4)
               gfr(2) = -cj*rr3 + sc(4)*rr5 - sc(6)*rr7
               gfr(3) = ci*rr3 + sc(3)*rr5 + sc(5)*rr7
               gfr(4) = 2.0d0 * rr5
               gfr(5) = 2.0d0 * (-cj*rr5+sc(4)*rr7-sc(6)*rr9)
               gfr(6) = 2.0d0 * (-ci*rr5-sc(3)*rr7-sc(5)*rr9)
               gfr(7) = 4.0d0 * rr7
c
c     get the induced-permanent term
c
               gfi(1) = bn(2)*(gli(1)+glip(1)+gli(6)+glip(6))*0.5d0
     &                + bn(2)* scip(2)*0.5d0
     &                + bn(3)*(gli(2)+glip(2)+gli(7)+glip(7))*0.5d0
     &                - bn(3)*0.5d0*(sci(3)*scip(4)+scip(3)*sci(4))
     &                + bn(4)*(gli(3)+glip(3))*0.5d0
               gfi(2) = -cj*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gfi(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gfi(4) = 2.0d0 * bn(2)
               gfi(5) = (sci(4) + scip(4))*bn(3)
               gfi(6) = -(sci(3) + scip(3))*bn(3)
               gfri(1) = rr5*0.5d0*((gli(1)+gli(6))*psc3
     &             + (glip(1)+glip(6))*dsc3 + scip(2)*ruscale3)
     &             + rr7*((gli(7)+gli(2))*psc5 + (glip(7)+glip(2))*dsc5
     &             - (sci(3)*scip(4)+scip(3)*sci(4))*ruscale5)*0.5d0
     &             + rr9*(gli(3)*psc7+glip(3)*dsc7)*0.5d0
               gfri(2) = -rr3*cj + rr5*sc(4) - rr7*sc(6)
               gfri(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfri(4) = 2.0d0 * rr5
               gfri(5) = rr7 * (sci(4)*psc7 + scip(4)*dsc7)
               gfri(6) = -rr7 * (sci(3)*psc7 + scip(3)*dsc7)
c
c       get the induced-induced terms for direct versus mutual
c
               gfdr = 0.5d0 * (rr5*scip(2)*ruscale3
     &                  - rr7*(scip(3)*sci(4)
     &                        +sci(3)*scip(4))*ruscale5)
               gfd = 0.5d0 * (bn(2)*scip(2)
     &                  - bn(3)*(scip(3)*sci(4)+sci(3)*scip(4)))
c
c     get the unscreened induced force components
c
               ftm2ri(1) = gfri(1)*xji + 0.5d0*
     &          (- rr3*cj*(uind(1,i)*psc3+uinp(1,i)*dsc3)
     &           + rr5*sc(4)*(uind(1,i)*psc5+uinp(1,i)*dsc5)
     &           - rr7*sc(6)*(uind(1,i)*psc7+uinp(1,i)*dsc7))
     &           +(rr3*ci*(uind(1,j)*psc3+uinp(1,j)*dsc3)
     &           + rr5*sc(3)*(uind(1,j)*psc5+uinp(1,j)*dsc5)
     &           + rr7*sc(5)*(uind(1,j)*psc7+uinp(1,j)*dsc7))*0.5d0
     &           + rr5*ruscale5*(sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &           + sci(3)*uinp(1,j)+scip(3)*uind(1,j))*0.5d0
     &           + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &           + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dj(1)
     &           + 0.5d0*gfri(4)*((qjui(1)-qiuj(1))*psc5
     &           + (qjuip(1)-qiujp(1))*dsc5)
     &           + gfri(5)*qir(1) + gfri(6)*qjr(1)
               ftm2ri(2) = gfri(1)*yji + 0.5d0*
     &          (- rr3*cj*(uind(2,i)*psc3+uinp(2,i)*dsc3)
     &           + rr5*sc(4)*(uind(2,i)*psc5+uinp(2,i)*dsc5)
     &           - rr7*sc(6)*(uind(2,i)*psc7+uinp(2,i)*dsc7))
     &           +(rr3*ci*(uind(2,j)*psc3+uinp(2,j)*dsc3)
     &           + rr5*sc(3)*(uind(2,j)*psc5+uinp(2,j)*dsc5)
     &           + rr7*sc(5)*(uind(2,j)*psc7+uinp(2,j)*dsc7))*0.5d0
     &           + rr5*ruscale5*(sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &           + sci(3)*uinp(2,j)+scip(3)*uind(2,j))*0.5d0
     &           + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &           + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dj(2)
     &           + 0.5d0*gfri(4)*((qjui(2)-qiuj(2))*psc5
     &           + (qjuip(2)-qiujp(2))*dsc5)
     &           + gfri(5)*qir(2) + gfri(6)*qjr(2)
               ftm2ri(3) = gfri(1)*zji  + 0.5d0*
     &          (- rr3*cj*(uind(3,i)*psc3+uinp(3,i)*dsc3)
     &           + rr5*sc(4)*(uind(3,i)*psc5+uinp(3,i)*dsc5)
     &           - rr7*sc(6)*(uind(3,i)*psc7+uinp(3,i)*dsc7))
     &           +(rr3*ci*(uind(3,j)*psc3+uinp(3,j)*dsc3)
     &           + rr5*sc(3)*(uind(3,j)*psc5+uinp(3,j)*dsc5)
     &           + rr7*sc(5)*(uind(3,j)*psc7+uinp(3,j)*dsc7))*0.5d0
     &           + rr5*ruscale5*(sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &           + sci(3)*uinp(3,j)+scip(3)*uind(3,j))*0.5d0
     &           + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &           + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dj(3)
     &           + 0.5d0*gfri(4)*((qjui(3)-qiuj(3))*psc5
     &           + (qjuip(3)-qiujp(3))*dsc5)
     &           + gfri(5)*qir(3) + gfri(6)*qjr(3)
c
c     get the real-space induced force
c
               ftm2i(1)=gfi(1)*xji + 0.5d0*(gfi(2)*(uind(1,i)+uinp(1,i))
     &              + bn(2)*(sci(4)*uinp(1,i)+scip(4)*uind(1,i))
     &              + gfi(3)*(uind(1,j)+uinp(1,j))
     &              + bn(2)*(sci(3)*uinp(1,j)+scip(3)*uind(1,j))
     &              + (sci(4)+scip(4))*bn(2)*di(1)
     &              + (sci(3)+scip(3))*bn(2)*dj(1)
     &              + gfi(4)*(qjui(1)+qjuip(1)-qiuj(1)-qiujp(1)))
     &              + gfi(5)*qir(1) + gfi(6)*qjr(1)
               ftm2i(2)=gfi(1)*yji + 0.5d0*(gfi(2)*(uind(2,i)+uinp(2,i))
     &              + bn(2)*(sci(4)*uinp(2,i)+scip(4)*uind(2,i))
     &              + gfi(3)*(uind(2,j)+uinp(2,j))
     &              + bn(2)*(sci(3)*uinp(2,j)+scip(3)*uind(2,j))
     &              + (sci(4)+scip(4))*bn(2)*di(2)
     &              + (sci(3)+scip(3))*bn(2)*dj(2)
     &              + gfi(4)*(qjui(2)+qjuip(2)-qiuj(2)-qiujp(2)))
     &              + gfi(5)*qir(2) + gfi(6)*qjr(2)
               ftm2i(3)=gfi(1)*zji + 0.5d0*(gfi(2)*(uind(3,i)+uinp(3,i))
     &              + bn(2)*(sci(4)*uinp(3,i)+scip(4)*uind(3,i))
     &              + gfi(3)*(uind(3,j)+uinp(3,j))
     &              + bn(2)*(sci(3)*uinp(3,j)+scip(3)*uind(3,j))
     &              + (sci(4)+scip(4))*bn(2)*di(3)
     &              + (sci(3)+scip(3))*bn(2)*dj(3)
     &              + gfi(4)*(qjui(3)+qjuip(3)-qiuj(3)-qiujp(3)))
     &              + gfi(5)*qir(3) + gfi(6)*qjr(3)
c
c     get the real permanent force without screening function
c
               ftm2r(1) = gfr(1)*xji + gfr(2)*di(1) + gfr(3)*dj(1)
     &                     + gfr(4)*(qjdi(1)-qidj(1)) + gfr(5)*qir(1)
     &                     + gfr(6)*qjr(1) + gfr(7)*(qiqjr(1)+qjqir(1))
               ftm2r(2) = gfr(1)*yji + gfr(2)*di(2) + gfr(3)*dj(2)
     &                     + gfr(4)*(qjdi(2)-qidj(2)) + gfr(5)*qir(2)
     &                     + gfr(6)*qjr(2) + gfr(7)*(qiqjr(2)+qjqir(2))
               ftm2r(3) = gfr(1)*zji + gfr(2)*di(3) + gfr(3)*dj(3)
     &                     + gfr(4)*(qjdi(3)-qidj(3)) + gfr(5)*qir(3)
     &                     + gfr(6)*qjr(3) + gfr(7)*(qiqjr(3)+qjqir(3))
c
c     get the real-space permanent force
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
               if (poltyp .eq. 'DIRECT') then
                  ftm2i(1) = ftm2i(1) - gfd*xji - 0.5d0*
     &               (bn(2)*(sci(4)*uinp(1,i)+scip(4)*uind(1,i))
     &               + bn(2)*(sci(3)*uinp(1,j)+scip(3)*uind(1,j)))
                  ftm2i(2) = ftm2i(2) - gfd*yji -0.5d0*
     &               (bn(2)*(sci(4)*uinp(2,i)+scip(4)*uind(2,i))
     &               + bn(2)*(sci(3)*uinp(2,j)+scip(3)*uind(2,j)))
                  ftm2i(3) = ftm2i(3) - gfd*zji -0.5d0*
     &               (bn(2)*(sci(4)*uinp(3,i)+scip(4)*uind(3,i))
     &               + bn(2)*(sci(3)*uinp(3,j)+scip(3)*uind(3,j)))
                  fdir(1) = gfdr*xji + 0.5d0*ruscale5*rr5*
     &               (sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &               + sci(3)*uinp(1,j)+scip(3)*uind(1,j))
                  fdir(2) = gfdr*yji + 0.5d0*ruscale5*rr5*
     &               (sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &               + sci(3)*uinp(2,j)+scip(3)*uind(2,j))
                  fdir(3) = gfdr*zji + 0.5d0*ruscale5*rr5*
     &               (sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &               + sci(3)*uinp(3,j)+scip(3)*uind(3,j))
               end if
c
c     handle of scaling for partially excluded interactions
c
               ftm2(1) = ftm2(1) - (1.0d0-mscale(jj))*ftm2r(1)
               ftm2(2) = ftm2(2) - (1.0d0-mscale(jj))*ftm2r(2)
               ftm2(3) = ftm2(3) - (1.0d0-mscale(jj))*ftm2r(3)
               ftm2i(1) = ftm2i(1) - ftm2ri(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - ftm2ri(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - ftm2ri(3) - fridmp(3) - findmp(3)
               if (poltyp .eq. 'DIRECT') then
                  ftm2i(1) = ftm2i(1) + fdir(1) + findmp(1)
                  ftm2i(2) = ftm2i(2) + fdir(2) + findmp(2)
                  ftm2i(3) = ftm2i(3) + fdir(3) + findmp(3)
               end if
c
c     now perform the torque calculation
c     intermediate terms for torque between multipoles i and j
c
               gtri(2) = 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5
               gtri(3) = 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5
               gtri(4) = gfri(4)
               gtri(5) = gfri(5)
               gtri(6) = gfri(6)
               gti(2) = 0.5d0*(sci(4)+scip(4))*bn(2)
               gti(3) = 0.5d0*(sci(3)+scip(3))*bn(2)
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
c
c     get the unscreened induced torque
c
               ttm2ri(1) = -rr3*(dixuj(1)*psc3+dixujp(1)*dsc3)*0.5d0
     &           + gtri(2)*dixr(1) + gtri(4)*((ujxqir(1)+rxqiuj(1))*psc5
     &           +(ujxqirp(1)+rxqiujp(1))*dsc5)*0.5d0 - gtri(5)*rxqir(1)
               ttm2ri(2) = -rr3*(dixuj(2)*psc3+dixujp(2)*dsc3)*0.5d0
     &           + gtri(2)*dixr(2) + gtri(4)*((ujxqir(2)+rxqiuj(2))*psc5
     &           +(ujxqirp(2)+rxqiujp(2))*dsc5)*0.5d0 - gtri(5)*rxqir(2)
               ttm2ri(3) = -rr3*(dixuj(3)*psc3+dixujp(3)*dsc3)*0.5d0
     &           + gtri(2)*dixr(3) + gtri(4)*((ujxqir(3)+rxqiuj(3))*psc5
     &           +(ujxqirp(3)+rxqiujp(3))*dsc5)*0.5d0 - gtri(5)*rxqir(3)
               ttm3ri(1) = -rr3*(djxui(1)*psc3+djxuip(1)*dsc3)*0.5d0
     &           + gtri(3)*djxr(1) - gtri(4)*((uixqjr(1)+rxqjui(1))*psc5
     &           +(uixqjrp(1)+rxqjuip(1))*dsc5)*0.5d0 - gtri(6)*rxqjr(1)
               ttm3ri(2) = -rr3*(djxui(2)*psc3+djxuip(2)*dsc3)*0.5d0
     &           + gtri(3)*djxr(2) - gtri(4)*((uixqjr(2)+rxqjui(2))*psc5
     &           +(uixqjrp(2)+rxqjuip(2))*dsc5)*0.5d0 - gtri(6)*rxqjr(2)
               ttm3ri(3) = -rr3*(djxui(3)*psc3+djxuip(3)*dsc3)*0.5d0
     &           + gtri(3)*djxr(3) - gtri(4)*((uixqjr(3)+rxqjui(3))*psc5
     &           +(uixqjrp(3)+rxqjuip(3))*dsc5)*0.5d0 - gtri(6)*rxqjr(3)
c
c     get the real-space induced torque
c
               ttm2i(1) = -bn(1)*(dixuj(1)+dixujp(1))*0.5d0
     &              + gti(2)*dixr(1) + gti(4)*(ujxqir(1)+rxqiuj(1)
     &              + ujxqirp(1)+rxqiujp(1))*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -bn(1)*(dixuj(2)+dixujp(2))*0.5d0
     &              + gti(2)*dixr(2) + gti(4)*(ujxqir(2)+rxqiuj(2)
     &              + ujxqirp(2)+rxqiujp(2))*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -bn(1)*(dixuj(3)+dixujp(3))*0.5d0
     &              + gti(2)*dixr(3) + gti(4)*(ujxqir(3)+rxqiuj(3)
     &              + ujxqirp(3)+rxqiujp(3))*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -bn(1)*(djxui(1)+djxuip(1))*0.5d0
     &              + gti(3)*djxr(1) - gti(4)*(uixqjr(1)+rxqjui(1)
     &              + uixqjrp(1)+rxqjuip(1))*0.5d0 - gti(6)*rxqjr(1)
               ttm3i(2) = -bn(1)*(djxui(2)+djxuip(2))*0.5d0
     &              + gti(3)*djxr(2) - gti(4)*(uixqjr(2)+rxqjui(2)
     &              + uixqjrp(2)+rxqjuip(2))*0.5d0 - gti(6)*rxqjr(2)
               ttm3i(3) = -bn(1)*(djxui(3)+djxuip(3))*0.5d0
     &              + gti(3)*djxr(3) - gti(4)*(uixqjr(3)+rxqjui(3)
     &              + uixqjrp(3)+rxqjuip(3))*0.5d0 - gti(6)*rxqjr(3)
c
c     get the real-space permanent torques
c
               ttm2(1) = -bn(1)*dixdj(1) + gf(2)*dixr(1)
     &            + gf(4)*(dixqjr(1)+djxqir(1)+rxqidj(1)-2.0d0*qixqj(1))
     &            - gf(5)*rxqir(1) - gf(7)*(rxqijr(1)+qjrxqir(1))
               ttm2(2) = -bn(1)*dixdj(2) + gf(2)*dixr(2)
     &            + gf(4)*(dixqjr(2)+djxqir(2)+rxqidj(2)-2.0d0*qixqj(2))
     &            - gf(5)*rxqir(2) - gf(7)*(rxqijr(2)+qjrxqir(2))
               ttm2(3) = -bn(1)*dixdj(3) + gf(2)*dixr(3)
     &            + gf(4)*(dixqjr(3)+djxqir(3)+rxqidj(3)-2.0d0*qixqj(3))
     &            - gf(5)*rxqir(3) - gf(7)*(rxqijr(3)+qjrxqir(3))
               ttm3(1) = bn(1)*dixdj(1) + gf(3)*djxr(1)
     &            - gf(4)*(dixqjr(1)+djxqir(1)+rxqjdi(1)-2.0d0*qixqj(1))
     &            - gf(6)*rxqjr(1) - gf(7)*(rxqjir(1)-qjrxqir(1))
               ttm3(2) = bn(1)*dixdj(2) + gf(3)*djxr(2)
     &            - gf(4)*(dixqjr(2)+djxqir(2)+rxqjdi(2)-2.0d0*qixqj(2))
     &            - gf(6)*rxqjr(2) - gf(7)*(rxqjir(2)-qjrxqir(2))
               ttm3(3) = bn(1)*dixdj(3) +gf(3)*djxr(3)
     &            - gf(4)*(dixqjr(3)+djxqir(3)+rxqjdi(3)-2.0d0*qixqj(3))
     &            - gf(6)*rxqjr(3) - gf(7)*(rxqjir(3)-qjrxqir(3))
c
c     get the real permanent torques
c
               ttm2r(1) = -rr3*dixdj(1) + gfr(2)*dixr(1)-gfr(5)*rxqir(1)
     &           + gfr(4)*(dixqjr(1)+djxqir(1)+rxqidj(1)-2.0d0*qixqj(1))
     &           - gfr(7)*(rxqijr(1)+qjrxqir(1))
               ttm2r(2) = -rr3*dixdj(2) + gfr(2)*dixr(2)-gfr(5)*rxqir(2)
     &           + gfr(4)*(dixqjr(2)+djxqir(2)+rxqidj(2)-2.0d0*qixqj(2))
     &           - gfr(7)*(rxqijr(2)+qjrxqir(2))
               ttm2r(3) = -rr3*dixdj(3) + gfr(2)*dixr(3)-gfr(5)*rxqir(3)
     &           + gfr(4)*(dixqjr(3)+djxqir(3)+rxqidj(3)-2.0d0*qixqj(3))
     &           - gfr(7)*(rxqijr(3)+qjrxqir(3))
               ttm3r(1) = rr3*dixdj(1) + gfr(3)*djxr(1) -gfr(6)*rxqjr(1)
     &           - gfr(4)*(dixqjr(1)+djxqir(1)+rxqjdi(1)-2.0d0*qixqj(1))
     &           - gfr(7)*(rxqjir(1)-qjrxqir(1))
               ttm3r(2) = rr3*dixdj(2) + gfr(3)*djxr(2) -gfr(6)*rxqjr(2)
     &           - gfr(4)*(dixqjr(2)+djxqir(2)+rxqjdi(2)-2.0d0*qixqj(2))
     &           - gfr(7)*(rxqjir(2)-qjrxqir(2))
               ttm3r(3) = rr3*dixdj(3) + gfr(3)*djxr(3) -gfr(6)*rxqjr(3)
     &           - gfr(4)*(dixqjr(3)+djxqir(3)+rxqjdi(3)-2.0d0*qixqj(3))
     &           - gfr(7)*(rxqjir(3)-qjrxqir(3))
c
c     handle the case where scaling is used
c
               ttm2(1) = ttm2(1) - ttm2r(1) * (1.0d0-mscale(jj))
               ttm2(2) = ttm2(2) - ttm2r(2) * (1.0d0-mscale(jj))
               ttm2(3) = ttm2(3) - ttm2r(3) * (1.0d0-mscale(jj))
               ttm3(1) = ttm3(1) - ttm3r(1) * (1.0d0-mscale(jj))
               ttm3(2) = ttm3(2) - ttm3r(2) * (1.0d0-mscale(jj))
               ttm3(3) = ttm3(3) - ttm3r(3) * (1.0d0-mscale(jj))
               ttm2i(1) = ttm2i(1) - ttm2ri(1)
               ttm2i(2) = ttm2i(2) - ttm2ri(2)
               ttm2i(3) = ttm2i(3) - ttm2ri(3)
               ttm3i(1) = ttm3i(1) - ttm3ri(1)
               ttm3i(2) = ttm3i(2) - ttm3ri(2)
               ttm3i(3) = ttm3i(3) - ttm3ri(3)
c
c     now compute the virial components for this interaction
c
               do k = 1, 3
                  frczi(k) = 0.0d0
                  frcxi(k) = 0.0d0
                  frczj(k) = 0.0d0
                  frcxj(k) = 0.0d0
               end do
               xiz = x(zaxis(i)) - x(ii)
               yiz = y(zaxis(i)) - y(ii)
               ziz = z(zaxis(i)) - z(ii)
               xix = x(xaxis(i)) - x(ii)
               yix = y(xaxis(i)) - y(ii)
               zix = z(xaxis(i)) - z(ii)
               xjz = x(zaxis(j)) - x(jj)
               yjz = y(zaxis(j)) - y(jj)
               zjz = z(zaxis(j)) - z(jj)
               xjx = x(xaxis(j)) - x(jj)
               yjx = y(xaxis(j)) - y(jj)
               zjx = z(xaxis(j)) - z(jj)
               do k = 1, 3
                  tq(k) = ttm2(k) + ttm2i(k)
               end do
               call torque1 (i,tq,frczi,frcxi)
               do k = 1, 3
                  tq(k) = ttm3(k) + ttm3i(k)
               end do
               call torque1 (j,tq,frczj,frcxj)
               vxx = -xji*(ftm2(1)+ftm2i(1)) + xix*frcxi(1)
     &                  + xjz*frczj(1) + xjx*frcxj(1) + xiz*frczi(1)
               vyx = -yji*(ftm2(1)+ftm2i(1)) + yix*frcxi(1)
     &                  + yjz*frczj(1) + yjx*frcxj(1) + yiz*frczi(1)
               vzx = -zji*(ftm2(1)+ftm2i(1)) + zix*frcxi(1)
     &                  + zjz*frczj(1) + zjx*frcxj(1) + ziz*frczi(1)
               vxy = -xji*(ftm2(2)+ftm2i(2)) + xix*frcxi(2)
     &                  + xjz*frczj(2) + xjx*frcxj(2) + xiz*frczi(2)
               vyy = -yji*(ftm2(2)+ftm2i(2)) + yix*frcxi(2)
     &                  + yjz*frczj(2) + yjx*frcxj(2) + yiz*frczi(2)
               vzy = -zji*(ftm2(2)+ftm2i(2)) + zix*frcxi(2)
     &                  + zjz*frczj(2) + zjx*frcxj(2) + ziz*frczi(2)
               vxz = -xji*(ftm2(3)+ftm2i(3)) + xix*frcxi(3)
     &                  + xjz*frczj(3) + xjx*frcxj(3) + xiz*frczi(3)
               vyz = -yji*(ftm2(3)+ftm2i(3)) + yix*frcxi(3)
     &                  + yjz*frczj(3) + yjx*frcxj(3) + yiz*frczi(3)
               vzz = -zji*(ftm2(3)+ftm2i(3)) + zix*frcxi(3)
     &                  + zjz*frczj(3) + zjx*frcxj(3) + ziz*frczi(3)
               vir(1,1) = vir(1,1) + f*vxx
               vir(2,1) = vir(2,1) + 0.5d0*f*(vyx+vxy)
               vir(3,1) = vir(3,1) + 0.5d0*f*(vzx+vxz)
               vir(1,2) = vir(1,2) + 0.5d0*f*(vxy+vyx)
               vir(2,2) = vir(2,2) + f*vyy
               vir(3,2) = vir(3,2) + 0.5d0*f*(vzy+vyz)
               vir(1,3) = vir(1,3) + 0.5d0*f*(vxz+vzx)
               vir(2,3) = vir(2,3) + 0.5d0*f*(vyz+vzy)
               vir(3,3) = vir(3,3) + f*vzz
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
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction with other unit cells
c
      do i = 1, npole
         ii = ipole(i)
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
         do j = 1, n
            mscale(j) = 1.0d0
            dscale(j) = 1.0d0
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
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
         end do
         do j = i, npole
            jj = ipole(j)
            do kn = 1, ncell
            if (.not. (use_polymer .and. r2.le.polycut2)) then
               mscale(jj) = 1.0d0
               pscale(jj) = 1.0d0
               dscale(jj) = 1.0d0
               uscale(jj) = 1.0d0
            end if
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
            call image (xji,yji,zji,kn)
            r2 = xji*xji + yji*yji + zji*zji
            if (r2 .le. off2) then
               r = sqrt(r2)
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
c
c     calculate the real space error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do m = 1, 5
                  bfac = dble(m+m-1)
                  alsq2n = alsq2 * alsq2n
                  bn(m) = (bfac*bn(m-1)+alsq2n*exp2a) / r2
               end do
c
c     zero out the temporary force and torque
c
               do k = 1, 3
                  ftm2(k) = 0.0d0
                  ftm2i(k) = 0.0d0
                  ttm2(k) = 0.0d0
                  ttm2i(k) = 0.0d0
                  ttm3(k) = 0.0d0
                  ttm3i(k) = 0.0d0
               end do
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               scale9 = 1.0d0
               do k = 1, 3
                  ddsc3(k) = 0.0d0
                  ddsc5(k) = 0.0d0
                  ddsc7(k) = 0.0d0
               end do
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
                     ddsc3(1) = -3.0d0*damp*exp(damp) * xji/r2
                     ddsc3(2) = -3.0d0*damp*exp(damp) * yji/r2
                     ddsc3(3) = -3.0d0*damp*exp(damp) * zji/r2
                     ddsc5(1) = -damp * ddsc3(1)
                     ddsc5(2) = -damp * ddsc3(2)
                     ddsc5(3) = -damp * ddsc3(3)
                     ddsc7(1) = (-0.2d0-0.6d0*damp) * ddsc5(1)
                     ddsc7(2) = (-0.2d0-0.6d0*damp) * ddsc5(2)
                     ddsc7(3) = (-0.2d0-0.6d0*damp) * ddsc5(3)
                  end if
               end if
               dsc3 = 1.0d0 - scale3*dscale(jj)
               dsc5 = 1.0d0 - scale5*dscale(jj)
               dsc7 = 1.0d0 - scale7*dscale(jj)
               psc3 = 1.0d0 - scale3*pscale(jj)
               psc5 = 1.0d0 - scale5*pscale(jj)
               psc7 = 1.0d0 - scale7*pscale(jj)
               ruscale3 = 1.0d0 - scale3*uscale(jj)
               ruscale5 = 1.0d0 - scale5*uscale(jj)
               ruscale7 = 1.0d0 - scale7*uscale(jj)
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
               dixujp(1) = di(2)*uinp(3,j) - di(3)*uinp(2,j)
               dixujp(2) = di(3)*uinp(1,j) - di(1)*uinp(3,j)
               dixujp(3) = di(1)*uinp(2,j) - di(2)*uinp(1,j)
               djxuip(1) = dj(2)*uinp(3,i) - dj(3)*uinp(2,i)
               djxuip(2) = dj(3)*uinp(1,i) - dj(1)*uinp(3,i)
               djxuip(3) = dj(1)*uinp(2,i) - dj(2)*uinp(1,i)
               dixr(1) = di(2)*zji - di(3)*yji
               dixr(2) = di(3)*xji - di(1)*zji
               dixr(3) = di(1)*yji - di(2)*xji
               djxr(1) = dj(2)*zji - dj(3)*yji
               djxr(2) = dj(3)*xji - dj(1)*zji
               djxr(3) = dj(1)*yji - dj(2)*xji
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
               qiujp(1)= qi(1)*uinp(1,j)+qi(4)*uinp(2,j)+qi(7)*uinp(3,j)
               qiujp(2)= qi(2)*uinp(1,j)+qi(5)*uinp(2,j)+qi(8)*uinp(3,j)
               qiujp(3)= qi(3)*uinp(1,j)+qi(6)*uinp(2,j)+qi(9)*uinp(3,j)
               qjuip(1)= qj(1)*uinp(1,i)+qj(4)*uinp(2,i)+qj(7)*uinp(3,i)
               qjuip(2)= qj(2)*uinp(1,i)+qj(5)*uinp(2,i)+qj(8)*uinp(3,i)
               qjuip(3)= qj(3)*uinp(1,i)+qj(6)*uinp(2,i)+qj(9)*uinp(3,i)
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
               uixqjrp(1) = uinp(2,i)*qjr(3) - uinp(3,i)*qjr(2)
               uixqjrp(2) = uinp(3,i)*qjr(1) - uinp(1,i)*qjr(3)
               uixqjrp(3) = uinp(1,i)*qjr(2) - uinp(2,i)*qjr(1)
               ujxqirp(1) = uinp(2,j)*qir(3) - uinp(3,j)*qir(2)
               ujxqirp(2) = uinp(3,j)*qir(1) - uinp(1,j)*qir(3)
               ujxqirp(3) = uinp(1,j)*qir(2) - uinp(2,j)*qir(1)
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
               rxqiujp(1) = yji*qiujp(3) - zji*qiujp(2)
               rxqiujp(2) = zji*qiujp(1) - xji*qiujp(3)
               rxqiujp(3) = xji*qiujp(2) - yji*qiujp(1)
               rxqjuip(1) = yji*qjuip(3) - zji*qjuip(2)
               rxqjuip(2) = zji*qjuip(1) - xji*qjuip(3)
               rxqjuip(3) = xji*qjuip(2) - yji*qjuip(1)
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
c     get intermediate variables for induction energy;
c     interaction between induce and permanent multipoles
c     final induction energy exclude induce-induce term
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
               scip(1) = uinp(1,i)*dj(1) + uinp(2,i)*dj(2)
     &                      + uinp(3,i)*dj(3) + di(1)*uinp(1,j)
     &                      + di(2)*uinp(2,j) + di(3)*uinp(3,j)
               scip(2) = uind(1,i)*uinp(1,j)+uind(2,i)*uinp(2,j)
     &                      + uind(3,i)*uinp(3,j)+uinp(1,i)*uind(1,j)
     &                      + uinp(2,i)*uind(2,j)+uinp(3,i)*uind(3,j)
               scip(3) = uinp(1,i)*xji + uinp(2,i)*yji + uinp(3,i)*zji
               scip(4) = uinp(1,j)*xji + uinp(2,j)*yji + uinp(3,j)*zji
               scip(7) = qir(1)*uinp(1,j) + qir(2)*uinp(2,j)
     &                      + qir(3)*uinp(3,j)
               scip(8) = qjr(1)*uinp(1,i) + qjr(2)*uinp(2,i)
     &                      + qjr(3)*uinp(3,i)
c
c     get the induced-induced term
c
               einin = sci(2)*rr3 - sci(3)*sci(4)*rr5
               findmp(1) = uscale(jj) * (scip(2)*rr3*ddsc3(1)
     &          -(sci(3)*scip(4)+scip(3)*sci(4))*rr5*ddsc5(1))*0.5d0
               findmp(2) = uscale(jj) * (scip(2)*rr3*ddsc3(2)
     &          -(sci(3)*scip(4)+scip(3)*sci(4))*rr5*ddsc5(2))*0.5d0
               findmp(3) = uscale(jj) * (scip(2)*rr3*ddsc3(3)
     &          -(sci(3)*scip(4)+scip(3)*sci(4))*rr5*ddsc5(3))*0.5d0
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
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = cj*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     calculate energy for this interaction and increment total
c
               e = bn(0)*gl(0) + bn(1)*(gl(1)+gl(6))
     &                + bn(2)*(gl(2)+gl(7)+gl(8))
     &                + bn(3)*(gl(3)+gl(5)) + bn(4)*gl(4)
               ei = 0.5d0 * (bn(1)*(gli(1)+gli(6))
     &                      + bn(2)*(gli(2)+gli(7)) + bn(3)*gli(3))
c
c    real energy without screening function
c
               erl = rr1*gl(0) + rr3*(gl(1)+gl(6))
     &                  + rr5*(gl(2)+gl(7)+gl(8))
     &                  + rr7*(gl(3)+gl(5)) + rr9*gl(4)
               erli = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &                   + rr5*(gli(2)+gli(7))*psc5
     &                   + rr7*gli(3)*psc7)
               e = e - (1.0d0-mscale(jj))*erl
               ei = ei - erli
               e = f * e
               ei = f * ei
               if (ii .eq. jj) then
                  e = 0.5d0 * e
                  ei = 0.5d0 * ei
               end if
               em = em + e
               ep = ep + ei
c
c
c
               fridmp(1) = 0.5d0*(rr3*((gli(1)+gli(6))*pscale(jj)
     &           + (glip(1)+glip(6))*dscale(jj))*ddsc3(1)
     &           + rr5*((gli(2)+gli(7))*pscale(jj)
     &           + (glip(2)+glip(7))*dscale(jj))*ddsc5(1)
     &           + rr7*(gli(3)*pscale(jj)+glip(3)*dscale(jj))*ddsc7(1))
               fridmp(2) = 0.5d0*(rr3*((gli(1)+gli(6))*pscale(jj)
     &           + (glip(1)+glip(6))*dscale(jj))*ddsc3(2)
     &           + rr5*((gli(2)+gli(7))*pscale(jj)
     &           + (glip(2)+glip(7))*dscale(jj))*ddsc5(2)
     &           + rr7*(gli(3)*pscale(jj)+glip(3)*dscale(jj))*ddsc7(2))
               fridmp(3) = 0.5d0*(rr3*((gli(1)+gli(6))*pscale(jj)
     &           + (glip(1)+glip(6))*dscale(jj))*ddsc3(3)
     &           + rr5*((gli(2)+gli(7))*pscale(jj)
     &           + (glip(2)+glip(7))*dscale(jj))*ddsc5(3)
     &           + rr7*(gli(3)*pscale(jj)+glip(3)*dscale(jj))*ddsc7(3))
c
c     increment the intramolecular energy if appropriate;
c     assumes intramolecular distances are less than half
c     of cell length and less than the ewald cutoff
c
               if (molcule(ii) .eq. molcule(jj)) then
                  eintra = eintra + mscale(jj)*erl*f
                  eintra = eintra + 0.5d0*pscale(jj)
     &                        * (rr3*(gli(1)+gli(6))*scale3
     &                              + rr5*(gli(2)+gli(7))*scale5
     &                              + rr7*gli(3)*scale7)
               end if
c
c     intermediate variables for permanent multipole force terms
c
               gf(1) = bn(1)*gl(0) + bn(2)*(gl(1)+gl(6))
     &                    + bn(3)*(gl(2)+gl(7)+gl(8))
     &                    + bn(4)*(gl(3)+gl(5)) + bn(5)*gl(4)
               gf(2) = -cj*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gf(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gf(4) = 2.0d0 * bn(2)
               gf(5) = 2.0d0 * (-cj*bn(2)+sc(4)*bn(3)-sc(6)*bn(4))
               gf(6) = 2.0d0 * (-ci*bn(2)-sc(3)*bn(3)-sc(5)*bn(4))
               gf(7) = 4.0d0 * bn(3)
               gfr(1) = rr3*gl(0) + rr5*(gl(1)+gl(6))
     &                     + rr7*(gl(2)+gl(7)+gl(8))
     &                     + rr9*(gl(3)+gl(5)) + rr11*gl(4)
               gfr(2) = -cj*rr3 + sc(4)*rr5 - sc(6)*rr7
               gfr(3) = ci*rr3 + sc(3)*rr5 + sc(5)*rr7
               gfr(4) = 2.0d0 * rr5
               gfr(5) = 2.0d0 * (-cj*rr5+sc(4)*rr7-sc(6)*rr9)
               gfr(6) = 2.0d0 * (-ci*rr5-sc(3)*rr7-sc(5)*rr9)
               gfr(7) = 4.0d0 * rr7
c
c     get the induced-permanent term
c
               gfi(1) = bn(2)*(gli(1)+glip(1)+gli(6)+glip(6))*0.5d0
     &                + bn(2)* scip(2)*0.5d0
     &                + bn(3)*(gli(2)+glip(2)+gli(7)+glip(7))*0.5d0
     &                - bn(3)*0.5d0*(sci(3)*scip(4)+scip(3)*sci(4))
     &                + bn(4)*(gli(3)+glip(3))*0.5d0
               gfi(2) = -cj*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gfi(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gfi(4) = 2.0d0 * bn(2)
               gfi(5) = (sci(4) + scip(4))*bn(3)
               gfi(6) = -(sci(3) + scip(3))*bn(3)
               gfri(1) = rr5*0.5d0*((gli(1)+gli(6))*psc3
     &             + (glip(1)+glip(6))*dsc3 + scip(2)*ruscale3)
     &             + rr7*((gli(7)+gli(2))*psc5 + (glip(7)+glip(2))*dsc5
     &             - (sci(3)*scip(4)+scip(3)*sci(4))*ruscale5)*0.5d0
     &             + rr9*(gli(3)*psc7+glip(3)*dsc7)*0.5d0
               gfri(2) = -rr3*cj + rr5*sc(4) - rr7*sc(6)
               gfri(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfri(4) = 2.0d0 * rr5
               gfri(5) = rr7 * (sci(4)*psc7 + scip(4)*dsc7)
               gfri(6) = -rr7 * (sci(3)*psc7 + scip(3)*dsc7)
c
c       get the induced-induced terms for direct versus mutual
c
               gfdr = 0.5d0 * (rr5*scip(2)*ruscale3
     &                  - rr7*(scip(3)*sci(4)
     &                        +sci(3)*scip(4))*ruscale5)
               gfd = 0.5d0 * (bn(2)*scip(2)
     &                  - bn(3)*(scip(3)*sci(4)+sci(3)*scip(4)))
c
c     get the unscreened induced force components
c
               ftm2ri(1) = gfri(1)*xji + 0.5d0*
     &          (- rr3*cj*(uind(1,i)*psc3+uinp(1,i)*dsc3)
     &           + rr5*sc(4)*(uind(1,i)*psc5+uinp(1,i)*dsc5)
     &           - rr7*sc(6)*(uind(1,i)*psc7+uinp(1,i)*dsc7))
     &           +(rr3*ci*(uind(1,j)*psc3+uinp(1,j)*dsc3)
     &           + rr5*sc(3)*(uind(1,j)*psc5+uinp(1,j)*dsc5)
     &           + rr7*sc(5)*(uind(1,j)*psc7+uinp(1,j)*dsc7))*0.5d0
     &           + rr5*ruscale5*(sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &           + sci(3)*uinp(1,j)+scip(3)*uind(1,j))*0.5d0
     &           + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &           + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dj(1)
     &           + 0.5d0*gfri(4)*((qjui(1)-qiuj(1))*psc5
     &           + (qjuip(1)-qiujp(1))*dsc5)
     &           + gfri(5)*qir(1) + gfri(6)*qjr(1)
               ftm2ri(2) = gfri(1)*yji + 0.5d0*
     &          (- rr3*cj*(uind(2,i)*psc3+uinp(2,i)*dsc3)
     &           + rr5*sc(4)*(uind(2,i)*psc5+uinp(2,i)*dsc5)
     &           - rr7*sc(6)*(uind(2,i)*psc7+uinp(2,i)*dsc7))
     &           +(rr3*ci*(uind(2,j)*psc3+uinp(2,j)*dsc3)
     &           + rr5*sc(3)*(uind(2,j)*psc5+uinp(2,j)*dsc5)
     &           + rr7*sc(5)*(uind(2,j)*psc7+uinp(2,j)*dsc7))*0.5d0
     &           + rr5*ruscale5*(sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &           + sci(3)*uinp(2,j)+scip(3)*uind(2,j))*0.5d0
     &           + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &           + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dj(2)
     &           + 0.5d0*gfri(4)*((qjui(2)-qiuj(2))*psc5
     &           + (qjuip(2)-qiujp(2))*dsc5)
     &           + gfri(5)*qir(2) + gfri(6)*qjr(2)
               ftm2ri(3) = gfri(1)*zji  + 0.5d0*
     &          (- rr3*cj*(uind(3,i)*psc3+uinp(3,i)*dsc3)
     &           + rr5*sc(4)*(uind(3,i)*psc5+uinp(3,i)*dsc5)
     &           - rr7*sc(6)*(uind(3,i)*psc7+uinp(3,i)*dsc7))
     &           +(rr3*ci*(uind(3,j)*psc3+uinp(3,j)*dsc3)
     &           + rr5*sc(3)*(uind(3,j)*psc5+uinp(3,j)*dsc5)
     &           + rr7*sc(5)*(uind(3,j)*psc7+uinp(3,j)*dsc7))*0.5d0
     &           + rr5*ruscale5*(sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &           + sci(3)*uinp(3,j)+scip(3)*uind(3,j))*0.5d0
     &           + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &           + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dj(3)
     &           + 0.5d0*gfri(4)*((qjui(3)-qiuj(3))*psc5
     &           + (qjuip(3)-qiujp(3))*dsc5)
     &           + gfri(5)*qir(3) + gfri(6)*qjr(3)
c
c     get the real-space induced force
c
               ftm2i(1)=gfi(1)*xji + 0.5d0*(gfi(2)*(uind(1,i)+uinp(1,i))
     &              + bn(2)*(sci(4)*uinp(1,i)+scip(4)*uind(1,i))
     &              + gfi(3)*(uind(1,j)+uinp(1,j))
     &              + bn(2)*(sci(3)*uinp(1,j)+scip(3)*uind(1,j))
     &              + (sci(4)+scip(4))*bn(2)*di(1)
     &              + (sci(3)+scip(3))*bn(2)*dj(1)
     &              + gfi(4)*(qjui(1)+qjuip(1)-qiuj(1)-qiujp(1)))
     &              + gfi(5)*qir(1) + gfi(6)*qjr(1)
               ftm2i(2)=gfi(1)*yji + 0.5d0*(gfi(2)*(uind(2,i)+uinp(2,i))
     &              + bn(2)*(sci(4)*uinp(2,i)+scip(4)*uind(2,i))
     &              + gfi(3)*(uind(2,j)+uinp(2,j))
     &              + bn(2)*(sci(3)*uinp(2,j)+scip(3)*uind(2,j))
     &              + (sci(4)+scip(4))*bn(2)*di(2)
     &              + (sci(3)+scip(3))*bn(2)*dj(2)
     &              + gfi(4)*(qjui(2)+qjuip(2)-qiuj(2)-qiujp(2)))
     &              + gfi(5)*qir(2) + gfi(6)*qjr(2)
               ftm2i(3)=gfi(1)*zji + 0.5d0*(gfi(2)*(uind(3,i)+uinp(3,i))
     &              + bn(2)*(sci(4)*uinp(3,i)+scip(4)*uind(3,i))
     &              + gfi(3)*(uind(3,j)+uinp(3,j))
     &              + bn(2)*(sci(3)*uinp(3,j)+scip(3)*uind(3,j))
     &              + (sci(4)+scip(4))*bn(2)*di(3)
     &              + (sci(3)+scip(3))*bn(2)*dj(3)
     &              + gfi(4)*(qjui(3)+qjuip(3)-qiuj(3)-qiujp(3)))
     &              + gfi(5)*qir(3) + gfi(6)*qjr(3)
c
c     get the real permanent force without screening function
c
               ftm2r(1) = gfr(1)*xji + gfr(2)*di(1) + gfr(3)*dj(1)
     &                     + gfr(4)*(qjdi(1)-qidj(1)) + gfr(5)*qir(1)
     &                     + gfr(6)*qjr(1) + gfr(7)*(qiqjr(1)+qjqir(1))
               ftm2r(2) = gfr(1)*yji + gfr(2)*di(2) + gfr(3)*dj(2)
     &                     + gfr(4)*(qjdi(2)-qidj(2)) + gfr(5)*qir(2)
     &                     + gfr(6)*qjr(2) + gfr(7)*(qiqjr(2)+qjqir(2))
               ftm2r(3) = gfr(1)*zji + gfr(2)*di(3) + gfr(3)*dj(3)
     &                     + gfr(4)*(qjdi(3)-qidj(3)) + gfr(5)*qir(3)
     &                     + gfr(6)*qjr(3) + gfr(7)*(qiqjr(3)+qjqir(3))
c
c     get the real-space permanent force
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
               if (poltyp .eq. 'DIRECT') then
                  ftm2i(1) = ftm2i(1) - gfd*xji - 0.5d0*
     &               (bn(2)*(sci(4)*uinp(1,i)+scip(4)*uind(1,i))
     &               + bn(2)*(sci(3)*uinp(1,j)+scip(3)*uind(1,j)))
                  ftm2i(2) = ftm2i(2) - gfd*yji -0.5d0*
     &               (bn(2)*(sci(4)*uinp(2,i)+scip(4)*uind(2,i))
     &               + bn(2)*(sci(3)*uinp(2,j)+scip(3)*uind(2,j)))
                  ftm2i(3) = ftm2i(3) - gfd*zji -0.5d0*
     &               (bn(2)*(sci(4)*uinp(3,i)+scip(4)*uind(3,i))
     &               + bn(2)*(sci(3)*uinp(3,j)+scip(3)*uind(3,j)))
                  fdir(1) = gfdr*xji + 0.5d0*ruscale5*rr5*
     &               (sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &               + sci(3)*uinp(1,j)+scip(3)*uind(1,j))
                  fdir(2) = gfdr*yji + 0.5d0*ruscale5*rr5*
     &               (sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &               + sci(3)*uinp(2,j)+scip(3)*uind(2,j))
                  fdir(3) = gfdr*zji + 0.5d0*ruscale5*rr5*
     &               (sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &               + sci(3)*uinp(3,j)+scip(3)*uind(3,j))
               end if
c
c     handle of scaling for partially excluded interactions
c
               ftm2(1) = ftm2(1) - (1.0d0-mscale(jj))*ftm2r(1)
               ftm2(2) = ftm2(2) - (1.0d0-mscale(jj))*ftm2r(2)
               ftm2(3) = ftm2(3) - (1.0d0-mscale(jj))*ftm2r(3)
               ftm2i(1) = ftm2i(1) - ftm2ri(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - ftm2ri(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - ftm2ri(3) - fridmp(3) - findmp(3)
               if (poltyp .eq. 'DIRECT') then
                  ftm2i(1) = ftm2i(1) + fdir(1) + findmp(1)
                  ftm2i(2) = ftm2i(2) + fdir(2) + findmp(2)
                  ftm2i(3) = ftm2i(3) + fdir(3) + findmp(3)
               end if
c
c     now perform the torque calculation
c     intermediate terms for torque between multipoles i and j
c
               gtri(2) = 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5
               gtri(3) = 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5
               gtri(4) = gfri(4)
               gtri(5) = gfri(5)
               gtri(6) = gfri(6)
               gti(2) = 0.5d0*(sci(4)+scip(4))*bn(2)
               gti(3) = 0.5d0*(sci(3)+scip(3))*bn(2)
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
c
c     get the unscreened induced torque
c
               ttm2ri(1) = -rr3*(dixuj(1)*psc3+dixujp(1)*dsc3)*0.5d0
     &           + gtri(2)*dixr(1) + gtri(4)*((ujxqir(1)+rxqiuj(1))*psc5
     &           +(ujxqirp(1)+rxqiujp(1))*dsc5)*0.5d0 - gtri(5)*rxqir(1)
               ttm2ri(2) = -rr3*(dixuj(2)*psc3+dixujp(2)*dsc3)*0.5d0
     &           + gtri(2)*dixr(2) + gtri(4)*((ujxqir(2)+rxqiuj(2))*psc5
     &           +(ujxqirp(2)+rxqiujp(2))*dsc5)*0.5d0 - gtri(5)*rxqir(2)
               ttm2ri(3) = -rr3*(dixuj(3)*psc3+dixujp(3)*dsc3)*0.5d0
     &           + gtri(2)*dixr(3) + gtri(4)*((ujxqir(3)+rxqiuj(3))*psc5
     &           +(ujxqirp(3)+rxqiujp(3))*dsc5)*0.5d0 - gtri(5)*rxqir(3)
               ttm3ri(1) = -rr3*(djxui(1)*psc3+djxuip(1)*dsc3)*0.5d0
     &           + gtri(3)*djxr(1) - gtri(4)*((uixqjr(1)+rxqjui(1))*psc5
     &           +(uixqjrp(1)+rxqjuip(1))*dsc5)*0.5d0 - gtri(6)*rxqjr(1)
               ttm3ri(2) = -rr3*(djxui(2)*psc3+djxuip(2)*dsc3)*0.5d0
     &           + gtri(3)*djxr(2) - gtri(4)*((uixqjr(2)+rxqjui(2))*psc5
     &           +(uixqjrp(2)+rxqjuip(2))*dsc5)*0.5d0 - gtri(6)*rxqjr(2)
               ttm3ri(3) = -rr3*(djxui(3)*psc3+djxuip(3)*dsc3)*0.5d0
     &           + gtri(3)*djxr(3) - gtri(4)*((uixqjr(3)+rxqjui(3))*psc5
     &           +(uixqjrp(3)+rxqjuip(3))*dsc5)*0.5d0 - gtri(6)*rxqjr(3)
c
c     get the real-space induced torque
c
               ttm2i(1) = -bn(1)*(dixuj(1)+dixujp(1))*0.5d0
     &              + gti(2)*dixr(1) + gti(4)*(ujxqir(1)+rxqiuj(1)
     &              + ujxqirp(1)+rxqiujp(1))*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -bn(1)*(dixuj(2)+dixujp(2))*0.5d0
     &              + gti(2)*dixr(2) + gti(4)*(ujxqir(2)+rxqiuj(2)
     &              + ujxqirp(2)+rxqiujp(2))*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -bn(1)*(dixuj(3)+dixujp(3))*0.5d0
     &              + gti(2)*dixr(3) + gti(4)*(ujxqir(3)+rxqiuj(3)
     &              + ujxqirp(3)+rxqiujp(3))*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -bn(1)*(djxui(1)+djxuip(1))*0.5d0
     &              + gti(3)*djxr(1) - gti(4)*(uixqjr(1)+rxqjui(1)
     &              + uixqjrp(1)+rxqjuip(1))*0.5d0 - gti(6)*rxqjr(1)
               ttm3i(2) = -bn(1)*(djxui(2)+djxuip(2))*0.5d0
     &              + gti(3)*djxr(2) - gti(4)*(uixqjr(2)+rxqjui(2)
     &              + uixqjrp(2)+rxqjuip(2))*0.5d0 - gti(6)*rxqjr(2)
               ttm3i(3) = -bn(1)*(djxui(3)+djxuip(3))*0.5d0
     &              + gti(3)*djxr(3) - gti(4)*(uixqjr(3)+rxqjui(3)
     &              + uixqjrp(3)+rxqjuip(3))*0.5d0 - gti(6)*rxqjr(3)
c
c     get the real-space permanent torques
c
               ttm2(1) = -bn(1)*dixdj(1) + gf(2)*dixr(1)
     &            + gf(4)*(dixqjr(1)+djxqir(1)+rxqidj(1)-2.0d0*qixqj(1))
     &            - gf(5)*rxqir(1) - gf(7)*(rxqijr(1)+qjrxqir(1))
               ttm2(2) = -bn(1)*dixdj(2) + gf(2)*dixr(2)
     &            + gf(4)*(dixqjr(2)+djxqir(2)+rxqidj(2)-2.0d0*qixqj(2))
     &            - gf(5)*rxqir(2) - gf(7)*(rxqijr(2)+qjrxqir(2))
               ttm2(3) = -bn(1)*dixdj(3) + gf(2)*dixr(3)
     &            + gf(4)*(dixqjr(3)+djxqir(3)+rxqidj(3)-2.0d0*qixqj(3))
     &            - gf(5)*rxqir(3) - gf(7)*(rxqijr(3)+qjrxqir(3))
               ttm3(1) = bn(1)*dixdj(1) + gf(3)*djxr(1)
     &            - gf(4)*(dixqjr(1)+djxqir(1)+rxqjdi(1)-2.0d0*qixqj(1))
     &            - gf(6)*rxqjr(1) - gf(7)*(rxqjir(1)-qjrxqir(1))
               ttm3(2) = bn(1)*dixdj(2) + gf(3)*djxr(2)
     &            - gf(4)*(dixqjr(2)+djxqir(2)+rxqjdi(2)-2.0d0*qixqj(2))
     &            - gf(6)*rxqjr(2) - gf(7)*(rxqjir(2)-qjrxqir(2))
               ttm3(3) = bn(1)*dixdj(3) +gf(3)*djxr(3)
     &            - gf(4)*(dixqjr(3)+djxqir(3)+rxqjdi(3)-2.0d0*qixqj(3))
     &            - gf(6)*rxqjr(3) - gf(7)*(rxqjir(3)-qjrxqir(3))
c
c     get the real permanent torques
c
               ttm2r(1) = -rr3*dixdj(1) + gfr(2)*dixr(1)-gfr(5)*rxqir(1)
     &           + gfr(4)*(dixqjr(1)+djxqir(1)+rxqidj(1)-2.0d0*qixqj(1))
     &           - gfr(7)*(rxqijr(1)+qjrxqir(1))
               ttm2r(2) = -rr3*dixdj(2) + gfr(2)*dixr(2)-gfr(5)*rxqir(2)
     &           + gfr(4)*(dixqjr(2)+djxqir(2)+rxqidj(2)-2.0d0*qixqj(2))
     &           - gfr(7)*(rxqijr(2)+qjrxqir(2))
               ttm2r(3) = -rr3*dixdj(3) + gfr(2)*dixr(3)-gfr(5)*rxqir(3)
     &           + gfr(4)*(dixqjr(3)+djxqir(3)+rxqidj(3)-2.0d0*qixqj(3))
     &           - gfr(7)*(rxqijr(3)+qjrxqir(3))
               ttm3r(1) = rr3*dixdj(1) + gfr(3)*djxr(1) -gfr(6)*rxqjr(1)
     &           - gfr(4)*(dixqjr(1)+djxqir(1)+rxqjdi(1)-2.0d0*qixqj(1))
     &           - gfr(7)*(rxqjir(1)-qjrxqir(1))
               ttm3r(2) = rr3*dixdj(2) + gfr(3)*djxr(2) -gfr(6)*rxqjr(2)
     &           - gfr(4)*(dixqjr(2)+djxqir(2)+rxqjdi(2)-2.0d0*qixqj(2))
     &           - gfr(7)*(rxqjir(2)-qjrxqir(2))
               ttm3r(3) = rr3*dixdj(3) + gfr(3)*djxr(3) -gfr(6)*rxqjr(3)
     &           - gfr(4)*(dixqjr(3)+djxqir(3)+rxqjdi(3)-2.0d0*qixqj(3))
     &           - gfr(7)*(rxqjir(3)-qjrxqir(3))
c
c     handle the case where scaling is used
c
               ttm2(1) = ttm2(1) - ttm2r(1) * (1.0d0-mscale(jj))
               ttm2(2) = ttm2(2) - ttm2r(2) * (1.0d0-mscale(jj))
               ttm2(3) = ttm2(3) - ttm2r(3) * (1.0d0-mscale(jj))
               ttm3(1) = ttm3(1) - ttm3r(1) * (1.0d0-mscale(jj))
               ttm3(2) = ttm3(2) - ttm3r(2) * (1.0d0-mscale(jj))
               ttm3(3) = ttm3(3) - ttm3r(3) * (1.0d0-mscale(jj))       
               ttm2i(1) = ttm2i(1) - ttm2ri(1)
               ttm2i(2) = ttm2i(2) - ttm2ri(2)
               ttm2i(3) = ttm2i(3) - ttm2ri(3)
               ttm3i(1) = ttm3i(1) - ttm3ri(1)
               ttm3i(2) = ttm3i(2) - ttm3ri(2)
               ttm3i(3) = ttm3i(3) - ttm3ri(3)
               if (ii .eq. jj) then
                  ttm2(1) = 0.5d0 * ttm2(1)
                  ttm2(2) = 0.5d0 * ttm2(2)
                  ttm2(3) = 0.5d0 * ttm2(3)
                  ttm3(1) = 0.5d0 * ttm3(1)
                  ttm3(2) = 0.5d0 * ttm3(2)
                  ttm3(3) = 0.5d0 * ttm3(3)
                  ttm2i(1) = 0.5d0 * ttm2i(1)
                  ttm2i(2) = 0.5d0 * ttm2i(2)
                  ttm2i(3) = 0.5d0 * ttm2i(3)
                  ttm3i(1) = 0.5d0 * ttm3i(1)
                  ttm3i(2) = 0.5d0 * ttm3i(2)
                  ttm3i(3) = 0.5d0 * ttm3i(3)
	       end if
c
c     now compute the virial components for this interaction
c
               do k = 1, 3
                  frczi(k) = 0.0d0
                  frcxi(k) = 0.0d0
                  frczj(k) = 0.0d0
                  frcxj(k) = 0.0d0
               end do
               xiz = x(zaxis(i)) - x(ii)
               yiz = y(zaxis(i)) - y(ii)
               ziz = z(zaxis(i)) - z(ii)
               xix = x(xaxis(i)) - x(ii)
               yix = y(xaxis(i)) - y(ii)
               zix = z(xaxis(i)) - z(ii)
               xjz = x(zaxis(j)) - x(jj)
               yjz = y(zaxis(j)) - y(jj)
               zjz = z(zaxis(j)) - z(jj)
               xjx = x(xaxis(j)) - x(jj)
               yjx = y(xaxis(j)) - y(jj)
               zjx = z(xaxis(j)) - z(jj)
               do k = 1, 3
                  tq(k) = ttm2(k) + ttm2i(k)
               end do
               call torque1 (i,tq,frczi,frcxi)
               do k = 1, 3
                  tq(k) = ttm3(k) + ttm3i(k)
               end do
               call torque1 (j,tq,frczj,frcxj)
               vxx = -xji*(ftm2(1)+ftm2i(1)) + xix*frcxi(1)
     &                  + xjz*frczj(1) + xjx*frcxj(1) + xiz*frczi(1)
               vyx = -yji*(ftm2(1)+ftm2i(1)) + yix*frcxi(1)
     &                  + yjz*frczj(1) + yjx*frcxj(1) + yiz*frczi(1)
               vzx = -zji*(ftm2(1)+ftm2i(1)) + zix*frcxi(1)
     &                  + zjz*frczj(1) + zjx*frcxj(1) + ziz*frczi(1)
               vxy = -xji*(ftm2(2)+ftm2i(2)) + xix*frcxi(2)
     &                  + xjz*frczj(2) + xjx*frcxj(2) + xiz*frczi(2)
               vyy = -yji*(ftm2(2)+ftm2i(2)) + yix*frcxi(2)
     &                  + yjz*frczj(2) + yjx*frcxj(2) + yiz*frczi(2)
               vzy = -zji*(ftm2(2)+ftm2i(2)) + zix*frcxi(2)
     &                  + zjz*frczj(2) + zjx*frcxj(2) + ziz*frczi(2)
               vxz = -xji*(ftm2(3)+ftm2i(3)) + xix*frcxi(3)
     &                  + xjz*frczj(3) + xjx*frcxj(3) + xiz*frczi(3)
               vyz = -yji*(ftm2(3)+ftm2i(3)) + yix*frcxi(3)
     &                  + yjz*frczj(3) + yjx*frcxj(3) + yiz*frczi(3)
               vzz = -zji*(ftm2(3)+ftm2i(3)) + zix*frcxi(3)
     &                  + zjz*frczj(3) + zjx*frcxj(3) + ziz*frczi(3)
               vir(1,1) = vir(1,1) + f*vxx
               vir(2,1) = vir(2,1) + 0.5d0*f*(vyx+vxy)
               vir(3,1) = vir(3,1) + 0.5d0*f*(vzx+vxz)
               vir(1,2) = vir(1,2) + 0.5d0*f*(vxy+vyx)
               vir(2,2) = vir(2,2) + f*vyy
               vir(3,2) = vir(3,2) + 0.5d0*f*(vzy+vyz)
               vir(1,3) = vir(1,3) + 0.5d0*f*(vxz+vzx)
               vir(2,3) = vir(2,3) + 0.5d0*f*(vyz+vzy)
               vir(3,3) = vir(3,3) + f*vzz
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
            end do
         end do
      end do
      end if
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
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine torque  --  convert multipole torque to force  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "torque" takes the torque values on sites defined by local
c     coordinate frames and distributes them to convert to forces
c     on the original sites and sites specifying the local frames
c
c     literature reference:
c
c     A. J. Stone and M. Alderton, "Distributed Multipole Analysis -
c     Methods and Applications", Molecular Physics, 56, 1047-1064 (1985)
c
c
      subroutine torque (trq,derivs)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      integer i,j,ia,ib,ic
      real*8 usiz,vsiz,wsiz
      real*8 upsiz,vpsiz
      real*8 dotdu,dotdv
      real*8 dphidu,dphidv,dphidw
      real*8 c,s,uvdis,vudis,du,dv
      real*8 u(3),v(3),w(3)
      real*8 up(3),vp(3),diff(3)
      real*8 trq(3,maxatm)
      real*8 derivs(3,maxatm)
c
c
c     coordinate frame motion described by rotation about u, v and w
c
      do i = 1, npole
         ia = zaxis(i)
         ib = ipole(i)
         ic = xaxis(i)
         u(1) = x(ia) - x(ib)
         u(2) = y(ia) - y(ib)
         u(3) = z(ia) - z(ib)
         usiz = sqrt(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))
         v(1) = x(ic) - x(ib)
         v(2) = y(ic) - y(ib)
         v(3) = z(ic) - z(ib)
         vsiz = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
         w(1) = u(2)*v(3) - u(3)*v(2)
         w(2) = u(3)*v(1) - u(1)*v(3)
         w(3) = u(1)*v(2) - u(2)*v(1)
         wsiz = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3))
         dotdu = 0.0d0
         dotdv = 0.0d0
         do j = 1, 3
            u(j) = u(j) / usiz
            v(j) = v(j) / vsiz
            w(j) = w(j) / wsiz
            diff(j) = v(j) - u(j)
            dotdu = dotdu + u(j)*diff(j)
            dotdv = dotdv + v(j)*diff(j)
         end do
c
c     get perpendiculars to u,v to get direction of motion of
c     u or v due to rotation about the cross product vector w
c
         upsiz = 0.0d0
         vpsiz = 0.0d0
         do j = 1, 3
            up(j) = diff(j) - dotdu*u(j)
            vp(j) = diff(j) - dotdv*v(j)
            upsiz = upsiz + up(j)*up(j)
            vpsiz = vpsiz + vp(j)*vp(j)
         end do
         upsiz = sqrt(upsiz)
         vpsiz = sqrt(vpsiz)
         do j = 1, 3
            up(j) = up(j) / upsiz
            vp(j) = vp(j) / vpsiz
         end do
c
c     negative of dot product of torque with unit vectors along u, v
c     and w give result of infinitesmal rotation along these vectors
c     i.e. dphi/dtheta = dot product, where dphi is work, dtheta is
c     angle; then dphi/dtheta is torque and the dot product is torque
c     component along unit vector
c
         dphidu = -trq(1,ib)*u(1) - trq(2,ib)*u(2) - trq(3,ib)*u(3)
         dphidv = -trq(1,ib)*v(1) - trq(2,ib)*v(2) - trq(3,ib)*v(3)
         dphidw = -trq(1,ib)*w(1) - trq(2,ib)*w(2) - trq(3,ib)*w(3)
c
c     get the projected distances between the vectors
c
         c = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
         s = sqrt(1.0d0 - c*c)
         uvdis = usiz * s
         vudis = vsiz * s
c
c     force distribution for the bisector local coordinate method
c
         if (polaxe(i) .eq. 'Bisector') then
            do j = 1, 3
               du = -w(j)*dphidv/uvdis + up(j)*dphidw/(2.0d0*usiz)
               dv = w(j)*dphidu/vudis + vp(j)*dphidw/(2.0d0*vsiz)
               derivs(j,ia) = derivs(j,ia) + du
               derivs(j,ic) = derivs(j,ic) + dv
               derivs(j,ib) = derivs(j,ib) - dv - du
            end do
c
c     force distribution for the Z-then-X local coordinate method
c
         else if (polaxe(i) .eq. 'Z-then-X') then
            do j = 1, 3
               du = -w(j)*dphidv/uvdis + up(j)*dphidw/usiz
               dv = w(j)*dphidu/vudis
               derivs(j,ia) = derivs(j,ia) + du
               derivs(j,ic) = derivs(j,ic) + dv
               derivs(j,ib) = derivs(j,ib) - dv - du
            end do
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine torque1  --  convert torque to force; 1-site  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "torque1" takes the torque value on a site defined by a local
c     coordinate frame and distributes it to convert to forces on
c     the original site and sites specifying the local frame
c
c     literature reference:
c
c     A. J. Stone and M. Alderton, "Distributed Multipole Analysis -
c     Methods and Applications", Molecular Physics, 56, 1047-1064 (1985)
c
c
      subroutine torque1 (i,trq,frcz,frcx)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      integer i,j,ia,ib,ic
      real*8 usiz,vsiz,wsiz
      real*8 upsiz,vpsiz
      real*8 dotdu,dotdv
      real*8 dphidu,dphidv,dphidw
      real*8 c,s,uvdis,vudis
      real*8 u(3),v(3),w(3)
      real*8 up(3),vp(3)
      real*8 trq(3),diff(3)
      real*8 frcz(3),frcx(3)
c
c
c     coordinate frame motion described by rotation about u, v and w
c
      ia = zaxis(i)
      ib = ipole(i)
      ic = xaxis(i)
      u(1) = x(ia) - x(ib)
      u(2) = y(ia) - y(ib)
      u(3) = z(ia) - z(ib)
      usiz = sqrt(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))
      v(1) = x(ic) - x(ib)
      v(2) = y(ic) - y(ib)
      v(3) = z(ic) - z(ib)
      vsiz = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
      w(1) = u(2)*v(3) - u(3)*v(2)
      w(2) = u(3)*v(1) - u(1)*v(3)
      w(3) = u(1)*v(2) - u(2)*v(1)
      wsiz = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3))
      dotdu = 0.0d0
      dotdv = 0.0d0
      do j = 1, 3
         u(j) = u(j) / usiz
         v(j) = v(j) / vsiz
         w(j) = w(j) / wsiz
         diff(j) = v(j) - u(j)
         dotdu = dotdu + u(j)*diff(j)
         dotdv = dotdv + v(j)*diff(j)
      end do
c
c     get perpendiculars to u,v to get direction of motion of
c     u or v due to rotation about the cross product vector w
c
      upsiz = 0.0d0
      vpsiz = 0.0d0
      do j = 1, 3
         up(j) = diff(j) - dotdu*u(j)
         vp(j) = diff(j) - dotdv*v(j)
         upsiz = upsiz + up(j)*up(j)
         vpsiz = vpsiz + vp(j)*vp(j)
      end do
      upsiz = sqrt(upsiz)
      vpsiz = sqrt(vpsiz)
      do j = 1, 3
         up(j) = up(j) / upsiz
         vp(j) = vp(j) / vpsiz
      end do
c
c     negative of dot product of torque with unit vectors along u, v
c     and w give result of infinitesmal rotation along these vectors
c     i.e. dphi/dtheta = dot product, where dphi is work, dtheta is
c     angle; then dphi/dtheta is torque and the dot product is torque
c     component along unit vector
c
      dphidu = -trq(1)*u(1) - trq(2)*u(2) - trq(3)*u(3)
      dphidv = -trq(1)*v(1) - trq(2)*v(2) - trq(3)*v(3)
      dphidw = -trq(1)*w(1) - trq(2)*w(2) - trq(3)*w(3)
c
c     get the projected distances between the vectors
c
      c = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
      s = sqrt(1.0d0 - c*c)
      uvdis = usiz * s
      vudis = vsiz * s
c
c     force distribution for the bisector local coordinate method
c
      if (polaxe(i) .eq. 'Bisector') then
         do j = 1, 3
            frcz(j) = -w(j)*dphidv/uvdis + up(j)*dphidw/(2.0d0*usiz)
            frcx(j) = w(j)*dphidu/vudis + vp(j)*dphidw/(2.0d0*vsiz)
         end do
c
c     force distribution for the Z-then-X local coordinate method
c
      else
         do j = 1, 3
            frcz(j) = -w(j)*dphidv/uvdis + up(j)*dphidw/usiz
            frcx(j) = w(j)*dphidu/vudis
         end do
      end if
      return
      end
