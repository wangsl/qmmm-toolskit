 
C     $Id: empole.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine empole  --  multipole & polarization energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "empole" calculates the electrostatic energy due to
c     atomic multipole interactions and dipole polarizability
c
c
      subroutine empole
      implicit none
      include 'cutoff.i'
c
c
c     choose pairwise double loop or Ewald summation version
c
      if (use_ewald) then
         call empole0b
      else
         call empole0a
      end if
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine empole0a  --  double loop multipole energy  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "empole0a" calculates the electrostatic energy due to
c     atomic multipole interactions and dipole polarizability
c     using a pairwise double loop
c
c
      subroutine empole0a
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'math.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k
      integer ii,kk
      integer ix,iz
      integer kx,kz
      real*8 e,ei,f
      real*8 fikm,fikp
      real*8 r,r2,xr,yr,zr
      real*8 fgrp,damp
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qix,qiy,qiz
      real*8 qkx,qky,qkz
      real*8 scale3,scale5
      real*8 scale7
      real*8 sc(10),sci(8)
      real*8 gl(0:4),gli(3)
      real*8 mscale(maxatm)
      real*8 pscale(maxatm)
      logical proceed,iuse,kuse
c
c
c     zero out the multipole and polarization energies
c
      em = 0.0d0
      ep = 0.0d0
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
c     set conversion factor, cutoff and scaling coefficients
c
      if (npole .eq. 0)  return
      f = electric / dielec
      call switch ('MPOLE')
c
c     calculate the multipole interaction energy term
c
      do i = 1, npole-1
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
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
         iuse = (use(ii) .or. use(iz) .or. use(ix))
         do j = i+1, npole
            mscale(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
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
c
c     decide whether to compute the current interaction
c
         do k = i+1, npole
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            kuse = (use(kk) .or. use(kz) .or. use(kx))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (iuse .or. kuse)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               call image (xr,yr,zr,0)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  ukx = uind(1,k)
                  uky = uind(2,k)
                  ukz = uind(3,k)
c
c     construct some intermediate quadrupole values
c
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
                  sc(2) = dix*dkx + diy*dky + diz*dkz
                  sc(3) = dix*xr + diy*yr + diz*zr
                  sc(4) = dkx*xr + dky*yr + dkz*zr
                  sc(5) = qix*xr + qiy*yr + qiz*zr
                  sc(6) = qkx*xr + qky*yr + qkz*zr
                  sc(7) = qix*dkx + qiy*dky + qiz*dkz
                  sc(8) = qkx*dix + qky*diy + qkz*diz
                  sc(9) = qix*qkx + qiy*qky + qiz*qkz
                  sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                        + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
                  sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                        + diy*uky + uiz*dkz + diz*ukz
                  sci(3) = uix*xr + uiy*yr + uiz*zr
                  sci(4) = ukx*xr + uky*yr + ukz*zr
                  sci(7) = qix*ukx + qiy*uky + qiz*ukz
                  sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
                  gl(0) = ci*ck
                  gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
                  gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                       + 2.0d0*(sc(7)-sc(8)+sc(10))
                  gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
                  gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
                  gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
                  gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                        - sc(3)*sci(4)
                  gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdamp(i) * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        scale3 = 1.0d0 - exp(damp)
                        scale5 = 1.0d0 - (1.0d0-damp)*exp(damp)
                        scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                          *exp(damp)
                     end if
                  end if
                  e = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                   + gl(3)*rr7 + gl(4)*rr9
                  ei = gli(1)*rr3*scale3 + gli(2)*rr5*scale5
     &                    + gli(3)*rr7*scale7
c
c     make the adjustment for scaled interactions
c
                  fikm = f * mscale(kk)
                  fikp = f * pscale(kk)
                  e = fikm * e
                  ei = 0.5d0 * fikp * ei
c
c     scale the interaction based on its group membership;
c     polarization cannot be group scaled as it is not pairwise
c
                  if (use_group) then
                     e = e * fgrp
c                    ei = ei * fgrp
                  end if
c
c     increment the overall multipole and polarization energies
c
                  em = em + e
                  ep = ep + ei
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
      do i = 1, npole
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
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
         iuse = (use(ii) .or. use(iz) .or. use(ix))
         do j = i, npole
            mscale(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
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
c
c     decide whether to compute the current interaction
c
         do k = i, npole
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            kuse = (use(kk) .or. use(kz) .or. use(kx))
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (iuse .or. kuse)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call image (xr,yr,zr,j)
                  r2 = xr*xr + yr* yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     ck = rpole(1,k)
                     dkx = rpole(2,k)
                     dky = rpole(3,k)
                     dkz = rpole(4,k)
                     qkxx = rpole(5,k)
                     qkxy = rpole(6,k)
                     qkxz = rpole(7,k)
                     qkyy = rpole(9,k)
                     qkyz = rpole(10,k)
                     qkzz = rpole(13,k)
                     ukx = uind(1,k)
                     uky = uind(2,k)
                     ukz = uind(3,k)
c
c     construct some intermediate quadrupole values
c
                     qix = qixx*xr + qixy*yr + qixz*zr
                     qiy = qixy*xr + qiyy*yr + qiyz*zr
                     qiz = qixz*xr + qiyz*yr + qizz*zr
                     qkx = qkxx*xr + qkxy*yr + qkxz*zr
                     qky = qkxy*xr + qkyy*yr + qkyz*zr
                     qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
                     sc(2) = dix*dkx + diy*dky + diz*dkz
                     sc(3) = dix*xr + diy*yr + diz*zr
                     sc(4) = dkx*xr + dky*yr + dkz*zr
                     sc(5) = qix*xr + qiy*yr + qiz*zr
                     sc(6) = qkx*xr + qky*yr + qkz*zr
                     sc(7) = qix*dkx + qiy*dky + qiz*dkz
                     sc(8) = qkx*dix + qky*diy + qkz*diz
                     sc(9) = qix*qkx + qiy*qky + qiz*qkz
                     sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                           + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
                     sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                           + diy*uky + uiz*dkz + diz*ukz
                     sci(3) = uix*xr + uiy*yr + uiz*zr
                     sci(4) = ukx*xr + uky*yr + ukz*zr
                     sci(7) = qix*ukx + qiy*uky + qiz*ukz
                     sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
                     gl(0) = ci*ck
                     gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
                     gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                          + 2.0d0*(sc(7)-sc(8)+sc(10))
                     gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
                     gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
                     gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
                     gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                           - sc(3)*sci(4)
                     gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
                     rr1 = 1.0d0 / r
                     rr3 = rr1 / r2
                     rr5 = 3.0d0 * rr3 / r2
                     rr7 = 5.0d0 * rr5 / r2
                     rr9 = 7.0d0 * rr7 / r2
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     scale7 = 1.0d0
                     damp = pdamp(i) * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           scale3 = 1.0d0 - exp(damp)
                           scale5 = 1.0d0 - (1.0d0-damp)*exp(damp)
                           scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                             *exp(damp)
                        end if
                     end if
                     e = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                      + gl(3)*rr7 + gl(4)*rr9
                     ei = gli(1)*rr3*scale3 + gli(2)*rr5*scale5
     &                       + gli(3)*rr7*scale7
                     fikm = f
                     fikp = f
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           fikm = fikm * mscale(kk)
                           fikp = fikp * pscale(kk)
                        end if
                     end if
                     e = fikm * e
                     ei = 0.5d0 * fikp * ei
c
c     scale the interaction based on its group membership;
c     polarization cannot be group scaled as it is not pairwise
c
                     if (use_group) then
                        e = e * fgrp
c                       ei = ei * fgrp
                     end if
c
c     increment the overall multipole and polarization energies
c
                     if (ii .eq. kk) then
                        e = 0.5d0 * e
                        ei = 0.5d0 * ei
                     end if
                     em = em + e
                     ep = ep + ei
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
c     ##  subroutine empole0b  --  Ewald summation multipole energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "empole0b" calculates the electrostatic energy due to
c     atomic multipole interactions and dipole polarizability
c     using a regular Ewald summation
c
c
      subroutine empole0b
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'energi.i'
      include 'ewald.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'units.i'
      integer i,ii
      real*8 e,ei,f
      real*8 term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the multipole and polarization energies
c
      em = 0.0d0
      ep = 0.0d0
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
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     compute the reciprocal space part of the Ewald summation
c
      call erecip
c
c     compute the real space portion of the Ewald summation
c
      call ereal
c
c     compute the self-energy portion of the Ewald summation
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
c     compute the cell dipole boundary correction term
c
      if (.not. tinfoil) then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xu = 0.0d0
         yu = 0.0d0
         zu = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            uix = uind(1,i)
            uiy = uind(2,i)
            uiz = uind(3,i)
            xd = xd + dix + rpole(1,i)*x(ii)
            yd = yd + diy + rpole(1,i)*y(ii)
            zd = zd + diz + rpole(1,i)*z(ii)
            xu = xu + uix
            yu = yu + uiy
            zu = zu + uiz
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         em = em + term*(xd*xd+yd*yd+zd*zd)
         ep = ep + term*(xd*xu+yd*yu+zd*zu)
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine erecip  --  multipole Ewald reciprocal energy  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "erecip" evaluates the reciprocal space portion of the regular
c     Ewald summation energy due to atomic multipole interactions
c     and dipole polarizability
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine erecip
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'energi.i'
      include 'ewald.i'
      include 'ewreg.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'units.i'
      integer i,j,k,l,ii
      integer jmin,jmax
      integer kmin,kmax
      integer lmin,lmax
      real*8 e,ei,f,term,cut
      real*8 eterm,expterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3,hsq
      real*8 t1,t2,t3,t4
      real*8 ck,dk,qk,uk
      real*8 q1,q2,q3
      real*8 ckr,skr
      real*8 cm(maxatm)
      real*8 dm(3,maxatm)
      real*8 qm(9,maxatm)
      real*8 um(3,maxatm)
      real*8 cjk(maxatm)
      real*8 sjk(maxatm)
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
                  do i = 1, npole
                     ckr = cjk(i)*elc(i,l) - sjk(i)*els(i,l)
                     skr = sjk(i)*elc(i,l) + cjk(i)*els(i,l)
                     ck = cm(i)
                     dk = h1*dm(1,i) + h2*dm(2,i) + h3*dm(3,i)
                     q1 = h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i)
                     q2 = h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i)
                     q3 = h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i)
                     qk = h1*q1 + h2*q2 + h3*q3
                     uk = h1*um(1,i) + h2*um(2,i) + h3*um(3,i)
                     t1 = t1 + (ck-qk)*skr + dk*ckr
                     t2 = t2 + (ck-qk)*ckr - dk*skr
                     t3 = t3 + uk*ckr
                     t4 = t4 - uk*skr
                  end do
                  expterm = eterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0 
                  end if
                  e = expterm * (t1*t1+t2*t2)
                  ei = expterm * (t1*t3+t2*t4)
                  em = em + e
                  ep = ep + ei
               end if
            end do
            lmin = -lmax
         end do
         kmin = -kmax
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ereal  --  multipole Ewald real space energy  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ereal" evaluates the real space portion of the regular Ewald
c     summation energy due to atomic multipole interactions and dipole
c     polarizability
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine ereal
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'ewald.i'
      include 'math.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'units.i'
      integer i,j,k,m
      integer ii,kk
      real*8 e,ei,f,erfc
      real*8 r,r2,xr,yr,zr
      real*8 bfac,exp2a
      real*8 efix,eifix
      real*8 ralpha,damp
      real*8 alsq2,alsq2n
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qix,qiy,qiz
      real*8 qkx,qky,qkz
      real*8 scale3,scale5
      real*8 scale7
      real*8 sc(10),sci(8)
      real*8 gl(0:4),gli(3)
      real*8 bn(0:4)
      real*8 mscale(maxatm)
      real*8 pscale(maxatm)
      external erfc
c
c
c     set conversion factor, cutoff and scaling coefficients
c
      if (npole .eq. 0)  return
      f = electric / dielec
      call switch ('EWALD')
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole-1
         ii = ipole(i)
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
         do j = i+1, npole
            mscale(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
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
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr,0)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. cut2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
               ukx = uind(1,k)
               uky = uind(2,k)
               ukz = uind(3,k)
c
c     calculate the error function damping terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)
     &            alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do m = 1, 4
                  bfac = dble(m+m-1)
                  alsq2n = alsq2 * alsq2n
                  bn(m) = (bfac*bn(m-1)+alsq2n*exp2a) / r2
               end do
c
c     construct some intermediate quadrupole values
c
               qix = qixx*xr + qixy*yr + qixz*zr
               qiy = qixy*xr + qiyy*yr + qiyz*zr
               qiz = qixz*xr + qiyz*yr + qizz*zr
               qkx = qkxx*xr + qkxy*yr + qkxz*zr
               qky = qkxy*xr + qkyy*yr + qkyz*zr
               qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
               sc(2) = dix*dkx + diy*dky + diz*dkz
               sc(3) = dix*xr + diy*yr + diz*zr
               sc(4) = dkx*xr + dky*yr + dkz*zr
               sc(5) = qix*xr + qiy*yr + qiz*zr
               sc(6) = qkx*xr + qky*yr + qkz*zr
               sc(7) = qix*dkx + qiy*dky + qiz*dkz
               sc(8) = qkx*dix + qky*diy + qkz*diz
               sc(9) = qix*qkx + qiy*qky + qiz*qkz
               sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                     + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
               sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                     + diy*uky + uiz*dkz + diz*ukz
               sci(3) = uix*xr + uiy*yr + uiz*zr
               sci(4) = ukx*xr + uky*yr + ukz*zr
               sci(7) = qix*ukx + qiy*uky + qiz*ukz
               sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
               gl(0) = ci*ck
               gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                    + 2.0d0*(sc(7)-sc(8)+sc(10))
               gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
               gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
               gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
               gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                     - sc(3)*sci(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
               e = gl(0)*bn(0) + gl(1)*bn(1) + gl(2)*bn(2)
     &                + gl(3)*bn(3) + gl(4)*bn(4)
               ei = gli(1)*bn(1) + gli(2)*bn(2) + gli(3)*bn(3)
c
c     full energies needeed for scaling and damping correction
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               damp = pdamp(i) * pdamp(k)
               scale3 = pscale(kk)
               scale5 = pscale(kk)
               scale7 = pscale(kk)
               damp = pdamp(i) * pdamp(k)
               if (damp .ne. 0.0d0) then
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     scale3 = scale3 * (1.0d0-exp(damp))
                     scale5 = scale5 * (1.0d0-(1.0d0-damp)*exp(damp))
                     scale7 = scale7 * (1.0d0-(1.0d0-damp+0.6d0*damp**2)
     &                                               *exp(damp))
                  end if
               end if
               efix = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                   + gl(3)*rr7 + gl(4)*rr9
               eifix = gli(1)*rr3*(1.0d0-scale3)
     &                    + gli(2)*rr5*(1.0d0-scale5)
     &                    + gli(3)*rr7*(1.0d0-scale7)
c
c     correction for scaled interactions and polarization damping
c
               e = e - efix*(1.0d0-mscale(kk))
               ei = ei - eifix
c
c     increment the overall multipole and polarization energies
c
               e = f * e
               ei = 0.5d0 * f * ei
               em = em + e
               ep = ep + ei
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
      do i = 1, npole
         ii = ipole(i)
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
         do j = i, npole
            mscale(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
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
         do k = i, npole
            kk = ipole(k)
            do j = 1, ncell
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               call image (xr,yr,zr,j)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. cut2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  ukx = uind(1,k)
                  uky = uind(2,k)
                  ukz = uind(3,k)
c
c     calculate the error function damping terms
c
                  ralpha = aewald * r
                  bn(0) = erfc(ralpha) / r
                  alsq2 = 2.0d0 * aewald**2
                  alsq2n = 0.0d0
                  if (aewald .gt. 0.0d0)
     &               alsq2n = 1.0d0 / (sqrtpi*aewald)
                  exp2a = exp(-ralpha**2)
                  do m = 1, 4
                     bfac = dble(m+m-1)
                     alsq2n = alsq2 * alsq2n
                     bn(m) = (bfac*bn(m-1)+alsq2n*exp2a) / r2
                  end do
c
c     construct some intermediate quadrupole values
c
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
                  sc(2) = dix*dkx + diy*dky + diz*dkz
                  sc(3) = dix*xr + diy*yr + diz*zr
                  sc(4) = dkx*xr + dky*yr + dkz*zr
                  sc(5) = qix*xr + qiy*yr + qiz*zr
                  sc(6) = qkx*xr + qky*yr + qkz*zr
                  sc(7) = qix*dkx + qiy*dky + qiz*dkz
                  sc(8) = qkx*dix + qky*diy + qkz*diz
                  sc(9) = qix*qkx + qiy*qky + qiz*qkz
                  sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                        + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
                  sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                        + diy*uky + uiz*dkz + diz*ukz
                  sci(3) = uix*xr + uiy*yr + uiz*zr
                  sci(4) = ukx*xr + uky*yr + ukz*zr
                  sci(7) = qix*ukx + qiy*uky + qiz*ukz
                  sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
                  gl(0) = ci*ck
                  gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
                  gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                       + 2.0d0*(sc(7)-sc(8)+sc(10))
                  gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
                  gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
                  gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
                  gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                        - sc(3)*sci(4)
                  gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
                  e = gl(0)*bn(0) + gl(1)*bn(1) + gl(2)*bn(2)
     &                   + gl(3)*bn(3) + gl(4)*bn(4)
                  ei = gli(1)*bn(1) + gli(2)*bn(2) + gli(3)*bn(3)
c
c     full energies needeed for scaling and damping correction
c
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
                  damp = pdamp(i) * pdamp(k)
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdamp(i) * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        scale3 = 1.0d0 - exp(damp)
                        scale5 = 1.0d0 - (1.0d0-damp)*exp(damp)
                        scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                          *exp(damp)
                     end if
                  end if
                  if (use_polymer .and. r2.le.polycut2) then
                     scale3 = scale3 * pscale(kk)
                     scale5 = scale5 * pscale(kk)
                     scale7 = scale7 * pscale(kk)
                  end if
                  efix = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                      + gl(3)*rr7 + gl(4)*rr9
                  eifix = gli(1)*rr3*(1.0d0-scale3)
     &                       + gli(2)*rr5*(1.0d0-scale5)
     &                       + gli(3)*rr7*(1.0d0-scale7)
c
c     correction for scaled interactions and polarization damping
c
                  if (use_polymer .and. r2.le.polycut2)
     &               e = e - efix*(1.0d0-mscale(kk))
                  ei = ei - eifix
c
c     increment the overall multipole and polarization energies
c
                  e = f * e
                  ei = 0.5d0 * f * ei
                  if (ii .eq. kk) then
                     e = 0.5d0 * e
                     ei = 0.5d0 * ei
                  end if
                  em = em + e
                  ep = ep + ei
               end if
            end do
         end do
      end do
      return
      end
