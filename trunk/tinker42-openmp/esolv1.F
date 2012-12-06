 
C     $Id: esolv1.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine esolv1  --  solvation energy and derivatives  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "esolv1" calculates the continuum solvation energy and
c     first derivatives with respect to Cartesian coordinates
c     using either the Eisenberg-McLachlan ASP, Ooi-Scheraga
c     SASA or various GB/SA solvation models
c
c
      subroutine esolv1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'energi.i'
      include 'math.i'
      include 'potent.i'
      include 'solute.i'
      include 'warp.i'
      integer i
      real*8 e,ai,ri,rb
      real*8 term,probe
      real*8 aes(maxatm)
c
c
c     zero out the continuum solvation energy and derivatives
c
      es = 0.0d0
      do i = 1, n
         drb(i) = 0.0d0
         des(1,i) = 0.0d0
         des(2,i) = 0.0d0
         des(3,i) = 0.0d0
      end do
c
c     set a value for the solvent molecule probe radius
c
      probe = 1.4d0
c
c     get the atomic Born radii values for the GB/SA models
c
      if (use_gbsa)  call born
c
c     compute the nonpolar solvation via ACE approximation
c
      if (use_gbsa .and. solvtyp.ne.'ONION') then
         term = 4.0d0 * pi
         do i = 1, n
            ai = asolv(i)
            ri = rsolv(i)
            rb = rborn(i)
            if (rb .ne. 0.0d0) then
               e = ai * term * (ri+probe)**2 * (ri/rb)**6
               es = es + e
               drb(i) = drb(i) - 6.0d0*e/rb
            end if
         end do
c
c     compute the nonpolar solvation via exact surface area
c
      else
         call surface (es,aes,des,rsolv,asolv,probe)
      end if
c
c     compute the generalized Born polarization solvation
c
      if (use_gbsa) then
         if (use_smooth) then
            call egbsa1b
         else
            call egbsa1a
         end if
      end if
c
c     increment derivatives due to Born radii chain rule terms
c
      if (use_gbsa)  call born1
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine egbsa1a  --  GB/SA energy and derivatives  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "egbsa1a" calculates the generalized Born energy and first
c     derivatives of the GB/SA solvation models
c
c     note application of distance cutoff scaling directly to
c     the Born radii chain rule term "derb" is an approximation
c
c
      subroutine egbsa1a
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'charge.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'molcul.i'
      include 'shunt.i'
      include 'solute.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i,k,ii,kk
      real*8 e,de,fgrp
      real*8 f,fi,fik
      real*8 fgb,fgb2
      real*8 rb2,rm2
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,r3,r4
      real*8 r5,r6,r7
      real*8 dwater,rbi,rbk
      real*8 dedx,dedy,dedz
      real*8 derb,drbi,drbk
      real*8 expterm,shift
      real*8 taper,dtaper
      real*8 trans,dtrans
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      logical proceed,iuse
c
c
c     set the solvent dielectric and energy conversion factor
c
      if (nion .eq. 0)  return
      dwater = 78.3d0
      f = -electric * (1.0d0 - 1.0d0/dwater)
c
c     set cutoff distances and switching function coefficients
c
      call switch ('CHARGE')
c
c     calculate GB/SA electrostatic polarization energy term
c
      do ii = 1, nion
         i = iion(ii)
         iuse = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         rbi = rborn(i)
         do kk = ii, nion
            k = iion(kk)
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
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rbk = rborn(k)
                  fik = fi * pchg(kk)
                  rb2 = rbi * rbk
                  expterm = exp(-0.25d0*r2/rb2)
                  fgb2 = r2 + rb2*expterm
                  fgb = sqrt(fgb2)
                  e = fik / fgb
                  de = -e * (r-0.25d0*r*expterm) / fgb2
                  derb = -e * expterm*(0.5d0+0.125d0*r2/rb2) / fgb2
c
c     use energy switching if near the cutoff distance
c
                  rm2 = (0.5d0 * (off+cut))**2
                  shift = fik / sqrt(rm2 + rb2*exp(-0.25d0*rm2/rb2))
                  e = e - shift
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     r6 = r3 * r3
                     r7 = r3 * r4
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                           + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                               + f3*r3 + f2*r2 + f1*r + f0)
                     dtrans = fik * (7.0d0*f7*r6 + 6.0d0*f6*r5
     &                               + 5.0d0*f5*r4 + 4.0d0*f4*r3
     &                             + 3.0d0*f3*r2 + 2.0d0*f2*r + f1)
                     derb = derb * taper
                     de = e*dtaper + de*taper + dtrans
                     e = e*taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     derb = derb * fgrp
                  end if
c
c     increment the overall energy and derivative expressions
c
                  if (i .eq. k) then
                     e = 0.5d0 * e
                     es = es + e
                     drbi = derb * rbk
                     drb(i) = drb(i) + drbi
                  else
                     es = es + e
                     de = de / r
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
                     des(1,i) = des(1,i) + dedx
                     des(2,i) = des(2,i) + dedy
                     des(3,i) = des(3,i) + dedz
                     des(1,k) = des(1,k) - dedx
                     des(2,k) = des(2,k) - dedy
                     des(3,k) = des(3,k) - dedz
                     drbi = derb * rbk
                     drbk = derb * rbi
                     drb(i) = drb(i) + drbi
                     drb(k) = drb(k) + drbk
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine egbsa1b  --  GB/SA energy/derivs for smoothing  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "egbsa1b" calculates the generalized Born energy and first
c     derivatives of the GB/SA solvation models for use with potential
c     smoothing methods
c
c
      subroutine egbsa1b
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'charge.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'math.i'
      include 'molcul.i'
      include 'solute.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      include 'warp.i'
      integer i,k,ii,kk
      real*8 e,de,fgrp
      real*8 f,fi,fik
      real*8 fgb,fgb2
      real*8 rb2,width
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,sterm
      real*8 expterm
      real*8 dwater,rbi,rbk
      real*8 dedx,dedy,dedz
      real*8 derb,drbi,drbk
      real*8 erf,erfterm,term
      real*8 wterm,rterm,bterm
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      logical proceed,iuse
      external erf
c
c
c     set the solvent dielectric and energy conversion factor
c
      if (nion .eq. 0)  return
      dwater = 78.3d0
      f = -electric * (1.0d0 - 1.0d0/dwater)
c
c     set the extent of smoothing to be performed
c
      sterm = 0.5d0 / sqrt(diffc)
c
c     calculate GB/SA electrostatic polarization energy term
c
      do ii = 1, nion
         i = iion(ii)
         iuse = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         rbi = rborn(i)
         do kk = ii, nion
            k = iion(kk)
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
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               r = sqrt(r2)
               rbk = rborn(k)
               fik = fi * pchg(kk)
               rb2 = rbi * rbk
               expterm = exp(-0.25d0*r2/rb2)
               fgb2 = r2 + rb2*expterm
               fgb = sqrt(fgb2)
               e = fik / fgb
               de = -e * (r-0.25d0*r*expterm) / fgb2
               derb = -e * expterm*(0.5d0+0.125d0*r2/rb2) / fgb2
c
c     use a smoothable GB/SA analogous to the Coulomb solution
c
               if (deform .gt. 0.0d0) then
                  wterm = exp(-0.006d0*rb2/deform)
                  width = sterm / sqrt(deform+0.15d0*rb2*wterm)
                  erfterm = erf(width*fgb)
                  term = width * exp(-(width*fgb)**2) / sqrtpi
                  rterm = term * (2.0d0*r-0.5d0*r*expterm)/fgb
                  bterm = term * ((expterm*(1.0d0+0.25d0*r2/rb2)/fgb)
     &                              - (fgb*(width/sterm)**2) * wterm
     &                                 * (0.15d0-0.0009d0*rb2/deform))
                  derb = derb*erfterm + e*bterm
                  de = de*erfterm + e*rterm
                  e = e * erfterm
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  de = de * fgrp
                  derb = derb * fgrp
               end if
c
c     increment the overall energy and derivative expressions
c
               if (i .eq. k) then
                  e = 0.5d0 * e
                  es = es + e
                  drbi = derb * rbk
                  drb(i) = drb(i) + drbi
               else
                  es = es + e
                  de = de / r
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  des(1,i) = des(1,i) + dedx
                  des(2,i) = des(2,i) + dedy
                  des(3,i) = des(3,i) + dedz
                  des(1,k) = des(1,k) - dedx
                  des(2,k) = des(2,k) - dedy
                  des(3,k) = des(3,k) - dedz
                  drbi = derb * rbk
                  drbk = derb * rbi
                  drb(i) = drb(i) + drbi
                  drb(k) = drb(k) + drbk
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
c     ##  subroutine born1  --  Born radii chain rule derivatives  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "born1" computes derivatives of the Born radii with respect
c     to atomic coordinates and increments total energy derivatives
c     and virial components for potentials involving Born radii
c
c
      subroutine born1
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'deriv.i'
      include 'math.i'
      include 'solute.i'
      include 'units.i'
      include 'virial.i'
      integer i,j,k,it,kt
      integer skip(maxatm)
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 de,de1,de2
      real*8 r,r2,r3,r4,r6
      real*8 p5inv,pip5
      real*8 gpi,vk,ratio
      real*8 ccf,cosq,dccf
      real*8 sinq,term,theta
      real*8 rb2,ri,rk,sk,sk2
      real*8 lik,lik2,lik3
      real*8 uik,uik2,uik3
      real*8 dlik,duik
      real*8 t1,t2,t3
      real*8 rbi,rbi2,vi
      real*8 ws2,s2ik,uik4
      real*8 expterm,rusum
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
c
c
c     get Born radius chain rule components for the Still method
c
      if (solvtyp .eq. 'STILL') then
         p5inv = 1.0d0 / p5
         pip5 = pi * p5
         do i = 1, n
            skip(i) = 0
         end do
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            skip(i) = i
            do j = 1, n12(i)
               skip(i12(j,i)) = i
            end do
            do j = 1, n13(i)
               skip(i13(j,i)) = i
            end do
            gpi = 2.0d0 * rborn(i)**2 / electric
            do k = 1, n
               if (skip(k) .ne. i) then
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  vk = vsolv(k)
                  r2 = xr**2 + yr**2 + zr**2
                  r =  sqrt(r2)
                  r6 = r2 * r2 * r2
                  ratio = r2 / (rsolv(i)+rsolv(k))**2
                  if (ratio .gt. p5inv) then
                     ccf = 1.0d0
                     dccf = 0.0d0
                  else
                     theta = ratio * pip5
                     cosq = cos(theta)
                     term = 0.5d0 * (1.0d0-cosq)
                     ccf = term * term
                     sinq = sin(theta)
                     dccf = 2.0d0 * term * sinq * pip5 * ratio
                  end if
                  de = drb(i) * p4 * gpi * vk * (4.0d0*ccf-dccf) / r6
c
c     increment the overall continuum solvation derivatives
c
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  des(1,i) = des(1,i) + dedx
                  des(2,i) = des(2,i) + dedy
                  des(3,i) = des(3,i) + dedz
                  des(1,k) = des(1,k) - dedx
                  des(2,k) = des(2,k) - dedy
                  des(3,k) = des(3,k) - dedz
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
         end do
c
c     get Born radius chain rule components for the HCT method
c
      else if (solvtyp .eq. 'HCT') then
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ri = rsolv(i) + doffset
            rb2 = rborn(i) * rborn(i)
            do k = 1, n
               if (k .ne. i) then
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  rk = rsolv(k) + doffset
                  sk = rk * shct(k)
                  sk2 = sk * sk
                  r2 = xr**2 + yr**2 + zr**2
                  r = sqrt(r2)
                  if (ri .lt. r+sk) then
                     lik = 1.0d0 / max(ri,r-sk)
                     uik = 1.0d0 / (r+sk)
                     lik2 = lik * lik
                     uik2 = uik * uik
                     lik3 = lik * lik2
                     uik3 = uik * uik2
                     dlik = 1.0d0
                     if (ri .ge. r-sk)  dlik = 0.0d0
                     duik = 1.0d0
                     t1 = 0.5d0*lik2 + 0.25d0*sk2*lik3/r
     &                       - 0.25d0*(lik/r+lik3*r)
                     t2 = -0.5d0*uik2 - 0.25d0*sk2*uik3/r
     &                       + 0.25d0*(uik/r+uik3*r)
                     t3 = 0.125d0*(1.0d0+sk2/r2)*(lik2-uik2)
     &                       + 0.25d0*log(uik/lik)/r2
                     de = drb(i) * rb2 * (dlik*t1+duik*t2+t3) / r
c
c     increment the overall continuum solvation derivatives
c
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
                     des(1,i) = des(1,i) + dedx
                     des(2,i) = des(2,i) + dedy
                     des(3,i) = des(3,i) + dedz
                     des(1,k) = des(1,k) - dedx
                     des(2,k) = des(2,k) - dedy
                     des(3,k) = des(3,k) - dedz
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
               end if
            end do
         end do
c
c     get Born radius chain rule components for the ACE method
c
      else if (solvtyp .eq. 'ACE') then
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            it = class(i)
            vi = vsolv(i)
            rbi = rborn(i)
            rbi2 = rbi * rbi
            do k = 1, n
               if (k .ne. i) then
                  xr = xi - x(k)
                  yr = yi - y(k)
                  zr = zi - z(k)
                  kt = class(k)
                  vk = vsolv(k)
                  s2ik = 1.0d0 / s2ace(it,kt)
                  ws2 = wace(it,kt) * s2ik
                  uik4 = uace(it,kt)**4
                  r2 = xr**2 + yr**2 + zr**2
                  r = sqrt(r2)
                  r3 = r2 * r
                  r4 = r2 * r2
                  r6 = r2 * r4
                  rusum = r4 + uik4
                  ratio = r3 / rusum
                  expterm = exp(-r2*s2ik)
                  de1 = -4.0d0 * r * ws2 * expterm
                  de2 = 3.0d0*r2/rusum - 4.0d0*r6/rusum**2
                  de = drb(i) * rbi2 * (de1+vk*ratio**3*de2/pi) / r
c
c     increment the overall continuum solvation derivatives
c
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  des(1,i) = des(1,i) + dedx
                  des(2,i) = des(2,i) + dedy
                  des(3,i) = des(3,i) + dedz
                  des(1,k) = des(1,k) - dedx
                  des(2,k) = des(2,k) - dedy
                  des(3,k) = des(3,k) - dedz
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
         end do
      end if
      return
      end
