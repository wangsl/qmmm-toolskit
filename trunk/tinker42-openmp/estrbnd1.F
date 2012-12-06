 
C     $Id: estrbnd1.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine estrbnd1   --  stretch-bend energy and derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "estrbnd1" calculates the stretch-bend potential energy and
c     first derivatives with respect to Cartesian coordinates
c
c
      subroutine estrbnd1
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atoms.i'
      include 'bond.i'
      include 'bound.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'math.i'
      include 'strbnd.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,k,istrbnd
      integer ia,ib,ic
      real*8 e,dt,dr,fgrp
      real*8 angle,force
      real*8 dot,cosine
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab,rab
      real*8 xcb,ycb,zcb,rcb
      real*8 xp,yp,zp,rp
      real*8 ddtdxia,ddtdyia,ddtdzia
      real*8 ddtdxic,ddtdyic,ddtdzic
      real*8 ddrdxia,ddrdyia,ddrdzia
      real*8 ddrdxic,ddrdyic,ddrdzic
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 term,terma,termc
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      logical proceed
c
c
c     zero out the energy and first derivative components
c
      eba = 0.0d0
      do i = 1, n
         deba(1,i) = 0.0d0
         deba(2,i) = 0.0d0
         deba(3,i) = 0.0d0
      end do
c
c     calculate the stretch-bend energy and first derivatives
c
      do istrbnd = 1, nstrbnd
         i = isb(1,istrbnd)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         force = ksb(istrbnd)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,0,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic))
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
c
c     compute the value of the bond angle
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            if (use_polymer) then
               call image (xab,yab,zab,0)
               call image (xcb,ycb,zcb,0)
            end if
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
            xp = ycb*zab - zcb*yab
            yp = zcb*xab - xcb*zab
            zp = xcb*yab - ycb*xab
            rp = sqrt(xp*xp + yp*yp + zp*zp)
            if (rp .ne. 0.0d0) then
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / (rab*rcb)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
c
c     find chain rule terms for the bond angle deviation
c
               dt = angle - anat(i)
               terma = -radian / (rab*rab*rp)
               termc = radian / (rcb*rcb*rp)
               ddtdxia = terma * (yab*zp-zab*yp)
               ddtdyia = terma * (zab*xp-xab*zp)
               ddtdzia = terma * (xab*yp-yab*xp)
               ddtdxic = termc * (ycb*zp-zcb*yp)
               ddtdyic = termc * (zcb*xp-xcb*zp)
               ddtdzic = termc * (xcb*yp-ycb*xp)
c
c     find chain rule terms for the bond length deviations
c
               dr = 0.0d0
               terma = 0.0d0
               termc = 0.0d0
               term = stbnunit * force
               j = isb(2,istrbnd)
               k = isb(3,istrbnd)
               if (j .ne. 0) then
                  dr = dr + rab - bl(j)
                  terma = 1.0d0 / rab
               end if
               if (k .ne. 0) then
                  dr = dr + rcb - bl(k)
                  termc = 1.0d0 / rcb
               end if
               ddrdxia = terma * xab
               ddrdyia = terma * yab
               ddrdzia = terma * zab
               ddrdxic = termc * xcb
               ddrdyic = termc * ycb
               ddrdzic = termc * zcb
c
c     scale the interaction based on its group membership
c
               if (use_group)  term = term * fgrp
c
c     get the energy and master chain rule terms for derivatives
c
               e = term * dt * dr
               dedxia = term * (dt*ddrdxia+ddtdxia*dr)
               dedyia = term * (dt*ddrdyia+ddtdyia*dr)
               dedzia = term * (dt*ddrdzia+ddtdzia*dr)
               dedxic = term * (dt*ddrdxic+ddtdxic*dr)
               dedyic = term * (dt*ddrdyic+ddtdyic*dr)
               dedzic = term * (dt*ddrdzic+ddtdzic*dr)
               dedxib = -dedxia - dedxic
               dedyib = -dedyia - dedyic
               dedzib = -dedzia - dedzic
c
c     increment the total stretch-bend energy and derivatives
c
               eba = eba + e
               deba(1,ia) = deba(1,ia) + dedxia
               deba(2,ia) = deba(2,ia) + dedyia
               deba(3,ia) = deba(3,ia) + dedzia
               deba(1,ib) = deba(1,ib) + dedxib
               deba(2,ib) = deba(2,ib) + dedyib
               deba(3,ib) = deba(3,ib) + dedzib
               deba(1,ic) = deba(1,ic) + dedxic
               deba(2,ic) = deba(2,ic) + dedyic
               deba(3,ic) = deba(3,ic) + dedzic
c
c     increment the internal virial tensor components
c
               vxx = xab*dedxia + xcb*dedxic
               vyx = yab*dedxia + ycb*dedxic
               vzx = zab*dedxia + zcb*dedxic
               vyy = yab*dedyia + ycb*dedyic
               vzy = zab*dedyia + zcb*dedyic
               vzz = zab*dedzia + zcb*dedzic
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
      return
      end
