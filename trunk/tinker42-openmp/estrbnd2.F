 
C     $Id: estrbnd2.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine estrbnd2  --  stretch-bend Hessian; analytical  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "estrbnd2" calculates the stretch-bend potential energy
c     second derivatives with respect to Cartesian coordinates
c
c
      subroutine estrbnd2 (iatom)
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atoms.i'
      include 'bond.i'
      include 'bound.i'
      include 'group.i'
      include 'hessn.i'
      include 'math.i'
      include 'strbnd.i'
      integer i,j,k,iatom
      integer ia,ib,ic,istrbnd
      real*8 angle,force
      real*8 dot,cosine
      real*8 dt,dr,fgrp
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab,rab
      real*8 xcb,ycb,zcb,rcb
      real*8 xp,yp,zp,rp,rp2
      real*8 term,terma,termc
      real*8 xrab,yrab,zrab,rab2
      real*8 xrcb,yrcb,zrcb,rcb2
      real*8 xabp,yabp,zabp
      real*8 xcbp,ycbp,zcbp
      real*8 ddtdxia,ddtdyia,ddtdzia
      real*8 ddtdxib,ddtdyib,ddtdzib
      real*8 ddtdxic,ddtdyic,ddtdzic
      real*8 ddrdxia,ddrdyia,ddrdzia
      real*8 ddrdxib,ddrdyib,ddrdzib
      real*8 ddrdxic,ddrdyic,ddrdzic
      real*8 dtxiaxia,dtxiayia,dtxiazia
      real*8 dtxibxib,dtxibyib,dtxibzib
      real*8 dtxicxic,dtxicyic,dtxiczic
      real*8 dtyiayia,dtyiazia,dtziazia
      real*8 dtyibyib,dtyibzib,dtzibzib
      real*8 dtyicyic,dtyiczic,dtziczic
      real*8 dtxibxia,dtxibyia,dtxibzia
      real*8 dtyibxia,dtyibyia,dtyibzia
      real*8 dtzibxia,dtzibyia,dtzibzia
      real*8 dtxibxic,dtxibyic,dtxibzic
      real*8 dtyibxic,dtyibyic,dtyibzic
      real*8 dtzibxic,dtzibyic,dtzibzic
      real*8 dtxiaxic,dtxiayic,dtxiazic
      real*8 dtyiaxic,dtyiayic,dtyiazic
      real*8 dtziaxic,dtziayic,dtziazic
      real*8 drxiaxia,drxiayia,drxiazia
      real*8 drxibxib,drxibyib,drxibzib
      real*8 drxicxic,drxicyic,drxiczic
      real*8 dryiayia,dryiazia,drziazia
      real*8 dryibyib,dryibzib,drzibzib
      real*8 dryicyic,dryiczic,drziczic
      logical proceed
c
c
c     compute the Hessian elements of the stretch-bends
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
         proceed = (iatom.eq.ia .or. iatom.eq.ib .or. iatom.eq.ic)
         if (proceed .and. use_group)
     &      call groups (proceed,fgrp,ia,ib,ic,0,0,0)
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
c     first derivatives of angle with respect to coordinates
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
               ddtdxib = -ddtdxia - ddtdxic
               ddtdyib = -ddtdyia - ddtdyic
               ddtdzib = -ddtdzia - ddtdzic
c
c     abbreviations used in defining chain rule terms
c
               rab2 = 2.0d0 / (rab*rab)
               xrab = xab * rab2
               yrab = yab * rab2
               zrab = zab * rab2
               rcb2 = 2.0d0 / (rcb*rcb)
               xrcb = xcb * rcb2
               yrcb = ycb * rcb2
               zrcb = zcb * rcb2
               rp2 = 1.0d0 / (rp*rp)
               xabp = (yab*zp-zab*yp) * rp2
               yabp = (zab*xp-xab*zp) * rp2
               zabp = (xab*yp-yab*xp) * rp2
               xcbp = (ycb*zp-zcb*yp) * rp2
               ycbp = (zcb*xp-xcb*zp) * rp2
               zcbp = (xcb*yp-ycb*xp) * rp2
c
c     second derivatives of angle with respect to coordinates
c
               dtxiaxia = terma*(xab*xcb-dot) + ddtdxia*(xcbp-xrab)
               dtxiayia = terma*(zp+yab*xcb) + ddtdxia*(ycbp-yrab)
               dtxiazia = terma*(zab*xcb-yp) + ddtdxia*(zcbp-zrab)
               dtyiayia = terma*(yab*ycb-dot) + ddtdyia*(ycbp-yrab)
               dtyiazia = terma*(xp+zab*ycb) + ddtdyia*(zcbp-zrab)
               dtziazia = terma*(zab*zcb-dot) + ddtdzia*(zcbp-zrab)
               dtxicxic = termc*(dot-xab*xcb) - ddtdxic*(xabp+xrcb)
               dtxicyic = termc*(zp-ycb*xab) - ddtdxic*(yabp+yrcb)
               dtxiczic = -termc*(yp+zcb*xab) - ddtdxic*(zabp+zrcb)
               dtyicyic = termc*(dot-yab*ycb) - ddtdyic*(yabp+yrcb)
               dtyiczic = termc*(xp-zcb*yab) - ddtdyic*(zabp+zrcb)
               dtziczic = termc*(dot-zab*zcb) - ddtdzic*(zabp+zrcb)
               dtxiaxic = terma*(yab*yab+zab*zab) - ddtdxia*xabp
               dtxiayic = -terma*xab*yab - ddtdxia*yabp
               dtxiazic = -terma*xab*zab - ddtdxia*zabp
               dtyiaxic = -terma*xab*yab - ddtdyia*xabp
               dtyiayic = terma*(xab*xab+zab*zab) - ddtdyia*yabp
               dtyiazic = -terma*yab*zab - ddtdyia*zabp
               dtziaxic = -terma*xab*zab - ddtdzia*xabp
               dtziayic = -terma*yab*zab - ddtdzia*yabp
               dtziazic = terma*(xab*xab+yab*yab) - ddtdzia*zabp
c
c     more angle deviation derivatives resulting from symmetry
c
               dtxibxia = -dtxiaxia - dtxiaxic
               dtxibyia = -dtxiayia - dtyiaxic
               dtxibzia = -dtxiazia - dtziaxic
               dtyibxia = -dtxiayia - dtxiayic
               dtyibyia = -dtyiayia - dtyiayic
               dtyibzia = -dtyiazia - dtziayic
               dtzibxia = -dtxiazia - dtxiazic
               dtzibyia = -dtyiazia - dtyiazic
               dtzibzia = -dtziazia - dtziazic
               dtxibxic = -dtxicxic - dtxiaxic
               dtxibyic = -dtxicyic - dtxiayic
               dtxibzic = -dtxiczic - dtxiazic
               dtyibxic = -dtxicyic - dtyiaxic
               dtyibyic = -dtyicyic - dtyiayic
               dtyibzic = -dtyiczic - dtyiazic
               dtzibxic = -dtxiczic - dtziaxic
               dtzibyic = -dtyiczic - dtziayic
               dtzibzic = -dtziczic - dtziazic
               dtxibxib = -dtxibxia - dtxibxic
               dtxibyib = -dtxibyia - dtxibyic
               dtxibzib = -dtxibzia - dtxibzic
               dtyibyib = -dtyibyia - dtyibyic
               dtyibzib = -dtyibzia - dtyibzic
               dtzibzib = -dtzibzia - dtzibzic
c
c     compute the values of the bond length deviations
c
               dr = 0.0d0
               terma = 0.0d0
               termc = 0.0d0
               term = stbnunit * force
               j = isb(2,istrbnd)
               k = isb(3,istrbnd)
               if (j .ne. 0) then
                  dr = dr + term*(rab-bl(j))
                  terma = term / rab
               end if
               if (k .ne. 0) then
                  dr = dr + term*(rcb-bl(k))
                  termc = term / rcb
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dr = dr * fgrp
                  terma = terma * fgrp
                  termc = termc * fgrp
               end if
c
c     first derivatives of bond length with respect to coordinates
c
               ddrdxia = terma * xab
               ddrdyia = terma * yab
               ddrdzia = terma * zab
               ddrdxic = termc * xcb
               ddrdyic = termc * ycb
               ddrdzic = termc * zcb
               ddrdxib = -ddrdxia - ddrdxic
               ddrdyib = -ddrdyia - ddrdyic
               ddrdzib = -ddrdzia - ddrdzic
c
c     abbreviations used in defining chain rule terms
c
               xab = xab / rab
               yab = yab / rab
               zab = zab / rab
               xcb = xcb / rcb
               ycb = ycb / rcb
               zcb = zcb / rcb
c
c     second derivatives of bond length with respect to coordinates
c
               drxiaxia = terma * (1.0d0-xab**2)
               drxiayia = -terma * xab*yab
               drxiazia = -terma * xab*zab
               dryiayia = terma * (1.0d0-yab**2)
               dryiazia = -terma * yab*zab
               drziazia = terma * (1.0d0-zab**2)
               drxicxic = termc * (1.0d0-xcb**2)
               drxicyic = -termc * xcb*ycb
               drxiczic = -termc * xcb*zcb
               dryicyic = termc * (1.0d0-ycb**2)
               dryiczic = -termc * ycb*zcb
               drziczic = termc * (1.0d0-zcb**2)
               drxibxib = drxiaxia + drxicxic
               drxibyib = drxiayia + drxicyic
               drxibzib = drxiazia + drxiczic
               dryibyib = dryiayia + dryicyic
               dryibzib = dryiazia + dryiczic
               drzibzib = drziazia + drziczic
c
c     increment diagonal and non-diagonal Hessian elements
c
               if (ia .eq. iatom) then
                  hessx(1,ia) = hessx(1,ia) + dt*drxiaxia + dr*dtxiaxia
     &                             + 2.0d0*ddtdxia*ddrdxia
                  hessx(2,ia) = hessx(2,ia) + dt*drxiayia + dr*dtxiayia
     &                             + ddtdxia*ddrdyia + ddtdyia*ddrdxia
                  hessx(3,ia) = hessx(3,ia) + dt*drxiazia + dr*dtxiazia
     &                             + ddtdxia*ddrdzia + ddtdzia*ddrdxia
                  hessy(1,ia) = hessy(1,ia) + dt*drxiayia + dr*dtxiayia
     &                             + ddtdyia*ddrdxia + ddtdxia*ddrdyia
                  hessy(2,ia) = hessy(2,ia) + dt*dryiayia + dr*dtyiayia
     &                             + 2.0d0*ddtdyia*ddrdyia
                  hessy(3,ia) = hessy(3,ia) + dt*dryiazia + dr*dtyiazia
     &                             + ddtdyia*ddrdzia + ddtdzia*ddrdyia
                  hessz(1,ia) = hessz(1,ia) + dt*drxiazia + dr*dtxiazia
     &                             + ddtdzia*ddrdxia + ddtdxia*ddrdzia
                  hessz(2,ia) = hessz(2,ia) + dt*dryiazia + dr*dtyiazia
     &                             + ddtdzia*ddrdyia + ddtdyia*ddrdzia
                  hessz(3,ia) = hessz(3,ia) + dt*drziazia + dr*dtziazia
     &                             + 2.0d0*ddtdzia*ddrdzia
                  hessx(1,ib) = hessx(1,ib) - dt*drxiaxia + dr*dtxibxia
     &                             + ddtdxia*ddrdxib + ddtdxib*ddrdxia
                  hessx(2,ib) = hessx(2,ib) - dt*drxiayia + dr*dtxibyia
     &                             + ddtdxia*ddrdyib + ddtdyib*ddrdxia
                  hessx(3,ib) = hessx(3,ib) - dt*drxiazia + dr*dtxibzia
     &                             + ddtdxia*ddrdzib + ddtdzib*ddrdxia
                  hessy(1,ib) = hessy(1,ib) - dt*drxiayia + dr*dtyibxia
     &                             + ddtdyia*ddrdxib + ddtdxib*ddrdyia
                  hessy(2,ib) = hessy(2,ib) - dt*dryiayia + dr*dtyibyia
     &                             + ddtdyia*ddrdyib + ddtdyib*ddrdyia
                  hessy(3,ib) = hessy(3,ib) - dt*dryiazia + dr*dtyibzia
     &                             + ddtdyia*ddrdzib + ddtdzib*ddrdyia
                  hessz(1,ib) = hessz(1,ib) - dt*drxiazia + dr*dtzibxia
     &                             + ddtdzia*ddrdxib + ddtdxib*ddrdzia
                  hessz(2,ib) = hessz(2,ib) - dt*dryiazia + dr*dtzibyia
     &                             + ddtdzia*ddrdyib + ddtdyib*ddrdzia
                  hessz(3,ib) = hessz(3,ib) - dt*drziazia + dr*dtzibzia
     &                             + ddtdzia*ddrdzib + ddtdzib*ddrdzia
                  hessx(1,ic) = hessx(1,ic) + dr*dtxiaxic
     &                             + ddtdxia*ddrdxic + ddtdxic*ddrdxia
                  hessx(2,ic) = hessx(2,ic) + dr*dtxiayic
     &                             + ddtdxia*ddrdyic + ddtdyic*ddrdxia
                  hessx(3,ic) = hessx(3,ic) + dr*dtxiazic
     &                             + ddtdxia*ddrdzic + ddtdzic*ddrdxia
                  hessy(1,ic) = hessy(1,ic) + dr*dtyiaxic
     &                             + ddtdyia*ddrdxic + ddtdxic*ddrdyia
                  hessy(2,ic) = hessy(2,ic) + dr*dtyiayic
     &                             + ddtdyia*ddrdyic + ddtdyic*ddrdyia
                  hessy(3,ic) = hessy(3,ic) + dr*dtyiazic
     &                             + ddtdyia*ddrdzic + ddtdzic*ddrdyia
                  hessz(1,ic) = hessz(1,ic) + dr*dtziaxic
     &                             + ddtdzia*ddrdxic + ddtdxic*ddrdzia
                  hessz(2,ic) = hessz(2,ic) + dr*dtziayic
     &                             + ddtdzia*ddrdyic + ddtdyic*ddrdzia
                  hessz(3,ic) = hessz(3,ic) + dr*dtziazic
     &                             + ddtdzia*ddrdzic + ddtdzic*ddrdzia
               else if (ib .eq. iatom) then
                  hessx(1,ib) = hessx(1,ib) + dt*drxibxib + dr*dtxibxib
     &                             + 2.0d0*ddtdxib*ddrdxib
                  hessx(2,ib) = hessx(2,ib) + dt*drxibyib + dr*dtxibyib
     &                             + ddtdxib*ddrdyib + ddtdyib*ddrdxib
                  hessx(3,ib) = hessx(3,ib) + dt*drxibzib + dr*dtxibzib
     &                             + ddtdxib*ddrdzib + ddtdzib*ddrdxib
                  hessy(1,ib) = hessy(1,ib) + dt*drxibyib + dr*dtxibyib
     &                             + ddtdyib*ddrdxib + ddtdxib*ddrdyib
                  hessy(2,ib) = hessy(2,ib) + dt*dryibyib + dr*dtyibyib
     &                             + 2.0d0*ddtdyib*ddrdyib
                  hessy(3,ib) = hessy(3,ib) + dt*dryibzib + dr*dtyibzib
     &                             + ddtdyib*ddrdzib + ddtdzib*ddrdyib
                  hessz(1,ib) = hessz(1,ib) + dt*drxibzib + dr*dtxibzib
     &                             + ddtdzib*ddrdxib + ddtdxib*ddrdzib
                  hessz(2,ib) = hessz(2,ib) + dt*dryibzib + dr*dtyibzib
     &                             + ddtdzib*ddrdyib + ddtdyib*ddrdzib
                  hessz(3,ib) = hessz(3,ib) + dt*drzibzib + dr*dtzibzib
     &                             + 2.0d0*ddtdzib*ddrdzib
                  hessx(1,ia) = hessx(1,ia) - dt*drxiaxia + dr*dtxibxia
     &                             + ddtdxib*ddrdxia + ddtdxia*ddrdxib
                  hessx(2,ia) = hessx(2,ia) - dt*drxiayia + dr*dtxibyia
     &                             + ddtdxib*ddrdyia + ddtdyia*ddrdxib
                  hessx(3,ia) = hessx(3,ia) - dt*drxiazia + dr*dtxibzia
     &                             + ddtdxib*ddrdzia + ddtdzia*ddrdxib
                  hessy(1,ia) = hessy(1,ia) - dt*drxiayia + dr*dtyibxia
     &                             + ddtdyib*ddrdxia + ddtdxia*ddrdyib
                  hessy(2,ia) = hessy(2,ia) - dt*dryiayia + dr*dtyibyia
     &                             + ddtdyib*ddrdyia + ddtdyia*ddrdyib
                  hessy(3,ia) = hessy(3,ia) - dt*dryiazia + dr*dtyibzia
     &                             + ddtdyib*ddrdzia + ddtdzia*ddrdyib
                  hessz(1,ia) = hessz(1,ia) - dt*drxiazia + dr*dtzibxia
     &                             + ddtdzib*ddrdxia + ddtdxia*ddrdzib
                  hessz(2,ia) = hessz(2,ia) - dt*dryiazia + dr*dtzibyia
     &                             + ddtdzib*ddrdyia + ddtdyia*ddrdzib
                  hessz(3,ia) = hessz(3,ia) - dt*drziazia + dr*dtzibzia
     &                             + ddtdzib*ddrdzia + ddtdzia*ddrdzib
                  hessx(1,ic) = hessx(1,ic) - dt*drxicxic + dr*dtxibxic
     &                             + ddtdxib*ddrdxic + ddtdxic*ddrdxib
                  hessx(2,ic) = hessx(2,ic) - dt*drxicyic + dr*dtxibyic
     &                             + ddtdxib*ddrdyic + ddtdyic*ddrdxib
                  hessx(3,ic) = hessx(3,ic) - dt*drxiczic + dr*dtxibzic
     &                             + ddtdxib*ddrdzic + ddtdzic*ddrdxib
                  hessy(1,ic) = hessy(1,ic) - dt*drxicyic + dr*dtyibxic
     &                             + ddtdyib*ddrdxic + ddtdxic*ddrdyib
                  hessy(2,ic) = hessy(2,ic) - dt*dryicyic + dr*dtyibyic
     &                             + ddtdyib*ddrdyic + ddtdyic*ddrdyib
                  hessy(3,ic) = hessy(3,ic) - dt*dryiczic + dr*dtyibzic
     &                             + ddtdyib*ddrdzic + ddtdzic*ddrdyib
                  hessz(1,ic) = hessz(1,ic) - dt*drxiczic + dr*dtzibxic
     &                             + ddtdzib*ddrdxic + ddtdxic*ddrdzib
                  hessz(2,ic) = hessz(2,ic) - dt*dryiczic + dr*dtzibyic
     &                             + ddtdzib*ddrdyic + ddtdyic*ddrdzib
                  hessz(3,ic) = hessz(3,ic) - dt*drziczic + dr*dtzibzic
     &                             + ddtdzib*ddrdzic + ddtdzic*ddrdzib
               else if (ic .eq. iatom) then
                  hessx(1,ic) = hessx(1,ic) + dt*drxicxic + dr*dtxicxic
     &                             + 2.0d0*ddtdxic*ddrdxic
                  hessx(2,ic) = hessx(2,ic) + dt*drxicyic + dr*dtxicyic
     &                             + ddtdxic*ddrdyic + ddtdyic*ddrdxic
                  hessx(3,ic) = hessx(3,ic) + dt*drxiczic + dr*dtxiczic
     &                             + ddtdxic*ddrdzic + ddtdzic*ddrdxic
                  hessy(1,ic) = hessy(1,ic) + dt*drxicyic + dr*dtxicyic
     &                             + ddtdyic*ddrdxic + ddtdxic*ddrdyic
                  hessy(2,ic) = hessy(2,ic) + dt*dryicyic + dr*dtyicyic
     &                             + 2.0d0*ddtdyic*ddrdyic
                  hessy(3,ic) = hessy(3,ic) + dt*dryiczic + dr*dtyiczic
     &                             + ddtdyic*ddrdzic + ddtdzic*ddrdyic
                  hessz(1,ic) = hessz(1,ic) + dt*drxiczic + dr*dtxiczic
     &                             + ddtdzic*ddrdxic + ddtdxic*ddrdzic
                  hessz(2,ic) = hessz(2,ic) + dt*dryiczic + dr*dtyiczic
     &                             + ddtdzic*ddrdyic + ddtdyic*ddrdzic
                  hessz(3,ic) = hessz(3,ic) + dt*drziczic + dr*dtziczic
     &                             + 2.0d0*ddtdzic*ddrdzic
                  hessx(1,ib) = hessx(1,ib) - dt*drxicxic + dr*dtxibxic
     &                             + ddtdxic*ddrdxib + ddtdxib*ddrdxic
                  hessx(2,ib) = hessx(2,ib) - dt*drxicyic + dr*dtxibyic
     &                             + ddtdxic*ddrdyib + ddtdyib*ddrdxic
                  hessx(3,ib) = hessx(3,ib) - dt*drxiczic + dr*dtxibzic
     &                             + ddtdxic*ddrdzib + ddtdzib*ddrdxic
                  hessy(1,ib) = hessy(1,ib) - dt*drxicyic + dr*dtyibxic
     &                             + ddtdyic*ddrdxib + ddtdxib*ddrdyic
                  hessy(2,ib) = hessy(2,ib) - dt*dryicyic + dr*dtyibyic
     &                             + ddtdyic*ddrdyib + ddtdyib*ddrdyic
                  hessy(3,ib) = hessy(3,ib) - dt*dryiczic + dr*dtyibzic
     &                             + ddtdyic*ddrdzib + ddtdzib*ddrdyic
                  hessz(1,ib) = hessz(1,ib) - dt*drxiczic + dr*dtzibxic
     &                             + ddtdzic*ddrdxib + ddtdxib*ddrdzic
                  hessz(2,ib) = hessz(2,ib) - dt*dryiczic + dr*dtzibyic
     &                             + ddtdzic*ddrdyib + ddtdyib*ddrdzic
                  hessz(3,ib) = hessz(3,ib) - dt*drziczic + dr*dtzibzic
     &                             + ddtdzic*ddrdzib + ddtdzib*ddrdzic
                  hessx(1,ia) = hessx(1,ia) + dr*dtxiaxic
     &                             + ddtdxic*ddrdxia + ddtdxia*ddrdxic
                  hessx(2,ia) = hessx(2,ia) + dr*dtyiaxic
     &                             + ddtdxic*ddrdyia + ddtdyia*ddrdxic
                  hessx(3,ia) = hessx(3,ia) + dr*dtziaxic
     &                             + ddtdxic*ddrdzia + ddtdzia*ddrdxic
                  hessy(1,ia) = hessy(1,ia) + dr*dtxiayic
     &                             + ddtdyic*ddrdxia + ddtdxia*ddrdyic
                  hessy(2,ia) = hessy(2,ia) + dr*dtyiayic
     &                             + ddtdyic*ddrdyia + ddtdyia*ddrdyic
                  hessy(3,ia) = hessy(3,ia) + dr*dtziayic
     &                             + ddtdyic*ddrdzia + ddtdzia*ddrdyic
                  hessz(1,ia) = hessz(1,ia) + dr*dtxiazic
     &                             + ddtdzic*ddrdxia + ddtdxia*ddrdzic
                  hessz(2,ia) = hessz(2,ia) + dr*dtyiazic
     &                             + ddtdzic*ddrdyia + ddtdyia*ddrdzic
                  hessz(3,ia) = hessz(3,ia) + dr*dtziazic
     &                             + ddtdzic*ddrdzia + ddtdzia*ddrdzic
               end if
            end if
         end if
      end do
      return
      end
