 
C     $Id: eurey1.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eurey1  --  bond stretch energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eurey1" calculates the Urey-Bradley interaction energy and
c     its first derivatives with respect to Cartesian coordinates
c
c
      subroutine eurey1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'urey.i'
      include 'urypot.i'
      include 'usage.i'
      include 'virial.i'
      integer i,ia,ib
      real*8 e,de,ideal,force
      real*8 dt,dt2,deddt,fgrp
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 xab,yab,zab,rab
      logical proceed
c
c
c     zero out the Urey-Bradley energy and first derivatives
c
      eub = 0.0d0
      do i = 1, n
         deub(1,i) = 0.0d0
         deub(2,i) = 0.0d0
         deub(3,i) = 0.0d0
      end do
c
c     calculate the Urey-Bradley 1-3 energy and first derivatives
c
      do i = 1, nurey
         ia = iury(1,i)
         ib = iury(2,i)
         ideal = ul(i)
         force = uk(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,0,0,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib))
c
c     compute the value of the 1-3 distance deviation
c
         if (proceed) then
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            if (use_polymer)  call image (xab,yab,zab,0)
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            dt = rab - ideal
            dt2 = dt * dt
            e = ureyunit * force * dt2 * (1.0d0+cury*dt+qury*dt2)
            deddt = 2.0d0 * ureyunit * force * dt
     &                 * (1.0d0+1.5d0*cury*dt+2.0d0*qury*dt2)
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               e = e * fgrp
               deddt = deddt * fgrp
            end if
c
c     compute chain rule terms needed for derivatives
c
            de = deddt / rab
            dedx = de * xab
            dedy = de * yab
            dedz = de * zab
c
c     increment the total Urey-Bradley energy and first derivatives
c
            eub = eub + e
            deub(1,ia) = deub(1,ia) + dedx
            deub(2,ia) = deub(2,ia) + dedy
            deub(3,ia) = deub(3,ia) + dedz
            deub(1,ib) = deub(1,ib) - dedx
            deub(2,ib) = deub(2,ib) - dedy
            deub(3,ib) = deub(3,ib) - dedz
c
c     increment the internal virial tensor components
c
            vxx = xab * dedx
            vyx = yab * dedx
            vzx = zab * dedx
            vyy = yab * dedy
            vzy = zab * dedy
            vzz = zab * dedz
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
      return
      end
