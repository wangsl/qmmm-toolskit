 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eurey2  --  atom-by-atom Urey-Bradley Hessian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eurey2" calculates second derivatives of the Urey-Bradley
c     interaction energy for a single atom at a time
c
c
      subroutine eurey2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'couple.i'
      include 'group.i'
      include 'hessn.i'
      include 'urey.i'
      include 'urypot.i'
      integer i,j,ia,ib,iurey
      real*8 ideal,force,fgrp
      real*8 xab,yab,zab
      real*8 rab,rab2
      real*8 dt,dt2,term
      real*8 termx,termy,termz
      real*8 de,deddt,d2eddt2
      real*8 d2e(3,3)
      logical proceed
c
c
c     compute the Hessian elements of the Urey-Bradley energy
c
      do iurey = 1, nurey
         ia = iury(1,iurey)
         ib = iury(2,iurey)
         ideal = ul(iurey)
         force = uk(iurey)
c
c     decide whether to compute the current interaction
c
         proceed = (i.eq.ia .or. i.eq.ib)
         if (proceed .and. use_group)
     &      call groups (proceed,fgrp,ia,ib,0,0,0,0)
c
c     compute the value of the 1-3 distance deviation
c
         if (proceed) then
            if (i .eq. ib) then
               ib = ia
               ia = i
            end if
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            if (use_polymer)  call image (xab,yab,zab,0)
            rab2 = xab*xab + yab*yab + zab*zab
            rab = sqrt(rab2)
            dt = rab - ideal
            dt2 = dt * dt
            deddt = 2.0d0 * ureyunit * force * dt
     &                 * (1.0d0+1.5d0*cury*dt+2.0d0*qury*dt2)
            d2eddt2 = 2.0d0 * ureyunit * force
     &                   * (1.0d0+3.0d0*cury*dt+6.0d0*qury*dt2)
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               deddt = deddt * fgrp
               d2eddt2 = d2eddt2 * fgrp
            end if
c
c     set the chain rule terms for the Hessian elements
c
            de = deddt / rab
            term = (d2eddt2-de) / rab2
            termx = term * xab
            termy = term * yab
            termz = term * zab
            d2e(1,1) = termx*xab + de
            d2e(1,2) = termx*yab
            d2e(1,3) = termx*zab
            d2e(2,1) = d2e(1,2)
            d2e(2,2) = termy*yab + de
            d2e(2,3) = termy*zab
            d2e(3,1) = d2e(1,3)
            d2e(3,2) = d2e(2,3)
            d2e(3,3) = termz*zab + de
c
c     increment diagonal and non-diagonal Hessian elements
c
            do j = 1, 3
               hessx(j,ia) = hessx(j,ia) + d2e(1,j)
               hessy(j,ia) = hessy(j,ia) + d2e(2,j)
               hessz(j,ia) = hessz(j,ia) + d2e(3,j)
               hessx(j,ib) = hessx(j,ib) - d2e(1,j)
               hessy(j,ib) = hessy(j,ib) - d2e(2,j)
               hessz(j,ib) = hessz(j,ib) - d2e(3,j)
            end do
         end if
      end do
      return
      end
