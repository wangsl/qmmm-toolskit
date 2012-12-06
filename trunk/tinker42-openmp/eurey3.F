 
C     $Id: eurey3.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine eurey3  --  Urey-Bradley energy & analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "eurey3" calculates the Urey-Bradley energy; also
c     partitions the energy among the atoms
c
c
      subroutine eurey3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'urey.i'
      include 'urypot.i'
      include 'usage.i'
      integer i,ia,ib
      real*8 e,ideal,force
      real*8 dt,dt2,fgrp
      real*8 xab,yab,zab,rab
      logical proceed
      logical header,huge
c
c
c     zero out the Urey-Bradley energy and partitioning terms
c
      neub = 0
      eub = 0.0d0
      do i = 1, n
         aeub(i) = 0.0d0
      end do
      header = .true.
c
c     calculate the Urey-Bradley 1-3 energy term
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
c
c     calculate the Urey-Bradley energy for this interaction
c
            e = ureyunit * force * dt2 * (1.0d0+cury*dt+qury*dt2)
c
c     scale the interaction based on its group membership
c
            if (use_group)  e = e * fgrp
c
c     increment the total Urey-Bradley energy
c
            neub = neub + 1
            eub = eub + e
            aeub(ia) = aeub(ia) + 0.5d0*e
            aeub(ib) = aeub(ib) + 0.5d0*e
c
c     print a message if the energy of this interaction is large
c
            huge = (e .gt. 5.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual Urey-Bradley Interactions :',
     &                    //,' Type',13x,'Atom Names',20x,'Ideal',
     &                       4x,'Actual',6x,'Energy',/)
               end if
               write (iout,20)  ia,name(ia),ib,name(ib),ideal,rab,e
   20          format (' UreyBrad',3x,i5,'-',a3,1x,i5,'-',a3,
     &                    12x,2f10.4,f12.4)
            end if
         end if
      end do
      return
      end
