 
C     $Id: eopbend3.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eopbend3  --  out-of-plane bending & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eopbend3" computes the out-of-plane bend potential energy at
c     trigonal centers via a Wilson-Decius-Cross angle bend; also
c     partitions the energy among the atoms
c
c
      subroutine eopbend3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'opbend.i'
      include 'usage.i'
      integer i,iopbend
      integer ia,ib,ic,id
      real*8 e,angle,force
      real*8 cosine,fgrp
      real*8 dt,dt2,dt3,dt4
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 xdb,ydb,zdb
      real*8 xad,yad,zad
      real*8 xcd,ycd,zcd
      real*8 rdb2,rad2,rcd2
      real*8 cc,ee,bkk2
      logical header,huge,proceed
c
c
c     zero out the out-of-plane bend energy and partitioning
c
      neopb = 0
      eopb = 0.0d0
      do i = 1, n
         aeopb(i) = 0.0d0
      end do
      header = .true.
c
c     calculate the out-of-plane bending energy term
c
      do iopbend = 1, nopbend
         i = iopb(iopbend)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(4,i)
         force = kopb(iopbend)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
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
            xid = x(id)
            yid = y(id)
            zid = z(id)
c
c     compute the out-of-plane bending angle
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdb = xid - xib
            ydb = yid - yib
            zdb = zid - zib
            xad = xia - xid
            yad = yia - yid
            zad = zia - zid
            xcd = xic - xid
            ycd = yic - yid
            zcd = zic - zid
            if (use_polymer) then
               call image (xab,yab,zab,0)
               call image (xcb,ycb,zcb,0)
               call image (xdb,ydb,zdb,0)
               call image (xad,yad,zad,0)
               call image (xcd,ycd,zcd,0)
            end if
            rdb2 = xdb*xdb + ydb*ydb + zdb*zdb
            rad2 = xad*xad + yad*yad + zad*zad
            rcd2 = xcd*xcd + ycd*ycd + zcd*zcd
            ee = xab*(ydb*zcb-zdb*ycb) +  yab*(zdb*xcb-xdb*zcb)
     &                      + zab*(xdb*ycb-ydb*xcb)
            cc = rad2*rcd2 - (xad*xcd+yad*ycd+zad*zcd)**2
            if (rdb2.ne.0.0d0 .and. cc.ne.0.0d0) then
               bkk2 = rdb2 - ee*ee/cc
               cosine = sqrt(bkk2/rdb2)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
c
c     find the out-of-plane angle bending energy
c
               dt = angle
               dt2 = dt * dt
               dt3 = dt2 * dt
               dt4 = dt2 * dt2
               e = opbunit * force * dt2
     &                * (1.0d0+cang*dt+qang*dt2+pang*dt3+sang*dt4)
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total out-of-plane bending energy
c
               neopb = neopb + 1
               eopb = eopb + e
               aeopb(ib) = aeopb(ib) + e
c
c     print a message if the energy of this interaction is large
c
               huge = (e .gt. 2.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Individual Out-of-Plane Bending',
     &                          ' Interactions :',
     &                       //,' Type',13x,'Atom Names',30x,'Angle',
     &                          6x,'Energy',/)
                  end if
                  write (iout,20)  ib,name(ib),id,name(id),angle,e
   20             format (' O-P-Bend',3x,i5,'-',a3,1x,i5,'-',a3,
     &                       22x,f10.4,f12.4)
               end if
            end if
         end if
      end do
      return
      end
