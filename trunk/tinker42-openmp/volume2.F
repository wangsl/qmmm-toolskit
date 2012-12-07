 
C     $Id$
c
c
c     ################################################################
c     ##  COPYRIGHT (C) 1990 by Craig Kundrot & Jay William Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine volume2  --  Cartesian excluded volume Hessian  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "volume2" calculates second derivatives of the total excluded
c     volume with respect to the Cartesian coordinates of the atoms
c
c     literature reference:
c
c     C. E. Kundrot, J. W. Ponder and F. M. Richards, "Algorithms for
c     Calculating Excluded Volume and Its Derivatives as a Function
c     of Molecular Conformation and Their Use in Energy Minimization",
c     Journal of Computational Chemistry, 12, 402-409 (1991)
c
c
      subroutine volume2 (iatom,radius,probe,xhess,yhess,zhess)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'math.i'
      integer maxarc
      parameter (maxarc=300)
      integer i,j,k,m
      integer in,iaa,ibb
      integer iatom,narc
      integer iblock,itemp
      integer idtemp,idfirst
      integer nnear,id(0:2)
      integer inear(maxarc)
      integer arciatom(maxarc)
      integer arcfatom(maxarc)
      real*8 probe,zstep,xr,yr,zr
      real*8 ztop,ztopshave,zstart
      real*8 aa,bb,temp,tempf
      real*8 phi1,phi2,phiold
      real*8 theta1,theta2,firsti
      real*8 zgrid,rsec2r,rsecr
      real*8 pix2,dist2,rcut2
      real*8 rr,rrx2,rrsq
      real*8 alpha,beta,gamma
      real*8 ti,tf,ri,s2,b,cosine
      real*8 rinsq,rsecn,rsec2n
      real*8 cos1,cos2,sin1,sin2
      real*8 phi_xy,phi_z
      real*8 delx(2),dely(2),delz(2)
      real*8 r_s(2),r_s2(2),u(2)
      real*8 r(0:2),r_r(0:2)
      real*8 duds(2),dudr(2)
      real*8 u_term(2)
      real*8 dfdtheta(3,2)
      real*8 dthetadx(2,3,0:2)
      real*8 dalphdx(2,3,0:2)
      real*8 dbetadx(2,2,0:2)
      real*8 dudx(2,3,0:2)
      real*8 dsdx(2,2,0:2)
      real*8 drdz(2,0:2)
      real*8 arci(maxarc)
      real*8 arcf(maxarc)
      real*8 dx(maxarc),dy(maxarc)
      real*8 dsq(maxarc),d(maxarc)
      real*8 radius(maxatm)
      real*8 vdwrad(maxatm)
      real*8 xhess(3,maxatm)
      real*8 yhess(3,maxatm)
      real*8 zhess(3,maxatm)
      logical covered
c
c
c     fix the stepsize in the z-direction; this value sets
c     the accuracy of the numerical derivatives; zstep=0.06
c     is a good balance between compute time and accuracy
c
      zstep = 0.0601d0
c
c     zero out the Hessian elements for current atom
c
      do i = 1, n
         do j = 1, 3
            xhess(j,i) = 0.0d0
            yhess(j,i) = 0.0d0
            zhess(j,i) = 0.0d0
         end do
      end do
      if (radius(iatom) .eq. 0.0d0)  return
      pix2 = 2.0d0 * pi
c
c     assign van der Waals radii to the atoms; note that
c     the radii are incremented by the size of the probe
c
      do i = 1, n
         vdwrad(i) = radius(i)
         if (vdwrad(i) .ne. 0.0d0)  vdwrad(i) = vdwrad(i) + probe
      end do
c
c     set the radius and coordinates for current atom
c
      rr = vdwrad(iatom)
      rrx2 = 2.0d0 * rr
      rrsq = rr**2
      xr = x(iatom)
      yr = y(iatom)
      zr = z(iatom)
c
c     select potential intersecting atoms
c
      nnear = 1
      do j = 1, n
         if (j.ne.iatom .and. vdwrad(j).ne.0.0d0) then
            dx(nnear) = x(j) - xr
            dy(nnear) = y(j) - yr
            dsq(nnear) = dx(nnear)**2 + dy(nnear)**2
            dist2 = dsq(nnear) + (z(j)-zr)**2
            rcut2 = (vdwrad(j) + rr)**2
            if (dist2 .lt. rcut2) then
               d(nnear) = sqrt(dsq(nnear))
               inear(nnear) = j
               nnear = nnear + 1
               if (nnear .gt. maxarc) then
                  write (iout,10)
   10             format (/,' VOLUME2  --  Increase',
     &                       ' the Value of MAXARC')
                  call fatal
               end if
            end if
         end if
      end do
      nnear = nnear - 1
c
c     determine the z resolution
c
      if (nnear .ne. 0) then
         ztop = zr + rr
         ztopshave = ztop - zstep
         zgrid = zr - rr
c
c     half of the part not covered by the planes
c
         zgrid = zgrid + (0.5d0*(rrx2-(int(rrx2/zstep)*zstep)))
         zstart = zgrid
c
c     section atom spheres perpendicular to the z axis
c
         dowhile (zgrid .le. ztop)
c
c     "rsecr" is radius of current atom sphere on the z-plane
c
            rsec2r = rrsq - (zgrid-zr)**2
            if (rsec2r .lt. 0.0d0) then
               rsec2r = 0.000001d0
            end if
            rsecr = sqrt(rsec2r)
            if (zgrid .ge. ztopshave) then
               phi1 = 0.0d0
            else
               phi1 = acos(((zgrid+0.5d0*zstep)-zr) / rr)
            end if
            if (zgrid .eq. zstart) then
               phi2 = pi
            else
               phi2 = phiold
            end if
c
c     check intersections of neighbor circles
c
            k = 0
            narc = 0
            covered = .false.
            dowhile (.not.covered .and. k.lt.nnear
     &                   .and. narc.lt.maxarc)
               k = k + 1
               in = inear(k)
               rinsq = vdwrad(in)**2
               rsec2n = rinsq - (zgrid-z(in))**2
               if (rsec2n .gt. 0.0d0) then
                  rsecn = sqrt(rsec2n)
                  if (d(k) .lt. rsecr+rsecn) then
                     b = rsecr - rsecn
                     if (d(k) .le. abs(b)) then
                        if (b .lt. 0.0d0) then
                           narc = 1
                           arci(narc) = 0.0d0
                           arcf(narc) = pix2
                           arciatom(narc) = in
                           arcfatom(narc) = in
                           covered = .true.
                        end if
                     else
                        narc = narc + 1
                        if (narc .gt. maxarc) then
                           write (iout,20)
   20                      format (/,' VOLUME2  -- Increase',
     &                                ' the Value of MAXARC')
                           call fatal
                        else
c
c     initial and final arc endpoints are found for intersection
c     of "ir" circle with another circle contained in same plane;
c     the initial endpoint of the enclosed arc is stored in "arci",
c     the final endpoint in "arcf"; get "cosine" via law of cosines
c
                           cosine = (dsq(k)+rsec2r-rsec2n) /
     &                                      (2.0d0*d(k)*rsecr)
                           cosine = min(1.0d0,max(-1.0d0,cosine))
c
c     "alpha" is the angle between a line containing either point
c     of intersection and the reference circle center and the
c     line containing both circle centers; "beta" is the angle
c     between the line containing both circle centers and x-axis
c
                           alpha = acos(cosine)
                           if (dx(k) .eq. 0.0d0) then
                              gamma = 0.5d0 * pi
                           else
                              gamma = atan(abs(dy(k)/dx(k)))
                           end if
                           if (dy(k) .gt. 0.0d0) then
                              if (dx(k) .gt. 0.0d0) then
                                 beta = gamma
                              else
                                 beta = pi - gamma
                              end if
                           else
                              if (dx(k) .gt. 0.0d0) then
                                 beta = pix2 - gamma
                              else
                                 beta = pi + gamma
                              end if
                           end if
c
c     finally, the arc endpoints
c
                           ti = beta - alpha
                           tf = beta + alpha
                           if (ti .lt. 0.0d0)  ti = ti + pix2
                           if (tf .gt. pix2)  tf = tf - pix2
                           arci(narc) = ti
                           arciatom(narc) = in
                           arcfatom(narc) = in
                           if (tf .lt. ti) then
                              arcf(narc) = pix2
                              narc = narc + 1
                              arci(narc) = 0.0d0
                              arciatom(narc) = in
                              arcfatom(narc) = in
                           end if
                           arcf(narc) = tf
                        end if
                     end if
                  end if
               end if
            end do
c
c     find the pre-area and pre-forces on this section (band)
c     through sphere "ir"; the "pre-" means a multiplicative
c     factor is yet to be applied
c
            if (narc .ne. 0) then
c
c     general case; sort arc endpoints
c
               k = 1
               dowhile (k .lt. narc)
                  aa = arci(k)
                  bb = arcf(k)
                  iaa = arciatom(k)
                  ibb = arcfatom(k)
                  temp = 10000000.0d0
                  do i = k, narc
                     if (arci(i) .le. temp) then
                        temp = arci(i)
                        itemp = i
                     end if
                  end do
                  arci(k) = arci(itemp)
                  arcf(k) = arcf(itemp)
                  arciatom(k) = arciatom(itemp)
                  arcfatom(k) = arcfatom(itemp)
                  arci(itemp) = aa
                  arcf(itemp) = bb
                  arciatom(itemp) = iaa
                  arcfatom(itemp) = ibb
                  k = k + 1
               end do
c
c     eliminate overlapping arc endpoints;
c     first, consolidate the occluded arcs
c
               m = 1
               tempf = arcf(1)
               idtemp = arcfatom(1)
               do k = 2, narc
                  if (tempf .lt. arci(k)) then
                     arcf(m) = tempf
                     arcfatom(m) = idtemp
                     m = m + 1
                     arci(m) = arci(k)
                     arciatom(m) = arciatom(k)
                     tempf = arcf(k)
                     idtemp = arcfatom(k)
                  else if (tempf .lt. arcf(k)) then
                     tempf = arcf(k)
                     idtemp = arcfatom(k)
                  end if
               end do
               arcf(m) = tempf
               arcfatom(m) = idtemp
               narc = m
c
c     change occluded arcs to accessible arcs
c
               if (narc .eq. 1) then
                  if (arci(1).eq.0.0d0 .and. arcf(1).eq.pix2) then
                     narc = 0
                  else
                     firsti = arci(1)
                     idfirst = arciatom(1)
                     arci(1) = arcf(1)
                     arciatom(1) = arcfatom(1)
                     arcf(1) = firsti + pix2
                     arcfatom(1) = idfirst
                  end if
               else
                  firsti = arci(1)
                  idfirst = arciatom(1)
                  do k = 1, narc-1
                     arci(k) = arcf(k)
                     arciatom(k) = arcfatom(k)
                     arcf(k) = arci(k+1)
                     arcfatom(k) = arciatom(k+1)
                  end do
c
c     check gap between first and last arcs; if the
c     occluded arc crossed zero, then no accessible arc
c
                  if (firsti.eq.0.0d0 .and. arcf(narc).eq.pix2) then
                     narc = narc - 1
                  else
                     arci(narc) = arcf(narc)
                     arciatom(narc) = arcfatom(narc)
                     arcf(narc) = firsti
                     arcfatom(narc) = idfirst
                  end if
               end if
c
c     setup prior to application of chain rule
c
               do k = 1, narc
                  ri = sqrt(rrsq - (zgrid-zr)**2)
                  do i = 1, 2
                     if (i .eq. 1) then
                        id(1) = arciatom(k)
                     else
                        id(2) = arcfatom(k)
                     end if
                     delx(i) = x(id(i)) - xr
                     dely(i) = y(id(i)) - yr
                     delz(i) = zgrid - z(id(i))
                     s2 = delx(i)**2 + dely(i)**2
                     r_s(i) = 1.0d0 / sqrt(s2)
                     r_s2(i) = r_s(i)**2
                     r(i) = sqrt(vdwrad(id(i))**2 - delz(i)**2)
                     r_r(i) = 1.0d0 / r(i)
                     u(i) = (ri**2+s2-r(i)**2) * (0.5d0*r_s(i) / ri)
                  end do
c
c     apply the chain rule repeatedly
c
                  theta1 = arci(k)
                  theta2 = arcf(k)
                  cos1 = cos(theta1)
                  cos2 = cos(theta2)
                  sin1 = sin(theta1)
                  sin2 = sin(theta2)
                  phi_xy = phi2 - phi1 - 0.5d0*(sin(2.0d0*phi2)
     &                                         -sin(2.0d0*phi1))
                  phi_z = sin(phi2)**2 - sin(phi1)**2
                  phi_xy = 0.5d0 * rrsq * phi_xy
                  phi_z = 0.5d0 * rrsq * phi_z
                  dfdtheta(1,1) = -cos1 * phi_xy
                  dfdtheta(2,1) = -sin1 * phi_xy
                  dfdtheta(3,1) = -phi_z
                  dfdtheta(1,2) =  cos2 * phi_xy
                  dfdtheta(2,2) =  sin2 * phi_xy
                  dfdtheta(3,2) =  phi_z
                  do i = 1, 2
                     dbetadx(i,1,0) = dely(i) * r_s2(i)
                     dbetadx(i,2,0) = -delx(i) * r_s2(i)
                     dbetadx(i,1,i) = -dbetadx(i,1,0)
                     dbetadx(i,2,i) = -dbetadx(i,2,0)
                  end do
                  do i = 1, 2
                     duds(i) = (1.0d0/ri) - (u(i)*r_s(i))
                     dsdx(i,1,i) = delx(i) * r_s(i)
                     dsdx(i,2,i) = dely(i) * r_s(i)
                     dsdx(i,1,0) = -dsdx(i,1,i)
                     dsdx(i,2,0) = -dsdx(i,2,i)
                     dudr(i) = -r(i) * r_s(i) / ri
                     drdz(i,i) = delz(i) * r_r(i)
                     drdz(i,0) = -drdz(i,i)
                  end do
                  do m = 0, 2
                     do i = 1, 2
                        dudx(i,1,m) = duds(i) * dsdx(i,1,m)
                        dudx(i,2,m) = duds(i) * dsdx(i,2,m)
                        dudx(i,3,m) = dudr(i) * drdz(i,m)
                     end do
                  end do
                  do i = 1, 2
                     u_term(i) = -1.0d0 / sqrt(1.0d0-u(i)**2)
                  end do
                  do j = 1, 3
                     do m = 0, 2
                        do i = 1, 2
                           dalphdx(i,j,m) = u_term(i) * dudx(i,j,m)
                        end do
                     end do
                  end do
                  do j = 1, 2
                     do m = 0, 2
                        dthetadx(1,j,m) = dbetadx(1,j,m)
     &                                       + dalphdx(1,j,m)
                        dthetadx(2,j,m) = dbetadx(2,j,m)
     &                                       - dalphdx(2,j,m)
                     end do
                  end do
                  do m = 0, 2
                     dthetadx(1,3,m) = dalphdx(1,3,m)
                     dthetadx(2,3,m) = -dalphdx(2,3,m)
                  end do
c
c     partials with respect to coordinates of serial atom id(m)
c
                  id(0) = iatom
                  do m = 0, 2
                     iblock = id(m)
                     do j = 1, 3
                        xhess(j,iblock) = xhess(j,iblock)
     &                     + dfdtheta(1,1) * dthetadx(1,j,m)
     &                     + dfdtheta(1,2) * dthetadx(2,j,m)
                        yhess(j,iblock) = yhess(j,iblock)
     &                     + dfdtheta(2,1) * dthetadx(1,j,m)
     &                     + dfdtheta(2,2) * dthetadx(2,j,m)
                        zhess(j,iblock) = zhess(j,iblock)
     &                     + dfdtheta(3,1) * dthetadx(1,j,m)
     &                     + dfdtheta(3,2) * dthetadx(2,j,m)
                     end do
                  end do
               end do
            end if
            zgrid = zgrid + zstep
            phiold = phi1
         end do
      end if
      return
      end
