 
C     $Id$
c
c
c     ################################################################
c     ##  COPYRIGHT (C) 1990 by Craig Kundrot & Jay William Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine volume1  --  Cartesian excluded volume derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "volume1" calculates first derivatives of the total excluded
c     volume with respect to the Cartesian coordinates of each atom
c
c     literature reference:
c
c     C. E. Kundrot, J. W. Ponder and F. M. Richards, "Algorithms for
c     Calculating Excluded Volume and Its Derivatives as a Function
c     of Molecular Conformation and Their Use in Energy Minimization",
c     Journal of Computational Chemistry, 12, 402-409 (1991)
c
c
      subroutine volume1 (radius,probe,dex)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'math.i'
      integer maxcube,maxarc
      parameter (maxcube=15)
      parameter (maxarc=300)
      integer i,j,k,m
      integer io,ir,in
      integer narc,nx,ny,nz
      integer istart,istop
      integer jstart,jstop
      integer kstart,kstop
      integer mstart,mstop
      integer isum,icube,itemp
      integer itab(maxatm)
      integer inov(maxarc)
      integer cube(2,maxcube,maxcube,maxcube)
      real*8 xmin,ymin,zmin
      real*8 xmax,ymax,zmax
      real*8 aa,bb,temp,phi_term
      real*8 theta1,theta2,dtheta
      real*8 seg_dx,seg_dy,seg_dz
      real*8 pre_dx,pre_dy,pre_dz
      real*8 rinsq,rdiff
      real*8 rsecn,rsec2n
      real*8 cosine,ti,tf
      real*8 alpha,beta
      real*8 ztop,zstart
      real*8 ztopshave
      real*8 phi1,cos_phi1
      real*8 phi2,cos_phi2
      real*8 zgrid,pix2
      real*8 rsec2r,rsecr
      real*8 rr,rrx2,rrsq
      real*8 rmax,edge
      real*8 xr,yr,zr
      real*8 dist2,vdwsum
      real*8 probe,zstep
      real*8 dex(3,maxatm)
      real*8 radius(maxatm)
      real*8 vdwrad(maxatm)
      real*8 arci(maxarc)
      real*8 arcf(maxarc)
      real*8 dx(maxarc),dy(maxarc)
      real*8 dsq(maxarc),d(maxarc)
      logical skip(maxatm)
c
c
c     fix the stepsize in the z-direction; this value sets
c     the accuracy of the numerical derivatives; zstep=0.06
c     is a good balance between compute time and accuracy
c
      zstep = 0.0601d0
c
c     initialize minimum and maximum ranges of atoms
c
      pix2 = 2.0d0 * pi
      rmax = 0.0d0
      xmin = x(1)
      xmax = x(1)
      ymin = y(1)
      ymax = y(1)
      zmin = z(1)
      zmax = z(1)
c
c     assign van der Waals radii to the atoms; note that
c     the radii are incremented by the size of the probe;
c     then get the maximum and minimum ranges of atoms
c
      do i = 1, n
         vdwrad(i) = radius(i)
         if (vdwrad(i) .eq. 0.0d0) then
            skip(i) = .true.
         else
            skip(i) = .false.
            vdwrad(i) = vdwrad(i) + probe
            if (vdwrad(i) .gt. rmax)  rmax = vdwrad(i)
            if (x(i) .lt. xmin)  xmin = x(i)
            if (x(i) .gt. xmax)  xmax = x(i)
            if (y(i) .lt. ymin)  ymin = y(i)
            if (y(i) .gt. ymax)  ymax = y(i)
            if (z(i) .lt. zmin)  zmin = z(i)
            if (z(i) .gt. zmax)  zmax = z(i)
         end if
      end do
c
c     load the cubes based on coarse lattice; first of all
c     set edge length to the maximum diameter of any atom
c
      edge = 2.0d0 * rmax
      nx = int((xmax-xmin)/edge) + 1
      ny = int((ymax-ymin)/edge) + 1
      nz = int((zmax-zmin)/edge) + 1
      if (max(nx,ny,nz) .gt. maxcube) then
         write (iout,10)
   10    format (/,' VOLUME1  --  Increase the Value of MAXCUBE')
         call fatal
      end if
c
c     initialize the coarse lattice of cubes
c
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz
               cube(1,i,j,k) = 0
               cube(2,i,j,k) = 0
            end do
         end do
      end do
c
c     find the number of atoms in each cube
c
      do m = 1, n
         if (.not. skip(m)) then
            i = int((x(m)-xmin)/edge) + 1
            j = int((y(m)-ymin)/edge) + 1
            k = int((z(m)-zmin)/edge) + 1
            cube(1,i,j,k) = cube(1,i,j,k) + 1
         end if
      end do
c
c     determine the highest index in the array "itab" for the
c     atoms that fall into each cube; the first cube that has
c     atoms defines the first index for "itab"; the final index
c     for the atoms in the present cube is the final index of
c     the last cube plus the number of atoms in the present cube
c
      isum = 0
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz
               icube = cube(1,i,j,k)
               if (icube .ne. 0) then
                  isum = isum + icube
                  cube(2,i,j,k) = isum
               end if
            end do
         end do
      end do
c
c     "cube(2,,,)" now contains a pointer to the array "itab"
c     giving the position of the last entry for the list of
c     atoms in that cube of total number equal to "cube(1,,,)"
c
      do m = 1, n
         if (.not. skip(m)) then
            i = int((x(m)-xmin)/edge) + 1
            j = int((y(m)-ymin)/edge) + 1
            k = int((z(m)-zmin)/edge) + 1
            icube = cube(2,i,j,k)
            itab(icube) = m
            cube(2,i,j,k) = icube - 1
         end if
      end do
c
c     set "cube(2,,,)" to be the starting index in "itab"
c     for atom list of that cube; and "cube(1,,,)" to be
c     the stop index
c
      isum = 0
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz
               icube = cube(1,i,j,k)
               if (icube .ne. 0) then
                  isum = isum + icube
                  cube(1,i,j,k) = isum
                  cube(2,i,j,k) = cube(2,i,j,k) + 1
               end if
            end do
         end do
      end do
c
c     process in turn each atom from the coordinate list;
c     first select the potential intersecting atoms
c
      do ir = 1, n
         pre_dx = 0.0d0
         pre_dy = 0.0d0
         pre_dz = 0.0d0
         if (skip(ir))  goto 50
         rr = vdwrad(ir)
         rrx2 = 2.0d0 * rr
         rrsq = rr * rr
         xr = x(ir)
         yr = y(ir)
         zr = z(ir)
c
c     find cubes to search for overlaps of current atom
c
         istart = int((xr-xmin)/edge)
         istop = min(istart+2,nx)
         istart = max(istart,1)
         jstart = int((yr-ymin)/edge)
         jstop = min(jstart+2,ny)
         jstart = max(jstart,1)
         kstart = int((zr-zmin)/edge)
         kstop = min(kstart+2,nz)
         kstart = max(kstart,1)
c
c     load all overlapping atoms into "inov"
c
         io = 0
         do i = istart, istop
            do j = jstart, jstop
               do k = kstart, kstop
                  mstart = cube(2,i,j,k)
                  if (mstart .ne. 0) then
                     mstop = cube(1,i,j,k)
                     do m = mstart, mstop
                        in = itab(m)
                        if (in .ne. ir) then
                           io = io + 1
                           if (io .gt. maxarc) then
                              write (iout,20)
   20                         format (/,' VOLUME1  --  Increase ',
     &                                   ' the Value of MAXARC')
                              call fatal
                           end if
                           dx(io) = x(in) - xr
                           dy(io) = y(in) - yr
                           dsq(io) = dx(io)**2 + dy(io)**2
                           dist2 = dsq(io) + (z(in)-zr)**2
                           vdwsum = (rr+vdwrad(in))**2
                           if (dist2.gt.vdwsum .or. dist2.eq.0.0d0) then
                              io = io - 1
                           else
                              d(io) = sqrt(dsq(io))
                              inov(io) = in
                           end if
                        end if
                     end do
                  end if
               end do
            end do
         end do
c
c     determine resolution along the z-axis
c
         if (io .ne. 0) then
            ztop = zr + rr
            ztopshave = ztop - zstep
            zgrid = zr - rr
c
c     half of the part not covered by the planes
c
            zgrid = zgrid + 0.5d0*(rrx2-(int(rrx2/zstep)*zstep))
            zstart = zgrid
c
c     section atom spheres perpendicular to the z axis
c
            dowhile (zgrid .le. ztop)
c
c     "rsecr" is radius of circle of intersection
c     of "ir" sphere on the current sphere
c
               rsec2r = rrsq - (zgrid-zr)**2
               if (rsec2r .lt. 0.0d0)  rsec2r = 0.000001d0
               rsecr = sqrt(rsec2r)
               if (zgrid .ge. ztopshave) then
                  cos_phi1 = 1.0d0
                  phi1 = 0.0d0
               else
                  cos_phi1 = (zgrid + 0.5d0*zstep - zr) / rr
                  phi1 = acos(cos_phi1)
               end if
               if (zgrid .eq. zstart) then
                  cos_phi2 = -1.0d0
                  phi2 = pi
               else
                  cos_phi2 = (zgrid - 0.5d0*zstep - zr) / rr
                  phi2 = acos(cos_phi2)
               end if
c
c     check intersections of neighbor circles
c
               narc = 0
               do k = 1, io
                  in = inov(k)
                  rinsq = vdwrad(in)**2
                  rsec2n = rinsq - (zgrid-z(in))**2
                  if (rsec2n .gt. 0.0d0) then
                     rsecn = sqrt(rsec2n)
                     if (d(k) .lt. rsecr+rsecn) then
                        rdiff = rsecr - rsecn
                        if (d(k) .le. abs(rdiff)) then
                           if (rdiff .lt. 0.0d0) then
                              narc = 1
                              arci(narc) = 0.0d0
                              arcf(narc) = pix2
                           end if
                           goto 40
                        end if
                        narc = narc + 1
                        if (narc .gt. maxarc) then
                           write (iout,30)
   30                      format (/,' VOLUME1  --  Increase',
     &                                ' the Value of MAXARC')
                           call fatal
                        end if
c
c     initial and final arc endpoints are found for intersection
c     of "ir" circle with another circle contained in same plane;
c     the initial endpoint of the enclosed arc is stored in "arci",
c     the final endpoint in "arcf"; get "cosine" via law of cosines
c
                        cosine = (dsq(k)+rsec2r-rsec2n)
     &                                     / (2.0d0*d(k)*rsecr)
                        cosine = min(1.0d0,max(-1.0d0,cosine))
c
c     "alpha" is the angle between a line containing either point
c     of intersection and the reference circle center and the
c     line containing both circle centers; "beta" is the angle
c     between the line containing both circle centers and x-axis
c
                        alpha = acos(cosine)
                        beta = atan2(dy(k),dx(k))
                        if (dy(k) .lt. 0.0d0)  beta = beta + pix2
                        ti = beta - alpha
                        tf = beta + alpha
                        if (ti .lt. 0.0d0)  ti = ti + pix2
                        if (tf .gt. pix2)  tf = tf - pix2
                        arci(narc) = ti
c
c     if the arc crosses zero, then it is broken into two segments;
c     the first ends at two pi and the second begins at zero
c
                        if (tf .lt. ti) then
                           arcf(narc) = pix2
                           narc = narc + 1
                           arci(narc) = 0.0d0
                        end if
                        arcf(narc) = tf
   40                   continue
                     end if
                  end if
               end do
c
c     find the pre-area and pre-forces on this section (band),
c     "pre-" means a multiplicative factor is yet to be applied
c
               if (narc .eq. 0) then
                  seg_dz = pix2 * (cos_phi1**2 - cos_phi2**2)
                  pre_dz = pre_dz + seg_dz
               else
c
c     sort the arc endpoint arrays, each with "narc" entries,
c     in order of increasing values of the arguments in "arci"
c
                  k = 1
                  dowhile (k .lt. narc)
                     aa = arci(k)
                     bb = arcf(k)
                     temp = 1000000.0d0
                     do i = k, narc
                        if (arci(i) .le. temp) then
                           temp = arci(i)
                           itemp = i
                        end if
                     end do
                     arci(k) = arci(itemp)
                     arcf(k) = arcf(itemp)
                     arci(itemp) = aa
                     arcf(itemp) = bb
                     k = k + 1
                  end do
c
c     consolidate arcs by removing overlapping arc endpoints
c
                  temp = arcf(1)
                  j = 1
                  do k = 2, narc
                     if (temp .lt. arci(k)) then
                        arcf(j) = temp
                        j = j + 1
                        arci(j) = arci(k)
                        temp = arcf(k)
                     else if (temp .lt. arcf(k)) then
                        temp = arcf(k)
                     end if
                  end do
                  arcf(j) = temp
                  narc = j
                  if (narc .eq. 1) then
                     narc = 2
                     arcf(2) = pix2
                     arci(2) = arcf(1)
                     arcf(1) = arci(1)
                     arci(1) = 0.0d0
                  else
                     temp = arci(1)
                     do k = 1, narc-1
                        arci(k) = arcf(k)
                        arcf(k) = arci(k+1)
                     end do
                     if (temp.eq.0.0d0 .and. arcf(narc).eq.pix2) then
                        narc = narc - 1
                     else
                        arci(narc) = arcf(narc)
                        arcf(narc) = temp
                     end if
                  end if
c
c     compute the numerical pre-derivative values
c
                  do k = 1, narc
                     theta1 = arci(k)
                     theta2 = arcf(k)
                     if (theta2 .ge. theta1) then
                        dtheta = theta2 - theta1
                     else
                        dtheta = (theta2+pix2) - theta1
                     end if
                     phi_term = phi2 - phi1 - 0.5d0*(sin(2.0d0*phi2)
     &                                              -sin(2.0d0*phi1))
                     seg_dx = (sin(theta2)-sin(theta1)) * phi_term
                     seg_dy = (cos(theta1)-cos(theta2)) * phi_term
                     seg_dz = dtheta * (cos_phi1**2 - cos_phi2**2)
                     pre_dx = pre_dx + seg_dx
                     pre_dy = pre_dy + seg_dy
                     pre_dz = pre_dz + seg_dz
                  end do
               end if
               zgrid = zgrid + zstep
            end do
         end if
   50    continue
         dex(1,ir) = 0.5d0 * rrsq * pre_dx
         dex(2,ir) = 0.5d0 * rrsq * pre_dy
         dex(3,ir) = 0.5d0 * rrsq * pre_dz
      end do
      return
      end
