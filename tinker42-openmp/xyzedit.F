 
C     $Id: xyzedit.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program xyzedit  --  editing of Cartesian coordinates  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "xyzedit" provides for modification and manipulation
c     of the contents of a Cartesian coordinates file
c
c
      program xyzedit
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'couple.i'
      include 'cutoff.i'
      include 'files.i'
      include 'iounit.i'
      include 'math.i'
      include 'refer.i'
      include 'titles.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,natom
      integer it,jt,ixyz
      integer nmode,mode
      integer offset,origin
      integer oldtype,newtype
      integer nlist,nmolecule
      integer freeunit
      integer list(maxatm)
      integer keep(maxatm)
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xran,yran,zran
      real*8 xbox,ybox,zbox
      real*8 xorig,yorig,zorig
      real*8 xcm,ycm,zcm
      real*8 phi,theta,psi
      real*8 cphi,ctheta,cpsi
      real*8 sphi,stheta,spsi
      real*8 dist2,cut2
      real*8 random,norm,weigh
      real*8 a(3,3)
      logical write
      character*120 coordfile
      character*120 record
      external merge
c
c
c     initialize various constants and the write flag
c
      call initial
      offset = 0
      write = .false.
c
c     read in the coordinates and force field definition
c
      call getxyz
      call field
      call katom
c
c     get the desired type of coordinate file modification
c
      write (iout,10)
   10 format (/,' The TINKER XYZ Edit Facility can Provide :',
     &        //,4x,'(1) Offset the Numbers of the Current Atoms',
     &        /,4x,'(2) Deletion of Individual Specified Atoms',
     &        /,4x,'(3) Deletion of Specified Types of Atoms',
     &        /,4x,'(4) Deletion of Atoms outside Cutoff Range',
     &        /,4x,'(5) Insertion of Individual Specified Atoms',
     &        /,4x,'(6) Replace Old Atom Type with a New Type',
     &        /,4x,'(7) Assign Connectivities based on Distance',
     &        /,4x,'(8) Convert Units from Bohrs to Angstroms',
     &        /,4x,'(9) Invert thru Origin to give Mirror Image',
     &        /,3x,'(10) Translate All Atoms by an X,Y,Z-Vector',
     &        /,3x,'(11) Translate Center of Mass to the Origin',
     &        /,3x,'(12) Translate a Specified Atom to the Origin',
     &        /,3x,'(13) Translate and Rotate to Inertial Frame',
     &        /,3x,'(14) Move to Specified Rigid Body Coordinates',
     &        /,3x,'(15) Move Stray Molecules into Periodic Box',
     &        /,3x,'(16) Create and Fill a Periodic Boundary Box',
     &        /,3x,'(17) Soak Current Molecule in Box of Solvent',
     &        /,3x,'(18) Append a Second XYZ file to Current One')
   20 continue
      nmode = 18
      mode = 0
      write (iout,30)
   30 format (/,' Enter the Number of the Desired Choice :  ',$)
      read (input,40,err=20)  mode
   40 format (i10)
      if (mode.lt.0 .or. mode.gt.nmode)  goto 20
c
c     get the offset value to be used in atom renumbering
c
      if (mode .eq. 1) then
   50    continue
         write (iout,60)
   60    format (/,' Offset used to Renumber the Current Atoms :  ',$)
         read (input,70,err=50)  offset
   70    format (i10)
         write = .true.
      end if
c
c     remove a specified list of individual atoms
c
      if (mode .eq. 2) then
         nlist = 0
         do i = 1, maxatm
            list(i) = 0
         end do
         write (iout,80)
   80    format (/,' Numbers of the Atoms to be Removed :  ',$)
         read (input,90)  record
   90    format (a120)
         read (record,*,err=100,end=100)  (list(i),i=1,maxatm)
  100    continue
         dowhile (list(nlist+1) .ne. 0)
            nlist = nlist + 1
         end do
         do i = 1, nlist
            if (list(i) .gt. n)  list(i) = n
            if (list(i) .lt. -n)  list(i) = -n
         end do
         call sort4 (nlist,list)
         do i = nlist, 1, -1
            if (i .gt. 1) then
               if (list(i-1) .lt. 0) then
                  do j = abs(list(i)), abs(list(i-1)), -1
                     call delete (j)
                  end do
               else if (list(i) .gt. 0) then
                  call delete (list(i))
               end if
            else if (list(i) .gt. 0) then
               call delete (list(i))
            end if
         end do
         write = .true.
         goto 20
      end if
c
c     remove all atoms with any of a specified list of atom types
c
      if (mode .eq. 3) then
         nlist = 0
         do i = 1, maxatm
            list(i) = 0
         end do
         write (iout,110)
  110    format (/,' Atom Types to be Removed :  ',$)
         read (input,120)  record
  120    format (a120)
         read (record,*,err=130,end=130)  (list(i),i=1,maxatm)
  130    continue
         dowhile (list(nlist+1) .ne. 0)
            nlist = nlist + 1
         end do
         natom = n
         do i = natom, 1, -1
            it = type(i)
            do j = 1, nlist
               if (list(j) .eq. it) then
                  call delete (i)
                  goto 140
               end if
            end do
  140       continue
         end do
         write = .true.
         goto 20
      end if
c
c     remove all atoms that are inactive and lie outside all cutoffs
c
      if (mode .eq. 4) then
         call active
         call cutoffs
         nlist = 0
         do i = 1, n
            keep(i) = 0
         end do
         cut2 = 0.0d0
         if (vdwcut .le. 1000.0d0)  cut2 = max(vdwcut**2,cut2)
         if (chgcut .le. 1000.0d0)  cut2 = max(chgcut**2,cut2)
         if (dplcut .le. 1000.0d0)  cut2 = max(dplcut**2,cut2)
         if (mpolecut .le. 1000.0d0)  cut2 = max(mpolecut**2,cut2)
         if (cut2 .eq. 0.0d0)  cut2 = 1.0d16
         do i = 1, n
            if (.not. use(i)) then
               do j = 1, n12(i)
                  keep(i12(j,i)) = i
               end do
               do j = 1, n13(i)
                  keep(i13(j,i)) = i
               end do
               do j = 1, n14(i)
                  keep(i14(j,i)) = i
               end do
               xi = x(i)
               yi = y(i)
               zi = z(i)
               do j = 1, n
                  if (use(j)) then
                     if (keep(j) .eq. i)  goto 150
                     dist2 = (x(j)-xi)**2 +(y(j)-yi)**2 +(z(j)-zi)**2
                     if (dist2 .le. cut2)  goto 150
                  end if
               end do
               nlist = nlist + 1
               list(nlist) = i
  150          continue
            end if
         end do
         do i = nlist, 1, -1
            call delete (list(i))
         end do
         write = .true.
         goto 20
      end if
c
c     insert a specified list of individual atoms
c
      if (mode .eq. 5) then
         nlist = 0
         do i = 1, maxatm
            list(i) = 0
         end do
         write (iout,160)
  160    format (/,' Numbers of the Atoms to be Inserted :  ',$)
         read (input,170)  record
  170    format (a120)
         read (record,*,err=180,end=180)  (list(i),i=1,maxatm)
  180    continue
         dowhile (list(nlist+1) .ne. 0)
            nlist = nlist + 1
         end do
         call sort4 (nlist,list)
         do i = nlist, 1, -1
            if (i .gt. 1) then
               if (list(i-1) .lt. 0) then
                  do j = abs(list(i-1)), abs(list(i))
                     call insert (j)
                  end do
               else if (list(i) .gt. 0) then
                  call insert (list(i))
               end if
            else if (list(i) .gt. 0) then
               call insert (list(i))
            end if
         end do
         write = .true.
         goto 20
      end if
c
c     get an old atom type and new atom type for replacement
c
      if (mode .eq. 6) then
  190    continue
         oldtype = 0
         newtype = 0
         write (iout,200)
  200    format (/,' Numbers of the Old and New Atom Types :  ',$)
         read (input,210)  record
  210    format (a120)
         read (record,*,err=190)  oldtype,newtype
         do i = 1, n
            if (type(i) .eq. oldtype)  type(i) = newtype
         end do
         call katom
         write = .true.
         goto 20
      end if
c
c     assign atom connectivities based on interatomic distances
c
      if (mode .eq. 7) then
         call unitcell
         call lattice
         do i = 1, n
            n12(i) = 0
         end do
         do i = 1, n-1
            xi = x(i)
            yi = y(i)
            zi = z(i)
            it = atomic(i)
            do j = i+1, n
               jt = atomic(j)
               xr = x(j) - xi
               yr = y(j) - yi
               zr = z(j) - zi
               if (use_image)  call image (xr,yr,zr,0)
               dist2 = xr*xr + yr*yr + zr*zr
               if (it.ne.1 .and. jt.ne.1) then
                  if (dist2 .lt. 4.0d0) then
                     n12(i) = n12(i) + 1
                     i12(n12(i),i) = j
                     n12(j) = n12(j) + 1
                     i12(n12(j),j) = i
                  end if
               else if (it.ne.1 .or. jt.ne.1) then
                  if (dist2 .lt. 2.25d0) then
                     n12(i) = n12(i) + 1
                     i12(n12(i),i) = j
                     n12(j) = n12(j) + 1
                     i12(n12(j),j) = i
                  end if
               end if
            end do
         end do
         do i = 1, n
            call sort (n12(i),i12(1,i))
         end do
         write = .true.
         goto 20
      end if
c
c     convert the coordinate units from Bohrs to Angstroms
c
      if (mode .eq. 8) then
         do i = 1, n
            x(i) = x(i) * bohr
            y(i) = y(i) * bohr
            z(i) = z(i) * bohr
         end do
         write = .true.
         goto 20
      end if
c
c     get mirror image by inverting coordinates through origin
c
      if (mode .eq. 9) then
         do i = 1, n
            x(i) = -x(i)
            y(i) = -y(i)
            z(i) = -z(i)
         end do
         write = .true.
         goto 20
      end if
c
c     translate the entire system by a specified x,y,z-vector
c
      if (mode .eq. 10) then
         xr = 0.0d0
         yr = 0.0d0
         zr = 0.0d0
         write (iout,220)
  220    format (/,' Enter Translation Vector Components :  ',$)
         read (input,230)  record
  230    format (a120)
         read (record,*,err=240,end=240)  xr,yr,zr
  240    continue
         do i = 1, n
            x(i) = x(i) + xr
            y(i) = y(i) + yr
            z(i) = z(i) + zr
         end do
         write = .true.
         goto 20
      end if
c
c     translate the center of mass to the coordinate origin
c
      if (mode .eq. 11) then
         xcm = 0.0d0
         ycm = 0.0d0
         zcm = 0.0d0
         norm = 0.0d0
         do i = 1, n
            weigh = mass(i)
            xcm = xcm + x(i)*weigh
            ycm = ycm + y(i)*weigh
            zcm = zcm + z(i)*weigh
            norm = norm + weigh
         end do
         xcm = xcm / norm
         ycm = ycm / norm
         zcm = zcm / norm
         do i = 1, n
            x(i) = x(i) - xcm
            y(i) = y(i) - ycm
            z(i) = z(i) - zcm
         end do
         write = .true.
         goto 20
      end if
c
c     translate to place a specified atom at the origin
c
      if (mode .eq. 12) then
         write (iout,250)
  250    format (/,' Number of the Atom to Move to the Origin :  ',$)
         read (input,260)  origin
  260    format (i10)
         xorig = x(origin)
         yorig = y(origin)
         zorig = z(origin)
         do i = 1, n
            x(i) = x(i) - xorig
            y(i) = y(i) - yorig
            z(i) = z(i) - zorig
         end do
         write = .true.
         goto 20
      end if
c
c     translate and rotate into standard orientation
c
      if (mode .eq. 13) then
         call inertia (2)
         write = .true.
         goto 20
      end if
c
c     translate and rotate to specified rigid body coordinates
c
      if (mode .eq. 14) then
         xcm = 0.0d0
         ycm = 0.0d0
         zcm = 0.0d0
         phi = 0.0d0
         theta = 0.0d0
         psi = 0.0d0
         write (iout,270)
  270    format (/,' Enter Rigid Body Coordinates :  ',$)
         read (input,280)  record
  280    format (a120)
         read (record,*,err=290,end=290)  xcm,ycm,zcm,phi,theta,psi
  290    continue
         call inertia (2)
         phi = phi / radian
         theta = theta / radian
         psi = psi / radian
         cphi = cos(phi)
         sphi = sin(phi)
         ctheta = cos(theta)
         stheta = sin(theta)
         cpsi = cos(psi)
         spsi = sin(psi)
         a(1,1) = ctheta * cphi
         a(2,1) = spsi*stheta*cphi - cpsi*sphi
         a(3,1) = cpsi*stheta*cphi + spsi*sphi
         a(1,2) = ctheta * sphi
         a(2,2) = spsi*stheta*sphi + cpsi*cphi
         a(3,2) = cpsi*stheta*sphi - spsi*cphi
         a(1,3) = -stheta
         a(2,3) = ctheta * spsi
         a(3,3) = ctheta * cpsi
         do i = 1, n
            xorig = x(i)
            yorig = y(i)
            zorig = z(i)
            x(i) = a(1,1)*xorig + a(2,1)*yorig + a(3,1)*zorig + xcm
            y(i) = a(1,2)*xorig + a(2,2)*yorig + a(3,2)*zorig + ycm
            z(i) = a(1,3)*xorig + a(2,3)*yorig + a(3,3)*zorig + zcm
         end do
         write = .true.
         goto 20
      end if
c
c     move stray molecules back into original periodic box
c
      if (mode .eq. 15) then
         call unitcell
         if (use_bounds) then
            call lattice
            call molecule
            call bounds
            write = .true.
         end if
         goto 20
      end if
c
c     create a random box full of the current coordinates file
c
      if (mode .eq. 16) then
         write (iout,300)
  300    format (/,' Enter Number of Molecules in Box :  ',$)
         read (input,310)  nmolecule
  310    format (i10)
         xbox = 0.0d0
         ybox = 0.0d0
         zbox = 0.0d0
         dowhile (xbox .eq. 0.0d0)
            write (iout,320)
  320       format (/,' Enter Periodic Box Dimensions (X,Y,Z) :  ',$)
            read (input,330)  record
  330       format (a120)
            read (record,*,err=340,end=340)  xbox,ybox,zbox
  340       continue
            if (ybox .eq. 0.0d0)  ybox = xbox
            if (zbox .eq. 0.0d0)  zbox = xbox
         end do
         ixyz = freeunit ()
         coordfile = filename(1:leng)//'.xyz'
         call version (coordfile,'new')
         open (unit=ixyz,file=coordfile,status='new')
         if (ltitle .eq. 0) then
            write (ixyz,350)  n*nmolecule
  350       format (i6)
         else
            write (ixyz,360)  n*nmolecule,title(1:ltitle)
  360       format (i6,2x,a)
         end if
         do k = 1, nmolecule
            offset = (k-1) * n
            xran = xbox * random ()
            yran = ybox * random ()
            zran = zbox * random ()
            do i = 1, n
               write (ixyz,370)  i+offset,name(i),x(i)+xran,
     &                           y(i)+yran,z(i)+zran,type(i),
     &                           (i12(j,i)+offset,j=1,n12(i))
  370          format (i6,2x,a3,3f12.6,5i6)
            end do
         end do
         close (unit=ixyz)
         write (iout,380)  coordfile
  380    format (/,' New Coordinates written to :  ',a)
         write = .false.
      end if
c
c     solvate the current system by insertion into a solvent box
c
      if (mode .eq. 17) then
         call soak
         write = .true.
      end if
c
c     append a second file to the current coordinates file
c
      if (mode .eq. 18) then
         call makeref
         call getxyz
         call merge
         write = .true.
         goto 20
      end if
c
c     write out the new coordinates file in its new format
c
      if (write) then
         ixyz = freeunit ()
         coordfile = filename(1:leng)//'.xyz'
         call version (coordfile,'new')
         open (unit=ixyz,file=coordfile,status='new')
         if (offset .eq. 0) then
            call prtxyz (ixyz)
         else
            if (ltitle .eq. 0) then
               write (ixyz,390)  n
  390          format (i6)
            else
               write (ixyz,400)  n,title(1:ltitle)
  400          format (i6,2x,a)
            end if
            do i = 1, n
               write (ixyz,410)  i+offset,name(i),x(i),y(i),z(i),
     &                           type(i),(i12(j,i)+offset,j=1,n12(i))
  410          format (i6,2x,a3,3f12.6,5i6)
            end do
         end if
         close (unit=ixyz)
         write (iout,420)  coordfile
  420    format (/,' New Coordinates written to File :  ',a)
      end if
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine soak  --  place a solute into a solvent box  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "soak" takes a currently defined solute system and places
c     it into a solvent box, with removal of any solvent molecules
c     that overlap the solute
c
c
      subroutine soak
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'iounit.i'
      include 'molcul.i'
      include 'refer.i'
      integer i,k,isolv
      integer ntot,freeunit
      real*8 xi,yi,zi
      real*8 xr,yr,zr,rik2
      real*8 close,close2
      logical remove(maxatm)
      character*120 solvfile
      external merge
c
c
c     make a copy of the solute coordinates and connectivities
c
      call makeref
c
c     read the coordinates for the solvent box
c
   10 continue
      write (iout,20)
   20 format (/,' Enter Name of Solvent Box Coordinates :  ',$)
      read (input,30)  solvfile
   30 format (a120)
      call suffix (solvfile,'xyz')
      call version (solvfile,'old')
      isolv = freeunit ()
      open (unit=isolv,file=solvfile,status='old',err=10)
      rewind (unit=isolv)
      call readxyz (isolv)
      close (unit=isolv)
c
c     combine solute and solvent into a single coordinate set
c
      call merge
c
c     count number of molecules and set lattice parameters
c
      call molecule
      call unitcell
      call lattice
c
c     initialize the list of solvent molecules to be deleted
c
      do i = 1, nmol
         remove(i) = .false.
      end do
c
c     search for close contacts between solute and solvent
c
      close = 1.5d0
      close2 = close * close
      do i = 1, nref
         xi = x(i)
         yi = y(i)
         zi = z(i)
         do k = nref+1, n
            xr = x(k) - xi
            yr = y(k) - yi
            zr = z(k) - zi
            if (use_image)  call image (xr,yr,zr,0)
            rik2 = xr*xr + yr*yr + zr*zr
            if (rik2 .lt. close2)  remove(molcule(k)) = .true.
         end do
      end do
c
c     remove solvent molecules that are too close to the solute
c
      ntot = n
      do i = ntot, nref+1, -1
         if (remove(molcule(i)))  call delete (i)
      end do
      return
      end
