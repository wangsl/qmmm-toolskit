 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program vibrate  --  vibrational analysis and normal modes  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "vibrate" performs a vibrational normal mode analysis; the
c     Hessian matrix of second derivatives is determined and then
c     diagonalized both directly and after mass weighting; output
c     consists of the eigenvalues of the force constant matrix as
c     well as the vibrational frequencies and displacements
c
c
      program vibrate
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'files.i'
      include 'hescut.i'
      include 'iounit.i'
      include 'math.i'
      include 'units.i'
      integer i,j,k,m
      integer ixyz,ihess
      integer lext,freeunit
      integer nfreq,ndummy
      integer nvib,ivib,nview
      integer hinit(3,maxatm)
      integer hstop(3,maxatm)
      integer hindex(maxhess)
      real*8 factor,vnorm,ratio
      real*8 h(maxhess),hdiag(3,maxatm)
      real*8 xref(maxatm),yref(maxatm)
      real*8 zref(maxatm)
      real*8 mass2(maxatm)
      real*8 matrix((maxvib+1)*maxvib/2)
      real*8 eigen(maxvib)
      real*8 vects(maxvib,maxvib)
      real*8 a(maxvib+1),b(maxvib+1)
      real*8 p(maxvib+1),w(maxvib+1)
      real*8 ta(maxvib+1),tb(maxvib+1)
      real*8 ty(maxvib+1)
      logical exist,query
      character*1 eigunits(maxvib)
      character*7 ext
      character*120 coordfile
      character*120 string
      character*120 record
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     initialize various things needed for vibrations
c
      nfreq = 3 * n
      ndummy = 0
      if (nfreq.gt.maxvib .or. nfreq**2.gt.maxhess) then
         write (iout,10)
   10    format (/,' VIBRATE  --  Too many Atoms in the Molecule')
         call fatal
      end if
      do i = 1, n
         if (atomic(i) .eq. 0) then
            ndummy = ndummy + 1
            mass(i) = 0.001d0
         end if
         mass2(i) = sqrt(mass(i))
      end do
      nvib = nfreq - 3*ndummy
c
c     calculate the Hessian matrix of second derivatives
c
      hesscut = 0.0d0
      call hessian (h,hinit,hstop,hindex,hdiag)
c
c     store upper triangle of the Hessian in "matrix"
c
      ihess = 0
      do i = 1, n
         do j = 1, 3
            ihess = ihess + 1
            matrix(ihess) = hdiag(j,i)
            do k = hinit(j,i), hstop(j,i)
               ihess = ihess + 1
               matrix(ihess) = h(k)
            end do
         end do
      end do
c
c     perform diagonalization to get Hessian eigenvalues
c
      call diagq (nfreq,maxvib,nfreq,matrix,eigen,vects,
     &                     a,b,p,w,ta,tb,ty)
      write (iout,20)
   20 format (/,' Eigenvalues of the Hessian Matrix :',/)
      write (iout,30)  (i,eigen(i),i=1,nvib)
   30 format (5(i5,f9.3,1x))
      write (iout,40)
   40 format ()
c
c     store upper triangle of the mass-weighted Hessian matrix
c
      ihess = 0
      do i = 1, n
         do j = 1, 3
            ihess = ihess + 1
            matrix(ihess) = hdiag(j,i) / mass(i)
            do k = hinit(j,i), hstop(j,i)
               m = (hindex(k)+2) / 3
               ihess = ihess + 1
               matrix(ihess) = h(k) / (mass2(i)*mass2(m))
            end do
         end do
      end do
c
c     diagonalize to get vibrational frequencies and normal modes
c
      call diagq (nfreq,maxvib,nfreq,matrix,eigen,vects,
     &                     a,b,p,w,ta,tb,ty)
      factor = sqrt(convert) / (2.0d0*pi*lightspd)
      do i = 1, nvib
         if (eigen(i) .lt. 0) then
            eigunits(i) = 'I'
         else
            eigunits(i) = ' '
         end if
         eigen(i) = factor * sqrt(abs(eigen(i)))
      end do
      write (iout,50)
   50 format (/,' Vibrational Frequencies (cm-1) :',/)
      write (iout,60)  (i,eigen(i),eigunits(i),i=1,nvib)
   60 format (5(i5,f9.3,a1))
c
c     form Cartesian coordinate displacements from normal modes
c
      do i = 1, nvib
         vnorm = 0.0d0
         do j = 1, nfreq
            k = (j+2) / 3
            vects(j,i) = vects(j,i) / mass2(k)
            vnorm = vnorm + vects(j,i)**2
         end do
         vnorm = sqrt(vnorm)
         do j = 1, nfreq
            vects(j,i) = vects(j,i) / vnorm
         end do
      end do
c
c     decide upon the vibrations to be output to a disk file
c
      query = .true.
   70 continue
      ivib = 0
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=80,end=80)  ivib
         query = .false.
      end if
   80 continue
      if (query) then
         write (iout,90)
   90    format (/,' Enter the Number of the Vibration',
     &              ' to be Output :  ',$)
         read (input,100)  record
  100    format (a120)
         read (record,*,err=110,end=110)  ivib
      end if
  110 continue
c
c     print the vibrational frequency and normal mode
c
      if (ivib.ge.1 .and. ivib.le.nvib) then
         write (iout,120)  ivib,eigen(ivib),eigunits(ivib)
  120    format (/,' Vibrational Normal Mode',i5,' with Frequency',
     &              f10.2,a1,' cm-1',
     &           //,5x,'Atom',5x,'Delta X',5x,'Delta Y',5x,'Delta Z',/)
         do i = 1, n
            j = 3 * (i-1)
            write (iout,130)  i,vects(j+1,ivib),vects(j+2,ivib),
     &                        vects(j+3,ivib)
  130       format (4x,i5,3f12.6)
         end do
c
c     create a name for the vibrational displacement file
c
         lext = 3
         call numeral (ivib,ext,lext)
         coordfile = filename(1:leng)//'.'//ext(1:lext)
         ixyz = freeunit ()
         call version (coordfile,'new')
         open (unit=ixyz,file=coordfile,status='new')
c
c     store the original atomic coordinates
c
         do i = 1, n
            xref(i) = x(i)
            yref(i) = y(i)
            zref(i) = z(i)
         end do
c
c     make file with plus and minus the current vibration
c
         nview = 3
         do i = -nview, nview
            ratio = i / dble(nview)
            do k = 1, n
               j = 3 * (k-1)
               x(k) = xref(k) + ratio*vects(j+1,ivib)
               y(k) = yref(k) + ratio*vects(j+2,ivib)
               z(k) = zref(k) + ratio*vects(j+3,ivib)
            end do
            call prtxyz (ixyz)
         end do
         close (unit=ixyz)
c
c     restore the original coordinates; get next vibration
c
         do i = 1, n
            x(i) = xref(i)
            y(i) = yref(i)
            z(i) = zref(i)
         end do
         goto 70
      end if
c
c     perform any final tasks before program exit
c
      call final
      end
