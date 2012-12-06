 
C     $Id: kewald.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine kewald  --  Ewald sum parameter assignment  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "kewald" assigns both regular Ewald summation and particle
c     mesh Ewald parameters for a periodic box
c
c
      subroutine kewald
      implicit none
      include 'sizes.i'
      include 'boxes.i'
      include 'cutoff.i'
      include 'ewald.i'
      include 'ewreg.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'pme.i'
      include 'potent.i'
      integer maxpower
      parameter (maxpower=50)
      integer i,k
      integer next,ntable
      integer ifft1,ifft2,ifft3
      integer jmax,kmax,lmax
      integer multi(maxpower)
      real*8 delta,eps
      logical use_pme
      character*20 keyword
      character*120 record
      character*120 string
      data multi  /   2,   3,   4,   5,   6,   8,   9,  10,  12,  15,
     &               16,  18,  20,  24,  25,  27,  30,  32,  36,  40,
     &               45,  48,  50,  54,  60,  64,  72,  75,  80,  81,
     &               90,  96, 100, 108, 120, 125, 128, 135, 144, 150,
     &              160, 162, 180, 192, 200, 216, 225, 240, 243, 250 /
c
c
c     decide between regular Ewald and particle mesh Ewald
c
      use_pme = .false.
      if (use_charge)  use_pme = .true.
c
c     set defaults for reciprocal fraction and boundary treatment
c
      frecip = 0.5d0
      tinfoil = .true.
c
c     estimate an optimal value for the Ewald coefficient
c
      eps = 1.0d-8
      call ewaldcof (aewald,ewaldcut,eps)
c
c     set defaults for PME B-spline order and charge grid size;
c     grid is a product of powers of 2, 3 and/or 5 for efficiency
c
      if (use_pme) then
         bsorder = 8
         delta = 1.0d-8
         nfft1 = 0
         nfft2 = 0
         nfft3 = 0
         ifft1 = int(1.5d0*xbox-delta) + 1
         ifft2 = int(1.5d0*ybox-delta) + 1
         ifft3 = int(1.5d0*zbox-delta) + 1
         do i = maxpower, 1, -1
            k = multi(i)
            if (k .le. maxfft) then
               if (k.ge.ifft1 .or. nfft1.eq.0)  nfft1 = k
               if (k.ge.ifft2 .or. nfft2.eq.0)  nfft2 = k
               if (k.ge.ifft3 .or. nfft3.eq.0)  nfft3 = k
            end if
         end do
      end if
c
c     search keywords for Ewald summation commands
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:12) .eq. 'EWALD-ALPHA ') then
            read (string,*,err=20)  aewald
         else if (keyword(1:15) .eq. 'EWALD-FRACTION ') then
            read (string,*,err=20)  frecip
         else if (keyword(1:15) .eq. 'EWALD-BOUNDARY ') then
            tinfoil = .false.
         else if (keyword(1:9) .eq. 'PME-GRID ') then
            nfft1 = 0
            nfft2 = 0
            nfft3 = 0
            read (string,*,err=10,end=10)  nfft1,nfft2,nfft3
   10       continue
            if (nfft2 .eq. 0)  nfft2 = nfft1
            if (nfft3 .eq. 0)  nfft3 = nfft1
         else if (keyword(1:10) .eq. 'PME-ORDER ') then
            read (string,*,err=20)  bsorder
         end if
   20    continue
      end do
c
c     check the B-spline order and charge grid dimension for PME
c
      if (use_pme) then
         if (bsorder .gt. maxorder) then
            write (iout,30)
   30       format (/,' KEWALD  --  B-Spline Order Too Large;',
     &                 ' Increase MAXORDER')
            call fatal
         end if
         if (max(nfft1,nfft2,nfft3) .gt. maxfft) then
            write (iout,40)
   40       format (/,' KEWALD  --  FFT Charge Grid Too Large;',
     &                 ' Increase MAXFFT')
            call fatal
         else if (nfft1.lt.ifft1 .or. nfft2.lt.ifft2
     &                .or. nfft3.lt.ifft3) then
            write (iout,50)
   50       format (/,' KEWALD  --  Warning, Small Charge Grid',
     &                 ' may give Poor Accuracy')
         end if
c
c     check the number of k-vectors for regular Ewald summation
c
      else
         jmax = int(frecip/recip(1,1))
         kmax = int(frecip/recip(2,2))
         lmax = int(frecip/recip(3,3))
         if (max(jmax,kmax,lmax) .gt. maxvec) then
            jmax = min(maxvec,jmax)
            kmax = min(maxvec,kmax)
            lmax = min(maxvec,lmax)
            write (iout,60)
   60       format (/,' KEWALD  --  Too many Reciprocal Space',
     &                 ' K-Vectors; Increase MAXVEC')
         end if
      end if
c
c     initialize the PME arrays that can be precomputed
c
      if (use_pme) then
         ntable = maxtable
         call moduli
         call fftsetup (nfft1,nfft2,nfft3,ntable,table)
      end if
c
c     print a message listing some of the Ewald parameters
c
      if (verbose) then
         if (use_pme) then
            write (iout,70)  aewald,nfft1,nfft2,nfft3,bsorder
   70       format (/,' Smooth Particle Mesh Ewald Parameters :',
     &              //,4x,'Ewald Coefficient',6x,'Charge Grid',
     &                 ' Dimensions',6x,'B-Spline Order',
     &              //,8x,f8.4,11x,3i6,12x,i6)
         else
            write (iout,80)  aewald,jmax,kmax,lmax,frecip
   80       format (/,' Regular Ewald Summation Parameters :',
     &              //,4x,'Ewald Coefficient',6x,'K-Vector',
     &                 ' Dimensions',5x,'Reciprocal Fraction',
     &              //,8x,f8.4,11x,3i5,14x,f8.4)
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ewaldcof  --  estimation of Ewald coefficient  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ewaldcof" finds a value of the Ewald coefficient such that
c     all terms beyond the specified cutoff distance will have an
c     value less than a specified tolerance
c
c
      subroutine ewaldcof (alpha,cutoff,eps)
      implicit none
      integer i,k
      real*8 alpha,cutoff,eps
      real*8 x,xlo,xhi,y
      real*8 ratio,erfc
      external erfc
c
c
c     get approximate value from cutoff and tolerance
c
      ratio = eps + 1.0d0
      x = 0.5d0
      i = 0
      dowhile (ratio .ge. eps)
         i = i + 1
         x = 2.0d0 * x
         y = x * cutoff
         ratio = erfc(y) / cutoff
      end do
c
c     use a binary search to refine the coefficient
c
      k = i + 60
      xlo = 0.0d0
      xhi = x
      do i = 1, k
         x = (xlo+xhi) / 2.0d0
         y = x * cutoff
         ratio = erfc(y) / cutoff
         if (ratio .ge. eps) then
            xlo = x
         else
            xhi = x
         end if
      end do
      alpha = x
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine moduli  --  store the inverse DFT moduli  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "moduli" sets the moduli of the inverse discrete Fourier
c     transform of the B-splines; bsmod[1-3] hold these values,
c     nfft[1-3] are the grid dimensions, bsorder is the order of
c     B-spline approximation
c
c
      subroutine moduli
      implicit none
      include 'sizes.i'
      include 'pme.i'
      integer i,big
      real*8 w,bsarray(maxfft)
      real*8 array(maxorder)
c
c
c     compute and load the moduli values
c
      w = 0.0d0
      call bspline (w,bsorder,array)
      big = max(nfft1,nfft2,nfft3)
      do i = 1, big
         bsarray(i) = 0.0d0
      end do
      do i = 2, bsorder+1
         bsarray(i) = array(i-1)
      end do
      call dftmod (bsmod1,bsarray,nfft1)
      call dftmod (bsmod2,bsarray,nfft2)
      call dftmod (bsmod3,bsarray,nfft3)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine dftmod  --  discrete Fourier transform modulus  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dftmod" computes the modulus of the discrete Fourier transform
c     of "bsarray", storing it into "bsmod"
c
c
      subroutine dftmod (bsmod,bsarray,nfft)
      implicit none
      include 'math.i'
      integer i,j,nfft
      real*8 eps,factor
      real*8 arg,sum1,sum2
      real*8 bsmod(*)
      real*8 bsarray(*)
c
c
c     get the modulus of the discrete Fourier transform
c
      eps = 1.0d-7
      factor = 2.0d0 * pi / dble(nfft)
      do i = 1, nfft
         sum1 = 0.0d0
         sum2 = 0.0d0
         do j = 1, nfft
            arg = factor * dble((i-1)*(j-1))
            sum1 = sum1 + bsarray(j)*cos(arg)
            sum2 = sum2 + bsarray(j)*sin(arg)
         end do
         bsmod(i) = sum1**2 + sum2**2
      end do
      do i = 1, nfft
         if (bsmod(i) .lt. eps)
     &      bsmod(i) = 0.5d0 * (bsmod(i-1)+bsmod(i+1))
      end do
      return
      end
