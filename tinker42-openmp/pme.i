 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  pme.i  --  parameters for particle mesh Ewald summation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     maxfft     maximum number of points along each FFT direction
c     maxorder   maximum order of the B-spline approximation
c     maxtable   maximum size of the FFT table array
c     maxgrid    maximum dimension of the PME charge grid array
c
c     bsmod1     B-spline moduli along the a-axis direction
c     bsmod2     B-spline moduli along the b-axis direction
c     bsmod3     B-spline moduli along the c-axis direction
c     table      intermediate array used by the FFT calculation
c     nfft1      number of grid points along the a-axis direction
c     nfft2      number of grid points along the b-axis direction
c     nfft3      number of grid points along the c-axis direction
c     bsorder    order of the PME B-spline approximation
c
c
      integer maxfft,maxorder,maxtable,maxgrid
      parameter (maxfft=100)
      parameter (maxorder=12)
      parameter (maxtable=4*maxfft+15)
      parameter (maxgrid=2*(maxfft/2)+1)
      integer nfft1,nfft2,nfft3,bsorder
      real*8 bsmod1,bsmod2,bsmod3
      real*8 table,qgrid
      common /pme/ bsmod1(maxfft),bsmod2(maxfft),bsmod3(maxfft),
     &             table(maxtable,3),qgrid(2,maxgrid,maxgrid,maxgrid),
     &             nfft1,nfft2,nfft3,bsorder
