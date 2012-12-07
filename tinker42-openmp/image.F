 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine image  --  pairwise distance of minimum image  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "image" takes the components of pairwise distance between
c     two points in the same or neighboring periodic boxes and
c     converts to the components of the minimum image distance
c
c
      subroutine image (xr,yr,zr,i)
      implicit none
      include 'sizes.i'
      include 'boxes.i'
      include 'cell.i'
      integer i
      real*8 xr,yr,zr
      real*8 xfrac,yfrac,zfrac
      real*8 xshift,yshift,zshift
      real*8 xcycle,ycycle,zcycle
      real*8 xcycle2,ycycle2,zcycle2
c
c
c     set dimensions for either single box or replicated cell
c
      
      if (i .ge. 0) then
         xcycle = xcell
         ycycle = ycell
         zcycle = zcell
         xcycle2 = xcell2
         ycycle2 = ycell2
         zcycle2 = zcell2
      else
         xcycle = xbox
         ycycle = ybox
         zcycle = zbox
         xcycle2 = xbox2
         ycycle2 = ybox2
         zcycle2 = zbox2
      end if
c
c     compute the distance to translate along each cell axis
c
      if (i .le. 0) then
         xshift = 0.0d0
         yshift = 0.0d0
         zshift = 0.0d0
      else
         xshift = icell(1,i) * xbox
         yshift = icell(2,i) * ybox
         zshift = icell(3,i) * zbox
      end if
c
c     for orthogonal lattice, find the desired image directly
c
      if (orthogonal) then
         xr = xr + xshift
         dowhile (xr .gt. xcycle2)
            xr = xr - xcycle
         end do
         dowhile (xr .lt. -xcycle2)
            xr = xr + xcycle
         end do
         yr = yr + yshift
         dowhile (yr .gt. ycycle2)
            yr = yr - ycycle
         end do
         dowhile (yr .lt. -ycycle2)
            yr = yr + ycycle
         end do
         zr = zr + zshift
         dowhile (zr .gt. zcycle2)
            zr = zr - zcycle
         end do
         dowhile (zr .lt. -zcycle2)
            zr = zr + zcycle
         end do
c
c     for monoclinic lattice, convert "xr" and "zr" to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
      else if (monoclinic) then
         zfrac = zr / beta_sin
         xfrac = xr - zfrac*beta_cos
         xfrac = xfrac + xshift
         dowhile (xfrac .gt. xcycle2)
            xfrac = xfrac - xcycle
         end do
         dowhile (xfrac .lt. -xcycle2)
            xfrac = xfrac + xcycle
         end do
         yr = yr + yshift
         dowhile (yr .gt. ycycle2)
            yr = yr - ycycle
         end do
         dowhile (yr .lt. -ycycle2)
            yr = yr + ycycle
         end do
         zfrac = zfrac + zshift
         dowhile (zfrac .gt. zcycle2)
            zfrac = zfrac - zcycle
         end do
         dowhile (zfrac .lt. -zcycle2)
            zfrac = zfrac + zcycle
         end do
         xr = xfrac + zfrac*beta_cos
         zr = zfrac * beta_sin
c
c     for triclinic lattice, convert pairwise components to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
      else if (triclinic) then
         zfrac = zr / gamma_term
         yfrac = (yr - zfrac*beta_term) / gamma_sin
         xfrac = xr - yfrac*gamma_cos - zfrac*beta_cos
         xfrac = xfrac + xshift
         dowhile (xfrac .gt. xcycle2)
            xfrac = xfrac - xcycle
         end do
         dowhile (xfrac .lt. -xcycle2)
            xfrac = xfrac + xcycle
         end do
         yfrac = yfrac + yshift
         dowhile (yfrac .gt. ycycle2)
            yfrac = yfrac - ycycle
         end do
         dowhile (yfrac .lt. -ycycle2)
            yfrac = yfrac + ycycle
         end do
         zfrac = zfrac + zshift
         dowhile (zfrac .gt. zcycle2)
            zfrac = zfrac - zcycle
         end do
         dowhile (zfrac .lt. -zcycle2)
            zfrac = zfrac + zcycle
         end do
         xr = xfrac + yfrac*gamma_cos + zfrac*beta_cos
         yr = yfrac*gamma_sin + zfrac*beta_term
         zr = zfrac * gamma_term
c
c     for truncated octahedron, use orthogonal box equations,
c     then perform extra tests to remove corner pieces
c
      else if (octahedron) then
         dowhile (xr .gt. xbox2)
            xr = xr - xbox
         end do
         dowhile (xr .lt. -xbox2)
            xr = xr + xbox
         end do
         dowhile (yr .gt. ybox2)
            yr = yr - ybox
         end do
         dowhile (yr .lt. -ybox2)
            yr = yr + ybox
         end do
         dowhile (zr .gt. zbox2)
            zr = zr - zbox
         end do
         dowhile (zr .lt. -zbox2)
            zr = zr + zbox
         end do
         if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
            xr = xr - sign(xbox2,xr)
            yr = yr - sign(ybox2,yr)
            zr = zr - sign(zbox2,zr)
         end if
      end if
      return
      end
