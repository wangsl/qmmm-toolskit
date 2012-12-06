 
C     $Id: born.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine born  --  find Born radii for GB/SA solvation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "born" computes the Born radius of each atom for use with
c     the various GB/SA solvation models
c
c     literature references:
c
c     W. C. Still, A. Tempczyk, R. C. Hawley and T. Hendrickson,
c     "A Semianalytical Treatment of Solvation for Molecular
c     Mechanics and Dynamics", J. Amer. Chem. Soc., 112, 6127-6129
c     (1990)  (Original Still; see supplimentary material)
c
c     D. Qiu, P.S. Shenkin, F.P. Hollinger and W.C. Still, "The
c     GB/SA Continuum Model for Solvation. A Fast Analytical Method
c     for the Calculation of Approximate Radii", J. Phys. Chem. A,
c     101, 3005-3014 (1997)  (Analytical Still Method)
c
c     G. D. Hawkins, C. J. Cramer and D. G. Truhlar, "Parametrized
c     Models of Aqueous Free Energies of Solvation Based on Pairwise
c     Descreening of Solute Atomic Charges from a Dielectric Medium",
c     J. Phys. Chem., 100, 19824-19839 (1996)  (HCT Method)
c
c     M. Schaefer C. Bartels and M. Karplus, "Solution Conformations
c     and Thermodynamics of Structured Peptides: Molecular Dynamics
c     Simulation with an Implicit Solvation Model", J. Mol. Biol.,
c     284, 835-848 (1998)  (Analytical Continuum Electrostatics, ACE)
c
c
      subroutine born
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'solute.i'
      include 'units.i'
      integer i,j,k,it,kt
      integer skip(maxatm)
      real*8 area,rold,t
      real*8 shell,fraction
      real*8 inner,outer,tinit
      real*8 ratio,total
      real*8 xi,yi,zi,ri
      real*8 rk,sk,sk2,sum
      real*8 lik,lik2,uik,uik2
      real*8 xr,yr,zr,rvdw
      real*8 r,r2,r3,r4
      real*8 gpi,pip5,p5inv
      real*8 theta,term,ccf
      real*8 expterm,rmu
      real*8 b0,gself,third
      real*8 roff(maxatm)
      logical done
c
c
c     compute atomic radii modified by the dielectric offset
c
      do i = 1, n
         roff(i) = rsolv(i) + doffset
      end do
c
c     get the Born radii via the numerical Still method
c
      if (solvtyp .eq. 'ONION') then
         tinit = 0.1d0
         ratio = 1.5d0
         do i = 1, n
            t = tinit
            rold = roff(i)
            total = 0.0d0
            done = .false.
            dowhile (.not. done)
               roff(i) = roff(i) + 0.5d0*t
               call surfatom (i,area,roff)
               fraction = area / (4.0d0*pi*roff(i)**2)
               if (fraction .lt. 0.99d0) then
                  inner = roff(i) - 0.5d0*t
                  outer = inner + t
                  shell = 1.0d0/inner - 1.0d0/outer
                  total = total + fraction*shell
                  roff(i) = roff(i) + 0.5d0*t
                  t = ratio * t
               else
                  inner = roff(i) - 0.5d0*t
                  total = total + 1.0d0/inner
                  done = .true.
               end if
            end do
            rborn(i) = 1.0d0 / total
            roff(i) = rold
         end do
c
c     get the Born radii via the analytical Still method;
c     note this code only loops over the variable parts
c
      else if (solvtyp .eq. 'STILL') then
         do i = 1, n
            skip(i) = 0
         end do
         p5inv = 1.0d0 / p5
         pip5 = pi * p5
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            gpi = gpol(i)
            skip(i) = i
            do j = 1, n12(i)
               skip(i12(j,i)) = i
            end do
            do j = 1, n13(i)
               skip(i13(j,i)) = i
            end do
            do k = 1, n
               if (skip(k) .ne. i) then
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  r2 = xr**2 + yr**2 + zr**2
                  r4 = r2 * r2
                  rvdw = rsolv(i) + rsolv(k)
                  ratio = r2 / (rvdw*rvdw)
                  if (ratio .gt. p5inv) then
                     ccf = 1.0d0
                  else
                     theta = ratio * pip5
                     term = 0.5d0 * (1.0d0-cos(theta))
                     ccf = term * term
                  end if
                  gpi = gpi + p4*ccf*vsolv(k)/r4
               end if
            end do
            rborn(i) = -0.5d0 * electric / gpi
         end do
c
c     get the Born radii via the Hawkins-Cramer-Truhlar method
c
      else if (solvtyp .eq. 'HCT') then
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ri = roff(i)
            sum = 1.0d0 / ri
            do k = 1, n
               if (i .ne. k) then
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  r2 = xr**2 + yr**2 + zr**2
                  r = sqrt(r2)
                  rk = roff(k)
                  sk = rk * shct(k)
                  sk2 = sk * sk
                  if (ri .lt. r+sk) then
                     lik = 1.0d0 / max(ri,r-sk)
                     uik = 1.0d0 / (r+sk)
                     lik2 = lik * lik
                     uik2 = uik * uik
                     term = lik - uik + 0.25d0*r*(uik2-lik2)
     &                         + (0.5d0/r)*log(uik/lik)
     &                         + (0.25d0*sk2/r)*(lik2-uik2)
                     sum = sum - 0.5d0*term
                  end if
               end if
            end do
            rborn(i) = 1.0d0 / sum
            rborn(i) = max(ri,rborn(i))
         end do
c
c     get the Born radii via analytical continuum electrostatics
c
      else if (solvtyp .eq. 'ACE') then
         third = 1.0d0 / 3.0d0
         b0 = 0.0d0
         do i = 1, n
            b0 = b0 + vsolv(i)
         end do
         b0 = (0.75d0*b0/pi)**third
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ri = rsolv(i)
            it = class(i)
            gself = 1.0d0/ri + 2.0d0*wace(it,it)
            do k = 1, n
               if (k .ne. i) then
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  kt = class(k)
                  r2 = xr**2 + yr**2 + zr**2
                  r3 = r2 * sqrt(r2)
                  r4 = r2 * r2
                  expterm = wace(it,kt) * exp(-r2/s2ace(it,kt))
                  rmu = r4 + uace(it,kt)**4
                  term = (vsolv(k)/(8.0d0*pi)) * (r3/rmu)**4
                  gself = gself - 2.0d0*(expterm+term)
               end if
            end do
            if (gself .ge. 0.5d0/b0) then
               rborn(i) = 1.0d0 / gself
            else
               rborn(i) = 2.0d0 * b0 * (1.0d0+b0*gself)
            end if
         end do
      end if
c
c     write out the final Born radius value for each atom
c
      if (debug) then
         write (iout,10)
   10    format (/,' Born Radii for Individual Atoms :',/)
         k = 1
         dowhile (k .le. n)
            write (iout,20)  (i,rborn(i),i=k,min(k+4,n))
   20       format (1x,5(i7,f8.3))
            k = k + 5
         end do
      end if
      return
      end
