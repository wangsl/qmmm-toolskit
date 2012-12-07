 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine esolv  --  continuum solvation energy  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "esolv" calculates the continuum solvation energy via
c     either the Eisenberg-McLachlan ASP model, Ooi-Scheraga
c     SASA model, various GB/SA methods or the ACE model
c
c
      subroutine esolv
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'energi.i'
      include 'math.i'
      include 'potent.i'
      include 'solute.i'
      include 'warp.i'
      integer i
      real*8 e,ai,ri,rb
      real*8 term,probe
      real*8 aes(maxatm)
      real*8 des(3,maxatm)
c
c
c     zero out the continuum solvation energy
c
      es = 0.0d0
c
c     set a value for the solvent molecule probe radius
c
      probe = 1.4d0
c
c     get the Born radius values for the GB/SA models
c
      if (use_gbsa)  call born
c
c     compute the nonpolar solvation via the ACE approximation
c
      if (use_gbsa .and. solvtyp.ne.'ONION') then
         term = 4.0d0 * pi
         do i = 1, n
            ai = asolv(i)
            ri = rsolv(i)
            rb = rborn(i)
            if (rb .ne. 0.0d0) then
               e = ai * term * (ri+probe)**2 * (ri/rb)**6
               es = es + e
            end if
         end do
c
c     compute the nonpolar solvation via exact surface area
c
      else
         call surface (es,aes,des,rsolv,asolv,probe)
      end if
c
c     get the polarization energy term for GB/SA solvation
c
      if (use_gbsa) then
         if (use_smooth) then
            call egbsa0b
         else
            call egbsa0a
         end if
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine egbsa0a  --  GB/SA polarization energy term  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "egbsa0a" calculates the generalized Born polarization energy
c     for the GB/SA solvation models
c
c
      subroutine egbsa0a
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'charge.i'
      include 'energi.i'
      include 'group.i'
      include 'shunt.i'
      include 'solute.i'
      include 'units.i'
      include 'usage.i'
      integer i,k,ii,kk
      real*8 e,f,fi,fik
      real*8 dwater,fgrp
      real*8 rb2,rm2,fgb,fgm
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,r3,r4
      real*8 r5,r6,r7
      real*8 shift,taper,trans
      logical proceed,iuse
c
c
c     set the solvent dielectric and energy conversion factor
c
      if (nion .eq. 0)  return
      dwater = 78.3d0
      f = -electric * (1.0d0 - 1.0d0/dwater)
c
c     set cutoff distances and switching function coefficients
c
      call switch ('CHARGE')
c
c     calculate GB/SA electrostatic polarization energy term
c
      do ii = 1, nion
         i = iion(ii)
         iuse = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         do kk = ii, nion
            k = iion(kk)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  fik = fi * pchg(kk)
                  rb2 = rborn(i) * rborn(k)
                  fgb = sqrt(r2 + rb2*exp(-0.25d0*r2/rb2))
                  e = fik / fgb
c
c     use shifted energy switching if near the cutoff distance
c
                  rm2 = (0.5d0 * (off+cut))**2
                  fgm = sqrt(rm2 + rb2*exp(-0.25d0*rm2/rb2))
                  shift = fik / fgm
                  e = e - shift
                  if (r2 .gt. cut2) then
                     r = sqrt(r2)
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     r6 = r3 * r3
                     r7 = r3 * r4
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                               + f3*r3 + f2*r2 + f1*r + f0)
                     e = e * taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall GB/SA solvation energy component
c
                  if (i .eq. k)  e = 0.5d0 * e
                  es = es + e
               end if
            end if
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine egbsa0b  --  GB/SA polarization for smoothing  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "egbsa0b" calculates the generalized Born polarization energy
c     for the GB/SA solvation models for use with potential smoothing
c     methods via analogy to the smoothing of Coulomb's law
c
c
      subroutine egbsa0b
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'charge.i'
      include 'energi.i'
      include 'group.i'
      include 'solute.i'
      include 'units.i'
      include 'usage.i'
      include 'warp.i'
      integer i,k,ii,kk
      real*8 e,fgrp
      real*8 f,fi,fik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 dwater,width
      real*8 erf,sterm
      real*8 r2,fgb,rb2
      logical proceed,iuse
      external erf
c
c
c     set the solvent dielectric and energy conversion factor
c
      if (nion .eq. 0)  return
      dwater = 78.3d0
      f = -electric * (1.0d0 - 1.0d0/dwater)
c
c     set the extent of smoothing to be performed
c
      sterm = 0.5d0 / sqrt(diffc)
c
c     calculate GB/SA electrostatic polarization energy term
c
      do ii = 1, nion
         i = iion(ii)
         iuse = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         do kk = ii, nion
            k = iion(kk)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               fik = fi * pchg(kk)
               rb2 = rborn(i) * rborn(k)
               fgb = sqrt(r2 + rb2*exp(-0.25d0*r2/rb2))
               e = fik / fgb
c
c     use a smoothable GB/SA analogous to Coulomb's law solution
c
               if (deform .gt. 0.0d0) then
                  width = deform + 0.15d0*rb2*exp(-0.006d0*rb2/deform)
                  width = sterm / sqrt(width)
                  e = e * erf(width*fgb)
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the overall GB/SA solvation energy component
c
               if (i .eq. k)  e = 0.5d0 * e
               es = es + e
            end if
         end do
      end do
      return
      end
