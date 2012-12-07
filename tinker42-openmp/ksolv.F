 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine ksolv  --  solvation parameter assignment  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "ksolv" assigns continuum solvation energy parameters for
c     the Eisenberg-McLachlan ASP, Ooi-Scheraga SASA or various
c     GB/SA solvation models
c
c     literature references:
c
c     L. Wesson and D. Eisenberg, "Atomic Solvation Parameters
c     Applied to Molecular Dynamics of Proteins in Solution",
c     Protein Science, 1, 227-235 (1992)  (Eisenberg-McLachlan ASP)
c
c     T. Ooi, M. Oobatake, G. Nemethy and H. A. Scheraga, "Accessible
c     Surface Areas as a Measure of the Thermodynamic Parameters of
c     Hydration of Peptides", PNAS, 84, 3086-3090 (1987)  (SASA)
c
c     M. Schaefer, C. Bartels, F. Leclerc and M. Karplus, "Effective
c     Atom Volumes for Implicit Solvent Models: Comparison between
c     Voronoi Volumes and Minimum Fluctuations Volumes", Journal of
c     Computational Chemistry, 22, 1857-1879 (2001)
c
c
      subroutine ksolv
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'couple.i'
      include 'keys.i'
      include 'math.i'
      include 'potent.i'
      include 'solute.i'
      include 'units.i'
      integer i,j,k,m,kc
      integer ia,ib,ic,id
      integer mm,nh,next
      integer atmnum,atmmas
      real*8 ri,ri2,rk,rk2
      real*8 r,r2,r4,rab,rbc
      real*8 cosine,factor
      real*8 h,ratio,term
      real*8 c1,c2,c3,pi2
      real*8 width,qterm,temp
      real*8 alpha,alpha2,alpha4
      real*8 vk,prod2,prod4
      real*8 fik,tik2,qik,uik
      real*8 s2ik,s3ik,omgik
      logical amide
      character*20 keyword
      character*20 value
      character*120 record
c
c
c     default is to not use the continuum solvation term
c
      use_solv = .false.
      use_gbsa = .false.
      solvtyp = '      '
c
c     search keywords for continuum solvation commands
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'SOLVATE ') then
            use_solv = .true.
            use_gbsa = .false.
            solvtyp = 'ASP'
            call getword (record,value,next)
            call upcase (value)
            if (value(1:3) .eq. 'ASP') then
               solvtyp = 'ASP'
            else if (value(1:4) .eq. 'SASA') then
               solvtyp = 'SASA'
            else if (value(1:5) .eq. 'ONION') then
               use_gbsa = .true.
               solvtyp = 'ONION'
            else if (value(1:5) .eq. 'STILL') then
               use_gbsa = .true.
               solvtyp = 'STILL'
            else if (value(1:3) .eq. 'HCT') then
               use_gbsa = .true.
               solvtyp = 'HCT'
            else if (value(1:3) .eq. 'ACE') then
               use_gbsa = .true.
               solvtyp = 'ACE'
            else if (value(1:4) .eq. 'GBSA') then
               use_gbsa = .true.
               solvtyp = 'STILL'
            end if
         end if
      end do
c
c     set the dielectric offset to the radii for GB/SA methods
c
      if (use_gbsa)  doffset = -0.09d0
c
c     set offset and scaling values for analytical Still method
c
      if (solvtyp .eq. 'STILL') then
         p1 = 0.073d0
         p2 = 0.921d0
         p3 = 6.211d0
         p4 = 15.236d0
         p5 = 1.254d0
         if (.not. use_bond)  call kbond
         if (.not. use_angle)  call kangle
      end if
c
c     set overlap scale factors for Hawkins-Cramer-Truhlar method
c
      if (solvtyp .eq. 'HCT') then
         do i = 1, n
            shct(i) = 0.80d0
            atmnum = atomic(i)
            if (atmnum .eq. 1)  shct(i) = 0.85d0
            if (atmnum .eq. 6)  shct(i) = 0.72d0
            if (atmnum .eq. 7)  shct(i) = 0.79d0
            if (atmnum .eq. 8)  shct(i) = 0.85d0
            if (atmnum .eq. 9)  shct(i) = 0.88d0
            if (atmnum .eq. 15)  shct(i) = 0.86d0
            if (atmnum .eq. 16)  shct(i) = 0.96d0
            if (atmnum .eq. 26)  shct(i) = 0.88d0
         end do
      end if
c
c     set the Gaussian width factor for the ACE method
c
      if (solvtyp .eq. 'ACE') then
         width = 1.2d0
      end if
c
c     assign surface area factors for various GB/SA methods
c
      if (solvtyp .eq. 'ONION') then
         do i = 1, n
            asolv(i) = 0.0072d0
         end do
      else if (solvtyp .eq. 'STILL') then
         do i = 1, n
            asolv(i) = 0.0049d0
         end do
      else if (solvtyp .eq. 'HCT') then
         do i = 1, n
            asolv(i) = 0.0054d0
         end do
      else if (solvtyp .eq. 'ACE') then
         do i = 1, n
            asolv(i) = 0.0030d0
         end do
      end if
c
c     assign the Eisenberg-McLachlan ASP solvation parameters;
c     parameters only available for protein-peptide groups
c
      if (solvtyp .eq. 'ASP') then
         do i = 1, n
            atmnum = atomic(i)
            if (atmnum .eq. 6) then
               rsolv(i) = 1.9d0
               asolv(i) = 0.004d0
            else if (atmnum .eq. 7) then
               rsolv(i) = 1.7d0
               asolv(i) = -0.113d0
               if (n12(i) .eq. 4) then
                  asolv(i) = -0.169d0
               end if
            else if (atmnum .eq. 8) then
               rsolv(i) = 1.4d0
               asolv(i) = -0.113d0
               if (n12(i).eq.1 .and. atomic(i12(1,i)).eq.6) then
                  do j = 1, n13(i)
                     k = i13(j,i)
                     if (n12(k).eq.1 .and. atomic(k).eq.8) then
                        asolv(i) = -0.166d0
                     end if
                  end do
               end if
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 15)  asolv(i) = -0.140d0
               end do
            else if (atmnum .eq. 15) then
               rsolv(i) = 1.9d0
               asolv(i) = -0.140d0
            else if (atmnum .eq. 16) then
               rsolv(i) = 1.8d0
               asolv(i) = -0.017d0
            else
               rsolv(i) = 0.0d0
               asolv(i) = 0.0d0
            end if
         end do
      end if
c
c     assign the Ooi-Scheraga SASA solvation parameters;
c     parameters only available for protein-peptide groups
c
      if (solvtyp .eq. 'SASA') then
         do i = 1, n
            atmnum = atomic(i)
            if (atmnum .eq. 6) then
               rsolv(i) = 2.0d0
               asolv(i) = 0.008d0
               if (n12(i) .eq. 3) then
                  rsolv(i) = 1.75d0
                  asolv(i) = -0.008d0
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k) .eq. 8) then
                        rsolv(i) = 1.55d0
                        asolv(i) = 0.427d0
                     end if
                  end do
               end if
            else if (atmnum .eq. 7) then
               rsolv(i) = 1.55d0
               asolv(i) = -0.132d0
               if (n12(i) .eq. 4)  asolv(i) = -1.212d0
            else if (atmnum .eq. 8) then
               rsolv(i) = 1.4d0
               if (n12(i) .eq. 1) then
                  asolv(i) = -0.038d0
                  if (atomic(i12(1,i)) .eq. 6) then
                     do j = 1, n13(i)
                        k = i13(j,i)
                        if (n12(k).eq.1 .and. atomic(k).eq.8) then
                           asolv(i) = -0.770d0
                        end if
                     end do
                  end if
               else if (n12(i) .eq. 2) then
                  asolv(i) = -0.172d0
               end if
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 15)  asolv(i) = -0.717d0
               end do
            else if (atmnum .eq. 15) then
               rsolv(i) = 2.1d0
               asolv(i) = 0.0d0
            else if (atmnum .eq. 16) then
               rsolv(i) = 2.0d0
               asolv(i) = -0.021d0
            else if (atmnum .eq. 17) then
               rsolv(i) = 2.0d0
               asolv(i) = 0.012d0
            else
               rsolv(i) = 0.0d0
               asolv(i) = 0.0d0
            end if
         end do
      end if
c
c     assign standard radii for GB/SA methods other than ACE;
c     taken from Macromodel and OPLS-AA, except for hydrogens
c
      if (use_gbsa .and. solvtyp.ne.'ACE') then
         do i = 1, n
            atmnum = atomic(i)
            if (atmnum .eq. 1) then
               rsolv(i) = 1.25d0
               k = i12(1,i)
               if (atomic(k) .eq. 7)  rsolv(i) = 1.15d0
               if (atomic(k) .eq. 8)  rsolv(i) = 1.05d0
            else if (atmnum .eq. 3) then
               rsolv(i) = 1.432d0
            else if (atmnum .eq. 6) then
               rsolv(i) = 1.90d0
               if (n12(i) .eq. 3)  rsolv(i) = 1.875d0
               if (n12(i) .eq. 2)  rsolv(i) = 1.825d0
            else if (atmnum .eq. 7) then
               rsolv(i) = 1.7063d0
               if (n12(i) .eq. 4)  rsolv(i) = 1.625d0
               if (n12(i) .eq. 1)  rsolv(i) = 1.60d0
            else if (atmnum .eq. 8) then
               rsolv(i) = 1.535d0
               if (n12(i) .eq. 1)  rsolv(i) = 1.48d0
            else if (atmnum .eq. 9) then
               rsolv(i) = 1.47d0
            else if (atmnum .eq. 10) then
               rsolv(i) = 1.39d0
            else if (atmnum .eq. 11) then
               rsolv(i) = 1.992d0
            else if (atmnum .eq. 12) then
               rsolv(i) = 1.70d0
            else if (atmnum .eq. 14) then
               rsolv(i) = 1.80d0
            else if (atmnum .eq. 15) then
               rsolv(i) = 1.87d0
            else if (atmnum .eq. 16) then
               rsolv(i) = 1.775d0
            else if (atmnum .eq. 17) then
               rsolv(i) = 1.735d0
            else if (atmnum .eq. 18) then
               rsolv(i) = 1.70d0
            else if (atmnum .eq. 19) then
               rsolv(i) = 2.123d0
            else if (atmnum .eq. 20) then
               rsolv(i) = 1.817d0
            else if (atmnum .eq. 35) then
               rsolv(i) = 1.90d0
            else if (atmnum .eq. 36) then
               rsolv(i) = 1.812d0
            else if (atmnum .eq. 37) then
               rsolv(i) = 2.26d0
            else if (atmnum .eq. 53) then
               rsolv(i) = 2.10d0
            else if (atmnum .eq. 54) then
               rsolv(i) = 1.967d0
            else if (atmnum .eq. 55) then
               rsolv(i) = 2.507d0
            else if (atmnum .eq. 56) then
               rsolv(i) = 2.188d0
            else
               rsolv(i) = 2.0d0
            end if
         end do
      end if
c
c     compute the atomic volumes for the analytical Still method
c
      if (solvtyp .eq. 'STILL') then
         do i = 1, n
            vsolv(i) = (4.0d0*pi/3.0d0) * rsolv(i)**3
            ri = rsolv(i)
            ri2 = ri * ri
            do j = 1, n12(i)
               k = i12(j,i)
               rk = rsolv(k)
               r = 1.01d0 * bl(bndlist(j,i))
               ratio = (rk*rk-ri2-r*r) / (2.0d0*ri*r)
               h = ri * (1.0d0+ratio)
               term = (pi/3.0d0) * h * h * (3.0d0*ri-h)
               vsolv(i) = vsolv(i) - term
            end do
         end do
c
c     get self-, 1-2 and 1-3 polarization for analytical Still method
c
         do i = 1, n
            gpol(i) = -0.5d0 * electric / (rsolv(i)+doffset+p1)
         end do
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            r = bl(i)
            r4 = r**4
            gpol(ia) = gpol(ia) + p2*vsolv(ib)/r4
            gpol(ib) = gpol(ib) + p2*vsolv(ia)/r4
         end do
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            factor = 1.0d0
            do j = 1, n12(ia)
               id = i12(j,ia)
               if (id .eq. ic) then
                  factor = 0.0d0
               else if (id .ne. ib) then
                  do k = 1, n12(ic)
                     if (i12(k,ic) .eq. id) then
                        factor = 0.5d0
                     end if
                  end do
               end if
            end do
            do j = 1, n12(ib)
               if (i12(j,ib) .eq. ia) then
                  rab = bl(bndlist(j,ib))
               else if (i12(j,ib) .eq. ic) then
                  rbc = bl(bndlist(j,ib))
               end if
            end do
            cosine = cos(anat(i)/radian)
            r2 = rab**2 + rbc**2 - 2.0d0*rab*rbc*cosine
            r4 = r2 * r2
            gpol(ia) = gpol(ia) + factor*p3*vsolv(ic)/r4
            gpol(ic) = gpol(ic) + factor*p3*vsolv(ia)/r4
         end do
      end if
c
c     assign the atomic radii and volumes for the ACE method;
c     volumes taken from average Voronoi values with hydrogens
c
      if (solvtyp .eq. 'ACE') then
         do i = 1, n
            atmnum = atomic(i)
            atmmas = nint(mass(i))
            if (atmnum .eq. 1) then
               rsolv(i) = 1.468d0
               vsolv(i) = 11.0d0
               k = i12(1,i)
               if (atomic(k).eq.6 .and. n12(k).eq.4) then
                  vsolv(i) = 11.895d0
               else if (atomic(k).eq.6 .and. n12(k).eq.3) then
                  vsolv(i) = 13.242d0
               else if (atomic(k).eq.7 .and. n12(k).eq.4) then
                  rsolv(i) = 0.60d0
                  vsolv(i) = 9.138d0
               else if (atomic(k).eq.7 .or. atomic(k).eq.8) then
                  rsolv(i) = 0.60d0
                  vsolv(i) = 9.901d0
               else if (atomic(k).ne.16) then
                  rsolv(i) = 1.468d0
                  vsolv(i) = 13.071d0
               end if
            else if (atmnum .eq. 6) then
               rsolv(i) = 2.49d0
               vsolv(i) = 7.0d0
               nh = 0
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 1)  nh = nh + 1
               end do
               if (n12(i) .eq. 4) then
                  if (nh .eq. 3) then
                     vsolv(i) = 3.042d0
                  else if (nh .eq. 2) then
                     vsolv(i) = 3.743d0
                  else if (nh .eq. 1) then
                     vsolv(i) = 4.380d0
                  end if
               else if (n12(i) .eq. 3) then
                  if (nh .eq. 1) then
                     rsolv(i) = 2.10d0
                     vsolv(i) = 7.482d0
                  else if (nh .eq. 0) then
                     rsolv(i) = 2.10d0
                     vsolv(i) = 8.288d0
                  end if
                  do j = 1, n12(i)
                     k = i12(1,j)
                     if (atomic(k).eq.8 .and. n12(k).eq.1) then
                        rsolv(i) = 2.10d0
                        vsolv(i) = 7.139d0
                     end if
                  end do
               end if
               if (atmmas .eq. 15) then
                  rsolv(i) = 2.165d0
                  vsolv(i) = 33.175d0
               else if (atmmas .eq. 14) then
                  rsolv(i) = 2.235d0
                  vsolv(i) = 20.862d0
               else if (atmmas.eq.13 .and. n12(i).eq.2) then
                  rsolv(i) = 2.10d0
                  vsolv(i) = 20.329d0
               else if (atmmas .eq. 13) then
                  rsolv(i) = 2.365d0
                  vsolv(i) = 11.784d0
               end if
            else if (atmnum .eq. 7) then
               rsolv(i) = 1.60d0
               vsolv(i) = 6.0d0
               nh = 0
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 1)  nh = nh + 1
               end do
               if (n12(i) .eq. 4) then
                  if (nh .eq. 3) then
                     vsolv(i) = 2.549d0
                  else if (nh .eq. 2) then
                     vsolv(i) = 3.304d0
                  end if
               else if (n12(i) .eq. 3) then
                  amide = .false.
                  do j = 1, n12(i)
                     m = i12(j,i)
                     if (atomic(m) .eq. 6) then
                        do k = 1, n12(m)
                           mm = i12(k,m)
                           if (atomic(mm).eq.8 .and. n12(mm).eq.1) then
                              amide = .true.
                           end if
                        end do
                     end if
                  end do
                  if (amide) then
                     if (nh .eq. 0) then
                        vsolv(i) = 7.189d0
                     else if (nh .eq. 1) then
                        vsolv(i) = 6.030d0
                     else if (nh .eq. 2) then
                        vsolv(i) = 5.693d0
                     end if
                  else
                     if (nh .eq. 2) then
                        vsolv(i) = 5.677d0
                     else if (nh .eq. 2) then
                        vsolv(i) = 6.498d0
                     end if
                  end if
               end if
            else if (atmnum .eq. 8) then
               rsolv(i) = 1.60d0
               vsolv(i) = 12.0d0
               if (n12(i) .eq. 1) then
                  vsolv(i) = 13.532d0
                  k = i12(1,i)
                  if (atomic(k) .eq. 15) then
                     vsolv(i) = 17.202d0
                  else
                     do j = 1, n13(i)
                        k = i13(j,i)
                        if (atomic(j).eq.8 .and. n12(j).eq.1) then
                           vsolv(i) = 15.400d0
                        end if
                     end do
                  end if
               else if (n12(i) .eq. 2) then
                  vsolv(i) = 10.642d0
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k) .eq. 15)  vsolv(i) = 11.416d0
                  end do
               end if
            else if (atmnum .eq. 12) then
               rsolv(i) = 1.0d0
               vsolv(i) = 15.235d0
           else if (atmnum .eq. 15) then
               rsolv(i) = 1.89d0
               vsolv(i) = 6.131d0
            else if (atmnum .eq. 16) then
               rsolv(i) = 1.89d0
               vsolv(i) = 17.232d0
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 16)  vsolv(i) = 18.465d0
               end do
            else if (atmnum .eq. 26) then
               rsolv(i) = 0.65d0
               vsolv(i) = 9.951d0
            else
               rsolv(i) = 0.0d0
               vsolv(i) = 0.0d0
            end if
         end do
c
c     calculate the pairwise parameters for the ACE method
c
         c1 = 4.0d0 / (3.0d0*pi)
         c2 = 77.0d0 * pi * sqrttwo / 512.0d0
         c3 = 2.0d0 * pi * sqrtpi
         pi2 = 1.0d0 / (pi*pi)
         do i = 1, n
            ic = class(i)
            ri = rsolv(i)
            ri2 = ri * ri
            do k = 1, n
               kc = class(k)
               rk = rsolv(kc)
               vk = vsolv(kc)
               rk2 = rk * rk
               alpha = max(width,ri/rk)
               alpha2 = alpha * alpha
               alpha4 = alpha2 * alpha2
               prod2 = alpha2 * rk2
               prod4 = prod2 * prod2
               ratio = alpha2 * rk2 / ri2
               tik2 = 0.5d0 * pi * ratio
               temp = 1.0d0 / (1.0d0+2.0d0*tik2)
               fik = 2.0d0/(1.0d0+tik2) - temp
               qik = tik2 * sqrt(temp)
               qterm = qik - atan(qik)
               if (k .ne. i) then
                  omgik = vk * qterm * pi2 / prod4
               else
                  omgik = c1 * qterm / (alpha4 * ri)
               end if
               s2ik = 3.0d0 * qterm * prod2
     &                   / ((3.0d0+fik)*qik-4.0d0*atan(qik))
               s3ik = s2ik * sqrt(s2ik)
               uik = c2 * ri / (1.0d0-(c3*s3ik*ri*omgik/vk))
               wace(ic,kc) = omgik
               s2ace(ic,kc) = s2ik
               uace(ic,kc) = uik
            end do
         end do
      end if
      return
      end
