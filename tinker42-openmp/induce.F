 
C     $Id$
c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine induce  --  evaluate induced dipole moments  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "induce" computes the induced dipole moment at each
c     polarizable site due to direct or mutual polarization;
c     assumes that multipole components have already been
c     rotated into the global coordinate frame
c
c
      subroutine induce
      implicit none
      include 'sizes.i'
      include 'cutoff.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'polar.i'
      include 'units.i'
      integer i,j,k
      real*8 norm
      logical header
c
c
c     choose pairwise double loop or Ewald summation version
c
      if (use_ewald) then
         call induce0b
      else
         call induce0a
      end if
c
c     print out a list of the final induced dipole moments
c
      header = .true.
      if (debug) then
         do i = 1, npole
            if (polarity(i) .ne. 0.0d0) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Induced Dipole Moments (Debyes) :')
                  if (digits .ge. 8) then
                     write (iout,20)
   20                format (/,4x,'Atom',14x,'X',15x,'Y',15x,'Z',
     &                          15x,'Total'/)
                  else if (digits .ge. 6) then
                     write (iout,30)
   30                format (/,4x,'Atom',14x,'X',13x,'Y',13x,'Z',
     &                          12x,'Total'/)
                  else
                     write (iout,40)
   40                format (/,4x,'Atom',14x,'X',11x,'Y',11x,'Z',
     &                          9x,'Total'/)
                  end if
               end if
               k = ipole(i)
               norm = sqrt(uind(1,i)**2+uind(2,i)**2+uind(3,i)**2)
               if (digits .ge. 8) then
                  write (iout,50)  k,(debye*uind(j,i),j=1,3),debye*norm
   50             format (i8,3x,4f16.8)
               else if (digits .ge. 6) then
                  write (iout,60)  k,(debye*uind(j,i),j=1,3),debye*norm
   60             format (i8,4x,4f14.6)
               else
                  write (iout,70)  k,(debye*uind(j,i),j=1,3),debye*norm
   70             format (i8,5x,4f12.4)
               end if
            end if
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine induce0a  --  induced dipoles via double loop  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "induce0a" computes the induced dipole moment at each
c     polarizable site using a pairwise double loop
c
c
      subroutine induce0a
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      include 'units.i'
      integer i,j,k,m
      integer iter,maxiter
      integer ii,ix,iz
      integer kk,kx,kz
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr3,rr5,rr7
      real*8 ci,dix,diy,diz
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,duir,puir
      real*8 dkr,dukr,pukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 scale3,scale5
      real*8 scale7,damp
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 field(3,maxatm)
      real*8 fieldp(3,maxatm)
      real*8 udir(3,maxatm)
      real*8 udirp(3,maxatm)
      real*8 uold(3,maxatm)
      real*8 uoldp(3,maxatm)
      real*8 dscale(maxatm)
      real*8 pscale(maxatm)
      logical proceed,done
c
c
c     zero out the induced dipole and the field at each site
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
      if (.not. use_polar)  return
c
c     set the switching function coefficients
c
      call switch ('MPOLE')
c
c     compute the direct induced dipole moment at each atom
c
      do i = 1, npole-1
         ii = ipole(i)
         iz = zaxis(ii)
         ix = xaxis(ii)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = i+1, npole
            dscale(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
         end do
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii)) 
     &            pscale(i14(j,ii)) = 0.5d0 * pscale(i14(j,ii))
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            kz = zaxis(kk)
            kx = xaxis(kk)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               call image (xr,yr,zr,0)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdamp(i) * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        scale3 = 1.0d0-exp(damp)
                        scale5 = 1.0d0-(1.0d0-damp)*exp(damp)
                        scale7 = 1.0d0-(1.0d0-damp+0.6d0*damp**2)
     &                                         *exp(damp)
                     end if
                  end if
                  rr3 = scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
                  dir = dix*xr + diy*yr + diz*zr
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qir = qix*xr + qiy*yr + qiz*zr
                  dkr = dkx*xr + dky*yr + dkz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qkr = qkx*xr + qky*yr + qkz*zr
                  fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkx + 2.0d0*rr5*qkx
                  fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dky + 2.0d0*rr5*qky
                  fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkz + 2.0d0*rr5*qkz
                  fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*dix - 2.0d0*rr5*qix
                  fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diy - 2.0d0*rr5*qiy
                  fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diz - 2.0d0*rr5*qiz
                  do j = 1, 3
                     field(j,i) = field(j,i) + fid(j)*dscale(kk)
                     field(j,k) = field(j,k) + fkd(j)*dscale(kk)
                     fieldp(j,i) = fieldp(j,i) + fid(j)*pscale(kk)
                     fieldp(j,k) = fieldp(j,k) + fkd(j)*pscale(kk)
                  end do
               end if
            end if
         end do
      end do
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
         do i = 1, npole
            ii = ipole(i)
            iz = zaxis(ii)
            ix = xaxis(ii)
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            qixx = rpole(5,i)
            qixy = rpole(6,i)
            qixz = rpole(7,i)
            qiyy = rpole(9,i)
            qiyz = rpole(10,i)
            qizz = rpole(13,i)
            do j = i, npole
               dscale(ipole(j)) = 1.0d0
               pscale(ipole(j)) = 1.0d0
            end do
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = p2scale
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = p3scale
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = p4scale
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = p5scale
            end do
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = d1scale
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = d2scale
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = d3scale
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = d4scale
            end do
            do k = i, npole
               kk = ipole(k)
               kz = zaxis(kk)
               kx = xaxis(kk)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
               do m = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call image (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     scale7 = 1.0d0
                     damp = pdamp(i) * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           scale3 = 1.0d0 - exp(damp)
                           scale5 = 1.0d0 - (1.0d0-damp)*exp(damp)
                           scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                             *exp(damp)
                        end if
                     end if
                     rr3 = scale3 / (r*r2)
                     rr5 = 3.0d0 * scale5 / (r*r2*r2)
                     rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
                     dir = dix*xr + diy*yr + diz*zr
                     qix = qixx*xr + qixy*yr + qixz*zr
                     qiy = qixy*xr + qiyy*yr + qiyz*zr
                     qiz = qixz*xr + qiyz*yr + qizz*zr
                     qir = qix*xr + qiy*yr + qiz*zr
                     dkr = dkx*xr + dky*yr + dkz*zr
                     qkx = qkxx*xr + qkxy*yr + qkxz*zr
                     qky = qkxy*xr + qkyy*yr + qkyz*zr
                     qkz = qkxz*xr + qkyz*yr + qkzz*zr
                     qkr = qkx*xr + qky*yr + qkz*zr
                     fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                           - rr3*dkx + 2.0d0*rr5*qkx
                     fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                           - rr3*dky + 2.0d0*rr5*qky
                     fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                           - rr3*dkz + 2.0d0*rr5*qkz
                     fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                           - rr3*dix - 2.0d0*rr5*qix
                     fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                           - rr3*diy - 2.0d0*rr5*qiy
                     fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                           - rr3*diz - 2.0d0*rr5*qiz
                     do j = 1, 3
                       fip(j) = fid(j)
                       fkp(j) = fkd(j)
                     end do
                     if (use_polymer .and. r2 .le. polycut2) then
                        do j = 1, 3
                           fid(j) = fid(j) * dscale(kk)
                           fip(j) = fip(j) * pscale(kk)
                           fkd(j) = fkd(j) * dscale(kk)
                           fkp(j) = fkp(j) * pscale(kk)
                        end do
                     end if
                     do j = 1, 3
                        field(j,i) = field(j,i) + fid(j)
                        fieldp(j,i) = fieldp(j,i) + fip(j)
                        if (ii .ne. kk) then
                           field(j,k) = field(j,k) + fkd(j)
                           fieldp(j,k) = fieldp(j,k) + fkp(j)
                        end if
                     end do
                  end if
               end do
            end do
         end do
      end if
c
c     set induced dipoles to polarizability times direct field
c
      do i = 1, npole
         do j = 1, 3
            udir(j,i) = polarity(i) * field(j,i)
            udirp(j,i) = polarity(i) * fieldp(j,i)
            uind(j,i) = udir(j,i)
            uinp(j,i) = udirp(j,i)
         end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         eps = 1.0d0
c
c     compute mutual induced dipole moments by an iterative method
c
         dowhile (.not. done)
            do i = 1, npole
               do j = 1, 3
                  field(j,i) = 0.0d0
                  fieldp(j,i) = 0.0d0
               end do
            end do
            do i = 1, npole-1
               ii = ipole(i)
               iz = zaxis(ii)
               ix = xaxis(ii)
               duix = uind(1,i)
               duiy = uind(2,i)
               duiz = uind(3,i)
               puix = uinp(1,i)
               puiy = uinp(2,i)
               puiz = uinp(3,i)
               do j = i+1, npole
                  dscale(ipole(j)) = 1.0d0
               end do
               do j = 1, np11(ii)
                  dscale(ip11(j,ii)) = u1scale
               end do
               do j = 1, np12(ii)
                  dscale(ip12(j,ii)) = u2scale
               end do
               do j = 1, np13(ii)
                  dscale(ip13(j,ii)) = u3scale
               end do
               do j = 1, np14(ii)
                  dscale(ip14(j,ii)) = u4scale
               end do
               do k = i+1, npole
                  kk = ipole(k)
                  kz = zaxis(kk)
                  kx = xaxis(kk)
                  proceed = .true.
                  if (use_intra)
     &               call groups (proceed,fgrp,ii,kk,0,0,0,0)
                  if (proceed) then
                     xr = x(kk) - x(ii)
                     yr = y(kk) - y(ii)
                     zr = z(kk) - z(ii)
                     call image (xr,yr,zr,0)
                     r2 = xr*xr + yr* yr + zr*zr
                     if (r2 .le. off2) then
                        r = sqrt(r2)
                        dukx = uind(1,k)
                        duky = uind(2,k)
                        dukz = uind(3,k)
                        pukx = uinp(1,k)
                        puky = uinp(2,k)
                        pukz = uinp(3,k)
                        scale3 = dscale(kk)
                        scale5 = dscale(kk)
                        damp = pdamp(i) * pdamp(k)
                        if (damp .ne. 0.0d0) then
                           damp = -pgamma * (r/damp)**3
                           if (damp .gt. -50.0d0) then
                              scale3 = scale3 * (1.0d0-exp(damp))
                              scale5 = scale5 * (1.0d0-(1.0d0-damp)
     &                                                     *exp(damp))
                           end if
                        end if
                        rr3 = scale3 / (r*r2)
                        rr5 = 3.0d0 * scale5 / (r*r2*r2)
                        duir = xr*duix + yr*duiy + zr*duiz
                        dukr = xr*dukx + yr*duky + zr*dukz
                        puir = xr*puix + yr*puiy + zr*puiz
                        pukr = xr*pukx + yr*puky + zr*pukz
                        fid(1) = -rr3*dukx + rr5*dukr*xr
                        fid(2) = -rr3*duky + rr5*dukr*yr
                        fid(3) = -rr3*dukz + rr5*dukr*zr
                        fkd(1) = -rr3*duix + rr5*duir*xr
                        fkd(2) = -rr3*duiy + rr5*duir*yr
                        fkd(3) = -rr3*duiz + rr5*duir*zr
                        fip(1) = -rr3*pukx + rr5*pukr*xr
                        fip(2) = -rr3*puky + rr5*pukr*yr
                        fip(3) = -rr3*pukz + rr5*pukr*zr
                        fkp(1) = -rr3*puix + rr5*puir*xr
                        fkp(2) = -rr3*puiy + rr5*puir*yr
                        fkp(3) = -rr3*puiz + rr5*puir*zr
                        do j = 1, 3
                           field(j,i) = field(j,i) + fid(j)
                           field(j,k) = field(j,k) + fkd(j)
                           fieldp(j,i) = fieldp(j,i) + fip(j)
                           fieldp(j,k) = fieldp(j,k) + fkp(j)
                        end do
                     end if
                  end if
               end do
            end do
c
c     periodic boundary for large cutoffs via replicates method
c
            if (use_replica) then
               do i = 1, npole
                  ii = ipole(i)
                  iz = zaxis(ii)
                  ix = xaxis(ii)
                  duix = uind(1,i)
                  duiy = uind(2,i)
                  duiz = uind(3,i)
                  puix = uinp(1,i)
                  puiy = uinp(2,i)
                  puiz = uinp(3,i)
                  do j = i, npole
                     dscale(ipole(j)) = 1.0d0
                  end do
                  do j = 1, np11(ii)
                     dscale(ip11(j,ii)) = u1scale
                  end do
                  do j = 1, np12(ii)
                     dscale(ip12(j,ii)) = u2scale
                  end do
                  do j = 1, np13(ii)
                     dscale(ip13(j,ii)) = u3scale
                  end do
                  do j = 1, np14(ii)
                     dscale(ip14(j,ii)) = u4scale
                  end do
                  do k = i, npole
                     kk = ipole(k)
                     kz = zaxis(kk)
                     kx = xaxis(kk)
                     dukx = uind(1,k)
                     duky = uind(2,k)
                     dukz = uind(3,k)
                     pukx = uinp(1,k)
                     puky = uinp(2,k)
                     pukz = uinp(3,k)
                     do m = 1, ncell
                        xr = x(kk) - x(ii)
                        yr = y(kk) - y(ii)
                        zr = z(kk) - z(ii)
                        call image (xr,yr,zr,m)
                        r2 = xr*xr + yr* yr + zr*zr
                        if (r2 .le. off2) then
                           r = sqrt(r2)
                           scale3 = 1.0d0
                           scale5 = 1.0d0
                           damp = pdamp(i) * pdamp(k)
                           if (damp .ne. 0.0d0) then
                              damp = -pgamma * (r/damp)**3
                              if (damp .gt. -50.0d0) then
                                 scale3 = 1.0d0 - exp(damp)
                                 scale5 = 1.0d0 - (1.0d0-damp)*exp(damp)
                              end if
                           end if
                           rr3 = scale3 / (r*r2)
                           rr5 = 3.0d0 * scale5 / (r*r2*r2)
                           duir = xr*duix + yr*duiy + zr*duiz
                           dukr = xr*dukx + yr*duky + zr*dukz
                           puir = xr*puix + yr*puiy + zr*puiz
                           pukr = xr*pukx + yr*puky + zr*pukz
                           fid(1) = -rr3*dukx + rr5*dukr*xr
                           fid(2) = -rr3*duky + rr5*dukr*yr
                           fid(3) = -rr3*dukz + rr5*dukr*zr
                           fkd(1) = -rr3*duix + rr5*duir*xr
                           fkd(2) = -rr3*duiy + rr5*duir*yr
                           fkd(3) = -rr3*duiz + rr5*duir*zr
                           fip(1) = -rr3*pukx + rr5*pukr*xr
                           fip(2) = -rr3*puky + rr5*pukr*yr
                           fip(3) = -rr3*pukz + rr5*pukr*zr
                           fkp(1) = -rr3*puix + rr5*puir*xr
                           fkp(2) = -rr3*puiy + rr5*puir*yr
                           fkp(3) = -rr3*puiz + rr5*puir*zr
                           if (use_polymer) then
                              if (r2 .le. polycut2) then
                                 do j = 1, 3
                                    fid(j) = fid(j) * dscale(kk)
                                    fkd(j) = fkd(j) * dscale(kk)
                                    fip(j) = fip(j) * dscale(kk)
                                    fkp(j) = fkp(j) * dscale(kk)
                                 end do
                              end if
                           end if
                           do j = 1, 3
                              field(j,i) = field(j,i) + fid(j)
                              fieldp(j,i) = fieldp(j,i) + fip(j)
                              if (ii .ne. kk) then
                                 field(j,k) = field(j,k) + fkd(j)
                                 fieldp(j,k) = fieldp(j,k) + fkp(j)
                              end if
                           end do
                        end if
                     end do
                  end do
               end do
            end if
c
c     check to see if the mutual induced dipoles have converged
c
            iter = iter + 1
            epsold = eps
            epsd = 0.0d0
            epsp = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  uold(j,i) = uind(j,i)
                  uoldp(j,i) = uinp(j,i)
                  uind(j,i) = udir(j,i) + polarity(i)*field(j,i)
                  uinp(j,i) = udirp(j,i) + polarity(i)*fieldp(j,i)
                  uind(j,i) = uold(j,i) + polsor*(uind(j,i)-uold(j,i))
                  uinp(j,i) = uoldp(j,i) + polsor*(uinp(j,i)-uoldp(j,i))
                  epsd = epsd + (uind(j,i)-uold(j,i))**2
                  epsp = epsp + (uinp(j,i)-uoldp(j,i))**2
               end do
            end do
            eps = max(epsd,epsp)
            eps = debye * sqrt(eps/dble(npolar))
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. maxiter)  done = .true.
         end do
c
c     terminate the calculation if dipoles failed to converge
c
         if (eps .gt. poleps) then
            write (iout,30)
   30       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine induce0b  --  Ewald summation induced dipoles  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "induce0b" computes the induced dipole moment at each
c     polarizable site using a regular Ewald summation
c
c
      subroutine induce0b
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'ewald.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'potent.i'
      include 'units.i'
      integer i,j,ii
      integer iter,maxiter
      real*8 eps,term
      real*8 epsd,epsp
      real*8 epsold
      real*8 ucell(3)
      real*8 ucellp(3)
      real*8 udir(3,maxatm)
      real*8 udirp(3,maxatm)
      real*8 uold(3,maxatm)
      real*8 uoldp(3,maxatm)
      real*8 field(3,maxatm)
      real*8 fieldp(3,maxatm)
      logical done
c
c
c     zero out the induced dipole and the field at each site
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
      if (.not. use_polar)  return
c
c     get the reciprical space part of the electrostatic field
c
      call udirect1 (field)
c
c     get the real space portion of the electrostatic field
c
      do i = 1, npole
         do j = 1, 3
            fieldp(j,i) = field(j,i)
         end do
      end do
      call udirect2 (field,fieldp)
c
c     get the self-energy portion of the electrostatic field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do i = 1, npole
         do j = 1, 3
            field(j,i) = field(j,i) + term*rpole(j+1,i)
            fieldp(j,i) = fieldp(j,i) + term*rpole(j+1,i)
         end do
      end do
c
c     compute the cell dipole boundary correction to field
c
      if (.not. tinfoil) then
         do i = 1, 3
            ucell(i) = 0.0d0
         end do
         do i = 1, npole
            ii = ipole(i)
            ucell(1) = ucell(1) + rpole(2,i) + rpole(1,i)*x(ii)
            ucell(2) = ucell(2) + rpole(3,i) + rpole(1,i)*y(ii)
            ucell(3) = ucell(3) + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do i = 1, npole
            do j = 1, 3
               field(j,i) = field(j,i) - term*ucell(j)
               fieldp(j,i) = fieldp(j,i) - term*ucell(j)
            end do
         end do
      end if
c
c     set induced dipoles to polarizability times direct field
c
      do i = 1, npole
         do j = 1, 3
            udir(j,i) = polarity(i) * field(j,i)
            udirp(j,i) = polarity(i) * fieldp(j,i)
            uind(j,i) = udir(j,i)
            uinp(j,i) = udirp(j,i)
         end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         eps = 1.0d0
c
c     compute mutual induced dipole moments by an iterative method
c
         dowhile (.not. done)
            do i = 1, npole
               do j = 1, 3
                  field(j,i) = 0.0d0
                  fieldp(j,i) = 0.0d0
               end do
            end do
            call umutual1 (field,fieldp)
            call umutual2 (field,fieldp)
            term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
            do i = 1, npole
               do j = 1, 3
                  field(j,i) = field(j,i) + term*uind(j,i)
                  fieldp(j,i) = fieldp(j,i) + term*uinp(j,i)
               end do
            end do
            if (.not. tinfoil) then
               do i = 1, 3
                  ucell(i) = 0.0d0
                  ucellp(i) = 0.0d0
               end do
               do i = 1, npole
                  do j = 1, 3
                     ucell(j) = ucell(j) + uind(j,i)
                     ucellp(j) = ucellp(j) + uinp(j,i)
                  end do
               end do
               term = (4.0d0/3.0d0) * pi/volbox
               do i = 1, npole
                  do j = 1, 3
                     field(j,i) = field(j,i) - term*ucell(j)
                     fieldp(j,i) = fieldp(j,i) - term*ucellp(j)
                  end do
               end do
            end if
c
c     check to see if the mutual induced dipoles have converged
c
            iter = iter + 1
            epsold = eps
            epsd = 0.0d0
            epsp = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  uold(j,i) = uind(j,i)
                  uoldp(j,i) = uinp(j,i)
                  uind(j,i) = udir(j,i) + polarity(i)*field(j,i)
                  uinp(j,i) = udirp(j,i) + polarity(i)*fieldp(j,i)
                  uind(j,i) = uold(j,i) + polsor*(uind(j,i)-uold(j,i))
                  uinp(j,i) = uoldp(j,i) + polsor*(uinp(j,i)-uoldp(j,i))
                  epsd = epsd + (uind(j,i)-uold(j,i))**2
                  epsp = epsp + (uinp(j,i)-uoldp(j,i))**2
               end do
            end do
            eps = max(epsd,epsp)
            eps = debye * sqrt(eps/dble(npolar))
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. maxiter)  done = .true.
         end do
c
c     terminate the calculation if dipoles failed to converge
c
         if (eps .gt. poleps) then
            write (iout,30)
   30       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine udirect1  --  Ewald recip direct induced field  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "udirect1" computes the reciprocal space contribution of the
c     permanent atomic multipole moments to the electrostatic field
c     for use in finding the direct induced dipole moments via a
c     regular Ewald summation
c
c
      subroutine udirect1 (field)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'ewald.i'
      include 'ewreg.i'
      include 'math.i'
      include 'mpole.i'
      include 'units.i'
      integer i,j,k,l,ii
      integer jmin,jmax
      integer kmin,kmax
      integer lmin,lmax
      real*8 expterm,cut
      real*8 term,fterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3,hsq
      real*8 qf,t1,t2
      real*8 ck,dk,qk
      real*8 q1,q2,q3
      real*8 ckr(maxatm)
      real*8 skr(maxatm)
      real*8 cjk(maxatm)
      real*8 sjk(maxatm)
      real*8 cm(maxatm)
      real*8 dm(3,maxatm)
      real*8 qm(9,maxatm)
      real*8 field(3,maxatm)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      term = -0.25d0 / aewald**2
      fterm = 8.0d0 * pi / volbox
c
c     set the number of vectors based on box dimensions
c
      cut = 4.0d0 * pi * pi * frecip
      jmin = 0
      kmin = 0
      lmin = 1
      jmax = min(maxvec,int(frecip/recip(1,1)))
      kmax = min(maxvec,int(frecip/recip(2,2)))
      lmax = min(maxvec,int(frecip/recip(3,3)))
c
c     copy the multipole moments into local storage areas
c
      do i = 1, npole
         cm(i) = rpole(1,i)
         dm(1,i) = rpole(2,i)
         dm(2,i) = rpole(3,i)
         dm(3,i) = rpole(4,i)
         qm(1,i) = rpole(5,i)
         qm(2,i) = rpole(6,i)
         qm(3,i) = rpole(7,i)
         qm(4,i) = rpole(8,i)
         qm(5,i) = rpole(9,i)
         qm(6,i) = rpole(10,i)
         qm(7,i) = rpole(11,i)
         qm(8,i) = rpole(12,i)
         qm(9,i) = rpole(13,i)
      end do
c
c     calculate and store the exponential factors
c
      do i = 1, npole
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
         ejc(i,0) = 1.0d0
         ejs(i,0) = 0.0d0
         ekc(i,0) = 1.0d0
         eks(i,0) = 0.0d0
         elc(i,0) = 1.0d0
         els(i,0) = 0.0d0
         ejc(i,1) = cos(xfr)
         ejs(i,1) = sin(xfr)
         ekc(i,1) = cos(yfr)
         eks(i,1) = sin(yfr)
         elc(i,1) = cos(zfr)
         els(i,1) = sin(zfr)
         ekc(i,-1) = ekc(i,1)
         eks(i,-1) = -eks(i,1)
         elc(i,-1) = elc(i,1)
         els(i,-1) = -els(i,1)
         do j = 2, jmax
            ejc(i,j) = ejc(i,j-1)*ejc(i,1) - ejs(i,j-1)*ejs(i,1)
            ejs(i,j) = ejs(i,j-1)*ejc(i,1) + ejc(i,j-1)*ejs(i,1)
         end do
         do j = 2, kmax
            ekc(i,j) = ekc(i,j-1)*ekc(i,1) - eks(i,j-1)*eks(i,1)
            eks(i,j) = eks(i,j-1)*ekc(i,1) + ekc(i,j-1)*eks(i,1)
            ekc(i,-j) = ekc(i,j)
            eks(i,-j) = -eks(i,j)
         end do
         do j = 2, lmax
            elc(i,j) = elc(i,j-1)*elc(i,1) - els(i,j-1)*els(i,1)
            els(i,j) = els(i,j-1)*elc(i,1) + elc(i,j-1)*els(i,1)
            elc(i,-j) = elc(i,j)
            els(i,-j) = -els(i,j)
         end do
      end do
c
c     loop over all k vectors from the reciprocal lattice
c
      do j = jmin, jmax
         rj = 2.0d0 * pi * dble(j)
         do k = kmin, kmax
            rk = 2.0d0 * pi * dble(k)
            do i = 1, npole
               cjk(i) = ejc(i,j)*ekc(i,k) - ejs(i,j)*eks(i,k)
               sjk(i) = ejs(i,j)*ekc(i,k) + ejc(i,j)*eks(i,k)
            end do
            do l = lmin, lmax
               rl = 2.0d0 * pi * dble(l)
               h1 = recip(1,1)*rj
               h2 = recip(2,1)*rj + recip(2,2)*rk
               h3 = recip(3,1)*rj + recip(3,2)*rk + recip(3,3)*rl
               hsq = h1*h1 + h2*h2 + h3*h3
               if (hsq .le. cut) then
                  t1 = 0.0d0
                  t2 = 0.0d0
                  do i = 1, npole
                     ckr(i) = cjk(i)*elc(i,l) - sjk(i)*els(i,l)
                     skr(i) = sjk(i)*elc(i,l) + cjk(i)*els(i,l)
                     ck = cm(i)
                     dk = h1*dm(1,i) + h2*dm(2,i) + h3*dm(3,i)
                     q1 = h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i)
                     q2 = h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i)
                     q3 = h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i)
                     qk = h1*q1 + h2*q2 + h3*q3
                     t1 = t1 + (ck-qk)*skr(i) + dk*ckr(i)
                     t2 = t2 + (ck-qk)*ckr(i) - dk*skr(i)
                  end do
                  expterm = fterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0 
                  end if
                  do i = 1, npole
                     qf = expterm * (skr(i)*t2-ckr(i)*t1)
                     field(1,i) = field(1,i) + h1*qf
                     field(2,i) = field(2,i) + h2*qf
                     field(3,i) = field(3,i) + h3*qf
                  end do
               end if
            end do
            lmin = -lmax
         end do
         kmin = -kmax
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine udirect2  --  Ewald real direct induced field  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "udirect2" computes the real space contribution of the permanent
c     atomic multipole moments to the electrostatic field for use in
c     finding the direct induced dipole moments via a regular Ewald
c     summation
c
c
      subroutine udirect2 (field,fieldp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'ewald.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'units.i'
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr,r,r2
      real*8 erfc,bfac,exp2a
      real*8 drr3,drr5,drr7
      real*8 prr3,prr5,prr7
      real*8 ci,dix,diy,diz
      real*8 qixx,qiyy,qizz
      real*8 qixy,qixz,qiyz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkyy,qkzz
      real*8 qkxy,qkxz,qkyz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 ralpha,damp
      real*8 alsq2,alsq2n
      real*8 scale3,scale5
      real*8 scale7
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 bn(0:3)
      real*8 fim(3),fkm(3)
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 field(3,maxatm)
      real*8 fieldp(3,maxatm)
      real*8 dscale(maxatm)
      real*8 pscale(maxatm)
      external erfc
c
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
      call switch ('EWALD')
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole-1
         ii = ipole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = i+1, npole
            dscale(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
         end do
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii)) 
     &            pscale(i14(j,ii)) = 0.5d0 * pscale(i14(j,ii))
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr,0)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. cut2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     calculate the error function damping terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)
     &            alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 3
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     compute the error function scaled and unscaled terms
c
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               damp = pdamp(i) * pdamp(k)
               if (damp .ne. 0.0d0) then
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     scale3 = 1.0d0-exp(damp)
                     scale5 = 1.0d0-(1.0d0-damp)*exp(damp)
                     scale7 = 1.0d0-(1.0d0-damp+0.6d0*damp**2)*exp(damp)
                  end if
               end if
               dsc3=scale3*dscale(kk)
               dsc5=scale5*dscale(kk)
               dsc7=scale7*dscale(kk)
               psc3=scale3*pscale(kk)
               psc5=scale5*pscale(kk)
               psc7=scale7*pscale(kk)
               drr3 = (1.0d0-dsc3) / (r*r2)
               drr5 = 3.0d0 * (1.0d0-dsc5) / (r*r2*r2)
               drr7 = 15.0d0 * (1.0d0-dsc7) / (r*r2*r2*r2)
               prr3 = (1.0d0-psc3) / (r*r2)
               prr5 = 3.0d0 * (1.0d0-psc5) / (r*r2*r2)
               prr7 = 15.0d0 * (1.0d0-psc7) / (r*r2*r2*r2)
               dir = dix*xr + diy*yr + diz*zr
               qix = qixx*xr + qixy*yr + qixz*zr
               qiy = qixy*xr + qiyy*yr + qiyz*zr
               qiz = qixz*xr + qiyz*yr + qizz*zr
               qir = qix*xr + qiy*yr + qiz*zr
               dkr = dkx*xr + dky*yr + dkz*zr
               qkx = qkxx*xr + qkxy*yr + qkxz*zr
               qky = qkxy*xr + qkyy*yr + qkyz*zr
               qkz = qkxz*xr + qkyz*yr + qkzz*zr
               qkr = qkx*xr + qky*yr + qkz*zr
               fim(1) = -xr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                     - bn(1)*dkx + 2.0d0*bn(2)*qkx
               fim(2) = -yr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                     - bn(1)*dky + 2.0d0*bn(2)*qky
               fim(3) = -zr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                     - bn(1)*dkz + 2.0d0*bn(2)*qkz
               fkm(1) = xr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                     - bn(1)*dix - 2.0d0*bn(2)*qix
               fkm(2) = yr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                     - bn(1)*diy - 2.0d0*bn(2)*qiy
               fkm(3) = zr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                     - bn(1)*diz - 2.0d0*bn(2)*qiz
               fid(1) = -xr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                     - drr3*dkx + 2.0d0*drr5*qkx
               fid(2) = -yr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                     - drr3*dky + 2.0d0*drr5*qky
               fid(3) = -zr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                     - drr3*dkz + 2.0d0*drr5*qkz
               fkd(1) = xr*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*dix - 2.0d0*drr5*qix
               fkd(2) = yr*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*diy - 2.0d0*drr5*qiy
               fkd(3) = zr*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*diz - 2.0d0*drr5*qiz
               fip(1) = -xr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                     - prr3*dkx + 2.0d0*prr5*qkx
               fip(2) = -yr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                     - prr3*dky + 2.0d0*prr5*qky
               fip(3) = -zr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                     - prr3*dkz + 2.0d0*prr5*qkz
               fkp(1) = xr*(prr3*ci+prr5*dir+prr7*qir)
     &                     - prr3*dix - 2.0d0*prr5*qix
               fkp(2) = yr*(prr3*ci+prr5*dir+prr7*qir)
     &                     - prr3*diy - 2.0d0*prr5*qiy
               fkp(3) = zr*(prr3*ci+prr5*dir+prr7*qir)
     &                     - prr3*diz - 2.0d0*prr5*qiz
c
c     increment the field at each site due to this interaction
c
               do j = 1, 3
                  field(j,i) = field(j,i) + fim(j) - fid(j)
                  field(j,k) = field(j,k) + fkm(j) - fkd(j)
                  fieldp(j,i) = fieldp(j,i) + fim(j) - fip(j)
                  fieldp(j,k) = fieldp(j,k) + fkm(j) - fkp(j)
               end do
            end if
         end do
      end do
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
         do i = 1, npole
            ii = ipole(i)
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            qixx = rpole(5,i)
            qixy = rpole(6,i)
            qixz = rpole(7,i)
            qiyy = rpole(9,i)
            qiyz = rpole(10,i)
            qizz = rpole(13,i)
            do j = i, npole
               dscale(ipole(j)) = 1.0d0
               pscale(ipole(j)) = 1.0d0
            end do
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = p2scale
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = p3scale
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = p4scale
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = p5scale
            end do
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = d1scale
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = d2scale
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = d3scale
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = d4scale
            end do
            do k = i, npole
               kk = ipole(k)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
               do m = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call image (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
c
c     calculate the error function damping terms
c
                  if (r2 .le. cut2) then
                     r = sqrt(r2)
                     ralpha = aewald * r
                     bn(0) = erfc(ralpha) / r
                     alsq2 = 2.0d0 * aewald**2
                     alsq2n = 0.0d0
                     if (aewald .gt. 0.0d0)
     &                  alsq2n = 1.0d0 / (sqrtpi*aewald)
                     exp2a = exp(-ralpha**2)
                     do j = 1, 3
                        bfac = dble(j+j-1)
                        alsq2n = alsq2 * alsq2n
                        bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
                     end do
c
c     compute the error function scaled and unscaled terms
c
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     scale7 = 1.0d0
                     damp = pdamp(i) * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           scale3 = 1.0d0 - exp(damp)
                           scale5 = 1.0d0 - (1.0d0-damp)*exp(damp)
                           scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                             *exp(damp)
                        end if
                     end if
                     dsc3 = scale3
                     dsc5 = scale5
                     dsc7 = scale7
                     psc3 = scale3
                     psc5 = scale5
                     psc7 = scale7
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           dsc3 = scale3 * dscale(kk)
                           dsc5 = scale5 * dscale(kk)
                           dsc7 = scale7 * dscale(kk)
                           psc3 = scale3 * pscale(kk)
                           psc5 = scale5 * pscale(kk)
                           psc7 = scale7 * pscale(kk)
                        end if
                     end if
                     drr3 = (1.0d0-dsc3) / (r*r2)
                     drr5 = 3.0d0 * (1.0d0-dsc5) / (r*r2*r2)
                     drr7 = 15.0d0 * (1.0d0-dsc7) / (r*r2*r2*r2)
                     prr3 = (1.0d0-psc3) / (r*r2)
                     prr5 = 3.0d0 * (1.0d0-psc5) / (r*r2*r2)
                     prr7 = 15.0d0 * (1.0d0-psc7) / (r*r2*r2*r2)
                     dir = dix*xr + diy*yr + diz*zr
                     qix = qixx*xr + qixy*yr + qixz*zr
                     qiy = qixy*xr + qiyy*yr + qiyz*zr
                     qiz = qixz*xr + qiyz*yr + qizz*zr
                     qir = qix*xr + qiy*yr + qiz*zr
                     dkr = dkx*xr + dky*yr + dkz*zr
                     qkx = qkxx*xr + qkxy*yr + qkxz*zr
                     qky = qkxy*xr + qkyy*yr + qkyz*zr
                     qkz = qkxz*xr + qkyz*yr + qkzz*zr
                     qkr = qkx*xr + qky*yr + qkz*zr
                     fim(1) = -xr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                           - bn(1)*dkx + 2.0d0*bn(2)*qkx
                     fim(2) = -yr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                           - bn(1)*dky + 2.0d0*bn(2)*qky
                     fim(3) = -zr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                           - bn(1)*dkz + 2.0d0*bn(2)*qkz
                     fkm(1) = xr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                           - bn(1)*dix - 2.0d0*bn(2)*qix
                     fkm(2) = yr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                           - bn(1)*diy - 2.0d0*bn(2)*qiy
                     fkm(3) = zr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                           - bn(1)*diz - 2.0d0*bn(2)*qiz
                     fid(1) = -xr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                           - drr3*dkx + 2.0d0*drr5*qkx
                     fid(2) = -yr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                           - drr3*dky + 2.0d0*drr5*qky
                     fid(3) = -zr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                           - drr3*dkz + 2.0d0*drr5*qkz
                     fkd(1) = xr*(drr3*ci+drr5*dir+drr7*qir)
     &                           - drr3*dix - 2.0d0*drr5*qix
                     fkd(2) = yr*(drr3*ci+drr5*dir+drr7*qir)
     &                           - drr3*diy - 2.0d0*drr5*qiy
                     fkd(3) = zr*(drr3*ci+drr5*dir+drr7*qir)
     &                           - drr3*diz - 2.0d0*drr5*qiz
                     fip(1) = -xr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                           - prr3*dkx + 2.0d0*prr5*qkx
                     fip(2) = -yr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                           - prr3*dky + 2.0d0*prr5*qky
                     fip(3) = -zr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                           - prr3*dkz + 2.0d0*prr5*qkz
                     fkp(1) = xr*(prr3*ci+prr5*dir+prr7*qir)
     &                           - prr3*dix - 2.0d0*prr5*qix
                     fkp(2) = yr*(prr3*ci+prr5*dir+prr7*qir)
     &                           - prr3*diy - 2.0d0*prr5*qiy
                     fkp(3) = zr*(prr3*ci+prr5*dir+prr7*qir)
     &                           - prr3*diz - 2.0d0*prr5*qiz
c
c     increment the field at each site due to this interaction
c
                     do j = 1, 3
                        field(j,i) = field(j,i) + fim(j) - fid(j)
                        fieldp(j,i) = fieldp(j,i) + fim(j) - fip(j)
                        if (ii .ne. kk) then
                           field(j,k) = field(j,k) + fkm(j) - fkd(j)
                           fieldp(j,k) = fieldp(j,k) + fkm(j) - fkp(j)
                        end if
                     end do
                  end if
               end do
            end do
         end do
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine umutual1  --  Ewald recip mutual induced field  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "umutual1" computes the reciprocal space contribution of the
c     induced atomic dipole moments to the electrostatic field for
c     use in iterative calculation of induced dipole moments via a
c     regular Ewald summation
c
c
      subroutine umutual1 (field,fieldp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'ewald.i'
      include 'ewreg.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'units.i'
      integer i,j,k,l,ii
      integer jmin,jmax
      integer kmin,kmax
      integer lmin,lmax
      real*8 expterm,cut
      real*8 term,fterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3,hsq
      real*8 duk,puk
      real*8 dqf,pqf
      real*8 t1,t2,t3,t4
      real*8 ckr(maxatm)
      real*8 skr(maxatm)
      real*8 cjk(maxatm)
      real*8 sjk(maxatm)
      real*8 field(3,maxatm)
      real*8 fieldp(3,maxatm)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      term = -0.25d0 / aewald**2
      fterm = 8.0d0 * pi / volbox
c
c     set the number of vectors based on box dimensions
c
      cut = 4.0d0 * pi * pi * frecip
      jmin = 0
      kmin = 0
      lmin = 1
      jmax = min(maxvec,int(frecip/recip(1,1)))
      kmax = min(maxvec,int(frecip/recip(2,2)))
      lmax = min(maxvec,int(frecip/recip(3,3)))
c
c     calculate and store the exponential factors
c
      do i = 1, npole
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
         ejc(i,0) = 1.0d0
         ejs(i,0) = 0.0d0
         ekc(i,0) = 1.0d0
         eks(i,0) = 0.0d0
         elc(i,0) = 1.0d0
         els(i,0) = 0.0d0
         ejc(i,1) = cos(xfr)
         ejs(i,1) = sin(xfr)
         ekc(i,1) = cos(yfr)
         eks(i,1) = sin(yfr)
         elc(i,1) = cos(zfr)
         els(i,1) = sin(zfr)
         ekc(i,-1) = ekc(i,1)
         eks(i,-1) = -eks(i,1)
         elc(i,-1) = elc(i,1)
         els(i,-1) = -els(i,1)
         do j = 2, jmax
            ejc(i,j) = ejc(i,j-1)*ejc(i,1) - ejs(i,j-1)*ejs(i,1)
            ejs(i,j) = ejs(i,j-1)*ejc(i,1) + ejc(i,j-1)*ejs(i,1)
         end do
         do j = 2, kmax
            ekc(i,j) = ekc(i,j-1)*ekc(i,1) - eks(i,j-1)*eks(i,1)
            eks(i,j) = eks(i,j-1)*ekc(i,1) + ekc(i,j-1)*eks(i,1)
            ekc(i,-j) = ekc(i,j)
            eks(i,-j) = -eks(i,j)
         end do
         do j = 2, lmax
            elc(i,j) = elc(i,j-1)*elc(i,1) - els(i,j-1)*els(i,1)
            els(i,j) = els(i,j-1)*elc(i,1) + elc(i,j-1)*els(i,1)
            elc(i,-j) = elc(i,j)
            els(i,-j) = -els(i,j)
         end do
      end do
c
c     loop over all k vectors from the reciprocal lattice
c
      do j = jmin, jmax
         rj = 2.0d0 * pi * dble(j)
         do k = kmin, kmax
            rk = 2.0d0 * pi * dble(k)
            do i = 1, npole
               cjk(i) = ejc(i,j)*ekc(i,k) - ejs(i,j)*eks(i,k)
               sjk(i) = ejs(i,j)*ekc(i,k) + ejc(i,j)*eks(i,k)
            end do
            do l = lmin, lmax
               rl = 2.0d0 * pi * dble(l)
               h1 = recip(1,1)*rj
               h2 = recip(2,1)*rj + recip(2,2)*rk
               h3 = recip(3,1)*rj + recip(3,2)*rk + recip(3,3)*rl
               hsq = h1*h1 + h2*h2 + h3*h3
               if (hsq .le. cut) then
                  t1 = 0.0d0
                  t2 = 0.0d0
                  t3 = 0.0d0
                  t4 = 0.0d0
                  do i = 1, npole
                     ckr(i) = cjk(i)*elc(i,l) - sjk(i)*els(i,l)
                     skr(i) = sjk(i)*elc(i,l) + cjk(i)*els(i,l)
                     duk = h1*uind(1,i) + h2*uind(2,i) + h3*uind(3,i)
                     puk = h1*uinp(1,i) + h2*uinp(2,i) + h3*uinp(3,i)
                     t1 = t1 + duk*ckr(i)
                     t2 = t2 - duk*skr(i)
                     t3 = t3 + puk*ckr(i)
                     t4 = t4 - puk*skr(i)
                  end do
                  expterm = fterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0 
                  end if
                  do i = 1, npole
                     dqf = expterm * (skr(i)*t2-ckr(i)*t1)
                     pqf = expterm * (skr(i)*t4-ckr(i)*t3)
                     field(1,i) = field(1,i) + h1*dqf
                     field(2,i) = field(2,i) + h2*dqf
                     field(3,i) = field(3,i) + h3*dqf
                     fieldp(1,i) = fieldp(1,i) + h1*pqf
                     fieldp(2,i) = fieldp(2,i) + h2*pqf
                     fieldp(3,i) = fieldp(3,i) + h3*pqf
                  end do
               end if
            end do
            lmin = -lmax
         end do
         kmin = -kmax
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine umutual2  --  Ewald real mutual induced field  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "umutual2" computes the real space contribution of the induced
c     atomic dipole moments to the electrostatic field for use in
c     iterative calculation of induced dipole moments via a regular
c     Ewald summation
c
c
      subroutine umutual2 (field,fieldp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'ewald.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'units.i'
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 erfc,bfac,exp2a
      real*8 duir,dukr
      real*8 puir,pukr
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 ralpha,damp
      real*8 alsq2,alsq2n
      real*8 scale3,scale5
      real*8 bn(0:2)
      real*8 fimd(3),fkmd(3)
      real*8 fimp(3),fkmp(3)
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 field(3,maxatm)
      real*8 fieldp(3,maxatm)
      real*8 dscale(maxatm)
      external erfc
c
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
      call switch ('EWALD')
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole-1
         ii = ipole(i)
         duix = uind(1,i)
         duiy = uind(2,i)
         duiz = uind(3,i)
         puix = uinp(1,i)
         puiy = uinp(2,i)
         puiz = uinp(3,i)
         do j = i+1, npole
            dscale(ipole(j)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = u4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr,0)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. cut2) then
               r = sqrt(r2)
               dukx = uind(1,k)
               duky = uind(2,k)
               dukz = uind(3,k)
               pukx = uinp(1,k)
               puky = uinp(2,k)
               pukz = uinp(3,k)
c
c     calculate the error function damping terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)
     &            alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 2
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     compute the error function scaled and unscaled terms
c
               scale3 = dscale(kk)
               scale5 = dscale(kk)
               damp = pdamp(i) * pdamp(k)
               if (damp .ne. 0.0d0) then
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     scale3 = scale3 * (1.0d0-exp(damp))
                     scale5 = scale5 * (1.0d0-(1.0d0-damp)*exp(damp))
                  end if
               end if
               rr3 = (1.0d0-scale3) / (r*r2)
               rr5 = 3.0d0 * (1.0d0-scale5) / (r*r2*r2)
               duir = xr*duix + yr*duiy + zr*duiz
               dukr = xr*dukx + yr*duky + zr*dukz
               puir = xr*puix + yr*puiy + zr*puiz
               pukr = xr*pukx + yr*puky + zr*pukz
               fimd(1) = -bn(1)*dukx + bn(2)*dukr*xr
               fimd(2) = -bn(1)*duky + bn(2)*dukr*yr
               fimd(3) = -bn(1)*dukz + bn(2)*dukr*zr
               fkmd(1) = -bn(1)*duix + bn(2)*duir*xr
               fkmd(2) = -bn(1)*duiy + bn(2)*duir*yr
               fkmd(3) = -bn(1)*duiz + bn(2)*duir*zr
               fimp(1) = -bn(1)*pukx + bn(2)*pukr*xr
               fimp(2) = -bn(1)*puky + bn(2)*pukr*yr
               fimp(3) = -bn(1)*pukz + bn(2)*pukr*zr
               fkmp(1) = -bn(1)*puix + bn(2)*puir*xr
               fkmp(2) = -bn(1)*puiy + bn(2)*puir*yr
               fkmp(3) = -bn(1)*puiz + bn(2)*puir*zr
               fid(1) = -rr3*dukx + rr5*dukr*xr
               fid(2) = -rr3*duky + rr5*dukr*yr
               fid(3) = -rr3*dukz + rr5*dukr*zr
               fkd(1) = -rr3*duix + rr5*duir*xr
               fkd(2) = -rr3*duiy + rr5*duir*yr
               fkd(3) = -rr3*duiz + rr5*duir*zr
               fip(1) = -rr3*pukx + rr5*pukr*xr
               fip(2) = -rr3*puky + rr5*pukr*yr
               fip(3) = -rr3*pukz + rr5*pukr*zr
               fkp(1) = -rr3*puix + rr5*puir*xr
               fkp(2) = -rr3*puiy + rr5*puir*yr
               fkp(3) = -rr3*puiz + rr5*puir*zr
c
c     increment the field at each site due to this interaction
c
               do j = 1, 3
                  field(j,i) = field(j,i) + fimd(j) - fid(j)
                  field(j,k) = field(j,k) + fkmd(j) - fkd(j)
                  fieldp(j,i) = fieldp(j,i) + fimp(j) - fip(j)
                  fieldp(j,k) = fieldp(j,k) + fkmp(j) - fkp(j)
               end do
            end if
         end do
      end do
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
         do i = 1, npole
            ii = ipole(i)
            duix = uind(1,i)
            duiy = uind(2,i)
            duiz = uind(3,i)
            puix = uinp(1,i)
            puiy = uinp(2,i)
            puiz = uinp(3,i)
            do j = i, npole
               dscale(ipole(j)) = 1.0d0
            end do
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = u1scale
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = u2scale
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = u3scale
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = u4scale
            end do
            do k = i, npole
               kk = ipole(k)
               dukx = uind(1,k)
               duky = uind(2,k)
               dukz = uind(3,k)
               pukx = uinp(1,k)
               puky = uinp(2,k)
               pukz = uinp(3,k)
               do m = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call image (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
c
c     calculate the error function damping terms
c
                  if (r2 .le. cut2) then
                     r = sqrt(r2)
                     ralpha = aewald * r
                     bn(0) = erfc(ralpha) / r
                     alsq2 = 2.0d0 * aewald**2
                     alsq2n = 0.0d0
                     if (aewald .gt. 0.0d0)
     &                  alsq2n = 1.0d0 / (sqrtpi*aewald)
                     exp2a = exp(-ralpha**2)
                     do j = 1, 2
                        bfac = dble(j+j-1)
                        alsq2n = alsq2 * alsq2n
                        bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
                     end do
c
c     compute the error function scaled and unscaled terms
c
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     damp = pdamp(i) * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           scale3 = 1.0d0 - exp(damp)
                           scale5 = 1.0d0 - (1.0d0-damp)*exp(damp)
                        end if
                     end if
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           scale3 = scale3 * dscale(kk)
                           scale5 = scale5 * dscale(kk)
                        end if
                     end if
                     rr3 = (1.0d0-scale3) / (r*r2)
                     rr5 = 3.0d0 * (1.0d0-scale5) / (r*r2*r2)
                     duir = xr*duix + yr*duiy + zr*duiz
                     dukr = xr*dukx + yr*duky + zr*dukz
                     puir = xr*puix + yr*puiy + zr*puiz
                     pukr = xr*pukx + yr*puky + zr*pukz
                     fimd(1) = -bn(1)*dukx + bn(2)*dukr*xr
                     fimd(2) = -bn(1)*duky + bn(2)*dukr*yr
                     fimd(3) = -bn(1)*dukz + bn(2)*dukr*zr
                     fkmd(1) = -bn(1)*duix + bn(2)*duir*xr
                     fkmd(2) = -bn(1)*duiy + bn(2)*duir*yr
                     fkmd(3) = -bn(1)*duiz + bn(2)*duir*zr
                     fimp(1) = -bn(1)*pukx + bn(2)*pukr*xr
                     fimp(2) = -bn(1)*puky + bn(2)*pukr*yr
                     fimp(3) = -bn(1)*pukz + bn(2)*pukr*zr
                     fkmp(1) = -bn(1)*puix + bn(2)*puir*xr
                     fkmp(2) = -bn(1)*puiy + bn(2)*puir*yr
                     fkmp(3) = -bn(1)*puiz + bn(2)*puir*zr
                     fid(1) = -rr3*dukx + rr5*dukr*xr
                     fid(2) = -rr3*duky + rr5*dukr*yr
                     fid(3) = -rr3*dukz + rr5*dukr*zr
                     fkd(1) = -rr3*duix + rr5*duir*xr
                     fkd(2) = -rr3*duiy + rr5*duir*yr
                     fkd(3) = -rr3*duiz + rr5*duir*zr
                     fip(1) = -rr3*pukx + rr5*pukr*xr
                     fip(2) = -rr3*puky + rr5*pukr*yr
                     fip(3) = -rr3*pukz + rr5*pukr*zr
                     fkp(1) = -rr3*puix + rr5*puir*xr
                     fkp(2) = -rr3*puiy + rr5*puir*yr
                     fkp(3) = -rr3*puiz + rr5*puir*zr
c
c     increment the field at each site due to this interaction
c
                     do j = 1, 3
                        field(j,i) = field(j,i) + fimd(j) - fid(j)
                        fieldp(j,i) = fieldp(j,i) + fimp(j) - fip(j)
                        if (ii .ne. kk) then
                           field(j,k) = field(j,k) + fkmd(j) -fkd(j)
                           fieldp(j,k) = fieldp(j,k) + fkmp(j) -fkd(j)
                        end if
                     end do
                  end if
               end do
            end do
         end do
      end if
      return
      end
