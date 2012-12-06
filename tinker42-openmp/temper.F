 
C     $Id: temper.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     #############################################################
c     ##  COPYRIGHT (C) 2003 by Alan Grossfield & Jay W. Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine temper  --  apply thermostat at half time step  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "temper" applies a velocity correction at the half time step
c     as needed for the Nose-Hoover extended system thermostat
c
c     literature references:
c
c     D. Frenkel and B. Smit, "Understanding Molecular Simulation,
c     2nd Edition", Academic Press, San Diego, CA, 2002; see Appendix
c     E.2 for implementation details
c
c     G. J. Martyna, M. E. Tuckerman, D. J. Tobias and M. L. Klein,
c     "Explicit Reversible Integrators for Extended Systems Dynamics",
c     Molecular Physics, 87, 1117-1157 (1996)
c
c
      subroutine temper (dt)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'mdstuf.i'
      include 'moldyn.i'
      include 'units.i'
      integer i,j
      real*8 dt,eksum
      real*8 ekt,scale
      real*8 dt2,dt4,dt8
      real*8 ekin(3,3)
c
c
c     make half-step velocity correction for Nose-Hoover system
c
      if (thermostat .eq. 'NOSE-HOOVER') then
         call kinetic (eksum,ekin)
         dt2 = dt / 2.0d0
         dt4 = dt / 4.0d0
         dt8 = dt / 8.0d0
         ekt = gasconst * kelvin
         gnh(2) = (qnh(1)*vnh(1)*vnh(1)-ekt) / qnh(2)
         vnh(2) = vnh(2) + gnh(2)*dt4
         vnh(1) = vnh(1) * exp(-vnh(2)*dt8)
         gnh(1) = (2.0d0*eksum-dble(nfree)*ekt) / qnh(1)
         vnh(1) = vnh(1) + gnh(1)*dt4
         vnh(1) = vnh(1) * exp(-vnh(2)*dt8)
         xnh(1) = xnh(1) + vnh(1)*dt2
         xnh(2) = xnh(2) + vnh(2)*dt2
         scale = exp(-vnh(1)*dt2)
C$OMP Parallel Do
C$OMP& Default(Shared), Private(I, J)
         do i = 1, n
            do j = 1, 3
               v(j,i) = scale * v(j,i)
            end do
         end do
C$OMP End Parallel Do
         eksum = eksum * scale * scale
         vnh(1) = vnh(1) * exp(-vnh(2)*dt8)
         gnh(1) = (2.0d0*eksum-dble(nfree)*ekt) / qnh(1)
         vnh(1) = vnh(1) + gnh(1)*dt4
         vnh(1) = vnh(1) * exp(-vnh(2)*dt8)
         gnh(2) = (qnh(1)*vnh(1)*vnh(1)-ekt) / qnh(2)
         vnh(2) = vnh(2) + gnh(2)*dt4
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine temper2  --  apply thermostat at full time step  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "temper2" computes the instantaneous temperature and applies a
c     thermostat via Berendsen velocity scaling, Andersen stochastic
c     collisions, Langevin piston or Nose-Hoover extended systems
c
c     literature references:
c
c     H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren,
c     A. DiNola and J. R. Hauk, "Molecular Dynamics with Coupling
c     to an External Bath", Journal of Chemical Physics, 81,
c     3684-3690 (1984)
c
c     H. C. Andersen, "Molecular Dynamics Simulations at Constant
c     Pressure and/or Temperature", Journal of Chemical Physics,
c     72, 2384-2393 (1980)
c
c
      subroutine temper2 (dt,eksum,temp)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'group.i'
      include 'mdstuf.i'
      include 'moldyn.i'
      include 'rgddyn.i'
      include 'units.i'
      include 'usage.i'
      integer i,j
      real*8 dt,eksum,temp
      real*8 scale,speed
      real*8 maxwell,random
      real*8 rate,trial
      real*8 dt2,dt4,dt8,ekt
      real*8 vec(3)
c
c
c     get the instantaneous temperature from the kinetic energy
c
      temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
      if (.not. isothermal)  return
c
c     couple to external temperature bath via Berendsen scaling
c
      if (thermostat .eq. 'BERENDSEN') then
         if (temp .eq. 0.0d0)  temp = 0.1d0
         scale = sqrt(1.0d0 + (dt/tautemp)*(kelvin/temp-1.0d0))
         if (integrate .eq. 'RIGIDBODY') then
            do i = 1, ngrp
               do j = 1, 3
                  vcm(j,i) = scale * vcm(j,i)
                  wcm(j,i) = scale * wcm(j,i)
               end do
            end do
         else
C$OMP Parallel Do
C$OMP& Default(Shared), Private(I, J)
            do i = 1, n
               if (use(i)) then
                  do j = 1, 3
                     v(j,i) = scale * v(j,i)
                  end do
               end if
            end do
C$OMP End Parallel Do
         end if
c
c     select random velocities via Andersen stochastic collisions
c
      else if (thermostat .eq. 'ANDERSEN') then
         rate = 1000.0d0 * collide * dt
         if (integrate .eq. 'RIGIDBODY') then
            rate = rate / dble(ngrp)**(2.0d0/3.0d0)
            do i = 1, ngrp
               trial = random ()
               if (trial .lt. rate) then
                  speed = maxwell (grpmass(i),kelvin)
                  call ranvec (vec)
                  do j = 1, 3
                     vcm(j,i) = speed * vec(j)
                  end do
               end if
            end do
         else
            rate = rate / dble(nuse)**(2.0d0/3.0d0)
            do i = 1, n
               if (use(i)) then
                  trial = random ()
                  if (trial .lt. rate) then
                     speed = maxwell (mass(i),kelvin)
                     call ranvec (vec)
                     do j = 1, 3
                        v(j,i) = speed * vec(j)
                     end do
                  end if
               end if
            end do
         end if
c
c     make full-step velocity correction for Nose-Hoover system
c
      else if (thermostat .eq. 'NOSE-HOOVER') then
         dt2 = dt / 2.0d0
         dt4 = dt / 4.0d0
         dt8 = dt / 8.0d0
         ekt = gasconst * kelvin
         gnh(2) = (qnh(1)*vnh(1)*vnh(1)-ekt) / qnh(2)
         vnh(2) = vnh(2) + gnh(2)*dt4
         vnh(1) = vnh(1) * exp(-vnh(2)*dt8)
         gnh(1) = (2.0d0*eksum-dble(nfree)*ekt) / qnh(1)
         vnh(1) = vnh(1) + gnh(1)*dt4
         vnh(1) = vnh(1) * exp(-vnh(2)*dt8)
         xnh(1) = xnh(1) + vnh(1)*dt2
         xnh(2) = xnh(2) + vnh(2)*dt2
         scale = exp(-vnh(1)*dt2)
         do i = 1, n
            do j = 1, 3
               v(j,i) = scale * v(j,i)
            end do
         end do
         eksum = eksum * scale * scale
         vnh(1) = vnh(1) * exp(-vnh(2)*dt8)
         gnh(1) = (2.0d0*eksum-dble(nfree)*ekt) / qnh(1)
         vnh(1) = vnh(1) + gnh(1)*dt4
         vnh(1) = vnh(1) * exp(-vnh(2)*dt8)
         gnh(2) = (qnh(1)*vnh(1)*vnh(1)-ekt) / qnh(2)
         vnh(2) = vnh(2) + gnh(2)*dt4
      end if
      return
      end
