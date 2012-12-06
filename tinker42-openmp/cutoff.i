 
C     $Id: cutoff.i 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  cutoff.i  --  cutoff distances for energy interactions  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     vdwcut      cutoff distance for van der Waals interactions
c     chgcut      cutoff distance for charge-charge interactions
c     dplcut      cutoff distance for dipole-dipole interactions
c     mpolecut    cutoff distance for atomic multipole interactions
c     vdwtaper    distance at which van der Waals switching begins
c     chgtaper    distance at which charge-charge switching begins
c     dpltaper    distance at which dipole-dipole switching begins
c     mpoletaper  distance at which atomic multipole switching begins
c     ewaldcut    cutoff distance for direct space Ewald summation
c     use_ewald   logical flag governing use of Ewald summation term
c     use_lights  logical flag to use method of lights neighbors
c
c
      real*8 vdwcut,chgcut
      real*8 dplcut,mpolecut
      real*8 vdwtaper,chgtaper
      real*8 dpltaper,mpoletaper
      real*8 ewaldcut
      logical use_ewald,use_lights
      common /cutoff/ vdwcut,chgcut,dplcut,mpolecut,vdwtaper,chgtaper,
     &                dpltaper,mpoletaper,ewaldcut,use_ewald,use_lights
