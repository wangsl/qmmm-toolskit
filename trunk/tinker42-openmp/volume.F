 
C     $Id: volume.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine volume  --  excluded volume term via Connolly  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "volume" calculates the excluded volume via the Connolly
c     analytical volume and surface area algorithm
c
c
      subroutine volume (volume_tot,radius,exclude)
      implicit none
      include 'sizes.i'
      real*8 exclude,probe
      real*8 volume_tot,area_tot
      real*8 radius(maxatm)
c
c
c     make call to the volume and surface area routine
c
      probe = 0.0d0
      call connolly (volume_tot,area_tot,radius,probe,exclude)
      return
      end
