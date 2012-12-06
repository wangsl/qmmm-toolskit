 
C     $Id: disgeo.i 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  disgeo.i  --  distance geometry bounds and parameters  ##
c     ##                                                         ##
c     #############################################################
c
c
c     bnd         distance geometry upper and lower bounds matrix
c     vdwrad      hard sphere radii for distance geometry atoms
c     vdwmax      maximum value of hard sphere sum for an atom pair
c     compact     index of local distance compaction on embedding
c     pathmax     maximum value of upper bound after smoothing
c     use_invert  flag to use enantiomer closest to input structure
c     use_anneal  flag to use simulated annealing refinement
c
c
      real*8 bnd,vdwrad,vdwmax
      real*8 compact,pathmax
      logical use_invert
      logical use_anneal
      common /disgeo/ bnd(maxgeo,maxgeo),vdwrad(maxatm),vdwmax,compact,
     &                pathmax,use_invert,use_anneal
