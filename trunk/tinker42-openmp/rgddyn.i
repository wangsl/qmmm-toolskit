 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  rgddyn.i  --  velocities and momenta for rigid body MD  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     vcm     current translational velocity of each rigid body
c     wcm     current angular velocity of each rigid body
c     lm      current angular momentum of each rigid body
c     linear  logical flag to mark group as linear or nonlinear
c
c
      real*8 vcm,wcm,lm
      logical linear
      common /rgddyn/ vcm(3,maxgrp),wcm(3,maxgrp),lm(3,maxgrp),
     &                linear(maxgrp)
