 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  angpot.i  --  specifics of bond angle functional forms  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     cang       cubic coefficient in angle bending potential
c     qang       quartic coefficient in angle bending potential
c     pang       quintic coefficient in angle bending potential
c     sang       sextic coefficient in angle bending potential
c     angunit    convert angle bending energy to kcal/mole
c     stbnunit   convert stretch-bend energy to kcal/mole
c     aaunit     convert angle-angle energy to kcal/mole
c     opbunit    convert out-of-plane bend energy to kcal/mole
c     opdunit    convert out-of-plane distance energy to kcal/mole
c     mm2stbn    logical flag governing use of MM2-style stretch-bend
c
c
      real*8 cang,qang,pang,sang
      real*8 angunit,stbnunit,aaunit
      real*8 opbunit,opdunit
      logical mm2stbn
      common /angpot/ cang,qang,pang,sang,angunit,stbnunit,aaunit,
     &                opbunit,opdunit,mm2stbn
