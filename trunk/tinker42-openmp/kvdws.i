 
C     $Id: kvdws.i 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  kvdws.i  --  forcefield parameters for van der Waals terms  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     rad      van der Waals radius parameter for each atom class
c     eps      van der Waals well depth parameter for each atom class
c     rad4     van der Waals radius parameter in 1-4 interactions
c     eps4     van der Waals well depth parameter in 1-4 interactions
c     reduct   van der Waals reduction factor for each atom class
c
c
      real*8 rad,eps,rad4,eps4,reduct
      common /kvdws/ rad(maxclass),eps(maxclass),rad4(maxclass),
     &               eps4(maxclass),reduct(maxclass)
