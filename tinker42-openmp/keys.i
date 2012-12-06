 
C     $Id: keys.i 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  keys.i  --  contents of current keyword parameter file  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     nkey      number of nonblank lines in the keyword file
c     keyline   contents of each individual keyword file line
c
c
      integer nkey
      character*120 keyline
      common /keys/ nkey,keyline(maxkey)
