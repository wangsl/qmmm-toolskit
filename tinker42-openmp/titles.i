 
C     $Id: titles.i 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  titles.i  --  title for the current molecular system  ##
c     ##                                                        ##
c     ############################################################
c
c
c     ltitle   length in characters of the nonblank title string
c     title    title used to describe the current structure
c
c
      integer ltitle
      character*120 title
      common /titles/ ltitle,title
