 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine lowcase  --  convert string to all lower case  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "lowcase" converts a text string to all lower case letters
c
c
      subroutine lowcase (string)
      implicit none
      integer i,length
      integer code,ichar
      character*1 char
      character*(*) string
c
c
c     convert upper case to lower case letters one at a time
c
      length = len(string)
      do i = 1, length
         code = ichar(string(i:i))
         if (code.ge.65 .and. code.le.90)
     &      string(i:i) = char(code+32)
      end do
      return
      end
