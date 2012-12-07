 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine upcase  --  convert string to all upper case  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "upcase" converts a text string to all upper case letters
c
c
      subroutine upcase (string)
      implicit none
      integer i,length
      integer code,ichar
      character*1 char
      character*(*) string
c
c
c     convert lower case to upper case letters one at a time
c
      length = len(string)
      do i = 1, length
         code = ichar(string(i:i))
         if (code.ge.97 .and. code.le.122)
     &      string(i:i) = char(code-32)
      end do
      return
      end
