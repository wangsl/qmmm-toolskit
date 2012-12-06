 
C     $Id: trimtext.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  function trimtext  --  find last non-blank character  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "trimtext" finds and returns the location of the last
c     non-blank character before the first null character in
c     an input text string; the function returns zero if no
c     such character is found
c
c
      function trimtext (string)
      implicit none
      integer i,size,len
      integer last,trimtext
      character*1 null
      character*(*) string
c
c
c     move forward through the string, one character
c     at a time, looking for first null character
c
      trimtext = 0
      null = char(0)
      size = len(string)
      last = size
      do i = 1, size
         if (string(i:i) .eq. null) then
            last = i - 1
            goto 10
         end if
      end do
   10 continue
c
c     move backward through the string, one character
c     at a time, looking for first non-blank character
c
      do i = last, 1, -1
         if (string(i:i) .gt. ' ') then
            trimtext = i
            goto 20
         end if
      end do
   20 continue
      return
      end
