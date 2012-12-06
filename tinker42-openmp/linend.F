 
C     $Id: linend.F 2 2008-02-28 22:52:52Z wangsl $

      integer function linend(cline)
      implicit real*8(a-h,o-z)
c
c     function which returns the length of a character string,
c     excluding trailing blanks, tabs and nulls.
c
      character*(*) cline, blank
      character*1 null, tab
      parameter (blank=' ')
c
      null = char(0)
      tab = char(9)
      do 5 i = len(cline), 1, -1
          if(cline(i:i).ne.blank.and.
     $    cline(i:i).ne.null.and.
     $    cline(i:i).ne.tab) goto 10
    5     continue
      linend = 0
      return
   10 linend = i
      return
      end
