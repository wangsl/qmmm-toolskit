 
C     $Id: gfloat.F 2 2008-02-28 22:52:52Z wangsl $
*Deck GFloat
      function gfloat (iarg)
      real*8 gfloat
c
c     integer to working-precision conversion
c
      gfloat=dfloat(iarg)
      return
      end
