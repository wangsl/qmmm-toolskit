 
C     $Id$
*Deck GFloat
      function gfloat (iarg)
      real*8 gfloat
c
c     integer to working-precision conversion
c
      gfloat=dfloat(iarg)
      return
      end
