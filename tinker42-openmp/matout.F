 
C     $Id$
*Deck MatOut
      Subroutine MatOut(X,M,N,MM,NN)
      Implicit Real*8(A-H,O-Z)
      Include 'iounit.i'
C
C     Matrix print routine.
C
C     Common/IO/In,IOut,IPunch
      Dimension X(M,N)
C
      Call OutMat(IOut,1,X,M,N,MM,NN)
      Return
      End
