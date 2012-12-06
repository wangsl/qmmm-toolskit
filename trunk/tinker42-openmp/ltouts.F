 
C     $Id: ltouts.F 2 2008-02-28 22:52:52Z wangsl $
*Deck LTOutS
      Subroutine LTOutS(IOut,String,Num,N,A,Key)
      Implicit Real*8(A-H,O-Z)
C
C     Print a heading, followed by a lower triangular matrix.
C     If Num is positive, it is also printed.
C
      Character*(*) String
      Dimension A(*)
 1000 Format(1X,A,I10)
 1010 Format(1X,A)
C
      LStr = Max0(LinEnd(String),1)
      If(Num.gt.0) then
        Write(IOut,1000) String(1:LStr), Num
      else
        Write(IOut,1010) String(1:LStr)
        endIf
      Call LTOut(IOut,N,A,Key)
      Return
      End
