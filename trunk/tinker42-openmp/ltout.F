 
C     $Id$
*Deck LTOut
      Subroutine LTOut(IOut,N,A,Key)
      Implicit Real*8(A-H,O-Z)
C
C     Working precision routine to print out the lower triangular part
C     of a symmetric matrix stored in compressed lower triangular form.
C
C        N         Dimension of matrix.
C        A         Array to be printed.
C        KEY    N<=0 ... Suppress elements with absolute values less
C                        than 10**(-6+N)
C                  1 ... Print complete matrix.
C
      Parameter (NumCol=5)
      Dimension A(*), S(NumCol)
 2001 Format(5(11X,I3))
 2003 Format(I4,9D14.6)
      LInd(I,J) = (Max0(I,J)*(Max0(I,J)-1)/2) + Min0(I,J)
C
      Zero = GFloat(0)
      If(Key.gt.0) then
        Thresh = Zero
      else
        Thresh = GFloat(10)**(Key-6)
        endIf
      Do 20 IStart = 1, N, NumCol
        IEnd = Min0(IStart+NumCol-1,N)
        Write(IOut,2001) (IR,IR=IStart,IEnd)
        Do 20 IRow = IStart, N
          IRange = Min0(IRow-IStart+1,NumCol)
          L = LInd(IRow,IStart)
          Do 10 I = 1, IRange
            S(I) = A(L+I-1)
            If(Abs(S(I)).lt.Thresh) S(I) = Zero
   10       Continue
          Write(IOut,2003) IRow, (S(I),I=1,IRange)
   20     Continue
      Return
      End
