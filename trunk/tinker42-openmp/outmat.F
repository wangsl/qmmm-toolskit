 
C     $Id$
*Deck OutMat
      Subroutine OutMat(IOut,Key,X,M,N,MM,NN)
      Implicit Real*8(A-H,O-Z)
C
C     Print matrix X.  ISmall is >0 to print all elements, <=0 to
C     print elements greater than 10**(Key-6)
C
      Parameter (NCol=5)
      Dimension X(M,N), Y(NCol)
      Save Zero
      Data Zero/0.0d0/
 1000 Format(8(10X,I3))
 1010 Format(I4,8D13.6)
C
      If(Key.gt.0) then
        Small = Zero
      else
        Small = GFloat(10)**(Key-6)
        endIf
      Do 40 ILower = 1, NN, NCol
        IUpper = Min0(ILower+NCol-1,NN)
        Num = IUpper - ILower + 1
        Write(IOut,1000) (J,J=ILower,IUpper)
        Do 30 I = 1, MM
          Do 10 J = ILower, IUpper
   10       Y(J-ILower+1) = X(I,J)
          Do 20 J = 1, Num
            If(Abs(Y(J)).lt.Small) Y(J) = Zero
   20       Continue
          Write(IOut,1010) I, (Y(J),J=1,Num)
   30     Continue
   40   Continue
      Return
      End
