 
C     $Id: squarelower.F 2 2008-02-28 22:52:52Z wangsl $
      
      SubRoutine SquareToLowerTriangular(A, N, L)
      Implicit None
      Integer N
      Real*8 A(N, *), L(*)
      Integer I, J, II

      II = 0
      Do I = 1, N
         Do J = 1, I
            II = II + 1
            L(II) = A(I, J)
         End Do
      End Do
      
      Return
      End

      SubRoutine LowerTriangularToSquare(L, A, N)
      Implicit None
      Integer N
      Real*8 L(*), A(N, *)
      Integer I, J, II

      II = 0
      Do I = 1, N
         Do J = 1, I
            II = II + 1
            A(I, J) = L(II)
            A(J, I) = A(I, J)
         End Do
      End Do

      Return
      End
