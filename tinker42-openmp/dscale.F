 
C     $Id: dscale.F 2 2008-02-28 22:52:52Z wangsl $


      SubRoutine DScale(N, S, A)
      Implicit None
      Integer N
      Real*8 S, A(*)
      Integer I

C$OMP Parallel Do
C$OMP& Default(Shared), Private(I)
      Do I = 1, N
         A(I) = S*A(I)
      End Do
C$OMP End Parallel Do

      Return
      End
