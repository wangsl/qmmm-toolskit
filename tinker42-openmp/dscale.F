 
C     $Id$


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
