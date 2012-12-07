 
C     $Id$

      SubRoutine DZero(N, A)
      Implicit None
      Integer N
      Real*8 A(*)
      Integer I

C$OMP Parallel Do
C$OMP& Default(Shared), Private(I)
      Do I = 1, N
         A(I) = 0.0D0
      End Do
C$OMP End Parallel Do
      
      Return
      End
      

      SubRoutine AClear(N, A)
      Implicit None
      Integer N
      Real*8 A(*)
      Call DZero(N, A)
      Return
      End
