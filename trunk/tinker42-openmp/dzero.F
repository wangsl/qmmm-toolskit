 
C     $Id: dzero.F 2 2008-02-28 22:52:52Z wangsl $

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
