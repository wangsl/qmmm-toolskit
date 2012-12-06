 
C     $Id: wrtegrads.F 2 2008-02-28 22:52:52Z wangsl $

      SubRoutine WriteEGradients(E, Einter, Grads, Vir)
      Implicit None
      Include 'sizes.i'
      Include 'atoms.i'
      Include 'usage.i'
      Real*8 E, EInter
      Real*8 Grads(3,*), Vir(3,3)
      Integer FreeUnit
      External FreeUnit
      Integer IUnit
      Integer I, J
      
      IUnit = FreeUnit()
      
      Open(IUnit, File = 'EGrads-OpenMP', Status = 'Unknown')
      
      Write(IUnit, '(2F20.8)') E, EInter
      Do I = 1, N
         If(Use(I))
     $        Write(IUnit, '(I8, 3F20.8)') I, (Grads(J,I), J = 1, 3)
      End Do
      Do I = 1, 3
         Write(IUnit, '(I8, 3F20.8)') I, (Vir(J,I), J = 1, 3)
      End Do
      
      Close(IUnit)

      Return
      End
