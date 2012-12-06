 
C     $Id: tinker-test.F 5 2008-03-08 02:30:43Z wangsl $

#include "memfort.h"

      Program TinkerTestMain
      Call TinkerTest
      End
      
      SubRoutine TinkerTest
      Implicit None
      Include 'iounit.i'
      Include 'sizes.i'
      Include 'atoms.i'
      Include 'linend.i'
      Include 'molcul.i'
      Include 'usage.i'
      Include 'atmtyp.i'
      Include 'couple.i'
      Include 'inter.i'
      Include 'virial.i'
      Integer I, J, IAtom, NH
      Logical FixedMol
      Real*8 E, DE(3,*)
      Pointer (pDE, DE)
      Data pDE/0/
      
      call initial
      call getxyz
      call mechanic
      
      E = 0.0D0
      Call FORTAllocReal8(pDE, 3*N)

      Call Gradient(E, DE)
      
      Call WriteEGRadients(E, EInter, DE, Vir)

      Write(IOut, *) 'E = ', E
      Write(IOut, *) 'EInter = ', EInter
      Write(IOut, *) 'Vir ', Vir

      Do I = 1, 50
         Write(IOut, '(I6, 3F20.10)') I, (DE(J,I), J = 1, 3)
      End Do

      If(pDE .ne. 0) Call FORTFree(pDE)
      
      Return
      End


