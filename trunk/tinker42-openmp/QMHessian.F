 
C     $Id$

#include "memfort.h"
      
      Program QMHessian
      Call CalculateQMHessian
      End

      SubRoutine CalculateQMHessian
      Implicit None
      Include 'iounit.i'
      Include 'sizes.i'
      Include 'atoms.i'
      Include 'linend.i'
      Include 'units.i'
      Include 'qmmm.i'
      
      Integer I, J, NAtom
      Real*8 H(*), H2(*)
      Pointer (pH, H), (pH2, H2)
      Data pH/0/, pH2/0/
      Integer NQM, LenH, LenH2
      Real*8 W(*), Work(*)
      Pointer (pW, W), (pWork, Work)
      Data pW, pWork/0, 0/
      Integer LWork, Info
      Integer FreeUnit, IUnit
      External FreeUnit

 1000 Format(/, ' QM Hessian matrix eigenvalues (Hartree/Bohr^2)', /)
      
      Call Initial
      Call GetXYZ
      Call Mechanic
      
      Write(IOut, *) IOut
      Write(IOut, *) QGrp, YGrp, M0Grp, MGrp,
     $     QMTotChg, QMSpMult

      NQM = QGrp + YGrp

      LenH = 3*NQM*(3*NQM+1)/2
      Call FORTAllocReal8(pH, LenH)
      Call DZero(LenH, H)
      
      Call HessianQM(H)

      Call LTOutS(IOut, 'Hessian', 3*NQM, 3*NQM, H, 1)

      LenH2 = 3*NQM * 3*NQM
      Call FORTAllocReal8(pH2, LenH2) 

      Call LowerTriangularToSquare(H, H2, 3*NQM)

      IUnit = FreeUnit()
      Open(IUnit, File = 'QM-Hessian', Status = 'Unknown')
      Write(IUnit, '(I8)') 3*NQM
      Do I = 1, 3*NQM*3*NQM
         Write(IUnit, '(E20.12)') H2(I)
      End Do
      Close(IUnit)      
      
      Call FORTAllocReal8(pW, 3*NQM)
      LWork = 3*3*NQM
      Call FORTAllocReal8(pWork, LWork)

      Call DSyEv('V', 'L', 3*NQM, H2, 3*NQM, W, Work, LWork, Info)

      If(pWork .ne. 0) Call FORTFree(pWork)

      If(Info .ne. 0) Call Crash('DSyeV error')

      Call DScale(3*NQM, Bohr**2/Hartree, W)

      Write(IOut, 1000)
      Do I = 1, 3*NQM
         Write(IOut, '(I4,F12.6)') I, W(I)
      End Do

      If(pW .ne. 0) Call FORTFree(pW)
      If(pH .ne. 0) Call FORTFree(pH)

      Return
      End



