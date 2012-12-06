
C     $Id: fixwater.F 3 2008-03-04 17:47:30Z wangsl $

      SubRoutine FixWater
      Implicit None
      Include 'iounit.i'
      Include 'sizes.i'
      Include 'atoms.i'
      Include 'usage.i'
      Include 'atmtyp.i'
      Include 'couple.i'
      Integer I, J, IAtom, NH
      Logical FixedMol
      
 1000 Format(' There are ', I6, ' active atoms after ',
     $     'fixed H2O molecules')
      
      Do I = 1, N
         If(Atomic(I) .ne. 8) Cycle
         
         If(Atomic(I) .eq. 8) then
            NH = 0
            Do J = 1, N12(I)
               If(Atomic(I12(J,I)) .eq. 1) NH = NH + 1
            End Do
         End If
         
         If(NH .ne. 2) Cycle
         
         FixedMol = .False.
         If(.not. Use(I)) FixedMol = .True.
         Do J = 1, N12(I)
            If(.not. Use(I12(J,I))) FixedMol = .True.
         End Do
         
         If(.not. FixedMol) Cycle
         
         If(Use(I)) then
            Use(I) = .False.
            NUse = NUse - 1
         End If
         
         Do J = 1, N12(I)
            IAtom = I12(J,I)
            If(Use(IAtom)) then
               Use(IAtom) = .False.
               NUse = NUse - 1
            End If
         End Do
      End Do
      
      Write(IOut, 1000) NUse

      Return
      End
