 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine hessian  --  atom-by-atom Hessian elements  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "hessian" calls subroutines to calculate the Hessian elements
c     for each atom in turn with respect to Cartesian coordinates
c

C     Modified by Shenglong Wang to calculate Hessian for QM atoms
C     from MM force field, which will be used as the initial guess
C     for QMMM iterative optimization
C     according to Yingkai's modification for Gaussian + Tinker

c
C     subroutine hessian (h,hinit,hstop,hindex,hdiag)

C     QMHess: lower triangular matrix

#include "memfort.h"

      SubRoutine HessianQM(QMHess)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'couple.i'
      include 'hescut.i'
      include 'hessn.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'potent.i'
      include 'rigid.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      Include 'linend.i'
      Include 'qmmm.i'
      integer i,j,k
      integer ii
      real*8 cutoff,rdn
      
      Real*8 XRed(*), YRed(*), ZRed(*)
      Pointer (pXRed, XRed), (pYRed, YRed), (pZRed, ZRed)
      Data pXRed, pYRed, pZRed/0, 0, 0/

      Real*8 Hx(1), Hy(1), Hz(1)
      Equivalence(Hx, HessX)
      Equivalence(Hy, HessY)
      Equivalence(Hz, HessZ)

      Real*8 QMHess(*)
      Integer Index, NQM, IHess

 1000 Format(' QM atom ', I4, ' is inactive')

      Index(I, J) = I*(I-1)/2 + J

      NQM = QGrp + YGrp

      Call AClear(Index(3*NQM, 3*NQM), QMHess)
      
c     maintain any periodic boundary conditions
c
      If(Use_bounds)
     $     Call QCrash('Sphere boundary only')

      if (use_bounds .and. .not.use_rigid)  call bounds
c     
c     alter bond and torsion constants for pisystem
c
      if (use_orbit)  call piscf
c
c     set the Born radii for use with GB/SA solvation
c
      if (use_gbsa)  call born
c
c     compute the induced dipoles at polarizable atoms
c
      if (use_polar) then
         call chkpole
         call rotpole
         call induce
      end if
c
c     calculate the "reduced" atomic coordinates
c     
      if (use_vdw) then
         Call FORTQAllocReal8(pXRed, N)
         Call FORTQAllocReal8(pYRed, N)
         Call FORTQAllocReal8(pZRed, N)
         do i = 1, n
            ii = ired(i)
            rdn = kred(i)
            xred(i) = rdn*(x(i)-x(ii)) + x(ii)
            yred(i) = rdn*(y(i)-y(ii)) + y(ii)
            zred(i) = rdn*(z(i)-z(ii)) + z(ii)
         end do
      end if
c
c     zero out the Hessian elements for the current atom

      do i = 1, NQM
         
         If(.not. Use(I)) then
            Write(IOut, 1000) I
            Call QCrash('QM atom is inactive')
         End If
         
         do k = 1, n
            do j = 1, 3
               hessx(j,k) = 0.0d0
               hessy(j,k) = 0.0d0
               hessz(j,k) = 0.0d0
            end do
         end do
c     
c     remove any previous use of the replicates method
c     
         cutoff = 0.0d0
         if (use_image)  call replica (cutoff)
c     
c     call the local geometry Hessian component routines
c     
         if (use_bond)  call ebond2 (i)
         if (use_angle)  call eangle2 (i)
         if (use_strbnd)  call estrbnd2 (i)
         if (use_urey)  call eurey2 (i)
         if (use_angang)  call eangang2 (i)
         if (use_opbend)  call eopbend2 (i)
         if (use_opdist)  call eopdist2 (i)
         if (use_improp)  call eimprop2 (i)
         if (use_imptor)  call eimptor2 (i)
         if (use_tors)  call etors2 (i)
         if (use_pitors)  call epitors2 (i)
         if (use_strtor)  call estrtor2 (i)
         if (use_tortor)  call etortor2 (i)
c
c     call the van der Waals Hessian component routines
c     
         if (use_vdw) then
            if (vdwtyp .eq. 'LENNARD-JONES') then
               call elj2 (i,xred,yred,zred)
            else if (vdwtyp .eq. 'BUCKINGHAM') then
               call ebuck2 (i,xred,yred,zred)
            else if (vdwtyp .eq. 'MM3-HBOND') then
               call emm3hb2 (i,xred,yred,zred)
            else if (vdwtyp .eq. 'BUFFERED-14-7') then
               call ehal2 (i,xred,yred,zred)
            else if (vdwtyp .eq. 'GAUSSIAN') then
               call egauss2 (i,xred,yred,zred)
            end if
         end if
c     
c     call the electrostatic Hessian component routines
c     
         if (use_charge) call echarge2 (i)
         if (use_chgdpl)  call echgdpl2 (i)
         if (use_dipole)  call edipole2 (i)
         if (use_mpole .or. use_polar)   call empole2 (i)
         if (use_rxnfld)   call erxnfld2 (i)
c     
c     call any miscellaneous Hessian component routines
c     
         if (use_solv)  call esolv2 (i)
         if (use_metal)  call emetal2 (i)
         if (use_geom)  call egeom2 (i)
         if (use_extra)  call extra2 (i)

C     Now copy HessX, HessY and HessZ to QMHess
         
         IHess = 3*(I-1)

C     Hess-X
         Do J = 1, IHess+1
            QMHess(Index(IHess+1, J)) = Hx(J)
         End Do
         
C     Hess-Y
         Do J = 1, IHess+2
            QMHess(Index(IHess+2, J)) = Hy(J)
         End Do

C     Hess-Z
         Do J = 1, IHess+3
            QMHess(Index(IHess+3, J)) = Hz(J)
         End Do

      End Do

      If(pXRed .ne. 0) Call FORTQFree(pXRed)
      If(pYRed .ne. 0) Call FORTQFree(pYRed)
      If(pZRed .ne. 0) Call FORTQFree(pZRed)

      Return
      End


