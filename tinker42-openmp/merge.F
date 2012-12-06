 
C     $Id: merge.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine merge  --  merge reference & current systems  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "merge" combines the reference and current structures into
c     a single new "current" structure containing the reference
c     atoms followed by the atoms of the current structure
c
c
      subroutine merge
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'iounit.i'
      include 'refer.i'
      integer i,j,k
      integer ntotal
c
c
c     check for too many total atoms in the combined system
c
      ntotal = n + nref
      if (ntotal .gt. maxatm) then
         write (iout,10)  maxatm
   10    format (/,' MERGE  --  The Maximum of',i6,' Atoms',
     &              ' has been Exceeded')
         call fatal
      end if
c
c     move the current structure to higher atom numbers
c
      do i = n, 1, -1
         k = i + nref
         x(k) = x(i)
         y(k) = y(i)
         z(k) = z(i)
         type(k) = type(i)
         name(k) = name(i)
         n12(k) = n12(i)
         do j = 1, n12(i)
            i12(j,k) = i12(j,i) + nref
         end do
      end do
c
c     place reference structure in the current structure
c
      call getref
      n = ntotal
      return
      end
