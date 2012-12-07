 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine getref  --  get structure from reference area  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "getref" copies structure information from the reference area
c     into the standard variables for the current system structure
c
c
      subroutine getref
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'files.i'
      include 'refer.i'
      include 'titles.i'
      integer i,j
c
c
c     retrieve the filename and title line for the structure
c
      filename = reffile
      leng = refleng
      title = reftitle
      ltitle = refltitle
c
c     retrieve the coordinates, type and connectivity of each atom
c
      n = nref
      do i = 1, n
         name(i) = refnam(i)
         x(i) = xref(i)
         y(i) = yref(i)
         z(i) = zref(i)
         type(i) = reftyp(i)
         n12(i) = n12ref(i)
         do j = 1, n12(i)
            i12(j,i) = i12ref(j,i)
         end do
      end do
      return
      end
