 
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
c     ##  subroutine makeref  --  copy structure to reference area  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "makeref" copies the information contained in the "xyz" file
c     of the current structure into corresponding reference areas
c
c
      subroutine makeref
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
c     copy the filename and title line for the structure
c
      reffile = filename
      refleng = leng
      reftitle = title
      refltitle = ltitle
c
c     copy the coordinates, type and connectivity of each atom
c
      nref = n
      do i = 1, n
         refnam(i) = name(i)
         xref(i) = x(i)
         yref(i) = y(i)
         zref(i) = z(i)
         reftyp(i) = type(i)
         n12ref(i) = n12(i)
         do j = 1, n12(i)
            i12ref(j,i) = i12(j,i)
         end do
      end do
      return
      end
