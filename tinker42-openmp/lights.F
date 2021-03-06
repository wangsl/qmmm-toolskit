 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine lights  --  get neighbors via method of lights  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "lights" computes the set of nearest neighbor interactions
c     using the method of lights algorithm
c
c     literature reference:
c
c     F. Sullivan, R. D. Mountain and J. O'Connell, "Molecular
c     Dynamics on Vector Computers", Journal of Computational
c     Physics, 61, 138-153 (1985)
c
c
      subroutine lights (nsite,map,xsort,ysort,zsort)
      implicit none
      include 'sizes.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'iounit.i'
      include 'light.i'
      include 'shunt.i'
      integer i,j,k,nsite
      integer map(maxlight)
      real*8 box,xcut,ycut,zcut
      real*8 xmove,ymove,zmove
      real*8 xsort(maxlight)
      real*8 ysort(maxlight)
      real*8 zsort(maxlight)
c
c
c     check that maximum number of replicates is not exceeded
c
      if (use_replica) then
         if (xcell2.gt.xbox .or. ycell2.gt.ybox
     &           .or. zcell2.gt.zbox) then
            write (iout,10)
   10       format (/,' LIGHTS  --  Cutoff Distance is Too Large',
     &                 ' for Method of Lights')
            call fatal
         end if
      end if
c
c     truncated octahedron periodicity is not handled at present
c
      if (use_bounds) then
         if (octahedron) then
            write (iout,20)
   20       format (/,' LIGHTS  --  Method of Lights not available',
     &                 ' for Truncated Octahedron')
            call fatal
         end if
      end if
c
c     when using images, move coordinates into periodic cell
c
      if (use_image) then
         do i = 1, nsite
            zsort(i) = zsort(i) / gamma_term
            ysort(i) = (ysort(i) - zsort(i)*beta_term) / gamma_sin
            xsort(i) = xsort(i) - ysort(i)*gamma_cos - zsort(i)*beta_cos
            dowhile (xsort(i) .gt. xcell2)
               xsort(i) = xsort(i) - xcell
            end do
            dowhile (xsort(i) .lt. -xcell2)
               xsort(i) = xsort(i) + xcell
            end do
            dowhile (ysort(i) .gt. ycell2)
               ysort(i) = ysort(i) - ycell
            end do
            dowhile (ysort(i) .lt. -ycell2)
               ysort(i) = ysort(i) + ycell
            end do
            dowhile (zsort(i) .gt. zcell2)
               zsort(i) = zsort(i) - zcell
            end do
            dowhile (zsort(i) .lt. -zcell2)
               zsort(i) = zsort(i) + zcell
            end do
         end do
      end if
c
c     generate replica coordinates for the sort arrays
c
      if (use_replica) then
         k = nsite
         do j = 1, ncell
            xmove = icell(1,j) * xbox
            ymove = icell(2,j) * ybox
            zmove = icell(3,j) * zbox
            do i = 1, nsite
               k = k + 1
               map(k) = i
               xsort(k) = xsort(i) + xmove
               ysort(k) = ysort(i) + ymove
               zsort(k) = zsort(i) + zmove
               dowhile (xsort(k) .gt. xcell2)
                  xsort(k) = xsort(k) - xcell
               end do
               dowhile (xsort(k) .lt. -xcell2)
                  xsort(k) = xsort(k) + xcell
               end do
               dowhile (ysort(k) .gt. ycell2)
                  ysort(k) = ysort(k) - ycell
               end do
               dowhile (ysort(k) .lt. -ycell2)
                  ysort(k) = ysort(k) + ycell
               end do
               dowhile (zsort(k) .gt. zcell2)
                  zsort(k) = zsort(k) - zcell
               end do
               dowhile (zsort(k) .lt. -zcell2)
                  zsort(k) = zsort(k) + zcell
               end do
            end do
         end do
      end if
c
c     map image and replicate sites onto their parent sites
c
      nlight = (ncell+1) * nsite
      do i = 0, ncell
         j = i * nsite
         do k = 1, nsite
            map(j+k) = k
         end do
      end do
c
c     sort the coordinate components into ascending order
c
      call sort2 (nlight,xsort,locx)
      call sort2 (nlight,ysort,locy)
      call sort2 (nlight,zsort,locz)
c
c     use of replicates requires secondary sorting along x-axis
c
      if (use_replica) then
         j = 1
         do i = 1, nlight-1
            if (xsort(i+1) .ne. xsort(i)) then
               call sort5 (i-j+1,locx(j),nsite)
               j = i + 1
            end if
         end do
         call sort5 (nlight-j+1,locx(j),nsite)
      end if
c
c     index the position of each atom in the sorted coordinates
c
      do i = 1, nlight
         rgx(locx(i)) = i
         rgy(locy(i)) = i
         rgz(locz(i)) = i
      end do
c
c     set the light width based on interaction cutoff
c
      xcut = off
      ycut = off
      zcut = off
      if (use_image) then
         if (monoclinic) then
            zcut = zcut / beta_sin
            xcut = xcut + zcut*abs(beta_cos)
         else if (triclinic) then
            zcut = zcut / gamma_term
            ycut = (ycut + zcut*abs(beta_term)) / gamma_sin
            xcut = xcut + ycut*abs(gamma_cos) + zcut*abs(beta_cos)
         end if
         xcut = min(xcut,xcell2)
         ycut = min(ycut,ycell2)
         zcut = min(zcut,zcell2)
      end if
c
c     find the negative x-coordinate boundary for each atom
c
      do i = nlight, 1, -1
         k = locx(i)
         if (k .le. nsite) then
            kbx(k) = i
         end if
      end do
c
c     find the positive x-coordinate boundary for each atom
c
      j = 1
      box = 0.0d0
      do i = 1, nlight
         k = locx(i)
         if (k .le. nsite) then
            dowhile (xsort(j)-xsort(i)+box .lt. xcut)
               if (j .eq. nlight) then
                  if (use_image) then
                     j = 0
                     box = xcell
                  end if
               end if
               j = j + 1
               if (j .gt. nlight)  goto 30
            end do
   30       continue
            j = j - 1
            if (j .lt. 1) then
               j = nlight
               box = 0.0d0
            end if
            kex(k) = j
         end if
      end do
c
c     find the negative y-coordinate boundary for each atom
c
      j = nlight
      box = 0.0d0
      do i = nlight, 1, -1
         k = locy(i)
         if (k .le. nsite) then
            dowhile (ysort(i)-ysort(j)+box .le. ycut)
               if (j .eq. 1) then
                  if (use_image) then
                     j = nlight + 1
                     box = ycell
                  end if
               end if
               j = j - 1
               if (j .lt. 1)  goto 40
            end do
   40       continue
            j = j + 1
            if (j .gt. nlight) then
               j = 1
               box = 0.0d0
            end if
            kby(k) = j
         end if
      end do
c
c     find the positive y-coordinate boundary for each atom
c
      j = 1
      box = 0.0d0
      do i = 1, nlight
         k = locy(i)
         if (k .le. nsite) then
            dowhile (ysort(j)-ysort(i)+box .lt. ycut)
               if (j .eq. nlight) then
                  if (use_image) then
                     j = 0
                     box = ycell
                  end if
               end if
               j = j + 1
               if (j .gt. nlight)  goto 50
            end do
   50       continue
            j = j - 1
            if (j .lt. 1) then
               j = nlight
               box = 0.0d0
            end if
            key(k) = j
         end if
      end do
c
c     find the negative z-coordinate boundary for each atom
c
      j = nlight
      box = 0.0d0
      do i = nlight, 1, -1
         k = locz(i)
         if (k .le. nsite) then
            dowhile (zsort(i)-zsort(j)+box .le. zcut)
               if (j .eq. 1) then
                  if (use_image) then
                     j = nlight + 1
                     box = zcell
                  end if
               end if
               j = j - 1
               if (j .lt. 1)  goto 60
            end do
   60       continue
            j = j + 1
            if (j .gt. nlight) then
               j = 1
               box = 0.0d0
            end if
            kbz(k) = j
         end if
      end do
c
c     find the positive z-coordinate boundary for each atom
c
      j = 1
      box = 0.0d0
      do i = 1, nlight
         k = locz(i)
         if (k .le. nsite) then
            dowhile (zsort(j)-zsort(i)+box .lt. zcut)
               if (j .eq. nlight) then
                  if (use_image) then
                     j = 0
                     box = zcell
                  end if
               end if
               j = j + 1
               if (j .gt. nlight)  goto 70
            end do
   70       continue
            j = j - 1
            if (j .lt. 1) then
               j = nlight
               box = 0.0d0
            end if
            kez(k) = j
         end if
      end do
      return
      end
