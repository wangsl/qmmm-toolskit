 
C     $Id: crystal.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program crystal  --  fractional coordinate manipulations  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "crystal" is a utility program which converts between
c     fractional and Cartesian coordinates, and can generate
c     full unit cells from asymmetric units
c
c
      program crystal
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'iounit.i'
      include 'files.i'
      include 'math.i'
      integer maxspace
      parameter (maxspace=24)
      integer i,ixyz,mode
      integer na,nb,nc
      integer next,freeunit
      logical exist,query
      character*1 answer
      character*10 sgroup(maxspace)
      character*120 coordfile
      character*120 record
      character*120 string
      data sgroup / 'P1        ', 'P1(-)     ', 'P21       ',
     &              'Cc        ', 'P21/a     ', 'P21/n     ',
     &              'P21/c     ', 'C2/c      ', 'P212121   ',
     &              'Pca21     ', 'Pna21     ', 'Pn21a     ',
     &              'Cmc21     ', 'Pccn      ', 'Pbcn      ',
     &              'Pbca      ', 'P41       ', 'I41/a     ',
     &              'P4(-)21c  ', 'P4(-)m2   ', 'R3c       ',
     &              'P6(3)/mmc ', 'Fm3(-)m   ', 'Im3(-)m   '/
c
c
c     get and read the Cartesian coordinates file
c
      call initial
      call getxyz
c
c     find out which unitcell manipulation is to be performed
c
      mode = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  mode
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' The TINKER Crystal Facility can Provide :',
     &           //,4x,'(1) Convert Fractional to Cartesian Coords',
     &           /,4x,'(2) Convert Cartesian to Fractional Coords',
     &           /,4x,'(3) Translate the Molecules into Unit Cell',
     &           /,4x,'(4) Make a Unit Cell from Asymmetric Unit',
     &           /,4x,'(5) Make a Big Block from Single Unit Cell')
         dowhile (mode.lt.1 .or. mode.gt.5)
            write (iout,30)
   30       format (/,' Enter the Number of the Desired Choice :  ',$)
            read (input,40)  mode
   40       format (i10)
         end do
      end if
c
c     get any cell dimensions found in the keyword list
c
      call unitcell
c
c     determine the space group if it will be needed later
c
      if (mode .eq. 4) then
         do i = 1, maxspace
            if (spacegrp .eq. sgroup(i))  goto 90
         end do
   50    continue
         write (iout,60)  (sgroup(i),i=1,maxspace)
   60    format (/,' Listing of the Available Space Groups :',
     &           //,'   (1) ',a10,'  (2) ',a10,'  (3) ',a10,
     &           /,'   (4) ',a10,'  (5) ',a10,'  (6) ',a10,
     &           /,'   (7) ',a10,'  (8) ',a10,'  (9) ',a10,
     &           /,'  (10) ',a10,' (11) ',a10,' (12) ',a10,
     &           /,'  (13) ',a10,' (14) ',a10,' (15) ',a10,
     &           /,'  (16) ',a10,' (17) ',a10,' (18) ',a10,
     &           /,'  (19) ',a10,' (20) ',a10,' (21) ',a10,
     &           /,'  (22) ',a10,' (23) ',a10,' (24) ',a10)
         write (iout,70)
   70    format (/,' Enter the Number of the Desired Choice :  ',$)
         read (input,80)  i
   80    format (i10)
         if (i.lt.1 .or. i.gt.maxspace)  goto 50
         spacegrp = sgroup(i)
   90    continue
      end if
c
c     if not in keyfile, get the unit cell axis lengths
c
      dowhile (xbox .eq. 0.0d0)
         write (iout,100)
  100    format (/,' Enter Unit Cell Axis Lengths :  ',$)
         read (input,110)  record
  110    format (a120)
         read (record,*,err=120,end=120)  xbox,ybox,zbox
  120    continue
         if (ybox .eq. 0.0d0)  ybox = xbox
         if (zbox .eq. 0.0d0)  zbox = xbox
         use_bounds = .true.
      end do
c
c     if not in keyfile, get the unit cell angle values
c
      dowhile (alpha .eq. 0.0d0)
         write (iout,130)
  130    format (/,' Enter Unit Cell Axis Angles :   ',$)
         read (input,140)  record
  140    format (a120)
         read (record,*,err=150,end=150)  alpha,beta,gamma
  150    continue
         if (alpha .eq. 0.0d0)  alpha = 90.0d0
         if (beta .eq. 0.0d0)  beta = alpha
         if (gamma .eq. 0.0d0)  gamma = alpha
         if (alpha.eq.90.0d0 .and. beta.eq.90.0d0
     &          .and. gamma.eq.90.0d0) then
            orthogonal = .true.
         else if (alpha.eq.90.0d0 .and. gamma.eq.90.0d0) then
            monoclinic = .true.
         else
            triclinic = .true.
         end if
      end do
c
c     find constants for coordinate interconversion
c
      call lattice
c
c     print out the initial cell dimensions to be used
c
      write (iout,160)  xbox,ybox,zbox,alpha,beta,gamma
  160 format (/,' Unit Cell Dimensions :      a    =',f10.4,
     &        /,'                             b    =',f10.4,
     &        /,'                             c    =',f10.4,
     &        /,'                            Alpha =',f10.4,
     &        /,'                            Beta  =',f10.4,
     &        /,'                            Gamma =',f10.4)
c
c     convert Cartesian to fractional coordinates
c
      if (mode.ne.1 .and. mode.ne.3) then
         do i = 1, n
            z(i) = (z(i)/gamma_term) / zbox
            y(i) = ((y(i)-z(i)*zbox*beta_term)/gamma_sin) / ybox
            x(i) = (x(i)-y(i)*ybox*gamma_cos-z(i)*zbox*beta_cos) / xbox
         end do
      end if
c
c     apply the appropriate space group symmetry operators
c
      if (mode .eq. 4) then
         write (iout,170)  spacegrp
  170    format (/,' Space Group Symbol :',12x,a10)
         call symmetry (spacegrp)
      end if
c
c     replicate the unit cell to make a block of unit cells
c
      if (mode .eq. 5) then
         na = 0
         nb = 0
         nc = 0
         write (iout,180)
  180    format (/,' Enter Number of Replicates along a-, b- and',
     &              ' c-Axes [1 1 1] :   ',$)
         read (input,190)  record
  190    format (a120)
         read (record,*,err=200,end=200)  na,nb,nc
  200    continue
         if (na .eq. 0)  na = 1
         if (nb .eq. 0)  nb = na
         if (nc .eq. 0)  nc = na
         if (na*nb*nc*n .gt. maxatm) then
            write (iout,210)  maxatm
  210       format (/,' CRYSTAL  --  The Maximum of',i6,' Atoms',
     &                 ' has been Exceeded')
            call fatal
         end if
         call bigblock (na,nb,nc)
         write (iout,220)  na,nb,nc,dble(na)*xbox,dble(nb)*ybox,
     &                     dble(nc)*zbox
  220    format (/,' Dimensions of the',i3,' x',i3,' x',i3,
     &              ' Cell Block :',
     &           //,' New Cell Dimensions :       a    =',f10.4,
     &            /,'                             b    =',f10.4,
     &            /,'                             c    =',f10.4)
      end if
c
c     convert fractional to Cartesian coordinates
c
      if (mode.ne.2 .and. mode.ne.3) then
         do i = 1, n
            x(i) = x(i)*xbox + y(i)*ybox*gamma_cos + z(i)*zbox*beta_cos
            y(i) = y(i)*ybox*gamma_sin + z(i)*zbox*beta_term
            z(i) = z(i)*zbox*gamma_term
         end do
      end if
c
c     translate any stray molecules back into the unit cell
c
      if (mode .eq. 3) then
         call field
         call katom
         call molecule
         call bounds
      end if
c
c     center unit cell at the origin and move molecules into cell
c
      if (mode .eq. 4) then
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,230)
  230       format (/,' Center on Origin and Move Molecules into',
     &                 ' Unit Cell [Y] :   ',$)
            read (input,240)  record
  240       format (a120)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .ne. 'N') then
            do i = 1, n
               z(i) = (z(i)/gamma_term) / zbox
               y(i) = ((y(i)-z(i)*zbox*beta_term)/gamma_sin) / ybox
               x(i) = (x(i)-y(i)*ybox*gamma_cos-z(i)*zbox*beta_cos)
     &                                  / xbox
            end do
            do i = 1, n
               x(i) = x(i) - 0.5d0
               y(i) = y(i) - 0.5d0
               z(i) = z(i) - 0.5d0
            end do
            do i = 1, n
               x(i) = x(i)*xbox + y(i)*ybox*gamma_cos
     &                   + z(i)*zbox*beta_cos
               y(i) = y(i)*ybox*gamma_sin + z(i)*zbox*beta_term
               z(i) = z(i)*zbox*gamma_term
            end do
            call field
            call katom
            call molecule
            call bounds
         end if
      end if
c
c     write out the new coordinates to a file
c
      ixyz = freeunit ()
      if (mode .eq. 2) then
         coordfile = filename(1:leng)//'.cel'
         call version (coordfile,'new')
         open (unit=ixyz,file=coordfile,status='new')
      else
         coordfile = filename(1:leng)//'.xyz'
         call version (coordfile,'new')
         open (unit=ixyz,file=coordfile,status='new')
      end if
      call prtxyz (ixyz)
      close (unit=ixyz)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine cellatom  --  add new atom to the unit cell  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "cellatom" completes the addition of a symmetry related atom
c     to a unit cell by updating the atom type and attachment arrays
c
c
      subroutine cellatom (jj,j)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      integer i,j,jj,delta
c
c
c     attachments of replicated atom are analogous to base atom
c
      delta = jj - j
      n12(jj) = n12(j)
      do i = 1, n12(j)
         i12(i,jj) = i12(i,j) + delta
      end do
      type(jj) = type(j)
      name(jj) = name(j)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine bigblock  --  create a block of unit cells  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "bigblock" replicates the coordinates of a single unit
c     cell to give a larger block of repeated units
c
c
      subroutine bigblock (na,nb,nc)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer i,j,k
      integer ii,jj,nsym
      integer na,nb,nc
      real*8 trans(3,maxcell)
c
c
c     construct translation offsets for the replicated cells
c
      nsym = 0
      do i = (1-na)/2, na/2
         do j = (1-nb)/2, nb/2
            do k = (1-nc)/2, nc/2
               nsym = nsym + 1
               trans(1,nsym) = i
               trans(2,nsym) = j
               trans(3,nsym) = k
            end do
         end do
      end do
c
c     put the original cell at the top of the replica list
c
      do i = 1, nsym
         if (trans(1,i).eq.0 .and. trans(2,i).eq.0
     &           .and. trans(3,i).eq.0)  k = i
      end do
      do i = k, 2, -1
         trans(1,i) = trans(1,i-1)
         trans(2,i) = trans(2,i-1)
         trans(3,i) = trans(3,i-1)
      end do
      trans(1,1) = 0
      trans(2,1) = 0
      trans(3,1) = 0
c
c     translate the original unit cell to make a block of cells
c
      do i = 2, nsym
         ii = (i-1) * n
         do j = 1, n
            jj = j + ii
            x(jj) = x(j) + trans(1,i)
            y(jj) = y(j) + trans(2,i)
            z(jj) = z(j) + trans(3,i)
            call cellatom (jj,j)
         end do
      end do
      n = nsym * n
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine symmetry  --  apply space group symmetry  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "symmetry" applies symmetry operators to the fractional
c     coordinates of the asymmetric unit in order to generate
c     the symmetry related atoms of the full unit cell
c
c
      subroutine symmetry (spacegrp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer i,j,ii,jj,nsym
      real*8 one3,two3
      character*10 spacegrp
c
c
c     P1 space group  (Intl. Tables 1)
c
      if (spacegrp .eq. 'P1        ') then
         nsym = 1
         do i = 1, n
            x(i) = x(i) - 0.5d0
            y(i) = y(i) - 0.5d0
            z(i) = z(i) - 0.5d0
         end do
c
c     P1(-) space group  (Intl. Tables 2)
c
      else if (spacegrp .eq. 'P1(-)     ') then
         nsym = 2
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P21 space group  (Intl. Tables 4)
c
      else if (spacegrp .eq. 'P21       ') then
         nsym = 2
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Cc space group  (Intl. Tables 9)
c
      else if (spacegrp .eq. 'Cc        ') then
         nsym = 4
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P21/a space group  (Intl. Tables 14)
c
      else if (spacegrp .eq. 'P21/a     ') then
         nsym = 4
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P21/n space group  (Intl. Tables 14)
c
      else if (spacegrp .eq. 'P21/n     ') then
         nsym = 4
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P21/c space group  (Intl. Tables 14)
c
      else if (spacegrp .eq. 'P21/c     ') then
         nsym = 4
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 4) then
                  x(jj) = x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     C2/c space group  (Intl. Tables 15)
c
      else if (spacegrp .eq. 'C2/c      ') then
         nsym = 8
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 5) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 7) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = -z(j)
               else if (i .eq. 8) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P212121 space group  (Intl. Tables 19)
c
      else if (spacegrp .eq. 'P212121   ') then
         nsym = 4
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pca21 space group  (Intl. Tables 29)
c
      else if (spacegrp .eq. 'Pca21     ') then
         nsym = 4
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pna21 space group  (Intl. Tables 33)
c
      else if (spacegrp .eq. 'Pna21     ') then
         nsym = 4
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pn21a space group  (Intl. Tables 33)
c
      else if (spacegrp .eq. 'Pn21a     ') then
         nsym = 4
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Cmc21 space group  (Intl. Tables 36)
c
      else if (spacegrp .eq. 'Cmc21     ') then
         nsym = 8
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 4) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 5) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 7) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 8) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pccn space group  (Intl. Tables 56)
c
      else if (spacegrp .eq. 'Pccn      ') then
         nsym = 8
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 4) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 5) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               else if (i .eq. 7) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 8) then
                  x(jj) = x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pbcn space group  (Intl. Tables 60)
c
      else if (spacegrp .eq. 'Pbcn      ') then
         nsym = 8
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 5) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 7) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = z(j)
               else if (i .eq. 8) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pbca space group  (Intl. Tables 61)
c
      else if (spacegrp .eq. 'Pbca      ') then
         nsym = 8
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = -z(j)
               else if (i .eq. 5) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 7) then
                  x(jj) = x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 8) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P41 space group  (Intl. Tables 76)
c
      else if (spacegrp .eq. 'P41       ') then
         nsym = 4
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = -y(j)
                  y(jj) = x(j)
                  z(jj) = 0.25d0 + z(j)
               else if (i .eq. 4) then
                  x(jj) = y(j)
                  y(jj) = -x(j)
                  z(jj) = 0.75d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     I41/a space group  (Intl. Tables 88)
c
      else if (spacegrp .eq. 'I41/a     ') then
         nsym = 16
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 4) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = z(j)
               else if (i .eq. 5) then
                  x(jj) = 0.75d0 - y(j)
                  y(jj) = 0.25d0 + x(j)
                  z(jj) = 0.25d0 + z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.25d0 - y(j)
                  y(jj) = 0.75d0 + x(j)
                  z(jj) = 0.75d0 + z(j)
               else if (i .eq. 7) then
                  x(jj) = 0.75d0 + y(j)
                  y(jj) = 0.75d0 - x(j)
                  z(jj) = 0.75d0 + z(j)
               else if (i .eq. 8) then
                  x(jj) = 0.25d0 + y(j)
                  y(jj) = 0.25d0 - x(j)
                  z(jj) = 0.25d0 + z(j)
               else if (i .eq. 9) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 10) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 11) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 12) then
                  x(jj) = x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               else if (i .eq. 13) then
                  x(jj) = 0.25d0 + y(j)
                  y(jj) = 0.75d0 - x(j)
                  z(jj) = 0.75d0 - z(j)
               else if (i .eq. 14) then
                  x(jj) = 0.75d0 + y(j)
                  y(jj) = 0.25d0 - x(j)
                  z(jj) = 0.25d0 - z(j)
               else if (i .eq. 15) then
                  x(jj) = 0.25d0 - y(j)
                  y(jj) = 0.25d0 + x(j)
                  z(jj) = 0.25d0 - z(j)
               else if (i .eq. 16) then
                  x(jj) = 0.75d0 - y(j)
                  y(jj) = 0.75d0 + x(j)
                  z(jj) = 0.75d0 - z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P4(-)21c space group  (Intl. Tables 114)
c
      else if (spacegrp .eq. 'P4(-)21c  ') then
         nsym = 8
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = y(j)
                  y(jj) = -x(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = -y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               else if (i .eq. 5) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 7) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 8) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P4(-)m2 space group  (Intl. Tables 115)
c
      else if (spacegrp .eq. 'P4(-)m2   ') then
         nsym = 8
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 4) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 5) then
                  x(jj) = -y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = y(j)
                  y(jj) = -x(j)
                  z(jj) = -z(j)
               else if (i .eq. 7) then
                  x(jj) = y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               else if (i .eq. 8) then
                  x(jj) = -y(j)
                  y(jj) = -x(j)
                  z(jj) = -z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     R3c space group  (Intl. Tables 161)
c
      else if (spacegrp .eq. 'R3c       ') then
         nsym = 18
         one3 = 1.0d0 / 3.0d0
         two3 = 2.0d0 / 3.0d0
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = -y(j)
                  y(jj) = x(j) - y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = y(j) - x(j)
                  y(jj) = -x(j)
                  z(jj) = z(j)
               else if (i .eq. 4) then
                  x(jj) = -y(j)
                  y(jj) = -x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 5) then
                  x(jj) = x(j)
                  y(jj) = x(j) - y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 6) then
                  x(jj) = y(j) - x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 7) then
                  x(jj) = x(j) + one3
                  y(jj) = y(j) + two3
                  z(jj) = z(j) + two3
               else if (i .eq. 8) then
                  x(jj) = -y(j) + one3
                  y(jj) = x(j) - y(j) + two3
                  z(jj) = z(j) + two3
               else if (i .eq. 9) then
                  x(jj) = y(j) - x(j) + one3
                  y(jj) = -x(j) + two3
                  z(jj) = z(j) + two3
               else if (i .eq. 10) then
                  x(jj) = -y(j) + one3
                  y(jj) = -x(j) + two3
                  z(jj) = 0.5d0 + z(j) + two3
               else if (i .eq. 11) then
                  x(jj) = x(j) + one3
                  y(jj) = x(j) - y(j) + two3
                  z(jj) = 0.5d0 + z(j) + two3
               else if (i .eq. 12) then
                  x(jj) = y(j) - x(j) + one3
                  y(jj) = y(j) + two3
                  z(jj) = 0.5d0 + z(j) + two3
               else if (i .eq. 13) then
                  x(jj) = x(j) + two3
                  y(jj) = y(j) + one3
                  z(jj) = z(j) + one3
               else if (i .eq. 14) then
                  x(jj) = -y(j) + two3
                  y(jj) = x(j) - y(j) + one3
                  z(jj) = z(j) + one3
               else if (i .eq. 15) then
                  x(jj) = y(j) - x(j) + two3
                  y(jj) = -x(j) + one3
                  z(jj) = z(j) + one3
               else if (i .eq. 16) then
                  x(jj) = -y(j) + two3
                  y(jj) = -x(j) + one3
                  z(jj) = 0.5d0 + z(j) + one3
               else if (i .eq. 17) then
                  x(jj) = x(j) + two3
                  y(jj) = x(j) - y(j) + one3
                  z(jj) = 0.5d0 + z(j) + one3
               else if (i .eq. 18) then
                  x(jj) = y(j) - x(j) + two3
                  y(jj) = y(j) + one3
                  z(jj) = 0.5d0 + z(j) + one3
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P6(3)/mmc space group  (Intl. Tables 194, Hexagonal Close Packed)
c
      else if (spacegrp .eq. 'P6(3)/mmc ') then
         nsym = 2
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = x(j) + 0.5d0
                  y(jj) = y(j)
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Fm3m space group  (Intl. Tables 225, Face Centered Cubic)
c
      else if (spacegrp .eq. 'Fm3(-)m   ') then
         nsym = 4
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = x(j)
                  y(jj) = y(j) + 0.5d0
                  z(jj) = z(j) + 0.5d0
               else if (i .eq. 3) then
                  x(jj) = x(j) + 0.5d0
                  y(jj) = y(j)
                  z(jj) = z(j) + 0.5d0
               else if (i .eq. 4) then
                  x(jj) = x(j) + 0.5d0
                  y(jj) = y(j) + 0.5d0
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Im3m space group  (Intl. Tables 229, Body Centered Cubic)
c
      else if (spacegrp .eq. 'Im3(-)m   ') then
         nsym = 2
         do i = 2, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 2) then
                  x(jj) = x(j) + 0.5d0
                  y(jj) = y(j) + 0.5d0
                  z(jj) = z(j) + 0.5d0
               end if
               call cellatom (jj,j)
            end do
         end do
      end if
c
c     set the total number of atoms to include full unitcell
c
      n = nsym * n
      return
      end
