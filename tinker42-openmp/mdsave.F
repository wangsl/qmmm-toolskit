 
C     $Id: mdsave.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine mdsave  --  save trajectory and restart files  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mdsave" writes molecular dynamics trajectory snapshots and
c     auxiliary files with velocity and induced dipole information;
c     also checks for user requested termination of a simulation
c
c
      subroutine mdsave (istep,dt,epot)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mdstuf.i'
      include 'moldyn.i'
      include 'mpole.i'
      include 'output.i'
      include 'polar.i'
      include 'potent.i'
      include 'socket.i'
      include 'titles.i'
      include 'units.i'
      integer i,j,k,istep
      integer ixyz,iind
      integer ivel,ifrc
      integer iend,idump,lext
      integer freeunit,trimtext
      integer moddump
      real*8 dt,epot,pico,wt
      logical exist
      character*7 ext
      character*120 endfile
      character*120 xyzfile
      character*120 velfile
      character*120 frcfile
      character*120 indfile
c
c
c     send data via external socket communication if desired
c
      if (.not.skt_init .or. use_socket)  call sktdyn (istep,dt,epot)
c
c     test for requested termination of the dynamics calculation
c
      endfile = 'tinker.end'
      inquire (file=endfile,exist=exist)
      if (.not. exist) then
         endfile = filename(1:leng)//'.end'
         inquire (file=endfile,exist=exist)
         if (exist) then
            iend = freeunit ()
            open (unit=iend,file=endfile,status='old')
            close (unit=iend,status='delete')
         end if
      end if
      if (exist) then
         write (iout,10)
   10    format (/,' MDSAVE  --  Dynamics Calculation Ending',
     &              ' due to User Request')
         call fatal
      end if
c
c     check number of steps between trajectory file dumps
c
      moddump = mod(istep,iwrite)
      if (moddump .ne. 0)  return
c
c     get the sequence number of the current trajectory frame
c
      idump = nprior + istep/iwrite
      lext = 3
      call numeral (idump,ext,lext)
c
c     print header for the instantaneous values at current step
c
      pico = dble(istep) * dt
      write (iout,20)  istep
   20 format (/,' Instantaneous Values for Frame saved at',
     &           i10,' Dynamics Steps')
      write (iout,30)  pico
   30 format (/,' Current Time',8x,f15.4,' Picosecond')
      write (iout,40)  epot
   40 format (' Current Potential',3x,f15.4,' Kcal/mole')
      if (use_bounds) then
         write (iout,50)  xbox,ybox,zbox
   50    format (' Lattice Lengths',6x,3f14.4)
         write (iout,60)  alpha,beta,gamma
   60    format (' Lattice Angles',7x,3f14.4)
      end if
      write (iout,70)  idump
   70 format (' Frame Number',13x,i10)
c
c     update the information needed to restart the trajectory
c
      call prtdyn
c
c     save coordinates to an archive or numbered structure file
c
      ixyz = freeunit ()
      if (archive) then
         xyzfile = filename(1:leng)
         call suffix (xyzfile,'arc')
         call version (xyzfile,'old')
         inquire (file=xyzfile,exist=exist)
         if (exist) then
            call openend (ixyz,xyzfile)
         else
            open (unit=ixyz,file=xyzfile,status='new')
         end if
      else
         xyzfile = filename(1:leng)//'.'//ext(1:lext)
         call version (xyzfile,'new')
         open (unit=ixyz,file=xyzfile,status='new')
      end if
      call prtxyz (ixyz)
      close (unit=ixyz)
      write (iout,80)  xyzfile(1:trimtext(xyzfile))
   80 format (' Coordinate File',12x,a)
c
c     save the velocity vector components at the current step
c
      if (velsave) then
         ivel = freeunit ()
         if (archive) then
            velfile = filename(1:leng)
            call suffix (velfile,'vel')
            call version (velfile,'old')
            inquire (file=velfile,exist=exist)
            if (exist) then
               call openend (ivel,velfile)
            else
               open (unit=ivel,file=velfile,status='new')
            end if
         else
            velfile = filename(1:leng)//'.'//ext(1:lext)//'v'
            call version (velfile,'new')
            open (unit=ivel,file=velfile,status='new')
         end if
         write (ivel,90)  n,title(1:ltitle)
   90    format (i6,2x,a)
         do i = 1, n
            write (ivel,100)  i,name(i),(v(j,i),j=1,3)
  100       format (i6,2x,a3,3x,d13.6,3x,d13.6,3x,d13.6)
         end do
         close (unit=ivel)
         write (iout,110)  velfile(1:trimtext(velfile))
  110    format (' Velocity File',15x,a)
      end if
c
c     save the force vector components for the current step
c
      if (frcsave) then
         ifrc = freeunit ()
         if (archive) then
            frcfile = filename(1:leng)
            call suffix (frcfile,'frc')
            call version (frcfile,'old')
            inquire (file=frcfile,exist=exist)
            if (exist) then
               call openend (ifrc,frcfile)
            else
               open (unit=ifrc,file=frcfile,status='new')
            end if
         else
            frcfile = filename(1:leng)//'.'//ext(1:lext)//'f'
            call version (frcfile,'new')
            open (unit=ifrc,file=frcfile,status='new')
         end if
         write (ifrc,120)  n,title(1:ltitle)
  120    format (i6,2x,a)
         do i = 1, n
            wt = mass(i) / convert
            write (ifrc,130)  i,name(i),(wt*a(j,i),j=1,3)
  130       format (i6,2x,a3,3x,d13.6,3x,d13.6,3x,d13.6)
         end do
         close (unit=ifrc)
         write (iout,140)  frcfile(1:trimtext(frcfile))
  140    format (' Force Vector File',11x,a)
      end if
c
c     save the current induced dipole moment at each site
c
      if (uindsave .and. use_polar) then
         iind = freeunit ()
         if (archive) then
            indfile = filename(1:leng)
            call suffix (indfile,'uind')
            call version (indfile,'old')
            inquire (file=indfile,exist=exist)
            if (exist) then
               call openend (iind,indfile)
            else
               open (unit=iind,file=indfile,status='new')
            end if
         else
            indfile = filename(1:leng)//'.'//ext(1:lext)//'u'
            call version (indfile,'new')
            open (unit=iind,file=indfile,status='new')
         end if
         write (iind,150)  n,title(1:ltitle)
  150    format (i6,2x,a)
         do i = 1, npole
            if (polarity(i) .ne. 0.0d0) then
               k = ipole(i)
               write (iind,160)  k,name(k),(debye*uind(j,i),j=1,3)
  160          format (i6,2x,a3,3f12.6)
            end if
         end do
         close (unit=iind)
         write (iout,170)  indfile(1:trimtext(indfile))
  170    format (' Induced Dipole File',10x,a)
      end if
c
c     skip an extra line to keep the output formating neat
c
      moddump = mod(istep,iprint)
      if (verbose .and. moddump.ne.0) then
         write (iout,180)
  180    format ()
      end if
      return
      end
