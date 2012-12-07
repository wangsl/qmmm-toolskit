 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2002  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  socket.i  --  control parameters for socket communication  ##
c     ##                                                             ##
c     #################################################################
c
c
c     runtyp      calculation type for passing socket information
c     cstep       current optimization or dynamics step number
c     cdt         current dynamics cumulative simulation time
c     cenergy     current potential energy from simulation
c     cdx         current gradient components along the x-axis
c     cdy         current gradient components along the y-axis
c     cdz         current gradient components along the z-axis
c     skt_init    logical flag set to true after socket initialization
c     use_socket  logical flag governing use of external sockets
c     use_gui     logical flag to show TINKER was invoked from GUI
c     closing     logical flag to indicate JVM and server shutdown
c
c
      integer runtyp,cstep
      real*8  cdt,cenergy
      real*8  cdx,cdy,cdz
      logical skt_init
      logical use_socket
      logical use_gui
      logical closing
      common /socket/ runtyp,cstep,cdt,cenergy,cdx(maxatm),cdy(maxatm),
     &                cdz(maxatm),skt_init,use_socket,use_gui,closing
