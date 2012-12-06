 
C     $Id: sockets.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2003 by Michael Schnieders & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  routines below implement dummy versions of the socket  ##
c     ##  communication calls required for the transmission of   ##
c     ##  information between TINKER and Force Field Explorer;   ##
c     ##  functional C code is in "sockets.c", while the dummy   ##
c     ##  calls in this file are written in standard Fortran     ##
c     ##                                                         ##
c     #############################################################
c
c     ############################
c     ##                        ##
c     ##  subroutine createjvm  ##
c     ##                        ##
c     ############################
c
c
      subroutine createjvm (flag)
      implicit none
      integer flag
      return
      end
c
c
c     #############################
c     ##                         ##
c     ##  subroutine destroyjvm  ##
c     ##                         ##
c     #############################
c
c
      subroutine destroyjvm ()
      implicit none
      return
      end
c
c
c     ###############################
c     ##                           ##
c     ##  subroutine createserver  ##
c     ##                           ##
c     ###############################
c
c
      subroutine createserver (flag)
      implicit none
      integer flag
      return
      end
c
c
c     ################################
c     ##                            ##
c     ##  subroutine destroyserver  ##
c     ##                            ##
c     ################################
c
c
      subroutine destroyserver ()
      implicit none
      return
      end
c
c
c     ###############################
c     ##                           ##
c     ##  subroutine createsystem  ##
c     ##                           ##
c     ###############################
c
c
      subroutine createsystem (n,nkey,flag)
      implicit none
      integer n,nkey,flag
      return
      end
c
c
c     #############################
c     ##                         ##
c     ##  subroutine getmonitor  ##
c     ##                         ##
c     #############################
c
c
      subroutine getmonitor ()
      implicit none
      return
      end
c
c
c     #################################
c     ##                             ##
c     ##  subroutine releasemonitor  ##
c     ##                             ##
c     #################################
c
c
      subroutine releasemonitor ()
      implicit none
      return
      end
c
c
c     ###############################
c     ##                           ##
c     ##  subroutine createupdate  ##
c     ##                           ##
c     ###############################
c
c
      subroutine createupdate (n,mode,amoeba,flag)
      implicit none
      integer n,mode,amoeba,flag
      return
      end
c
c
c     #############################
c     ##                         ##
c     ##  subroutine needupdate  ##
c     ##                         ##
c     #############################
c
c
      subroutine needupdate (flag)
      implicit none
      integer flag
      return
      end
c
c
c     #############################
c     ##                         ##
c     ##  subroutine setupdated  ##
c     ##                         ##
c     #############################
c
c
      subroutine setupdated ()
      implicit none
      return
      end
c
c
c     ##########################
c     ##                      ##
c     ##  subroutine setfile  ##
c     ##                      ##
c     ##########################
c
c
      subroutine setfile (filename)
      implicit none
      character*(*) filename
      return
      end
c
c
c     ################################
c     ##                            ##
c     ##  subroutine setforcefield  ##
c     ##                            ##
c     ################################
c
c
      subroutine setforcefield (forcefield)
      implicit none
      character*(*) forcefield
      return
      end
c
c
c     #############################
c     ##                         ##
c     ##  subroutine setkeyword  ##
c     ##                         ##
c     #############################
c
c
      subroutine setkeyword (i,keyline)
      implicit none
      integer i
      character*(*) keyline
      return
      end
c
c
c     ###############################
c     ##                           ##
c     ##  subroutine setatomtypes  ##
c     ##                           ##
c     ###############################
c
c
      subroutine setatomtypes (n,type)
      implicit none
      integer n
      integer type(*)
      return
      end
c
c
c     ############################
c     ##                        ##
c     ##  subroutine setatomic  ##
c     ##                        ##
c     ############################
c
c
      subroutine setatomic (n,atomic)
      implicit none
      integer n
      integer atomic(*)
      return
      end
c
c
c     ##########################
c     ##                      ##
c     ##  subroutine setmass  ##
c     ##                      ##
c     ##########################
c
c
      subroutine setmass (n,mass)
      implicit none
      integer n
      real*8 mass(*)
      return
      end
c
c
c     ############################
c     ##                        ##
c     ##  subroutine setcharge  ##
c     ##                        ##
c     ############################
c
c
      subroutine setcharge (n,charge)
      implicit none
      integer n
      real*8 charge(*)
      return
      end
c
c
c     ##################################
c     ##                              ##
c     ##  subroutine setconnectivity  ##
c     ##                              ##
c     ##################################
c
c
      subroutine setconnectivity (n,b1,b2,b3,b4)
      implicit none
      integer n
      integer b1(*),b2(*),b3(*),b4(*)
      return
      end
c
c
c     ##########################
c     ##                      ##
c     ##  subroutine setname  ##
c     ##                      ##
c     ##########################
c
c
      subroutine setname (i,name)
      implicit none
      integer i
      character*(*) name
      return
      end
c
c
c     ###########################
c     ##                       ##
c     ##  subroutine setstory  ##
c     ##                       ##
c     ###########################
c
c
      subroutine setstory (i,story)
      implicit none
      integer i
      character*(*) story
      return
      end
c
c
c     #################################
c     ##                             ##
c     ##  subroutine setcoordinates  ##
c     ##                             ##
c     #################################
c
c
      subroutine setcoordinates (n,x,y,z)
      implicit none
      integer n
      real*8 x(*),y(*),z(*)
      return
      end
c
c
c     ##########################
c     ##                      ##
c     ##  subroutine setstep  ##
c     ##                      ##
c     ##########################
c
c
      subroutine setstep (ncycle)
      implicit none
      integer ncycle
      return
      end
c
c
c     ##########################
c     ##                      ##
c     ##  subroutine settime  ##
c     ##                      ##
c     ##########################
c
c
      subroutine settime (time)
      implicit none
      real*8 time
      return
      end
c
c
c     ############################
c     ##                        ##
c     ##  subroutine setenergy  ##
c     ##                        ##
c     ############################
c
c
      subroutine setenergy (energy)
      implicit none
      real*8 energy
      return
      end
c
c
c     ###############################
c     ##                           ##
c     ##  subroutine setgradients  ##
c     ##                           ##
c     ###############################
c
c
      subroutine setgradients (n,dx,dy,dz)
      implicit none
      integer n
      real*8 dx(*),dy(*),dz(*)
      return
      end
c
c
c     ##############################
c     ##                          ##
c     ##  subroutine setvelocity  ##
c     ##                          ##
c     ##############################
c
c
      subroutine setvelocity (n,vx,vy,vz)
      implicit none
      integer n
      real*8 vx(*),vy(*),vz(*)
      return
      end
c
c
c     ##################################
c     ##                              ##
c     ##  subroutine setacceleration  ##
c     ##                              ##
c     ##################################
c
c
      subroutine setacceleration (n,ax,ay,az)
      implicit none
      integer n
      real*8 ax(*),ay(*),az(*)
      return
      end
c
c
c     #############################
c     ##                         ##
c     ##  subroutine setinduced  ##
c     ##                         ##
c     #############################
c
c
      subroutine setinduced (n,ux,uy,uz)
      implicit none
      integer n
      real*8 ux(*),uy(*),uz(*)
      return
      end
