 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  bath.i  --  temperature and pressure control parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     maxnose     maximum length of the Nose-Hoover chain
c
c     kelvin0     target value for the system temperature (K)
c     kelvin      variable target temperature for thermostat (K)
c     atmsph      target value for the system pressure (atm)
c     tautemp     time constant for Berendsen thermostat (psec)
c     taupres     time constant for Berendsen barostat (psec)
c     compress    isothermal compressibility of medium (atm-1)
c     collide     collision frequency for Andersen thermostat
c     xnh         position of each chained Nose-Hoover thermostat
c     vnh         velocity of each chained Nose-Hoover thermostat
c     qnh         mass for each chained Nose-Hoover thermostat
c     gnh         coupling between chained Nose-Hoover thermostats
c     isothermal  logical flag governing use of temperature control
c     isobaric    logical flag governing use of pressure control
c     tempvary    logical flag to enable variable target thermostat
c     thermostat  choice of temperature control method to be used
c     barostat    choice of pressure control method to be used
c
c
      integer maxnose
      parameter (maxnose=2)
      real*8 kelvin0,kelvin
      real*8 atmsph
      real*8 tautemp,taupres
      real*8 compress,collide
      real*8 xnh,vnh,qnh,gnh
      logical isothermal
      logical isobaric
      logical tempvary
      character*9 barostat
      character*11 thermostat
      common /bath/ kelvin0,kelvin,atmsph,tautemp,taupres,compress,
     &              collide,xnh(maxnose),vnh(maxnose),qnh(maxnose),
     &              gnh(maxnose),isothermal,isobaric,tempvary,
     &              thermostat,barostat
