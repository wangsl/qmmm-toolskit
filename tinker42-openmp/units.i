 
C     $Id$
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  units.i  --  physical constants and unit conversions  ##
c     ##                                                        ##
c     ############################################################
c
c
c     avogadro    Avogadro's number (N) in particles/mole
c     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/K/mole
c     gasconst    ideal gas constant (R) in kcal/mole/K
c     lightspd    speed of light in vacuum (c) in cm/ps
c     bohr        conversion from Bohrs to Angstroms
c     joule       conversion from calories to joules
c     evolt       conversion from Hartree to electron-volts
c     hartree     conversion from Hartree to kcal/mole
c     electric    conversion from electron**2/Ang to kcal/mole
c     debye       conversion from electron-Ang to Debyes
c     prescon     conversion from kcal/mole/Ang**3 to Atm
c     convert     conversion from kcal to g*Ang**2/ps**2
c
c
      real*8 avogadro,boltzmann
      real*8 gasconst,lightspd
      real*8 bohr,joule,evolt
      real*8 hartree,electric
      real*8 debye,prescon,convert
      parameter (avogadro=6.02214199d+23)
      parameter (boltzmann=0.83143435d0)
      parameter (gasconst=1.9872065d-3)
      parameter (lightspd=2.99792458d-2)
      parameter (bohr=0.5291772083d0)
      parameter (joule=4.184d0)
      parameter (evolt=27.2113834d0)
      parameter (hartree=627.5094709d0)
      parameter (electric=332.05382d0)
      parameter (debye=4.8033324d0)
      parameter (prescon=6.85695d+4)
      parameter (convert=4.184d+2)
