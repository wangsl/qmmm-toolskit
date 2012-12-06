 
C     $Id: mutate.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mutate  --  set parameters for hybrid system  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mutate" constructs the hybrid hamiltonian for a specified
c     initial state, final state and mutation parameter "lambda"
c
c
      subroutine mutate
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'katoms.i'
      include 'keys.i'
      include 'mutant.i'
      integer i,ihyb
      integer it0,it1,next
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     zero number of hybrid atoms and hybrid atom list
c
      nhybrid = 0
      do i = 1, n
         alter(i) = .false.
      end do
c
c     search keywords for free energy perturbation options
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'LAMBDA ') then
            nhybrid = 0
            lambda = 0.0d0
            string = record(next:120)
            read (string,*,err=10)  lambda
         else if (keyword(1:7) .eq. 'MUTATE ') then
            string = record(next:120)
            read (string,*,err=10)  ihyb,it0,it1
            nhybrid = nhybrid + 1
            ihybrid(nhybrid) = ihyb
            type0(nhybrid) = it0
            type1(nhybrid) = it1
            class0(nhybrid) = atmcls(it0)
            class1(nhybrid) = atmcls(it1)
            alter(ihyb) = .true.
         end if
   10    continue
      end do
      return
      end
