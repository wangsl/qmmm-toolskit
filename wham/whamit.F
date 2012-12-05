
c     Note WHAM iteration by Shenglong Wang Oct. 6, 2005
c     Reference: 
c     (1) Marc Souaille and Benoit Roux  Comput. Phys. Commun. 135, 40 (2001)
c     (2) Benoit Roux  Comput. Phys. Commun. 91, 275 (1995)
c     (3) Erik M. Boczko and Charles L. Brooks, III  J. Phys. Chem. 97, 4509 (1993)
c     (4) Alan  Grossfield  http://dasher.wustl.edu/alan/wham/
c     (5) Shankar Kumar et.al  J. Comput. Chem. 13, 1011 (1992)

      subroutine wham_iteration_with_fortran(ebw, Nwind, Niter, 
     $     tolerance, nt, maxnt, kBT, ebf, ebf2,
     $     IterationOutput, LenIterationOutput)
      implicit none
      integer Nwind, Niter, maxnt
      Real*8 tolerance
      Integer nt(Nwind)
c     Please note ebw is saved in row form (c)
c     I use the pointer from wham_iteration function directly
c     without any modification
c     Nov. 16, 2005 Shenglong Wang
      Real*8 ebw(Nwind, maxnt, Nwind)

      Real*8 ebf(Nwind)
      Real*8 ebf2(Nwind)
      Real*8 kBT(Nwind)

      Character*(*) IterationOutput
      Integer LenIterationOutput

      Real*8 fact(Nwind)
      integer i, j, k, l, n, max_n
      Real*8 ebfk, bottom, delta
      logical converged
      Real*8 cputime, start_time, end_time, elapsed_time
      
 100  Format('Start to run WHAM iterations in Fortran ...',/)
 900  Format('Maxnt error')
 1000 Format(/,'Fortran WHAM iteration step: ', i6)
 1100 Format('Window index: ', i4, '  ebf: ', 1pe16.8, '  Error: ',
     $     1pe16.8)
 1200 Format(/, '*** Congratulations: WHAM finished succefully ',
     $     'from Fortran ***',/)
 1300 Format(/, 75('='), /,
     $     'Sorry, WHAM failed with maximum iterations',
     $     i4, ' and tolerance', 1pe12.4 /,
     $     'Please restart WHAM iterations ',
     $     'with the output ebf from this iterations', //,
     $     '** The exp(beta*f) have been written to file "', A, '"',
     $     /, 75('='), /)
 1400 Format('Finished WHAM iteration with elapsed CPU time ', 
     $     f8.2, ' seconds', /)

      start_time = cputime()
      
      write(6, 100)
      call write_header
      
      max_n = -100000
      do k = 1, Nwind
         if(nt(k) .gt. max_n) max_n = nt(k)
      end do
      
      if(maxnt .lt. max_n) then
         write(6, 900)
         stop
      end if

      do k = 1, Nwind
         fact(k) = nt(k)*ebf(k)
         ebf2(k) = 1.0d0/ebf(k)
      end do

      converged = .false.
      n = 0
      
      do while(.not. converged .and. n .lt. Niter)

         do k = 1, Nwind
            ebfk = 0.0d0
c$omp parallel do schedule(static, 1)
c$omp& default(shared)
c$omp& private(i, j, l, bottom)
c$omp& firstprivate(k)
c$omp& reduction(+:ebfk)
            do i = 1, Nwind
               do l = 1, nt(i)
                  bottom = 0.0d0
                  do j = 1, Nwind
                     bottom = bottom + ebw(j,l,i)*fact(j)
                  end do
                  ebfk = ebfk + ebw(k,l,i)/bottom
               end do
            end do
c$omp end parallel do
            ebf2(k) = ebfk
            ebf(k) = ebf(k)/ebf(1)
            fact(k) = nt(k)/ebf2(k)*ebf2(1)
         end do

         n = n+1
         write(6, 1000) n
         
         converged = .true.
         do k = 1, Nwind
            delta = abs(kBT(k)*log(ebf(k)*ebf2(k)))
            if (delta .ge. tolerance) converged = .false.
            ebf(k) = ebf2(1)/ebf2(k)
            write(6, 1100) k, ebf(k), delta
         end do
         
      end do
      
      if(converged) then
         do k = 1, Nwind
            ebf2(k) = 1.0d0/ebf(k)
         end do
         write(6, 1200)
      else 
         write(6, 1300)  Niter, tolerance,
     $        IterationOutput(1:LenIterationOutput)
      end if
      
      Open(Unit = 20, File = IterationOutput(1:LenIterationOutput), 
     $     Status = 'Unknown')
      Do k = 1, Nwind
         write(20, '(1pe15.8)') ebf(k)
      End Do
      Close(20)

      end_time = cputime()
      elapsed_time = end_time - start_time
      write(6, 1400) elapsed_time

      Call Flush(6)

      return
      end
            
