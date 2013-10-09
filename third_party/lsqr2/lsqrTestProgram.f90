!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File lsqrTestProgram.f90
!
!    lsqrTestProgram
!
! Main program for testing LSQR via subroutine lsqrtest in lsqrTestModule.
!
! Maintained by Michael Saunders <saunders@stanford.edu>.
!
! 24 Oct 2007: Use real(8) instead of double precision or -r8.
! 26 Oct 2012: Main program outputs some helpful info to the screen.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program lsqrTestProgram

  use   lsqrDataModule, only : dp
  use   lsqrTestModule, only : lsqrtest
  implicit none

  !---------------------------------------------------------------------
  ! This program calls lsqrtest(...) to generate a series of test problems
  ! Ax = b or Ax ~= b and solve them with LSQR.
  ! The matrix A is m x n.  It is defined by routines in lsqrTestModule.
  !
  ! 23 Sep 2007: First version of lsqrTestProgram.f90.
  ! 24 Oct 2007: Use real(dp) instead of compiler option -r8.
  ! 26 Oct 2012: Add date and cpu output.
  !---------------------------------------------------------------------

  intrinsic     :: date_and_time, cpu_time

  ! Local variables
  integer       :: ios,m,n,nbar,ndamp,nduplc,npower,nout
  real(dp)      :: time1, time2
  real(dp)      :: damp

  character(8)  :: date
  character(10) :: time
  character(80) :: output_file


  nout   = 10
  output_file = 'LSQR.txt'
  open(nout, file=output_file, status='unknown', iostat=ios)

  if (ios /= 0) then
     write(*,*)
     write(*,*) "Error opening file ", trim(output_file)
     stop
  end if

  call date_and_time( date, time )
  write(*,*)
  write(*,*) 'Date: ', date, '        Time: ', time
  call cpu_time( time1 )

  nbar   = 1000
  nduplc = 40

  m = 2*nbar        ! Over-determined systems
  n = nbar
  do ndamp = 2,7
     npower = ndamp
     damp   = 0.0_dp
     if (ndamp > 2) damp   = 10.0_dp**(-ndamp)
     call lsqrtest(m,n,nduplc,npower,damp,nout)
  end do

  m = nbar          ! Square systems
  n = nbar
  do ndamp = 2,7
     npower = ndamp
     damp   = 0.0_dp
     if (ndamp > 2) damp   = 10.0_dp**(-ndamp-6)
     call lsqrtest(m,n,nduplc,npower,damp,nout)
  end do

  m = nbar          ! Under-determined systems
  n = 2*nbar
  do ndamp = 2,6
     npower = ndamp
     damp   = 0.0_dp
     if (ndamp > 2) damp   = 10.0_dp**(-ndamp-6)
     call lsqrtest(m,n,nduplc,npower,damp,nout)
  end do

  close(nout)
  call cpu_time( time2 )

  write(*,'(a, f13.3)') " Total CPU time (seconds) ", time2-time1
  write(*,*) "Results are in output file   ", trim(output_file)
  write(*,*) "Search the file for 'appears'"
  write(*,*) "For example,    grep appears ", trim(output_file)
  write(*,*)

end program lsqrTestProgram
