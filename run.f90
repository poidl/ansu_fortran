!This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
!Partially modified by P. Barker (2010-13)
!Partially modified by S. Riha (2013)
!Principal investigator: Trevor McDougall
!
!Translated to Fortran by S. Riha (2013)

program run
    use ncutils_mod
    use ans_mod
    implicit none

    integer, parameter :: rk = selected_real_kind(14,30)
    integer, parameter :: NX = 90, NY = 43, NZ=101
    real(rk), dimension(NX,NY) :: sns, ctns, pns, cut_off_choice
    real(rk), dimension(NX,NY) :: e1t, e2t, ex, ey
    real(rk), dimension(NX,NY,NZ) :: s, ct, p
    real(rk) :: zero, nan
    real(rk), dimension(3,4) :: test
    integer :: i
    type(cell), dimension(:), pointer :: regions


    call ncread(sns,ctns,pns,s,ct,p,e1t,e2t)
    !call ncread_debug(ctns)
!    write(*,'(A, F20.16)') 'hoit: ', sns(90,43)
!    call ncwrite_debug(pack(ctns,.true.),'ctns.nc','ctns',2)
    call mld(s,ct,p,cut_off_choice)
    call epsilon_(sns,ctns,pns,e1t,e2t,ex,ey)
    call ncwrite(pack(ex,.true.),'ex.nc','ex',2)
!    call ncwrite(pack(ey,.true.),'ey.nc','ey',2)
    zero=0.0d0
    nan = zero/ zero
    ex=merge(nan,ex,pns<=cut_off_choice)
    ey=merge(nan,ey,pns<=cut_off_choice)
    pns=merge(nan,pns,pns<=cut_off_choice)



    !test=reshape( (/1.0d0, 1.0d0, nan, nan, nan, nan, nan, nan, nan,nan,nan,nan/), shape = (/3,4/))
    !test=reshape( (/1.0d0, 1.0d0, nan, nan, nan, nan, 1.0d0,1.0d0, nan,nan,1.0d0,1.0d0/), shape = (/3,4/))
    !test=reshape( (/1.0d0, nan, 1.0d0, nan, nan, nan, nan, nan, nan, 1.0d0, nan, 1.0d0/), shape = (/3,4/))
    !test=reshape( (/1.0d0, nan, 1.0d0, nan, nan, nan, nan, nan, nan, 1.0d0, nan, 1.0d0/), shape = (/3,4/))
    !call find_regions(test,regions);

    call find_regions(pns,regions);

!    write(*,*) 'size(regions): ', size(regions)
!    do i=1,size(regions)
!        write(*,*) 'size(regions(',i,'): ',size(regions(i)%points)
!    enddo

    !write(*,*) 'size(regions): ', size(regions(1)%points)
    !write(*,*) 'regions(5)%points(4): ',  regions(5)%points(4)
    !write(*,'(A, F20.16)') 'hoit: ', sns(50,6)

end program run

