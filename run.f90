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

    real(rk), dimension(NX,NY) :: sns, ctns, pns, cut_off_choice
    real(rk), dimension(NX,NY) :: e1t, e2t, ex, ey, drho
    real(rk), dimension(NX*NY) :: x
    real(rk), dimension(NX,NY,NZ) :: s, ct, p
    real(rk) :: nan
    real(rk), dimension(3,4) :: test
    integer :: i, setnan
    type(cell), dimension(:), pointer :: regions_test

    setnan=0 ! hide division by 0 at compile time
    nan=0d0/setnan

    call ncread(sns,ctns,pns,s,ct,p)
    !call ncread_debug(ctns)
!    write(*,'(A, F20.16)') 'hoit: ', sns(90,43)
!    call ncwrite_debug(pack(ctns,.true.),'ctns.nc','ctns',2)
    call mld(s,ct,p,cut_off_choice)
    call ncwrite(pack(cut_off_choice,.true.),'cut_off_choice.nc','cutoff',2)

    call delta_tilde_rho(sns,ctns,pns)
    call ncwrite(pack(drhox,.true.),'drhox.nc','drhox',2)
    call ncwrite(pack(drhoy,.true.),'drhoy.nc','drhoy',2)

!    call ncwrite(pack(ey,.true.),'ey.nc','ey',2)

!    ex=merge(nan,ex,pns<=cut_off_choice)
!    ey=merge(nan,ey,pns<=cut_off_choice)
!    pns=merge(nan,pns,pns<=cut_off_choice)

    !test=reshape( (/1.0d0, 1.0d0, nan, nan, nan, nan, nan, nan, nan,nan,nan,nan/), shape = (/3,4/))
    !test=reshape( (/1.0d0, 1.0d0, nan, nan, nan, nan, 1.0d0,1.0d0, nan,nan,1.0d0,1.0d0/), shape = (/3,4/))
    !test=reshape( (/1.0d0, nan, 1.0d0, nan, nan, nan, nan, nan, nan, 1.0d0, nan, 1.0d0/), shape = (/3,4/))
    !test=reshape( (/1.0d0, nan, 1.0d0, nan, nan, nan, nan, nan, nan, 1.0d0, nan, 1.0d0/), shape = (/3,4/))
    !call find_regions(test,regions);

    call find_regions(pns);

    allocate(regions_test(1))
    allocate(regions_test(1)%points(size(regions(1)%points)))
    regions_test(1)%points=regions(1)%points


    call solve_lsqr(drho)

    write(*,*) 'size(regions_test): ', size(regions_test)
    do i=1,size(regions_test)
        write(*,*) 'size(regions_test(',i,'): ',size(regions_test(i)%points)
    enddo

    !write(*,*) 'size(regions): ', size(regions(1)%points)
    !write(*,*) 'regions(5)%points(4): ',  regions(5)%points(4)
    write(*,'(A, F20.16)') 'hoit: ', sns(50,6)

end program run

