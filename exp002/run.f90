!This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
!Partially modified by P. Barker (2010-13)
!Partially modified by S. Riha (2013)
!Principal investigator: Trevor McDougall
!
!Translated to Fortran by S. Riha (2013)

program run
    use ncutils
    use ansu
    implicit none


    real(rk), dimension(nx,ny,nz) :: s, ct, p
    real(rk), dimension(nx,ny) :: sns, ctns, pns


    call ncread(sns,ctns,pns,s,ct,p)

    call optimize_surface(sns,ctns,pns,s,ct,p)

    call ncwrite(pack(sns,.true.),'sns.nc','sns',2)

    write(*,*) '****END****'

end program run


!    call ncwrite(pack(ey,.true.),'ey.nc','ey',2)

!    ex=merge(nan,ex,pns<=cut_off_choice)
!    ey=merge(nan,ey,pns<=cut_off_choice)
!    pns=merge(nan,pns,pns<=cut_off_choice)

    !test=reshape( (/1.0d0, 1.0d0, nan, nan, nan, nan, nan, nan, nan,nan,nan,nan/), shape = (/3,4/))
    !test=reshape( (/1.0d0, 1.0d0, nan, nan, nan, nan, 1.0d0,1.0d0, nan,nan,1.0d0,1.0d0/), shape = (/3,4/))
    !test=reshape( (/1.0d0, nan, 1.0d0, nan, nan, nan, nan, nan, nan, 1.0d0, nan, 1.0d0/), shape = (/3,4/))
    !test=reshape( (/1.0d0, nan, 1.0d0, nan, nan, nan, nan, nan, nan, 1.0d0, nan, 1.0d0/), shape = (/3,4/))
    !call find_regions(test,regions);

!    write(*,*) 'size(regions_test): ', size(regions_test)
!    do i=1,size(regions)
!        write(*,*) 'size(regions(',i,'): ',size(regions(i)%points)
!    enddo

    !write(*,*) 'size(regions): ', size(regions(1)%points)
!    do i=1,size(regions(1)%points)
!        write(*,*) 'regions(1)%point(', i, '): ',  int(regions(1)%points(i))
!    enddo

