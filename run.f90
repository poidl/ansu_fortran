! Main program. Reads initial surface and 3-d data set from netcdf
! file, and calculates NIT iterations.

program run
    use ncutils
    use ansu
    implicit none

    integer, parameter :: nx = 90, ny = 43, nz=101

    real(rk), dimension(nx,ny,nz) :: s, ct, p
    real(rk), dimension(nx,ny) :: sns, ctns, pns
    integer :: it, nit

    nit=15

    call ncread(sns,ctns,pns,s,ct,p)

    call ncwrite(pack(pns,.true.),shape(pns),'pns0.nc','pns')

    do it=1,nit

        write(*,*) 'Iteration', it

        call optimize_surface(sns,ctns,pns,s,ct,p)

        if (it.eq.1) then
            call ncwrite(pack(pns,.true.),shape(pns),'pns1.nc','pns')
        else if (it.eq.2) then
            call ncwrite(pack(pns,.true.),shape(pns),'pns2.nc','pns')
        else if (it.eq.3) then
            call ncwrite(pack(pns,.true.),shape(pns),'pns3.nc','pns')
        endif
    enddo


!    call ncwrite(pack(sns,.true.),shape(sns),'sns.nc','sns')
!    call ncwrite(pack(ctns,.true.),shape(ctns),'ctns.nc','ctns')
!    call ncwrite(pack(pns,.true.),shape(pns),'pns.nc','pns')


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

