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

    nit=6

    ! read input data: sns, ctns, pns are 2-d and
    ! represent sal., temp. and pressure on the initial surface.
    ! s, ct, p are the 3-d ocean hydrography
    call ncread(sns,ctns,pns,s,ct,p)

    ! write the initial surface to a netcdf file
    !call ncwrite(pack(pns,.true.),shape(pns),'pns0.nc','pns')

    ! In each iteration, the following will be written to standard out.
    ! A region is a 4-connected neighbourhood on the current surface
    ! (http://en.wikipedia.org/wiki/4-connected_neighborhood).
    write(*,*), ''
    write(*,*), 'Col. 1: Region Number'
    write(*,*), 'Col. 2: Region Size'
    write(*,*), 'Col. 3: Number of LSQR iterations'

    ! Main loop
    do it=1,nit

        write(*,*) ''
        write (*,'(a20,2x,i3,1x,a11)') &
            '******** Iteration ', it,' *********'

        ! sns, ctns and pns are input and output
        call optimize_surface(sns,ctns,pns,s,ct,p)

        ! write surfaces in each iteration
!        if (it.eq.1) then
!            call ncwrite(pack(pns,.true.),shape(pns),'pns1.nc','pns')
!        else if (it.eq.2) then
!            call ncwrite(pack(pns,.true.),shape(pns),'pns2.nc','pns')
!        else if (it.eq.3) then
!            call ncwrite(pack(pns,.true.),shape(pns),'pns3.nc','pns')
!        endif
    enddo

    ! write final surface
    call ncwrite(pack(sns,.true.),shape(sns),'sns.nc','sns')
    call ncwrite(pack(ctns,.true.),shape(ctns),'ctns.nc','ctns')
    call ncwrite(pack(pns,.true.),shape(pns),'pns.nc','pns')

    write(*,*) ''
    write(*,*) '****END****'

end program run
