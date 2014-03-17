! Quick hack to read/write stuff from/to nc files.
! For debugging only.

module ncutils
    use definitions ! real kind
    use netcdf

    implicit none

    public ncwrite
    public ncread
    private check

contains

    subroutine ncread(sns,ctns,pns,s,ct,p)

        !character (len = *), parameter :: FILE_NAME = "/home/z3439823/mymatlab/omega/ansu_utils/exp219/data/os_input.nc"
        character (len = *), parameter :: FILE_NAME = "/home/z3439823/mymatlab/omega/ansu_utils/exp229/data/os_input.nc"
        real(rk), dimension(:,:), intent(inout) :: sns, ctns, pns
        real(rk), dimension(:,:,:), intent(inout) :: s, ct, p
        integer :: ncid, varid

        call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )

        call check( nf90_inq_varid(ncid, "sns", varid) )
        call check( nf90_get_var(ncid, varid, sns) )

        call check( nf90_inq_varid(ncid, "ctns", varid) )
        call check( nf90_get_var(ncid, varid, ctns) )
        call check( nf90_inq_varid(ncid, "pns", varid) )
        call check( nf90_get_var(ncid, varid, pns) )

        call check( nf90_inq_varid(ncid, "sa", varid) )
        call check( nf90_get_var(ncid, varid, s) )
        call check( nf90_inq_varid(ncid, "ct", varid) )
        call check( nf90_get_var(ncid, varid, ct) )
        call check( nf90_inq_varid(ncid, "p", varid) )
        call check( nf90_get_var(ncid, varid, p) )

        ! Close
        call check( nf90_close(ncid) )

    end subroutine ncread


    subroutine ncwrite(va,myshape,fname,vname)
        real(rk), dimension(:), intent(in) :: va
        character (len = *), intent(in) :: fname
        character (len = *), intent(in) :: vname
        integer, dimension(:), intent(in) :: myshape
        integer :: ndims, ncid, varid, x_dimid, y_dimid, z_dimid
        integer, dimension(:), allocatable :: dimids
        integer :: n1,n2,n3
        real(rk) :: zero, fillval

        ndims=size(myshape)
        !write(*,*) 'ndims: ', ndims


        allocate(dimids(ndims))

        call check( nf90_create(fname,ior(nf90_netcdf4, nf90_classic_model),ncid))

        n1=myshape(1)
        call check( nf90_def_dim(ncid,"x",n1,x_dimid))

        if (ndims.eq.1) then
            dimids =  [x_dimid]
        endif
        if (ndims.gt.1) then
            n2=myshape(2)
            call check( nf90_def_dim(ncid,"y",n2,y_dimid))
            if (ndims.eq.2) then
                dimids =  [x_dimid, y_dimid]
            endif
            if (ndims.eq.3) then
                n3=myshape(3)
                call check( nf90_def_dim(ncid,"z",n3,z_dimid))
                dimids =  [x_dimid, y_dimid, z_dimid]
            endif
        endif

        call check( nf90_def_var(ncid,vname,NF90_DOUBLE,dimids,varid))
        call check( nf90_enddef(ncid) )


        if (ndims.eq.1) then
            call check( nf90_put_var(ncid,varid,va))
        else if (ndims.eq.2) then
            call check( nf90_put_var(ncid,varid,reshape(va,[n1,n2])))
        else if (ndims.eq.3) then
            call check( nf90_put_var(ncid,varid,reshape(va,[n1,n2,n3])))
        endif

        call check( nf90_close(ncid) )

    end subroutine ncwrite


    subroutine check(status)
        integer, intent ( in) :: status

        if(status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if
    end subroutine check
end module ncutils
