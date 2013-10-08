!This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
!Partially modified by P. Barker (2010-13)
!Partially modified by S. Riha (2013)
!Principal investigator: Trevor McDougall
!
!Translated to Fortran by S. Riha (2013)

module ans_mod
    use ncutils_mod
    implicit none
    public mld

    type, public :: cell
        real, dimension(:), pointer :: points
    end type cell


    integer, parameter, public :: rk = selected_real_kind(14,30)


    type(cell), dimension(:), pointer, public :: regions
    integer, parameter, public :: NX = 90, NY = 43, NZ=101

    real(rk),  dimension(NX,NY), public :: drhox, drhoy
    real(rk),  allocatable, dimension(:) :: y
    integer,  allocatable, dimension(:) ::  j1, j2



    interface
        function gsw_rho(sa,ct,p)
            integer, parameter :: rk = selected_real_kind(14,30)
            real (rk)  :: gsw_rho
            real (rk) :: sa, ct, p
        end function gsw_rho
        function gsw_alpha(sa,ct,p)
            integer, parameter :: rk = selected_real_kind(14,30)
            real (rk)  :: gsw_alpha
            real (rk) :: sa, ct, p
        end function gsw_alpha
        function gsw_beta(sa,ct,p)
            integer, parameter :: rk = selected_real_kind(14,30)
            real (rk)  :: gsw_beta
            real (rk) :: sa, ct, p
        end function gsw_beta
    end interface


contains
    subroutine mld(s,ct,p,cut_off_choice)


        real(rk), dimension(:,:,:), intent(in) :: s, ct, p
        real(rk), dimension(:,:), intent(out) :: cut_off_choice
        real(rk), dimension(size(s,1)*size(s,2),size(s,3)) :: rho, surf_dens_, thresh
        real(rk), dimension(size(s,1)*size(s,2)*size(s,3)) :: pflat
        real(rk), dimension(size(s,1)*size(s,2)) :: surf_dens, p_, cut_off_choice_
        logical, dimension(size(s,1)*size(s,2),size(s,3)) :: pos
        integer, dimension(size(s,1)*size(s,2)) :: ip
        integer, dimension(4) :: C3
        integer :: nxy, i, j, k, setnan

        nxy=nx*ny

        do k=1,nz
            do j=1,ny
                do i=1,nx
                    rho(i+(j-1)*nx,k)=gsw_rho(s(i,j,k),ct(i,j,k),0d0*p(i,j,k))
                enddo
            enddo
        enddo

        surf_dens=rho(:,1)

!        call ncwrite(pack(surf_dens,.true.),'surf_dens.nc','surf_dens',2)

        surf_dens_= spread(surf_dens,2,nz)

!        call ncwrite(pack(surf_dens_,.true.),'surf_dens_.nc','surf_dens_',3)

        thresh= surf_dens_+0.3d0

        pos=(thresh-rho)>0.0d0

        ip=count(pos,2)

        pflat=pack(p,.true.)
        do i=1,nxy
            p_(i)=pflat(i+(ip(i)-1)*nxy)
        enddo

!        call ncwrite(p_,'p_.nc','p_',2)

        setnan=0 ! hide division by 0 at compile time
        cut_off_choice_=0d0/setnan

        cut_off_choice_=merge(p_,cut_off_choice_,ip>0)

        cut_off_choice=reshape(cut_off_choice_,[nx,ny])

!        call ncwrite(pack(cut_off_choice,.true.),'cut_off_choice.nc','cutoff',2)


    end subroutine mld

    subroutine delta_tilde_rho(sns,ctns,pns)
        real(rk), dimension(:,:), intent(in) :: sns, ctns, pns
        real(rk), dimension(size(sns,1), size(sns,2)) :: gradx_s, grady_s, gradx_ct, grady_ct
        real(rk), dimension(size(sns,1),size(sns,2)) :: r1, r2, sns_, ctns_, pmid
        real(rk), dimension(size(sns,1),size(sns,2)) :: debug
        integer :: i, j, nxy, setnan
        real(rk) :: nan
        logical :: zonally_periodic

        namelist /user_input/ zonally_periodic

        open(1,file='user_input.nml')
        read(1,user_input)
        close(1)

        setnan=0 ! hide division by 0 at compile time
        nan=0d0/setnan

        ! drhox
        pmid=0.5d0*(pns+cshift(pns,1,1))

        do j=1,ny
            do i=1,nx
                r1(i,j)=gsw_rho(sns(i,j),ctns(i,j),pmid(i,j))
            enddo
        enddo

        sns_=cshift(sns,1,1)
        ctns_=cshift(ctns,1,1)
        do j=1,ny
            do i=1,nx
                r2(i,j)=gsw_rho(sns_(i,j),ctns_(i,j),pmid(i,j))
            enddo
        enddo

        drhox=r2-r1

        if (.not.(zonally_periodic)) then
            drhox(nx,:)=nan
        endif

        ! drhoy
        pmid=0.5d0*(pns+cshift(pns,1,2))

        do j=1,ny
            do i=1,nx
                r1(i,j)=gsw_rho(sns(i,j),ctns(i,j),pmid(i,j))
            enddo
        enddo

        sns_=cshift(sns,1,2)
        ctns_=cshift(ctns,1,2)
        do j=1,ny
            do i=1,nx
                r2(i,j)=gsw_rho(sns_(i,j),ctns_(i,j),pmid(i,j))
            enddo
        enddo

        drhoy=r2-r1
        drhoy(:,ny)=nan

    end subroutine delta_tilde_rho



    subroutine epsilon_(sns,ctns,pns,e1t,e2t,ex,ey)
        real(rk), dimension(:,:), intent(in) :: sns, ctns, pns, e1t, e2t
        real(rk), dimension(:,:), intent(out) :: ex, ey
        real(rk), dimension(size(sns,1), size(sns,2)) :: gradx_s, grady_s, gradx_ct, grady_ct
        real(rk), dimension(size(sns,1),size(sns,2)) :: alpha, beta, alphax, alphay, betax, betay
        real(rk), dimension(size(sns,1),size(sns,2)) :: debug
        integer :: i, j, nxy, setnan
        real(rk) :: nan
        logical :: zonally_periodic

        namelist /user_input/ zonally_periodic

        open(1,file='user_input.nml')
        read(1,user_input)
        close(1)

!        call ncwrite(pack(sns,.true.),'e1t.nc','e1t',2)
        call grad_surf(ctns,e1t,e2t,gradx_ct,grady_ct)
        call grad_surf(sns,e1t,e2t,gradx_s,grady_s)
!        call ncwrite(pack(ctns,.true.),'e1t.nc','e1t',2)
!        call ncwrite(pack(grady_ct,.true.),'grady_ct.nc','grady_ct',2)
!        call ncwrite(pack(gradx_s,.true.),'gradx_s.nc','gradx_s',2)

        nxy=nx*ny

        do j=1,ny
            do i=1,nx
                alpha(i,j)=gsw_alpha(sns(i,j),ctns(i,j),pns(i,j))
            enddo
        enddo
        do j=1,ny
            do i=1,nx
                beta(i,j)=gsw_beta(sns(i,j),ctns(i,j),pns(i,j))
            enddo
        enddo

        alphax=0.5d0*(alpha+cshift(alpha,1,1))
        alphay=0.5d0*(alpha+cshift(alpha,1,2))
        betax=0.5d0*(beta+cshift(beta,1,1))
        betay=0.5d0*(beta+cshift(beta,1,2))

        setnan=0 ! hide division by 0 at compile time
        nan=0d0/setnan

        alphay(:,ny)=nan
        betay(:,ny)=nan

        if (.not.(zonally_periodic)) then
            alphax(nx,:)=nan
            betax(nx,:)=nan
        endif

        ex=betax*gradx_s-alphax*gradx_ct
        ey=betay*grady_s-alphay*grady_ct


    end subroutine epsilon_

    subroutine get_j()
        real(rk),  dimension(nx*ny) :: drhox_, drhoy_
        logical, dimension(nx*ny) :: reg
        integer, dimension(nx*ny) :: sreg, sreg_en, sreg_nn
        integer, allocatable, dimension(:) :: region, en, nn, j1_ew, j2_ew, j1_ns, j2_ns
        integer :: setnan, i, j
        real(rk) :: nan
        logical :: zonally_periodic

        namelist /user_input/ zonally_periodic

        open(1,file='user_input.nml')
        read(1,user_input)
        close(1)

        setnan=0 ! hide division by 0 at compile time
        nan=0d0/setnan

        drhox_=pack(drhox,.true.)
        drhoy_=pack(drhoy,.true.)

        allocate(region(size(regions(1)%points)))
        region=regions(1)%points

        reg=.false.
        do i=1,size(region)
            reg(region(i))=.true.
        enddo

        en=reg.and.cshift(reg,1)
        if (.not.zonally_periodic) then
            en(nx:ny*nx:nx)=.false.
        endif

        sreg(1)=reg(1) ! cumulative sum
        do i=2,size(reg)
            j=reg(i)
            sreg(i)=sreg(i-1)+j
        enddo
        sreg_en=cshift(sreg,1)

        allocate(j1_ew(sum(en)))
        allocate(j2_ew(sum(en)))
        allocate(y(sum(en)+sum(nn)+1)) ! A*x=y
        j=1
        do i=1,size(en)
            if (en(i)==1) then
                j1_ew(j)=sreg_en(i) ! j1 are j-indices for matrix coefficient 1
                j2_ew(j)=sreg(i) ! j2 are j-indices for matrix coefficient -1
                y(j)=drhox_(i)
                j=j+1
            endif
        enddo


        ! north-south
        nn=reg.and.cshift(reg,nx)
        nn(nx*(ny-1)+1:ny*nx:nx)=.false.

        sreg_nn=cshift(sreg,nx)


        allocate(j1_ns(sum(nn)))
        allocate(j2_ns(sum(nn)))
        j=1
        do i=1,size(nn)
            if (nn(i)==1) then
                j1_ns(j)=sreg_nn(i) ! j1 are j-indices for matrix coefficient 1
                j2_ns(j)=sreg(i) ! j2 are j-indices for matrix coefficient -1
                y(sum(en)+j)=drhoy_(i)
                j=j+1
            endif
        enddo

        allocate(j1(size(j1_ew)+size(j1_ns)))
        j1(1:size(j1_ew))=j1_ew
        j1(size(j1_ew)+1: size(j1_ew)+size(j1_ns))=j1_ns

        allocate(j2(size(j2_ew)+size(j2_ns)))
        j2(1:size(j2_ew))=j2_ew
        j2(size(j2_ew)+1: size(j2_ew)+size(j2_ns))=j2_ns

        !condition
        y(size(y))=0

    end subroutine get_j


    subroutine solve_lsqr( drho)
        real(rk),  dimension(:,:), intent(out) :: drho
        real(rk),  allocatable, dimension(:,:) :: y
        integer, allocatable, dimension(:) :: j1, j2
        real(rk), allocatable, dimension(:) :: x
        integer :: setnan, i, j
        real(rk) :: nan
        logical :: zonally_periodic

        namelist /user_input/ zonally_periodic

        open(1,file='user_input.nml')
        read(1,user_input)
        close(1)

        setnan=0 ! hide division by 0 at compile time
        nan=0d0/setnan

        call get_j()

        allocate(x(size(regions(1)%points)))

        call Aprod1(size(y),size(regions(1)%points), x,y)

!        allocate(y(size(drhox)+size(drhoy)+1))
!
!        call get_j
!        call Aprod1(size(j1)+1,size(regions(1)%points),drho,y)
!        LSQR  ( m, n, Aprod1, Aprod2, b, damp, wantse,         &
!                     x, se,                                         &
!                     atol, btol, conlim, itnlim, nout,              &
!                     istop, itn, Anorm, Acond, rnorm, Arnorm, xnorm )

    end subroutine solve_lsqr

    subroutine Aprod1(m,n,x,y)
        integer,  intent(in)    :: m,n
        real(rk), intent(in)    :: x(n)
        real(rk), intent(inout) :: y(m)
        integer :: i

        do i=1,m-1
            y(i)=y(i)+x(j1(i))-x(j2(i))
        enddo
    end subroutine Aprod1

    subroutine Aprod2(m,n,x,y)
        integer,  intent(in)    :: m,n
        real(rk), intent(inout)    :: x(n)
        real(rk), intent(in) :: y(m)
        integer, allocatable, dimension(:) :: j1, j2
        integer :: j
        real(rk), dimension(size(y)) :: tmp

        do j=1,n
            tmp=0
            where (j1==j)
                tmp=y
            end where
                x(j)=x(j)+sum(tmp)
                tmp=0
            where (j2==j)
                tmp=y
            end where
                x(j)=x(j)-sum(tmp)
                tmp=0
        enddo
    end subroutine Aprod2

    subroutine grad_surf(f,e1t,e2t,fx,fy)
        real(rk),  dimension(:,:), intent(in) :: f, e1t, e2t
        real(rk),  dimension(:,:), intent(out) :: fx, fy
        integer :: setnan
        real(rk) :: nan
        logical :: zonally_periodic

        namelist /user_input/ zonally_periodic

        open(1,file='user_input.nml')
        read(1,user_input)
        close(1)

        fx=(cshift(f,1,1)-f)/e1t


        setnan=0 ! hide division by 0 at compile time
        nan=0d0/setnan
        if (.not.(zonally_periodic)) then
            fx(size(f,1),:)=nan
        endif

        fy=(cshift(f,1,2)-f)/e2t

!        call ncwrite(pack(fy,.true.),'fy.nc','fy',2)


    end subroutine grad_surf


    subroutine find_regions(pns)
        real(rk),  dimension(:,:), intent(in) :: pns
        logical,  dimension(size(pns,1),size(pns,2)) :: wet
        logical,  dimension(size(pns,1)*size(pns,2)) :: bool, wet_
        integer,  dimension(size(pns,1)*size(pns,2)) :: L_
        integer, dimension(size(pns,2)) :: east, west
        integer :: k, i, j, ii, iregion,kk,kkk
        integer, allocatable, dimension(:) :: idx, neighbours, neighbours_tmp
        logical :: zonally_periodic

        namelist /user_input/ zonally_periodic

        open(1,file='user_input.nml')
        read(1,user_input)
        close(1)


        wet=.not.(isnan(pns))
        wet_=pack(wet,.true.)

        L_=0 ! label matrix
        iregion=1 ! region label

        do while (.true.)
            if (allocated(idx)) then
                deallocate(idx)
            endif
            allocate(idx(1))
            ! find index of first wet point
            do i = 1, size(wet_)
                if (wet_(i)) then
                    ii = i
                    exit
                endif
            end do

            if (i==size(wet_)+1) then
                exit
            endif

            idx(1)=ii ! linear indices of the pixels that were just labeled
            wet_(idx(1))=.false. ! set the pixel to 0
            L_(idx(1)) = iregion ! label first point

            do while (size(idx)/=0) ! find neighbours
                if (allocated(neighbours)) then
                    deallocate(neighbours)
                endif
                allocate(neighbours(size(idx)*4)) ! each point has 4 neighbours
                neighbours(1:size(idx))=idx+1
                neighbours(1*size(idx)+1:2*size(idx))=idx-1
                neighbours(2*size(idx)+1:3*size(idx))=idx+nx
                neighbours(3*size(idx)+1:4*size(idx))=idx-nx

                ! zonal boundaries
                east= (/(i, i=0   , (ny-1)*nx, nx)/) ! eastern neighbours of eastern bdy
                west= (/(i, i=nx+1, ny*nx+1  , nx)/) ! western neighbours of western bdy
                allocate(neighbours_tmp(size(idx)))
                if (zonally_periodic) then
                    do i=1,ny
                        where (neighbours(1*size(idx)+1:2*size(idx)) == east(ny-i+1))
                            neighbours(1*size(idx)+1:2*size(idx))=neighbours(1*size(idx)+1:2*size(idx))+nx
                        end where
                        where (neighbours(1:size(idx)) == west(i))
                            neighbours(1:size(idx))=neighbours(1:size(idx))-nx
                        end where
                    enddo
                else
                    do i=1,ny
                        where (neighbours(1*size(idx)+1:2*size(idx)) == east(i))
                            neighbours(1*size(idx)+1:2*size(idx))=-2*nx
                        end where
                        where (neighbours(1:size(idx)) == west(i))
                            neighbours(1:size(idx))=-2*nx
                        end where
                    enddo
                endif
                deallocate(neighbours_tmp)

                ! meridional boundaries
                where (neighbours<1)
                    neighbours=-2*nx ! flagging as invalid
                end where
                where (neighbours>nx*ny)
                    neighbours=-2*nx
                end where

                allocate(neighbours_tmp(count(.not.(neighbours==-2*nx))))
                j=1
                do i=1,size(neighbours)
                    if (neighbours(i)/=-2*nx) then
                        neighbours_tmp(j)=neighbours(i)
                        j=j+1
                    endif
                enddo
                deallocate(neighbours)
                allocate(neighbours(size(neighbours_tmp)))
                neighbours=neighbours_tmp
                deallocate(neighbours_tmp)

                ! Remove duplicate entries. This is copied from
                ! http://rosettacode.org/wiki/Remove_duplicate_elements#Fortran
                allocate(neighbours_tmp(size(neighbours)))

                if (size(neighbours)>0) then
                    neighbours_tmp(1) = neighbours(1)

                    k = 1
                    outer: do i=2,size(neighbours)
                        do j=1,k
                            if (neighbours_tmp(j) == neighbours(i)) then
                              ! Found a match so start looking again
                                cycle outer
                            end if
                        end do
                        ! No match found so add it to the output
                        k = k + 1
                        neighbours_tmp(k) = neighbours(i)
                    end do outer
                    deallocate(neighbours)
                    allocate(neighbours(k))
                    neighbours=neighbours_tmp(1:k)
                endif
                deallocate(neighbours_tmp)
                ! end remove duplicate entries

                ! keep only wet neighbours
                do i=1,size(neighbours)
                    if (.not.(wet_(neighbours(i)))) then
                        neighbours(i)=-2*nx ! flagging as invalid
                    endif
                enddo

                allocate(neighbours_tmp(count(.not.(neighbours==-2*nx))))
                    j=1
                    do i=1,size(neighbours)
                        if (neighbours(i)/=-2*nx) then
                            neighbours_tmp(j)=neighbours(i)
                            j=j+1
                        endif
                    enddo
                deallocate(neighbours)
                allocate(neighbours(size(neighbours_tmp)))
                neighbours=neighbours_tmp
                deallocate(neighbours_tmp)
                deallocate(idx)
                allocate(idx(size(neighbours)))
                idx=neighbours
                deallocate(neighbours)

                ! done?
                if (size(idx)==0) then
                    exit
                endif

                do i=1,size(idx)
                    L_(idx(i))=iregion
                    wet_(idx(i))=.false. ! set the pixels that were just labeled to 0
                enddo
            enddo

            iregion=iregion+1
        enddo


        allocate(regions(iregion-1))
        bool=.false.

        do k=1,iregion-1
            j=1
            bool=.false.
            where (L_==k)
                bool=.true.
            end where
            allocate(regions(k)%points(1:count(bool)))
            do i=1,size(L_)
                if (L_(i)==k) then
                    regions(k)%points(j)=i
                    j=j+1
                endif
            enddo
        enddo
    end subroutine find_regions


end module ans_mod
