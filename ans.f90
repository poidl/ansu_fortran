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


    interface
        function gsw_rho(sa,ct,p)
            integer, parameter :: r14 = selected_real_kind(14,30)
            real (r14)  :: gsw_rho
            real (r14) :: sa, ct, p
        end function gsw_rho
        function gsw_alpha(sa,ct,p)
            integer, parameter :: r14 = selected_real_kind(14,30)
            real (r14)  :: gsw_alpha
            real (r14) :: sa, ct, p
        end function gsw_alpha
        function gsw_beta(sa,ct,p)
            integer, parameter :: r14 = selected_real_kind(14,30)
            real (r14)  :: gsw_beta
            real (r14) :: sa, ct, p
        end function gsw_beta
    end interface

    type cell
        real, dimension(:), pointer :: points
    end type cell
contains
    subroutine mld(s,ct,p,cut_off_choice)

        integer, parameter :: rk = selected_real_kind(14,30)
        real(rk), dimension(:,:,:), intent(in) :: s, ct, p
        real(rk), dimension(:,:), intent(out) :: cut_off_choice
        real(rk), dimension(size(s,1)*size(s,2),size(s,3)) :: rho, surf_dens_, thresh
        real(rk), dimension(size(s,1)*size(s,2)*size(s,3)) :: pflat
        real(rk), dimension(size(s,1)*size(s,2)) :: surf_dens, p_, cut_off_choice_
        logical, dimension(size(s,1)*size(s,2),size(s,3)) :: pos
        integer, dimension(size(s,1)*size(s,2)) :: ip
        integer, dimension(4) :: C3
        integer :: nx, ny, nz, nxy, i, j, k, setnan

        nx=size(s,1)
        ny=size(s,2)
        nz=size(s,3)
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

        call ncwrite(pack(cut_off_choice,.true.),'cut_off_choice.nc','cutoff',2)


    end subroutine mld



    subroutine epsilon_(sns,ctns,pns,e1t,e2t,ex,ey)
        integer, parameter :: rk = selected_real_kind(14,30)
        real(rk), dimension(:,:), intent(in) :: sns, ctns, pns, e1t, e2t
        real(rk), dimension(:,:), intent(out) :: ex, ey
        real(rk), dimension(size(sns,1), size(sns,2)) :: gradx_s, grady_s, gradx_ct, grady_ct
        real(rk), dimension(size(sns,1),size(sns,2)) :: alpha, beta, alphax, alphay, betax, betay
        real(rk), dimension(size(sns,1),size(sns,2)) :: debug
        integer :: i, j, nx, ny, nxy, setnan
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

        nx=size(sns,1)
        ny=size(sns,2)
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



    subroutine grad_surf(f,e1t,e2t,fx,fy)
        integer, parameter :: rk = selected_real_kind(14,30)
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


    subroutine find_regions(pns, regions)
        integer, parameter :: rk = selected_real_kind(14,30)
        real(rk),  dimension(:,:), intent(in) :: pns
        logical,  dimension(size(pns,1),size(pns,2)) :: wet
        logical,  dimension(size(pns,1)*size(pns,2)) :: bool, wet_
        integer,  dimension(size(pns,1)*size(pns,2)) :: L_
        integer, dimension(size(pns,2)) :: east, west
        integer :: k, i, j, ii, iregion,kk,kkk, xn,yn
        integer, allocatable, dimension(:) :: idx, neighbours, neighbours_tmp
        type(cell), dimension(:), pointer, intent(out) :: regions
        logical :: zonally_periodic

        namelist /user_input/ zonally_periodic

        open(1,file='user_input.nml')
        read(1,user_input)
        close(1)

        xn=size(pns,1)
        yn=size(pns,2)

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
                neighbours(2*size(idx)+1:3*size(idx))=idx+xn
                neighbours(3*size(idx)+1:4*size(idx))=idx-xn

                ! zonal boundaries
                east= (/(i, i=0   , (yn-1)*xn, xn)/) ! eastern neighbours of eastern bdy
                west= (/(i, i=xn+1, yn*xn+1  , xn)/) ! western neighbours of western bdy
                allocate(neighbours_tmp(size(idx)))
                if (zonally_periodic) then
                    do i=1,yn
                        where (neighbours(1*size(idx)+1:2*size(idx)) == east(yn-i+1))
                            neighbours(1*size(idx)+1:2*size(idx))=neighbours(1*size(idx)+1:2*size(idx))+xn
                        end where
                        where (neighbours(1:size(idx)) == west(i))
                            neighbours(1:size(idx))=neighbours(1:size(idx))-xn
                        end where
                    enddo
                else
                    do i=1,yn
                        where (neighbours(1*size(idx)+1:2*size(idx)) == east(i))
                            neighbours(1*size(idx)+1:2*size(idx))=-2*xn
                        end where
                        where (neighbours(1:size(idx)) == west(i))
                            neighbours(1:size(idx))=-2*xn
                        end where
                    enddo
                endif
                deallocate(neighbours_tmp)

                ! meridional boundaries
                where (neighbours<1)
                    neighbours=-2*xn ! flagging as invalid
                end where
                where (neighbours>xn*yn)
                    neighbours=-2*xn
                end where

                allocate(neighbours_tmp(count(.not.(neighbours==-2*xn))))
                j=1
                do i=1,size(neighbours)
                    if (neighbours(i)/=-2*xn) then
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
                if (allocated(neighbours_tmp)) then
                    deallocate(neighbours_tmp)
                endif
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
                ! end remove duplicate entries

                ! keep only wet neighbours
                do i=1,size(neighbours)
                    if (.not.(wet_(neighbours(i)))) then
                        neighbours(i)=-2*xn ! flagging as invalid
                    endif
                enddo

                if (allocated(neighbours_tmp)) then
                    deallocate(neighbours_tmp)
                endif
                allocate(neighbours_tmp(count(.not.(neighbours==-2*xn))))
                    j=1
                    do i=1,size(neighbours)
                        if (neighbours(i)/=-2*xn) then
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
