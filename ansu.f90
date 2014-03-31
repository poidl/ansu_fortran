!This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
!Modified by P. Barker (2010-2013)
!Modified by S. Riha (2013-2014)
!Principal investigator: Trevor McDougall
!
!Translated to Fortran by S. Riha (2013-2014)

module ansu
    use definitions
    use lsqrModule, only : LSQR
    use ncutils
    implicit none

!=============================================================
    public :: optimize_surface, mld, delta_tilde_rho, find_regions, &
               solve_lsqr, dz_from_drho, wetting

!=============================================================
    private

    real(rk),  allocatable, dimension(:) :: y
    integer,  allocatable, dimension(:) ::  j1, j2

!=============================================================

    interface
        function gsw_rho(sa,ct,p)
            integer, parameter :: rk = selected_real_kind(14,30)
            real (rk)  :: gsw_rho
            real (rk) :: sa, ct, p
        end function gsw_rho
!        function gsw_rho_vectorized(sa,ct,p)
!            integer, parameter :: rk = selected_real_kind(14,30)
!            real (rk), dimension(:) :: gsw_rho_vectorized
!            real (rk), dimension(:) :: sa, ct, p
!        end function gsw_rho_vectorized
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

    subroutine optimize_surface(sns,ctns,pns,s,ct,p)
        ! Main function. All other functions in this module could, in principle,
        ! be private (but are not for debugging purposes).

        real(rk), dimension(:,:), intent(inout) :: sns, ctns, pns
        real(rk), dimension(:,:,:), intent(in) :: s, ct, p
        type(region_type), dimension(:), allocatable :: regions
        real(rk), dimension(:,:), allocatable :: cut_off_choice, drhox, drhoy, drho
        integer :: nneighbours, nxy, nx, ny

        nx=size(s,1)
        ny=size(s,2)
        nxy=nx*ny
        allocate(cut_off_choice(nx,ny))
        allocate(drhox(nx,ny))
        allocate(drhoy(nx,ny))
        allocate(drho(nx,ny))

        ! disregard data above mixed layer depth
        !call mld(s,ct,p,cut_off_choice)

        ! lateral extension of outcropping/undercropping surface
        call wetting_simple(sns,ctns,pns,s,ct,p,nneighbours)
        write(*,'(a34,1x,i4)')  &
            'Wetting...number of added points:', nneighbours

        ! compute lateral gradient
        call delta_tilde_rho(sns,ctns,pns,drhox,drhoy)

        ! find regions
        call find_regions(pns,regions)

        ! lateral integration
        call solve_lsqr(regions,drhox,drhoy,drho)

        ! vertical integration
        call dz_from_drho(sns, ctns, pns, s, ct, p, drho)

    end subroutine optimize_surface


    subroutine dz_from_drho(sns, ctns, pns, s, ct, p, drho)
        real(rk), dimension(:,:), intent(inout) :: sns, ctns, pns
        real(rk), dimension(:,:,:), intent(in) :: s, ct, p
        real(rk), dimension(:,:), intent(in) :: drho
        real(rk), dimension(:,:), allocatable :: s_, ct_, p_
        real(rk), dimension(:), allocatable :: s0, ct0, p0, drho0, sns_,ctns_,pns_
        integer :: nxy, nx, ny, nz

        nx=size(s,1)
        ny=size(s,2)
        nz=size(s,3)
        nxy=nx*ny

        allocate(s_(nxy,nz))
        allocate(ct_(nxy,nz))
        allocate(p_(nxy,nz))

        allocate(s0(nxy))
        allocate(ct0(nxy))
        allocate(p0(nxy))
        allocate(drho0(nxy))

        allocate(sns_(nxy))
        allocate(ctns_(nxy))
        allocate(pns_(nxy))

        s0=pack(sns,.true.)
        ct0=pack(ctns,.true.)
        p0=pack(pns,.true.)
        drho0=pack(drho,.true.)

        s_=reshape(s,[nx*ny,nz])
        ct_=reshape(ct,[nx*ny,nz])
        p_=reshape(p,[nx*ny,nz])

        !call depth_ntp_iter_drho(s0,ct0,p0,s_,ct_,p_,drho0,sns_,ctns_,pns_)
        call depth_ntp_simple(s0,ct0,p0,s_,ct_,p_,drho0,sns_,ctns_,pns_)

        sns=reshape(sns_,[nx,ny])
        ctns=reshape(ctns_,[nx,ny])
        pns=reshape(pns_,[nx,ny])

    end subroutine dz_from_drho

    subroutine wetting(sns,ctns,pns,s,ct,p,nneighbours)
        real(rk), dimension(:,:), intent(inout) :: sns, ctns, pns
        real(rk), dimension(:,:,:), intent(in) :: s, ct, p
        integer, intent(out) :: nneighbours

        real(rk), dimension(:), allocatable :: sns_, ctns_, pns_, s1
        real(rk), dimension(:,:), allocatable :: s_, ct_, p_
        real(rk), dimension(:), allocatable :: sns_new, ctns_new, pns_new
        logical, dimension(:), allocatable :: wet, wets, nn, sn, en, wn, nbr
        integer, dimension(:), allocatable :: inds, inds_nbr
        integer, dimension(:), allocatable :: inb, isb
        integer, dimension(:), allocatable :: ieb, iwb
        integer, dimension(:), allocatable :: inbr, iwo
        integer :: i, nxy, nx, ny, nz
        logical :: zonally_periodic

        namelist /domain/ zonally_periodic

        open(1,file='user_input.nml')
        read(1,domain)
        close(1)

        nx=size(s,1)
        ny=size(s,2)
        nz=size(s,3)
        nxy=nx*ny

        allocate(sns_(nxy))
        allocate(ctns_(nxy))
        allocate(pns_(nxy))
        allocate(s1(nxy))

        allocate(s_(nxy,nz))
        allocate(ct_(nxy,nz))
        allocate(p_(nxy,nz))

        allocate(wet(nxy))
        allocate(wets(nxy))
        allocate(nn(nxy))
        allocate(sn(nxy))
        allocate(en(nxy))
        allocate(wn(nxy))
        allocate(nbr(nxy))

        allocate(inb(nx))
        allocate(isb(nx))
        allocate(ieb(ny))
        allocate(iwb(ny))

        sns_=pack(sns,.true.)
        ctns_=pack(ctns,.true.)
        pns_=pack(pns,.true.)

        s_=reshape(s,[nxy,nz])
        ct_=reshape(ct,[nxy,nz])
        p_=reshape(p,[nxy,nz])

        s1=pack(s(:,:,1),.true.)

        wet=(.not.(isnan(s1))).and.(isnan(sns_)) ! wet points at ocean surface excluding ans
        wets=(.not.isnan(sns_)) ! wet points on ans

        inb=[(i, i= (ny-1)*nx+1, ny*nx)] ! indices of northern boundary
        isb=[(i, i= 1, nx)]
        ieb=[(i, i= nx, nxy, nx)]
        iwb=[(i, i= 1, (ny-1)*nx+1, nx)]

        nn=(wet).and.(cshift(wets,nx)) ! wet points with northern neighbour on ans
        nn(inb)=.false.
        sn=(wet).and.(cshift(wets,-nx))
        sn(isb)=.false.
        en=(wet).and.(cshift(wets,1))
        wn=(wet).and.(cshift(wets,-1))
        if (zonally_periodic) then
            en(ieb)=(wet(ieb)).and.(wets(iwb))
            wn(iwb)=(wet(iwb)).and.(wets(ieb))
        else
            en(ieb)=.false.
            wn(iwb)=.false.
        endif

        ! if a point adjacent to ans boundary has multiple neighbours, just do one neutral
        ! calculation
        wn= (wn) .and. (.not.(en))
        nn= (nn) .and. (.not.(wn)) .and. (.not.(en))
        sn= (sn) .and. (.not.(nn)) .and. (.not.(wn)) .and. (.not.(en))

        ! expand in western direction
        inds=[(i,i=1,nxy)]
        inds_nbr=inds+1
        inds_nbr(ieb)=iwb

        call find(en,iwo)
        inbr=inds_nbr(iwo)

        allocate(sns_new(size(iwo)))
        allocate(ctns_new(size(iwo)))
        allocate(pns_new(size(iwo)))

        call depth_ntp_iter_drho(sns_(inbr),ctns_(inbr),pns_(inbr),  &
                            s_(iwo,:),ct_(iwo,:),p_(iwo,:), 0*sns_(inbr), &
                            sns_new,ctns_new,pns_new)
!        call depth_ntp_simple(sns_(inbr),ctns_(inbr),pns_(inbr),  &
!                            s_(iwo,:),ct_(iwo,:),p_(iwo,:), 0*sns_(inbr), &
!                            sns_new,ctns_new,pns_new)

        sns_(iwo)=sns_new
        ctns_(iwo)=ctns_new
        pns_(iwo)=pns_new

        nneighbours=count(.not.(isnan(sns_new)))

        deallocate(sns_new)
        deallocate(ctns_new)
        deallocate(pns_new)


        ! expand in eastern direction
        inds=[(i,i=1,nxy)]
        inds_nbr=inds-1
        inds_nbr(iwb)=ieb

        call find(wn,iwo)
        inbr=inds_nbr(iwo)

        allocate(sns_new(size(iwo)))
        allocate(ctns_new(size(iwo)))
        allocate(pns_new(size(iwo)))

        call depth_ntp_iter_drho(sns_(inbr),ctns_(inbr),pns_(inbr),  &
                            s_(iwo,:),ct_(iwo,:),p_(iwo,:), 0*sns_(inbr), &
                            sns_new,ctns_new,pns_new)
!        call depth_ntp_simple(sns_(inbr),ctns_(inbr),pns_(inbr),  &
!                            s_(iwo,:),ct_(iwo,:),p_(iwo,:), 0*sns_(inbr), &
!                            sns_new,ctns_new,pns_new)

        sns_(iwo)=sns_new
        ctns_(iwo)=ctns_new
        pns_(iwo)=pns_new

        nneighbours=nneighbours+count(.not.(isnan(sns_new)))

        deallocate(sns_new)
        deallocate(ctns_new)
        deallocate(pns_new)


        ! expand in southern direction
        inds=[(i,i=1,nxy)]
        inds_nbr=inds+nx
        inds_nbr(inb)=isb

        call find(nn,iwo)
        inbr=inds_nbr(iwo)

        allocate(sns_new(size(iwo)))
        allocate(ctns_new(size(iwo)))
        allocate(pns_new(size(iwo)))

        call depth_ntp_iter_drho(sns_(inbr),ctns_(inbr),pns_(inbr),  &
                            s_(iwo,:),ct_(iwo,:),p_(iwo,:), 0*sns_(inbr), &
                            sns_new,ctns_new,pns_new)
!       call depth_ntp_simple(sns_(inbr),ctns_(inbr),pns_(inbr),  &
!                           s_(iwo,:),ct_(iwo,:),p_(iwo,:), 0*sns_(inbr), &
!                           sns_new,ctns_new,pns_new)

        sns_(iwo)=sns_new
        ctns_(iwo)=ctns_new
        pns_(iwo)=pns_new

        nneighbours=nneighbours+count(.not.(isnan(sns_new)))

        deallocate(sns_new)
        deallocate(ctns_new)
        deallocate(pns_new)

        ! expand in northern direction
        inds=[(i,i=1,nxy)]
        inds_nbr=inds-nx
        inds_nbr(isb)=inb

        call find(sn,iwo)
        inbr=inds_nbr(iwo)

        allocate(sns_new(size(iwo)))
        allocate(ctns_new(size(iwo)))
        allocate(pns_new(size(iwo)))

        call depth_ntp_iter_drho(sns_(inbr),ctns_(inbr),pns_(inbr),  &
                            s_(iwo,:),ct_(iwo,:),p_(iwo,:), 0*sns_(inbr), &
                            sns_new,ctns_new,pns_new)
!       call depth_ntp_simple(sns_(inbr),ctns_(inbr),pns_(inbr),  &
!                           s_(iwo,:),ct_(iwo,:),p_(iwo,:), 0*sns_(inbr), &
!                           sns_new,ctns_new,pns_new)

        sns_(iwo)=sns_new
        ctns_(iwo)=ctns_new
        pns_(iwo)=pns_new

        nneighbours=nneighbours+count(.not.(isnan(sns_new)))

        sns=reshape(sns_, [nx,ny])
        ctns=reshape(ctns_, [nx,ny])
        pns=reshape(pns_, [nx,ny])

    end subroutine wetting


    subroutine wetting_simple(sns,ctns,pns,s,ct,p,nneighbours)
        real(rk), dimension(:,:), intent(inout) :: sns, ctns, pns
        real(rk), dimension(:,:,:), intent(in) :: s, ct, p
        integer, intent(out) :: nneighbours

        real(rk), dimension(:), allocatable :: sns_, ctns_, pns_, s1
        real(rk), dimension(:,:), allocatable :: s_, ct_, p_
        real(rk), dimension(:), allocatable :: sns_new, ctns_new, pns_new
        logical, dimension(:), allocatable :: wet, wets, nn, sn, en, wn, nbr
        integer, dimension(:), allocatable :: inds, inds_nbr
        integer, dimension(:), allocatable :: inb, isb
        integer, dimension(:), allocatable :: ieb, iwb
        integer, dimension(:), allocatable :: inbr, iwo
        integer :: i, nxy, nx, ny, nz
        logical :: zonally_periodic

        namelist /domain/ zonally_periodic

        open(1,file='user_input.nml')
        read(1,domain)
        close(1)

        nx=size(s,1)
        ny=size(s,2)
        nz=size(s,3)
        nxy=nx*ny

        allocate(sns_(nxy))
        allocate(ctns_(nxy))
        allocate(pns_(nxy))
        allocate(s1(nxy))

        allocate(s_(nxy,nz))
        allocate(ct_(nxy,nz))
        allocate(p_(nxy,nz))

        allocate(wet(nxy))
        allocate(wets(nxy))
        allocate(nn(nxy))
        allocate(sn(nxy))
        allocate(en(nxy))
        allocate(wn(nxy))
        allocate(nbr(nxy))

        allocate(inb(nx))
        allocate(isb(nx))
        allocate(ieb(ny))
        allocate(iwb(ny))

        sns_=pack(sns,.true.)
        ctns_=pack(ctns,.true.)
        pns_=pack(pns,.true.)

        s_=reshape(s,[nxy,nz])
        ct_=reshape(ct,[nxy,nz])
        p_=reshape(p,[nxy,nz])

        s1=pack(s(:,:,1),.true.)

        wet=(.not.(isnan(s1))).and.(isnan(sns_)) ! wet points at ocean surface excluding ans
        wets=(.not.isnan(sns_)) ! wet points on ans

        inb=[(i, i= (ny-1)*nx+1, ny*nx)] ! indices of northern boundary
        isb=[(i, i= 1, nx)]
        ieb=[(i, i= nx, nxy, nx)]
        iwb=[(i, i= 1, (ny-1)*nx+1, nx)]

        nn=(wet).and.(cshift(wets,nx)) ! wet points with northern neighbour on ans
        nn(inb)=.false.
        sn=(wet).and.(cshift(wets,-nx))
        sn(isb)=.false.
        en=(wet).and.(cshift(wets,1))
        wn=(wet).and.(cshift(wets,-1))
        if (zonally_periodic) then
            en(ieb)=(wet(ieb)).and.(wets(iwb))
            wn(iwb)=(wet(iwb)).and.(wets(ieb))
        else
            en(ieb)=.false.
            wn(iwb)=.false.
        endif

        ! if a point adjacent to ans boundary has multiple neighbours, just do one neutral
        ! calculation
        wn= (wn) .and. (.not.(en))
        nn= (nn) .and. (.not.(wn)) .and. (.not.(en))
        sn= (sn) .and. (.not.(nn)) .and. (.not.(wn)) .and. (.not.(en))

        ! expand in western direction
        inds=[(i,i=1,nxy)]
        inds_nbr=inds+1
        inds_nbr(ieb)=iwb

        call find(en,iwo)
        inbr=inds_nbr(iwo)

        allocate(sns_new(size(iwo)))
        allocate(ctns_new(size(iwo)))
        allocate(pns_new(size(iwo)))

!        call depth_ntp_iter_drho(sns_(inbr),ctns_(inbr),pns_(inbr),  &
!                            s_(iwo,:),ct_(iwo,:),p_(iwo,:), 0*sns_(inbr), &
!                            sns_new,ctns_new,pns_new)
        call depth_ntp_simple(sns_(inbr),ctns_(inbr),pns_(inbr),  &
                            s_(iwo,:),ct_(iwo,:),p_(iwo,:), 0*sns_(inbr), &
                            sns_new,ctns_new,pns_new)

        sns_(iwo)=sns_new
        ctns_(iwo)=ctns_new
        pns_(iwo)=pns_new

        nneighbours=count(.not.(isnan(sns_new)))

        deallocate(sns_new)
        deallocate(ctns_new)
        deallocate(pns_new)


        ! expand in eastern direction
        inds=[(i,i=1,nxy)]
        inds_nbr=inds-1
        inds_nbr(iwb)=ieb

        call find(wn,iwo)
        inbr=inds_nbr(iwo)

        allocate(sns_new(size(iwo)))
        allocate(ctns_new(size(iwo)))
        allocate(pns_new(size(iwo)))

!        call depth_ntp_iter_drho(sns_(inbr),ctns_(inbr),pns_(inbr),  &
!                            s_(iwo,:),ct_(iwo,:),p_(iwo,:), 0*sns_(inbr), &
!                            sns_new,ctns_new,pns_new)
        call depth_ntp_simple(sns_(inbr),ctns_(inbr),pns_(inbr),  &
                            s_(iwo,:),ct_(iwo,:),p_(iwo,:), 0*sns_(inbr), &
                            sns_new,ctns_new,pns_new)

        sns_(iwo)=sns_new
        ctns_(iwo)=ctns_new
        pns_(iwo)=pns_new

        nneighbours=nneighbours+count(.not.(isnan(sns_new)))

        deallocate(sns_new)
        deallocate(ctns_new)
        deallocate(pns_new)


        ! expand in southern direction
        inds=[(i,i=1,nxy)]
        inds_nbr=inds+nx
        inds_nbr(inb)=isb

        call find(nn,iwo)
        inbr=inds_nbr(iwo)

        allocate(sns_new(size(iwo)))
        allocate(ctns_new(size(iwo)))
        allocate(pns_new(size(iwo)))

!        call depth_ntp_iter_drho(sns_(inbr),ctns_(inbr),pns_(inbr),  &
!                            s_(iwo,:),ct_(iwo,:),p_(iwo,:), 0*sns_(inbr), &
!                            sns_new,ctns_new,pns_new)
        call depth_ntp_simple(sns_(inbr),ctns_(inbr),pns_(inbr),  &
                            s_(iwo,:),ct_(iwo,:),p_(iwo,:), 0*sns_(inbr), &
                            sns_new,ctns_new,pns_new)

        sns_(iwo)=sns_new
        ctns_(iwo)=ctns_new
        pns_(iwo)=pns_new

        nneighbours=nneighbours+count(.not.(isnan(sns_new)))

        deallocate(sns_new)
        deallocate(ctns_new)
        deallocate(pns_new)

        ! expand in northern direction
        inds=[(i,i=1,nxy)]
        inds_nbr=inds-nx
        inds_nbr(isb)=inb

        call find(sn,iwo)
        inbr=inds_nbr(iwo)

        allocate(sns_new(size(iwo)))
        allocate(ctns_new(size(iwo)))
        allocate(pns_new(size(iwo)))

!        call depth_ntp_iter_drho(sns_(inbr),ctns_(inbr),pns_(inbr),  &
!                            s_(iwo,:),ct_(iwo,:),p_(iwo,:), 0*sns_(inbr), &
!                            sns_new,ctns_new,pns_new)
        call depth_ntp_simple(sns_(inbr),ctns_(inbr),pns_(inbr),  &
                            s_(iwo,:),ct_(iwo,:),p_(iwo,:), 0*sns_(inbr), &
                            sns_new,ctns_new,pns_new)

        sns_(iwo)=sns_new
        ctns_(iwo)=ctns_new
        pns_(iwo)=pns_new

        nneighbours=nneighbours+count(.not.(isnan(sns_new)))

        sns=reshape(sns_, [nx,ny])
        ctns=reshape(ctns_, [nx,ny])
        pns=reshape(pns_, [nx,ny])

    end subroutine wetting_simple

    subroutine find(boolvec,itrue)
        logical, dimension(:), intent(in) :: boolvec
        integer, dimension(:), allocatable, intent(out) :: itrue

        integer :: i,j

        if (allocated(itrue)) deallocate(itrue)
        allocate(itrue(count(boolvec)))
        j=1
        do i=1,size(boolvec)
            if (boolvec(i)) then
                itrue(j)=i
                j=j+1
            endif
        enddo

    end subroutine find



    subroutine depth_ntp_iter(s0,ct0,p0,s,ct,p,sns,ctns,pns)

        real(rk), dimension(:), intent(in) :: s0, ct0, p0
        real(rk), dimension(:,:), intent(in) :: s, ct, p
        real(rk), dimension(:), intent(out) :: sns, ctns, pns

        real(rk), dimension(:), allocatable :: s0_, ct0_, p0_
!!!!!!!!!!!!!!!!!!!!!!!!
!        real(rk), dimension(:), allocatable :: mycast_flat, mybottle_flat
!        real(rk), dimension(:,:), allocatable :: pmid_, s0_stack, ct0_stack
!        real(rk), dimension(:,:), allocatable :: mycast, mybottle
!!!!!!!!!!!!!!!!!!!!!!!!
        real(rk), dimension(:,:), allocatable :: s_, ct_, p_
        real(rk), dimension(:,:), allocatable :: s0_stacked, ct0_stacked, p0_stacked
        integer :: stack, nxy, nz, refine_ints, cnt, i,j,k
        integer, dimension(:), allocatable :: inds
        logical, dimension(:), allocatable :: fr
        real(rk), dimension(:), allocatable :: s0_old, ct0_old, p0_old
        real(rk), dimension(:,:), allocatable :: F
        real(rk) :: cast_, bottle_

        call getnan(nan)

        nxy=size(s,1)
        nz=size(s,2)

        allocate(s0_(nxy))
        allocate(ct0_(nxy))
        allocate(p0_(nxy))
        s0_=s0
        ct0_=ct0
        p0_=p0

        allocate(s_(nxy,nz))
        allocate(ct_(nxy,nz))
        allocate(p_(nxy,nz))
        s_=s
        ct_=ct
        p_=p

        allocate(inds(nxy))
        inds=(/(i, i=1,nxy)/)
        allocate(fr(nxy))
        fr=.true.

        sns=nan
        ctns=nan
        pns=nan

        stack=nz
        refine_ints=2

        cnt=0
        do while (.true.)
            cnt=cnt+1

            allocate(F(count(fr),stack))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            allocate(pmid_(count(fr),stack))
!            allocate(s0_stack(count(fr),stack))
!            allocate(ct0_stack(count(fr),stack))
!            allocate(mycast(count(fr),stack))
!            allocate(mybottle(count(fr),stack))
!            allocate(mycast_flat(count(fr)*stack))
!            allocate(mybottle_flat(count(fr)*stack))
!
!            do k=1,stack
!               do i=1,size(F,1)
!                    pmid_(i,k)=0.5d0*(p_(i,k)+p0_(i))
!                    s0_stack(i,k)=s0_(i)
!                    ct0_stack(i,k)=ct0_(i)
!               enddo
!            enddo
!
!            call gsw_rho_vectorized(pack(s_,.true.),pack(ct_,.true.),pack(pmid_,.true.),mycast_flat)
!            call gsw_rho_vectorized(pack(s0_stack,.true.),pack(ct0_stack,.true.),pack(pmid_,.true.),mybottle_flat)
!
!            mycast=reshape(mycast_flat,[count(fr),stack])
!            mybottle=reshape(mybottle_flat,[count(fr),stack])
!
!            F=mycast-mybottle
!
!            deallocate(pmid_)
!            deallocate(s0_stack)
!            deallocate(ct0_stack)
!            deallocate(mycast)
!            deallocate(mybottle)
!            deallocate(mycast_flat)
!            deallocate(mybottle_flat)

            do k=1,stack
               do i=1,size(F,1)
                    cast_=gsw_rho(s_(i,k),ct_(i,k),0.5d0*(p_(i,k)+p0_(i)))
                    bottle_=gsw_rho(s0_(i),ct0_(i),0.5d0*(p_(i,k)+p0_(i)))

                    F(i,k)=cast_-bottle_
               enddo
            enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            deallocate(fr)
            call root_core(F,inds,refine_ints, &
                            s_,ct_,p_,sns,ctns,pns,fr)

            if (all(.not.(fr))) then
                exit
            endif

            stack=refine_ints+1

            deallocate(F)

            allocate(s0_old(size(fr)))
            allocate(ct0_old(size(fr)))
            allocate(p0_old(size(fr)))

            s0_old=s0_
            ct0_old=ct0_
            p0_old=p0_

            deallocate(s0_)
            deallocate(ct0_)
            deallocate(p0_)

            allocate(s0_(count(fr)))
            allocate(ct0_(count(fr)))
            allocate(p0_(count(fr)))

            j=1
            do i=1,size(fr)
                if (fr(i)) then
                    s0_(j)=s0_old(i)
                    ct0_(j)=ct0_old(i)
                    p0_(j)=p0_old(i)
                    j=j+1
                endif
            enddo

            deallocate(s0_old)
            deallocate(ct0_old)
            deallocate(p0_old)

        enddo

    end subroutine depth_ntp_iter


    subroutine depth_ntp_iter_drho(s0,ct0,p0,s,ct,p,drho,sns,ctns,pns)

        real(rk), dimension(:), intent(in) :: s0, ct0, p0
        real(rk), dimension(:,:), intent(in) :: s, ct, p
        real(rk), dimension(:), intent(in) :: drho
        real(rk), dimension(:), intent(out) :: sns, ctns, pns

        real(rk), dimension(:), allocatable :: s0_, ct0_, p0_, drho_
        real(rk), dimension(:,:), allocatable :: s_, ct_, p_
!!!!!!!!!!!!!!!!!!!!!!!!
!        real(rk), dimension(:), allocatable :: mycast_flat, mybottle_flat
!        real(rk), dimension(:,:), allocatable :: pmid_, s0_stack, ct0_stack, drho_stack
!        real(rk), dimension(:,:), allocatable :: mycast, mybottle
!!!!!!!!!!!!!!!!!!!!!!!!
        real(rk), dimension(:,:), allocatable :: s0_stacked, ct0_stacked, p0_stacked
        integer :: stack, nxy, nz, refine_ints, cnt, i,j,k
        integer, dimension(:), allocatable :: inds
        logical, dimension(:), allocatable :: fr
        real(rk), dimension(:), allocatable :: s0_old, ct0_old, p0_old, drho_old
        real(rk), dimension(:,:), allocatable :: F
        real(rk) :: cast_, bottle_

        call getnan(nan)

        nxy=size(s,1)
        nz=size(s,2)

        allocate(s0_(nxy))
        allocate(ct0_(nxy))
        allocate(p0_(nxy))
        allocate(drho_(nxy))
        s0_=s0
        ct0_=ct0
        p0_=p0
        drho_=drho

        allocate(s_(nxy,nz))
        allocate(ct_(nxy,nz))
        allocate(p_(nxy,nz))
        s_=s
        ct_=ct
        p_=p


        allocate(inds(nxy))
        inds=(/(i, i=1,nxy)/)
        allocate(fr(nxy))
        fr=.true.

        sns=nan
        ctns=nan
        pns=nan

!        ! discard land
!        nn=~isnan(s0);
!        iwet=~(sum(nn,1)==0);
!
!        s0=s0(iwet);
!        ct0=ct0(iwet);
!        p0=p0(iwet);
!        drho=drho(iwet);
!        s=s(:,iwet);
!        ct=ct(:,iwet);
!        p=p(:,iwet);
!        inds=inds(iwet);

        stack=nz
        refine_ints=2

        cnt=0
        do while (.true.)
            cnt=cnt+1

            allocate(F(count(fr),stack))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            allocate(pmid_(count(fr),stack))
!            allocate(s0_stack(count(fr),stack))
!            allocate(ct0_stack(count(fr),stack))
!            allocate(drho_stack(count(fr),stack))
!            allocate(mycast(count(fr),stack))
!            allocate(mybottle(count(fr),stack))
!            allocate(mycast_flat(count(fr)*stack))
!            allocate(mybottle_flat(count(fr)*stack))
!
!            do k=1,stack
!               do i=1,size(F,1)
!                    pmid_(i,k)=0.5d0*(p_(i,k)+p0_(i))
!                    s0_stack(i,k)=s0_(i)
!                    ct0_stack(i,k)=ct0_(i)
!                    drho_stack(i,k)=drho_(i)
!               enddo
!            enddo
!
!            call gsw_rho_vectorized(pack(s_,.true.),pack(ct_,.true.),pack(pmid_,.true.),mycast_flat)
!            call gsw_rho_vectorized(pack(s0_stack,.true.),pack(ct0_stack,.true.),pack(pmid_,.true.),mybottle_flat)
!
!            mycast=reshape(mycast_flat,[count(fr),stack])
!            mybottle=reshape(mybottle_flat,[count(fr),stack])
!
!            F=mycast-mybottle+drho_stack
!
!            deallocate(pmid_)
!            deallocate(s0_stack)
!            deallocate(ct0_stack)
!            deallocate(drho_stack)
!            deallocate(mycast)
!            deallocate(mybottle)
!            deallocate(mycast_flat)
!            deallocate(mybottle_flat)

            do k=1,stack
               do i=1,size(F,1)
                    cast_=gsw_rho(s_(i,k),ct_(i,k),0.5d0*(p_(i,k)+p0_(i)))
                    bottle_=gsw_rho(s0_(i),ct0_(i),0.5d0*(p_(i,k)+p0_(i)))

                    F(i,k)=cast_-bottle_+drho_(i)
               enddo
            enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            deallocate(fr)
            call root_core(F,inds,refine_ints, &
                            s_,ct_,p_,sns,ctns,pns,fr)

            if (all(.not.(fr))) then
                exit
            endif

            stack=refine_ints+1

            deallocate(F)

            allocate(s0_old(size(fr)))
            allocate(ct0_old(size(fr)))
            allocate(p0_old(size(fr)))
            allocate(drho_old(size(fr)))

            s0_old=s0_
            ct0_old=ct0_
            p0_old=p0_
            drho_old=drho_

            deallocate(s0_)
            deallocate(ct0_)
            deallocate(p0_)
            deallocate(drho_)

            allocate(s0_(count(fr)))
            allocate(ct0_(count(fr)))
            allocate(p0_(count(fr)))
            allocate(drho_(count(fr)))

            j=1
            do i=1,size(fr)
                if (fr(i)) then
                    s0_(j)=s0_old(i)
                    ct0_(j)=ct0_old(i)
                    p0_(j)=p0_old(i)
                    drho_(j)=drho_old(i)
                    j=j+1
                endif
            enddo

            deallocate(s0_old)
            deallocate(ct0_old)
            deallocate(p0_old)
            deallocate(drho_old)

        enddo

    end subroutine depth_ntp_iter_drho


    subroutine depth_ntp_simple(s0_,ct0_,p0_,s_,ct_,p_,drho_,sns,ctns,pns)

        real(rk), dimension(:), intent(in) :: s0_, ct0_, p0_
        real(rk), dimension(:,:), intent(in) :: s_, ct_, p_
        real(rk), dimension(:), intent(in) :: drho_
        real(rk), dimension(:), intent(out) :: sns, ctns, pns

        real(rk), dimension(:), allocatable :: s0, ct0, p0, drho
        real(rk), dimension(:,:), allocatable :: s, ct, p
        real(rk), dimension(:), allocatable :: F1, F2, s1, s2, ct1, ct2, p1, p2


        integer, dimension(:), allocatable :: inds, k1
        logical, dimension(:), allocatable :: wet, go_deep, go_shallow
        logical, dimension(:), allocatable :: search_initial, done
        logical, dimension(:), allocatable :: zoom, found_1, found_2, found
        real(rk) :: pmid, cast, bottle, F, ct_mid, s_mid, p_mid, F_mid, delta_root
        integer :: nxy, nz, nxy_wet, cnt, i, k

        real(rk), dimension(:), allocatable ::  debu
        integer, dimension(:), allocatable ::  debu_int

        namelist /root_finding/ delta_root

        open(1,file='user_input.nml')
        read(1,root_finding)
        close(1)

        call getnan(nan)
        sns=nan
        ctns=nan
        pns=nan

        nxy=size(s_,1)
        nz=size(s_,2)

        allocate(wet(nxy))

        ! discard land
        wet=.not.(isnan(s0_))
        nxy_wet=count(wet)

        allocate(s0(nxy_wet))
        allocate(ct0(nxy_wet))
        allocate(p0(nxy_wet))
        allocate(drho(nxy_wet))

        allocate(s(nxy_wet,nz))
        allocate(ct(nxy_wet,nz))
        allocate(p(nxy_wet,nz))

        s0=nan
        ct0=nan
        p0=nan
        drho=nan
        s=nan
        ct=nan
        p=nan

        allocate(inds(nxy_wet))

        cnt=1
        do i=1,nxy
            if (wet(i)) then
                s0(cnt)=s0_(i)
                ct0(cnt)=ct0_(i)
                p0(cnt)=p0_(i)
                drho(cnt)=drho_(i)
                inds(cnt)=i
                s(cnt,:)=s_(i,:)
                ct(cnt,:)=ct_(i,:)
                p(cnt,:)=p_(i,:)
                cnt=cnt+1
            endif
        enddo

        allocate(k1(nxy_wet))
        allocate(F1(nxy_wet))
        allocate(F2(nxy_wet))
        allocate(go_shallow(nxy_wet))
        allocate(go_deep(nxy_wet))
        allocate(search_initial(nxy_wet))
        allocate(zoom(nxy_wet))
        allocate(done(nxy_wet))
        allocate(found(nxy_wet))
        allocate(found_1(nxy_wet))
        allocate(found_2(nxy_wet))
        allocate(s1(nxy_wet))
        allocate(s2(nxy_wet))
        allocate(ct1(nxy_wet))
        allocate(ct2(nxy_wet))
        allocate(p1(nxy_wet))
        allocate(p2(nxy_wet))

        F1=nan
        F2=nan
        s1=nan
        s2=nan
        ct1=nan
        ct2=nan
        p1=nan
        p2=nan
        k1=-99 ! k1 of type integer
        found=.false.
        done=.false.

        do k=1,nz-1
            do i=1,nxy_wet
                if ((p0(i)>=p(i,k)).and.(p0(i)<p(i,k+1))) then
                    k1(i)=k
                    pmid=0.5d0*(p(i,k)+p0(i))
                    cast=gsw_rho(s(i,k),ct(i,k),pmid)
                    bottle=gsw_rho(s0(i),ct0(i),pmid)
                    F1(i)=cast-bottle+drho(i)

                    pmid=0.5d0*(p(i,k+1)+p0(i))
                    cast=gsw_rho(s(i,k+1),ct(i,k+1),pmid)
                    bottle=gsw_rho(s0(i),ct0(i),pmid)
                    F2(i)=cast-bottle+drho(i)
                endif
            enddo
        enddo

        go_shallow = F1 >= 0.0d0+delta_root ! go shallower to find first bottle pair
        go_deep    = F2 <= 0.0d0-delta_root ! go deeper to find first bottle pair

        search_initial = (go_deep).or.(go_shallow)

        zoom = (F1<0.0d0+delta_root).and.(F2>0.0d0-delta_root) ! bisect where true

        done = ( .not.(go_shallow.or.go_deep) ).and.(.not.(zoom)) ! nans

        do i=1,nxy_wet
            if (zoom(i)) then
                s1(i)=s(i,k1(i))
                s2(i)=s(i,k1(i)+1)
                ct1(i)=ct(i,k1(i))
                ct2(i)=ct(i,k1(i)+1)
                p1(i)=p(i,k1(i))
                p2(i)=p(i,k1(i)+1)
            endif
        enddo

!!!!!!!!!!!!!!!!!!
!        allocate(debu(nxy))
!        allocate(debu_int(nxy))
!        debu=nan
!        write(*,*) 'nsn: ', nan
!        do i=1,nxy_wet
!!            if (k1(i).ne.-99) then
!!                debu(inds(i))=dble(k1(i))
!!            else
!!                 debu(inds(i))=nan
!!            endif
!
!            if (zoom(i)) then
!                debu(inds(i))=1.0d0
!            else
!                debu(inds(i))=0.0d0
!            endif
!!
!!            debu(inds(i))=s1(i)
!!            write(*,*) 'k1(i): ', k1(i)
!        enddo
!        call ncwrite(pack(debu,.true.),(/90,43/),'debu.nc','debu')
!
!!!!!!!!!!!!!!!!!!!

        do while (any(.not.(done)))
            do i=1,nxy_wet
                if (.not.(done(i))) then
                    if (search_initial(i)) then

                        if (go_deep(i)) then
                            k1(i)=k1(i)+1 ! upper bottle moves one deeper
                            k=k1(i)+1 ! only need to evaluate the lower bottle
                            if (k.le.nz) then
                                pmid=0.5d0*(p(i,k)+p0(i))
                                cast=gsw_rho(s(i,k),ct(i,k),pmid)
                                bottle=gsw_rho(s0(i),ct0(i),pmid)
                                F=cast-bottle+drho(i)
                                if (.not.(isnan(F))) then
                                    F1(i)=F2(i)
                                    F2(i)=F
                                else
                                    done(i)=.true. ! bottom
                                endif
                            else
                                done(i)=.true. ! below model domain
                            endif
                            if (F>0.0d0-delta_root) then
                                search_initial(i)=.false.
                                zoom(i)=.true. ! start bisection
                            endif
                        elseif (go_shallow(i)) then
                            k1(i)=k1(i)-1 ! upper bottle moves one up
                            k=k1(i) ! only need to evaluate upper bottle
                            if (k.ge.1) then
                                pmid=0.5d0*(p(i,k)+p0(i))
                                cast=gsw_rho(s(i,k),ct(i,k),pmid)
                                bottle=gsw_rho(s0(i),ct0(i),pmid)
                                F=cast-bottle+drho(i)
                                if (.not.(isnan(F))) then
                                    F2(i)=F1(i)
                                    F1(i)=F
                                else
                                    done(i)=.true. ! in sea ice
                                endif
                            else
                                done(i)=.true. ! above free surface
                            endif
                            if (F<0.0d0+delta_root) then
                                search_initial(i)=.false.
                                zoom(i)=.true. ! start bisection
                            endif
                        endif
                        if (zoom(i)) then
                            s1(i)=s(i,k1(i))
                            s2(i)=s(i,k1(i)+1)
                            ct1(i)=ct(i,k1(i))
                            ct2(i)=ct(i,k1(i)+1)
                            p1(i)=p(i,k1(i))
                            p2(i)=p(i,k1(i)+1)
                        endif

                    elseif (zoom(i)) then
                        ! check if we are done
                        found_1(i) = (abs(F1(i))).lt.(delta_root)
                        found_2(i) = (abs(F2(i))).lt.(delta_root)
                        if (found_1(i)) then
                            sns(inds(i))=s1(i)
                            ctns(inds(i))=ct1(i)
                            pns(inds(i))=p1(i)
                        elseif (found_2(i)) then
                            sns(inds(i))=s2(i)
                            ctns(inds(i))=ct2(i)
                            pns(inds(i))=p2(i)
                        endif
                        done(i) = (found_1(i)).or.(found_2(i)) ! found zero crossing

                        if (.not.(done(i))) then ! bisect
                            s_mid=0.5d0*(s1(i)+s2(i))
                            ct_mid=0.5d0*(ct1(i)+ct2(i))
                            p_mid=0.5d0*(p1(i)+p2(i))
                            pmid=0.5d0*(p_mid+p0(i))
                            cast=gsw_rho(s_mid,ct_mid,pmid)
                            bottle=gsw_rho(s0(i),ct0(i),pmid)
                            F_mid=cast-bottle+drho(i)
                            if (F_mid>0) then
                                F2(i)=F_mid
                                s2(i)=s_mid
                                ct2(i)=ct_mid
                                p2(i)=p_mid
                            else
                                F1(i)=F_mid
                                s1(i)=s_mid
                                ct1(i)=ct_mid
                                p1(i)=p_mid
                            endif
                        endif
                    endif ! search initial or zoom
                endif

            enddo
        enddo


    end subroutine depth_ntp_simple


    subroutine root_core(F,inds,refine_ints, &
                            s,ct,p,sns,ctns,pns,fr)

        real(rk), allocatable, dimension(:,:), intent(in) :: F
        integer, intent(in) :: refine_ints
        integer, allocatable, dimension(:), intent(inout) :: inds
        logical, allocatable, dimension(:), intent(out) :: fr
        real(rk), allocatable, dimension(:,:), intent(inout) :: s,ct,p
        real(rk), dimension(:), intent(inout) :: sns, ctns,pns

        logical, allocatable, dimension(:,:) :: F_p, F_n, zc_F_stable
        logical, allocatable, dimension(:) :: myfinal, cond1, any_zc_F_stable
        real(rk), allocatable, dimension(:,:) :: s_,ct_,p_
        real(rk), allocatable, dimension(:) :: F_neg
        integer, allocatable, dimension(:,:) :: cs
        integer, allocatable, dimension(:) :: k_zc,k_, inds_tmp, inds_local

        real(rk) :: delta_root, d
        integer :: i, j, stack, nxy, nxy_new


        namelist /root_finding/ delta_root

        open(1,file='user_input.nml')
        read(1,root_finding)
        close(1)

        stack=size(F,2)
        nxy=size(F,1)

        allocate(F_p(nxy,stack))
        allocate(F_n(nxy,stack))

        F_p= F>=0
        F_n= F<0

        allocate(zc_F_stable(nxy,stack))
        zc_F_stable= (F_n).and.(cshift(F_p,1,2))
!        deallocate(F_p)
!        deallocate(F_n)
        zc_F_stable(:,stack)=.false.

        allocate(cs(nxy,stack))

        do j=1,stack
            do i=1,nxy
                cs(i,j)=count(zc_F_stable(i,1:j))+1;
            enddo
        enddo

        cs=merge(0,cs,cs/=1)

        allocate(k_zc(nxy))
        k_zc=sum(cs,2)+1;! vertical index of shallowest stable zero crossing

        allocate(any_zc_F_stable(nxy))
        any_zc_F_stable=any(zc_F_stable,2)

!        deallocate(zc_F_stable)

        allocate(F_neg(nxy))
        F_neg=nan
        do i=1,nxy
            if (any_zc_F_stable(i)) then
                F_neg(i)=F(i,k_zc(i))
            endif
        enddo

        allocate(myfinal(nxy))
        allocate(cond1(nxy))
        allocate(fr(nxy))
        myfinal=(abs(F_neg)<=delta_root)
        cond1=abs(F_neg)>delta_root;
        fr= (any_zc_F_stable).and.(cond1)
!        deallocate(cond1)
!        deallocate(F_neg)

        do i=1,nxy
            if (myfinal(i)) then
                sns(inds(i))=s(i,k_zc(i))
                ctns(inds(i))=ct(i,k_zc(i))
                pns(inds(i))=p(i,k_zc(i))
            endif
        enddo

        if (all(.not.(fr))) then
            return
        endif

!        deallocate(myfinal)

        nxy_new=count(fr)

        allocate(inds_tmp(nxy_new))
        allocate(k_(nxy_new))
        allocate(inds_local(nxy_new))
        j=1
        do i=1,nxy
            if (fr(i)) then
                inds_tmp(j)=inds(i)
                inds_local(j)=i
                k_(j)=k_zc(i)
                j=j+1
            endif
        enddo

        deallocate(inds)
        allocate(inds(nxy_new))
        inds=inds_tmp

        allocate(s_(nxy_new,refine_ints+1))
        allocate(ct_(nxy_new,refine_ints+1))
        allocate(p_(nxy_new,refine_ints+1))

        do i=1,nxy_new
            d=(s(inds_local(i),k_(i)+1)-s(inds_local(i),k_(i)))/refine_ints
            s_(i,:)=s(inds_local(i),k_(i))+d*[(i, i=0,refine_ints,1)]
            d=(ct(inds_local(i),k_(i)+1)-ct(inds_local(i),k_(i)))/refine_ints
            ct_(i,:)=ct(inds_local(i),k_(i))+d*[(i, i=0,refine_ints,1)]
            d=(p(inds_local(i),k_(i)+1)-p(inds_local(i),k_(i)))/refine_ints
            p_(i,:)=p(inds_local(i),k_(i))+d*[(i, i=0,refine_ints,1)]
        enddo

        deallocate(s)
        deallocate(ct)
        deallocate(p)
        allocate(s(nxy_new,refine_ints+1))
        allocate(ct(nxy_new,refine_ints+1))
        allocate(p(nxy_new,refine_ints+1))
        s=s_
        ct=ct_
        p=p_

    end subroutine root_core


    subroutine mld(s,ct,p,cut_off_choice)

        real(rk), dimension(:,:,:), intent(in) :: s, ct, p
        real(rk), dimension(:,:), intent(out) :: cut_off_choice

        real(rk), dimension(:,:), allocatable :: rho, surf_dens_, thresh
        real(rk), dimension(:), allocatable :: pflat

        real(rk), dimension(:), allocatable :: surf_dens, p_, cut_off_choice_
        logical, dimension(:,:), allocatable :: pos
        integer, dimension(:), allocatable :: ip
        integer, dimension(4) :: C3
        integer :: i, j, k, nxy, nx, ny, nz

        nx=size(s,1)
        ny=size(s,2)
        nz=size(s,3)
        nxy=nx*ny

        allocate(rho(nx*ny,nz))
        allocate(surf_dens_(nx*ny,nz))
        allocate(thresh(nx*ny,nz))
        allocate(pflat(nx*ny*nz))

        allocate(surf_dens(nx*ny))
        allocate(p_(nx*ny))
        allocate(cut_off_choice_(nx*ny))
        allocate(pos(nx*ny,nz))
        allocate(ip(nx*ny))

        do k=1,nz
            do j=1,ny
                do i=1,nx
                    rho(i+(j-1)*nx,k)=gsw_rho(s(i,j,k),ct(i,j,k),0d0*p(i,j,k))
                enddo
            enddo
        enddo

        surf_dens=rho(:,1)
        surf_dens_= spread(surf_dens,2,nz)
        thresh= surf_dens_+0.3d0
        pos=(thresh-rho)>0.0d0
        ip=count(pos,2)

        pflat=pack(p,.true.)
        do i=1,nxy
            p_(i)=pflat(i+(ip(i)-1)*nxy)
        enddo

        call getnan(nan)

        cut_off_choice_=nan
        cut_off_choice_=merge(p_,cut_off_choice_,ip>0)
        cut_off_choice=reshape(cut_off_choice_,[nx,ny])

    end subroutine mld


    subroutine delta_tilde_rho(sns,ctns,pns,drhox,drhoy)

        real(rk), dimension(:,:), intent(in) :: sns, ctns, pns
        real(rk), dimension(:,:), intent(out) :: drhox,drhoy
        real(rk), dimension(:,:), allocatable :: gradx_s, grady_s, gradx_ct, grady_ct
        real(rk), dimension(:,:), allocatable :: r1, r2, sns_, ctns_, pmid
        integer :: i, j, nxy, nx, ny
        logical :: zonally_periodic

        namelist /domain/ zonally_periodic

        open(1,file='user_input.nml')
        read(1,domain)
        close(1)

        nx=size(sns,1)
        ny=size(sns,2)

        allocate(gradx_s(nx,ny))
        allocate(grady_s(nx,ny))
        allocate(gradx_ct(nx,ny))
        allocate(grady_ct(nx,ny))
        allocate(r1(nx,ny))
        allocate(r2(nx,ny))
        allocate(sns_(nx,ny))
        allocate(ctns_(nx,ny))
        allocate(pmid(nx,ny))

        call getnan(nan)

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


    subroutine lsqr_Ay(region,drhox,drhoy)
        !type(region_type), dimension(:), allocatable, intent(in) :: regions
        integer, dimension(:), intent(in) :: region
        real(rk), dimension(:,:), intent(in) :: drhox,drhoy

        real(rk),  dimension(:), allocatable :: drhox_, drhoy_
        logical, dimension(:), allocatable :: reg, en, nn
        integer, dimension(:), allocatable :: sreg, sreg_en, sreg_nn
        integer, allocatable, dimension(:) :: j1_ew, j2_ew, j1_ns, j2_ns
        integer :: i, j, nxy, nx, ny
        logical :: zonally_periodic
        integer, allocatable, dimension(:) :: pts

        namelist /domain/ zonally_periodic

        open(1,file='user_input.nml')
        read(1,domain)
        close(1)

        call getnan(nan)

        nx=size(drhox,1)
        ny=size(drhox,2)
        nxy=nx*ny

        allocate(drhox_(nxy))
        allocate(drhoy_(nxy))
        allocate(reg(nxy))
        allocate(en(nxy))
        allocate(nn(nxy))
        allocate(sreg(nxy))
        allocate(sreg_en(nxy))
        allocate(sreg_nn(nxy))

        drhox_=pack(drhox,.true.)
        drhoy_=pack(drhoy,.true.)

        !allocate(region(size(regions(1)%points)))
        !region=regions(1)%points

        reg=.false.
        do i=1,size(region)
            reg(region(i))=.true.
        enddo

        en=reg.and.cshift(reg,1) ! en is true at a point if its eastward neighbor is in the region
        do i=nx,nx*ny,nx ! correct eastern boundary
            en(i)=reg(i).and.reg(i-nx+1)
        enddo

        if (.not.zonally_periodic) then
            en(nx:ny*nx:nx)=.false.
        endif

        nn=reg.and.cshift(reg,nx)
        nn(nx*(ny-1)+1:ny*nx)=.false.

        sreg(1)=merge(1,0,reg(1)) ! cumulative sum
        ! sparse indices of points forming the region (points of non-region are indexed with dummy)
        do i=2,size(reg)
            sreg(i)=sreg(i-1)+merge(1,0,reg(i))
        enddo
        sreg_en=cshift(sreg,1) ! sparse indices of eastward neighbours

        do i=nx,nx*ny,nx ! correct eastern boundary
            sreg_en(i)=sreg(i-nx+1)
        enddo

        allocate(j1_ew(count(en)))
        allocate(j2_ew(count(en)))
        allocate(y(count(en)+count(nn)+1)) ! A*x=y

        j=1
        do i=1,size(en)
            if (en(i)) then
                j1_ew(j)=sreg_en(i) ! j1 are j-indices for matrix coefficient 1
                j2_ew(j)=sreg(i) ! j2 are j-indices for matrix coefficient -1
                y(j)=drhox_(i)
                j=j+1
            endif
        enddo

        sreg_nn=cshift(sreg,nx)

        allocate(j1_ns(count(nn)))
        allocate(j2_ns(count(nn)))
        j=1
        do i=1,size(nn)
            if (nn(i)) then
                j1_ns(j)=sreg_nn(i) ! j1 are j-indices for matrix coefficient 1
                j2_ns(j)=sreg(i) ! j2 are j-indices for matrix coefficient -1
                y(count(en)+j)=drhoy_(i)
                j=j+1
            endif
        enddo

        if (allocated(j1)) then
            deallocate(j1)
        endif
        allocate(j1(size(j1_ew)+size(j1_ns)))
        j1(1:size(j1_ew))=j1_ew
        j1(size(j1_ew)+1: size(j1_ew)+size(j1_ns))=j1_ns

        if (allocated(j2)) then
            deallocate(j2)
        endif
        allocate(j2(size(j2_ew)+size(j2_ns)))
        j2(1:size(j2_ew))=j2_ew
        j2(size(j2_ew)+1: size(j2_ew)+size(j2_ns))=j2_ns

        !condition
        y(size(y))=0.0d0


    end subroutine lsqr_Ay


    subroutine solve_lsqr(regions,drhox,drhoy,drho)
        type(region_type), dimension(:), allocatable, intent(in) :: regions
        real(rk), dimension(:,:), intent(in) :: drhox,drhoy
        real(rk),  dimension(:,:), allocatable, intent(out) :: drho

        real(rk),  dimension(:), allocatable :: drho_
        integer, allocatable, dimension(:) :: region, j1, j2
        real(rk), allocatable, dimension(:) :: x
        integer :: i, j, nx, ny, ireg, nregions
        logical :: zonally_periodic
        namelist /domain/ zonally_periodic

        ! begin LSQR arguments, as described in lsqrModule.f90
        integer :: m, n
        real(rk) :: damp=0.0d0 ! damping parameter
        logical :: wantse=.false. ! standard error estimates
        real(rk), dimension(1) :: se=(/0.0d0/)
        !real(rk) :: atol=1.0d-4
        real(rk) :: atol=1.0d-15
        real(rk) :: btol=0.0d0
        real(rk) :: conlim=0.0d0
        !integer :: itnlim=450 ! max. number of iterations
        integer :: itnlim=50000 ! max. number of iterations
        integer :: nout=-99 ! output file
        integer :: istop=-99 ! reason for termination
        integer :: itn=-99 ! number of iterations performed
        real(rk) :: Anorm, Acond, rnorm, Arnorm, xnorm
        ! end LSQR arguments

        open(1,file='user_input.nml')
        read(1,domain)
        close(1)

        call getnan(nan)

        nx=size(drhox,1)
        ny=size(drhox,2)
        allocate(drho_(nx*ny))
        drho_(:)=nan

        nregions=size(regions)

        do ireg=1,nregions
            allocate(region(size(regions(ireg)%points)))
            region=regions(ireg)%points
            !write(*,*) 'size(regions): ', size(regions)

            call lsqr_Ay(region,drhox,drhoy)

            allocate(x(size(region)))
            x=0.0d0
            m=size(y)
            n=size(x)

            call LSQR  ( m, n, Aprod1, Aprod2, y, damp, wantse,         &
                         x, se,                                         &
                         atol, btol, conlim, itnlim, nout,              &
                         istop, itn, Anorm, Acond, rnorm, Arnorm, xnorm )
            deallocate(y)

            !write(*,*), 'Region Number'
            !write(*,*), 'Region Size'
            !write(*,*), 'Number of LSQR iterations'
            write(*,'(i3,1x,i5,1x,i5)') &
                ireg,size(region),itn

            do i=1,size(region)
                drho_(region(i))= x(i)
            enddo

            deallocate(x)
            deallocate(region)

        enddo

        allocate(drho(nx,ny))
        drho=reshape( drho_,(/nx,ny/) )

    end subroutine solve_lsqr


    subroutine Aprod1(m,n,x,y)

        integer,  intent(in)    :: m,n
        real(rk), intent(in)    :: x(n)
        real(rk), intent(inout) :: y(m)
        integer :: i

        do i=1,m-1
            y(i)=y(i)+x(j1(i))-x(j2(i))
        enddo
        ! condition
        y(m)=y(m)+sum(x)

    end subroutine Aprod1


    subroutine Aprod2(m,n,x,y)

        integer,  intent(in)    :: m,n
        real(rk), intent(inout)    :: x(n)
        real(rk), intent(in) :: y(m)
        integer :: i,j
        real(rk), dimension(m) :: tmp

        do i=1,m-1
            x(j1(i))=x(j1(i))+y(i)
            x(j2(i))=x(j2(i))-y(i)
        enddo

        ! condition
        do j=1,n
            x(j)=x(j)+y(m)
        enddo

    end subroutine Aprod2


!    subroutine grad_surf(f,e1t,e2t,fx,fy)
!
!        real(rk),  dimension(:,:), intent(in) :: f, e1t, e2t
!        real(rk),  dimension(:,:), intent(out) :: fx, fy
!        logical :: zonally_periodic
!
!        namelist /domain/ zonally_periodic
!
!        open(1,file='user_input.nml')
!        read(1,domain)
!        close(1)
!
!        fx=(cshift(f,1,1)-f)/e1t
!
!        call getnan(nan)
!
!        if (.not.(zonally_periodic)) then
!            fx(nx,:)=nan
!        endif
!
!        fy=(cshift(f,1,2)-f)/e2t
!
!    end subroutine grad_surf


    subroutine find_regions(pns,regions)

        type(region_type), dimension(:), allocatable, intent(out) :: regions
        real(rk),  dimension(:,:), intent(in) :: pns
        logical,  dimension(:,:), allocatable :: wet
        logical,  dimension(:), allocatable :: bool, wet_
        integer,  dimension(:), allocatable :: L_
        integer, dimension(:), allocatable :: east, west
        integer :: k, i, j, ii, iregion, nxy, nx, ny
        integer, allocatable, dimension(:) :: idx, neighbours, neighbours_tmp
        logical :: zonally_periodic

        namelist /domain/ zonally_periodic

        open(1,file='user_input.nml')
        read(1,domain)
        close(1)

        nx=size(pns,1)
        ny=size(pns,2)
        nxy=nx*ny

        allocate(wet(nx,ny))
        allocate(bool(nxy))
        allocate(wet_(nxy))
        allocate(L_(nxy))
        allocate(east(ny))
        allocate(west(ny))

        wet=.not.(isnan(pns))
        wet_=pack(wet,.true.)

        L_=0 ! label matrix
        iregion=1 ! region label

        do while (.true.)
            if (allocated(idx)) deallocate(idx)
            allocate(idx(1))
            ! find index of first wet point
            do i = 1, size(wet_)
                if (wet_(i)) then
                    ii = i
                    exit
                endif
            end do

            if (i==size(wet_)+1) exit

            idx(1)=ii ! linear indices of the pixels that were just labeled
            wet_(idx(1))=.false. ! set the pixel to 0
            L_(idx(1)) = iregion ! label first point

            do while (size(idx)/=0) ! find neighbours
                if (allocated(neighbours)) deallocate(neighbours)
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


    subroutine gsw_rho_vectorized(sa,ct,p,rho)

        real(rk),  dimension(:), intent(in) :: sa,ct,p
        real(rk),  dimension(:), allocatable, intent(out) :: rho
        real(rk),  dimension(:), allocatable :: sqrtsa, v_hat_denominator, v_hat_numerator

        real (rk), parameter :: v01 =  9.998420897506056d2, v02 =  2.839940833161907d0
        real (rk), parameter :: v03 = -3.147759265588511d-2, v04 =  1.181805545074306d-3
        real (rk), parameter :: v05 = -6.698001071123802d0, v06 = -2.986498947203215d-2
        real (rk), parameter :: v07 =  2.327859407479162d-4, v08 = -3.988822378968490d-2
        real (rk), parameter :: v09 =  5.095422573880500d-4, v10 = -1.426984671633621d-5
        real (rk), parameter :: v11 =  1.645039373682922d-7, v12 = -2.233269627352527d-2
        real (rk), parameter :: v13 = -3.436090079851880d-4, v14 =  3.726050720345733d-6
        real (rk), parameter :: v15 = -1.806789763745328d-4, v16 =  6.876837219536232d-7
        real (rk), parameter :: v17 = -3.087032500374211d-7, v18 = -1.988366587925593d-8
        real (rk), parameter :: v19 = -1.061519070296458d-11, v20 =  1.550932729220080d-10
        real (rk), parameter :: v21 =  1.0d0, v22 =  2.775927747785646d-3, v23 = -2.349607444135925d-5
        real (rk), parameter :: v24 =  1.119513357486743d-6, v25 =  6.743689325042773d-10
        real (rk), parameter :: v26 = -7.521448093615448d-3, v27 = -2.764306979894411d-5
        real (rk), parameter :: v28 =  1.262937315098546d-7, v29 =  9.527875081696435d-10
        real (rk), parameter :: v30 = -1.811147201949891d-11, v31 = -3.303308871386421d-5
        real (rk), parameter :: v32 =  3.801564588876298d-7, v33 = -7.672876869259043d-9
        real (rk), parameter :: v34 = -4.634182341116144d-11, v35 =  2.681097235569143d-12
        real (rk), parameter :: v36 =  5.419326551148740d-6, v37 = -2.742185394906099d-5
        real (rk), parameter :: v38 = -3.212746477974189d-7, v39 =  3.191413910561627d-9
        real (rk), parameter :: v40 = -1.931012931541776d-12, v41 = -1.105097577149576d-7
        real (rk), parameter :: v42 =  6.211426728363857d-10, v43 = -1.119011592875110d-10
        real (rk), parameter :: v44 = -1.941660213148725d-11, v45 = -1.864826425365600d-14
        real (rk), parameter :: v46 =  1.119522344879478d-14, v47 = -1.200507748551599d-15
        real (rk), parameter :: v48 =  6.057902487546866d-17


        allocate(rho(size(sa)))
        allocate(sqrtsa(size(sa)))
        allocate(v_hat_denominator(size(sa)))
        allocate(v_hat_numerator(size(sa)))

        sqrtsa = sqrt(sa)

        v_hat_denominator = v01 + ct*(v02 + ct*(v03 + v04*ct))  &
                     + sa*(v05 + ct*(v06 + v07*ct) &
                 + sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct)))) &
                      + p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct) &
                      + p*(v17 + ct*(v18 + v19*ct) + v20*sa))

        v_hat_numerator = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct))) &
                   + sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa &
               + sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))  &
                    + p*(v37 + ct*(v38 + ct*(v39 + v40*ct))  &
                   + sa*(v41 + v42*ct) &
                    + p*(v43 + ct*(v44 + v45*ct + v46*sa) &
                    + p*(v47 + v48*ct)))

        rho = v_hat_denominator/v_hat_numerator


    end subroutine gsw_rho_vectorized

end module ansu
