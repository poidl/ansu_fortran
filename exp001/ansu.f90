!This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
!Partially modified by P. Barker (2010-13)
!Partially modified by S. Riha (2013)
!Principal investigator: Trevor McDougall
!
!Translated to Fortran by S. Riha (2013)

module ansu
    use definitions
    use grid_params
    use lsqrModule, only : LSQR
!    use ncutils
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

        real(rk), dimension(:,:), intent(inout) :: sns, ctns, pns
        real(rk), dimension(:,:,:), intent(in) :: s, ct, p
        type(region_type), dimension(:), allocatable :: regions
        real(rk), dimension(nx,ny) :: cut_off_choice, drhox, drhoy, drho
        integer :: nneighbours

        !call mld(s,ct,p,cut_off_choice) ! disregard data above mixed layer depth
        call wetting(sns,ctns,pns,s,ct,p,nneighbours) ! lateral extension of outcropping/undercropping surface
        call delta_tilde_rho(sns,ctns,pns,drhox,drhoy) ! compute lateral gradient
        call find_regions(pns,regions); ! find regions
        call solve_lsqr(regions,drhox,drhoy,drho) ! lateral integration
        call dz_from_drho(sns, ctns, pns, s, ct, p, drho) ! vertical integration

    end subroutine optimize_surface


    subroutine dz_from_drho(sns, ctns, pns, s, ct, p, drho)
        real(rk), dimension(:,:), intent(inout) :: sns, ctns, pns
        real(rk), dimension(:,:,:), intent(in) :: s, ct, p
        real(rk), dimension(:,:), intent(in) :: drho
        real(rk), dimension(nx*ny) :: sns_out, ctns_out, pns_out
        real(rk), dimension(:,:), allocatable :: s_, ct_, p_
        real(rk), dimension(nx,ny) :: rho_surf, t2
        real(rk), dimension(:), allocatable :: pns_, t2_, pns_old,t2_old
        real(rk), allocatable, dimension(:,:) ::  F
        integer, dimension(:), allocatable :: inds
        logical, dimension(:), allocatable :: fr

        integer :: i,j,k, refine_ints, cnt, stack

        call getnan(nan)

        sns_out=nan;
        ctns_out=nan;
        pns_out=nan;

        allocate(pns_(nx*ny))
        pns_=pack(pns,.true.)

        allocate(inds(nx*ny))
        inds=(/(i, i=1,nx*ny)/)
        allocate(fr(nx*ny))
        fr=.true.

        refine_ints=100

        allocate(s_(nx*ny,nz))
        allocate(ct_(nx*ny,nz))
        allocate(p_(nx*ny,nz))
        s_=reshape(s,[nx*ny,nz])
        ct_=reshape(ct,[nx*ny,nz])
        p_=reshape(p,[nx*ny,nz])

        stack=nz

        rho_surf=nan
        do i=1,nx
            do j=1,ny
                rho_surf(i,j)=gsw_rho(sns(i,j),ctns(i,j),pns(i,j))
            enddo
        enddo

        t2=rho_surf-drho
        allocate(t2_(nx*ny))
        t2_=pack(t2,.true.)

        cnt=0
        do while (.true.)
            cnt=cnt+1

            allocate(F(count(fr),stack))

            do k=1,stack
                do i=1,size(F,1)
                    F(i,k)=gsw_rho(s_(i,k),ct_(i,k),pns_(i)) &
                                    -t2_(i)
                enddo
            enddo

            deallocate(fr)
            call root_core(F,inds,refine_ints, &
                            s_,ct_,p_,sns_out,ctns_out,pns_out,fr)

            if (all(.not.(fr))) then
                exit
            endif

            stack=refine_ints+1

            allocate(pns_old(size(fr)))
            allocate(t2_old(size(fr)))

            pns_old=pns_
            t2_old=t2_

            deallocate(pns_)
            deallocate(t2_)
            allocate(pns_(count(fr)))
            allocate(t2_(count(fr)))

            j=1
            do i=1,size(fr)
                if (fr(i)) then
                    pns_(j)=pns_old(i)
                    t2_(j)=t2_old(i)
                    j=j+1
                endif
            enddo

            deallocate(pns_old)
            deallocate(t2_old)

            deallocate(F)

        enddo

        sns=reshape(sns_out,[nx,ny])
        ctns=reshape(ctns_out,[nx,ny])
        pns=reshape(pns_out,[nx,ny])

    end subroutine dz_from_drho


    subroutine wetting(sns,ctns,pns,s,ct,p,nneighbours)
        real(rk), dimension(:,:), intent(inout) :: sns, ctns, pns
        real(rk), dimension(:,:,:), intent(in) :: s, ct, p
        integer, intent(out) :: nneighbours

        real(rk), dimension(nx*ny) :: sns_, ctns_, pns_, s1
        real(rk), dimension(nx*ny,nz) :: s_, ct_, p_
        real(rk), dimension(:), allocatable :: sns_new, ctns_new, pns_new
        logical, dimension(nx*ny) :: wet, wets, nn, sn, en, wn, nbr
        integer, dimension(nx*ny) :: inds, inds_nbr
        integer, dimension(nx) :: inb, isb
        integer, dimension(ny) :: ieb, iwb
        integer :: i
        integer, dimension(:), allocatable :: inbr, iwo
        logical :: zonally_periodic
        !!! debug
        real(rk), dimension(:), allocatable :: snstmp,ctnstmp,pnstmp
        real(rk), dimension(:,:), allocatable :: stmp,cttmp,ptmp
        real(rk), dimension(nx,ny) :: ma
        real(rk), dimension(4) :: test=[1,2,3,4]
        real(rk), dimension(2) :: jj=[0,0]
        !!! debug

        namelist /domain/ zonally_periodic

        open(1,file='user_input.nml')
        read(1,domain)
        close(1)

        sns_=pack(sns,.true.)
        ctns_=pack(ctns,.true.)
        pns_=pack(pns,.true.)

        s_=reshape(s,[nx*ny,nz])
        ct_=reshape(ct,[nx*ny,nz])
        p_=reshape(p,[nx*ny,nz])

        s1=pack(s(:,:,1),.true.)

        wet=(.not.(isnan(s1))).and.(isnan(sns_)) ! wet points at ocean surface excluding ans
        wets=(.not.isnan(sns_)) ! wet points on ans

        inb=[(i, i= (ny-1)*nx+1, ny*nx)] ! indices of northern boundary
        isb=[(i, i= 1, nx)]
        ieb=[(i, i= nx, nx*ny, nx)]
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
        inds=[(i,i=1,nx*ny)]
        inds_nbr=inds+1
        inds_nbr(ieb)=iwb

        call find(en,iwo)
        inbr=inds_nbr(iwo)

        allocate(sns_new(size(iwo)))
        allocate(ctns_new(size(iwo)))
        allocate(pns_new(size(iwo)))

        call depth_ntp_iter(sns_(inbr),ctns_(inbr),pns_(inbr),  &
                            s_(iwo,:),ct_(iwo,:),p_(iwo,:),  &
                            sns_new,ctns_new,pns_new)

        sns_(iwo)=sns_new
        ctns_(iwo)=ctns_new
        pns_(iwo)=pns_new

        nneighbours=count(.not.(isnan(sns_new)))

        deallocate(sns_new)
        deallocate(ctns_new)
        deallocate(pns_new)


        ! expand in eastern direction
        inds=[(i,i=1,nx*ny)]
        inds_nbr=inds-1
        inds_nbr(iwb)=ieb

        call find(wn,iwo)
        inbr=inds_nbr(iwo)

        allocate(sns_new(size(iwo)))
        allocate(ctns_new(size(iwo)))
        allocate(pns_new(size(iwo)))

        call depth_ntp_iter(sns_(inbr),ctns_(inbr),pns_(inbr),  &
                            s_(iwo,:),ct_(iwo,:),p_(iwo,:),  &
                            sns_new,ctns_new,pns_new)

        sns_(iwo)=sns_new
        ctns_(iwo)=ctns_new
        pns_(iwo)=pns_new

        nneighbours=nneighbours+count(.not.(isnan(sns_new)))

        deallocate(sns_new)
        deallocate(ctns_new)
        deallocate(pns_new)


        ! expand in southern direction
        inds=[(i,i=1,nx*ny)]
        inds_nbr=inds+nx
        inds_nbr(inb)=isb

        call find(nn,iwo)
        inbr=inds_nbr(iwo)

        allocate(sns_new(size(iwo)))
        allocate(ctns_new(size(iwo)))
        allocate(pns_new(size(iwo)))

        call depth_ntp_iter(sns_(inbr),ctns_(inbr),pns_(inbr),  &
                            s_(iwo,:),ct_(iwo,:),p_(iwo,:),  &
                            sns_new,ctns_new,pns_new)

        sns_(iwo)=sns_new
        ctns_(iwo)=ctns_new
        pns_(iwo)=pns_new

        nneighbours=nneighbours+count(.not.(isnan(sns_new)))

        deallocate(sns_new)
        deallocate(ctns_new)
        deallocate(pns_new)

        ! expand in northern direction
        inds=[(i,i=1,nx*ny)]
        inds_nbr=inds-nx
        inds_nbr(isb)=inb

        call find(sn,iwo)
        inbr=inds_nbr(iwo)

        allocate(sns_new(size(iwo)))
        allocate(ctns_new(size(iwo)))
        allocate(pns_new(size(iwo)))

        call depth_ntp_iter(sns_(inbr),ctns_(inbr),pns_(inbr),  &
                            s_(iwo,:),ct_(iwo,:),p_(iwo,:),  &
                            sns_new,ctns_new,pns_new)

        sns_(iwo)=sns_new
        ctns_(iwo)=ctns_new
        pns_(iwo)=pns_new

        nneighbours=nneighbours+count(.not.(isnan(sns_new)))

        sns=reshape(sns_, [nx,ny])
        ctns=reshape(ctns_, [nx,ny])
        pns=reshape(pns_, [nx,ny])

    end subroutine wetting

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
        real(rk), dimension(:,:), allocatable :: s_, ct_, p_
        real(rk), dimension(:,:), allocatable :: s0_stacked, ct0_stacked, p0_stacked
        integer :: stack, nxy, refine_ints, cnt, i,j,k
        integer, dimension(:), allocatable :: inds
        logical, dimension(:), allocatable :: fr
        real(rk), dimension(:), allocatable :: s0_old, ct0_old, p0_old
        real(rk), dimension(:,:), allocatable :: F,cast,bottle
        real(rk) :: cast_, bottle_

        call getnan(nan)

        nxy=size(s,1)

        allocate(s0_(nxy))
        allocate(ct0_(nxy))
        allocate(p0_(nxy))
        s0_=s0
        ct0_=ct0
        p0_=p0
!         call ncwrite(s0_,'s0_.nc','fort',1)
        allocate(s_(nxy,nz))
        allocate(ct_(nxy,nz))
        allocate(p_(nxy,nz))
        s_=s
        ct_=ct
        p_=p
!        call ncwrite_new(s_,'s_.nc','fort')

        allocate(inds(nxy))
        inds=(/(i, i=1,nxy)/)
        allocate(fr(nxy))
        fr=.true.

        sns=nan
        ctns=nan
        pns=nan

        stack=nz
        refine_ints=100

        cnt=0
        do while (.true.)
            cnt=cnt+1

            allocate(cast(count(fr),stack))
            allocate(bottle(count(fr),stack))
            allocate(F(count(fr),stack))

            do k=1,stack
               do i=1,size(F,1)
                    cast_=gsw_rho(s_(i,k),ct_(i,k),0.5d0*(p_(i,k)+p0_(i)))
                    bottle_=gsw_rho(s0_(i),ct0_(i),0.5d0*(p_(i,k)+p0_(i)))

                    F(i,k)=cast_-bottle_
               enddo
            enddo

            deallocate(fr)
            call root_core(F,inds,refine_ints, &
                            s_,ct_,p_,sns,ctns,pns,fr)

            if (all(.not.(fr))) then
                exit
            endif

            stack=refine_ints+1

            deallocate(cast)
            deallocate(bottle)
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


!    subroutine depth_ntp_iter_drho(s0,ct0,p0,s,ct,p,sns,ctns,pns)
!
!        real(rk), dimension(:), intent(in) :: s0, ct0, p0
!        real(rk), dimension(:,:), intent(in) :: s, ct, p
!        real(rk), dimension(:), intent(out) :: sns, ctns, pns
!
!        real(rk), dimension(:), allocatable :: s0_, ct0_, p0_
!        real(rk), dimension(:,:), allocatable :: s_, ct_, p_
!        real(rk), dimension(:,:), allocatable :: s0_stacked, ct0_stacked, p0_stacked
!        integer :: stack, nxy, refine_ints, cnt, i,j,k
!        integer, dimension(:), allocatable :: inds
!        logical, dimension(:), allocatable :: fr
!        real(rk), dimension(:), allocatable :: s0_old, ct0_old, p0_old
!        real(rk), dimension(:,:), allocatable :: F,cast,bottle
!        real(rk) :: cast_, bottle_
!
!        call getnan(nan)
!
!        nxy=size(s,1)
!
!        allocate(s0_(nxy))
!        allocate(ct0_(nxy))
!        allocate(p0_(nxy))
!        s0_=s0
!        ct0_=ct0
!        p0_=p0
!!         call ncwrite(s0_,'s0_.nc','fort',1)
!        allocate(s_(nxy,nz))
!        allocate(ct_(nxy,nz))
!        allocate(p_(nxy,nz))
!        s_=s
!        ct_=ct
!        p_=p
!!        call ncwrite_new(s_,'s_.nc','fort')
!
!        allocate(inds(nxy))
!        inds=(/(i, i=1,nxy)/)
!        allocate(fr(nxy))
!        fr=.true.
!
!        sns=nan
!        ctns=nan
!        pns=nan
!
!        stack=nz
!        refine_ints=100
!
!        cnt=0
!        do while (.true.)
!            cnt=cnt+1
!
!            allocate(cast(count(fr),stack))
!            allocate(bottle(count(fr),stack))
!            allocate(F(count(fr),stack))
!
!            do k=1,stack
!               do i=1,size(F,1)
!                    cast_=gsw_rho(s_(i,k),ct_(i,k),0.5d0*(p_(i,k)+p0_(i)))
!                    bottle_=gsw_rho(s0_(i),ct0_(i),0.5d0*(p_(i,k)+p0_(i)))
!
!                    F(i,k)=cast_-bottle_
!               enddo
!            enddo
!
!            deallocate(fr)
!            call root_core(F,inds,refine_ints, &
!                            s_,ct_,p_,sns,ctns,pns,fr)
!
!            if (all(.not.(fr))) then
!                exit
!            endif
!
!            stack=refine_ints+1
!
!            deallocate(cast)
!            deallocate(bottle)
!            deallocate(F)
!
!            allocate(s0_old(size(fr)))
!            allocate(ct0_old(size(fr)))
!            allocate(p0_old(size(fr)))
!
!            s0_old=s0_
!            ct0_old=ct0_
!            p0_old=p0_
!
!            deallocate(s0_)
!            deallocate(ct0_)
!            deallocate(p0_)
!
!            allocate(s0_(count(fr)))
!            allocate(ct0_(count(fr)))
!            allocate(p0_(count(fr)))
!
!            j=1
!            do i=1,size(fr)
!                if (fr(i)) then
!                    s0_(j)=s0_old(i)
!                    ct0_(j)=ct0_old(i)
!                    p0_(j)=p0_old(i)
!                    j=j+1
!                endif
!            enddo
!
!            deallocate(s0_old)
!            deallocate(ct0_old)
!            deallocate(p0_old)
!
!        enddo
!
!    end subroutine depth_ntp_iter_drho


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
        deallocate(F_p)
        deallocate(F_n)
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

        deallocate(zc_F_stable)

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
        deallocate(cond1)
        deallocate(F_neg)

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

        deallocate(myfinal)

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

        real(rk), dimension(nx*ny,nz) :: rho, surf_dens_, thresh
        real(rk), dimension(nx*ny*nz) :: pflat
        real(rk), dimension(nx*ny) :: surf_dens, p_, cut_off_choice_
        logical, dimension(nx*ny,nz) :: pos
        integer, dimension(nx*ny) :: ip
        integer, dimension(4) :: C3
        integer :: nxy, i, j, k

        nxy=nx*ny

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
        real(rk), dimension(nx, ny) :: gradx_s, grady_s, gradx_ct, grady_ct
        real(rk), dimension(nx,ny) :: r1, r2, sns_, ctns_, pmid
        real(rk), dimension(nx,ny) :: debug
        integer :: i, j, nxy
        logical :: zonally_periodic

        namelist /domain/ zonally_periodic

        open(1,file='user_input.nml')
        read(1,domain)
        close(1)

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


    subroutine lsqr_Ay(regions,drhox,drhoy)
        type(region_type), dimension(:), allocatable, intent(in) :: regions
        real(rk), dimension(:,:), intent(in) :: drhox,drhoy
        real(rk),  dimension(nx*ny) :: drhox_, drhoy_
        logical, dimension(nx*ny) :: reg, en, nn
        integer, dimension(nx*ny) :: sreg, sreg_en, sreg_nn
        integer, allocatable, dimension(:) :: region, j1_ew, j2_ew, j1_ns, j2_ns
        integer :: i, j
        logical :: zonally_periodic
        integer, allocatable, dimension(:) :: pts

        namelist /domain/ zonally_periodic

        open(1,file='user_input.nml')
        read(1,domain)
        close(1)

        call getnan(nan)

        drhox_=pack(drhox,.true.)
        drhoy_=pack(drhoy,.true.)

        allocate(region(size(regions(1)%points)))
        region=regions(1)%points

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

        allocate(j1(size(j1_ew)+size(j1_ns)))
        j1(1:size(j1_ew))=j1_ew
        j1(size(j1_ew)+1: size(j1_ew)+size(j1_ns))=j1_ns

        allocate(j2(size(j2_ew)+size(j2_ns)))
        j2(1:size(j2_ew))=j2_ew
        j2(size(j2_ew)+1: size(j2_ew)+size(j2_ns))=j2_ns

        !condition
        y(size(y))=0.0d0

!        call ncwrite(pack(dble(j1),.true.),'j1.nc','j1',1)
!        call ncwrite(pack(dble(j2),.true.),'j2.nc','j2',1)
!        call ncwrite(pack(y,.true.),'y.nc','y',1)

    end subroutine lsqr_Ay


    subroutine solve_lsqr(regions,drhox,drhoy,drho)
        type(region_type), dimension(:), allocatable, intent(in) :: regions
        real(rk), dimension(:,:), intent(in) :: drhox,drhoy
        real(rk),  dimension(nx,ny), intent(out) :: drho
        real(rk),  dimension(nx*ny) :: drho_
        integer, allocatable, dimension(:) :: j1, j2
        real(rk), allocatable, dimension(:) :: x,b
        integer :: i, j
        logical :: zonally_periodic
        namelist /domain/ zonally_periodic

        ! begin LSQR arguments, as described in lsqrModule.f90
        integer :: m, n
        real(rk) :: damp=0.0d0 ! damping parameter
        logical :: wantse=.false. ! standard error estimates
        real(rk), dimension(1) :: se=(/0.0d0/)
        real(rk) :: atol=1.0d-11
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

        call lsqr_Ay(regions,drhox,drhoy)

        allocate(x(size(regions(1)%points)))
        allocate(b(size(y)))
        x=0.0d0
        b=y
        m=size(y)
        n=size(x)

        write(*,*) 'calling LSQR...'
        call LSQR  ( m, n, Aprod1, Aprod2, b, damp, wantse,         &
                     x, se,                                         &
                     atol, btol, conlim, itnlim, nout,              &
                     istop, itn, Anorm, Acond, rnorm, Arnorm, xnorm )
        write(*,*) 'istop: ',istop
        write(*,*) 'rnorm: ',rnorm
        write(*,*) 'Arnorm/(Anorm*rnorm): ', Arnorm/(Anorm*rnorm)

        drho_(:)=nan
!        call ncwrite(pack(drho_,.true.),'drho_.nc','drho_',2)
        do i=1,size(regions(1)%points)
            drho_(regions(1)%points(i))= x(i)
        enddo
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


    subroutine grad_surf(f,e1t,e2t,fx,fy)

        real(rk),  dimension(:,:), intent(in) :: f, e1t, e2t
        real(rk),  dimension(:,:), intent(out) :: fx, fy
        logical :: zonally_periodic

        namelist /domain/ zonally_periodic

        open(1,file='user_input.nml')
        read(1,domain)
        close(1)

        fx=(cshift(f,1,1)-f)/e1t

        call getnan(nan)

        if (.not.(zonally_periodic)) then
            fx(nx,:)=nan
        endif

        fy=(cshift(f,1,2)-f)/e2t

    end subroutine grad_surf


    subroutine find_regions(pns,regions)

        type(region_type), dimension(:), allocatable, intent(out) :: regions
        real(rk),  dimension(:,:), intent(in) :: pns
        logical,  dimension(nx,ny) :: wet
        logical,  dimension(nx*ny) :: bool, wet_
        integer,  dimension(nx*ny) :: L_
        integer, dimension(ny) :: east, west
        integer :: k, i, j, ii, iregion,kk,kkk
        integer, allocatable, dimension(:) :: idx, neighbours, neighbours_tmp
        logical :: zonally_periodic

        namelist /domain/ zonally_periodic

        open(1,file='user_input.nml')
        read(1,domain)
        close(1)

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


!    subroutine epsilon_(sns,ctns,pns,e1t,e2t,ex,ey)
!        real(rk), dimension(:,:), intent(in) :: sns, ctns, pns, e1t, e2t
!        real(rk), dimension(:,:), intent(out) :: ex, ey
!        real(rk), dimension(nx, ny) :: gradx_s, grady_s, gradx_ct, grady_ct
!        real(rk), dimension(nx,ny) :: alpha, beta, alphax, alphay, betax, betay
!        real(rk), dimension(nx,ny) :: debug
!        integer :: i, j, nxy, setnan
!        real(rk) :: nan
!        logical :: zonally_periodic
!
!        namelist /user_input/ zonally_periodic
!
!        open(1,file='user_input.nml')
!        read(1,user_input)
!        close(1)
!
!        call ncwrite(pack(sns,.true.),'e1t.nc','e1t',2)
!        call grad_surf(ctns,e1t,e2t,gradx_ct,grady_ct)
!        call grad_surf(sns,e1t,e2t,gradx_s,grady_s)
!        call ncwrite(pack(ctns,.true.),'e1t.nc','e1t',2)
!        call ncwrite(pack(grady_ct,.true.),'grady_ct.nc','grady_ct',2)
!        call ncwrite(pack(gradx_s,.true.),'gradx_s.nc','gradx_s',2)
!
!        nxy=nx*ny
!
!        do j=1,ny
!            do i=1,nx
!                alpha(i,j)=gsw_alpha(sns(i,j),ctns(i,j),pns(i,j))
!            enddo
!        enddo
!        do j=1,ny
!            do i=1,nx
!                beta(i,j)=gsw_beta(sns(i,j),ctns(i,j),pns(i,j))
!            enddo
!        enddo
!
!        alphax=0.5d0*(alpha+cshift(alpha,1,1))
!        alphay=0.5d0*(alpha+cshift(alpha,1,2))
!        betax=0.5d0*(beta+cshift(beta,1,1))
!        betay=0.5d0*(beta+cshift(beta,1,2))
!
!        setnan=0 ! hide division by 0 at compile time
!        nan=0d0/setnan
!
!        alphay(:,ny)=nan
!        betay(:,ny)=nan
!
!        if (.not.(zonally_periodic)) then
!            alphax(nx,:)=nan
!            betax(nx,:)=nan
!        endif
!
!        ex=betax*gradx_s-alphax*gradx_ct
!        ey=betay*grady_s-alphay*grady_ct
!
!
!    end subroutine epsilon_

end module ansu
