!This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
!Partially modified by P. Barker (2010-13)
!Partially modified by S. Riha (2013)
!Principal investigator: Trevor McDougall
!
!Translated to Fortran by S. Riha (2013)

module definitions
    integer, parameter, public :: rk = selected_real_kind(14,30)

    type, public :: region_type
        integer, dimension(:), allocatable :: points
    end type region_type

    public :: getnan
    public :: nan
    real(rk) :: nan

contains
    subroutine getnan(nan)
        real(rk), intent(out) :: nan
        integer :: setnan
        setnan=0 ! hide division by 0 at compile time
        nan=0d0/setnan
    end subroutine getnan

end module definitions
