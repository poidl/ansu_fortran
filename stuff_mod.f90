!This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
!Partially modified by P. Barker (2010-13)
!Partially modified by S. Riha (2013)
!Principal investigator: Trevor McDougall
!
!Translated to Fortran by S. Riha (2013)

module stuff_mod
    integer, parameter, public :: rk = selected_real_kind(14,30)
    !integer, parameter, public :: NX = 2, NY = 2, NZ=101
    integer, parameter, public :: NX = 90, NY = 43, NZ=101
end module stuff_mod
