### ansu_fortran
A Fortran program for calculating approximately neutral surfaces. 

This program transforms a single poorly adjusted surface (such as a 
potential density surface) into an 
approximately neutral surface with minimal error during a number of iterations. 

References: 
Klocker, A., McDougall, T., Jackett, D., 2009. A new method for forming approximately neutral
surfaces. Ocean Science 5, 155–172.
Riha, S., McDougall, T. J., Barker, P.M. (unpublished manuscript, 2014): Improvements of an algorithm for 
forming approximately neutral surfaces. http://www.hoitaus.com/drupal/files/publications/paper_syd1_draft.pdf


#### DIRECTORY STRUCTURE AND BUILD:

##### ansu.f90: Fortran module providing the subroutine optimize_surface.f90. Depends on the 
GSW Oceanographic Toolbox and the LSQR algorithm (see below). To build the module, edit the 
Makefile to point to the dependencies, and type
	make ansu.o
	
##### run.f90: An executable program that illustrates the use of optimize_surface.f90. Additionaly
to the GSW Toolbox and LSQR, the NetCDF library must be available to read/write data. The input
data set (10 MB) can be downloaded from http://www.hoitaus.com/drupal/files/data/os_input.nc
To build the executable, edit the Makefile to point to the dependencies, and type
	make run


#### DEPENDENCIES:
ansu.f90 depends on two third-party libraries.

1) Gibbs-SeaWater (GSW) Oceanographic Toolbox (Fortran library)
This library is not distributed with ansu. It is available at
http://www.teos-10.org/software.htm

2) LSQR: Sparse Equations and Least Squares
Provided by Systems Optimization Laboratory, Stanford 
University the terms of the OSI Common Public License (CPL)
Reference: C. C. Paige and M. A. Saunders, LSQR: An algorithm for sparse
linear equations and sparse least squares, TOMS 8(1), 43-71 (1982). 
This library is distributed with ansu.

3) run.f90 additionaly requires NetCDF.
http://www.unidata.ucar.edu/software/netcdf/


#### REPOSITORY:
https://github.com/poidl/ansu_fortran

[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/1b6b47d26b067861a6dbf1387417841f "githalytics.com")](http://githalytics.com/poidl/ansu_fortran.git)
