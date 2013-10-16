FORT=gfortran
NCDIR = /home/z3439823/programs/netcdf
GSWDIR= /home/z3439823/fortran/gsw_fortran_v3_01
LSQRDIR= /home/z3439823/eclipse/workspace/ansu/third_party/lsqr2

#FLAGS= -g -O 
#FLAGS= -O3
FLAGS= -pg

# LLFLAGS: options preceeding "-Wl' are for linking at compile time, 
# succeeding ones are for linking at load time (or runtime?). In the 
# case where
# a netcdf library is installed in a default location, but I want
# to use one located in my home directory, then omitting the succeeding
# ones will cause the program to call the default library (check with ldd).
# This is surprisingly tricky...am I missing something?
LLFLAGS = -L$(NCDIR)/lib -lnetcdff -Wl,-rpath,$(NCDIR)/lib
INCFLAGS = -I$(NCDIR)/include -I$(LSQRDIR)


$(LSQRDIR)/%.o : $(LSQRDIR)/%.f90 ; ${FORT} ${FLAGS} -c -o $@ $< -J $(LSQRDIR)
%.o : %.f90 ; ${FORT} ${FLAGS} -c -o $@ $< $(INCFLAGS)

lsqrfiles= $(LSQRDIR)/lsqrDataModule.o $(LSQRDIR)/lsqrblas.o \
$(LSQRDIR)/lsqrblasInterface.o $(LSQRDIR)/lsqrDataModule.o \
$(LSQRDIR)/lsqrModule.o 

files = $(lsqrfiles) stuff.o  ncutils.o ans.o run.o 

run_exe: $(files)
	$(FORT) $(FLAGS) -o $@  $(files) \
	$(GSWDIR)/gsw_oceanographic_toolbox.o $(LLFLAGS)

clean:
	rm -f *.o *.mod
	rm -f $(LSQRDIR)/*.o $(LSQRDIR)/*.mod
	
