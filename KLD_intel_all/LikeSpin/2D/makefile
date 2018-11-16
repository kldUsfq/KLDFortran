# Start makefile
# Defining variables
objects = var_mod_2D.o atoms_par_2D.o lm_mod_2D.o sphere_lebedev_rule_2D.o wfn_mod_2D.o cube_mod_2D.o rho_mod_2D.o becke88_mod_2D.o KLD_2D.o
modules = var_mod_2D.mod atoms_par_2D.mod lm_mod_2D.mod sphere_lebedev_rule_2D.mod wfn_mod_2D.mod cube_mod_2D.mod rho_mod_2D.mod becke88_mod_2D.mod
f90comp = mpiifort
switch = -ipo -O3 -parallel -no-prec-div -qopenmp -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -mkl=sequential
#MKL Libraries
MKLROOT = /opt/intel/mkl

#LAPACK
lapack = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a

#MPI Libraries
mpilibs=/opt/intel/compilers_and_libraries_2018.3.222/linux/compiler/lib/intel64
mpilibsiomp5 = $(mpilibs)/libiomp5.a

# libraries
LIBS=/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64
# executable name
EXE = KLD2D
# Makefile
$(EXE):	 		$(objects)
			$(f90comp) $(objects)  $(lapack) $(switch) -o$(EXE) -L$(MKLROOT)/lib/intel64 -L$(mpilibs) 
var_mod_2D.o:		var_mod_2D.f90
			$(f90comp) -c $(switch) var_mod_2D.f90
atoms_par_2D.o:		atoms_par_2D.f90
			$(f90comp) -c $(switch) atoms_par_2D.f90
lm_mod_2D.o:		lm_mod_2D.f90
			$(f90comp) -c $(switch) lm_mod_2D.f90
sphere_lebedev_rule_2D.o:	sphere_lebedev_rule_2D.f90 
			$(f90comp) -c $(switch) sphere_lebedev_rule_2D.f90
wfn_mod_2D.o:		wfn_mod_2D.f90
			$(f90comp) -c $(switch) wfn_mod_2D.f90
cube_mod_2D.o:		cube_mod_2D.f90
			$(f90comp) -c $(switch) cube_mod_2D.f90
rho_mod_2D.o:		rho_mod_2D.f90
			$(f90comp) -c $(switch) rho_mod_2D.f90
becke88_mod_2D.o:	becke88_mod_2D.f90
			$(f90comp) -c $(switch) becke88_mod_2D.f90
KLD_2D.o:		KLD_2D.f90
			$(f90comp) -c $(switch) KLD_2D.f90 
# Cleaning everything
clean:
			rm $(modules) 
			rm $(objects)
			rm $(EXE)
# End makefile
