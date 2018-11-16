
1d:		../1D/*.*90
		cd ../1D/ && make
2d:		../2D/*.*90
		cd ../2D/ && make
3d:		../3D/*.*90
		cd ../3D/ && make
all:		
		cd ../1D/ && make
		cd ../2D/ && make
		cd ../3D/ && make

#####

objects=../1D/var_mod_1D.o ../1D/atoms_par_1D.o ../1D/lm_mod_1D.o ../1D/sphere_lebedev_rule_1D.o ../1D/wfn_mod_1D.o ../1D/cube_mod_1D.o ../1D/rho_mod_1D.o ../1D/becke88_mod_1D.o ../2D/var_mod_2D.o ../2D/atoms_par_2D.o ../2D/lm_mod_2D.o ../2D/sphere_lebedev_rule_2D.o ../2D/wfn_mod_2D.o ../2D/cube_mod_2D.o ../2D/rho_mod_2D.o ../2D/becke88_mod_2D.o ../3D/var_mod_3D.o ../3D/atoms_par_3D.o ../3D/lm_mod_3D.o ../3D/sphere_lebedev_rule_3D.o ../3D/wfn_mod_3D.o ../3D/cube_mod_3D.o ../3D/rho_mod_3D.o ../3D/becke88_mod_3D.o ./main.o 
modules=../1D/var_mod_1D.mod ../1D/atoms_par_1D.mod ../1D/lm_mod_1D.mod ../1D/sphere_lebedev_rule_1D.mod ../1D/wfn_mod_1D.mod ../1D/cube_mod_1D.mod ../1D/rho_mod_1D.mod ../1D/becke88_mod_1D.mod ../2D/var_mod_2D.mod ../2D/atoms_par_2D.mod ../2D/lm_mod_2D.mod ../2D/sphere_lebedev_rule_2D.mod ../2D/wfn_mod_2D.mod ../2D/cube_mod_2D.mod ../2D/rho_mod_2D.mod ../2D/becke88_mod_2D.mod ../3D/var_mod_3D.mod ../3D/atoms_par_3D.mod ../3D/lm_mod_3D.mod ../3D/sphere_lebedev_rule_3D.mod ../3D/wfn_mod_3D.mod ../3D/cube_mod_3D.mod ../3D/rho_mod_3D.mod ../3D/becke88_mod_3D.mod 
f90comp= mpiifort
switch= -ipo -O3 -parallel -no-prec-div -qopenmp -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -mkl=sequential
#switch= -O3 -parallel -no-prec-div -qopenmp -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -mkl=sequential
LIBS = /opt/intel/compilers_and_libraries/linux/mkl/lib/intel64
#CPATH=../1D/,../2D/,../3D/
#CPATH=../1D/
#CPATH=$(CPATH):../2D/

KLD_ALL:	$(objects) 
		$(f90comp) $(objects) $(lapack) $(switch) -o KLD_ALL -L$(MKLROOT)/lib/intel64 -L$(mpilibs) 
main.o: 	main.f90 
		$(f90comp) -I../1D/ -I../2D/ -I../3D/ -c $(switch) main.f90 

clean:		
		cd ../1D/ && make clean
		cd ../2D/ && make clean
		cd ../3D/ && make clean
