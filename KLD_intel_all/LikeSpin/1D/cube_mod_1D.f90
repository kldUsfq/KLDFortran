MODULE cube_mod_1D

USE var_mod_1D

implicit none

CONTAINS

     SUBROUTINE cube_input_1D

     IMPLICIT NONE

     OPEN(UNIT=5,FILE='cube.dat',STATUS='OLD')

         READ(5,*) output_file 
         READ(5,*) XX0,YY0,ZZ0
         READ(5,*) INX,IXA,IXB,IXC
         READ(5,*) INY,IYA,IYB,IYC
         READ(5,*) INZ,IZA,IZB,IZC

     CLOSE(5)

     END SUBROUTINE cube_input_1D

END MODULE cube_mod_1D
