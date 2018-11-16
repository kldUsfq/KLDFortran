MODULE cube_mod_3D

USE var_mod_3D

implicit none

CONTAINS

     SUBROUTINE cube_input_3D

     IMPLICIT NONE

     OPEN(UNIT=5,FILE='cube.dat',STATUS='OLD')

!        READ(5,*) cube_file 
         READ(5,*) output_file 
         READ(5,*) XX0,YY0,ZZ0
         READ(5,*) INX,IXA,IXB,IXC
         READ(5,*) INY,IYA,IYB,IYC
         READ(5,*) INZ,IZA,IZB,IZC

     CLOSE(5)

     END SUBROUTINE cube_input_3D

     SUBROUTINE cube_output_3D

     IMPLICIT NONE

     INTEGER :: i,k,kmin,kmax
     INTEGER :: III,JJJ,KKK
     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: fint

     ALLOCATE(fint(INX*INY*INZ))

     k= 0

     DO III = 1,INX
        DO JJJ = 1,INY
           DO KKK = 1,INZ
              k = k+1
              fint(k) = int(III,JJJ,KKK)
           END DO
        END DO
     END DO 

! construccion del archivo cube

     OPEN(UNIT=6,FILE=cube_file)

     write(6,*) 'CUBE FILE'
     write(6,*)
     write(6,'(I5,3(F12.6))') wfn(1)%NAT,XX0,YY0,ZZ0
     write(6,'(I5,3(F12.6))') INX,IXA,IXB,IXC
     write(6,'(I5,3(F12.6))') INY,IYA,IYB,IYC
     write(6,'(I5,3(F12.6))') INZ,IZA,IZB,IZC
 
     DO i=1,wfn(1)%NAT
 
        WRITE(6,'(I5,4(F12.6))') NINT(wfn(1)%CHAR(I)),wfn(1)%CHAR(I),wfn(1)%X0(I),wfn(1)%Y0(I),wfn(1)%Z0(I)
 
     END DO
 
     kmax = 0
 
     DO while(kmax<INX*INY*INZ)
 
        kmin = kmax + 1
        kmax = kmin + 5
        if(kmax>INX*INY*INZ) kmax=INX*INY*INZ
        write(6,'(6E13.5)') (fint(k),k=kmin,kmax)
 
     END DO

     DEALLOCATE(fint)

     CLOSE(6)

     END SUBROUTINE cube_output_3D

END MODULE cube_mod_3D
