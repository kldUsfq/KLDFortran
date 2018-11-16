!
! KLD CONSOLIDATED VERSION 
!

PROGRAM main_KLD 

IMPLICIT NONE

INTEGER :: VERSION 

! READING ARGUMENT 

OPEN (UNIT=1, FILE='DIM')
READ(1,*) VERSION
!READ(*,*) VERSION
WRITE( *, '(a)' ) ' '
WRITE( *, '(a)' ) '  THIS IS VERSION ' 
PRINT *, VERSION
WRITE( *, '(a)' ) ' '

IF (VERSION.EQ.1) THEN
   CALL KLD_1D
ELSE IF (VERSION.EQ.2) THEN
   CALL KLD_2D
ELSE
   CALL KLD_3D
END IF
END PROGRAM main_KLD

!
! KLD_1D: calculate the Kullback-Leibler divergence in a 1D-grid
! Version : 1.0.0  (Fortran 90/MPI/OpenMP)
! Author: Luis Rincon, USFQ, Aug 27, 2015
! Revision: Luis Rincon, USFQ, Sep 10,2017
!

SUBROUTINE KLD_1D

USE wfn_mod_1D
USE cube_mod_1D
USE becke88_mod_1D

IMPLICIT NONE

INTEGER :: count_0,count_1,count_rate,count_max
REAL    :: start_time,finish_time 
REAL    :: t1,t2


! call cpu time

    CALL cpu_time(t1)

! call system_clock (10/09/2017)

    CALL system_clock(count_0,count_rate,count_max)

     start_time = count_0/count_rate

! HEADER

write ( *, '(a)' ) ' '
write ( *, '(a)' ) '  PROGRAM KLD 2D                   '
write ( *, '(a)' ) '  FORTRAN90/OpenMP/MPI version     '
write ( *, '(a)' ) ' '

! read the wfn input

    CALL wfn_input_1D

! read cube parameters

    CALL cube_input_1D

! cube calculation and output

    CALL becke88_1D


    DEALLOCATE(wfn)   

! Print elapsed CPU time

    CALL cpu_time(t2)

    write(*,*) ' Elapsed CPU time =  ', t2 - t1,' sec' 

! Print elapsed real time

    CALL system_clock(count_1,count_rate,count_max)

    finish_time = count_1/count_rate 

    write(*,*) ' Elapsed Real time = ', finish_time-start_time,' sec' 

    write(*,*) ' End of the program '

END SUBROUTINE KLD_1D

!
! KLD_2D: calculate the Kullback-Leibler divergence in a 2D-grid
! Version : 1.0.0  (Fortran 90/MPI/OpenMP)
! Author: Luis Rincon, USFQ, Aug 27, 2015
! Revision: Luis Rincon, USFQ, Sep 10,2017
!

SUBROUTINE KLD_2D

USE wfn_mod_2D
USE cube_mod_2D
USE becke88_mod_2D

IMPLICIT NONE

INTEGER :: count_0,count_1,count_rate,count_max
REAL    :: start_time,finish_time 
REAL    :: t1,t2


! call cpu time

    CALL cpu_time(t1)

! call system_clock (10/09/2017)

    CALL system_clock(count_0,count_rate,count_max)

     start_time = count_0/count_rate

! HEADER

write ( *, '(a)' ) ' '
write ( *, '(a)' ) '  PROGRAM KLD 2D                   '
write ( *, '(a)' ) '  FORTRAN90/OpenMP/MPI version     '
write ( *, '(a)' ) ' '

! read the wfn input

    CALL wfn_input_2D

! read cube parameters

    CALL cube_input_2D

! cube calculation and output

    CALL becke88_2D


    DEALLOCATE(wfn)   

! Print elapsed CPU time

    CALL cpu_time(t2)

    write(*,*) ' Elapsed CPU time =  ', t2 - t1,' sec' 

! Print elapsed real time

    CALL system_clock(count_1,count_rate,count_max)

    finish_time = count_1/count_rate 

    write(*,*) ' Elapsed Real time = ', finish_time-start_time,' sec' 

    write(*,*) ' End of the program '

END SUBROUTINE KLD_2D

!
! KLD: calculate the Kullback-Leibler divergence in a 3D-grid cube
! Version : 1.0.0  (Fortran 90/MPI/OpenMP)
! Author: Luis Rincon, USFQ, Aug 27, 2015
! Revision: Luis Rincon, USFQ, Sep 10,2017
!

SUBROUTINE KLD_3D

USE var_mod_3D
USE wfn_mod_3D
USE cube_mod_3D
USE becke88_mod_3D

IMPLICIT NONE

INTEGER :: count_0,count_1,count_rate,count_max
REAL    :: start_time,finish_time 
REAL    :: t1,t2


! call cpu time

    CALL cpu_time(t1)

! call system_clock (10/09/2017)

    CALL system_clock(count_0,count_rate,count_max)

     start_time = count_0/count_rate

! HEADER

write ( *, '(a)' ) ' '
write ( *, '(a)' ) '  PROGRAM KLD                      '
write ( *, '(a)' ) '  FORTRAN90/OpenMP/MPI version     '
write ( *, '(a)' ) ' '

! read the wfn input

    CALL wfn_input_3D

! read cube parameters

    CALL cube_input_3D

! cube calculation

    CALL becke88_3D

! OUTPUT

     CALL cube_output_3D


    DEALLOCATE(wfn,int)   

! Print elapsed CPU time

    CALL cpu_time(t2)

    print *, ' Elapsed CPU time =  ', t2 - t1,' sec' 

! Print elapsed real time

    CALL system_clock(count_1,count_rate,count_max)

    finish_time = count_1/count_rate 

    print *, ' Elapsed Real time = ', finish_time-start_time,' sec' 

    print *, ' End of the program '

END SUBROUTINE KLD_3D

