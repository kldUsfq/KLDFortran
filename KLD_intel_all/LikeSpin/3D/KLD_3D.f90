!
! KLD: calculate the Kullback-Leibler divergence in a 3D-grid cube
! Version : 1.0.0  (Fortran 90/MPI/OpenMP)
! Author: Luis Rincon, USFQ, Aug 27, 2015
! Revision: Luis Rincon, USFQ, Sep 10,2017
!

program KLD_3D

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

!     CALL cube_output_3D


!   DEALLOCATE(wfn,int)   
    DEALLOCATE(wfn)   

! Print elapsed CPU time

    CALL cpu_time(t2)

    print *, ' Elapsed CPU time =  ', t2 - t1,' sec' 

! Print elapsed real time

    CALL system_clock(count_1,count_rate,count_max)

    finish_time = count_1/count_rate 

    print *, ' Elapsed Real time = ', finish_time-start_time,' sec' 

    print *, ' End of the program '

end program KLD_3D

