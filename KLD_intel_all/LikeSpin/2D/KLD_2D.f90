!
! KLD_2D: calculate the Kullback-Leibler divergence in a 2D-grid
! Version : 1.0.0  (Fortran 90/MPI/OpenMP)
! Author: Luis Rincon, USFQ, Aug 27, 2015
! Revision: Luis Rincon, USFQ, Sep 10,2017
!

program KLD_2D

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

end program KLD_2D

