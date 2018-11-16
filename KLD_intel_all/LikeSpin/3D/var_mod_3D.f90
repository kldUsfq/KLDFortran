MODULE var_mod_3D

IMPLICIT NONE

    CHARACTER(LEN=40), PUBLIC :: cube_file
    CHARACTER(LEN=40), PUBLIC :: output_file 

    TYPE wfn_type
         CHARACTER(LEN=40), PUBLIC :: wfn_file
         CHARACTER(LEN=2), PUBLIC, ALLOCATABLE, DIMENSION(:) :: AT
         INTEGER, PUBLIC :: NMO,NPO,NAT
         INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: IX,NX
         DOUBLE PRECISION, PUBLIC :: NEL
         DOUBLE PRECISION, PUBLIC, ALLOCATABLE, DIMENSION(:) :: X0,Y0,Z0,CHAR,ALPHA,NOCC,NA
         DOUBLE PRECISION, PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: C
    END TYPE wfn_type

    INTEGER, PUBLIC :: NWFN
    TYPE(wfn_type), PUBLIC, ALLOCATABLE, DIMENSION(:) :: wfn 

    INTEGER, PUBLIC :: INX,INY,INZ
    DOUBLE PRECISION, PUBLIC :: XX0,YY0,ZZ0
    DOUBLE PRECISION, PUBLIC :: IXA,IXB,IXC
    DOUBLE PRECISION, PUBLIC :: IYA,IYB,IYC
    DOUBLE PRECISION, PUBLIC :: IZA,IZB,IZC

    DOUBLE PRECISION, PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: int

END MODULE var_mod_3D
