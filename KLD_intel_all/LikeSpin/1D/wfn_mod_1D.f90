MODULE wfn_mod_1D

!
! this module read the wfn input file from Gaussian or Gamess 
! Author: Luis Rincon, ULA, 03.2014 (Include wfn_type :: 24/07/2017)
!

USE var_mod_1D
USE lm_mod_1D

IMPLICIT NONE

CONTAINS

SUBROUTINE wfn_input_1D

IMPLICIT NONE

INTEGER :: II,i,j,k,kmin,kmax

    OPEN(UNIT=5,FILE='wfn.dat',STATUS='OLD')

    READ(5,*) NWFN

    ALLOCATE(wfn(NWFN)) 

    DO i=1,NWFN
       READ(5,*) wfn(i)%wfn_file
    END DO     

    CLOSE(UNIT=5)

    DO II=1,NWFN 

       OPEN(UNIT=5,FILE=wfn(II)%wfn_file,STATUS='OLD')

       READ(5,*)

!   NMO = Number of Molecular Orbitals
!   NPO = Number of Primitives Gaussian Orbitals
!   NAT = Number of Atoms

       READ(5,'(19X,I4,16X,I4,16X,I4)') wfn(II)%NMO,wfn(II)%NPO,wfn(II)%NAT

       ALLOCATE(wfn(II)%AT(wfn(II)%NAT),wfn(II)%X0(wfn(II)%NAT),wfn(II)%Y0(wfn(II)%NAT),wfn(II)%Z0(wfn(II)%NAT))
       ALLOCATE(wfn(II)%CHAR(wfn(II)%NAT),wfn(II)%IX(wfn(II)%NPO),wfn(II)%NX(wfn(II)%NPO))
       ALLOCATE(wfn(II)%ALPHA(wfn(II)%NPO),wfn(II)%C(wfn(II)%NPO,wfn(II)%NMO),wfn(II)%NOCC(wfn(II)%NMO))

!   read the atomic coordinates in atomic units


       DO i=1,wfn(II)%NAT
       READ(5,'(2X,A2,20X,3F12.8,11X,F4.2)') wfn(II)%AT(I),wfn(II)%X0(I),wfn(II)%Y0(I),wfn(II)%Z0(I),wfn(II)%CHAR(I)
       END DO

!   read the nuclear center associated with each primitive

       kmin = 1
       DO WHILE(kmin<=wfn(II)%NPO)
          kmax = kmin+19
          IF(kmax>wfn(II)%NPO) kmax=wfn(II)%NPO
          READ(5,'(21X,20(I2,X))') (wfn(II)%IX(I),I=kmin,kmax)
          kmin = kmin+20
       END DO

!   read the type of orbital for each primitive

       kmin = 1
       DO WHILE(kmin<=wfn(II)%NPO)
          kmax = kmin+19
          IF(kmax>wfn(II)%NPO) kmax=wfn(II)%NPO
          READ(5,'(21X,20(I2,X))') (wfn(II)%NX(I),I=kmin,kmax)
          kmin = kmin+20
       END DO 

!   read the orbital exponent for each primitive

       kmin = 1
       DO WHILE(kmin<=wfn(II)%NPO)
          kmax = kmin+4
          IF(kmax>wfn(II)%NPO) kmax=wfn(II)%NPO
          READ(5,'(11X,5(E13.7,X))') (wfn(II)%ALPHA(I),I=kmin,kmax)
          kmin = kmin+5
       END DO

!   read the molecular orbital coefficients

       DO I=1,wfn(II)%NMO
          READ(5,'(38X,F11.8)') wfn(II)%NOCC(I)
          kmin = 1
          DO WHILE(kmin<=wfn(II)%NPO)
             kmax = kmin+4
             IF(kmax>wfn(II)%NPO) kmax=wfn(II)%NPO
             READ(5,'(5(1X,F15.8))') (wfn(II)%C(K,I),K=kmin,kmax)
             kmin = kmin+5
          END DO
       END DO

!  

       CLOSE(UNIT=5)

       ALLOCATE(wfn(II)%NA(wfn(II)%NMO))
       wfn(II)%NEL = 0.d0
        do i=1,wfn(II)%NMO
           wfn(II)%NOCC(I) = 0.5d0*wfn(II)%NOCC(I)
           if(wfn(II)%NOCC(I)>1.d0) wfn(II)%NOCC(I)=1.d0
           if(wfn(II)%NOCC(I)<0.d0) wfn(II)%NOCC(I)=0.d0 
           wfn(II)%NEL = wfn(II)%NEL + wfn(II)%NOCC(I)
           wfn(II)%NA(I) = wfn(II)%NOCC(I)
           if(wfn(II)%NOCC(I)==1.d0) wfn(II)%NA(I) = 0.d0
        end do

        CALL lm_opt_1D(wfn(II)%NMO,wfn(II)%NOCC,wfn(II)%NA)

    END DO

! end subroutine wfn_input

END SUBROUTINE wfn_input_1D

END MODULE wfn_mod_1D
