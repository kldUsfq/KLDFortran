!*** lm_mod: Levenberg Marquardt optimization of chi**2
!*** DATE: 20.07.2017
!*** LUIS RINCON (USFQ)
MODULE lm_mod_2D 

IMPLICIT NONE

!*** Levenberg Marqueardt parameters
INTEGER, PRIVATE :: MaxIter
DOUBLE PRECISION, PRIVATE, PARAMETER :: epsilon1 = 1D-3
DOUBLE PRECISION, PRIVATE, PARAMETER :: epsilon2 = 1D-3
DOUBLE PRECISION, PRIVATE, PARAMETER :: epsilon3 = 1D-7
DOUBLE PRECISION, PRIVATE, PARAMETER :: epsilon4 = 1D-1
DOUBLE PRECISION, PRIVATE, PARAMETER :: lambda0  = 1D-2
DOUBLE PRECISION, PRIVATE, PARAMETER :: lambda_up = 11.D0
DOUBLE PRECISION, PRIVATE, PARAMETER :: lambda_dn = 9.D0

CONTAINS

    SUBROUTINE lm_opt_2D(NNMO,NNOCC,NNA)

    INTEGER, INTENT(IN) :: NNMO
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: NNOCC
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:) :: NNA   

    LOGICAL :: Iconv 
    INTEGER :: i,j,k,Niter 
    DOUBLE PRECISION :: chi2n,chi2o,lambda,YYJ,rhoi,den1,den2,max_grad,grad,max_par,par
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: YI,YJ
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: JACB,HESS

    INTEGER :: N,NRHS,LDA,LDB,INFO
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A,B
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: IPIV

    MaxIter = 10*NNMO

    ALLOCATE(YI(NNMO),YJ(NNMO),JACB(NNMO,NNMO),HESS(NNMO,NNMO))

! initialization of the lm parameters to the occupations numbers

        DO i=1,NNMO
           YI(i) = NNOCC(i)*(NNOCC(i)-1.d0)
        END DO
        chi2n = 0.d0
        DO i=1,NNMO
           YJ(i) = 0.d0
           DO j=1,NNMO
              IF(i/=j) then
                YJ(i) = YJ(i) - NNA(i)*NNA(j)
              END IF
           END DO
           chi2n = chi2n + (YI(i)-YJ(i))**2
        END DO

        chi2o  = 0.d0
        lambda = lambda0
        Iconv  = .TRUE.
        if(chi2n<epsilon3) Iconv = .FALSE.
        Niter = 0

! main loop

        DO WHILE(Iconv)

           Niter = Niter + 1

           write(*,*) 'LM Iteration: ',Niter,' CHI2: ',chi2n,' DIF: ',chi2n-chi2o


           chi2o = chi2n

! Jacobian

          DO i=1,NNMO
             DO j=1,NNMO
                if(i==j) then
                  JACB(i,i) = 0.d0
                  DO k=1,NNMO
                     if(i/=k) then
                       JACB(i,i) = JACB(i,i) - NNA(k)
                     end if
                  end do
                else
                  JACB(i,j) = - NNA(i)
                end if
             END DO
          END DO 

! Hessian

         DO i=1,NNMO
            DO j=1,NNMO
               HESS(i,j) = 0.d0
               DO k=1,NNMO
                  HESS(i,j) = HESS(i,j) + JACB(k,i)*JACB(k,j)
               END DO
            END DO
         END DO
! compute the solution of the real system of linear equations AX=B using DGESV
! A = (HESS + lambda*DIAG(HESS))
! B = J(YI-YJ)

        N = NNMO
        NRHS = 1 
        LDA = NNMO 
        LDB = NNMO

        ALLOCATE(A(LDA,N),B(LDB,NRHS),IPIV(N))

        DO i=1,NNMO
           DO j=1,NNMO
              if(i==j) then
                A(i,i) = (1.d0+lambda)*HESS(i,i)
              else
                A(i,j) = HESS(i,j)
              end if
           END DO
        END DO

        DO i=1,NNMO
           B(i,1) = 0.d0
           DO j=1,NNMO
              B(i,1) = B(i,1) + JACB(j,i)*(YI(j)-YJ(j))
           END DO
        END DO

        CALL DGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)

        IF(INFO/=0) THEN
          write(*,*) '**** DGESV ERROR. FORTRAN STOP ****'
          STOP
        END IF

! NEW CHI2

        chi2n = 0.d0
        DO i=1,NNMO
           YYJ = 0.d0
           DO j=1,NNMO
              if(i/=j) then
                YYJ = YYJ - (NNA(i)+B(i,1))*(NNA(j)+B(j,1))
              end if
           END DO
           chi2n = chi2n + (YI(i)-YYJ)**2
        END DO

! RHOI

          den1 = 0.d0
          do i=1,NNMO
             den1 = den1 + lambda*B(i,1)*HESS(i,i)*B(i,1)
          end do 
          den2 = 0.d0
          do i=1,NNMO
             do j=1,NNMO
                den2 = den2 + B(i,1)*JACB(j,i)*(YI(j)-YJ(j))
             end do
          end do

          rhoi = (chi2o-chi2n)/(den1+den2)

! ACEPTANCE CRITERIA

          if(rhoi>epsilon4) then
            lambda = MAX(lambda/lambda_dn,1D-7)
            do i=1,NNMO
               NNA(i) = NNA(i)+B(i,1)
            end do
            do i=1,NNMO
               YJ(i) = 0.d0
               do j=1,NNMO
                  if(i/=j) then
                    YJ(i) = YJ(i) - NNA(i)*NNA(j)
                  end if
               end do
            end do
          else
            lambda = MIN(lambda*lambda_up,1D7)
          end if 

! convergence criteria
! gradient
          max_grad = 0.d0
          do i=1,NNMO
             grad = 0.d0
             do j=1,NNMO
                grad = grad + JACB(j,i)*(YI(j)-YJ(j))
             end do
             if(abs(grad)>max_grad) max_grad = abs(grad)
          end do
! parameters
          max_par = 0.d0
          do i=1,NNMO
             par = B(i,1)/YJ(i)
             if(abs(par)>max_par) max_par = abs(par)
          end do
          if(max_grad<epsilon1) Iconv=.FALSE.
          if(max_par<epsilon2)  IConv=.FALSE.
          if(chi2n<epsilon3)    Iconv=.FALSE.
          if(Niter>MaxIter) Iconv=.FALSE. 

        DEALLOCATE(A,B,IPIV)

        END DO

        write(*,*) ' chi**2  =  ',chi2n

! end main loop

    END SUBROUTINE lm_opt_2D

END MODULE lm_mod_2D
