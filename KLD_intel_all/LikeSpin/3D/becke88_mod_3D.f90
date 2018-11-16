MODULE becke88_mod_3D

!
! NUMERICAL INTEGRATION PROGRAM (MPI/OPENMP)
! Author: Luis Rincon, NIST, Dec 15, 2010
! Revision: Luis Rincon, ULA, Jun 09, 2015
!           Luis Rincon, USFQ, Sep 01, 2017
! Reference: AD Becke, JCP, 88, 2547-2553 (1988)
! 

!!!!! call to MPI libraries
USE mpi

!!!!! call to openMP libraries
USE omp_lib

USE var_mod_3D
USE sphere_lebedev_rule_3D
USE atoms_par_3D
USE rho_mod_3D

IMPLICIT NONE

integer, PUBLIC, parameter :: rule_lebedev=7

CONTAINS

SUBROUTINE becke88_3D

IMPLICIT NONE

integer :: ia,ja,ka,npoints,kpoint
double precision :: r_bs,r_gc
double precision, allocatable, dimension(:) :: x,y,z,w_becke
integer (kind=4) :: available_lebedev, order_lebedev, precision_lebedev
real (kind=8),  allocatable, dimension(:) :: x_lebedev,y_lebedev,z_lebedev,w_lebedev
double precision :: x_gc,w_gc,wgc
double precision, parameter :: PI = 3.1415926535897932384626433832795
integer :: IJK,III,JJJ,KKK
double precision :: X1,Y1,Z1,KLD
integer, parameter :: root_process = 0 ! let process 0 be the root process.
integer :: num_procs,my_id,ierr
double precision, allocatable, dimension(:,:,:) :: partial_int
double precision, allocatable, dimension(:,:,:) :: total_int
double precision :: R12
double precision,parameter :: DomeSigma = 4.5
 
!
!  call to the subroutines that generate Lebedev grids for integration on a sphere.
!  the subroutines are included in the file SPHERE_LEBEDEV_RULE
!
!  Author:
!
!    Dmitri Laikov
!
!  Reference:
!
!    Vyacheslav Lebedev, Dmitri Laikov,
!    A quadrature formula for the sphere of the 131st
!    algebraic order of accuracy,
!    Russian Academy of Sciences Doklady Mathematics,
!    Volume 59, Number 3, 1999, pages 477-481.
!
     available_lebedev = available_table ( rule_lebedev )
     if ( available_lebedev == 1 ) then
         order_lebedev = order_table ( rule_lebedev )
         write(*,*) ' order_lebedev = ', order_lebedev
         precision_lebedev = precision_table ( rule_lebedev )
         write(*,*) ' precision_lebedev ',precision_lebedev
         ALLOCATE(x_lebedev(order_lebedev),y_lebedev(order_lebedev),z_lebedev(order_lebedev),w_lebedev(order_lebedev))
         CALL ld_by_order ( order_lebedev, x_lebedev, y_lebedev, z_lebedev, w_lebedev )
     end if

     npoints = 0 
     DO ia = 1,wfn(1)%NAT
        npoints = npoints + get_nbc(wfn(1)%AT(ia))*order_lebedev 
     END DO
     write(*,*) '   Number of points  =  ',npoints

     ALLOCATE(x(npoints),y(npoints),z(npoints),w_becke(npoints))
     kpoint = 0
     DO ia = 1,wfn(1)%NAT
! Bragg-Slater radius of atom ia
        r_bs = get_bsr(wfn(1)%AT(ia))
! radial loop: Gauss-Chevichev
        DO ja = 1,get_nbc(wfn(1)%AT(ia))
!    Gauss-Chevichev quadrature grid
        x_gc = DCOS((dfloat(ja)/dfloat(get_nbc(wfn(1)%AT(ia))+1))*PI)
        w_gc = (PI/dfloat(get_nbc(wfn(1)%AT(ia))+1)) * DSIN((dfloat(ja)/dfloat(get_nbc(wfn(1)%AT(ia))+1))*PI)**2
        wgc = 2.0 * r_bs**3 * w_gc * sqrt(((1.0+x_gc)**3)/((1.0-x_gc)**9))
! Gauss-Chevichev radius
           r_gc = r_bs * ((1.0+x_gc)/(1.0-x_gc))
! angular loop: Levedev's quadrature
              DO ka = 1,order_lebedev
                 kpoint = kpoint + 1 
                 x(kpoint) = wfn(1)%X0(ia) + x_lebedev(ka) * r_gc
                 y(kpoint) = wfn(1)%Y0(ia) + y_lebedev(ka) * r_gc
                 z(kpoint) = wfn(1)%Z0(ia) + z_lebedev(ka) * r_gc
                 w_becke(kpoint) = 4.0 * PI * w_lebedev(ka) * wv(ia,x(kpoint),y(kpoint),z(kpoint)) * wgc
              END DO
! end angular loop
        END DO
! end radial loop
     END DO

     DEALLOCATE(x_lebedev,y_lebedev,z_lebedev,w_lebedev)

! start the MPI paralelization

     call MPI_INIT(ierr)

! find out MY process ID, and how many processes were started.

     call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)

     write ( *, '(a)' ) ' '
     write ( *, '(a,i8)' ) '  The number of processors available = ', omp_get_num_procs ( )
     write ( *, '(a,i8)' ) '  The number of threads available    = ', omp_get_max_threads ( )
     write ( *, '(a)' ) ' '

     ALLOCATE(total_int(INX,INY,INZ))
     ALLOCATE(partial_int(INX,INY,INZ))

! the integration information is sent from the root to all others process

     call MPI_BCAST(npoints,1,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(x,npoints,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(y,npoints,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(z,npoints,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(w_becke,npoints,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)

! the cube information is sent from the root to all others process   

     call MPI_BCAST(INX,1,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr) 
     call MPI_BCAST(INY,1,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(INZ,1,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(XX0,1,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(YY0,1,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(ZZ0,1,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(IXA,1,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(IXB,1,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(IXC,1,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(IYA,1,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(IYB,1,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(IYC,1,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(IZA,1,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(IZB,1,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(IZC,1,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)

! the WFN information is sent from the root to all others process 

     call MPI_BCAST(wfn(1)%NMO,1,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr) 
     call MPI_BCAST(wfn(1)%NPO,1,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(wfn(1)%NAT,1,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(wfn(1)%NEL,1,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(wfn(1)%NOCC,wfn(1)%NMO,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(wfn(1)%NA,wfn(1)%NMO,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(wfn(1)%X0,wfn(1)%NAT,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(wfn(1)%Y0,wfn(1)%NAT,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(wfn(1)%Z0,wfn(1)%NAT,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(wfn(1)%IX,wfn(1)%NPO,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(wfn(1)%NX,wfn(1)%NPO,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(wfn(1)%ALPHA,wfn(1)%NPO,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(wfn(1)%C,wfn(1)%NPO*wfn(1)%NMO,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)

     IF (NWFN==2) THEN

        call MPI_BCAST(wfn(2)%NMO,1,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(wfn(2)%NPO,1,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(wfn(2)%NAT,1,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(wfn(2)%NEL,1,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(wfn(2)%NOCC,wfn(2)%NMO,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(wfn(2)%NA,wfn(2)%NMO,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(wfn(2)%X0,wfn(2)%NAT,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(wfn(2)%Y0,wfn(2)%NAT,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(wfn(2)%Z0,wfn(2)%NAT,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(wfn(2)%IX,wfn(2)%NPO,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(wfn(2)%NX,wfn(2)%NPO,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(wfn(2)%ALPHA,wfn(2)%NPO,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(wfn(2)%C,wfn(2)%NPO*wfn(2)%NMO,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)

     END IF
 
     partial_int = 0.d0
     DO IJK = (my_id+1),npoints,num_procs
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(X1,Y1,Z1,III,JJJ,KKK,KLD)
        DO III = 1,INX
           DO JJJ = 1,INY
              DO KKK = 1,INZ
                 X1 = XX0 + float(III-1)*IXA+float(JJJ-1)*IXB+float(KKK-1)*IXC
                 Y1 = YY0 + float(III-1)*IYA+float(JJJ-1)*IYB+float(KKK-1)*IYC
                 Z1 = ZZ0 + float(III-1)*IZA+float(JJJ-1)*IZB+float(KKK-1)*IZC
                 R12 = DSQRT((x(IJK)-X1)**2+(y(IJK)-Y1)**2+(z(IJK)-Z1)**2)
                 IF (R12 <= DomeSigma) THEN
                 IF(NWFN==1) THEN
                   KLD = KLD1_FUNC(x(IJK),y(IJK),z(IJK),X1,Y1,Z1,wfn(1)%NMO,wfn(1)%NEL,wfn(1)%NOCC,wfn(1)%NA,wfn(1)%NPO,&
&                         wfn(1)%X0,wfn(1)%Y0,wfn(1)%Z0,wfn(1)%IX,wfn(1)%NX,wfn(1)%ALPHA,wfn(1)%C)
                 ELSE IF(NWFN==2) THEN
                    KLD = KLD2_FUNC(x(IJK),y(IJK),z(IJK),X1,Y1,Z1,wfn(1)%NMO,wfn(1)%NEL,wfn(1)%NOCC,wfn(1)%NA,wfn(1)%NPO,&
&                   wfn(1)%X0,wfn(1)%Y0,wfn(1)%Z0,wfn(1)%IX,wfn(1)%NX,wfn(1)%ALPHA,wfn(1)%C, &
&                   wfn(2)%NMO,wfn(2)%NOCC,wfn(2)%NEL,wfn(2)%NA,wfn(2)%NPO,             &
&                   wfn(2)%X0,wfn(2)%Y0,wfn(2)%Z0,wfn(2)%IX,wfn(2)%NX,wfn(2)%ALPHA,wfn(2)%C)
                 END IF
                 ELSE
                 KLD = 0.0
                 END IF
                 partial_int(III,JJJ,KKK) = partial_int(III,JJJ,KKK) + w_becke(IJK)*KLD
              END DO
          END DO
        END DO
!$OMP END PARALLEL DO
     END DO

! engage in a reduction in which all partials_int are combined, and the
! total int appears in the root process

     call MPI_REDUCE(partial_int,total_int,INX*INY*INZ,MPI_DOUBLE_PRECISION,MPI_SUM,root_process,MPI_COMM_WORLD,ierr) 

! end MPI paralelization
     call MPI_FINALIZE(ierr)


     if(my_id==0) then
        OPEN(UNIT=6,FILE=output_file)
        DO III = 1,INX
           DO JJJ = 1,INY
              DO KKK = 1,INZ
                 X1 = XX0 + float(III-1)*IXA+float(JJJ-1)*IXB+float(KKK-1)*IXC
                 Y1 = YY0 + float(III-1)*IYA+float(JJJ-1)*IYB+float(KKK-1)*IYC
                 Z1 = ZZ0 + float(III-1)*IZA+float(JJJ-1)*IZB+float(KKK-1)*IZC
!                write(*,*) X1,Y1,Z1,total_int(III,JJJ,KKK)
                 write(6,'(4(F12.6,4X))') X1,Y1,Z1,total_int(III,JJJ,KKK)
             END DO
          END DO
       END DO
       CLOSE(UNIT=6)
     end if

     DEALLOCATE(total_int)
     DEALLOCATE(x,y,z,w_becke)

END SUBROUTINE becke88_3D

DOUBLE PRECISION FUNCTION wv(ia,x,y,z) 

IMPLICIT NONE

INTEGER, INTENT(IN) :: ia
DOUBLE PRECISION, INTENT(IN) :: x,y,z

integer, parameter :: mu_max = 3
integer :: i,li,lj
double precision :: ri,rj,Rij,muij,rbsi,rbsj,chij,uij,aij
DOUBLE PRECISION, allocatable, dimension(:) :: ww

ALLOCATE(ww(wfn(1)%NAT))
ww = 1.0
DO li=2,wfn(1)%NAT
   DO lj=1,li-1
      ri = dsqrt((x-wfn(1)%X0(li))**2+(y-wfn(1)%Y0(li))**2+(z-wfn(1)%Z0(li))**2)
      rj = dsqrt((x-wfn(1)%X0(lj))**2+(y-wfn(1)%Y0(lj))**2+(z-wfn(1)%Z0(lj))**2)
      Rij = dsqrt((wfn(1)%X0(li)-wfn(1)%X0(lj))**2+(wfn(1)%Y0(li)-wfn(1)%Y0(lj))**2 &
&     +(wfn(1)%Z0(li)-wfn(1)%Z0(lj))**2)
      muij = (ri - rj)/Rij
! estimation of aij
      rbsi = get_bsr(wfn(1)%AT(li))
      if(wfn(1)%AT(li)/='H ') rbsi = 2.0*rbsi
      rbsj = get_bsr(wfn(1)%AT(lj))
      if(wfn(1)%AT(lj)/='H ') rbsj = 2.0*rbsi
      chij = rbsi/rbsj
      uij = (chij-1.0)/(chij+1.0)
      aij = uij/(uij**2 - 1.0)
      if(aij>0.50) aij=0.50
      if(aij<-0.50) aij=-0.50
      muij = muij + aij * (1.0 - muij**2)
      DO i = 1,mu_max
          muij = 1.50 * muij - 0.50 * muij**3
      END DO
      ww(li) = 0.50 * ww(li) * (1.0 - muij)
      ww(lj) = 0.50 * ww(lj) * (1.0 + muij)
   END DO
END DO
wv = ww(ia)/SUM(ww)
DEALLOCATE(ww)
END FUNCTION wv
END MODULE becke88_mod_3D
