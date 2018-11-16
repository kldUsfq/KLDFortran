MODULE rho_mod_2D

IMPLICIT NONE

CONTAINS

     DOUBLE PRECISION FUNCTION gaussian(NX,xx,yy,zz,ALPHA)

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: NX
     DOUBLE PRECISION, INTENT(IN) :: xx,yy,zz,ALPHA

     DOUBLE PRECISION :: pro,ri

        if(NX==1) pro = 1.0
        if(NX==2) pro = xx
        if(NX==3) pro = yy
        if(NX==4) pro = zz
        if(NX==5) pro = xx**2
        if(NX==6) pro = yy**2
        if(NX==7) pro = zz**2
        if(NX==8) pro = xx*yy
        if(NX==9) pro = xx*zz
        if(NX==10) pro = yy*zz
        if(NX==11) pro = xx**3
        if(NX==12) pro = yy**3
        if(NX==13) pro = zz**3
        if(NX==14) pro = (xx**2)*yy
        if(NX==15) pro = (xx**2)*zz
        if(NX==16) pro = (yy**2)*zz
        if(NX==17) pro = xx*(yy**2)
        if(NX==18) pro = xx*(zz**2)
        if(NX==19) pro = yy*(zz**2)
        if(NX==20) pro = xx*yy*zz
        if(NX==21) pro = xx**4
        if(NX==22) pro = yy**4
        if(NX==23) pro = zz**4
        if(NX==24) pro = (xx**3)*yy
        if(NX==25) pro = (xx**3)*zz
        if(NX==26) pro = xx*(yy**3)
        if(NX==27) pro = (yy**3)*zz
        if(NX==28) pro = xx*(zz**3)
        if(NX==29) pro = yy*(zz**3)
        if(NX==30) pro = (xx**2)*(yy**2)
        if(NX==31) pro = (xx**2)*(zz**2)
        if(NX==32) pro = (yy**2)*(zz**2)
        if(NX==33) pro = (xx**2)*yy*zz
        if(NX==34) pro = xx*(yy**2)*zz
        if(NX==35) pro = xx*yy*(zz**2)

        ri = sqrt(xx**2+yy**2+zz**2)

        gaussian = pro*exp(-ALPHA*(ri**2))

     END FUNCTION gaussian

     DOUBLE PRECISION FUNCTION orb(II,xx,yy,zz,NPO,X0,Y0,Z0,IX,NX,ALPHA,C)

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: II,NPO
     INTEGER, INTENT(IN) , DIMENSION(:) :: IX,NX
     DOUBLE PRECISION, INTENT(IN) :: xx,yy,zz
     DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: X0,Y0,Z0,ALPHA
     DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:) :: C

     INTEGER :: i

     orb = 0.0

     do i=1,NPO

        orb = orb + c(i,ii)*gaussian(NX(i),xx-X0(IX(i)),yy-Y0(IX(i)),zz-Z0(IX(i)),ALPHA(i))

     end do

     END FUNCTION orb

     DOUBLE PRECISION FUNCTION rho(xx,yy,zz,NMO,NOCC,NPO,X0,Y0,Z0,IX,NX,ALPHA,C)

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: NMO,NPO
     INTEGER, INTENT(IN) , DIMENSION(:) :: IX,NX
     DOUBLE PRECISION, INTENT(IN) :: xx,yy,zz
     DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: X0,Y0,Z0,ALPHA,NOCC
     DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:) :: C

     INTEGER :: i
     DOUBLE PRECISION :: orbi

     rho = 0.d0
     do i=1,NMO
        if(nocc(i)>0.d0) then
          orbi = orb(i,xx,yy,zz,NPO,X0,Y0,Z0,IX,NX,ALPHA,C)
          rho = rho + nocc(i)*(orbi**2)
        end if
     end do

     END FUNCTION rho

     DOUBLE PRECISION FUNCTION rho2(x1,y1,z1,x2,y2,z2,NMO,NOCC,NA,NPO,X0,Y0,Z0,IX,NX,ALPHA,C) 

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: NMO,NPO
     INTEGER, INTENT(IN) , DIMENSION(:) :: IX,NX
     DOUBLE PRECISION, INTENT(IN) :: x1,y1,z1,x2,y2,z2
     DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: X0,Y0,Z0,ALPHA,NOCC,NA
     DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:) :: C

     INTEGER :: i,j
     DOUBLE PRECISION :: orbi1,orbi2,orbj1,orbj2

     rho2 = 0.d0
     do i=1,NMO
     if(nocc(i)>0.d0) then
        orbi1 = orb(i,x1,y1,z1,NPO,X0,Y0,Z0,IX,NX,ALPHA,C)
        orbi2 = orb(i,x2,y2,z2,NPO,X0,Y0,Z0,IX,NX,ALPHA,C)
        do j=1,NMO
        if(nocc(j)>0.d0) then
           orbj1 = orb(j,x1,y1,z1,NPO,X0,Y0,Z0,IX,NX,ALPHA,C)
           orbj2 = orb(j,x2,y2,z2,NPO,X0,Y0,Z0,IX,NX,ALPHA,C)
           rho2  = rho2 + (nocc(i)*nocc(j)-na(i)*na(j))*(orbi1*orbj2*orbi1*orbj2 - orbi1*orbj2*orbj1*orbi2)
        end if
        end do
     end if
     end do

     END FUNCTION rho2

     DOUBLE PRECISION FUNCTION rho_cond(x1,y1,z1,x2,y2,z2,NMO,NOCC,NA,NPO,X0,Y0,Z0,IX,NX,ALPHA,C)

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: NMO,NPO
     INTEGER, INTENT(IN) , DIMENSION(:) :: IX,NX
     DOUBLE PRECISION, INTENT(IN) :: x1,y1,z1,x2,y2,z2
     DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: X0,Y0,Z0,ALPHA,NOCC,NA
     DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:) :: C

     rho_cond = rho2(x1,y1,z1,x2,y2,z2,NMO,NOCC,NA,NPO,X0,Y0,Z0,IX,NX,ALPHA,C)/rho(x2,y2,z2,NMO,NOCC,NPO,X0,Y0,Z0,IX,NX,ALPHA,C)

     END FUNCTION rho_cond
 
     DOUBLE PRECISION FUNCTION f_cutoff(rho,rho_cut)

     IMPLICIT NONE

     DOUBLE PRECISION, INTENT(IN) :: rho,rho_cut

            f_cutoff = 0.5d0*(1.D0+DERF(0.5D0*DLOG(rho/rho_cut)))

     END FUNCTION f_cutoff

     DOUBLE PRECISION FUNCTION KLD1_FUNC(x1,y1,z1,x2,y2,z2,NMO,NEL,NOCC,NA,NPO,X0,Y0,Z0,IX,NX,ALPHA,C,rho_cut)

      IMPLICIT NONE

     INTEGER, INTENT(IN) :: NMO,NPO
     INTEGER, INTENT(IN) , DIMENSION(:) :: IX,NX
     DOUBLE PRECISION, INTENT(IN) :: x1,y1,z1,x2,y2,z2
     DOUBLE PRECISION, INTENT(IN) :: NEL
     DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: X0,Y0,Z0,ALPHA,NOCC,NA
     DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:) :: C
     DOUBLE PRECISION, INTENT(IN) :: rho_cut

     DOUBLE PRECISION :: rcd,rx1,fco

     rcd = rho_cond(x1,y1,z1,x2,y2,z2,NMO,NOCC,NA,NPO,X0,Y0,Z0,IX,NX,ALPHA,C)/(NEL-1.d0)
     rx1 = rho(x1,y1,z1,NMO,NOCC,NPO,X0,Y0,Z0,IX,NX,ALPHA,C)/NEL 
     fco = f_cutoff(rho(x2,y2,z2,NMO,NOCC,NPO,X0,Y0,Z0,IX,NX,ALPHA,C),rho_cut) 

     IF((rcd<=0.d0).OR.(rx1<=0.d0)) THEN
       KLD1_FUNC = 0.d0
     ELSE
       KLD1_FUNC = (NEL-1.d0)*fco*rcd*(DLOG(rcd/rx1)/DLOG(2.d0))
     END IF

     END FUNCTION KLD1_FUNC

     DOUBLE PRECISION FUNCTION KLD2_FUNC(x1,y1,z1,x2,y2,z2,NMO1,NEL1,NOCC1,NA1,NPO1,X01,Y01,Z01,IX1,NX1,ALPHA1,C1, &
&                              NMO2,NOCC2,NEL2,NA2,NPO2,X02,Y02,Z02,IX2,NX2,ALPHA2,C2,rho_cut)

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: NMO1,NPO1
     INTEGER, INTENT(IN) , DIMENSION(:) :: IX1,NX1
     INTEGER, INTENT(IN) :: NMO2,NPO2
     INTEGER, INTENT(IN) , DIMENSION(:) :: IX2,NX2
     DOUBLE PRECISION, INTENT(IN) :: x1,y1,z1,x2,y2,z2
     DOUBLE PRECISION, INTENT(IN) :: NEL1,NEL2
     DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: X01,Y01,Z01,ALPHA1,NOCC1,NA1
     DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: X02,Y02,Z02,ALPHA2,NOCC2,NA2
     DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:) :: C1,C2
     DOUBLE PRECISION, INTENT(IN) :: rho_cut

     DOUBLE PRECISION :: rcd1,rcd2,fco

     rcd1 = rho_cond(x1,y1,z1,x2,y2,z2,NMO1,NOCC1,NA1,NPO1,X01,Y01,Z01,IX1,NX1,ALPHA1,C1)/(NEL1-1.d0)
     rcd2 = rho_cond(x1,y1,z1,x2,y2,z2,NMO2,NOCC2,NA2,NPO2,X02,Y02,Z02,IX2,NX2,ALPHA2,C2)/(NEL2-1.d0)
     fco = f_cutoff(rho(x2,y2,z2,NMO1,NOCC1,NPO1,X01,Y01,Z01,IX1,NX1,ALPHA1,C1),rho_cut)

     IF((rcd1<=0.d0).OR.(rcd2<=0.d0)) THEN
       KLD2_FUNC = 0.d0
     ELSE
       KLD2_FUNC = NEL1*fco*rcd1*(DLOG(rcd1/rcd2)/DLOG(2.d0))
     END IF

     END FUNCTION KLD2_FUNC

END MODULE rho_mod_2D
