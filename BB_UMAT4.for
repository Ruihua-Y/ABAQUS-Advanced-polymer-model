      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1 DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     2 CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     3 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     1  DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     2  TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),DROT(3,3),
     3  DFGRD0(3,3),DFGRD1(3,3)
      INTEGER I,J,K,L,KL,IERR
      DOUBLE PRECISION F(3,3),SIG(3,3),SIGP(3,3),DDS(3,3,3,3)
      DOUBLE PRECISION STN(18),STP(18),FPERT(3,3)
      DOUBLE PRECISION DT,THETA,PERT,DETFP,DETF
      DOUBLE PRECISION PAR(18)
      DO I = 1,18
        PAR(I) = PROPS(I)
      END DO
      DT    = DTIME
      THETA = TEMP
      IF (DT .LE. 0.0D0) DT = 1.0D-6
      DO I = 1,3
        DO J = 1,3
          F(I,J) = DFGRD1(I,J)
        END DO
      END DO
      DO I = 1,18
        STN(I) = STATEV(I)
      END DO
      CALL STRESS_UPDATE(F,STN,PAR,DT,THETA,SIG,STP,IERR)
      IF (IERR .NE. 0) THEN
        PNEWDT = 0.25D0
        DO I = 1,NTENS
          STRESS(I) = 0.0D0
          DO J = 1,NTENS
            DDSDDE(I,J) = 0.0D0
          END DO
        END DO
        RETURN
      END IF
      DO I = 1,18
        STATEV(I) = STP(I)
      END DO
      CALL TENS2VOIGT(SIG,STRESS)
      PERT = 1.0D-7
      DETF = DET3(F)
      DO I = 1,3
        DO J = 1,3
          DO K = 1,3
            DO L = 1,3
              DDS(I,J,K,L) = 0.0D0
            END DO
          END DO
        END DO
      END DO
      DO KL = 1,6
        CALL PERTURB_F(F,KL,PERT,FPERT)
        CALL STRESS_UPDATE(FPERT,STN,PAR,DT,THETA,SIGP,STP,IERR)
        IF (IERR .NE. 0) THEN
          PNEWDT = 0.25D0
          RETURN
        END IF
        DETFP = DET3(FPERT)
        DO I = 1,3
          DO J = 1,3
            DDS(I,J,KDIR1(KL),KDIR2(KL)) =
     1        (DETFP/DETF*SIGP(I,J) - SIG(I,J))/PERT
          END DO
        END DO
      END DO
      CALL DDS_TO_VOIGT(DDS,DDSDDE)
      CALL GEO_TANGENT(SIG,DDSDDE)
      RETURN
      END

      INTEGER FUNCTION KDIR1(KL)
      INTEGER KL
      IF (KL .EQ. 1) KDIR1 = 1
      IF (KL .EQ. 2) KDIR1 = 2
      IF (KL .EQ. 3) KDIR1 = 3
      IF (KL .EQ. 4) KDIR1 = 1
      IF (KL .EQ. 5) KDIR1 = 1
      IF (KL .EQ. 6) KDIR1 = 2
      RETURN
      END

      INTEGER FUNCTION KDIR2(KL)
      INTEGER KL
      IF (KL .EQ. 1) KDIR2 = 1
      IF (KL .EQ. 2) KDIR2 = 2
      IF (KL .EQ. 3) KDIR2 = 3
      IF (KL .EQ. 4) KDIR2 = 2
      IF (KL .EQ. 5) KDIR2 = 3
      IF (KL .EQ. 6) KDIR2 = 3
      RETURN
      END

      SUBROUTINE PERTURB_F(F,KL,PERT,FP)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER KL,I,J,K,L
      DOUBLE PRECISION F(3,3),FP(3,3),PERT
      K = KDIR1(KL)
      L = KDIR2(KL)
      DO I = 1,3
        DO J = 1,3
          FP(I,J) = F(I,J) + 0.5D0*PERT*
     1      (DELTA(I,K)*F(L,J) + DELTA(I,L)*F(K,J))
        END DO
      END DO
      RETURN
      END

      DOUBLE PRECISION FUNCTION DELTA(I,J)
      INTEGER I,J
      IF (I .EQ. J) THEN
        DELTA = 1.0D0
      ELSE
        DELTA = 0.0D0
      END IF
      RETURN
      END

      DOUBLE PRECISION FUNCTION DET3(A)
      DOUBLE PRECISION A(3,3)
      DET3 = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
     1     - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
     2     + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      RETURN
      END

      SUBROUTINE TENS2VOIGT(A,AV)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION A(3,3),AV(6)
      AV(1) = A(1,1)
      AV(2) = A(2,2)
      AV(3) = A(3,3)
      AV(4) = A(1,2)
      AV(5) = A(1,3)
      AV(6) = A(2,3)
      RETURN
      END

      SUBROUTINE VOIGT2TENS(AV,A)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION A(3,3),AV(6)
      A(1,1) = AV(1)
      A(2,2) = AV(2)
      A(3,3) = AV(3)
      A(1,2) = AV(4)
      A(2,1) = AV(4)
      A(1,3) = AV(5)
      A(3,1) = AV(5)
      A(2,3) = AV(6)
      A(3,2) = AV(6)
      RETURN
      END

      SUBROUTINE DDS_TO_VOIGT(DDS,DDSDDE)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION DDS(3,3,3,3),DDSDDE(6,6)
      INTEGER I,J,K,L,II,JJ
      INTEGER MAP(3,3)
      DATA MAP /1,4,5, 4,2,6, 5,6,3/
      DO II = 1,6
        DO JJ = 1,6
          DDSDDE(II,JJ) = 0.0D0
        END DO
      END DO
      DO I = 1,3
        DO J = I,3
          II = MAP(I,J)
          DO K = 1,3
            DO L = 1,3
              JJ = MAP(K,L)
              DDSDDE(II,JJ) = DDSDDE(II,JJ) +
     1          0.5D0*(DDS(I,J,K,L)+DDS(J,I,K,L))
            END DO
          END DO
        END DO
      END DO
      RETURN
      END

      SUBROUTINE GEO_TANGENT(SIG,DDSDDE)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION SIG(3,3),DDSDDE(6,6)
      INTEGER II,JJ,I,J,K,L
      DO II = 1,6
        I = KDIR1(II)
        J = KDIR2(II)
        DO JJ = 1,6
          K = KDIR1(JJ)
          L = KDIR2(JJ)
          DDSDDE(II,JJ) = DDSDDE(II,JJ) - 0.5D0*(DELTA(I,K)*SIG(J,L)
     1      + DELTA(I,L)*SIG(J,K) + SIG(I,K)*DELTA(J,L)
     2      + SIG(I,L)*DELTA(J,K))
        END DO
      END DO
      RETURN
      END

      SUBROUTINE STRESS_UPDATE(F,STN,PAR,DT,THETA,SIG,STP,IERR)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION F(3,3),STN(18),PAR(18),SIG(3,3),STP(18),DT,THETA
      INTEGER IERR
      DOUBLE PRECISION C(3,3),LAM2(3),NN(3,3),LAM(3),EPSP(3)
      DOUBLE PRECISION EPS(3,3),EPSN(3,3),RMAT(3,3),UINV(3,3)
      DOUBLE PRECISION TAU(3,3),TAUB(3,3),TAUC(3,3)
      DOUBLE PRECISION TAUBDEVN(3,3),TAUBDEV(3,3)
      DOUBLE PRECISION ENN(3),EN(3),EPSTRI(3),EPSNP(3)
      DOUBLE PRECISION TAUC_P(3)
      DOUBLE PRECISION JJ,OMEGA
      INTEGER I,J,K,IT
      DOUBLE PRECISION D,EE,TAUBM0,MU0,GAMMA0,LAMC
      DOUBLE PRECISION XI,TMP
      DOUBLE PRECISION TB_OCT,MU,ZETA_DEV(3,3),ZETA_OCT,DEPDEV(3,3)
      DOUBLE PRECISION TB_LO,TB_HI,TB_NEW
      DOUBLE PRECISION LNJ
      DOUBLE PRECISION A2,RES,KT
      DOUBLE PRECISION R,KB,GB,KAPPA,MUSTAR,CVISC,DH0,VS,VP,NS,GS,CS,
     1  TSTAR,TINF,ALPHA,ETA,ACRIT,BCRIT
      IERR = 0
      R      = PAR(1)
      KB     = PAR(2)
      GB     = PAR(3)
      KAPPA  = PAR(4)
      MUSTAR = PAR(5)
      CVISC  = PAR(6)
      DH0    = PAR(7)
      VS     = PAR(8)
      VP     = PAR(9)
      NS     = PAR(10)
      GS     = PAR(11)
      CS     = PAR(12)
      TSTAR  = PAR(13)
      TINF   = PAR(14)
      ALPHA  = PAR(15)
      ETA    = PAR(16)
      ACRIT  = PAR(17)
      BCRIT  = PAR(18)
      A2 = ALPHA*ALPHA
      D  = NS*KB*THETA/2.0D0
      EE = (1.0D0+ETA)*(1.0D0-A2)
      TAUBM0 = -NS*KB*THETA*((1.0D0+A2)/(1.0D0-3.0D0*A2)
     1          + ETA/(1.0D0+ETA))
      MU0 = MUSTAR*EXP(CVISC/(THETA-TINF) - CVISC/(TSTAR-TINF)
     1      + DH0/(R*THETA) - DH0/(R*TSTAR))
      GAMMA0 = GS*EXP(CS/(THETA-TINF) - CS/(TSTAR-TINF))
      LAMC = ACRIT + BCRIT*(THETA-237.15D0)
      DO I = 1,3
        DO J = 1,3
          C(I,J) = 0.0D0
          DO K = 1,3
            C(I,J) = C(I,J) + F(K,I)*F(K,J)
          END DO
        END DO
      END DO
      CALL JACOBI(C,LAM2,NN,IERR)
      IF (IERR .NE. 0) RETURN
      DO I = 1,3
        IF (LAM2(I) .LE. 0.0D0) THEN
          IERR = 2
          RETURN
        END IF
        LAM(I)  = SQRT(LAM2(I))
        EPSP(I) = LOG(LAM(I))
      END DO
      JJ = LAM(1)*LAM(2)*LAM(3)
      LNJ = LOG(JJ)
      DO I = 1,3
        DO J = 1,3
          UINV(I,J) = 0.0D0
          DO K = 1,3
            UINV(I,J) = UINV(I,J) + (1.0D0/LAM(K))*NN(I,K)*NN(J,K)
          END DO
        END DO
      END DO
      DO I = 1,3
        DO J = 1,3
          RMAT(I,J) = 0.0D0
          DO K = 1,3
            RMAT(I,J) = RMAT(I,J) + F(I,K)*UINV(K,J)
          END DO
        END DO
      END DO
      DO I = 1,3
        DO J = 1,3
          EPS(I,J) = 0.0D0
          DO K = 1,3
            EPS(I,J) = EPS(I,J) + EPSP(K)*NN(I,K)*NN(J,K)
          END DO
        END DO
      END DO
      CALL VOIGT2TENS(STN(1),EPSN)
      CALL VOIGT2TENS(STN(7),TAUBDEVN)
      DO I = 1,3
        ENN(I)   = STN(12+I)
        EPSNP(I) = STN(15+I)
      END DO
      CALL DEVIATORIC(EPS,EPSN,DEPDEV)
      DO I = 1,3
        DO J = 1,3
          ZETA_DEV(I,J) = TAUBDEVN(I,J) + 2.0D0*GB*DEPDEV(I,J)
        END DO
      END DO
      CALL OCTAHEDRAL(ZETA_DEV,ZETA_OCT)
      OMEGA = EXP(-VP*KAPPA*LNJ/(R*THETA))
      TB_OCT = ZETA_OCT
      TB_LO = 0.0D0
      TB_HI = ZETA_OCT
      DO IT = 1,80
        XI  = VS*TB_OCT/(2.0D0*R*THETA)
        MU  = MU0*OMEGA*XOVSINH(XI)
        RES = TB_OCT - ZETA_OCT/(1.0D0 + 2.0D0*GB*DT/MU)
        IF (ABS(RES) .LT. 1.0D-10*MAX(1.0D0,ABS(TB_OCT))) GOTO 110
        IF (RES .GT. 0.0D0) THEN
          TB_HI = TB_OCT
        ELSE
          TB_LO = TB_OCT
        END IF
        KT = 1.0D0 - 2.0D0*DVARPI(MU0,OMEGA,XI,VS,R,THETA)
     1       *ZETA_OCT*GB*DT/(2.0D0*GB*DT+MU)**2
        TB_NEW = TB_OCT - RES/KT
        IF (TB_NEW .LE. TB_LO .OR. TB_NEW .GE. TB_HI
     1      .OR. KT .LE. 0.0D0) THEN
          TB_OCT = 0.5D0*(TB_LO + TB_HI)
        ELSE
          TB_OCT = TB_NEW
        END IF
      END DO
      IERR = 1
      RETURN
  110 CONTINUE
      XI = VS*TB_OCT/(2.0D0*R*THETA)
      MU = MU0*OMEGA*XOVSINH(XI)
      TMP = MU/(MU+2.0D0*GB*DT)
      DO I = 1,3
        DO J = 1,3
          TAUBDEV(I,J) = TMP*ZETA_DEV(I,J)
        END DO
      END DO
      TMP = KAPPA*LNJ + TAUBM0
      DO I = 1,3
        DO J = 1,3
          TAUB(I,J) = TAUBDEV(I,J)
        END DO
        TAUB(I,I) = TAUB(I,I) + TMP
      END DO
      DO I = 1,3
        EPSTRI(I) = ENN(I) + (EPSP(I) - EPSNP(I))
      END DO
      CALL CONFORMATIONAL(DT,EPSTRI,EPSNP,ENN,D,EE,A2,ETA,GAMMA0,LAMC,
     1                    EN,TAUC_P,IERR)
      IF (IERR .NE. 0) RETURN
      DO I = 1,3
        DO J = 1,3
          TAUC(I,J) = 0.0D0
          DO K = 1,3
            TAUC(I,J) = TAUC(I,J) + TAUC_P(K)*NN(I,K)*NN(J,K)
          END DO
        END DO
      END DO
      DO I = 1,3
        DO J = 1,3
          TAU(I,J) = TAUB(I,J) + TAUC(I,J)
        END DO
      END DO
      CALL ROTATE_STRESS(RMAT,TAU,SIG)
      DO I = 1,3
        DO J = 1,3
          SIG(I,J) = SIG(I,J)/JJ
        END DO
      END DO
      CALL TENS2VOIGT(EPS,STP(1))
      CALL TENS2VOIGT(TAUBDEV,STP(7))
      DO I = 1,3
        STP(12+I) = EN(I)
        STP(15+I) = EPSP(I)
      END DO
      RETURN
      END

      SUBROUTINE ROTATE_STRESS(R,T,A)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION R(3,3),T(3,3),A(3,3),TMP(3,3)
      INTEGER I,J,K
      DO I = 1,3
        DO J = 1,3
          TMP(I,J) = 0.0D0
          DO K = 1,3
            TMP(I,J) = TMP(I,J) + T(I,K)*R(J,K)
          END DO
        END DO
      END DO
      DO I = 1,3
        DO J = 1,3
          A(I,J) = 0.0D0
          DO K = 1,3
            A(I,J) = A(I,J) + R(I,K)*TMP(K,J)
          END DO
        END DO
      END DO
      RETURN
      END

      SUBROUTINE OCTAHEDRAL(A,OCT)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION A(3,3),OCT
      INTEGER I,J
      DOUBLE PRECISION S
      S = 0.0D0
      DO I = 1,3
        DO J = 1,3
          S = S + A(I,J)*A(I,J)
        END DO
      END DO
      OCT = SQRT(S/3.0D0)
      RETURN
      END

      SUBROUTINE DEVIATORIC(EPS,EPSN,DEPDEV)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION EPS(3,3),EPSN(3,3),DEPDEV(3,3)
      INTEGER I,J
      DOUBLE PRECISION TR
      TR = (EPS(1,1)-EPSN(1,1)) + (EPS(2,2)-EPSN(2,2))
     1   + (EPS(3,3)-EPSN(3,3))
      TR = TR/3.0D0
      DO I = 1,3
        DO J = 1,3
          DEPDEV(I,J) = EPS(I,J) - EPSN(I,J)
        END DO
        DEPDEV(I,I) = DEPDEV(I,I) - TR
      END DO
      RETURN
      END

      DOUBLE PRECISION FUNCTION XOVSINH(X)
      DOUBLE PRECISION X
      IF (ABS(X) .LT. 1.0D-6) THEN
        XOVSINH = 1.0D0 - X*X/6.0D0
      ELSE
        XOVSINH = X/SINH(X)
      END IF
      RETURN
      END

      DOUBLE PRECISION FUNCTION DVARPI(MU0,OMEGA,XI,VS,R,THETA)
      DOUBLE PRECISION MU0,OMEGA,XI,VS,R,THETA
      IF (ABS(XI) .LT. 1.0D-6) THEN
        DVARPI = MU0*OMEGA*(-XI/3.0D0)*VS/(2.0D0*R*THETA)
      ELSE
        DVARPI = MU0*OMEGA*(1.0D0/SINH(XI))*(1.0D0-XI/TANH(XI))*VS
     1           /(2.0D0*R*THETA)
      END IF
      RETURN
      END

      SUBROUTINE CONFORMATIONAL(DT,EPSTRI,EPSNP,ENN,D,EE,A2,ETA,GAMMA0,
     1                          LAMC,EN,TAUC_P,IERR)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION DT,EPSTRI(3),EPSNP(3),ENN(3),D,EE,A2,ETA,GAMMA0,
     1  LAMC,EN(3),TAUC_P(3)
      INTEGER IT,I,J,K,IERR
      DOUBLE PRECISION LAM(3),L2(3),TAU(3),R(3),KMAT(3,3),Q(3,3)
      DOUBLE PRECISION DTILDE(3,3),XI,ZETA,TMP,LAMMAX,GAMMA,DG_DLMAX
      DOUBLE PRECISION T(3),DXI_DLAM(3),DZETA_DLAM(3)
      DOUBLE PRECISION DTAU_DXI(3),DTAU_DZETA(3),CHI
      DOUBLE PRECISION SOL(3,3),DEN(3),RMAX
      IERR = 0
      DO I = 1,3
        DO J = 1,3
          DTILDE(I,J) = -1.0D0/3.0D0
        END DO
        DTILDE(I,I) = 2.0D0/3.0D0
      END DO
      DO I = 1,3
        EN(I) = EPSTRI(I)
      END DO
      DO IT = 1,80
        DO I = 1,3
          LAM(I) = EXP(EN(I))
          L2(I) = LAM(I)*LAM(I)
        END DO
        CALL TAU_C_PRINCIPAL(LAM,D,EE,A2,ETA,TAU)
        LAMMAX = MAX(LAM(1),LAM(2),LAM(3))
        IF (LAMMAX .LT. LAMC) THEN
          GAMMA = GAMMA0*(LAMC-1.0D0)/(LAMC-LAMMAX)
          DG_DLMAX = GAMMA0*(LAMC-1.0D0)/(LAMC-LAMMAX)**2
        ELSE
          GAMMA = 1.0D30
          DG_DLMAX = 0.0D0
        END IF
        DO I = 1,3
          R(I) = EN(I) - EPSTRI(I)
          DO J = 1,3
            R(I) = R(I) + (DT/GAMMA)*DTILDE(I,J)*TAU(J)
          END DO
        END DO
        CALL CONFORMATIONAL_Q(LAM,D,EE,A2,ETA,Q)
        DO I = 1,3
          T(I) = 0.0D0
          IF (ABS(LAM(I)-LAMMAX) .LT. 1.0D-12) T(I) = 1.0D0
        END DO
        DO I = 1,3
          DO J = 1,3
            KMAT(I,J) = DTILDE(I,J)
          END DO
          KMAT(I,I) = KMAT(I,I) + 1.0D0
        END DO
        CALL MATMUL3(DTILDE,Q,SOL)
        DO I = 1,3
          DO J = 1,3
            KMAT(I,J) = KMAT(I,J) + (DT/GAMMA)*SOL(I,J)
          END DO
        END DO
        TMP = DG_DLMAX*LAMMAX*DT/(GAMMA*GAMMA)
        DO I = 1,3
          DEN(I) = 0.0D0
          DO J = 1,3
            DEN(I) = DEN(I) + DTILDE(I,J)*TAU(J)
          END DO
        END DO
        DO I = 1,3
          DO J = 1,3
            KMAT(I,J) = KMAT(I,J) - TMP*DEN(I)*T(J)
          END DO
        END DO
        CALL SOLVE3(KMAT,R,DEN)
        DO I = 1,3
          EN(I) = EN(I) - DEN(I)
        END DO
        RMAX = MAX(ABS(R(1)),ABS(R(2)),ABS(R(3)))
        IF (RMAX .LT. 1.0D-10*MAX(1.0D0,ABS(EN(1)),
     1      ABS(EN(2)),ABS(EN(3)))) GOTO 210
      END DO
      IERR = 4
      RETURN
  210 CONTINUE
      DO I = 1,3
        LAM(I) = EXP(EN(I))
      END DO
      CALL TAU_C_PRINCIPAL(LAM,D,EE,A2,ETA,TAUC_P)
      RETURN
      END

      SUBROUTINE TAU_C_PRINCIPAL(LAM,D,EE,A2,ETA,TAU)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION LAM(3),D,EE,A2,ETA,TAU(3)
      INTEGER I
      DOUBLE PRECISION XI,ZETA,L2
      XI = 0.0D0
      ZETA = 0.0D0
      DO I = 1,3
        L2 = LAM(I)*LAM(I)
        XI = XI + L2/(1.0D0+ETA*L2)
        ZETA = ZETA + L2
      END DO
      ZETA = 1.0D0 - A2*ZETA
      DO I = 1,3
        L2 = LAM(I)*LAM(I)
        TAU(I) = 2.0D0*D*A2*L2*(EE*XI/(ZETA*ZETA) - 1.0D0/ZETA)
     1         + D*EE*2.0D0*L2/(ZETA*(1.0D0+ETA*L2)**2)
     2         + D*2.0D0*ETA*L2/(1.0D0+ETA*L2)
      END DO
      RETURN
      END

      SUBROUTINE CONFORMATIONAL_Q(LAM,D,EE,A2,ETA,Q)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION LAM(3),D,EE,A2,ETA,Q(3,3)
      INTEGER A,B
      DOUBLE PRECISION XI,ZETA,L2(3)
      DOUBLE PRECISION DXI_DLAM(3),DZETA_DLAM(3)
      DOUBLE PRECISION DTAU_DXI(3),DTAU_DZETA(3),CHI,DTDL
      DO A = 1,3
        L2(A) = LAM(A)*LAM(A)
      END DO
      XI = 0.0D0
      ZETA = 0.0D0
      DO A = 1,3
        XI = XI + L2(A)/(1.0D0+ETA*L2(A))
        ZETA = ZETA + L2(A)
      END DO
      ZETA = 1.0D0 - A2*ZETA
      DO A = 1,3
        DXI_DLAM(A)   = 2.0D0*LAM(A)/(1.0D0+ETA*L2(A))**2
        DZETA_DLAM(A) = -2.0D0*A2*LAM(A)
        DTAU_DXI(A)   = 2.0D0*D*EE*A2*L2(A)/(ZETA*ZETA)
        DTAU_DZETA(A) = 2.0D0*D*A2*L2(A)*(1.0D0/(ZETA*ZETA)
     1      - 2.0D0*EE*XI/(ZETA**3))
     2      - D*EE*2.0D0*L2(A)/(ZETA*ZETA*(1.0D0+ETA*L2(A))**2)
      END DO
      DO B = 1,3
        DO A = 1,3
          DTDL = DTAU_DXI(A)*DXI_DLAM(B)
     1         + DTAU_DZETA(A)*DZETA_DLAM(B)
          IF (A .EQ. B) THEN
            CHI = 4.0D0*D*A2*LAM(A)*(EE*XI/(ZETA*ZETA)
     1            -1.0D0/ZETA)
     2          + D*EE/ZETA*4.0D0*LAM(A)*(1.0D0-ETA*L2(A))
     3            /(1.0D0+ETA*L2(A))**3
     4          + D*4.0D0*ETA*LAM(A)/(1.0D0+ETA*L2(A))**2
            DTDL = DTDL + CHI
          END IF
          Q(A,B) = DTDL*LAM(B)
        END DO
      END DO
      RETURN
      END

      SUBROUTINE MATMUL3(A,B,C)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION A(3,3),B(3,3),C(3,3)
      INTEGER I,J,K
      DO I = 1,3
        DO J = 1,3
          C(I,J) = 0.0D0
          DO K = 1,3
            C(I,J) = C(I,J) + A(I,K)*B(K,J)
          END DO
        END DO
      END DO
      RETURN
      END

      SUBROUTINE SOLVE3(A,B,X)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION A(3,3),B(3),X(3)
      DOUBLE PRECISION DET,T1,T2,T3
      DET = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
     1    - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
     2    + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      IF (ABS(DET) .LT. 1.0D-30) DET = 1.0D-30
      T1 = B(1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
     1   - A(1,2)*(B(2)*A(3,3)-A(2,3)*B(3))
     2   + A(1,3)*(B(2)*A(3,2)-A(2,2)*B(3))
      T2 = A(1,1)*(B(2)*A(3,3)-A(2,3)*B(3))
     1   - B(1)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
     2   + A(1,3)*(A(2,1)*B(3)-B(2)*A(3,1))
      T3 = A(1,1)*(A(2,2)*B(3)-B(2)*A(3,2))
     1   - A(1,2)*(A(2,1)*B(3)-B(2)*A(3,1))
     2   + B(1)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      X(1) = T1/DET
      X(2) = T2/DET
      X(3) = T3/DET
      RETURN
      END

      SUBROUTINE JACOBI(A,LAM,V,IERR)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION A(3,3),LAM(3),V(3,3),S(3,3),TMP(3,3)
      INTEGER I,J,P,Q,ITER,IERR
      DOUBLE PRECISION THETA,T,C,SS,APP,AQQ,APQ,AMAX
      IERR = 0
      DO I = 1,3
        DO J = 1,3
          S(I,J) = A(I,J)
          V(I,J) = 0.0D0
        END DO
        V(I,I) = 1.0D0
      END DO
      DO ITER = 1,60
        P = 1
        Q = 2
        AMAX = ABS(S(1,2))
        IF (ABS(S(1,3)) .GT. AMAX) THEN
          AMAX = ABS(S(1,3))
          P = 1
          Q = 3
        END IF
        IF (ABS(S(2,3)) .GT. AMAX) THEN
          AMAX = ABS(S(2,3))
          P = 2
          Q = 3
        END IF
        IF (AMAX .LT. 1.0D-14) GOTO 410
        APP = S(P,P)
        AQQ = S(Q,Q)
        APQ = S(P,Q)
        THETA = 0.5D0*(AQQ-APP)/APQ
        T = SIGN(1.0D0,THETA)/(ABS(THETA)+SQRT(THETA*THETA+1.0D0))
        C = 1.0D0/SQRT(T*T+1.0D0)
        SS = T*C
        DO I = 1,3
          IF (I .NE. P .AND. I .NE. Q) THEN
            TMP(I,P) = C*S(I,P) - SS*S(I,Q)
            TMP(I,Q) = SS*S(I,P) + C*S(I,Q)
          END IF
        END DO
        DO I = 1,3
          IF (I .NE. P .AND. I .NE. Q) THEN
            S(I,P) = TMP(I,P)
            S(I,Q) = TMP(I,Q)
            S(P,I) = S(I,P)
            S(Q,I) = S(I,Q)
          END IF
        END DO
        S(P,P) = C*C*APP - 2.0D0*SS*C*APQ + SS*SS*AQQ
        S(Q,Q) = SS*SS*APP + 2.0D0*SS*C*APQ + C*C*AQQ
        S(P,Q) = 0.0D0
        S(Q,P) = 0.0D0
        DO I = 1,3
          TMP(I,P) = C*V(I,P) - SS*V(I,Q)
          TMP(I,Q) = SS*V(I,P) + C*V(I,Q)
        END DO
        DO I = 1,3
          V(I,P) = TMP(I,P)
          V(I,Q) = TMP(I,Q)
        END DO
      END DO
      IERR = 3
      RETURN
  410 CONTINUE
      LAM(1) = S(1,1)
      LAM(2) = S(2,2)
      LAM(3) = S(3,3)
      CALL SORT3(LAM,V)
      RETURN
      END

      SUBROUTINE SORT3(LAM,V)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION LAM(3),V(3,3),TMP,VTMP(3)
      INTEGER I,J
      DO I = 1,2
        DO J = I+1,3
          IF (LAM(J) .GT. LAM(I)) THEN
            TMP = LAM(I)
            LAM(I) = LAM(J)
            LAM(J) = TMP
            VTMP(1) = V(1,I)
            VTMP(2) = V(2,I)
            VTMP(3) = V(3,I)
            V(1,I) = V(1,J)
            V(2,I) = V(2,J)
            V(3,I) = V(3,J)
            V(1,J) = VTMP(1)
            V(2,J) = VTMP(2)
            V(3,J) = VTMP(3)
          END IF
        END DO
      END DO
      RETURN
      END
