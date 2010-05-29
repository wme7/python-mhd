C     ------------
CN    NAME: R I E M A N N _ V T 
C     ------------

CP    PURPOSE:
CP    THIS PROGRAM COMPUTES THE SOLUTION OF A 1D   
CP    RELATIVISTIC RIEMANN PROBLEM WITH ARBITRARY TANGENTIAL VELOCITIES
CP    (FOR CONSTANT-GAMMA IDEAL GASES) 
CP    WITH INITIAL DATA UL IF X<R0 AND UR IF X>R0   
CP    IN THE WHOLE SPATIAL DOMAIN [R0 - 0.5,R0 + 0.5]  
C

CC    COMMENTS:
CC    SEE PONS, MARTI AND MUELLER, JFM, 2000
CC
CC    WRITTEN BY:     Jose-Maria Marti
CC                    Departamento de Astronomia y Astrofisica 
CC                    Universidad de Valencia 
CC                    46100 Burjassot (Valencia), Spain
CC                    jose-maria.marti@uv.es
CC    AND
CC                    Ewald Mueller
CC                    Max-Planck-Institut fuer Astrophysik
CC                    Karl-Schwarzschild-Str. 1
CC                    85741 Garching, Germany
CC                    emueller@mpa-garching.mpg.de
C

      SUBROUTINE RIEMANN_VT( IT, IGB, IPL, IRHOL, IVELL, IVELTL,
     &                                IPR, IRHOR, IVELR, IVELTR )

      IMPLICIT NONE

      INCLUDE 'npoints'

C     ----- 
C     ARGUMENTS 
C     -----

      DOUBLE PRECISION IT, IGB
      DOUBLE PRECISION IPL, IRHOL, IVELL, IVELTL
      DOUBLE PRECISION IPR, IRHOR, IVELR, IVELTR

C     ------
C     COMMON BLOCKS
C     ------

      DOUBLE PRECISION RHOL, PL, UL, HL, CSL, VELL, VELTL, WL, 
     &                 RHOR, PR, UR, HR, CSR, VELR, VELTR, WR 
      COMMON /STATES/  RHOL, PL, UL, HL, CSL, VELL, VELTL, WL, 
     &                 RHOR, PR, UR, HR, CSR, VELR, VELTR, WR

      DOUBLE PRECISION RHOLS, ULS, HLS, CSLS, VELLS, VELTLS, VSHOCKL 
      COMMON /LS/      RHOLS, ULS, HLS, CSLS, VELLS, VELTLS, VSHOCKL

      DOUBLE PRECISION RHORS, URS, HRS, CSRS, VELRS, VELTRS, VSHOCKR 
      COMMON /RS/      RHORS, URS, HRS, CSRS, VELRS, VELTRS, VSHOCKR

      DOUBLE PRECISION GB
      COMMON /GB/      GB

      DOUBLE PRECISION QX(NG), QW(NG)
      COMMON /INT/     QX, QW
 
C     -----------
C     INTERNAL VARIABLES
C     -----------

      INTEGER          MN, N, I, ILOOP 
      PARAMETER        (MN = 400)

      DOUBLE PRECISION TOL, PMIN, PMAX, DVEL1, DVEL2, CHECK, R0
      PARAMETER        (R0=0.D0)
 
      DOUBLE PRECISION PS, VELS

      DOUBLE PRECISION RHOA(MN), PA(MN), VELA(MN), VELTA(MN),
     &                 UA(MN), HA(MN)
  
      DOUBLE PRECISION A
 
      DOUBLE PRECISION RAD(MN), X1, X2, X3, X4, X5, T, LAM, LAMS

      DOUBLE PRECISION KL, KR, CONST
 
C     ------- 
C     INITIAL STATES 
C     -------

      GB    = IGB
      T     = IT

C     -----
C     LEFT STATE
C     -----
      PL    = IPL
      RHOL  = IRHOL
      VELL  = IVELL
      VELTL = IVELTL

C     -----
C     RIGHT STATE
C     -----
      PR    = IPR
      RHOR  = IRHOR
      VELR  = IVELR
      VELTR = IVELTR

C     ---------------------------------------
C     SPECIFIC INTERNAL ENERGY, SPECIFIC ENTHALPY, SOUND SPEED, 
C     ADIABATIC CONSTANT AND FLOW LORENTZ FACTORS IN THE INITIAL STATES 
C     ---------------------------------------

      UL = PL/(GB - 1.D0)/RHOL
      UR = PR/(GB - 1.D0)/RHOR

      HL = 1.D0 + GB*UL
      HR = 1.D0 + GB*UR

      CSL= DSQRT((GB - 1.D0)*(HL - 1.D0)/HL)
      CSR= DSQRT((GB - 1.D0)*(HR - 1.D0)/HR)

      KL = PL/RHOL**GB
      KR = PR/RHOR**GB

      WL = 1.D0/DSQRT(1.D0 - VELL*VELL - VELTL*VELTL)
      WR = 1.D0/DSQRT(1.D0 - VELR*VELR - VELTR*VELTR)

C     ------
C     COEFFICIENTS FOR NUMERICAL INTEGRATION IN RAREFACTIONS
C     ------

      CALL GAULEG(-1.D0,1.D0,QX,QW,NG)

C     -------- 
C     NUMBER OF POINTS 
C     --------

      N   = 400

C     ------------- 
C     TOLERANCE FOR THE SOLUTION 
C     -------------

      TOL = 0.D0

C

      ILOOP = 0

      IF ((PL.EQ.PR).AND.(VELL.EQ.VELR)) THEN
 
         PS      = PL
         VELS    = VELL

         VSHOCKL = VELL
         RHOLS   = (PS/KL)**(1.D0/GB)
         ULS     = PS/(GB - 1.D0)/RHOLS
         
         VSHOCKR = VELL
         RHORS   = (PS/KR)**(1.D0/GB)
         URS     = PS/(GB - 1.D0)/RHORS

      ELSE

        PMIN  = (PL + PR)/2.D0
        PMAX  = PMIN
 
5       ILOOP = ILOOP + 1

        PMIN  = 0.5D0*MAX(PMIN,0.D0)
        PMAX  = 2.D0*PMAX
 

        CALL GETDVEL2(PMIN, DVEL1)

        CALL GETDVEL2(PMAX, DVEL2)
        
        CHECK = DVEL1*DVEL2
        IF (CHECK.GT.0.D0) GOTO 5


C     --------------------------- 
C     PRESSURE AND FLOW VELOCITY IN THE INTERMEDIATE STATES 
C     ---------------------------

        CALL GETP2(PMIN, PMAX, TOL, PS)

        VELS = 0.5D0*(VELLS + VELRS)

C        WRITE(*,*) 'VELS = ', VELS, 'PS = ', PS

      ENDIF
      
C     -------
C     SOLUTION ON THE NUMERICAL MESH 
C     -------

C     -----------
C     POSITIONS OF THE WAVES
C     -----------

      IF (PL.GE.PS) THEN

        CONST = HL*WL*VELTL

        CALL FLAMB(KL, CONST, PL, VELL, 'L', LAM )

        CALL FLAMB(KL, CONST, PS, VELS, 'L', LAMS)

        X1 = R0 + LAM *T
        X2 = R0 + LAMS*T 

      ELSE

        X1 = R0 + VSHOCKL*T 
        X2 = X1

      END IF
 
      X3 = R0 + VELS*T 
 
      IF (PR.GE.PS) THEN

        CONST = HR*WR*VELTR

        CALL FLAMB(KR, CONST, PS, VELS, 'R', LAMS)

        CALL FLAMB(KR, CONST, PR, VELR, 'R', LAM )

        X4 = R0 + LAMS*T 
        X5 = R0 + LAM *T
 
      ELSE

        X4 = R0 + VSHOCKR*T 
        X5 = X4

      END IF

C     ---------- 
C     SOLUTION ON THE MESH 
C     ----------

      DO 10 I=1,N

        RAD(I) = R0 + DFLOAT(I)/DFLOAT(N) - 0.5D0

 10     CONTINUE
 
      DO 120 I=1,N

        IF (RAD(I).LE.X1) THEN

          PA(I)    = PL
          RHOA(I)  = RHOL
          VELA(I)  = VELL
          VELTA(I) = VELTL  
          UA(I)    = UL
          HA(I)    = 1.D0 + GB*UA(I)

        ELSE IF (RAD(I).LE.X2) THEN

          A = (RAD(I) - R0)/T 

          CALL RAREF2(A, PS, RHOL, PL, UL, CSL, VELL, VELTL,  
     &                'L', RHOA(I), PA(I), UA(I), VELA(I), VELTA(I))

        ELSE IF (RAD(I).LE.X3) THEN

          PA(I)    = PS
          RHOA(I)  = RHOLS
          VELA(I)  = VELS
          UA(I)    = ULS
          HA(I)    = 1.D0 + GB*UA(I)
          CONST    = HL*WL*VELTL
          VELTA(I) = CONST*DSQRT((1.D0 - VELS*VELS)/
     &               (CONST**2 + HA(I)**2))

        ELSE IF (RAD(I).LE.X4) THEN

          PA(I)    = PS
          RHOA(I)  = RHORS
          VELA(I)  = VELS
          UA(I)    = URS
          HA(I)    = 1.D0 + GB*UA(I)
          CONST    = HR*WR*VELTR
          VELTA(I) = CONST*DSQRT((1.D0 - VELS*VELS)/
     &               (CONST**2 + HA(I)**2))

        ELSE IF (RAD(I).LE.X5) THEN

          A = (RAD(I) - R0)/T 
 
          CALL RAREF2(A, PS, RHOR, PR, UR, CSR, VELR, VELTR,  
     &                'R', RHOA(I), PA(I), UA(I), VELA(I), VELTA(I))

        ELSE

          PA(I)    = PR
          RHOA(I)  = RHOR
          VELA(I)  = VELR
          VELTA(I) = VELTR
          UA(I)    = UR
          HA(I)    = 1.D0 + GB*UA(I)

        END IF

120     CONTINUE

      OPEN (3,FILE='solution.dat',FORM='FORMATTED',STATUS='NEW')

C      WRITE(3,150) N, T 
 150  FORMAT(I5,1X,F10.5)

      DO 60 I=1,N
        WRITE(3,200) RAD(I),PA(I),RHOA(I),VELA(I),VELTA(I)
60      CONTINUE

200   FORMAT(5(E14.8,1X))
 
      CLOSE(3)
      END
 
C     ---------- 
CN    NAME: G E T D V E L 2
C     ----------

CP    PURPOSE: 
CP    COMPUTE THE DIFFERENCE IN FLOW SPEED BETWEEN LEFT AND RIGHT INTERMEDIATE 
CP    STATES FOR GIVEN LEFT AND RIGHT STATES AND PRESSURE 
C 

CC    COMMENTS:
CC    NONE
 
      SUBROUTINE GETDVEL2( P, DVEL )

      IMPLICIT NONE

      INCLUDE 'npoints'

C     ------
C     ARGUMENTS
C     ------
 
      DOUBLE PRECISION P, DVEL

C     ------
C     COMMON BLOCKS
C     ------

      DOUBLE PRECISION RHOL, PL, UL, HL, CSL, VELL, VELTL, WL, 
     &                 RHOR, PR, UR, HR, CSR, VELR, VELTR, WR 
      COMMON /STATES/  RHOL, PL, UL, HL, CSL, VELL, VELTL, WL, 
     &                 RHOR, PR, UR, HR, CSR, VELR, VELTR, WR
 
      DOUBLE PRECISION RHOLS, ULS, HLS, CSLS, VELLS, VELTLS, VSHOCKL 
      COMMON /LS/      RHOLS, ULS, HLS, CSLS, VELLS, VELTLS, VSHOCKL

      DOUBLE PRECISION RHORS, URS, HRS, CSRS, VELRS, VELTRS, VSHOCKR 
      COMMON /RS/      RHORS, URS, HRS, CSRS, VELRS, VELTRS, VSHOCKR

      DOUBLE PRECISION GB
      COMMON /GB/      GB

      DOUBLE PRECISION QX(NG), QW(NG)
      COMMON /INT/     QX, QW

C     -----
C     LEFT WAVE 
C     -----

      CALL GETVEL2(P, RHOL, PL, UL,  HL,  CSL,  VELL,  VELTL,  WL, 'L',
     &                RHOLS,    ULS, HLS, CSLS, VELLS, VELTLS, VSHOCKL)

C     -----
C     RIGHT WAVE
C     -----

      CALL GETVEL2(P, RHOR, PR, UR,  HR,  CSR,  VELR,  VELTR,  WR, 'R',
     &                RHORS,    URS, HRS, CSRS, VELRS, VELTRS, VSHOCKR)

      DVEL = VELLS - VELRS

      RETURN
      END
C     ------- 
CN    NAME: G E T P 2
C     -------

CP    PURPOSE: 
CP    FIND THE PRESSURE IN THE INTERMEDIATE STATE OF A RIEMANN PROBLEM IN 
CP    RELATIVISTIC HYDRODYNAMICS 
C 

CC    COMMENTS: 
CC    THIS ROUTINE USES A COMBINATION OF INTERVAL BISECTION AND INVERSE 
CC    QUADRATIC INTERPOLATION TO FIND THE ROOT IN A SPECIFIED INTERVAL. 
CC    IT IS ASSUMED THAT DVEL(PMIN) AND DVEL(PMAX) HAVE OPPOSITE SIGNS WITHOUT 
CC    A CHECK. 
CC    ADAPTED FROM "COMPUTER METHODS FOR MATHEMATICAL COMPUTATION", 
CC    BY G. E. FORSYTHE, M. A. MALCOLM, AND C. B. MOLER, 
CC    PRENTICE-HALL, ENGLEWOOD CLIFFS N.J. 

      SUBROUTINE GETP2( PMIN, PMAX, TOL, PS )

      IMPLICIT NONE

C     ----- 
C     ARGUMENTS 
C     -----

      DOUBLEPRECISION PMIN, PMAX, TOL, PS

C     ------- 
C     COMMON BLOCKS 
C     -------

      DOUBLE PRECISION GB
      COMMON /GB/      GB

      DOUBLE PRECISION RHOL, PL, UL, HL, CSL, VELL, VELTL, WL, 
     &                 RHOR, PR, UR, HR, CSR, VELR, VELTR, WR 
      COMMON /STATES/  RHOL, PL, UL, HL, CSL, VELL, VELTL, WL, 
     &                 RHOR, PR, UR, HR, CSR, VELR, VELTR, WR

C     --------- 
C     INTERNAL VARIABLES 
C     ---------

      DOUBLEPRECISION A, B, C, D, E, EPS, FA, FB, FC, TOL1, 
     &                XM, P, Q, R, S

C     ------------- 
C     COMPUTE MACHINE PRECISION 
C     -------------

      EPS  = 1.D0 
10    EPS  = EPS/2.D0 
      TOL1 = 1.D0 + EPS 
      IF( TOL1 .GT. 1.D0 ) GO TO 10

C     ------- 
C     INITIALIZATION 
C     -------

      A  = PMIN 
      B  = PMAX 
      CALL GETDVEL2(A,FA) 
      CALL GETDVEL2(B,FB)

C     ----- 
C     BEGIN STEP 
C     -----

20    C  = A 
      FC = FA 
      D  = B - A 
      E  = D 
30    IF( DABS(FC) .GE. DABS(FB) )GO TO 40 
      A  = B 
      B  = C 
      C  = A 
      FA = FB 
      FB = FC 
      FC = FA

C     -------- 
C     CONVERGENCE TEST 
C     --------

40    TOL1 = 2.D0*EPS*DABS(B) + 0.5D0*TOL 
      XM   = 0.5D0*(C - B) 
      IF( DABS(XM) .LE. TOL1 ) GO TO 90 
      IF( FB .EQ. 0.D0 ) GO TO 90

C     ------------ 
C     IS BISECTION NECESSARY? 
C     ------------

      IF( DABS(E) .LT. TOL1 ) GO TO 70 
      IF( DABS(FA) .LE. DABS(FB) ) GO TO 70

C     ------------------ 
C     IS QUADRATIC INTERPOLATION POSSIBLE? 
C     ------------------

      IF( A .NE. C ) GO TO 50

C     ---------- 
C     LINEAR INTERPOLATION 
C     ----------

      S = FB/FA 
      P = 2.D0*XM*S 
      Q = 1.D0 - S 
      GO TO 60

C     ---------------- 
C     INVERSE QUADRATIC INTERPOLATION 
C     ----------------

50    Q = FA/FC 
      R = FB/FC 
      S = FB/FA 
      P = S*(2.D0*XM*Q*(Q - R) - (B - A)*(R - 1.D0)) 
      Q = (Q - 1.D0)*(R - 1.D0)*(S - 1.D0)

C     ------ 
C     ADJUST SIGNS 
C     ------

60    IF( P .GT. 0.D0 ) Q = -Q 
      P = DABS(P)

C     -------------- 
C     IS INTERPOLATION ACCEPTABLE? 
C     --------------

      IF( (2.D0*P) .GE. (3.D0*XM*Q-DABS(TOL1*Q)) ) GO TO 70 
      IF( P .GE. DABS(0.5D0*E*Q) ) GO TO 70 
      E = D 
      D = P/Q 
      GO TO 80

C     ----- 
C     BISECTION 
C     -----

70    D = XM 
      E = D

C     ------- 
C     COMPLETE STEP 
C     -------

80    A  = B 
      FA = FB 
      IF( DABS(D) .GT. TOL1 ) B = B+D 
      IF( DABS(D) .LE. TOL1 ) B = B+DSIGN(TOL1,XM) 
      CALL GETDVEL2(B,FB) 
      IF( (FB*(FC/DABS(FC))) .GT. 0.D0) GO TO 20 
      GO TO 30

C     -- 
C     DONE 
C     --

90    PS = B

      RETURN 
      END

C     ------
CN    NAME: F L A M B
C     ------

CP    PURPOSE:
CP    COMPUTE THE VALUE OF THE SELF-SIMILARITY VARIABLE INSIDE A RAREFACTION 
CP    CONNECTED TO A SPECIFIED LEFT / RIGHT STATE
C

CC    COMMENTS:
CC    NONE

      SUBROUTINE FLAMB(K, A, P, VEL, S, XI)

      IMPLICIT NONE

C     --------
C     ARGUMENTS
C     --------

      DOUBLE PRECISION K, A, P, VEL

      CHARACTER*1      S

      DOUBLE PRECISION XI

C     -------
C     COMMON BLOCKS
C     -------

      DOUBLE PRECISION G
      COMMON /GB/      G

C     --------------
C     INTERNAL VARIABLES
C     --------------
      
      DOUBLE PRECISION SIGN
      DOUBLE PRECISION RHO, H, CS2, VELT2, V2, BETA, DISC

      IF (S.EQ.'L') SIGN = -1.D0

      IF (S.EQ.'R') SIGN =  1.D0
  
      RHO   = (P/K)**(1.D0/G)
      CS2   = G*(G - 1.D0)*P/(G*P +(G - 1.D0)*RHO)
      H     = 1.D0/(1.D0 - CS2/(G - 1.D0))

      VELT2 = (1.D0 - VEL*VEL)*A*A/(H*H + A*A)
      V2    = VELT2 + VEL*VEL
     
      BETA  = (1.D0 - V2)*CS2/(1.D0 - CS2)
      DISC  = DSQRT(BETA*(1.D0 + BETA - VEL*VEL))
      
      XI = (VEL + SIGN*DISC)/(1.D0 + BETA)

      RETURN
      END

C     -------- 
CN    NAME: R A R E F 2
C     --------

CP    PURPOSE: 
CP    COMPUTE THE FLOW STATE IN A RAREFACTION FOR GIVEN PRE-WAVE STATE 
C
 
CC    COMMENTS: 
CC    THE VELOCITY IN THE RAREFACTION IS WRITTEN IN TERMS OF THE PRESCRIBED 
CC    LEFT / RIGHT STATE  AND PRESSURE ACCORDING TO EXPRESSIONS (3.25) AND 
CC    (3.26) OF REZZOLLA, ZANOTTI AND PONS, JFM, 2002.
CC    THE INTEGRAL IN THE VELOCITY EXPRESSION IS COMPUTED THROUGH A GAUSSIAN 
CC    QUADRATURE 

      SUBROUTINE RAREF2(XI, PS, RHOA, PA, UA, CSA, VELA, VELTA, S, 
     &                      RHO,  P,  U,       VEL,  VELT)

      IMPLICIT NONE

      INCLUDE 'npoints'

C     ------
C     ARGUMENTS
C     ------

      DOUBLE PRECISION XI, PS, RHOA, PA, UA, CSA, VELA, VELTA

      CHARACTER*1      S

      DOUBLE PRECISION RHO, P, U, VEL, VELT

C     ------
C     COMMON BLOCKS
C     ------

      DOUBLE PRECISION GB
      COMMON /GB/      GB

      DOUBLE PRECISION QX(NG), QW(NG)
      COMMON /INT/     QX, QW

C     --------
C     INTERNAL VARIABLES
C     --------

      INTEGER          I

      DOUBLE PRECISION HA, WA, SIGN

      DOUBLE PRECISION CONST, K, XIO, XIP, PO, H

      DOUBLE PRECISION SUMW, DIFW, INTEGRAL, XX, RRHO, CCS2, HH, 
     &                 FUNR, A, FP, DFDP

      HA = 1.D0 + UA + PA/RHOA

      WA = 1.D0/DSQRT(1.D0 - VELA*VELA - VELTA*VELTA)

      CONST = HA*WA*VELTA

      K     = PA/RHOA**GB

      IF (S.EQ.'L') SIGN = -1.D0

      IF (S.EQ.'R') SIGN =  1.D0

      CALL FLAMB(K, CONST, PA, VELA, S, XIO)

      PO = PA
      P  = 0.95D0*PA

 20   CONTINUE

      SUMW = 0.5D0*(P + PA)
      DIFW = 0.5D0*(P - PA)
      INTEGRAL = 0.D0

      DO 10 I = 1, NG
         
         XX   = DIFW*QX(I) + SUMW
         RRHO = (XX/K)**(1.D0/GB)
         CCS2 = GB*(GB - 1.D0)*XX/(GB*XX + (GB - 1.0)*RRHO) 
         HH   = 1.D0/(1.D0 - CCS2/(GB - 1.D0)) 

         FUNR = DSQRT(HH*HH + CONST*CONST*(1.D0 - CCS2))/
     &         (HH*HH + CONST*CONST)/(RRHO*DSQRT(CCS2))
 
         INTEGRAL = INTEGRAL + DIFW*QW(I)*FUNR
        
 10      CONTINUE

      A    = SIGN*INTEGRAL + 0.5D0*DLOG((1.D0 + VELA)/(1.D0 - VELA))
      VEL  = DTANH(A)
 
      CALL FLAMB(K, CONST, P, VEL, S, XIP)

      FP   = XIP - XI

      DFDP = (XIP - XIO)/(P - PO)

      PO  = P
      XIO = XIP 

      P = P - FP/DFDP
      P = DMAX1(P,PS)

      IF (DABS(FP).GT.1.D-10) GOTO 20

      RHO = (P/K)**(1.D0/GB)

      U   = P/RHO/(GB -1.D0)
      
      H   = 1.D0 + U + P/RHO

      VELT = CONST*DSQRT((1.D0 - VEL*VEL)/(CONST*CONST + H*H))

      RETURN
      END

C     --------- 
CN    NAME: G E T V E L 2
C     ---------

CP    PURPOSE: 
CP    COMPUTE THE FLOW VELOCITY BEHIND A RAREFACTION OR SHOCK IN TERMS OF THE 
CP    POST-WAVE PRESSURE FOR A GIVEN STATE AHEAD THE WAVE IN A RELATIVISTIC 
CP    FLOW 
C 

CC    COMMENTS: 
CC    THIS ROUTINE CLOSELY FOLLOWS THE EXPRESSIONS IN PONS, MARTI AND MUELLER, 
CC    JFM, 2002

      SUBROUTINE GETVEL2(P, RHOA, PA, UA, HA, CSA, VELA, VELTA, WA, S,
     &                      RHO,      U,  H,  CS,  VEL,  VELT,  VSHOCK)

      IMPLICIT NONE

      INCLUDE 'npoints'

C     --------
C     ARGUMENTS
C     --------

      DOUBLE PRECISION P, RHOA, PA, UA, HA, CSA, VELA, VELTA, WA
      CHARACTER*1      S
      DOUBLE PRECISION RHO, U, H, CS, VEL, VELT, VSHOCK

C     -----
C     COMMON BLOCKS
C     -----

      DOUBLE PRECISION GB
      COMMON /GB/      GB

      DOUBLE PRECISION QX(NG), QW(NG)
      COMMON /INT/     QX, QW

C     --------
C     INTERNAL VARIABLES
C     --------

      INTEGER          I

      DOUBLE PRECISION A, B, C, SIGN
      DOUBLE PRECISION J, WSHOCK
      DOUBLE PRECISION K, CONST, SUMW, DIFW, INTEGRAL, XX, RRHO,
     &                 HH, CCS2, FUNR

C     ------------
C     LEFT OR RIGHT PROPAGATING WAVE
C     ------------

      IF (S.EQ.'L') SIGN = -1.D0

      IF (S.EQ.'R') SIGN =  1.D0

      IF (P.GE.PA) THEN

C       ---
C       SHOCK 
C       ---

        A  = 1.D0 - (GB - 1.D0)*(P - PA)/GB/P
        B  = 1.D0 - A
        C  = HA*(PA - P)/RHOA - HA**2
 
C       ---------------- 
C       CHECK FOR UNPHYSICAL ENTHALPIES 
C       ----------------

        IF (C.GT.(B**2/4.D0/A)) STOP
     & 'GETVEL2: UNPHYSICAL SPECIFIC ENTHALPY IN INTERMEDIATE STATE'
 
C       ---------------------------------
C       SPECIFIC ENTHALPY AT THE LEFT OF THE CONTACT DISCONTINUITY 
C       (OBTAINED FROM THE EQUATION OF STATE AND THE TAUB ADIABAT)
C       ---------------------------------

        H = (-B + DSQRT(B**2 - 4.D0*A*C))/2.D0/A
 
C       -----------------------------------
C       DENSITY AT THE LEFT OF THE CONTACT DISCONTINUITY
C       (OBTAINED FROM SPECIFIC ENTHALPY AND THE EQUATION OF STATE)
C       -----------------------------------

        RHO = GB*P/(GB - 1.D0)/(H - 1.D0)

C       ----------------------------------
C       SPECIFIC INT. ENERGY AT THE LEFT OF THE CONTACT DISCONTINUITY 
C       (OBTAINED FROM THE EQUATION OF STATE)
C       ----------------------------------

        U = P/(GB - 1.D0)/RHO

C       -----------------------------
C       MASS FLUX ACROSS LEFT WAVE
C       (OBTAINED FROM THE RANKINE-HUGONIOT RELATIONS)
C       -----------------------------

        J = SIGN*DSQRT((P - PA)/(HA/RHOA - H/RHO))
 
C       ------------------------------
C       SHOCK VELOCITY 
C       (OBTAINED FROM THE DEFINITION OF MASS FLUX)
C       ------------------------------

        A      = J**2 + (RHOA*WA)**2
        B      = -2.D0*VELA*RHOA**2*WA**2
        C      = (RHOA*WA*VELA)**2 - J**2

        VSHOCK = (-B + SIGN*DSQRT(B*B - 4.D0*A*C))/(2.D0*A)
        WSHOCK = 1.D0/DSQRT(1.D0 - VSHOCK**2)
 
C       --------------------------
C       VELOCITY AT THE LEFT OF THE CONTACT DISCONTINUITY 
C       --------------------------

        A    = WSHOCK*(P - PA)/J + HA*WA*VELA
        B    = HA*WA + (P - PA)*(WSHOCK*VELA/J + 1.D0/RHOA/WA)
 
        VEL  = A/B
 
        A    = HA*WA*VELTA

        VELT = A*DSQRT((1.D0 - VEL*VEL)/(H*H + A*A))
 
C       -----------------------------
C       SOUND SPEED AT THE LEFT OF THE CONTACT DISCONTINUITY
C       -----------------------------

        CS = DSQRT(GB*P/RHO/H)
 
      ELSE

C       ------
C       RAREFACTION 
C       ------

        CONST = HA*WA*VELTA

        K     = PA/RHOA**GB

C       --------------------------
C       DENSITY AT THE LEFT SIDE OF THE CONTACT DISCONTINUITY
C       --------------------------

        RHO = (P/K)**(1.D0/GB)
 
C       ----------------------------------
C       SPECIFIC INT. ENERGY AT THE LEFT OF THE CONTACT DISCONTINUITY 
C       (OBTAINED FROM THE EQUATION OF STATE)
C       ----------------------------------

        U = P/(GB - 1.D0)/RHO
        H = 1.D0 + GB*U
 
C       ----------------------------------
C       SOUND SPEED AT THE LEFT OF THE CONTACT DISCONTINUITY 
C       ----------------------------------

        CS = DSQRT(GB*P/(RHO + GB*P/(GB - 1.D0)))
 
C       ----------------------------------
C       VELOCITY AT THE LEFT OF THE CONTACT DISCONTINUITY 
C       ----------------------------------

C       ------
C       INTEGRAL 
C       ------

        SUMW = 0.5D0*(P + PA)

        DIFW = 0.5D0*(P - PA)

        INTEGRAL = 0.D0

        DO 110 I = 1, NG

           XX  = DIFW*QX(I) + SUMW
           RRHO = (XX/K)**(1.D0/GB)
           CCS2 = GB*(GB - 1.D0)*XX/(GB*XX + (GB - 1.0)*RRHO) 
           HH   = 1.D0/(1.D0 - CCS2/(GB - 1.D0)) 

           FUNR = DSQRT(HH*HH + CONST*CONST*(1.D0 - CCS2))/
     &           (HH*HH + CONST*CONST)/(RRHO*DSQRT(CCS2))
 
           INTEGRAL = INTEGRAL + DIFW*QW(I)*FUNR
        
 110       CONTINUE

        A    = SIGN*INTEGRAL + 0.5D0*DLOG((1.D0 + VELA)/(1.D0 - VELA))
        VEL  = DTANH(A)
        VELT = CONST*DSQRT((1.D0 - VEL*VEL)/(CONST**2 + H**2))
        
      END IF
 
      RETURN
      END

C     ------
CN    NAME: G A U L E G 
C     ------

CP    PURPOSE:
CP    COMPUTE ABCISSAS AND WEIGHTS FOR GAUSS-LEGENDRE QUADRATURE INTEGRATION
C

CC    COMMENTS:
CC    ADAPTED FROM PRESS ET AL., "NUMERICAL RECIPES", CAMBRIDGE, 1988

      SUBROUTINE GAULEG(X1,X2,X,W,N)

      IMPLICIT NONE

C     --------
C     ARGUMENTS
C     --------

      INTEGER          N
      DOUBLE PRECISION X1,X2,X(N),W(N)

C     ---------
C     INTERNAL VARIABLES
C     ---------

      DOUBLE PRECISION EPS
      PARAMETER        (EPS=3.D-14)
      INTEGER          I, J, M
      DOUBLE PRECISION P1, P2, P3, PP, XL, XM, Z, Z1

      M  = (N + 1)/2
      XM = 0.5D0*(X2 + X1)
      XL = 0.5D0*(X2 - X1)

      DO 12 I = 1, M

        Z = COS(3.141592654D0*(I - .25D0)/(N + .5D0))
1       CONTINUE

        P1 = 1.D0
        P2 = 0.D0
          
        DO 11 J = 1, N

          P3 = P2
          P2 = P1
          P1 = ((2.D0*J - 1.D0)*Z*P2 - (J - 1.D0)*P3)/J
11        CONTINUE

        PP = N*(Z*P1 - P2)/(Z*Z - 1.D0)
        Z1 = Z
        Z  = Z1 - P1/PP
        
        IF (ABS(Z -Z1).GT.EPS) GOTO 1

        X(I)     = XM - XL*Z
        X(N+1-I) = XM + XL*Z
        W(I)     = 2.D0*XL/((1.D0 - Z*Z)*PP*PP)
        W(N+1-I) = W(I)
12      CONTINUE

      RETURN
      END
