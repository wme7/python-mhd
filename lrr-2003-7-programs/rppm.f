C     -------- 
CN    NAME: R P P M
C     --------

CP    PURPOSE: 
CP    SOLVE THE 1D EQUATIONS OF RELATIVISTIC HYDRODYNAMICS IN PLANAR SYMMETRY
CP    FOR AN IDEAL GAS EQUATION OF STATE.
C

CC    COMMENTS: 
CC    THIS PROGRAM IS DESCRIBED IN THE PAPER BY MARTI & MUELLER, JCP, 1996.
CC    IT USES AN EXACT RIEMANN SOLVER (MARTI & MUELLER, JFM, 1994), PPM SPATIAL
CC    RECONSTRUCTION AND AVERAGING IN THE DOMAIN OF DEPENDENCE OF THE 
CC    CELL INTERFACES FOR TIME ADVANCE.
CC    LIGHT SPEED IN CODE UNITS IS EQUAL TO 1.
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

      PROGRAM RPPM

      IMPLICIT NONE

      INCLUDE 'size'

C     ----------
C     COMMON BLOCKS
C     ----------

      INTEGER         BNDMNX,BNDMXX
      COMMON /BOUN/   BNDMNX,BNDMXX

      INTEGER         NEND,NOUT,ITSTP,NX
      COMMON /INPTI/  NEND,NOUT,ITSTP,NX

      INTEGER         NSTEP
      COMMON /NSTEP/  NSTEP

      INTEGER         NOUT1
      COMMON /OUTI/   NOUT1

      DOUBLEPRECISION GAMMA
      COMMON /ADIND/  GAMMA

      DOUBLEPRECISION X(-4:MN5),XL(-4:MN5),XR(-4:MN5),DX(-4:MN5)              
      COMMON /GRD/    X,XL,XR,DX

      DOUBLEPRECISION P(-4:MN5),RHO(-4:MN5),VEL(-4:MN5),W(-4:MN5),
     &                U(-4:MN5),CS(-4:MN5),H(-4:MN5),DPDRH(-4:MN5),
     &                DPDU(-4:MN5),R(-4:MN5),M(-4:MN5),E(-4:MN5)
      COMMON /HYDRO/  P,RHO,VEL,W,U,CS,H,DPDRH,DPDU,R,M,E

      DOUBLEPRECISION TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,SMALLU,
     &                GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2
      COMMON /INPTF/  TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,SMALLU,
     &                GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2

      DOUBLEPRECISION TOUT1
      COMMON /OUTF/   TOUT1

      DOUBLEPRECISION TIME,DT
      COMMON /ZEIT/   TIME,DT

      CHARACTER*7     OUTFIL
      CHARACTER*8     LABEL
      CHARACTER*4     BASENM
      CHARACTER*2     SUFFIX
      COMMON /CHRC/   LABEL,OUTFIL,BASENM,SUFFIX

C     --------------------
C     READ INITIAL PARAMETERS
C     --------------------

      CALL INPUT

C     --------------------
C     CONSTRUCT A NEW MODEL
C     --------------------

      IF (SUFFIX(2:2).NE.'A') THEN

        WRITE(6,2000)
 2000   FORMAT('RPPM: CHECK INPUT FILE SUFFIX')
        STOP

      END IF

      WRITE(6,2100)
 2100 FORMAT('RPPM: CONSTRUCTING NEW INITIAL MODEL')

      CALL GRID

      CALL INIT

      CALL TSTEP

      DT = MIN(DT, DTINI)

      OUTFIL = BASENM//'O'//SUFFIX

      NOUT1 = 0
      TOUT1 = 0.D0

C     -------------
C     START TIME LOOP
C     -------------

      DO 100 NSTEP = 1, NEND

         TIME  = TIME  + DT
         NOUT1 = NOUT1 + 1
         TOUT1 = TOUT1 + DT

         IF (TIME.GT.TMAX) GOTO 200

         CALL BNDRY

         CALL HYDROW

         IF ( (NOUT1.GE.NOUT) .OR. (TOUT1.GE.TOUT) ) CALL PLTOUT

         CALL TSTEP

 100     CONTINUE

 200  CONTINUE

      STOP 'RPPM: NORMAL TERMINATION'
      END

C     -------- 
CN    NAME: I N P U T
C     --------

CP    PURPOSE: 
CP    READS THE INITIAL PARAMETERS FROM FILE inpt.dat
C

CC    COMMENTS: 
CC    SEE INSERTED COMMENTS

      SUBROUTINE INPUT

      IMPLICIT NONE

      INCLUDE 'size'

C     ---------
C     COMMON BLOCKS
C     ---------

      INTEGER         NEND,NOUT,ITSTP,NX
      COMMON /INPTI/  NEND,NOUT,ITSTP,NX

      DOUBLEPRECISION TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,SMALLU,
     &                GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2
      COMMON /INPTF/  TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,SMALLU,
     &                GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2

      CHARACTER*7     OUTFIL
      CHARACTER*8     LABEL,LABEL1
      CHARACTER*4     BASENM
      CHARACTER*2     SUFFIX
      COMMON /CHRC/   LABEL,OUTFIL,BASENM,SUFFIX

C     ---------
C     INTERNAL VARIABLES
C     ---------

      CHARACTER*72 TEXT
      CHARACTER*2  TXT
      CHARACTER*15 TXTXT
      DATA TXTXT /'.............  '/

      OPEN(1,FILE='inpt.dat',FORM='FORMATTED',STATUS='OLD')

      READ (1,*) TEXT
      WRITE(6,*) TEXT
      READ (1,*) TEXT
      WRITE(6,*) TEXT
      READ (1,*) TEXT
      WRITE(6,*) TEXT
      PRINT*, '-------------------------------------------------------'

C     ----------------------------
C     BASENM IS THE ROOT FOR THE OUTPUT, PLOT AND RESTART FILE NAMES (ROOTS
C     RST_, RBW_, RSR_, RBWI STAND FOR SPECIAL TESTS)
C     ----------------------------

      READ (1,*) TXT,LABEL,BASENM
      WRITE(6,*) TXT,LABEL,TXTXT,BASENM
      LABEL1 = 'basenm'
      IF (LABEL.NE.LABEL1) GOTO 10

C     ----------------------------
C     NEND IS THE TOTAL NUMBER OF TIMESTEPS
C     ----------------------------

      READ (1,*) TXT,LABEL,NEND
      WRITE(6,*) TXT,LABEL,TXTXT,NEND
      LABEL1 = 'nend'
      IF (LABEL.NE.LABEL1) GOTO 10

C     ----------------------------
C     THE PROGRAM STOPS WHEN TIME IS .GE. THAN TMAX
C     ----------------------------
 
      READ (1,*) TXT,LABEL,TMAX
      WRITE(6,*) TXT,LABEL,TXTXT,TMAX
      LABEL1 = 'tmax'
      IF (LABEL.NE.LABEL1) GOTO 10

C     ----------------------------
C     SUFFIX IS THE SUFFIX FOR THE OUTPUT AND RESTART FILE NAMES
C     ----------------------------

      READ (1,*) TXT,LABEL,SUFFIX
      WRITE(6,*) TXT,LABEL,TXTXT,SUFFIX
      LABEL1 = 'suffix'
      IF (LABEL.NE.LABEL1) GOTO 10

C     ----------------------------
C     AN OUTPUT FILE IS WRITEN EVERY NOUT TIMESTEPS
C     ----------------------------

      READ (1,*) TXT,LABEL,NOUT
      WRITE(6,*) TXT,LABEL,TXTXT,NOUT
      LABEL1 = 'nout'
      IF (LABEL.NE.LABEL1) GOTO 10

C     ----------------------------
C     AN OUTPUT FILE IS WRITEN EVERY TOUT TIME UNITS
C     ----------------------------

      READ (1,*) TXT,LABEL,TOUT
      WRITE(6,*) TXT,LABEL,TXTXT,TOUT
      LABEL1 = 'tout'
      IF (LABEL.NE.LABEL1) GOTO 10

C     -----------------------------
C     DT IS WRITEN ON THE SCREEN EVERY ITSTP TIMESTEPS
C     -----------------------------

      READ (1,*) TXT,LABEL,ITSTP
      WRITE(6,*) TXT,LABEL,TXTXT,ITSTP
      LABEL1 = 'itstp'
      IF (LABEL.NE.LABEL1) GOTO 10

C     -----------------------------
C     CFL IS THE TIME-LIMITING FACTOR (.LT. 1)
C     -----------------------------

      READ (1,*) TXT,LABEL,CFL
      WRITE(6,*) TXT,LABEL,TXTXT,CFL
      LABEL1 = 'cfl'
      IF (LABEL.NE.LABEL1) GOTO 10

C     -------------------
C     DTINI IS THE INITIAL DT
C     -------------------

      READ (1,*) TXT,LABEL,DTINI
      WRITE(6,*) TXT,LABEL,TXTXT,DTINI
      LABEL1 = 'dtini'
      IF (LABEL.NE.LABEL1) GOTO 10

C     -------------------------------------------
C     SMALL IS THE THRESHOLD FOR NONDIMENSIONAL NUMBERS (I.E. VELOCITY)
C     -------------------------------------------

      READ (1,*) TXT,LABEL,SMALL
      WRITE(6,*) TXT,LABEL,TXTXT,SMALL
      LABEL1 = 'small'
      IF (LABEL.NE.LABEL1) GOTO 10 

C     -----------------------------------------
C     SMLRHO IS THE THRESHOLD FOR DENSITIES (RHO, R)
C     -----------------------------------------

      READ (1,*) TXT,LABEL,SMLRHO
      WRITE(6,*) TXT,LABEL,TXTXT,SMLRHO
      LABEL1 = 'smlrho'
      IF (LABEL.NE.LABEL1) GOTO 10

C     --------------------------------
C     SMALLP IS THE THRESHOLD FOR PRESSURE
C     --------------------------------

      READ (1,*) TXT,LABEL,SMALLP
      WRITE(6,*) TXT,LABEL,TXTXT,SMALLP
      LABEL1 = 'smallp'
      IF (LABEL.NE.LABEL1) GOTO 10

C     ---------------------------------
C     SMALLU IS THE THRESHOLD FOR INTERNAL ENERGY
C     ---------------------------------

      READ (1,*) TXT,LABEL,SMALLU
      WRITE(6,*) TXT,LABEL,TXTXT,SMALLU
      LABEL1 = 'smallu'
      IF (LABEL.NE.LABEL1) GOTO 10

C     ------------------------
C     GRIDLX IS THE LENGTH OF THE GRID
C     ------------------------

      READ (1,*) TXT,LABEL,GRIDLX
      WRITE(6,*) TXT,LABEL,TXTXT,GRIDLX
      LABEL1 = 'gridlx'
      IF (LABEL.NE.LABEL1) GOTO 10

      IF ((BASENM.EQ.'RST_'.OR.BASENM.EQ.'RBW_'.OR.BASENM.EQ.'RSR_'.OR.
     &     BASENM.EQ.'RBWI').AND.
     &     GRIDLX.NE.1.D0) THEN

         GRIDLX = 1.D0
         WRITE(6,1200)
 1200    FORMAT('INPUT: GRIDLX RESET TO 1.D0')

      END IF

C     ---------------------------
C     NX IS THE NUMBER OF GRID POINTS
C     ---------------------------

      READ (1,*) TXT,LABEL,NX
      WRITE(6,*) TXT,LABEL,TXTXT,NX
      LABEL1 = 'nx'
      IF (LABEL.NE.LABEL1) GOTO 10

      IF (NX.LT.4.OR.NX.GT.MNX) STOP 'INPUT: UNSUITABLE NX'

C     ----------------------------------------
C     ETA1 IS USED IN SUBROUTINE DETECT (PPM RECONSTRUCTION)
C     TYPICAL VALUE: 5.D0
C     ----------------------------------------

      READ (1,*) TXT,LABEL,ETA1
      WRITE(6,*) TXT,LABEL,TXTXT,ETA1
      LABEL1 = 'eta1'
      IF (LABEL.NE.LABEL1) GOTO 10

C     -------------------------------------------
C     ETA2 IS USED IN SUBROUTINE DETECT (PPM RECONSTRUCTION)
C     TYPICAL VALUE: 5.D-2
C     -------------------------------------------

      READ (1,*) TXT,LABEL,ETA2
      WRITE(6,*) TXT,LABEL,TXTXT,ETA2
      LABEL1 = 'eta2'
      IF (LABEL.NE.LABEL1) GOTO 10

C     ---------------------------------------
C     EPSLN IS USED IN SUBROUTINE DETECT (PPM RECONSTRUCTION)
C     TYPICAL VALUE: 1.D-1
C     ---------------------------------------

      READ (1,*) TXT,LABEL,EPSLN
      WRITE(6,*) TXT,LABEL,TXTXT,EPSLN
      LABEL1 = 'epsln'
      IF (LABEL.NE.LABEL1) GOTO 10

C     -----------------------------------------
C     AK0 IS USED IN SUBROUTINE DETECT (PPM RECONSTRUCTION)
C     TYPICAL VALUE: 1.D0
C     -----------------------------------------
      READ (1,*) TXT,LABEL,AK0
      WRITE(6,*) TXT,LABEL,TXTXT,AK0
      LABEL1 = 'ak0'
      IF (LABEL.NE.LABEL1) GOTO 10

C     ----------------------------------------
C     EPSILN IS USED IN SUBROUTINE FLATEN (PPM RECONSTRUCTION)
C     TYPICAL VALUE: 1.D0
C     ----------------------------------------

      READ (1,*) TXT,LABEL,EPSILN
      WRITE(6,*) TXT,LABEL,TXTXT,EPSILN
      LABEL1 = 'epsiln'
      IF (LABEL.NE.LABEL1) GOTO 10

C     ---------------------------------------
C     OMG1 IS USED IN SUBROUTINE FLATEN (PPM RECONSTRUCTION)
C     TYPICAL VALUE: 5.2D-1
C     ---------------------------------------

      READ (1,*) TXT,LABEL,OMG1
      WRITE(6,*) TXT,LABEL,TXTXT,OMG1
      LABEL1 = 'omg1'
      IF (LABEL.NE.LABEL1) GOTO 10

C     ----------------------------------------
C     OMG2 IS USED IN SUBROUTINE FLATEN (PPM RECONSTRUCTION)
C     TYPICAL VALUE: 1.D1
C     ----------------------------------------

      READ (1,*) TXT,LABEL,OMG2
      WRITE(6,*) TXT,LABEL,TXTXT,OMG2
      LABEL1 = 'omg2'
      IF (LABEL.NE.LABEL1) GOTO 10

      PRINT*, '-------------------------------------------------------'

      RETURN

 10   CONTINUE
      PRINT*, ' '
      WRITE(6,1020)
 1020 FORMAT('INPUT: INCORRECT INPUT DECK')
      WRITE(6,1001) LABEL, LABEL1
 1001 FORMAT('       LABEL = ',A6,'  EXPECTED LABEL = ',A6)
      STOP

      END

C     -------- 
CN    NAME: G R I D
C     --------

CP    PURPOSE: 
CP    ESTABLISHES THE BOUNDARY CONDITIONS AND SETS UP THE NUMERICAL GRID
C
 
CC    COMMENTS: 
CC    BOUNDARY CONDITIONS ARE SPECIFIED FOR A SERIES OF 1D TESTS. BOUNDARIES.
CC    APPROPRIATE BOUNDARIES MUST BE CHOSEN BY THE USER FOR SPECIFIC PROBLEMS.
CC    AN EQUIDISTANT X-GRID IS GENERATED BY DEFAULT, HOWEVER THE CODE IS
CC    SUITED FOR NON-EQUIDISTANT GRIDS

      SUBROUTINE GRID

      IMPLICIT NONE

      INCLUDE 'size'

C     ----------
C     COMMON BLOCKS
C     ----------

      INTEGER         BNDMNX,BNDMXX
      COMMON /BOUN/   BNDMNX,BNDMXX

      INTEGER         NEND,NOUT,ITSTP,NX
      COMMON /INPTI/  NEND,NOUT,ITSTP,NX

      DOUBLEPRECISION X(-4:MN5),XL(-4:MN5),XR(-4:MN5),DX(-4:MN5)
      COMMON /GRD/    X,XL,XR,DX

      DOUBLEPRECISION TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,SMALLU,
     &                GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2
      COMMON /INPTF/  TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,SMALLU,
     &                GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2

      CHARACTER*7     OUTFIL
      CHARACTER*8     LABEL
      CHARACTER*4     BASENM
      CHARACTER*2     SUFFIX
      COMMON /CHRC/   LABEL,OUTFIL,BASENM,SUFFIX

C     -------------
C     INTERNAL VARIABLES
C     -------------

      INTEGER         I

      DOUBLEPRECISION DELX

C     ---------------------------------
C     BOUNDARY CONDITIONS
C     BNDM.. = 1 ===> REFLECTING BOUNDARY
C     BNDM.. = 2 ===> FLOW OUT BOUNDARY
C     BNDM.. = 3 ===> FLOW IN BOUNDARY
C     BNDM.. = 4 ===> PERIODIC BOUNDARY
C     BNDM.. = 5 ===> ANY OTHER BOUNDARY
C     ---------------------------------

      IF (BASENM.EQ.'RST_'.OR.BASENM.EQ.'RBW_'.OR.BASENM.EQ.'RBWI') THEN
         BNDMNX = 2
         BNDMXX = 2
      ELSE IF (BASENM.EQ.'RSR_') THEN
         BNDMNX = 2
         BNDMXX = 1
      ELSE
         BNDMNX = 2
         BNDMXX = 2
      END IF

      IF (BNDMNX.EQ.4.AND.BNDMXX.NE.4) STOP 'GRID: INCORRECT BOUNDARIES'

C     ----------
C     SET UP X-GRID
C     ----------

      X(1) = 0.D0

      DELX = GRIDLX/DFLOAT(NX)

      DO 10 I=2,NX+1
         XL(I) = XL(I-1) + DELX
 10      CONTINUE


      DO 20 I=1,NX
         XR(I) = XL(I+1)
 20      CONTINUE

      DO 30 I=1,NX
         X(I) = 0.5D0*(XL(I) + XR(I))
 30      CONTINUE

      DO 40 I=1,NX
         DX(I) = XR(I) - XL(I)
 40      CONTINUE

      RETURN
      END

C     -------- 
CN    NAME: I N I T
C     --------

CP    PURPOSE: 
CP    DEFINES THE INITIAL MODEL
C
 
CC    COMMENTS: 
CC    DEFINES INITIAL DATA FOR A SERIES OF STANDARD 1D TEST PROBLEMS.
CC    APPROPRIATE INITIAL DATA MUST BE DEFINED BY THE USER FOR SPECIFIC 
CC    PROBLEMS.
   
      SUBROUTINE INIT

      IMPLICIT NONE

      INCLUDE 'size'

C     -----------
C     COMMON BLOCKS
C     -----------

      INTEGER         NEND,NOUT,ITSTP,NX
      COMMON /INPTI/  NEND,NOUT,ITSTP,NX

      DOUBLEPRECISION GAMMA
      COMMON /ADIND/  GAMMA

      DOUBLEPRECISION X(-4:MN5),XL(-4:MN5),XR(-4:MN5),DX(-4:MN5)
      COMMON /GRD/    X,XL,XR,DX

      DOUBLEPRECISION P(-4:MN5),RHO(-4:MN5),VEL(-4:MN5),W(-4:MN5),
     &                U(-4:MN5),CS(-4:MN5),H(-4:MN5),DPDRH(-4:MN5),
     &                DPDU(-4:MN5),R(-4:MN5),M(-4:MN5),E(-4:MN5)
      COMMON /HYDRO/  P,RHO,VEL,W,U,CS,H,DPDRH,DPDU,R,M,E

      DOUBLEPRECISION TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,SMALLU,
     &                GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2
      COMMON /INPTF/  TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,SMALLU,
     &                GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2

      DOUBLEPRECISION TIME,DT
      COMMON /ZEIT/   TIME,DT

      CHARACTER*7     OUTFIL
      CHARACTER*8     LABEL
      CHARACTER*4     BASENM
      CHARACTER*2     SUFFIX
      COMMON /CHRC/   LABEL,OUTFIL,BASENM,SUFFIX

C     -----------
C     INTERNAL VARIABLES
C     -----------

      INTEGER         I

C     ---------
C     INITIAL TIME
C     ---------

      TIME = 0.D0

C     ---------------------
C     RELATIVISTIC SOD'S TUBE
C     ---------------------

      IF (BASENM.EQ.'RST_') THEN

         GAMMA = 1.4D0
         
         DO 100 I=1,NX

            IF (X(I).LE.GRIDLX/2.D0) THEN
               VEL(I) = 0.D0
               RHO(I) = 1.D0
               U(I)   = 2.5D0
            ELSE
               VEL(I) = 0.D0
               RHO(I) = 0.125D0
               U(I)   = 2.D0
            END IF

 100        CONTINUE

         GOTO 450

      END IF

C     ----------------------
C     SCHNEIDER ET AL.'S TEST
C     ----------------------

      IF (BASENM.EQ.'SCHN') THEN

         GAMMA = 5.D0/3.D0
         
         DO 125 I=1,NX

            IF (X(I).LE.GRIDLX/2.D0) THEN
               VEL(I) = 0.D0
               RHO(I) = 10.D0
               U(I)   = 2.D0
            ELSE
               VEL(I) = 0.D0
               RHO(I) = 1.D0
               U(I)   = 1.D-6
            END IF

 125        CONTINUE

         GOTO 450

      END IF

C     --------------------
C     RELATIVISTIC BLAST WAVE
C     --------------------

      IF (BASENM.EQ.'RBW_') THEN

         GAMMA = 5.D0/3.D0

         DO 200 I=1,NX

            IF (X(I).LE.GRIDLX/2.D0) THEN
               VEL(I) = 0.D0
               RHO(I) = 1.D0
               U(I)   = 1.5D3
            ELSE
               VEL(I) = 0.D0
               RHO(I) = 1.D0
               U(I)   = 1.5D-2
            END IF

 200        CONTINUE

         GOTO 450

      END IF

C     -----------------------
C     RELATIVISTIC SHOCK REFLECTION
C     -----------------------

      IF (BASENM.EQ.'RSR_') THEN

         GAMMA = 4.D0/3.D0
         
         DO 300 I=1,NX

            VEL(I) = 0.99999D0
            RHO(I) = 1.D0
            U(I)   = 1.D-7/DSQRT(1.D0 - VEL(I)*VEL(I))

 300        CONTINUE

         GOTO 450

      END IF

C     ------------------------
C     RELATIVISTIC BLAST WAVE INTERACTION
C     ------------------------

      IF (BASENM.EQ.'RBWI') THEN

         GAMMA = 1.4D0
         
         DO 400 I=1,NX

            IF (X(I).LE.0.1D0*GRIDLX) THEN
               VEL(I) = 0.D0
               RHO(I) = 1.D0
               U(I)   = 2.5D3
            ELSE IF (X(I).LE.0.9D0*GRIDLX) THEN
               VEL(I) = 0.D0
               RHO(I) = 1.D0
               U(I)   = 2.5D-2
            ELSE
               VEL(I) = 0.D0
               RHO(I) = 1.D0
               U(I)   = 2.5D2
            END IF

 400        CONTINUE

        GOTO 450
      END IF
      
      STOP 'INIT: NO INITIAL DATA SPECIFIED'

 450  CONTINUE

      CALL EOS (NX, RHO, U, GAMMA, P, H, CS, DPDRH, DPDU)

      DO 500 I=1,NX

            W(I)     = 1.D0/DSQRT(1.D0 - VEL(I)*VEL(I))
      
            R(I)     = RHO(I)*W(I)
            M(I)     = R(I)*H(I)*W(I)*VEL(I)
            E(I)     = R(I)*H(I)*W(I) - P(I) - R(I)

 500        CONTINUE

      RETURN
      END

C     -------- 
CN    NAME: T S T E P 
C     --------

CP    PURPOSE: 
CP    COMPUTES THE NEW TIMESTEP VALUE FROM COURANT CONDITION
C
 
CC    COMMENTS: 
CC    NONE

      SUBROUTINE TSTEP

      IMPLICIT NONE

      INCLUDE 'size'

C     -----------
C     COMMON BLOCKS
C     -----------

      INTEGER         NEND,NOUT,ITSTP,NX
      COMMON /INPTI/  NEND,NOUT,ITSTP,NX

      INTEGER         NSTEP
      COMMON /NSTEP/  NSTEP

      DOUBLEPRECISION X(-4:MN5),XL(-4:MN5),XR(-4:MN5),DX(-4:MN5)
      COMMON /GRD/    X,XL,XR,DX

      DOUBLEPRECISION P(-4:MN5),RHO(-4:MN5),VEL(-4:MN5),W(-4:MN5),
     &                U(-4:MN5),CS(-4:MN5),H(-4:MN5),DPDRH(-4:MN5),
     &                DPDU(-4:MN5),R(-4:MN5),M(-4:MN5),E(-4:MN5)
      COMMON /HYDRO/  P,RHO,VEL,W,U,CS,H,DPDRH,DPDU,R,M,E

      DOUBLEPRECISION TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,SMALLU,
     &                GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2
      COMMON /INPTF/  TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,SMALLU,
     &                GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2

      DOUBLEPRECISION TIME,DT
      COMMON /ZEIT/   TIME,DT

C     ------------
C     INTERNAL VARIABLES
C     ------------

      INTEGER         I,IC

      DOUBLEPRECISION DTCC,DTEST(MN)

      DOUBLEPRECISION LAMBD1,LAMBD4,LAMBDA

      DOUBLEPRECISION V

      IC   = 0
      DTCC = 0.D0

      DO 10 I=1,NX

         LAMBD1   = (VEL(I) - CS(I))/(1.D0 - VEL(I)*CS(I))
         LAMBD4   = (VEL(I) + CS(I))/(1.D0 + VEL(I)*CS(I))
         LAMBDA   = DMAX1(DABS(LAMBD1),DABS(LAMBD4))
         DTEST(I) = LAMBDA/(XR(I) - XL(I))

 10      CONTINUE

      DO 13 I=1,NX

         IF (DTEST(I).GT.DTCC) THEN
            IC   = I
            DTCC = DTEST(I)
         END IF

 13      CONTINUE

      DT = CFL/DTCC

      V  = DABS(VEL(IC))

      IF (MOD(NSTEP,ITSTP).NE.0) RETURN

      WRITE(6,1001) NSTEP,DT,IC,CS(IC),V

 1001 FORMAT(I5,2X,1PE8.1,2X,I5,2X,1P2E11.2,2X,1P2E11.2)

      RETURN
      END

C     -------- 
CN    NAME: B N D R Y
C     --------

CP    PURPOSE: 
CP    PROVIDES DIFFERENT TYPES OF BOUNDARY CONDITIONS
C

CC    COMMENTS: 
CC    NEW BOUNDARY CONDITIONS CAN BE SPECIFIED AT THE USER'S WILL

      SUBROUTINE BNDRY

      IMPLICIT NONE

      INCLUDE 'size'

C     --------
C     COMMON BLOCKS
C     --------

      INTEGER         BNDMNX,BNDMXX
      COMMON /BOUN/   BNDMNX,BNDMXX

      INTEGER         NEND,NOUT,ITSTP,NX
      COMMON /INPTI/  NEND,NOUT,ITSTP,NX

      DOUBLEPRECISION X(-4:MN5),XL(-4:MN5),XR(-4:MN5),DX(-4:MN5)
      COMMON /GRD/    X,XL,XR,DX

      DOUBLEPRECISION P(-4:MN5),RHO(-4:MN5),VEL(-4:MN5),W(-4:MN5),
     &                U(-4:MN5),CS(-4:MN5),H(-4:MN5),DPDRH(-4:MN5),
     &                DPDU(-4:MN5),R(-4:MN5),M(-4:MN5),E(-4:MN5)
      COMMON /HYDRO/  P,RHO,VEL,W,U,CS,H,DPDRH,DPDU,R,M,E

      DOUBLEPRECISION TIME,DT
      COMMON /ZEIT/   TIME,DT
      
C     ---------
C     INTERNAL VARIABLES
C     ---------

      INTEGER         I

C     ------------
C     LEFT BOUNDARY
C     ------------

      GOTO (410, 420, 430, 440, 450), BNDMNX

C     -----------------
C     REFLECTING BOUNDARY
C     -----------------

 410  CONTINUE

      DO 415 I=-4,0

         RHO(I)  = RHO(1-I)
         VEL(I)  = -VEL(1-I)
         U(I)    = U(1-I)
         P(I)    = P(1-I)
         CS(I)   = CS(1-I)
         H(I)    = H(1-I)
         W(I)    = W(1-I)
         DX(I)   = DX(1-I)

 415     CONTINUE

      GOTO 500

C     ---------------
C     FLOW OUT BOUNDARY
C     ---------------

 420  CONTINUE

      DO 425 I=-4,0

         RHO(I)  = RHO(1)
         VEL(I)  = VEL(1)
         U(I)    = U(1)
         P(I)    = P(1)
         CS(I)   = CS(1)
         H(I)    = H(1)
         W(I)    = W(1)
         DX(I)   = DX(1)

 425     CONTINUE

      GOTO 500

C     ---------------
C     FLOW IN BOUNDARY
C     ---------------

 430  CONTINUE

      STOP 'BNDRY: INFLOW BOUNDARY MUST BE SUPPLIED BY THE USER'

      GOTO 500

C     --------------
C     PERIODIC BOUNDARY
C     --------------

 440  CONTINUE

      DO 445 I=-4,0

         RHO(I)  = RHO(NX+I)
         VEL(I)  = VEL(NX+I)
         U(I)    = U(NX+I)
         P(I)    = P(NX+I)
         CS(I)   = CS(NX+I)
         H(I)    = H(NX+I)
         W(I)    = W(NX+I)
         DX(I)   = DX(NX+I)

 445     CONTINUE

      GOTO 500

C     ---------------------------------------------------------
C     SPECIAL BOUNDARY (ADD ANY NONSTANDARD BOUNDARY CONDITION HERE)
C     INFLOW JET BOUNDARY
C     ---------------------------------------------------------

 450  CONTINUE

      STOP 'BNDRY: NON STANDARD BOUNDARY. TO BE SUPPLIED BY THE USER'
            
 500  CONTINUE

      DO 505 I=0,-4,-1

         XL(I) = XL(I+1)-DX(I)
         XR(I) = XR(I+1)-DX(I+1)
         X(I)  = 0.5*(XL(I)+XR(I))

 505     CONTINUE

C     ------------
C     RIGHT BOUNDARY
C     ------------

      GOTO (510, 520, 530, 540, 550), BNDMXX

C     -----------------
C     REFLECTING BOUNDARY
C     -----------------

 510  CONTINUE

      DO 515 I=1,5

         RHO(NX+I)  = RHO(NX+1-I)
         VEL(NX+I)  = -VEL(NX+1-I)
         U(NX+I)    = U(NX+1-I)
         P(NX+I)    = P(NX+1-I)
         CS(NX+I)   = CS(NX+1-I)
         H(NX+I)    = H(NX+1-I)
         W(NX+I)    = W(NX+1-I)
         DX(NX+I)   = DX(NX+1-I)

 515     CONTINUE

      GOTO 600

C     --------------
C     FLOW OUT BOUNDARY
C     --------------

 520  CONTINUE

      DO 525 I=NX+1,NX+5

         RHO(I)  = RHO(NX)
         VEL(I)  = VEL(NX)
         U(I)    = U(NX)
         P(I)    = P(NX)
         CS(I)   = CS(NX)
         H(I)    = H(NX)
         W(I)    = W(NX)
         DX(I)   = DX(NX)

 525     CONTINUE

      GOTO 600

C     --------------
C     FLOW IN BOUNDARY
C     --------------

 530  CONTINUE

      STOP 'BNDRY: INFLOW BOUNDARY MUST BE SUPPLIED BY THE USER'

      GOTO 600

C     ---------------
C     PERIODIC BOUNDARY
C     ---------------

 540  CONTINUE

      DO 545 I=1,5

         RHO(NX+I)  = RHO(I)
         VEL(NX+I)  = VEL(I)
         U(NX+I)    = U(I)
         P(NX+I)    = P(I)
         CS(NX+I)   = CS(I)
         H(NX+I)    = H(I)
         W(NX+I)    = W(I)
         DX(NX+I)   = DX(I)

 545     CONTINUE

      GOTO 600

C     ---------------------------------------
C     SPECIAL BOUNDARY
C     ADD ANY NONSTANDARD BOUNDARY CONDITION HERE
C     ---------------------------------------

 550  CONTINUE

      DO 555 I=NX+1,NX+5

         DX(I)   = DX(NX)
         XL(I)   = XL(I-1) + DX(I-1)
         XR(I)   = XR(I-1) + DX(I)
         X(I)    = 0.5*(XL(I)+XR(I))
         RHO(I)  = 1.D0+DABS(VEL(NX))*TIME/X(I)
         VEL(I)  = VEL(NX)
         U(I)    = U(NX)
         P(I)    = P(NX)
         CS(I)   = CS(NX)
         H(I)    = H(NX)
         W(I)    = W(NX)

 555     CONTINUE

 600  CONTINUE

      DO 650 I=NX+1,NX+5

         XL(I) = XL(I-1) + DX(I-1)
         XR(I) = XR(I-1) + DX(I)
         X(I)  = 0.5*(XL(I)+XR(I))

 650     CONTINUE

      RETURN
      END

C     -------- 
CN    NAME: H Y D R O W
C     --------

CP    PURPOSE: 
CP    ADVANCE IN TIME THE  1D EQUATIONS OF RELATIVISTIC HYDRODYNAMICS 
CP    (IN CONSERVATION FORM) IN PLANAR COORDINATES.
C

CC    COMMENTS: 
CC    NONE

      SUBROUTINE HYDROW

      IMPLICIT NONE

      INCLUDE 'size'

C     ---------
C     COMMON BLOCKS
C     ---------

      INTEGER         NEND,NOUT,ITSTP,NX
      COMMON /INPTI/  NEND,NOUT,ITSTP,NX

      DOUBLEPRECISION COEFF1(-4:MN5),COEFF2(-4:MN5),COEFF3(-4:MN5),
     &                COEFF4(-4:MN5),COEFF5(-4:MN5)
      COMMON /COEFF/  COEFF1,COEFF2,COEFF3,COEFF4,COEFF5

      DOUBLEPRECISION DELP(-4:MN5),DELRHO(-4:MN5),DELVEL(-4:MN5),
     &                DELU(-4:MN5)
      COMMON /DELU/   DELP,DELRHO,DELVEL,DELU

      DOUBLEPRECISION FICT(-4:MN5)
      COMMON /FICT/   FICT

      DOUBLEPRECISION FLATN(-4:MN5),FLATN1(-4:MN5)
      COMMON /FLAT/   FLATN,FLATN1

      DOUBLEPRECISION GAMMA
      COMMON /ADIND/  GAMMA
      
      DOUBLEPRECISION X(-4:MN5),XL(-4:MN5),XR(-4:MN5),DX(-4:MN5)
      COMMON /GRD/    X,XL,XR,DX

      DOUBLEPRECISION P(-4:MN5),RHO(-4:MN5),VEL(-4:MN5),
     &                W(-4:MN5),U(-4:MN5),CS(-4:MN5),H(-4:MN5),
     &                DPDRH(-4:MN5),DPDU(-4:MN5),R(-4:MN5),M(-4:MN5),
     &                E(-4:MN5)
      COMMON /HYDRO/  P,RHO,VEL,W,U,CS,H,DPDRH,DPDU,R,M,E

      DOUBLEPRECISION TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,
     &                SMALLU,GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2
      COMMON /INPTF/  TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,
     &                SMALLU,GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2

      DOUBLEPRECISION PL(-4:MN6),PR(-4:MN6),RHOL(-4:MN6),RHOR(-4:MN6),
     &                VELL(-4:MN6),VELR(-4:MN6),
     &                UL(-4:MN6),UR(-4:MN6),CSL(-4:MN6),
     &                CSR(-4:MN6),RL(-4:MN6),RR(-4:MN6),ML(-4:MN6),
     &                MR(-4:MN6),EL(-4:MN6),
     &                ER(-4:MN6)
      COMMON /INTERF/ PL,PR,RHOL,RHOR,VELL,VELR,UL,UR,CSL,
     &                CSR,RL,RR,ML,MR,EL,ER

      DOUBLEPRECISION TIME,DT
      COMMON /ZEIT/   TIME,DT

      DOUBLEPRECISION DP(-4:MN5),P6(-4:MN5),DRHO(-4:MN5),RHO6(-4:MN5),
     &                DVEL(-4:MN5),VEL6(-4:MN5),
     &                DU(-4:MN5),U6(-4:MN5)
      COMMON /U6/     DP,P6,DRHO,RHO6,DVEL,VEL6,DU,U6

      DOUBLEPRECISION PM(-4:MN5),PP(-4:MN5),RHOM(-4:MN5),RHOP(-4:MN5),
     &                VELM(-4:MN5),VELP(-4:MN5),
     &                UM(-4:MN5),UP(-4:MN5)
      COMMON /UMP/    PM,PP,RHOM,RHOP,VELM,VELP,UM,UP

C     ------------
C     INTERNAL VARIABLES
C     ------------

      INTEGER         I

      DOUBLEPRECISION AUX1

      DOUBLEPRECISION DTDX(-4:MN5)

      DOUBLEPRECISION RFLX(-4:MN6),MFLX(-4:MN6),EFLX(-4:MN6)

C     -------------
C     SPATIAL RECONSTRUCTION (PPM)
C     -------------

      CALL COEF(NX,COEFF1,COEFF2,COEFF3,COEFF4,COEFF5)

      CALL INTERP(NX,PM   ,P   ,PP   ,DELP  )

      CALL INTERP(NX,RHOM ,RHO ,RHOP ,DELRHO)

      CALL DETECT(NX,RHOM ,RHO ,RHOP ,DELRHO)

      CALL INTERP(NX,VELM ,VEL ,VELP ,DELVEL)

      CALL FLATEN


      DO 17 I=0,NX+1
         RHOM(I)  = FLATN(I)*RHO(I) + FLATN1(I)*RHOM(I)
         RHOP(I)  = FLATN(I)*RHO(I) + FLATN1(I)*RHOP(I)
         VELM(I)  = FLATN(I)*VEL(I) + FLATN1(I)*VELM(I)
         VELP(I)  = FLATN(I)*VEL(I) + FLATN1(I)*VELP(I)
         PM(I)    = FLATN(I)*P(I)   + FLATN1(I)*PM(I)
         PP(I)    = FLATN(I)*P(I)   + FLATN1(I)*PP(I)
 17      CONTINUE

      CALL MONOT(NX,PM   ,P   ,PP   ,DP   ,P6   )

      CALL MONOT(NX,RHOM ,RHO ,RHOP ,DRHO ,RHO6 )

      CALL MONOT(NX,VELM ,VEL ,VELP ,DVEL ,VEL6 )

C     -------------------
C     AVERAGED STATES FOR TIME ADVANCE
C     -------------------

      CALL STAT1D

C     ----------------
C     COMPUTATION OF NUMERICAL FLUXES (WITH AN EXACT RIEMANN SOLVER)
C     ----------------

      DO 23 I=1,NX+1 

         CALL NFLUX(RHOL(I),RHOR(I),PL(I),PR(I),VELL(I),VELR(I),
     &              UL(I),UR(I),CSL(I),CSR(I),RFLX(I),MFLX(I),EFLX(I))
 23      CONTINUE
        
C     ----------
C     TIME ADVANCE
C     ----------

      DO 50 I=1,NX
         DTDX(I) = DT/DX(I)
         R(I)    = R(I) - DTDX(I)*(RFLX(I+1) - RFLX(I))
         R(I)    = DMAX1(SMLRHO,R(I))
 50      CONTINUE

      DO 60 I=1,NX
         AUX1  = -DTDX(I)*(MFLX(I+1) - MFLX(I) )
         M(I)  = M(I) + AUX1 + DT*FICT(I)
         AUX1  = -DTDX(I)*(EFLX(I+1) - EFLX(I) )
         E(I)  = E(I) + AUX1
 60      CONTINUE

C     --------------
C     PRIMITIVE VARIABLES
C     --------------

      CALL GETPRFQ(NX,R,M,E,VEL,W,RHO,U,P,H,CS,DPDRH,DPDU)

      RETURN
      END

C     -------- 
CN    NAME: P L T O U T
C     --------

CP    PURPOSE: 
CP    COMPUTES THE ANALYTICAL SOLUTION OF STANDARD 1D PROBLEMS AND WRITES
CP    THE RESULTS IN THE OUTPUT FILES
C
 
CC    COMMENTS: 
CC    NONE

      SUBROUTINE PLTOUT

      IMPLICIT NONE  

      INCLUDE 'size'  

C     -----------
C     COMMON BLOCKS
C     -----------

      INTEGER         NSTEP  
      COMMON /NSTEP/  NSTEP  

      INTEGER         NOUT1
      COMMON /OUTI/   NOUT1  

      INTEGER         NEND,NOUT,ITSTP,NX
      COMMON /INPTI/  NEND,NOUT,ITSTP,NX
 
      DOUBLEPRECISION P(-4:MN5),RHO(-4:MN5),VEL(-4:MN5),W(-4:MN5),
     &                U(-4:MN5),CS(-4:MN5),H(-4:MN5),DPDRH(-4:MN5),
     &                DPDU(-4:MN5),R(-4:MN5),M(-4:MN5),E(-4:MN5)
      COMMON /HYDRO/  P,RHO,VEL,W,U,CS,H,DPDRH,DPDU,R,M,E

      DOUBLEPRECISION X(-4:MN5),XL(-4:MN5),XR(-4:MN5),DX(-4:MN5)
      COMMON /GRD/    X,XL,XR,DX

      DOUBLEPRECISION TOUT1
      COMMON /OUTF/   TOUT1  

      DOUBLEPRECISION TIME,DT  
      COMMON /ZEIT/   TIME,DT  

      DOUBLEPRECISION G  
      COMMON /ADIND/  G  

      CHARACTER*7     OUTFIL  
      CHARACTER*8     LABEL  
      CHARACTER*4     BASENM  
      CHARACTER*2     SUFFIX  
      COMMON /CHRC/   LABEL,OUTFIL,BASENM,SUFFIX  

C     ------------
C     INTERNAL VARIABLES
C     ------------

      INTEGER         I  

      DOUBLEPRECISION P1,RHO1,VEL1,U1,CS1  

      DOUBLEPRECISION OCS2,CS2,FCS2,DFDCS2,CS9  

      DOUBLEPRECISION P2,RHO2,VEL2,U2  

      DOUBLEPRECISION P3,RHO3,VEL3,U3,CS3  

      DOUBLEPRECISION P4,RHO4,VEL4,U4,CS4  

      DOUBLEPRECISION P5,RHO5,VEL5,U5,CS5  

      DOUBLEPRECISION P6,RHO6,VEL6,U6,CS6  

      DOUBLEPRECISION P7,RHO7,VEL7,U7,CS7  

      DOUBLEPRECISION P8,RHO8,VEL8,U8,CS8  

      DOUBLEPRECISION P10,RHO10,VEL10,U10,CS10  

      DOUBLEPRECISION X1,X2,X3,X4,X5,X6,X7,X8,X9  

      DOUBLEPRECISION A,B,C,D,K,L  

      DOUBLEPRECISION PA(MN),RHOA(MN),VELA(MN),UA(MN)  

      DOUBLEPRECISION WW,VS,XS  

C     ----------------------
C     ANALYTICAL SOLUTION FOR STANDARD TESTS
C     ----------------------

      IF (BASENM.EQ.'RSR_') THEN  

        WW   = 1.D0/SQRT(1.D0 - VEL(1)*VEL(1))  
        U1   = 0.D0  
        U2   = WW - 1.D0  
        RHO1 = 1.D0  
        RHO2 = ((G + 1.D0)/(G - 1.D0) + G*U2/(G - 1.D0))*RHO1  
        VEL1 = VEL(1)  
        VEL2 = 0.D0  
        P1   = (G - 1.D0)*RHO1*U1  
        P2   = (G - 1.D0)*RHO2*U2  
        VS   = -VEL1/(RHO2/RHO1/WW - 1.D0)  
        XS   = XR(NX) + VS*TIME
  
        DO 30 I = 1,NX  

          IF (X(I).LT.XS) THEN  
            UA(I)   = U1  
            RHOA(I) = RHO1  
            VELA(I) = VEL1  
            PA(I)   = P1  
          ELSE  
            UA(I)   = U2  
            RHOA(I) = RHO2  
            VELA(I) = VEL2  
            PA(I)   = P2  
          END IF  

 30       CONTINUE  

      END IF  
 
      IF (BASENM.EQ.'RST_') THEN
  
        P1   = 1.D0
        RHO1 = 1.D0  
        VEL1 = 0.D0  
        U1   = P1/(G - 1.D0)/RHO1  
        CS1  = DSQRT(G*P1*(G - 1.D0)/(RHO1*(G - 1.D0) + G*P1))  
 
        P3   = 0.3115D0  
        RHO3 = 0.4345D0  
        VEL3 = 0.4262D0  
        U3   = P3/(G - 1.D0)/RHO3  
        CS3  = DSQRT(G*P3*(G - 1.D0)/(RHO3*(G - 1.D0) + G*P3))  

        P4   = P3  
        RHO4 = 0.273D0  
        VEL4 = VEL3  
        U4   = P4/(G - 1.D0)/RHO4  
        CS4  = DSQRT(G*P4*(G - 1.D0)/(RHO4*(G - 1.D0) + G*P4))  
 
        P5   = 0.1D0  
        RHO5 = 0.125D0  
        VEL5 = 0.D0  
        U5   = P5/(G - 1.D0)/RHO5  
        CS5  = DSQRT(G*P5*(G - 1.D0)/(RHO5*(G - 1.D0) + G*P5))  
 
        X1 = 0.5D0 + (VEL1 - CS1)*TIME/(1.D0 - VEL1*CS1)  
 
        X2 = 0.5D0 + (VEL3 - CS3)*TIME/(1.D0 - VEL3*CS3)  
 
        X3 = 0.5D0 + VEL3*TIME  
 
        X4 = 0.5D0 + VEL4*TIME/(1.D0 - RHO5*DSQRT(1.D0 - VEL4**2)/RHO4)  
 
      ELSE IF(BASENM.EQ.'RBW_') THEN  
 
        P1   = 1000.D0  
        RHO1 = 1.D0  
        VEL1 = 0.D0  
        U1   = P1/(G-1.D0)/RHO1  
        CS1  = DSQRT(G*P1*(G - 1.D0)/(RHO1*(G - 1.D0) + G*P1))  
 
        P3   = 18.6D0  
        RHO3 = 9.15D-2  
        VEL3 = 0.960D0  
        U3   = P3/(G-1.D0)/RHO3  
        CS3  = DSQRT(G*P3*(G - 1.D0)/(RHO3*(G - 1.D0) + G*P3))  
 
        P4   = P3  
        RHO4 = 10.75D0  
        VEL4 = VEL3  
        U4   = P4/(G - 1.D0)/RHO4  
        CS4  = DSQRT(G*P4*(G - 1.D0)/(RHO4*(G - 1.D0) + G*P4))  
 
        P5   = 1.D-2  
        RHO5 = 1.D0  
        VEL5 = 0.D0  
        U5   = P5/(G - 1.D0)/RHO5  
        CS5  = DSQRT(G*P5*(G - 1.D0)/(RHO5*(G - 1.D0) + G*P5))  
 
        X1 = 0.5D0 + (VEL1 - CS1)*TIME/(1.D0 - VEL1*CS1)  
 
        X2 = 0.5D0 + (VEL3 - CS3)*TIME/(1.D0 - VEL3*CS3)  
 
        X3 = 0.5D0 + VEL3*TIME  
 
        X4 = 0.5D0 + VEL4*TIME/(1.D0 - RHO5*DSQRT(1.D0 - VEL4**2)/RHO4)  

      ELSE IF(BASENM.EQ.'SCHN') THEN  

        P1   = 13.33333333D0  
        RHO1 = 10.D0  
        VEL1 = 0.D0  
        U1   = P1/(G - 1.D0)/RHO1  
        CS1  = DSQRT(G*P1*(G - 1.D0)/(RHO1*(G - 1.D0) + G*P1))  
 
        P3   = 1.448D0  
        RHO3 = 2.639D0  
        VEL3 = 0.714D0  
        U3   = P3/(G - 1.D0)/RHO3  
        CS3  = DSQRT(G*P3*(G - 1.D0)/(RHO3*(G - 1.D0) + G*P3))  
 
        P4   = P3  
        RHO4 = 5.071D0  
        VEL4 = VEL3  
        U4   = P4/(G - 1.D0)/RHO4  
        CS4  = DSQRT(G*P4*(G - 1.D0)/(RHO4*(G - 1.D0)+G*P4))  
 
        P5   = 0.666666666666D-6  
        RHO5 = 1.D0  
        VEL5 = 0.D0  
        U5   = P5/(G - 1.D0)/RHO5  
        CS5  = DSQRT(G*P5*(G - 1.D0)/(RHO5*(G - 1.D0) + G*P5))  
 
        X1 = 0.5D0 + (VEL1 - CS1)*TIME/(1.D0 - VEL1*CS1)  
 
        X2 = 0.5D0 + (VEL3 - CS3)*TIME/(1.D0 - VEL3*CS3)  
 
        X3 = 0.5D0 + VEL3*TIME  
 
        X4 = 0.5D0 + VEL4*TIME/(1.D0 - RHO5*DSQRT(1.D0 - VEL4**2)/RHO4)  

      ELSE IF (BASENM.EQ.'RBWI') THEN  

         U1   = 2.5D3  
         RHO1 = 1.D0  
         VEL1 = 0.D0  
         P1   = (G - 1.D0)*U1*RHO1  
         CS1  = DSQRT(G*P1*(G - 1.D0)/(RHO1*(G - 1.D0) + G*P1))  
           
         RHO3 = 0.0491D0  
         VEL3 = 0.957D0  
         P3   = 14.71D0  
         U3   = P3/(G - 1.D0)/RHO3  
         CS3  = DSQRT(G*P3*(G - 1.D0)/(RHO3*(G - 1.D0) + G*P3))  
            
         RHO4 = 14.39D0
         VEL4 = VEL3  
         P4   = P3  
         U4   = P4/(G - 1.D0)/RHO4  
         CS4  = DSQRT(G*P4*(G - 1.D0)/(RHO4*(G - 1.D0) + G*P4))  
           
         IF (TIME.LT.0.4200) THEN  

            U5   = 2.5D-2  
            RHO5 = 1.D0  
            VEL5 = 0.D0  
            P5   = (G - 1.D0)*RHO5*U5  
            CS5  = DSQRT(G*P5*(G - 1.D0)/(RHO5*(G - 1.D0) + G*P5))  
              
            U6   = U5  
            RHO6 = RHO5  
            VEL6 = VEL5  
            P6   = P5  
            CS6  = CS5  

         ELSE  

            RHO5 = 104.41D0  
            VEL5 = 0.456D0  
            P5   = 369.84D0  
            U5   = P5/(G - 1.D0)/RHO5  
            CS5  = DSQRT(G*P5*(G - 1.D0)/(RHO5*(G - 1.D0) + G*P5))  
              
            RHO6 = 117.25D0  
            VEL6 = VEL5  
            P6   = P5  
            U6   = P6/(G - 1.D0)/RHO6  
            CS6  = DSQRT(G*P6*(G - 1.D0)/(RHO6*(G - 1.D0) + G*P6))  

         END IF  

         RHO7 = 9.72D0
         VEL7 = -0.882D0  
         P7   = 4.639D0  
         U7   = P7/(G - 1.D0)/RHO7  
         CS7  = DSQRT(G*P7*(G - 1.D0)/(RHO7*(G - 1.D0) + G*P7))  
           
         RHO8 = 0.112D0  
         VEL8 = VEL7  
         P8   = P7  
         U8   = P8/(G - 1.D0)/RHO8  
         CS8  = DSQRT(G*P8*(G - 1.D0)/(RHO8*(G - 1.D0) + G*P8))  
           
         U10  = 2.5D2  
         RHO10= 1.D0  
         VEL10= 0.D0  
         P10  = (G - 1.D0)*U10*RHO10  
         CS10 = DSQRT(G*P10*(G - 1.D0)/(RHO10*(G - 1.D0) + G*P10))  

         X1 = 0.1D0 - TIME*CS1  
           
         X2 = 0.1D0 + TIME*(VEL3 - CS3)/(1.D0 - VEL3*CS3)  
           
         X3 = 0.1D0 + TIME*VEL3  
           
         IF (TIME.LT.0.4200D0) THEN  
            X4 = 0.1D0 + TIME*0.9776D0  
              
            X5 = 0.9D0 - TIME*0.9274D0  
              
            X6 = X5  
         ELSE  
            X4 = 0.5106D0 + (TIME - 0.4200D0)*0.088D0  
              
            X5 = 0.5106D0 + (TIME - 0.4200D0)*VEL5  
              
            X6 = 0.5106D0 + (TIME - 0.4200D0)*0.703D0  
         END IF   
            
         X7 = 0.9D0 + TIME*VEL7  
           
         X8 = 0.9D0 + TIME*(VEL8 + CS8)/(1.D0 + VEL8*CS8)  
           
         X9 = 0.9D0 + TIME*CS10  

      END IF  
 
      IF (BASENM.EQ.'RST_'.OR.BASENM.EQ.'RBW_'.OR.  
     &    BASENM.EQ.'SCHN') THEN  

        DO 70 I=1,NX  

          IF (X(I).LT.X1) THEN  

            PA(I)   = P1  
            RHOA(I) = RHO1  
            VELA(I) = VEL1  
            UA(I)   = U1  
      
          ELSE IF (X(I).LT.X2) THEN  
         
            A = (X(I) - 0.5D0)/TIME  
            B = DSQRT(G - 1.D0)  
            C = (B + CS1)/(B - CS1)  
            D = -B/2.D0  
            K = (1.D0 + A)/(1.D0 - A)  
            L = C*K**D  
            OCS2 = CS1  
 50         FCS2 = L*(1.D0 + OCS2)**D*(OCS2 - B) +
     &             (1.D0 - OCS2)**D*(OCS2+B)  
            DFDCS2 = L*(1.D0 + OCS2)**D*(1.D0 + D*(OCS2 - B)/
     &               (1.D0 + OCS2)) +  
     &                 (1.D0 - OCS2)**D*(1.D0 - D*(OCS2 + B)/
     &               (1.D0 - OCS2))  
            CS2 = OCS2 - FCS2/DFDCS2  
            
            IF (DABS(CS2 - OCS2)/OCS2.GT.5.D-10) THEN  
               OCS2 = CS2  
               GOTO 50  
            END IF  
         
            VELA(I) = (A + CS2)/(1.D0 + A*CS2)  
            RHOA(I) = RHO1*((CS2**2*(G - 1.D0 - CS1**2))/  
     &                (CS1**2*(G - 1.D0 - CS2**2.)))**(1.D0/(G - 1.D0))  
            PA(I)   = CS2**2*(G - 1.D0)*RHOA(I)/(G - 1.D0 - CS2**2)/G  
            UA(I)   = PA(I)/(G - 1.D0)/RHOA(I)
  
          ELSE IF (X(I).LT.X3) THEN  

            PA(I)   = P3  
            RHOA(I) = RHO3  
            VELA(I) = VEL3  
            UA(I)   = U3  

          ELSE IF (X(I).LT.X4) THEN  
  
            PA(I)   = P4  
            RHOA(I) = RHO4  
            VELA(I) = VEL4  
            UA(I)   = U4  
        
          ELSE
  
            PA(I)   = P5  
            RHOA(I) = RHO5  
            VELA(I) = VEL5  
            UA(I)   = U5  
          END IF  
 
70        CONTINUE  
 
      END IF  

      IF (BASENM.EQ.'RBWI') THEN  

        DO 80 I=1,NX  
      
          IF (X(I).LT.X1) THEN  

            RHOA(I) = RHO1  
            VELA(I) = VEL1  
            PA(I)   = P1  
      
          ELSE IF (X(I).LT.X2) THEN  
 
            A = (X(I) - 0.1D0)/TIME  
            B = DSQRT(G - 1.D0)  
            C = (B + CS1)/(B - CS1)  
            D = -B/2.D0  
            K = (1.D0 + A)/(1.D0 - A)  
            L = C*K**D  
            OCS2 = CS1  
 52         FCS2 = L*(1.D0 + OCS2)**D*(OCS2 - B) +
     &               (1.D0 - OCS2)**D*(OCS2 + B)  
            DFDCS2 = L*(1.D0 + OCS2)**D*(1.D0 + D*(OCS2 - B)/
     &               (1.D0 + OCS2)) +  
     &                 (1.D0 - OCS2)**D*(1.D0 - D*(OCS2 + B)/
     &               (1.D0 - OCS2))  
            CS2 = OCS2 - FCS2/DFDCS2  
         
            IF (ABS(CS2-OCS2)/OCS2.GT.5.E-10)THEN  
               OCS2 = CS2  
               GOTO 52  
            END IF 
 
            VELA(I) = (A+CS2)/(1.D0+A*CS2)  
            RHOA(I) = RHO1*((CS2**2.*(G-1.D0-CS1**2))/  
     &                (CS1**2*(G-1.D0-CS2**2)))**(1.D0/(G-1.D0))  
            PA(I)   = CS2**2*(G-1.D0)*RHOA(I)/(G-1.D0-CS2**2)/G  
            UA(I)   = PA(I)/(G-1.D0)/RHOA(I)  

          ELSE IF (X(I).LT.X3) THEN  

            RHOA(I) = RHO3  
            VELA(I) = VEL3  
            PA(I)   = P3  

          ELSE IF (X(I).LT.X4) THEN  

            RHOA(I) = RHO4  
            VELA(I) = VEL4  
            PA(I)   = P4  

          ELSE IF (X(I).LT.X5) THEN  

            RHOA(I) = RHO5  
            VELA(I) = VEL5  
            PA(I)   = P5  

          ELSE IF (X(I).LT.X6) THEN  

            RHOA(I) = RHO6  
            VELA(I) = VEL6  
            PA(I)   = P6  

          ELSE IF (X(I).LT.X7) THEN  

            RHOA(I) = RHO7  
            VELA(I) = VEL7  
            PA(I)   = P7  

          ELSE IF (X(I).LT.X8) THEN  

            RHOA(I) = RHO8  
            VELA(I) = VEL8  
            PA(I)   = P8  
          
          ELSE IF (X(I).LT.X9) THEN  
         
            A = (X(I) - 0.9D0)/TIME  
            B = SQRT(G - 1.D0)  
            C = (B + CS10)/(B - CS10)  
            D = B/2.D0  
            K = (1.D0 + A)/(1.D0 - A)  
            L = C*K**D  
            OCS2 = CS10  
 54         FCS2 = L*(1.D0 - OCS2)**D*(OCS2 - B) +
     &               (1.D0 + OCS2)**D*(OCS2 + B)  
            DFDCS2 = L*(1.D0 - OCS2)**D*(1.D0 + D*(OCS2 - B)/
     &               (1.D0 - OCS2))+  
     &                 (1.D0 + OCS2)**D*(1.D0 - D*(OCS2 + B)/
     &               (1.D0 + OCS2))  
            CS9 = OCS2 - FCS2/DFDCS2  
         
            IF (DABS(CS9-OCS2)/OCS2.GT.5.D-10)THEN  
               OCS2 = CS9  
               GOTO 54  
            END IF
  
            VELA(I) = (A - CS9)/(1.D0 - A*CS9)  
            RHOA(I) = RHO10*((CS9 **2*(G - 1.D0 - CS10**2))/  
     &               (CS10**2*(G - 1.D0 - CS9 **2)))**(1.D0/(G - 1.D0))  
            PA(I)   = CS9**2*(G - 1.D0)*RHOA(I)/(G - 1.D0 - CS9**2)/G  
            UA(I)   = PA(I)/(G - 1.D0)/RHOA(I)  
 
          ELSE   
    
            RHOA(I) = RHO10  
            VELA(I) = VEL10  
            PA(I)   = P10  

         END IF  

 80      CONTINUE
  
      END IF  

C     ----
C     OUTPUT 
C     ----

      OPEN(10,FILE='DATA/'//OUTFIL,FORM='FORMATTED',STATUS='NEW')  
      WRITE(10,111) NSTEP,TIME
 111  FORMAT('N =', I6, 3X, 'TIME = ', 1PE10.3)   

      DO 85 I=1,NX  
        WRITE(10,200) X(I),P(I),PA(I),RHO(I),RHOA(I),  
     &                VEL(I),VELA(I)  
 85     CONTINUE  

 200  FORMAT(F6.4,1X,2(F11.4,1X),2(F9.4,1X),2(F7.5,1X))  
      CLOSE(10)  

      NOUT1 = 0  
      TOUT1 = 0.D0  

      CALL FILNAM  

      RETURN  
      END  

C     -------- 
CN    NAME: E O S 
C     --------

CP    PURPOSE: 
CP    COMPUTES DERIVED THERMODYNAMICAL QUANTITIES
C

CC    COMMENTS: 
CC    GAMMA-LAW EOS

      SUBROUTINE EOS(N,RHO,U,G,P,H,CS,DPDRH,DPDU)

      IMPLICIT NONE

      INCLUDE 'size'

C     -------
C     ARGUMENTS
C     -------

      INTEGER         N

      DOUBLEPRECISION U(-4:MN5),RHO(-4:MN5),P(-4:MN5),H(-4:MN5),
     &                CS(-4:MN5),DPDRH(-4:MN5),DPDU(-4:MN5)
      
      DOUBLEPRECISION G

C     -----------
C     INTERNAL VARIABLES
C     -----------

      INTEGER         I

      DOUBLEPRECISION GAM1

      DO 10 I=1,N
         GAM1     = G -1.D0
         P(I)     = GAM1*RHO(I)*U(I)
         DPDRH(I) = GAM1*U(I)
         DPDU(I)  = GAM1*RHO(I)
         H(I)     = 1.D0 + U(I) + P(I)/RHO(I)
         CS(I)    = DSQRT((DPDRH(I) + P(I)*DPDU(I)/RHO(I)/RHO(I))/H(I))
 10      CONTINUE

      RETURN
      END

C     -------- 
CN    NAME: C O E F
C     --------

CP    PURPOSE: 
CP    COMPUTES THE COEFFICIENTS FOR INTERPOLATED VALUES 
C

CC    COMMENTS: 
CC    COMPUTES DE GRID-DEPENDENT COEFFICIENTS APPEARING IN EQS.59-61 OF MARTI
CC    AND MUELLER (1996), JCP, VOL. 123, 1-14

      SUBROUTINE COEF(N,COEFF1,COEFF2,COEFF3,COEFF4,COEFF5)
 
      IMPLICIT NONE
 
      INCLUDE 'size'
 
C     --------
C     ARGUMENTS
C     --------

      INTEGER N
 
      DOUBLEPRECISION COEFF1(-4:MN5),COEFF2(-4:MN5),COEFF3(-4:MN5)
      DOUBLEPRECISION COEFF4(-4:MN5),COEFF5(-4:MN5)

C     ----------
C     COMMON BLOCKS
C     ----------

      DOUBLEPRECISION X(-4:MN5),XL(-4:MN5),XR(-4:MN5),DX(-4:MN5)
      COMMON /GRD/    X,XL,XR,DX

C     --------------
C     INTERNAL VARIABLES
C     --------------
      
      INTEGER         I

      DOUBLEPRECISION SCRCH1(-4:MN6),SCRCH2(-4:MN6),SCRCH3(-4:MN6),
     &                SCRCH4(-4:MN6)
 
 
      DO 10 I=0,N+2
 
        SCRCH1(I) = DX(I) + DX(I-1)
        SCRCH2(I) = SCRCH1(I) + DX(I)
        SCRCH3(I) = SCRCH1(I) + DX(I-1)

10      CONTINUE
 
      DO 20 I=0,N+1

        SCRCH4(I) = DX(I)/(SCRCH1(I) + DX(I+1))
        COEFF1(I) = SCRCH4(I)*SCRCH3(I)/SCRCH1(I+1)
        COEFF2(I) = SCRCH4(I)*SCRCH2(I+1)/SCRCH1(I)

20      CONTINUE
 
      DO 30 I=0,N
        
        SCRCH4(I) = 1.D0/(SCRCH1(I) + SCRCH1(I+2))
        COEFF3(I) = -SCRCH4(I)*DX(I)*SCRCH1(I)/SCRCH3(I+1)
        COEFF4(I) = SCRCH4(I)*DX(I+1)*SCRCH1(I+2)/SCRCH2(I+1)
        COEFF5(I) = DX(I) - 2.D0*(DX(I+1)*COEFF3(I) + DX(I)*COEFF4(I))
        COEFF5(I) = COEFF5(I)/SCRCH1(I+1)

30      CONTINUE
  
      RETURN
      END
 
C     -------- 
CN    NAME: I N T E R P
C     --------

CP    PURPOSE: 
CP    COMPUTES INTERPOLATED VALUES AT INTERFACES
C

CC    COMMENTS: 
CC    STEP 1 IN THE RECONSTRUCTION PROCEDURE (SEE APPENDIX I IN MARTI 
CC    & MUELLER 1996)

      SUBROUTINE INTERP(N,UM,U,UP,DELU)
 
      IMPLICIT NONE
 
      INCLUDE 'size'

C     --------
C     ARGUMENTS
C     --------

      INTEGER N
 
      DOUBLEPRECISION UM(-4:MN5),U(-4:MN5),UP(-4:MN5),DELU(-4:MN5)

C     ---------
C     COMMON BLOCKS
C     ---------
 
      DOUBLEPRECISION COEFF1(-4:MN5),COEFF2(-4:MN5),COEFF3(-4:MN5)
      DOUBLEPRECISION COEFF4(-4:MN5),COEFF5(-4:MN5)
      COMMON /COEFF/  COEFF1,COEFF2,COEFF3,COEFF4,COEFF5
 
C     -----------
C     INTERNAL VARIABLES
C     -----------

      INTEGER         I

      DOUBLEPRECISION SCRCH1(-4:MN6),SCRCH2(-4:MN6)
 
      DOUBLEPRECISION SDELU
 
      DO 10 I=-2,N+3

        SCRCH1(I) = U(I) - U(I-1)

10    CONTINUE
 
C     ------------------------------------------------------
C     DELU(I) AS IN EQ.61 OF MARTI AND MUELLER (1996), JCP, VOL. 123, 1-14
C     ------------------------------------------------------

      DO 20 I=0,N+1

        DELU(I) = COEFF1(I)*SCRCH1(I+1) + COEFF2(I)*SCRCH1(I)

20    CONTINUE
 
C     -------------------------------------------------------
C     DELU(I) AS IN EQ.60 OF MARTI AND MUELLER (1996), JCP, VOL. 123, 1-14
C     -------------------------------------------------------

      DO 30 I=0,N+1

        IF (SCRCH1(I+1)*SCRCH1(I).GT.0.D0) THEN
          SDELU     = DELU(I)/DABS(DELU(I))
          SCRCH2(I) = MIN(DABS(SCRCH1(I)),DABS(SCRCH1(I+1)))
          DELU(I)   = MIN(DABS(DELU(I)),2.D0*SCRCH2(I))*SDELU
        ELSE
          DELU(I)   = 0.D0
        END IF

30    CONTINUE

C     ----------------------------------------------
C     INTERFACE VALUES AS IN  EQ.59 OF MARTI AND MUELLER (1996), JCP, 
C     VOL. 123, 1-14
C     ----------------------------------------------

      DO 40 I=0,N

        UP(I)   = U(I)  + COEFF5(I)*SCRCH1(I+1) + COEFF3(I)*DELU(I+1)
        UP(I)   = UP(I) + COEFF4(I)*DELU(I)
        UM(I+1) = UP(I)

40    CONTINUE
 
      RETURN
      END
 
C     -------- 
CN    NAME: D E T E C T
C     --------

CP    PURPOSE: 
CP    DETECTS CONTACT DISCONTINUITIES AND STEEPENS THE CORRESPONDING
CP    RECONSTRUCTED VALUES AT INTERFACES
C

CC    COMMENTS: 
CC    STEP 2 IN THE RECONSTRUCTION PROCEDURE (SEE APPENDIX I IN MARTI 
CC    & MUELLER 1996)

      SUBROUTINE DETECT(N,UM,U,UP,DELU)
 
      IMPLICIT NONE
 
      INCLUDE 'size'
 
C     -------
C     ARGUMENTS
C     -------

      INTEGER         N
 
      DOUBLEPRECISION UM(-4:MN5),U(-4:MN5),UP(-4:MN5),DELU(-4:MN5)

C     -----------
C     COMMON BLOCKS
C     -----------
 
      DOUBLEPRECISION X(-4:MN5),XL(-4:MN5),XR(-4:MN5),DX(-4:MN5)
      COMMON /GRD/    X,XL,XR,DX

      DOUBLEPRECISION P(-4:MN5),RHO(-4:MN5),VEL(-4:MN5),W(-4:MN5),
     &                UU(-4:MN5),CS(-4:MN5),H(-4:MN5),DPDRH(-4:MN5),
     &                DPDU(-4:MN5),R(-4:MN5),M(-4:MN5),E(-4:MN5)
      COMMON /HYDRO/  P,RHO,VEL,W,UU,CS,H,DPDRH,DPDU,R,M,E
 
      DOUBLEPRECISION TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,SMALLU,
     &                GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2
      COMMON /INPTF/  TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,SMALLU,
     &                GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2
 
      DOUBLEPRECISION GB
      COMMON /ADIND/  GB
 
C     ---------
C     INTERNAL VARIABLES
C     ---------

      INTEGER         I

      DOUBLEPRECISION D2U(-4:MN6),ETA(-4:MN6),ETATIL(-4:MN6)
      
      DOUBLEPRECISION SCRCH1(-4:MN6),SCRCH2(-4:MN6),SCRCH3(-4:MN6),
     &                SCRCH4(-4:MN6)

 
      DO 10 I=-1,N+3

        SCRCH1(I) = DX(I) + DX(I-1)
        SCRCH2(I) = SCRCH1(I) + DX(I+1)
        SCRCH3(I) = U(I) - U(I-1)
        SCRCH1(I) = SCRCH3(I)/SCRCH1(I)

10      CONTINUE
 
      DO 20 I=-1,N+2
        
        D2U(I)    = (SCRCH1(I+1) - SCRCH1(I))/SCRCH2(I)
        SCRCH4(I) = X(I) - X(I-1)
        SCRCH4(I) = SCRCH4(I)*SCRCH4(I)*SCRCH4(I)

20      CONTINUE
 
      DO 30 I=0,N+1

        SCRCH1(I) = D2U(I+1)*D2U(I-1)
        SCRCH3(I) = DABS(U(I+1) - U(I-1))
        SCRCH3(I) = SCRCH3(I) - EPSLN*MIN(DABS(U(I+1)),DABS(U(I-1)))

30      CONTINUE
 
C     ------------------------------------------------------
C     ETATIL(I) AS IN EQ.67 OF MARTI AND MUELLER (1996), JCP, VOL. 123, 1-14
C     -------------------------------------------------------

      DO 40 I=0,N+1
      
        IF ((U(I+1) - U(I-1)).EQ.0.D0) THEN
           SCRCH2(I) = SMALL*SMLRHO
        ELSE
           SCRCH2(I) = U(I+1) - U(I-1)
        END IF

        IF ((SCRCH1(I).GT.0.D0).OR.(SCRCH3(I).LT.0.D0)) THEN
           ETATIL(I) = 0.D0
        ELSE
           ETATIL(I) = (D2U(I-1) - D2U(I+1))*(SCRCH4(I) + SCRCH4(I+1))
           ETATIL(I) = ETATIL(I)/(X(I+1) - X(I-1))/SCRCH2(I)
        END IF

40      CONTINUE
 
C     -------------------------------------------------------
C     ETA(I) AS IN EQ.66 OF MARTI AND MUELLER (1996), JCP, VOL. 123, 1-14
C     ONLY FOR ZONES VERIFYING EQ.63 (OTHERWISE, ZERO)
C     -------------------------------------------------------

      DO 50 I=0,N+1

        ETA(I)    = MAX(0.D0,MIN(ETA1*(ETATIL(I) - ETA2),1.D0))
        SCRCH1(I) = DABS(P  (I+1) - P  (I-1))/MIN(P  (I+1),P  (I-1))
        SCRCH2(I) = DABS(RHO(I+1) - RHO(I-1))/MIN(RHO(I+1),RHO(I-1))

50      CONTINUE
 
      DO 60 I=0,N+1

        IF (GB*AK0*SCRCH2(I).LT.SCRCH1(I)) THEN
          ETA(I) = 0.D0
        END IF

60      CONTINUE
 
C     ----------------------------
C     NEW RECONSTRUCTED VALUES (EQ.65)
C     ----------------------------

      DO 70 I=0,N+1

        SCRCH1(I) = U(I-1) + 0.5D0*DELU(I-1)
        SCRCH2(I) = U(I+1) - 0.5D0*DELU(I+1)
        UM(I)     = UM(I)*(1.D0-ETA(I)) + SCRCH1(I)*ETA(I)
        UP(I)     = UP(I)*(1.D0-ETA(I)) + SCRCH2(I)*ETA(I)

70    CONTINUE
 
      RETURN
      END
 
C     -------- 
CN    NAME: M O N O T 
C     --------

CP    PURPOSE: 
CP    RECONSTRUCTED VALUES OBTAINED IN STEPS 1 TO 3 ARE MODIFIED SUCH THAT 
CP    THE INTERPOLATION PARABOLA IN EACH ZONE IS MONOTONE
C

CC    COMMENTS: 
CC    STEP 4 IN THE RECONSTRUCTION PROCEDURE (SEE APPENDIX I IN MARTI 
CC    & MUELLER 1996)

      SUBROUTINE MONOT(N,UM,U,UP,DU,U6)
 
      IMPLICIT NONE
 
      INCLUDE 'size'

C     --------
C     ARGUMENTS
C     --------

      INTEGER         N
 
      DOUBLEPRECISION UM(-4:MN5),U(-4:MN5),UP(-4:MN5),DU(-4:MN5)
      DOUBLEPRECISION U6(-4:MN5)
 
C     --------
C     INTERNAL VARIABLES
C     --------

      INTEGER         I

      DOUBLEPRECISION SCRCH1(-4:MN6),SCRCH2(-4:MN6),SCRCH3(-4:MN6)
  
C     -----------------------------------------------------
C     NEW RECONSRUCTED VALUES IF CONDITION IN EQ.73 OF MARTI & MUELLER
C     1996 HOLDS
C     -----------------------------------------------------

      DO 10 I=0,N+1

        DU(I)     = UP(I) - UM(I)
        SCRCH1(I) = UP(I) - U(I)
        SCRCH1(I) = SCRCH1(I)*(UM(I)-U(I))

10      CONTINUE
 
      DO 20 I=0,N+1

        IF (SCRCH1(I).GE.0.D0) THEN
          UM(I) = U(I)
          UP(I) = U(I)
        END IF

20      CONTINUE

C     --------------------------------------------------------
C     NEW RECONSTRUCTED VALUES IF CONDITION IN EQ.74 OR EQ.75 OF MARTI 
C     & MUELLER 1996 HOLDS
C     --------------------------------------------------------
 
      DO 30 I=0,N+1

        DU(I)     = UP(I) - UM(I)
        SCRCH1(I) = (UP(I) - U(I))*(UM(I) - U(I))

        IF (SCRCH1(I).EQ.0.D0) THEN
          SCRCH2(I) = UM(I)
          SCRCH3(I) = UP(I)
        ELSE
          SCRCH2(I) = 3.D0*U(I) - 2.D0*UP(I)
          SCRCH3(I) = 3.D0*U(I) - 2.D0*UM(I)
        END IF

30      CONTINUE
 
      DO 40 I=0,N+1

        IF (DU(I)*(UM(I) - SCRCH2(I)).LT.0.D0) THEN
          UM(I) = SCRCH2(I)
        END IF

        IF (DU(I)*(SCRCH3(I) - UP(I)).LT.0.D0) THEN
          UP(I) = SCRCH3(I)
        END IF

40      CONTINUE
 
C     -------------------------------------------
C     COMPUTATION OF AUXILIAR VARIABLES FOR TIME ADVANCE
C     -------------------------------------------
 
      DO 50 I=0,N+1

        DU(I) = UP(I) - UM(I)
        U6(I) = 6.D0*U(I) - 3.D0*(UM(I) + UP(I))

50      CONTINUE
 
      RETURN
      END

C     -------- 
CN    NAME: F L A T E N  
C     --------

CP    PURPOSE: 
CP    THIS SUBROUTINE FLATTENS ZONE STRUCTURE IN REGIONS WHERE SHOCKS
CP    ARE TOO THIN    
C

CC    COMMENTS: 
CC    STEP 3 IN THE RECONSTRUCTION PROCEDURE (SEE APPENDIX I IN MARTI 
CC    & MUELLER 1996)

      SUBROUTINE FLATEN

      IMPLICIT NONE
 
      INCLUDE 'size'

C     -------------
C     COMMON BLOCKS
C     -------------

      INTEGER         NEND,NOUT,ITSTP,NX
      COMMON /INPTI/  NEND,NOUT,ITSTP,NX

      DOUBLEPRECISION FLATN(-4:MN5),FLATN1(-4:MN5)
      COMMON /FLAT/   FLATN,FLATN1

      DOUBLEPRECISION P(-4:MN5),RHO(-4:MN5),VEL(-4:MN5),W(-4:MN5),
     &                U(-4:MN5),CS(-4:MN5),H(-4:MN5),DPDRH(-4:MN5),
     &                DPDU(-4:MN5),R(-4:MN5),M(-4:MN5),E(-4:MN5)
      COMMON /HYDRO/  P,RHO,VEL,W,U,CS,H,DPDRH,DPDU,R,M,E

      DOUBLEPRECISION TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,SMALLU,
     &                GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2
      COMMON /INPTF/  TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,SMALLU,
     &                GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2
 
C     ------------
C     INTERNAL VARIABLES
C     ------------

      INTEGER         I

      DOUBLEPRECISION SCRCH1(-4:MN6),SCRCH2(-4:MN6),SCRCH3(-4:MN6)
 
      DOUBLEPRECISION DP(-4:MN5),DVEL(-4:MN5),DP2
 
      DO 20 I=-2,NX+3

        DP(I)     = P(I+1)   - P(I-1)
        DVEL(I)   = VEL(I+1) - VEL(I-1)
        SCRCH1(I) = EPSILN*DMIN1(P(I+1),P(I-1)) - DABS(DP(I))

        IF (SCRCH1(I).LT.0.D0.AND.DVEL(I).LT.0.D0) THEN
          SCRCH1(I) = 1.D0
        ELSE
          SCRCH1(I) = 0.D0
        END IF

20      CONTINUE
 
      DO 30 I=-1,NX+2

        DP2 = P(I+2) - P(I-2)

          IF (DP2.EQ.0.D0) THEN

            IF (DP(I).EQ.0.D0) THEN
              SCRCH2(I) = -OMG1
            ELSE
              SCRCH2(I) = 1.D0 - OMG1
            END IF

          ELSE
            SCRCH2(I) = DP(I)/DP2 - OMG1
        END IF

        SCRCH3(I) = SCRCH1(I)*DMAX1(0.D0,SCRCH2(I)*OMG2)

30    CONTINUE
 
      DO 40 I=0,NX+1

        IF (DP(I).LT.0.D0) THEN
          SCRCH2(I) = SCRCH3(I+1)
        ELSE
          SCRCH2(I) = SCRCH3(I-1)
        END IF

40    CONTINUE
 
      DO 45 I=0,NX+1

        FLATN(I)  = DMAX1(SCRCH2(I),SCRCH3(I))
        FLATN(I)  = DMAX1(0.D0,DMIN1(1.D0,FLATN(I)))
        FLATN1(I) = 1.D0 - FLATN(I)

45    CONTINUE
 
      RETURN
      END


C     -------- 
CN    NAME: S T A T 1 D 
C     --------

CP    PURPOSE: 
CP    THIS SUBROUTINE CALCULATES EFFECTIVE SECOND-ORDER-ACCURATE LEFT 
CP    AND RIGHT STATES FOR RIEMANN PROBLEMS IN ONE DIMENSIONAL 
CP    CALCULATIONS.
C

CC    COMMENTS: 
CC    THIS ROUTINE CLOSELY FOLLOWS THE ANALYTICAL DEVELOPMENTS DESCRIBED IN 
CC    MARTI & MUELLER, JCP, 1996

      SUBROUTINE STAT1D

      IMPLICIT NONE
 
      INCLUDE 'size'

C     ---------
C     COMMON BLOCKS
C     ---------

      INTEGER         NEND,NOUT,ITSTP,NX
      COMMON /INPTI/  NEND,NOUT,ITSTP,NX

      DOUBLEPRECISION TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,SMALLU,
     &                GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2
      COMMON /INPTF/  TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,SMALLU,
     &                GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2

      DOUBLEPRECISION P(-4:MN5),RHO(-4:MN5),VEL(-4:MN5),W(-4:MN5),
     &                U(-4:MN5),CS(-4:MN5),H(-4:MN5),DPDRH(-4:MN5),
     &                DPDU(-4:MN5),R(-4:MN5),M(-4:MN5),E(-4:MN5)
      COMMON /HYDRO/  P,RHO,VEL,W,U,CS,H,DPDRH,DPDU,R,M,E

      DOUBLEPRECISION PM(-4:MN5),PP(-4:MN5),RHOM(-4:MN5),RHOP(-4:MN5),
     &                VELM(-4:MN5),VELP(-4:MN5),UM(-4:MN5),UP(-4:MN5)
      COMMON /UMP/    PM,PP,RHOM,RHOP,VELM,VELP,UM,UP

      DOUBLEPRECISION DP(-4:MN5),P6(-4:MN5),DRHO(-4:MN5),RHO6(-4:MN5),
     &                DVEL(-4:MN5),VEL6(-4:MN5),DU(-4:MN5),U6(-4:MN5)
      COMMON /U6/     DP,P6,DRHO,RHO6,DVEL,VEL6,DU,U6

      DOUBLEPRECISION PL(-4:MN6),PR(-4:MN6),RHOL(-4:MN6),RHOR(-4:MN6),
     &                VELL(-4:MN6),VELR(-4:MN6),UL(-4:MN6),UR(-4:MN6),
     &                CSL(-4:MN6),CSR(-4:MN6),RL(-4:MN6),RR(-4:MN6),
     &                ML(-4:MN6),MR(-4:MN6),EL(-4:MN6),ER(-4:MN6)
      COMMON /INTERF/ PL,PR,RHOL,RHOR,VELL,VELR,UL,UR,CSL,CSR,RL,RR,ML,
     &                MR,EL,ER

      DOUBLEPRECISION FICT(-4:MN5)
      COMMON /FICT/   FICT

      DOUBLEPRECISION TIME,DT
      COMMON /ZEIT/   TIME,DT

      DOUBLEPRECISION X(-4:MN5),XL(-4:MN5),XR(-4:MN5),DX(-4:MN5)
      COMMON /GRD/    X,XL,XR,DX

      DOUBLEPRECISION GB
      COMMON /ADIND/  GB
 
C     ---------
C     INTERNAL VARIABLES
C     ---------

      INTEGER I
 
      DOUBLEPRECISION LAMB1(-4:MN5),LAMB2(-4:MN5),LAMB3(-4:MN5)
 
      DOUBLEPRECISION P1L(-4:MN6),P1R(-4:MN6),P2L(-4:MN6),P2R(-4:MN6),
     &                P3L(-4:MN6),P3R(-4:MN6),RHO1L(-4:MN6),
     &                RHO1R(-4:MN6),RHO2L(-4:MN6),RHO2R(-4:MN6),
     &                RHO3L(-4:MN6),RHO3R(-4:MN6),VEL1L(-4:MN6),
     &                VEL1R(-4:MN6),VEL2L(-4:MN6),VEL2R(-4:MN6),
     &                VEL3L(-4:MN6),VEL3R(-4:MN6)
 
      DOUBLEPRECISION BETA1L(-4:MN6),BETA2L(-4:MN6),BETA3L(-4:MN6),
     &                BETA1R(-4:MN6),BETA2R(-4:MN6),BETA3R(-4:MN6)

 
      DOUBLEPRECISION SCRCH1(-4:MN6),SCRCH2(-4:MN6),SCRCH3(-4:MN6),
     &                SCRCH4(-4:MN6),SCRCH5(-4:MN6)

      DO 10 I=0,NX+1

        LAMB1(I) = (VEL(I) - CS(I))/(1.D0 - VEL(I)*CS(I))
        LAMB2(I) =  VEL(I)
        LAMB3(I) = (VEL(I) + CS(I))/(1.D0 + VEL(I)*CS(I))

10      CONTINUE

      CALL AVRG1D(NX,DT,LAMB1,PM,PP,DP,P6,RHOM,RHOP,DRHO,RHO6,
     &            VELM,VELP,DVEL,VEL6,
     &            P1L,P1R,RHO1L,RHO1R,VEL1L,VEL1R)
 
      CALL AVRG1D(NX,DT,LAMB2,PM,PP,DP,P6,RHOM,RHOP,DRHO,RHO6,
     &            VELM,VELP,DVEL,VEL6,
     &            P2L,P2R,RHO2L,RHO2R,VEL2L,VEL2R)
 
      CALL AVRG1D(NX,DT,LAMB3,PM,PP,DP,P6,RHOM,RHOP,DRHO,RHO6,
     &            VELM,VELP,DVEL,VEL6,
     &            P3L,P3R,RHO3L,RHO3R,VEL3L,VEL3R)
 
C     -------------
C     EFFECTIVE LEFT STATES
C     -------------
 
      DO 20 I=1,NX+1

        SCRCH1(I) = P3L(I)/(GB - 1.D0)/RHO3L(I)
        SCRCH1(I) = 1.D0 + SCRCH1(I) + P3L(I)/RHO3L(I)
        SCRCH2(I) = DSQRT(GB*P3L(I)/RHO3L(I)/SCRCH1(I))
        SCRCH1(I) = SCRCH1(I)*SCRCH2(I)
        SCRCH3(I) = 1.D0/(1.D0 - VEL3L(I)**2)
        SCRCH3(I) = RHO3L(I)*SCRCH3(I)
        SCRCH4(I) = 0.D0
        SCRCH5(I) = 0.D0
 
        BETA1L(I) = 0.5D0*(VEL3L(I) - VEL1L(I) - (P3L(I) - P1L(I))/
     &              SCRCH3(I)/SCRCH1(I) - DT*SCRCH4(I))
        BETA2L(I) = RHO3L(I) - RHO2L(I) - (P3L(I) - P2L(I))/
     &              SCRCH1(I)/SCRCH2(I)
        BETA3L(I) = -0.5D0*DT*SCRCH5(I)

20      CONTINUE
 
      DO 25 I=1,NX+1

        IF (LAMB3(I-1).LE.0.D0) THEN
          BETA3L(I) = 0.D0
        END IF

        IF (LAMB2(I-1).LE.0.D0) THEN
          BETA2L(I) = 0.D0
        END IF

        IF (LAMB1(I-1).LE.0.D0) THEN
          BETA1L(I) = 0.D0
        END IF

25      CONTINUE
 
      DO 30 I=1,NX+1
        
        PL(I)     = P3L(I) + SCRCH3(I)*SCRCH1(I)*(BETA1L(I) + BETA3L(I))
        PL(I)     = DMAX1(SMALLP,PL(I))
 
        RHOL(I)   = RHO3L(I) + SCRCH3(I)/SCRCH2(I)*(BETA1L(I) + 
     &              BETA3L(I)) - BETA2L(I)
        RHOL(I)   = DMAX1(SMLRHO,RHOL(I))
 
        VELL(I)   = VEL3L(I) + BETA3L(I) - BETA1L(I)
 
        UL(I)     = PL(I)/(GB - 1.D0)/RHOL(I)
 
        SCRCH1(I) = 1.D0 + UL(I) + PL(I)/RHOL(I)
 
        CSL(I)    = DSQRT(GB*PL(I)/RHOL(I)/SCRCH1(I))
 
30      CONTINUE
 
C     ----------
C     EFFECTIVE RIGHT STATES
C     ----------

      DO 35 I=1,NX+1

        SCRCH1(I) = P1R(I)/(GB - 1.D0)/RHO1R(I)
        SCRCH1(I) = 1.D0 + SCRCH1(I) + P1R(I)/RHO1R(I)
        SCRCH2(I) = DSQRT(GB*P1R(I)/RHO1R(I)/SCRCH1(I))
        SCRCH1(I) = SCRCH1(I)*SCRCH2(I)
        SCRCH3(I) = 1.D0/(1.D0 - VEL1R(I)**2)
        SCRCH3(I) = RHO1R(I)*SCRCH3(I)
        SCRCH4(I) = 0.D0
        SCRCH5(I) = 0.D0
 
        BETA1R(I) = -0.5D0*DT*SCRCH4(I)
        BETA2R(I) = RHO1R(I) - RHO2R(I) - (P1R(I)-P2R(I))/
     &              SCRCH1(I)/SCRCH2(I)
        BETA3R(I) = -0.5D0*(VEL1R(I) - VEL3R(I) + (P1R(I) - P3R(I))/
     &              SCRCH3(I)/SCRCH1(I) + DT*SCRCH5(I))
35      CONTINUE
 
      DO 40 I=1,NX+1

        IF (LAMB3(I).GE.0.D0) THEN
          BETA3R(I) = 0.D0
        END IF

        IF (LAMB2(I).GE.0.D0) THEN
          BETA2R(I) = 0.D0
        END IF

        IF (LAMB1(I).GE.0.D0) THEN
          BETA1R(I) = 0.D0
        END IF

40    CONTINUE
 
      DO 45 I=1,NX+1

        PR(I)     = P1R(I) + SCRCH3(I)*SCRCH1(I)*(BETA1R(I) + BETA3R(I))
        PR(I)     = DMAX1(SMALLP,PR(I))
 
        RHOR(I)   = RHO1R(I) + SCRCH3(I)/SCRCH2(I)*(BETA1R(I) + 
     &              BETA3R(I)) - BETA2R(I)
        RHOR(I)   = DMAX1(SMLRHO,RHOR(I))
 
        VELR(I)   = VEL1R(I) + BETA3R(I) - BETA1R(I)
 
        UR(I)     = PR(I)/(GB - 1.D0)/RHOR(I)
 
        SCRCH1(I) = 1.D0 + UR(I) + PR(I)/RHOR(I)
 
        CSR(I)    = DSQRT(GB*PR(I)/RHOR(I)/SCRCH1(I))
 
45      CONTINUE
 
C     --------------
C     CONSERVED VARIABLES
C     --------------

      DO 60 I=1,NX+1

         SCRCH1(I) = 1.D0/DSQRT(1.D0 - VELL(I)*VELL(I))
         SCRCH2(I) = 1.D0 + UL(I) + PL(I)/RHOL(I)
         RL(I)     = RHOL(I)*SCRCH1(I)
         ML(I)     = RL(I)*SCRCH2(I)*SCRCH1(I)*VELL(I)
         EL(I)     = RL(I)*SCRCH2(I)*SCRCH1(I) - PL(I) - RL(I)

60       CONTINUE
         
      DO 70 I=1,NX+1

         SCRCH1(I) = 1.D0/DSQRT(1.D0 - VELR(I)*VELR(I))
         SCRCH2(I) = 1.D0 + UR(I) + PR(I)/RHOR(I)
         RR(I)     = RHOR(I)*SCRCH1(I)
         MR(I)     = RR(I)*SCRCH2(I)*SCRCH1(I)*VELR(I)
         ER(I)     = RR(I)*SCRCH2(I)*SCRCH1(I) - PR(I) - RR(I)

70       CONTINUE
 
      RETURN
      END

C     -------- 
CN    NAME: A V R G 1 D 
C     --------

CP    PURPOSE: 
CP    THIS SUBROUTINE CALCULATES AVERAGES OF QUANTITIES P,RHO,VEL, OVER
CP    THE PART OF THE DOMAIN OF DEPENDENCE FOR THE LAMBDA 
CP    CHARACTERISTIC OF RADM(I) FOR THE TIME INTERVAL (T(N),T(N+1)).
C

CC    COMMENTS: 
CC    THIS ROUTINE CLOSELY FOLLOWS THE ANALYTICAL DEVELOPMENTS DESCRIBED IN 
CC    MARTI & MUELLER, JCP, 1996

      SUBROUTINE AVRG1D(N,DT,LAMB,PM,PP,DP,P6,RHOM,RHOP,DRHO,
     &                  RHO6,VELM,VELP,DVEL,VEL6,PL,PR,RHOL,RHOR,
     &                  VELL,VELR)
 
      IMPLICIT NONE
 
      INCLUDE 'size'
 
C     --------
C     ARGUMENTS
C     --------

      INTEGER         N
 
      DOUBLEPRECISION DT
 
      DOUBLEPRECISION LAMB(-4:MN5)

      DOUBLEPRECISION PM(-4:MN5),PP(-4:MN5),DP(-4:MN5),P6(-4:MN5),
     &                RHOM(-4:MN5),RHOP(-4:MN5),DRHO(-4:MN5),
     &                RHO6(-4:MN5),VELM(-4:MN5),VELP(-4:MN5),
     &                DVEL(-4:MN5),VEL6(-4:MN5)

      DOUBLEPRECISION PL(-4:MN6),PR(-4:MN6),RHOL(-4:MN6),RHOR(-4:MN6),
     &                VELL(-4:MN6),VELR(-4:MN6)

C     ------
C     COMMON BLOCKS
C     ------

      DOUBLEPRECISION X(-4:MN5),XL(-4:MN5),XR(-4:MN5),DX(-4:MN5)
      COMMON /GRD/    X,XL,XR,DX

C     ----------
C     INTERNAL VARIABLES
C     ----------

      INTEGER         I
 
      DOUBLEPRECISION SCRCH1(-4:MN6)
 

      DO 10 I=0,N+1

        SCRCH1(I) = DMAX1(0.D0,DT*LAMB(I)/DX(I))

10      CONTINUE
 
      DO 20 I=1,N+1
        
        PL(I)   = PP(I-1)   - SCRCH1(I-1)/2.D0*
     &           (DP(I-1)   - (1.D0 - 2.D0*SCRCH1(I-1)/3.D0)*P6(I-1)  )
 
        RHOL(I) = RHOP(I-1) - SCRCH1(I-1)/2.D0*
     &           (DRHO(I-1) - (1.D0 - 2.D0*SCRCH1(I-1)/3.D0)*RHO6(I-1))
 
        VELL(I) = VELP(I-1) - SCRCH1(I-1)/2.D0*
     &           (DVEL(I-1) - (1.D0 - 2.D0*SCRCH1(I-1)/3.D0)*VEL6(I-1))

20      CONTINUE
 
      DO 30 I=0,N+1

        SCRCH1(I) = DMAX1(0.D0,-DT*LAMB(I)/DX(I))

30      CONTINUE
 
      DO 40 I=1,N+1
        
        PR(I)   = PM(I)   + SCRCH1(I)/2.D0*
     &           (DP(I)   + (1.D0 - 2.D0*SCRCH1(I)/3.D0)*P6(I)  )
 
        RHOR(I) = RHOM(I) + SCRCH1(I)/2.D0*
     &           (DRHO(I) + (1.D0 - 2.D0*SCRCH1(I)/3.D0)*RHO6(I))
 
        VELR(I) = VELM(I) + SCRCH1(I)/2.D0*
     &           (DVEL(I) + (1.D0 - 2.D0*SCRCH1(I)/3.D0)*VEL6(I))

40      CONTINUE

      RETURN
      END

C     -------- 
CN    NAME: N F L U X 
C     --------

CP    PURPOSE: 
CP    COMPUTES THE NUMERICAL FLUXES 
C

CC    COMMENTS: 
CC    COMPUTES THE NUMERICAL FLUXES FROM THE EXACT SOLUTION OF THE 
CC    RELATIVISTIC RIEMANN PROBLEM AS DESCRIBED IN MARTI AND MUELLER, JFM, 
CC    1994
      
      SUBROUTINE NFLUX(RHOL1,RHOR1,PL1,PR1,VELL1,VELR1,UL1,UR1,
     &                 CSL1,CSR1,FR,FM,FE)

      IMPLICIT NONE

C     -----------
C     ARGUMENTS
C     -----------

      DOUBLE PRECISION RHOL1, PL1, UL1, CSL1, VELL1, 
     &                 RHOR1, PR1, UR1, CSR1, VELR1 

      DOUBLE PRECISION FR, FM, FE

C     ------- 
C     COMMON BLOCKS 
C     -------

      DOUBLE PRECISION RHOL, PL, UL, HL, CSL, VELL, WL, 
     &                 RHOR, PR, UR, HR, CSR, VELR, WR  
      COMMON /STATES/  RHOL, PL, UL, HL, CSL, VELL, WL, 
     &                 RHOR, PR, UR, HR, CSR, VELR, WR

      DOUBLE PRECISION RHOLS, ULS, HLS, CSLS, VELLS, VSHOCKL 
      COMMON /LS/      RHOLS, ULS, HLS, CSLS, VELLS, VSHOCKL

      DOUBLE PRECISION RHORS, URS, HRS, CSRS, VELRS, VSHOCKR 
      COMMON /RS/      RHORS, URS, HRS, CSRS, VELRS, VSHOCKR

      DOUBLE PRECISION GAMMA 
      COMMON /ADIND/   GAMMA

C     --------- 
C     INTERNAL VARIABLES 
C     ---------

      INTEGER          ILOOP

      DOUBLE PRECISION TOL, PMIN, PMAX, DVEL1, DVEL2, CHECK

      DOUBLE PRECISION PS, VELS

      DOUBLE PRECISION RHOA, PA, VELA, UA

      DOUBLE PRECISION XI, XI1, XI2, XI3, XI4, XI5

C     ------------------------------ 
C     SPECIFIC ENTHALPY AND  
C     FLOW LORENTZ FACTORS IN THE INITIAL STATES 
C     ------------------------------ 

      RHOL = RHOL1
      RHOR = RHOR1

      PL   = PL1
      PR   = PR1

      UL   = UL1
      UR   = UR1

      VELL = VELL1
      VELR = VELR1

      CSL  = CSL1
      CSR  = CSR1

      HL  = 1.D0+UL+PL/RHOL 
      HR  = 1.D0+UR+PR/RHOR

      WL  = 1.D0/DSQRT(1.D0-VELL**2) 
      WR  = 1.D0/DSQRT(1.D0-VELR**2)

C     ------------- 
C     TOLERANCE FOR THE SOLUTION 
C     -------------

      TOL = 1.D-8

C

      ILOOP = 0

      PMIN  = (PL + PR)/2.D0 
      PMAX  = PMIN

 5    ILOOP = ILOOP + 1

      PMIN  = 0.5D0*MAX(PMIN,0.D0) 
      PMAX  = 2.D0*PMAX

      CALL GETDVEL(PMIN, DVEL1)

      CALL GETDVEL(PMAX, DVEL2)

      CHECK = DVEL1*DVEL2 
      IF (CHECK.GT.0.D0) GOTO 5

C     --------------------------- 
C     PRESSURE AND FLOW VELOCITY IN THE INTERMEDIATE STATES 
C     ---------------------------

      CALL GETP(PMIN, PMAX, TOL, PS)

      VELS = 0.5D0*(VELLS + VELRS)

C     --------------- 
C     SOLUTION ON THE NUMERICAL INTERFACE 
C     ---------------

C     ----------- 
C     POSITIONS OF THE WAVES 
C     -----------

      IF (PL.GE.PS) THEN

        XI1 = (VELL - CSL )/(1.D0 - VELL*CSL ) 
        XI2 = (VELS - CSLS)/(1.D0 - VELS*CSLS)

      ELSE

        XI1 = VSHOCKL
        XI2 = XI1

      END IF

        XI3 = VELS

      IF (PR.GE.PS) THEN

        XI4 = (VELS + CSRS)/(1.D0 + VELS*CSRS) 
        XI5 = (VELR + CSR )/(1.D0 + VELR*CSR )

      ELSE

        XI4 = VSHOCKR
        XI5 = XI4

      END IF

C     ---------- 
C     SOLUTION ON THE INTERFACE AT X = 0 (XI = 0)
C     ----------

      XI = 0.D0

      IF (XI1.GE.XI) THEN

        PA   = PL 
        RHOA = RHOL 
        VELA = VELL 
        UA   = UL

      ELSE IF (XI2.GE.XI) THEN

        CALL RAREF(XI,RHOL,CSL,VELL,'L',RHOA,PA,UA,VELA)

      ELSE IF (XI3.GE.XI) THEN

        PA   = PS 
        RHOA = RHOLS 
        VELA = VELS 
        UA   = ULS

      ELSE IF (XI4.GE.XI) THEN

        PA   = PS 
        RHOA = RHORS 
        VELA = VELS 
        UA   = URS

      ELSE IF (XI5.GE.XI) THEN

        CALL RAREF(XI,RHOR,CSR,VELR,'R',RHOA,PA,UA,VELA)

      ELSE

        PA   = PR 
        RHOA = RHOR 
        VELA = VELR 
        UA   = UR

      END IF

C     -----------
C     NUMERICAL FLUXES 
C     -----------

      FR = RHOA*VELA/DSQRT(1.D0 - VELA**2)  
      FM = RHOA*(1.D0 + UA + PA/RHOA)*VELA**2/(1.D0 - VELA**2) + PA  
      FE = RHOA*(1.D0 + UA + PA/RHOA)*VELA/(1.D0 - VELA**2) -
     &     RHOA*VELA/DSQRT(1.D0 - VELA**2)
   
      RETURN
      END

C     ---------- 
CN    NAME: G E T D V E L 
C     ----------

CP    PURPOSE: 
CP    COMPUTE THE DIFFERENCE IN FLOW SPEED BETWEEN LEFT AND RIGHT INTERMEDIATE 
CP    STATES FOR GIVEN LEFT AND RIGHT STATES AND PRESSURE 
C 

CC    COMMENTS:
CC    NONE
C 
      SUBROUTINE GETDVEL( P, DVEL )

      IMPLICIT NONE

C     ----- 
C     ARGUMENTS 
C     -----

      DOUBLEPRECISION P, DVEL

C     ------- 
C     COMMON BLOCKS 
C     -------

      DOUBLE PRECISION RHOLS,ULS,HLS,CSLS,VELLS,VSHOCKL 
      COMMON /LS/      RHOLS,ULS,HLS,CSLS,VELLS,VSHOCKL

      DOUBLE PRECISION RHORS,URS,HRS,CSRS,VELRS,VSHOCKR 
      COMMON /RS/      RHORS,URS,HRS,CSRS,VELRS,VSHOCKR

      DOUBLE PRECISION RHOL, PL, UL, HL, CSL, VELL, WL, 
     & RHOR, PR, UR, HR, CSR, VELR, WR 
      COMMON /STATES/  RHOL, PL, UL, HL, CSL, VELL, WL, 
     & RHOR, PR, UR, HR, CSR, VELR, WR

      DOUBLE PRECISION GAMMA 
      COMMON /ADIND/   GAMMA

C     ----- 
C     LEFT WAVE 
C     -----

      CALL GETVEL(P, RHOL, PL, UL,  HL,  CSL,  VELL,  WL, 'L', 
     & RHOLS,    ULS, HLS, CSLS, VELLS, VSHOCKL )

C     ----- 
C     RIGHT WAVE 
C     -----

      CALL GETVEL(P, RHOR, PR, UR,  HR,  CSR,  VELR,  WR, 'R', 
     & RHORS,    URS, HRS, CSRS, VELRS, VSHOCKR )

      DVEL = VELLS - VELRS

      RETURN 
      END

C     ------- 
CN    NAME: G E T P 
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
C
      SUBROUTINE GETP( PMIN, PMAX, TOL, PS )

      IMPLICIT NONE

C     ----- 
C     ARGUMENTS 
C     -----

      DOUBLEPRECISION PMIN, PMAX, TOL, PS

C     ------- 
C     COMMON BLOCKS 
C     -------

      DOUBLEPRECISION GAMMA 
      COMMON /ADIND/  GAMMA

      DOUBLEPRECISION RHOL, PL, UL, HL, CSL, VELL, WL, 
     & RHOR, PR, UR, HR, CSR, VELR, WR 
      COMMON /STATES/ RHOL, PL, UL, HL, CSL, VELL, WL, 
     & RHOR, PR, UR, HR, CSR, VELR, WR

C     --------- 
C     INTERNAL VARIABLES 
C     ---------

      DOUBLEPRECISION A, B, C, D, E, EPS, FA, FB, FC, TOL1, 
     & XM, P, Q, R, S

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
      CALL GETDVEL(A,FA) 
      CALL GETDVEL(B,FB)

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
      CALL GETDVEL(B,FB) 
      IF( (FB*(FC/DABS(FC))) .GT. 0.D0) GO TO 20 
      GO TO 30

C     -- 
C     DONE 
C     --

90    PS = B

      RETURN 
      END

C     --------- 
CN    NAME: G E T V E L 
C     ---------

CP    PURPOSE: 
CP    COMPUTE THE FLOW VELOCITY BEHIND A RAREFACTION OR SHOCK IN TERMS OF THE 
CP    POST-WAVE PRESSURE FOR A GIVEN STATE AHEAD THE WAVE IN A RELATIVISTIC 
CP    FLOW 
C 

CC    COMMENTS: 
CC    THIS ROUTINE CLOSELY FOLLOWS THE EXPRESSIONS IN MARTI AND MUELLER, 
CC    J. FLUID MECH., (1994)

      SUBROUTINE GETVEL( P, RHOA, PA, UA, HA, CSA, VELA, WA, S, 
     & RHO,      U,  H,  CS,  VEL,  VSHOCK )

      IMPLICIT NONE

C     ----- 
C     ARGUMENTS 
C     -----

      DOUBLE PRECISION P, RHOA, PA, UA, HA, CSA, VELA, WA  
      CHARACTER*1      S 
      DOUBLE PRECISION RHO, U, H, CS, VEL, VSHOCK

C     ------- 
C     COMMON BLOCKS 
C     -------

      DOUBLE PRECISION GAMMA 
      COMMON /ADIND/   GAMMA

C     --------- 
C     INTERNAL VARIABLES 
C     ---------

      DOUBLE PRECISION A, B, C, SIGN 
      DOUBLE PRECISION J, WSHOCK 
      DOUBLE PRECISION K, SQGL1

C     --------------- 
C     LEFT OR RIGHT PROPAGATING WAVE 
C     ---------------

      IF (S.EQ.'L') SIGN = -1.D0

      IF (S.EQ.'R') SIGN =  1.D0

C

      IF (P/PA - 1.D0.GT.1.D-10) THEN

C       --- 
C       SHOCK 
C       ---

        A  = 1.D0+(GAMMA-1.D0)*(PA-P)/GAMMA/P 
        B  = 1.D0-A 
        C  = HA*(PA-P)/RHOA-HA**2

C       ---------------- 
C       CHECK FOR UNPHYSICAL ENTHALPIES 
C       ----------------

        IF (C.GT.(B**2/4.D0/A)) STOP 
     & 'GETVEL: UNPHYSICAL SPECIFIC ENTHALPY IN INTERMEDIATE STATE'

C       ----------------------------- 
C       SPECIFIC ENTHALPY IN THE POST-WAVE STATE 
C       (FROM THE EQUATION OF STATE AND THE TAUB ADIABAT, 
C       EQ.(74), MM94) 
C       -----------------------------

        H = (-B+DSQRT(B**2-4.D0*A*C))/2.D0/A

C       --------------- 
C       DENSITY IN THE POST-WAVE STATE 
C       (FROM EQ.(73), MM94) 
C       ---------------

        RHO = GAMMA*P/(GAMMA-1.D0)/(H-1.D0)

C       ------------------------ 
C       SPECIFIC INTERNAL ENERGY IN THE POST-WAVE STATE 
C       (FROM THE EQUATION OF STATE) 
C       ------------------------

        U = P/(GAMMA-1.D0)/RHO

C       -------------------------- 
C       MASS FLUX ACROSS THE WAVE  
C       (FROM THE RANKINE-HUGONIOT RELATIONS, EQ.(71), MM94) 
C       --------------------------

        J = SIGN*DSQRT((P-PA)/(HA/RHOA-H/RHO))

C       ---------- 
C       SHOCK VELOCITY 
C       (FROM EQ.(86), MM94 
C       ----------

        A      = J**2+(RHOA*WA)**2 
        B      = -VELA*RHOA**2*WA**2 
        VSHOCK = (-B+SIGN*J**2*DSQRT(1.D0+RHOA**2/J**2))/A 
        WSHOCK = 1.D0/DSQRT(1.D0-VSHOCK**2)

C       ------------------- 
C       FLOW VELOCITY IN THE POST-SHOCK STATE 
C       (FROM EQ.(67), MM94) 
C       -------------------

        A = WSHOCK*(P-PA)/J+HA*WA*VELA 
        B = HA*WA+(P-PA)*(WSHOCK*VELA/J+1.D0/RHOA/WA)

        VEL = A/B

C       --------------------- 
C       LOCAL SOUND SPEED IN THE POST-SHOCK STATE 
C       (FROM THE EQUATION OF STATE) 
C       ---------------------

        CS = DSQRT(GAMMA*P/RHO/H)

      ELSE IF (P/PA - 1.D0.GT.0.D0) THEN
       
C       --------------
C       ALMOST CONSTANT INTERMEDIATE STATE
C       --------------

        RHO    = RHOA
        U      = UA
        H      = 1.D0 + U + P/RHO
        CS     = DSQRT(GAMMA*P/RHO/H)
        VEL    = VELA
        VSHOCK = VELA

      ELSE

C       ------ 
C       RAREFACTION 
C       ------

C       --------------------------- 
C       POLITROPIC CONSTANT OF THE GAS ACROSS THE RAREFACTION 
C       ---------------------------

        K = PA/RHOA**GAMMA

C       --------------- 
C       DENSITY BEHIND THE RAREFACTION 
C       ---------------

        RHO = (P/K)**(1.D0/GAMMA)

C       ------------------------ 
C       SPECIFIC INTERNAL ENERGY BEHIND THE RAREFACTION 
C       (FROM THE EQUATION OF STATE) 
C       ------------------------

        U = P/(GAMMA-1.D0)/RHO

C       -----------
C       SPECIFIC ENTHALPY
C       -----------

        H = 1.D0 + U + P/RHO
C       -------------------- 
C       LOCAL SOUND SPEED BEHIND THE RAREFACTION 
C       (FROM THE EQUATION OF STATE) 
C       --------------------

        CS = DSQRT(GAMMA*P/RHO/H)

C       ------------------ 
C       FLOW VELOCITY BEHIND THE RAREFACTION 
C       ------------------

        SQGL1 = DSQRT(GAMMA-1.D0) 
        A   = (1.D0+VELA)/(1.D0-VELA)* 
     & ((SQGL1+CSA)/(SQGL1-CSA)* 
     & (SQGL1-CS )/(SQGL1+CS ))**(-SIGN*2.D0/SQGL1)

        VEL = (A-1.D0)/(A+1.D0)

      END IF
      END

C     -------- 
CN    NAME: R A R E F 
C     --------

CP    PURPOSE: 
CP    COMPUTE THE FLOW STATE IN A RAREFACTION FOR GIVEN PRE-WAVE STATE 
C

CC    COMMENTS: 
CC    THIS ROUTINE CLOSELY FOLLOWS THE EXPRESSIONS IN MARTI AND MUELLER, 
CC    J. FLUID MECH., (1994)

      SUBROUTINE RAREF( XI, RHOA, CSA, VELA, S, RHO, P, U, VEL )

      IMPLICIT NONE

C     ----- 
C     ARGUMENTS 
C     -----

      DOUBLE PRECISION XI

      DOUBLE PRECISION RHOA, CSA, VELA

      CHARACTER        S

      DOUBLE PRECISION RHO, P, U, VEL

C     ------- 
C     COMMON BLOCKS 
C     -------

      DOUBLE PRECISION GAMMA 
      COMMON /ADIND/   GAMMA

C     --------- 
C     INTERNAL VARIABLES 
C     ---------

      DOUBLE PRECISION B, C, D, K, L, V, OCS2, FCS2, DFDCS2, CS2, SIGN

C     --------------- 
C     LEFT OR RIGHT PROPAGATING WAVE 
C     ---------------

      IF (S.EQ.'L') SIGN =  1.D0

      IF (S.EQ.'R') SIGN = -1.D0

      B    = DSQRT(GAMMA - 1.D0) 
      C    = (B + CSA)/(B - CSA) 
      D    = -SIGN*B/2.D0 
      K    = (1.D0 + XI)/(1.D0 - XI) 
      L    = C*K**D 
      V    = ((1.D0 - VELA)/(1.D0 + VELA))**D

      OCS2 = CSA

25    FCS2   = L*V*(1.D0 + SIGN*OCS2)**D*(OCS2 - B) + 
     & (1.D0 - SIGN*OCS2)**D*(OCS2 + B)

      DFDCS2 = L*V*(1.D0 + SIGN*OCS2)**D* 
     & (1.D0 + SIGN*D*(OCS2 - B)/(1.D0 + SIGN*OCS2)) + 
     & (1.D0 - SIGN*OCS2)**D*
     & (1.D0 - SIGN*D*(OCS2 + B)/(1.D0 - SIGN*OCS2))

      CS2 = OCS2 - FCS2/DFDCS2

      IF (ABS(CS2 - OCS2)/OCS2.GT.5.E-7)THEN 
        OCS2 = CS2 
        GOTO 25 
      END IF

      VEL = (XI + SIGN*CS2)/(1.D0 + SIGN*XI*CS2)

      RHO = RHOA*((CS2**2*(GAMMA - 1.D0 - CSA**2))/ 
     & (CSA**2*(GAMMA - 1.D0 - CS2**2))) 
     & **(1.D0/(GAMMA - 1.D0))

      P   = CS2**2*(GAMMA - 1.D0)*RHO/(GAMMA - 1.D0 - CS2**2)/GAMMA

      U   = P/(GAMMA - 1.D0)/RHO

      RETURN  
      END 

C     -------- 
CN    NAME: F I L N A M
C     --------

CP    PURPOSE: 
CP    CONSTRUCTS NEW FILENAMES FOR OUTPUT AND RESTART FILES
C

CC    COMMENTS: 
CC    NONE

      SUBROUTINE FILNAM

      IMPLICIT NONE

C     ---------
C     COMMON BLOCKS
C     ---------

      CHARACTER*7   OUTFIL
      CHARACTER*8   LABEL
      CHARACTER*4   BASENM
      CHARACTER*2   SUFFIX
      CHARACTER*1   SF1,SF2

      COMMON /CHRC/ LABEL,OUTFIL,BASENM,SUFFIX

C     ---------
C     INTERNAL VARIABLES
C     ---------

      INTEGER       ISF1,ISF2

      IF (SUFFIX(2:2).EQ.'z'.OR.SUFFIX(2:2).EQ.'Z') THEN
         SF1  = SUFFIX(1:1)
         ISF1 = ICHAR(SF1)
         SF2  = SUFFIX(2:2)
         ISF2 = ICHAR(SF2)

         ISF1 = ISF1 + 1
         ISF2 = ISF2 - 25
         SUFFIX(1:1) = CHAR(ISF1)
         SUFFIX(2:2) = CHAR(ISF2)
      ELSE
         SF2  = SUFFIX(2:2)
         ISF2 = ICHAR(SF2)
         ISF2 = ISF2 + 1
         SUFFIX(2:2) = CHAR(ISF2)
      END IF

      OUTFIL = BASENM // 'O' // SUFFIX

      RETURN
      END

C     -------- 
CN    NAME: G E T P R F Q
C     --------

CP    PURPOSE: 
CP    COMPUTE THE PRIMITIVE QUANTITIES FROM THE CONSERVED ONES  
C

CC    COMMENTS: 
CC    PRIMITIVE VARIABLES ARE OBTAINED BY SOLVING AN IMPLICIT EQUATION FOR 
CC    THE PRESSURE BY MEANS OF A NEWTON-RAPHSON METHOD. HARDWIRED FOR A 
CC    CONSTANT-GAMMA IDEAL GAS

      SUBROUTINE GETPRFQ(N,R,M,E,
     &                   VEL,W,RHO,U,P,H,CS,DPDRH,DPDU)
 
      IMPLICIT NONE

      INCLUDE 'size'

C     -------
C     ARGUMENTS
C     -------

      INTEGER N

      DOUBLEPRECISION P(-4:MN5),RHO(-4:MN5),VEL(-4:MN5),
     &                U(-4:MN5),CS(-4:MN5),W(-4:MN5),H(-4:MN5),
     &                DPDRH(-4:MN5),DPDU(-4:MN5)

      DOUBLEPRECISION R(-4:MN5),M(-4:MN5),E(-4:MN5)

C     ---------
C     COMMON BLOCKS
C     ---------

      INTEGER         NEND,NOUT,ITSTP,NX
      COMMON /INPTI/  NEND,NOUT,ITSTP,NX

      DOUBLEPRECISION TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,
     &                SMALLU,GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2
      COMMON /INPTF/  TMAX,TOUT,CFL,DTINI,SMALL,SMLRHO,SMALLP,
     &                SMALLU,GRIDLX,EPSILN,ETA1,ETA2,EPSLN,AK0,OMG1,OMG2

      DOUBLEPRECISION GB
      COMMON /ADIND/  GB

      DOUBLEPRECISION TIME,DT
      COMMON /ZEIT/   TIME,DT

C     ---------
C     INTERNAL VARIABLES
C     ---------

      INTEGER I,COUNT

      DOUBLEPRECISION MOMENT,VELCTY

      DOUBLEPRECISION PMIN(-4:MN5),PMAX,OP(-4:MN5)

      DOUBLEPRECISION FP,DFDP,ERRP

      DO 5 I=1,N

         R(I) = DMAX1(R(I),SMLRHO)
         E(I) = DMAX1(E(I),SMALLU)

 5       CONTINUE

      DO 9 I=1,N

         COUNT   = 0

         MOMENT  = DABS(M(I))
         PMIN(I) = DMAX1(MOMENT - E(I) - R(I) + MOMENT*SMALL,SMALLP)
         PMAX    = (GB-1.D0)*E(I)
         IF (PMIN(I).GT.PMAX) GOTO 990

         OP(I)      = 0.5D0*(PMIN(I)+PMAX)

 8       COUNT   = COUNT + 1
         OP(I)   = DMAX1(OP(I),PMIN(I))
         VELCTY  = MOMENT/(E(I) + R(I) + OP(I))
         W(I)    = 1.D0/DSQRT(1.D0 - VELCTY*VELCTY)
         FP      = (GB - 1.D0)*(E(I) + R(I)*(1.D0 - W(I))+
     &             OP(I)*(1.D0 - W(I)*W(I)))/W(I)/W(I) - OP(I)
         DFDP    = (GB - 1.D0)*VELCTY*VELCTY*
     &             (E(I) + R(I)*(1.D0 - W(I))+OP(I))/
     &             (E(I) + R(I) + OP(I)) - 1.D0

         P(I)    = DMAX1(OP(I) - FP/DFDP,PMIN(I))

         ERRP    = DABS(1.D0 - P(I)/OP(I))

         OP(I)      = P(I)

         IF (COUNT.GE.10000) GOTO 999

         IF (ERRP.GT.1.D-8) GOTO 8

         VEL(I)  = M(I)/(E(I)+R(I)+OP(I))
         IF (DABS(VEL(I)) .LT.SMALL*SMALL) VEL(I)  = 0.D0
 
         RHO(I)  = R(I)/W(I)

 9       CONTINUE

      DO 30 I=1,N

	 U(I) = P(I)/(GB - 1.D0)/RHO(I)

         IF (P(I).EQ.PMIN(I)) THEN
            WRITE(6,*) 'GETPRFQ: MINIMUM PRESSURE REACHED AT POINT:'
            WRITE(6,*) '         I = ', I,' T = ', TIME
         END IF

 30      CONTINUE

      CALL EOS(N,RHO,U,GB,P,H,CS,DPDRH,DPDU)

      GOTO 1000

 990  WRITE(6,*) 'GETPRFQ: NO PHYSICAL PRESSURE AVAILABLE'
      WRITE(6,*) '         T          = ', TIME
      WRITE(6,*) '         I          = ', I
      WRITE(6,*) '         R          = ', R(I),   ' MOMENT = ', MOMENT 
      WRITE(6,*) '         E          = ', E(I)
      WRITE(6,*) '         MOMENT-E-D = ', MOMENT - R(I) - E(I)
      WRITE(6,*) '         (GB-1)E    = ', (GB - 1.D0)*E(I) 
      STOP

 999  WRITE(6,*) 'GETPRFQ: NON CONVERGENCE IN PRESSURE'
      WRITE(6,*) '         T   = ', TIME
      WRITE(6,*) '         I   = ', I
      WRITE(6,*) '         P   = ', P(I),   ' PMIN = ', PMIN(I)
      WRITE(6,*) '         R   = ', R(I),   ' M    = ', M(I)
      WRITE(6,*) '         E   = ', E(I)
      WRITE(6,*) '         VEL = ', VEL(I)
      STOP

 1000 CONTINUE

      RETURN
      END
