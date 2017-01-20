MODULE IRKC_M

  IMPLICIT NONE

! Declare everything in the module to be PRIVATE unless specifically
! declared as PUBLIC.

  PRIVATE

  TYPE, PUBLIC :: IRKC_SOL
    INTEGER :: NEQN,NPDES
    DOUBLE PRECISION :: T,TEND,HLAST,RTOL,ATOL,HMAX,H0
    DOUBLE PRECISION, DIMENSION(:), POINTER :: Y,YLAST,YP,YPLAST
    INTEGER :: NFE,NFI,NSTEPS,NACCPT,NREJCT,NFESIG,MAXM,ERR_FLAG
    LOGICAL :: DONE,ONE_STEP,CONSTANT_J,STOP_ON_ERROR
  END TYPE IRKC_SOL

  ! A structure SOL of type IRKC_SOL contains information
  ! about the numerical solution.  Of most direct interest
  ! are the fields that hold the current value t of the
  ! independent variable and the NEQN components of the
  ! solution y(t).
  !
  !   SOL%T     -- t
  !
  !   SOL%Y(*)  -- y(t)  (vector of NEQN components)
  !
  ! The logical quantity SOL%DONE indicates whether
  ! the integration is complete.
  !
  ! After a step by IRKC to SOL%T of size SOL%HLAST, the
  ! solution structure SOL contains all the information
  ! needed by the IRKC_VAL function to evaluate an
  ! approximate solution YARG(*) for any argument ARG in
  ! the interval [SOL%T-SOL%HLAST, SOL%T] that has the
  ! same order of accuracy as SOL%Y(*).
  !
  ! Statistics about the integration are available as
  ! the integer fields
  !
  !   SOL%NFE    -- number of evaluations of F_E
  !   SOL%NFI    -- number of evaluations of F_I
  !   SOL%NSTEPS -- number of steps
  !   SOL%NACCPT -- number of accepted steps
  !   SOL%NREJCT -- number of rejected steps
  !   SOL%NFESIG -- number of evaluations of F_E used
  !                 in estimating the spectral radius
  !   SOL%MAXM   -- maximum number of stages used
  !
  ! They can be displayed conveniently with the auxiliary
  ! subroutine IRKC_STATS.

  INTERFACE IRKC_SET

    MODULE PROCEDURE SET1,SET2

  END INTERFACE


  ! Global variables within module
  DOUBLE PRECISION :: UROUND=EPSILON(1D0),SQRTU
  LOGICAL :: HAVE_SPCRAD,STOP_ON_ERROR
  INTEGER :: NFE,NFI,NSTEPS,NACCPT,NREJCT,NFESIG,MAXM,MMAX
  DOUBLE PRECISION :: TDIR,HMAX,HLAST,RTOL,ATOL,HMUS1

  ! Externally visible subroutines and functions:
  PUBLIC :: IRKC, IRKC_SET, IRKC_VAL, IRKC_STATS


CONTAINS

!*******************************************************
!***          PUBLIC SUBROUTINES/FUNCTIONS           ***
!*******************************************************

  SUBROUTINE IRKC(SOL,F_E,F_I,SPCRAD)
    TYPE(IRKC_SOL) :: SOL
    DOUBLE PRECISION, OPTIONAL :: SPCRAD
    EXTERNAL F_E,F_I,SPCRAD

    INTEGER :: NEQN,NPDES
    DOUBLE PRECISION :: T,TEND
    DOUBLE PRECISION, DIMENSION(SOL%NEQN) :: YNEW

    LOGICAL :: RHO_CONVG_FAIL,IT_CONVG_FAIL,ONE_STEP,CONSTANT_J,&
               LAST_STEP
    INTEGER :: M,I,IER,ERR_FLAG
    DOUBLE PRECISION :: JACNRM,EST,ERR,FAC,H,HMIN,TEMP1,TEMP2,TEMP3

    ! In the one-step mode, a few quantities need to be kept
    ! from one call to the next:
    INTEGER, SAVE :: NSTSIG
    DOUBLE PRECISION, SAVE :: HOLD,ERROLD,ABSH,SPRAD
    DOUBLE PRECISION, ALLOCATABLE, SAVE :: EV(:)
    LOGICAL, SAVE :: NEWSPC,JACATT

    IF (SOL%DONE) THEN
      PRINT *,' Cannot call IRKC after integration is done.'
      STOP
    END IF

    ERR_FLAG = SOL%ERR_FLAG
    IF ( (SOL%NFE > 0) .AND. (ERR_FLAG > 0) ) THEN
      PRINT *,' Cannot call IRKC after a fatal error.'
      STOP
    END IF

    ! Local variables for quantities in SOL
    NEQN = SOL%NEQN
    NPDES = SOL%NPDES
    T = SOL%T
    TEND = SOL%TEND
    TDIR = SIGN(1D0,TEND - T)
    HMAX = SOL%HMAX
    RTOL = SOL%RTOL
    ATOL = SOL%ATOL
    ONE_STEP = SOL%ONE_STEP
    CONSTANT_J = SOL%CONSTANT_J
    STOP_ON_ERROR = SOL%STOP_ON_ERROR
    HLAST = SOL%HLAST
    NFE = SOL%NFE
    NFI = SOL%NFI
    NSTEPS = SOL%NSTEPS
    NACCPT = SOL%NACCPT
    NREJCT = SOL%NREJCT
    NFESIG = SOL%NFESIG
    MAXM = SOL%MAXM

    ! Load a work array.
    YNEW = SOL%Y

    ! Initialize on the first call.
    IF (NFE == 0) THEN

      IF ((RTOL > 0.1D0) .OR. (RTOL < 10D0*UROUND)) THEN
        IF (STOP_ON_ERROR) THEN
          sol%t=-1d0;return    !CALL REPORT_ERRORS(1)
        ELSE
          SOL%ERR_FLAG = 1
          RETURN
        END IF
      END IF

      IF (ATOL <= 0D0 ) THEN
        IF (STOP_ON_ERROR) THEN
          sol%t=-1d0;return    !CALL REPORT_ERRORS(2)
        ELSE
          SOL%ERR_FLAG = 2
          RETURN
        END IF
      END IF

      HAVE_SPCRAD = PRESENT(SPCRAD)
      IF (.NOT. HAVE_SPCRAD) THEN
        ALLOCATE(EV(NEQN),STAT=IER)
        IF (IER /= 0) THEN
          IF (STOP_ON_ERROR) THEN
            sol%t=-1d0;return    !CALL REPORT_ERRORS(7)
          ELSE
            SOL%ERR_FLAG = 7
            RETURN
          END IF
        END IF
      END IF

      SQRTU = SQRT(UROUND)
      MMAX = NINT(SQRT(RTOL/(10D0*UROUND)))
      MMAX = MAX(MMAX,2)
      JACATT = .FALSE.
      NSTSIG = 0
      NEWSPC = .TRUE.

      CALL F_E(NEQN,T,YNEW,SOL%YP)
      NFE = NFE + 1
      ! Get maximum norm of Jacobian of F_I for computation
      ! of initial step size.
      CALL ADD_F_I(NEQN,T,YNEW,SOL%YP,NPDES,F_I,JACNRM)

      HMIN = 10D0*UROUND*MAX(ABS(T),ABS(TEND))

    END IF

    ! Start of loop for taking a step:
    TAKE_A_STEP:  DO

      !  Estimate the spectral radius of the Jacobian
      !  when NEWSPC = .TRUE. A convergence failure
      !  in RHO is reported by RHO_CONVG_FAIL.
      IF (NEWSPC) THEN
        IF (HAVE_SPCRAD) THEN
          SPRAD = SPCRAD(NEQN,T,YNEW)
        ELSE
          CALL RHO(F_E,NEQN,T,YNEW,SOL%YP,SOL%YLAST,SOL%YPLAST,EV,&
                   SPRAD,RHO_CONVG_FAIL)
          IF (RHO_CONVG_FAIL) THEN
            IF (STOP_ON_ERROR) THEN
              sol%t=-1d0;return    !CALL REPORT_ERRORS(4)
            ELSE
              SOL%ERR_FLAG = 4
              RETURN
            END IF
          END IF
          ! The value SOL%YP returned by RHO corresponds to F_E only.
          CALL ADD_F_I(NEQN,T,YNEW,SOL%YP,NPDES,F_I)
        END IF
        JACATT = .TRUE.
      END IF

      !  Compute an initial step size.
      IF (NSTEPS == 0) THEN
        IF (SOL%H0 > 0D0) THEN
          ABSH = MIN(SOL%H0,HMAX)
        ELSE
          ABSH = HMAX
          IF (MAX(SPRAD,JACNRM)*ABSH > 1D0) THEN
            ABSH = 1D0/MAX(SPRAD,JACNRM)
          END IF
        END IF
        ABSH = MAX(ABSH,HMIN)
        SOL%YLAST = YNEW + ABSH*SOL%YP
        CALL F_E(NEQN,T+ABSH,SOL%YLAST,SOL%YPLAST)
        NFE = NFE + 1
        CALL ADD_F_I(NEQN,T+ABSH,SOL%YLAST,SOL%YPLAST,NPDES,F_I)

        EST = RMSNORM( (SOL%YPLAST - SOL%YP) / (ATOL + RTOL*ABS(YNEW)) )
        EST = ABSH*SQRT(EST/NEQN)

        IF (0.1D0*ABSH < HMAX*SQRT(EST)) THEN
          ABSH = MAX(0.1D0*ABSH/SQRT(EST),HMIN)
        ELSE
          ABSH = HMAX
        END IF

      END IF

      !  Adjust step size and determine number of stages M.
      LAST_STEP = .FALSE.
      IF (1.1D0*ABSH >= ABS(TEND - T)) THEN
        ABSH = ABS(TEND - T)
        LAST_STEP = .TRUE.
      END IF
      M = 1 + INT(SQRT(1.54D0*ABSH*SPRAD + 1D0))
      !  Limit M to MMAX to control the growth of roundoff error.
      IF (M > MMAX) THEN
        M = MMAX
        ABSH = (M**2 - 1)/(1.54D0*SPRAD)
        LAST_STEP = .FALSE.
      END IF
      MAXM = MAX(M,MAXM)

      ! A tentative solution at T+H is returned in
      ! SOL%Y and its slope is evaluated in SOL%YLAST.
      H = TDIR*ABSH
      HMIN = 10D0*UROUND*MAX(ABS(T),ABS(T+H))
      CALL STEP(NEQN,F_E,NPDES,F_I,T,YNEW,&
                H,M,SOL%Y,SOL%YLAST,SOL%YPLAST,&
                IT_CONVG_FAIL,ERR_FLAG)
      IF (ERR_FLAG > 0) THEN
        IF (STOP_ON_ERROR) THEN
          sol%t=-1d0;return    !CALL REPORT_ERRORS(ERR_FLAG)
        ELSE
          SOL%ERR_FLAG = ERR_FLAG
          RETURN
        END IF
      END IF

      IF (IT_CONVG_FAIL) THEN
        ERR = 10D0   ! Dummy value to trigger failed step.
      ELSE
        ! Here YNEW is y(T) and SOL%Y is the tentative y(T+H).
        ! Reuse of storage has y'(T) stored in SOL%YP and
        ! y'(T+H) is to be stored in SOL%YLAST.
        CALL F_E(NEQN,T+H,SOL%Y,SOL%YLAST)
        NFE = NFE + 1
        CALL ERR_EST(NEQN,NPDES,T,H,YNEW,&
                     SOL%YP,SOL%Y,SOL%YLAST,F_I,ERR,ERR_FLAG)
        IF (ERR_FLAG > 0) THEN
          IF (STOP_ON_ERROR) THEN
            sol%t=-1d0;return    !CALL REPORT_ERRORS(ERR_FLAG)
          ELSE
            SOL%ERR_FLAG = ERR_FLAG
            RETURN
          END IF
        END IF
      END IF
      NSTEPS = NSTEPS + 1

      IF (ERR > 1D0) THEN
        ! Step is rejected.
        NREJCT = NREJCT + 1
        IF (IT_CONVG_FAIL) THEN
          ABSH = ABSH/2D0
        ELSE
          ABSH = 0.8D0*ABSH/SQRT(ERR)
        END IF
        IF (ABSH < HMIN) THEN
          IF (STOP_ON_ERROR) THEN
            sol%t=-1d0;return    !CALL REPORT_ERRORS(3)
          ELSE
            SOL%ERR_FLAG = 3
            RETURN
          END IF
        ELSE
          NEWSPC = .NOT. JACATT
          CYCLE TAKE_A_STEP
        END IF
      END IF

      !  Step is accepted.
      NACCPT = NACCPT + 1
      T = T + H
      JACATT = CONSTANT_J
      NSTSIG = MOD(NSTSIG+1,25)
      NEWSPC = .FALSE.
      IF (HAVE_SPCRAD .OR. (NSTSIG == 0)) THEN
        NEWSPC = .NOT. JACATT
      END IF

      ! Update stored quantities.
      HLAST = H
      DO I = 1,NEQN
        TEMP3 = SOL%YLAST(I)
        SOL%YLAST(I) = YNEW(I)
        YNEW(I) = SOL%Y(I)
        SOL%YPLAST(I) = SOL%YP(I)
        SOL%YP(I) = TEMP3
      END DO

      ! Compute new step size.
      FAC = 10D0
      IF (NACCPT == 1) THEN
        TEMP2 = SQRT(ERR)
        IF (0.8D0 < FAC*TEMP2) THEN
          FAC = 0.8D0/TEMP2
        END IF
      ELSE
        TEMP1 = 0.8D0*ABSH*SQRT(ERROLD)
        TEMP2 = ABS(HOLD)*ERR
        IF (TEMP1 < FAC*TEMP2) THEN
          FAC = TEMP1/TEMP2
        END IF
      END IF
      ABSH = MAX(0.1D0,FAC)*ABSH
      ABSH = MAX(HMIN,MIN(HMAX,ABSH))
      ERROLD = ERR
      HOLD = H
      H = TDIR*ABSH

      IF (LAST_STEP .OR. ONE_STEP) THEN
        SOL%DONE = LAST_STEP
        EXIT TAKE_A_STEP
      END IF

      ! End of loop for taking a step:
    END DO TAKE_A_STEP


    ! Load working variables into SOL for return
    SOL%T      = T
    SOL%HLAST  = HLAST
    SOL%NFE    = NFE
    SOL%NFI    = NFI
    SOL%NSTEPS = NSTEPS
    SOL%NACCPT = NACCPT
    SOL%NREJCT = NREJCT
    SOL%NFESIG = NFESIG
    SOL%MAXM   = MAXM

    IF (SOL%DONE .AND. (.NOT. HAVE_SPCRAD) ) THEN
      DEALLOCATE(EV,STAT=IER)
      IF (IER /= 0) THEN
        IF (STOP_ON_ERROR) THEN
          sol%t=-1d0;return    !CALL REPORT_ERRORS(7)
        ELSE
          SOL%ERR_FLAG = 7
          RETURN
        END IF
      END IF
    END IF

  END SUBROUTINE IRKC


!--------------BEGIN GENERIC FUNCTION IRKC_SET--------------

  FUNCTION SET1(T0,Y0,TEND,RE,AE,ONE_STEP,CONSTANT_J,&
                    STOP_ON_ERROR,NPDES,HMAX,H0) RESULT(SOL)

    DOUBLE PRECISION :: T0,Y0(:),TEND
    INTEGER, OPTIONAL :: NPDES
    DOUBLE PRECISION, OPTIONAL :: RE,AE,HMAX,H0
    LOGICAL, OPTIONAL :: ONE_STEP,CONSTANT_J,STOP_ON_ERROR

    TYPE(IRKC_SOL) :: SOL
    INTEGER :: NEQN,IER

    SOL%T = T0
    SOL%TEND = TEND

    NEQN = SIZE(Y0)
    SOL%NEQN = NEQN
    ALLOCATE(SOL%Y(NEQN),SOL%YP(NEQN),SOL%YLAST(NEQN),&
             SOL%YPLAST(NEQN),STAT=IER)
    IF (IER /= 0) THEN
      PRINT *,'A storage allocation error occurred in IRKC_SET.'
      STOP
    END IF
    SOL%Y = Y0
    SOL%YP = 0D0
    SOL%YLAST = 0D0
    SOL%YPLAST = 0D0

    IF (PRESENT(RE)) THEN
      SOL%RTOL = RE
    ELSE
      SOL%RTOL = 1D-2
    END IF

    IF (PRESENT(AE)) THEN
      SOL%ATOL = AE
    ELSE
      SOL%ATOL = 1D-3
    END IF

    IF (PRESENT(ONE_STEP)) THEN
      SOL%ONE_STEP = ONE_STEP
    ELSE
      SOL%ONE_STEP = .TRUE.
    END IF

    IF (PRESENT(CONSTANT_J)) THEN
      SOL%CONSTANT_J = CONSTANT_J
    ELSE
      SOL%CONSTANT_J = .FALSE.
    END IF

    IF (PRESENT(STOP_ON_ERROR)) THEN
      SOL%STOP_ON_ERROR = STOP_ON_ERROR
    ELSE
      SOL%STOP_ON_ERROR = .TRUE.
    END IF
    SOL%ERR_FLAG = 0

    IF (PRESENT(NPDES)) THEN
      SOL%NPDES = NPDES
    ELSE
      SOL%NPDES = 1
    END IF

    IF (PRESENT(HMAX)) THEN
      SOL%HMAX = HMAX
    ELSE
      SOL%HMAX = ABS(TEND - T0)
    END IF

    IF (PRESENT(H0)) THEN
      SOL%H0 = ABS(H0)
    ELSE
      SOL%H0 = -1D0
    END IF

    SOL%NFE = 0
    SOL%NFI = 0
    SOL%NSTEPS = 0
    SOL%NACCPT = 0
    SOL%NREJCT = 0
    SOL%NFESIG = 0
    SOL%MAXM = 0
    SOL%HLAST = 0D0
    SOL%DONE = .FALSE.

  END FUNCTION SET1

  FUNCTION SET2(SOLIN,TEND,RE,AE,ONE_STEP,CONSTANT_J,&
                    STOP_ON_ERROR,NPDES,HMAX,H0) RESULT(SOL)

    TYPE(IRKC_SOL) :: SOLIN,SOL
    DOUBLE PRECISION :: TEND
    DOUBLE PRECISION, OPTIONAL :: RE,AE,HMAX,H0
    LOGICAL, OPTIONAL :: ONE_STEP,CONSTANT_J,STOP_ON_ERROR
    INTEGER, OPTIONAL :: NPDES

    IF (PRESENT(NPDES)) THEN
      IF (NPDES /= SOLIN%NPDES) THEN
        PRINT *,' Cannot change NPDES.'
        STOP
      END IF
    END IF

    ! Is this a new TEND?
    IF (ABS(TEND - SOLIN%TEND) > 0D0) THEN
      ! Is this a change of direction?
      IF (SIGN(1D0,TEND - SOLIN%TEND)*SOLIN%HLAST < 0) THEN
        PRINT *,' Cannot change direction without restarting.'
        STOP
      END IF
      SOLIN%TEND = TEND ! Reset in SOLIN for later copy to SOL.
    END IF

    SOL = SOLIN

    IF (PRESENT(RE)) SOL%RTOL = RE

    IF (PRESENT(AE)) SOL%ATOL = AE

    IF (PRESENT(ONE_STEP)) SOL%ONE_STEP = ONE_STEP

    IF (PRESENT(CONSTANT_J)) SOL%CONSTANT_J = CONSTANT_J

    IF (PRESENT(STOP_ON_ERROR)) SOL%STOP_ON_ERROR = STOP_ON_ERROR

    IF (PRESENT(HMAX)) SOL%HMAX = HMAX

    IF (PRESENT(H0)) SOL%H0 = ABS(H0)

    SOL%DONE = .FALSE.

  END FUNCTION SET2

!----------------END GENERIC FUNCTION IRKC_SET--------------


  FUNCTION IRKC_VAL(SOL,ARG) RESULT(YARG)
  !  IRKC_VAL is used to compute approximate solutions at specific t
  !  and to compute cheaply the large number of approximations that
  !  may be needed for plotting or locating when events occur.

  !  After a step by IRKC to SOL%T of size SOL%HLAST, the solution
  !  structure SOL contains all the information needed by IRKC_VAL to
  !  evaluate an approximate solution YARG(*) for any argument ARG
  !  in the interval [SOL%T-SOL%HLAST, SOL%T] that has the same order
  !  of accuracy as SOL%Y(*).
    TYPE(IRKC_SOL) :: SOL
    DOUBLE PRECISION :: ARG
    DOUBLE PRECISION, DIMENSION(SOL%NEQN) :: YARG
    DOUBLE PRECISION :: TLAST,A1,A2,B1,B2,S

    TLAST = SOL%T - SOL%HLAST
    S  = (ARG - TLAST) / SOL%HLAST
    A1 = (1D0 + 2D0*S)*(S - 1D0)**2
    A2 = (3D0 - 2D0*S)*S**2
    B1 = SOL%HLAST*S*(S - 1D0)**2
    B2 = SOL%HLAST*(S - 1D0)*S**2

    YARG = A1*SOL%YLAST + A2*SOL%Y + B1*SOL%YPLAST + B2*SOL%YP

  END FUNCTION IRKC_VAL


  SUBROUTINE IRKC_STATS(SOL)
    TYPE(IRKC_SOL) :: SOL
    INTEGER :: ITEMP

    ITEMP = DBLE(SOL%NFI)/( DBLE(SOL%NEQN)/DBLE(SOL%NPDES) )

    PRINT *,' '
    PRINT *,' Integration statistics: '
    PRINT *,' '
    PRINT *,' Evaluations of F_E                      = ',  &
            SOL%NFE
    PRINT *,' Evaluations of F_I / (NEQN/NPDES)       = ',  &
            ITEMP
    PRINT *,' Number of steps                         = ',  &
            SOL%NSTEPS
    PRINT *,' Number of accepted steps                = ',  &
            SOL%NACCPT
    PRINT *,' Number of rejected steps                = ',  &
            SOL%NREJCT
    PRINT *,' Evaluations of F_E for spectral radius  = ',  &
            SOL%NFESIG
    PRINT *,' Maximum number of stages used           = ',  &
            SOL%MAXM
    PRINT *,' '

  END SUBROUTINE IRKC_STATS


!*******************************************************
!***         PRIVATE SUBROUTINES/FUNCTIONS           ***
!*******************************************************

  SUBROUTINE ADD_F_I(NEQN,T,Y,YP,NPDES,F_I,JACNRM)
  ! Evaluate the terms to be handled implicitly
  ! and add in blocks of size NPDEs. Optionally
  ! compute the maximum norm of the Jacobian.
    INTEGER :: NEQN,NPDES,I
    DOUBLE PRECISION :: T,Y(NEQN),YP(NEQN)
    DOUBLE PRECISION, OPTIONAL :: JACNRM
    DOUBLE PRECISION :: TEMP(NPDES),JAC(NPDES,NPDES)
    EXTERNAL F_I

    IF (.NOT. PRESENT(JACNRM)) THEN
      DO I = 1,NEQN,NPDES
        CALL F_I(I,NPDES,T,Y(I:I+NPDES-1),TEMP,.FALSE.,JAC)
        NFI = NFI + 1
        YP(I:I+NPDES-1) = YP(I:I+NPDES-1) + TEMP
      END DO
    ELSE
      JACNRM = 0D0
      DO I = 1,NEQN,NPDES
        CALL F_I(I,NPDES,T,Y(I:I+NPDES-1),TEMP,.TRUE.,JAC)
        NFI = NFI + 1
        YP(I:I+NPDES-1) = YP(I:I+NPDES-1) + TEMP
        ! Compute maximum norm of the matrix JAC.
        JACNRM = MAX(JACNRM,MAXVAL( SUM(ABS(JAC),DIM=2) ))
      END DO
    END IF

  END SUBROUTINE ADD_F_I


  SUBROUTINE STEP(NEQN,F_E,NPDES,F_I,T,YN,H,M,Y,YJM1,YJM2,&
                  IT_CONVG_FAIL,ERR_FLAG)
  ! Take an implicit step of size H from T to T+H
  ! to get Y(*).
    INTEGER :: NEQN,M,NPDES,ERR_FLAG
    DOUBLE PRECISION :: T,H
    DOUBLE PRECISION, DIMENSION(NEQN) :: YN,Y,YJM1,YJM2
    LOGICAL :: IT_CONVG_FAIL
    EXTERNAL F_E,F_I

    INTEGER :: J,I
    DOUBLE PRECISION :: AJM1,ARG,BJ,BJM1,BJM2,DZJ,DZJM1,&
          DZJM2,D2ZJ,D2ZJM1,D2ZJM2,MU,MUS,NU,TEMP1,TEMP2,THJ,&
          THJM1,THJM2,W0,W1,ZJ,ZJM1,ZJM2,GAMMA,MUS1
    DOUBLE PRECISION :: RHS(NEQN),DYE_0(NEQN),DYE_JM1(NEQN),&
                        DYI_0(NEQN),DYI(NEQN,3),JAC(NPDES,NPDES)
    ! Swap columns of DYI by swapping these pointers:
    INTEGER :: PY,PY_JM1,PY_JM2,ITEMP
    PY     = 1
    PY_JM1 = 2
    PY_JM2 = 3

    ! Flag to report convergence failure.
    IT_CONVG_FAIL = .FALSE.

    W0 = 1D0 + 2D0/(13D0*M**2)
    TEMP1 = W0**2 - 1D0
    TEMP2 = SQRT(TEMP1)
    ARG = M*LOG(W0 + TEMP2)
    W1 = SINH(ARG)*TEMP1 / (COSH(ARG)*M*TEMP2 - W0*SINH(ARG))
    BJM1 = 1D0/W0
    BJM2 = 1D0/(2D0*W0)**2

    ! Evaluate the first stage.
    ! FN of explicit step is recomputed because
    ! in calling program it includes the implicit
    ! part which we need here in its own right.
    CALL F_E(NEQN,T,YN,DYE_0)
    NFE = NFE + 1
    DO I = 1,NEQN,NPDES
      CALL F_I(I,NPDES,T,YN(I:I+NPDES-1),&
               DYI_0(I:I+NPDES-1),.FALSE.,JAC)
      NFI = NFI + 1
    END DO
    DYI(:,PY_JM2) = DYI_0        ! Initialize for j = 2

    YJM2 = YN
    MUS  = W1*BJM1
    MUS1 = MUS
    HMUS1 = H*MUS1         ! Global variable
    RHS = YN + HMUS1*DYE_0
    YJM1 = YN
    CALL STAGE(NEQN,NPDES,F_I,T+HMUS1,RHS,YJM1,&
               DYI(:,PY_JM1),IT_CONVG_FAIL,ERR_FLAG)
    IF (IT_CONVG_FAIL .OR. ERR_FLAG > 0) RETURN

    THJM2  = 0D0
    THJM1  = MUS
    ZJM1   = W0
    ZJM2   = 1D0
    DZJM1  = 1D0
    DZJM2  = 0D0
    D2ZJM1 = 0D0
    D2ZJM2 = 0D0

    ! Evaluate stages j = 2,...,m.
    DO J = 2, M
      ZJ    =   2D0*W0*ZJM1 - ZJM2
      DZJ   =   2D0*W0*DZJM1 - DZJM2 + 2D0*ZJM1
      D2ZJ  =   2D0*W0*D2ZJM1 - D2ZJM2 + 4D0*DZJM1
      BJ    =   D2ZJ/DZJ**2
      AJM1  =   1D0 - ZJM1*BJM1
      MU    =   2D0*W0*BJ/BJM1
      NU    = - BJ/BJM2
      MUS   =   MU*W1/W0
      GAMMA = - AJM1*MUS

      CALL F_E(NEQN,T + H*THJM1,YJM1,DYE_JM1)
      NFE = NFE + 1
      RHS = (1D0 - MU - NU)*YN + MU*YJM1 + NU*YJM2 + &
             H*MUS*DYE_JM1 + H*GAMMA*DYE_0 + &
             (GAMMA - (1D0 - MU - NU)*MUS1)*H*DYI_0 &
             - NU*HMUS1*DYI(:,PY_JM2)
      Y = YJM1   ! Guess for next stage
      THJ = MU*THJM1 + NU*THJM2 + MUS*(1D0 - AJM1)
      CALL STAGE(NEQN,NPDES,F_I,T+H*THJ,&
                 RHS,Y,DYI(:,PY),IT_CONVG_FAIL,ERR_FLAG)
      IF (IT_CONVG_FAIL .OR. ERR_FLAG > 0) RETURN

      ! Shift the data for the next stage.
      IF (J < M) THEN
        YJM2   = YJM1
        YJM1   = Y

        ! Swap column pointers.
        ITEMP  = PY_JM2
        PY_JM2 = PY_JM1
        PY_JM1 = PY
        PY     = ITEMP

        THJM2  = THJM1
        THJM1  = THJ
        BJM2   = BJM1
        BJM1   = BJ
        ZJM2   = ZJM1
        ZJM1   = ZJ
        DZJM2  = DZJM1
        DZJM1  = DZJ
        D2ZJM2 = D2ZJM1
        D2ZJM1 = D2ZJ
      END IF

    END DO

  END SUBROUTINE STEP


  SUBROUTINE STAGE(NEQN,NPDES,F_I,TNJ,RHS,YS,DYS,&
                   IT_CONVG_FAIL,ERR_FLAG)
    INTEGER :: NEQN,NPDES,ERR_FLAG
    DOUBLE PRECISION :: TNJ,RHS(NEQN),YS(NEQN),DYS(NEQN)
    LOGICAL :: IT_CONVG_FAIL
    EXTERNAL  F_I

    INTEGER, PARAMETER :: ITMAX=10
    INTEGER :: I,M,INFO,IPVT(NPDES)
    DOUBLE PRECISION :: ERR,W(NPDES),RES(NPDES),WT(NPDES),&
                        ITMAT(NPDES,NPDES)

    IT_CONVG_FAIL = .FALSE.
    ! Compute implicit stage a block at at time
    DO I = 1,NEQN,NPDES
      ! Initialize guess for block
      W = YS(I:I+NPDES-1)

      ! Form and factor iteration matrix.
      CALL F_I(I,NPDES,TNJ,W,DYS(I:I+NPDES-1),.TRUE.,ITMAT)
      NFI = NFI + 1
      ! HMUS1 = H*MUS1 is a global variable.
      ITMAT = EYE(NPDES) - HMUS1*ITMAT
      CALL FACTOR(ITMAT,NPDES,INFO,IPVT)
      IF (INFO > 0) THEN
        ERR_FLAG = 5
        RETURN
      END IF

      DO M = 1,ITMAX
        ! Compute the residual of current iterate W.
        RES = RHS(I:I+NPDES-1) + HMUS1*DYS(I:I+NPDES-1) - W
        ! Compute the correction.
        RES = SOLVE(ITMAT,NPDES,IPVT,RES)

        WT = ABS(W)
        W = W + RES
        WT = MAX(WT,ABS(W))
        ! Prepare for another iteration;
        ! place result for slope in DYS.
        CALL F_I(I,NPDES,TNJ,W,DYS(I:I+NPDES-1),.FALSE.,ITMAT)
        NFI = NFI + 1

        ! Test for convergence
        ERR = RMSNORM( RES/(ATOL + RTOL*WT) )
        IF (ERR <= 0.5D0) THEN
          EXIT
        ELSEIF (M == ITMAX) THEN
          IT_CONVG_FAIL = .TRUE.
          RETURN
        END IF

      END DO

      ! Place result for solution in YS.
      YS(I:I+NPDES-1) = W

    END DO

  END SUBROUTINE STAGE


  SUBROUTINE ERR_EST(NEQN,NPDES,T,H,YNEW,YP,Y,YLAST,&
                                   F_I,ERR,ERR_FLAG)
  ! Here YNEW is y(T) and Y is the tentative y(T+H). Reuse of
  ! storage has y'(T) stored in YP and y'(T+H) is to be returned
  ! in YLAST.
    INTEGER :: NEQN,NPDES,ERR_FLAG,I,INFO
    DOUBLE PRECISION :: T,H,ERR
    DOUBLE PRECISION, DIMENSION(NEQN) :: YNEW,YP,Y,YLAST
    EXTERNAL F_I

    DOUBLE PRECISION :: TEMP(NPDES),EST(NPDES),FMAT(NPDES,NPDES)
    INTEGER :: IPVT(NPDES)

    ! The error estimate is formed, filtered, weighted, and
    ! the squares of the components are summed in blocks.
    ! After all blocks are processed, the RMS norm is formed.

    ! Form y'(T+H) in YLAST. The input YLAST contains F_E.
    ! Here F_I is evaluated in blocks and added.  The input
    ! YP holds y'(T).  However, we also need F_I at y(T),
    ! so it is evaluated in blocks and used in the error
    ! estimate in blocks.

    ERR = 0D0
    DO I = 1,NEQN,NPDES

      CALL F_I(I,NPDES,T+H,Y(I:I+NPDES-1),TEMP,.FALSE.,FMAT)
      NFI = NFI + 1
      ! Form a block of y'(T+H) in YLAST.
      YLAST(I:I+NPDES-1) = YLAST(I:I+NPDES-1) + TEMP
      EST = 0.5D0*H*(YLAST(I:I+NPDES-1) - YP(I:I+NPDES-1)) &
            + HMUS1*TEMP
      ! Get Jacobian along with F_I.
      CALL F_I(I,NPDES,T,YNEW(I:I+NPDES-1),TEMP,.TRUE.,FMAT)
      NFI = NFI + 1
      ! Complete computation of the unfiltered estimate.
      EST = EST - HMUS1*TEMP

      ! Form and factor the filter.
      FMAT = EYE(NPDES) - H*FMAT
      CALL FACTOR(FMAT,NPDES,INFO,IPVT)
      IF (INFO > 0) THEN
        ERR_FLAG = 6
        RETURN
      END IF
      ! Filter the estimate.
      EST = SOLVE(FMAT,NPDES,IPVT,EST)

      ! Weight the components of the estimated error.
      TEMP = EST / ( ATOL + &
             RTOL*MAX(ABS(Y(I:I+NPDES-1)),ABS(YNEW(I:I+NPDES-1))) )
      ! Sum the squares of the weighted error estimates.
      ERR = ERR + DOT_PRODUCT(TEMP,TEMP)

    END DO
    ! Form the RMS norm
    ERR = SQRT(ERR/NEQN)

  END SUBROUTINE ERR_EST


  SUBROUTINE RHO(F_E,NEQN,T,YN,FN,V,FV,EV,SPCRAD,RHO_CONVG_FAIL)
  ! RHO attempts to compute a close upper bound, SPCRAD, on
  ! the spectral radius of the Jacobian matrix using a nonlinear
  ! power method.  A convergence failure is reported by
  ! RHO_CONVG_FAIL = .TRUE.
    EXTERNAL F_E
    INTEGER :: NEQN
    DOUBLE PRECISION :: T,SPCRAD
    DOUBLE PRECISION, DIMENSION(NEQN) :: YN,FN,V,FV,EV
    LOGICAL :: RHO_CONVG_FAIL

    INTEGER, PARAMETER :: ITMAX=50
    INTEGER :: ITER,INDEX
    DOUBLE PRECISION :: YNRM,SIGMA,SIGMAL,DYNRM,DFNRM,VNRM,SMALL

    ! HMAX is a global variable.  SPCRAD smaller than 1/HMAX are
    ! not interesting because they do not constrain the step size.
    SMALL = 1D0/HMAX
    RHO_CONVG_FAIL = .FALSE.

    ! The initial slope is used as guess when NSTEPS = 0 and
    ! thereafter the last computed eigenvector.  Some care
    ! is needed to deal with special cases. Approximations to
    ! the eigenvector are normalized so that their Euclidean
    ! norm has the constant value DYNRM.

    ! It is the spectral radius of the Jacobian of F_E that is
    ! computed, but the input value FN includes the effect of
    ! the implicit term, so FN must be recomputed here.
    CALL F_E(NEQN,T,YN,FN)
    NFESIG = NFESIG + 1

    YNRM = ENORM(YN)
    IF (NSTEPS == 0) THEN
      V = FN
    ELSE
      V = EV
    END IF
    VNRM = ENORM(V)
    IF (ABS(YNRM) > 0D0 .AND. ABS(VNRM) > 0D0) THEN
      DYNRM = YNRM*SQRTU
      V = YN + V*(DYNRM/VNRM)
    ELSEIF (ABS(YNRM) > 0D0) THEN
      DYNRM = YNRM*SQRTU
      V = YN*(1D0 + SQRTU)
    ELSEIF (ABS(VNRM) > 0D0) THEN
      DYNRM = UROUND
      V = V*(DYNRM/VNRM)
    ELSE
      DYNRM = UROUND
      V = DYNRM
    END IF

    ! Now iterate with a nonlinear power method.
    SIGMA = 0D0
    DO ITER = 1,ITMAX

      CALL F_E(NEQN,T,V,FV)
      NFESIG = NFESIG + 1
      DFNRM = ENORM(FV - FN)
      SIGMAL = SIGMA
      SIGMA = DFNRM/DYNRM
      ! SPCRAD is a little bigger than the estimate SIGMA of the
      ! spectral radius, so is more likely to be an upper bound.
      SPCRAD = 1.2D0*SIGMA
      IF (ITER >= 2 .AND. &
          ABS(SIGMA - SIGMAL) <= MAX(SIGMA,SMALL)*0.01D0) THEN

        EV = V - YN
        RETURN
      END IF

      ! The next V(*) is the change in f
      ! scaled so that norm(V - YN) = DYNRM.
      IF (ABS(DFNRM) > 0D0) THEN
        V = YN + (FV - FN)*(DYNRM/DFNRM)
      ELSE
        ! The new V(*) degenerated to YN(*)--"randomly" perturb
        ! current approximation to the eigenvector by changing
        ! the sign of one component.
        INDEX = 1 + MOD(ITER,NEQN)
        V(INDEX) = YN(INDEX) - (V(INDEX) - YN(INDEX))
      END IF

    END DO

    ! Set flag to report a convergence failure.
    RHO_CONVG_FAIL = .TRUE.

  END SUBROUTINE RHO


  SUBROUTINE REPORT_ERRORS(ERR_FLAG)
    INTEGER :: ERR_FLAG

    SELECT CASE (ERR_FLAG)
      CASE (1)
        PRINT *,' Must have 10*UROUND <= RE <= 0.1.'
      CASE (2)
        PRINT *,' AE must be positive.'
      CASE (3)
        PRINT *,' Unable to achieve the desired accuracy with the',&
              ' precision available.'
      CASE (4)
        PRINT *,' Convergence failure in estimating the spectral radius.'
      CASE (5)
        PRINT *,' Singular iteration matrix.'
      CASE (6)
        PRINT *,' Singular filter matrix.'
      CASE (7)
        PRINT *,' A storage allocation error occurred in IRKC.'
    END SELECT
    STOP

  END SUBROUTINE REPORT_ERRORS


  FUNCTION ENORM(V)
    ! Euclidean norm of a vector.
    DOUBLE PRECISION :: V(:),ENORM

    ENORM = SQRT( DOT_PRODUCT(V,V) )

  END FUNCTION ENORM


  FUNCTION RMSNORM(V)
    ! RMS norm of a vector.
    DOUBLE PRECISION :: V(:),RMSNORM

    RMSNORM = SQRT( DOT_PRODUCT(V,V) / SIZE(V) )

  END FUNCTION RMSNORM


  FUNCTION EYE(N) RESULT(IDENTITY)
    INTEGER :: N,I,J
    DOUBLE PRECISION :: IDENTITY(N,N)

    IDENTITY = RESHAPE( (/ ((0D0,J=1,I-1),1D0,(0D0,J=I+1,N), &
                                     I=1,N) /), (/ N,N /) )

  END FUNCTION EYE


  SUBROUTINE FACTOR(A,NEQ,FLAG,PIVOTS)
  ! FACTOR decomposes the matrix A using Gaussian elimination. It
  ! is used in conjunction with SOLVE to solve A*X = B, a system
  ! of NEQ linear equations in NEQ unknowns.  It is sometimes
  ! convenient for the matrix to be defined as the portion
  ! A(1:NEQ,1:NEQ) of a larger matrix, and similarly for the
  ! vector B. FACTOR/SOLVE are coded so that this is possible.
  !
  ! Input arguments:
  !    A       = real matrix A(1:NEQ,1:NEQ) to be factored.
  !    NEQ     = the number of equations to be solved, an integer.
  !
  ! Output arguments:
  !    A       = contains the upper triangular matrix U in its upper
  !              portion (by rows) and a permuted version of a lower
  !              triangular matrix (I - L).  The factorization is
  !              such that (permutation matrix)*A = L*U.
  !   FLAG     = an integer that reports success or the reason for
  !              failure.  FLAG = 0 indicates success. If FLAG > 0,
  !              a zero pivot occurred at equation number FLAG and
  !              the computation was terminated.  FLAG = -1 means
  !              there was an input error: either NEQ <= 0 or A does
  !              not have at least NEQ rows and columns.
  !   PIVOTS   = an integer vector of at least NEQ entries that
  !              records row interchanges. The entry PIVOTS(NEQ) =
  !              (-1)**(number of row interchanges).
  !
  !  When FLAG > 0, the determinant of A is 0 and when FLAG = 0,
  !              det(A) = pivots(NEQ) * A(1,1) *  ... * A(NEQ,NEQ).

    INTEGER :: NEQ, FLAG, PIVOTS(NEQ)
    DOUBLE PRECISION :: A(:,:)

    INTEGER :: I, K, M, OCCURRED(1)
    DOUBLE PRECISION :: T, TEMP(NEQ)

    IF ((NEQ <= 0) .OR. &
        (SIZE(A,DIM=1) < NEQ) .OR. (SIZE(A,DIM=2) < NEQ)) THEN
      FLAG = -1
      RETURN
    END IF
    FLAG = 0
    PIVOTS(NEQ) = 1

    IF (NEQ == 1) THEN
      IF (ABS(A(1,1)) <= 0D0) THEN
        FLAG = 1
        RETURN
      END IF
    END IF

    ! Gaussian elimination with partial pivoting.
    DO K = 1,NEQ-1

      ! Determine the row M containing the largest element in
      ! magnitude to be used as a pivot.
      OCCURRED = MAXLOC(ABS(A(K:NEQ,K)))
      M = OCCURRED(1) + K - 1

      ! If all possible pivots are 0, A is numerically singular.
      IF (ABS(A(M,K)) <= 0D0) THEN
        FLAG = K
        RETURN
      END IF
      PIVOTS(K) = M
      IF (M /= K) THEN
        ! Interchange the current row K with the pivot row M.
        TEMP(K:NEQ) = A(K,K:NEQ)
        A(K,K:NEQ)  = A(M,K:NEQ)
        A(M,K:NEQ)  = TEMP(K:NEQ)
        PIVOTS(NEQ) = - PIVOTS(NEQ)
      END IF

      ! Eliminate subdiagonal entries of column K.
      DO I = K+1,NEQ
        T = A(I,K)/A(K,K)
        A(I,K) = - T
        IF (ABS(T) > 0D0) THEN
          A(I,K+1:NEQ) = A(I,K+1:NEQ) - T*A(K,K+1:NEQ)
        END IF
      END DO

    END DO

    IF (ABS(A(NEQ,NEQ)) <= 0D0) THEN
      FLAG = NEQ
      RETURN
    END IF

  END SUBROUTINE FACTOR


  FUNCTION SOLVE(A,NEQ,PIVOTS,B) RESULT(X)
  ! SOLVE solves A*X = B, a system of NEQ linear equations in NEQ
  ! unknowns using the decomposition obtained from a successful call
  ! to FACTOR.
  !
  ! Input arguments:
  !    A           = output of FACTOR.  Triangular decomposition
  !                  of the coefficient matrix.
  !    NEQ         = the number of equations to be solved.
  !    PIVOTS      = output of FACTOR. Record of row interchanges.
  !    B           = right hand side vector B, a real vector of at
  !                  least NEQ entries.
  !
  ! Output argument:
  !    X           = solution vector of the same size and type as B.

    INTEGER :: NEQ,PIVOTS(:)
    DOUBLE PRECISION :: A(:,:),B(:)
    DOUBLE PRECISION, DIMENSION(SIZE(B)) :: X

    INTEGER :: I, K, M
    DOUBLE PRECISION :: TEMP

    X = B

    IF (NEQ == 1) THEN
      X(1) = X(1)/A(1,1)
    ELSE
      ! Forward elimination.
      DO K = 1,NEQ-1
        M = PIVOTS(K)
        TEMP = X(M)
        X(M) = X(K)
        X(K) = TEMP
        X(K+1:NEQ) = X(K+1:NEQ) + A(K+1:NEQ,K)*X(K)
      END DO

      ! Back substitution.
      X(NEQ) = X(NEQ)/A(NEQ,NEQ)
      DO I = NEQ-1,1,-1
        X(I) = (X(I) - DOT_PRODUCT(A(I,I+1:NEQ),X(I+1:NEQ)))/A(I,I)
      END DO

    END IF

  END FUNCTION SOLVE


END MODULE IRKC_M
