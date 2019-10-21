
! -------------------------------------------------------------------------------------------------
SUBROUTINE gpg_descent(&
    nTasks, nFeatures, nObs, nGroups,&
    obsBoundByTaskIn, obsBoundByTaskOut,&
    groupBoundByTaskIn, groupBoundByTaskOut, &
    obsBoundByGroupIn, obsBoundByGroupOut, tiesTotalWeight,&
    weight, features, censored,&
    regParameter, penaltyFactor, regPower, regMixing, iStrongRuleSet, &
    algorithm, threshold, backtrackingFraction, maxIteration,&
    beta, linearPredictor, logLik,&
    nCycles, nUpdates, iError&
)
! This function performs proximal gradient descent on an active set
! either with stepsizes minorized by majorizing hessian
! or by backtracking (eventually add IRLS)
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE

! INPUT
INTEGER, INTENT(IN)                     :: nTasks
INTEGER, INTENT(IN)                     :: nFeatures
INTEGER, INTENT(IN)                     :: nObs
INTEGER, INTENT(IN)                     :: nGroups

INTEGER, INTENT(IN)                     :: obsBoundByTaskIn(nTasks)
INTEGER, INTENT(IN)                     :: obsBoundByTaskOut(nTasks)
INTEGER, INTENT(IN)                     :: groupBoundByTaskIn(nTasks)
INTEGER, INTENT(IN)                     :: groupBoundByTaskOut(nTasks)
INTEGER, INTENT(IN)                     :: obsBoundByGroupIn(nGroups)
INTEGER, INTENT(IN)                     :: obsBoundByGroupOut(nGroups)
DOUBLE PRECISION, INTENT(IN)            :: tiesTotalWeight(nGroups)

DOUBLE PRECISION, INTENT(IN)            :: weight(nObs)
DOUBLE PRECISION, INTENT(IN)            :: features(nObs, nFeatures)
INTEGER, INTENT(IN)                     :: censored(nObs)

DOUBLE PRECISION, INTENT(IN)            :: penaltyFactor(nFeatures)
DOUBLE PRECISION, INTENT(IN)            :: regParameter
INTEGER, INTENT(IN)                     :: regPower
DOUBLE PRECISION, INTENT(IN)            :: regMixing
INTEGER, INTENT(INOUT)                  :: iStrongRuleSet(nFeatures)

INTEGER, INTENT(IN)                     :: algorithm
DOUBLE PRECISION, INTENT(IN)            :: threshold
DOUBLE PRECISION, INTENT(IN)            :: backtrackingFraction
INTEGER, INTENT(IN)                     :: maxIteration

! INPUT/OUTPUT
DOUBLE PRECISION, INTENT(INOUT)         :: linearPredictor(nObs)
DOUBLE PRECISION, INTENT(INOUT)         :: beta(nFeatures, nTasks)
DOUBLE PRECISION, INTENT(INOUT)         :: logLik(nTasks)

! OUTPUT
INTEGER, INTENT(OUT)                    :: nCycles
INTEGER, INTENT(OUT)                    :: nUpdates
INTEGER, INTENT(OUT)                    :: iError(5)

! - LOCAL
DOUBLE PRECISION, PARAMETER             :: small = 1.0D-16
INTEGER                                 :: j
DOUBLE PRECISION                        :: stepsize(nFeatures)
DOUBLE PRECISION                        :: beta_tmp(nFeatures,nTasks)
INTEGER                                 :: converged
DOUBLE PRECISION                        :: gradient(nFeatures,nTasks)
DOUBLE PRECISION                        :: hessian(nFeatures,nTasks)
DOUBLE PRECISION                        :: logLik_tmp(nTasks)
INTEGER                                 :: nCycles_local
INTEGER                                 :: iActive(nFeatures)


! - INITIALIZATIONS
converged = 0
nCycles_local = 0
beta_tmp = 0.0D0
gradient = 0.0D0
hessian = 0.0D0
logLik_tmp = 0.0D0
iActive = iStrongRuleSet
DO WHILE (converged == 0)
    stepsize = 0.0D0
    nCycles_local = nCycles_local + 1
    nCycles = nCycles + 1
    IF(nCycles > maxIteration) THEN
        iError = (/ 2,3,0,0,0 /)
        EXIT
    ENDIF
    ! Save previous estimates in current variable to update
    beta_tmp = beta
    logLik_tmp = logLik
    SELECT CASE (algorithm)
        CASE (2)
            CALL gpg_cycle_backtracking(&
            nTasks, nFeatures, nObs, nGroups,&
            obsBoundByTaskIn, obsBoundByTaskOut,&
            groupBoundByTaskIn, groupBoundByTaskOut, &
            obsBoundByGroupIn, obsBoundByGroupOut, tiesTotalWeight,&
            weight, features, censored, iActive,&
            regParameter, penaltyFactor, regPower, regMixing, backtrackingFraction, &
            beta_tmp, linearPredictor, logLik_tmp, stepsize,&
            nUpdates, gradient, hessian, iError&
            )
        CASE DEFAULT
            CALL gpg_cycle(&
                nTasks, nFeatures, nObs, nGroups,&
                obsBoundByTaskIn, obsBoundByTaskOut,&
                groupBoundByTaskIn, groupBoundByTaskOut, &
                obsBoundByGroupIn, obsBoundByGroupOut, tiesTotalWeight,&
                weight, features, censored, iActive,&
                regParameter, penaltyFactor, regPower, regMixing,&
                beta_tmp, linearPredictor, logLik_tmp, stepsize,&
                nUpdates, gradient, hessian, iError&
            )
    END SELECT
    ! Check for fatal errors
    IF(iError(1) == 1) RETURN
    ! Now beta_tmp contains the new estimates and llk_tmp the new likelihood
    ! Compute difference in estimates
    IF(maxval(abs(beta_tmp-beta)) < threshold)THEN
        converged = 1
    ENDIF
    beta = beta_tmp
    logLik = logLik_tmp
    ! Update the active set
    IF(nCycles_local == 1)THEN
        DO j=1,nFeatures
            IF(maxval(abs(beta(j,:))) < small) iActive(j) = 1
        ENDDO
    ENDIF
    ! Non convergence
    IF(nCycles_local >= maxIteration / 10) THEN
        iError = (/ 2,4,0,0,0 /)
        EXIT
    ENDIF
ENDDO

! -------------------------------------------------------------------------------------------------
END SUBROUTINE gpg_descent
! -------------------------------------------------------------------------------------------------




! -------------------------------------------------------------------------------------------------
SUBROUTINE gpg_cycle(&
    nTasks, nFeatures, nObs, nGroups,&
    obsBoundByTaskIn, obsBoundByTaskOut,&
    groupBoundByTaskIn, groupBoundByTaskOut, &
    obsBoundByGroupIn, obsBoundByGroupOut, tiesTotalWeight,&
    weight, features, censored, iActive, &
    regParameter, penaltyFactor, regPower, regMixing,&
    beta, linearPredictor, logLik, stepsize,&
    nUpdates, gradient, hessian, iError&
)
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE
! This function performs one cycle of proximal gradient descent on an active set
! with stepsizes minorized by majorizing hessian

! INPUT
INTEGER, INTENT(IN)                     :: nTasks
INTEGER, INTENT(IN)                     :: nFeatures
INTEGER, INTENT(IN)                     :: nObs
INTEGER, INTENT(IN)                     :: nGroups

INTEGER, INTENT(IN)                     :: obsBoundByTaskIn(nTasks)
INTEGER, INTENT(IN)                     :: obsBoundByTaskOut(nTasks)
INTEGER, INTENT(IN)                     :: groupBoundByTaskIn(nTasks)
INTEGER, INTENT(IN)                     :: groupBoundByTaskOut(nTasks)
INTEGER, INTENT(IN)                     :: obsBoundByGroupIn(nGroups)
INTEGER, INTENT(IN)                     :: obsBoundByGroupOut(nGroups)
DOUBLE PRECISION, INTENT(IN)            :: tiesTotalWeight(nGroups)

DOUBLE PRECISION, INTENT(IN)            :: weight(nObs)
DOUBLE PRECISION, INTENT(IN)            :: features(nObs, nFeatures)
INTEGER, INTENT(IN)                     :: censored(nObs)
INTEGER, INTENT (IN)                    :: iActive(nFeatures)

DOUBLE PRECISION, INTENT(IN)            :: regParameter
DOUBLE PRECISION, INTENT(IN)            :: penaltyFactor(nFeatures)
INTEGER, INTENT(IN)                     :: regPower
DOUBLE PRECISION, INTENT(IN)            :: regMixing

! INPUT/OUTPUT
DOUBLE PRECISION, INTENT(INOUT)         :: linearPredictor(nObs)
DOUBLE PRECISION, INTENT(INOUT)         :: beta(nFeatures, nTasks)
DOUBLE PRECISION, INTENT(INOUT)         :: logLik(nTasks)
DOUBLE PRECISION, INTENT(INOUT)         :: stepsize(nFeatures)
INTEGER, INTENT(INOUT)                  :: nUpdates

! OUTPUT
DOUBLE PRECISION, INTENT(OUT)           :: gradient(nFeatures, nTasks)
DOUBLE PRECISION, INTENT(OUT)           :: hessian(nFeatures, nTasks)
INTEGER, INTENT(OUT)                    :: iError(5)

! - LOCAL
DOUBLE PRECISION, PARAMETER             :: small = 1.0D-16
INTEGER                                 :: i, j, k
DOUBLE PRECISION                        :: beta_previous(nFeatures,nTasks)
DOUBLE PRECISION                        :: gradient_tmp(nTasks)
DOUBLE PRECISION                        :: hessian_tmp(nTasks)
DOUBLE PRECISION                        :: beta_tmp(nTasks)
! for ranges 
INTEGER, ALLOCATABLE                    :: ind(:)
INTEGER                                 :: length

! INITIALIZATION
gradient = 0.0D0
hessian = 0.0D0
beta_previous = beta
! GPG CYCLE
DO j=1,nFeatures
    IF(iActive(j)==1) CYCLE
    gradient_tmp = 0.0D0
    hessian_tmp = 0.0D0
    ! GRADIENT AND HESSIAN
    CALL PartialDerivatives(&
        nTasks, nObs, nGroups, obsBoundByGroupIn, obsBoundByGroupOut,&
        groupBoundByTaskIn, groupBoundByTaskOut, &
        weight, features(:,j), censored, tiesTotalWeight, linearPredictor, gradient_tmp, hessian_tmp &
    )
    gradient(j,:) = gradient_tmp
    hessian(j,:) = hessian_tmp
    ! MAJORIZATION IF STEPSIZE NOT SUPPLIED
    IF(stepsize(j) == 0.0D0) THEN
        stepsize(j) = maxval(hessian_tmp)
        IF(stepsize(j) < small) THEN
            iError = (/2,5,j,0,0/)
            RETURN
        ENDIF
        stepsize(j) = 1.0D0 / stepsize(j)
    ENDIF
    ! PROXIMAL GRADIENT DESCENT STEP
    beta_tmp = beta(j,:)
    CALL Proximal(regParameter, penaltyFactor(j), stepsize(j), regMixing, regPower, nTasks, beta_tmp, gradient_tmp)
    beta(j,:) = beta_tmp
    ! UPDATE LINEAR PREDICTOR
    DO k=1,nTasks
        length = obsBoundByTaskOut(k) - obsBoundByTaskIn(k) + 1
        ALLOCATE(ind(length))
        ind = (/ (i,i=obsBoundByTaskIn(k),obsBoundByTaskOut(k),1) /)
        linearPredictor(ind) = linearPredictor(ind) + features(ind,j) * (beta(j,k) - beta_previous(j,k))
        DEALLOCATE(ind)
    ENDDO
    nUpdates = nUpdates + 1
    ! UPDATE LIKELIHOOD
    CALL LogLikelihood(&
        nTasks, nObs, nGroups, obsBoundByGroupIn, obsBoundByGroupOut,&
        groupBoundByTaskIn, groupBoundByTaskOut, &
        weight, censored, tiesTotalWeight, linearPredictor, logLik&
    )
ENDDO

! -------------------------------------------------------------------------------------------------
END SUBROUTINE gpg_cycle
! -------------------------------------------------------------------------------------------------




! -------------------------------------------------------------------------------------------------
SUBROUTINE gpg_cycle_backtracking(&
    nTasks, nFeatures, nObs, nGroups,&
    obsBoundByTaskIn, obsBoundByTaskOut,&
    groupBoundByTaskIn, groupBoundByTaskOut, &
    obsBoundByGroupIn, obsBoundByGroupOut, tiesTotalWeight,&
    weight, features, censored, iActive, &
    regParameter, penaltyFactor, regPower, regMixing, backtrackingFraction, &
    beta, linearPredictor, logLik, stepsize,&
    nUpdates, gradient, hessian, iError&
)
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE
! This function performs one cycle of proximal gradient descent on an active set
! with stepsizes minorized by majorizing hessian

! INPUT
INTEGER, INTENT(IN)                     :: nTasks
INTEGER, INTENT(IN)                     :: nFeatures
INTEGER, INTENT(IN)                     :: nObs
INTEGER, INTENT(IN)                     :: nGroups

INTEGER, INTENT(IN)                     :: obsBoundByTaskIn(nTasks)
INTEGER, INTENT(IN)                     :: obsBoundByTaskOut(nTasks)
INTEGER, INTENT(IN)                     :: groupBoundByTaskIn(nTasks)
INTEGER, INTENT(IN)                     :: groupBoundByTaskOut(nTasks)
INTEGER, INTENT(IN)                     :: obsBoundByGroupIn(nGroups)
INTEGER, INTENT(IN)                     :: obsBoundByGroupOut(nGroups)
DOUBLE PRECISION, INTENT(IN)            :: tiesTotalWeight(nGroups)

DOUBLE PRECISION, INTENT(IN)            :: weight(nObs)
DOUBLE PRECISION, INTENT(IN)            :: features(nObs, nFeatures)
INTEGER, INTENT(IN)                     :: censored(nObs)
INTEGER, INTENT (IN)                    :: iActive(nFeatures)

DOUBLE PRECISION, INTENT(IN)            :: regParameter
DOUBLE PRECISION, INTENT(IN)            :: penaltyFactor(nFeatures)
INTEGER, INTENT(IN)                     :: regPower
DOUBLE PRECISION, INTENT(IN)            :: regMixing
DOUBLE PRECISION, INTENT(IN)            :: backtrackingFraction

! INPUT/OUTPUT
DOUBLE PRECISION, INTENT(INOUT)         :: linearPredictor(nObs)
DOUBLE PRECISION, INTENT(INOUT)         :: beta(nFeatures, nTasks)
DOUBLE PRECISION, INTENT(INOUT)         :: logLik(nTasks)
DOUBLE PRECISION, INTENT(INOUT)         :: stepsize(nFeatures)
INTEGER, INTENT(INOUT)                  :: nUpdates

! OUTPUT
DOUBLE PRECISION, INTENT(OUT)           :: gradient(nFeatures, nTasks)
DOUBLE PRECISION, INTENT(OUT)           :: hessian(nFeatures, nTasks)
INTEGER, INTENT(OUT)                    :: iError(5)

! - LOCAL
DOUBLE PRECISION, PARAMETER             :: small = 1.0D-16
INTEGER                                 :: i, j, k
INTEGER                                 :: armijoFlag
DOUBLE PRECISION                        :: beta_previous(nFeatures,nTasks)
DOUBLE PRECISION                        :: condition, nTries
DOUBLE PRECISION                        :: linearPredictor_tmp(nObs)
DOUBLE PRECISION                        :: logLik_tmp(nTasks)
DOUBLE PRECISION                        :: gradient_tmp(nTasks)
DOUBLE PRECISION                        :: gradient_pseudo(nTasks)
DOUBLE PRECISION                        :: hessian_tmp(nTasks)
DOUBLE PRECISION                        :: beta_tmp(nTasks)
! for ranges 
INTEGER, ALLOCATABLE                    :: ind(:)
INTEGER                                 :: length

! INITIALIZATION
gradient = 0.0D0
hessian = 0.0D0
beta_previous = beta
CALL LogLikelihood(&
    nTasks, nObs, nGroups, obsBoundByGroupIn, obsBoundByGroupOut,&
    groupBoundByTaskIn, groupBoundByTaskOut, &
    weight, censored, tiesTotalWeight, linearPredictor, logLik&
)
! GPG CYCLE
DO j=1,nFeatures
    IF(iActive(j)==1) CYCLE
    gradient_tmp = 0.0D0
    gradient_pseudo = 0.0D0
    hessian_tmp = 0.0D0
    condition = 0.0D0
    ! GRADIENT AND HESSIAN
    CALL PartialDerivatives(&
        nTasks, nObs, nGroups, obsBoundByGroupIn, obsBoundByGroupOut,&
        groupBoundByTaskIn, groupBoundByTaskOut, &
        weight, features(:,j), censored, tiesTotalWeight, linearPredictor, gradient_tmp, hessian_tmp &
    )
    gradient(j,:) = gradient_tmp
    hessian(j,:) = hessian_tmp
    ! COMPUTE INITIAL STEPSIZE
    IF(stepsize(j) == 0.0D0) THEN
        stepsize(j) = minval(hessian_tmp)
        IF(stepsize(j) < small) THEN
            iError = (/2,5,j,0,0/)
            RETURN
        ENDIF
        stepsize(j) = 1.0D0 / stepsize(j)
    ENDIF
    ! BACKTRACKING LOOP
    armijoFlag = 0
    nTries = 0 
    DO WHILE (armijoFlag == 0 .AND. nTries <= 20)
        linearPredictor_tmp = linearPredictor
        nUpdates = nUpdates + 1
        nTries = nTries + 1
        ! TEMPORARY UPDATE
        beta_tmp = beta(j,:)
        CALL Proximal(regParameter, penaltyFactor(j), stepsize(j), regMixing, regPower, nTasks, beta_tmp, gradient_tmp)
        ! TEMPORARY LINEAR PREDICTOR
        DO k=1,nTasks
            length = obsBoundByTaskOut(k) - obsBoundByTaskIn(k) + 1
            ALLOCATE(ind(length))
            ind = (/ (i,i=obsBoundByTaskIn(k),obsBoundByTaskOut(k),1) /)
            linearPredictor_tmp(ind) = linearPredictor_tmp(ind) + features(ind,j) * (beta_tmp(k) - beta_previous(j,k))
            DEALLOCATE(ind)
        ENDDO
        ! PSEUDO GRADIENT
        gradient_pseudo = (beta_tmp - beta_previous(j,:)) / stepsize(j)
        ! TEMPORARY LOG LIKELIHOOD
        CALL LogLikelihood(&
            nTasks, nObs, nGroups, obsBoundByGroupIn, obsBoundByGroupOut,&
            groupBoundByTaskIn, groupBoundByTaskOut, &
            weight, censored, tiesTotalWeight, linearPredictor_tmp, logLik_tmp&
        )
        ! COMPUTE CONDITION
        condition = sum(logLik) - stepsize(j) * sum(gradient_tmp * gradient_pseudo) +&
            0.5D0 * stepsize(j) * sum(gradient_pseudo * gradient_pseudo)
        ! CHECK ARMIJO CONDITION
        IF (condition + 1e-3 < sum(logLik_tmp)) THEN
            ! FAIL: REDUCE STEPSIZE
            stepsize(j) = stepsize(j) * backtrackingFraction
        ELSE
            ! PASS: EXIT LOOP
            armijoFlag = 1
        ENDIF
    ENDDO
    beta(j,:) = beta_tmp
    linearPredictor = linearPredictor_tmp
    logLik = logLik_tmp
ENDDO

! -------------------------------------------------------------------------------------------------
END SUBROUTINE gpg_cycle_backtracking
! -------------------------------------------------------------------------------------------------






