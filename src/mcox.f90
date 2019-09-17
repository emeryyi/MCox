! -------------------------------------------------------------------------------------------------
SUBROUTINE mcox_solution_path(&
    nTasks, nObs, nFeatures, &
    time, iTask, censored, weight, features, &
    penaltyFactor, regPower, regMixing, iExclude, df_max, df_max_ever, &
    nLambda, lambda, lambdaFraction, &
    algorithm, threshold, backtrackingFraction, maxIteration, &
    nLambda_done, beta_path, betaNorm_path, logLik_path, baselineHazard_path, nBeta, nBeta_ever, iEntered, &
    nUpdates, nCycles, iError &
)

! -------------------------------------------------------------------------------------------------
IMPLICIT NONE
! -------------------------------------------------------------------------------------------------
! INPUT

! DIMENSIONS
INTEGER, INTENT(IN)                     :: nTasks
INTEGER, INTENT(IN)                     :: nFeatures
INTEGER, INTENT(IN)                     :: nObs
! DATA INPUT
DOUBLE PRECISION, INTENT(IN)            :: time(nObs)
INTEGER, INTENT(IN)                     :: iTask(nObs)
INTEGER, INTENT(IN)                     :: censored(nObs)
DOUBLE PRECISION, INTENT(INOUT)         :: weight(nObs)
DOUBLE PRECISION, INTENT(INOUT)         :: features(nObs, nFeatures)
! REGULARIZATION INPUT
DOUBLE PRECISION, INTENT(IN)            :: penaltyFactor(nFeatures)
INTEGER, INTENT(IN)                     :: regPower
DOUBLE PRECISION, INTENT(IN)            :: regMixing
INTEGER, INTENT(IN)                     :: iExclude(nFeatures)
INTEGER, INTENT(IN)                     :: df_max
INTEGER, INTENT(IN)                     :: df_max_ever
! REGULARIZATION PATH
INTEGER, INTENT(IN)                     :: nLambda
DOUBLE PRECISION, INTENT(INOUT)         :: lambda(nLambda)
DOUBLE PRECISION, INTENT(IN)            :: lambdaFraction
! ALGORITHM INPUT
INTEGER, INTENT(IN)                     :: algorithm
DOUBLE PRECISION, INTENT(IN)            :: threshold
DOUBLE PRECISION, INTENT(IN)            :: backtrackingFraction
INTEGER, INTENT(IN)                     :: maxIteration

! -------------------------------------------------------------------------------------------------
! OUTPUT

! SOLUTION PATH OUTPUT
INTEGER, INTENT(OUT)                    :: nLambda_done
DOUBLE PRECISION, INTENT(OUT)           :: beta_path(nTasks,nFeatures,nLambda)
DOUBLE PRECISION, INTENT(OUT)           :: betaNorm_path(nTasks,nLambda)
DOUBLE PRECISION, INTENT(OUT)           :: logLik_path(nTasks,nLambda)
DOUBLE PRECISION, ALLOCATABLE, INTENT(OUT):: baselineHazard_path(:,:)
INTEGER, INTENT(OUT)                    :: nBeta(nLambda)
INTEGER, INTENT(OUT)                    :: nBeta_ever(nLambda)
INTEGER, INTENT(OUT)                    :: iEntered(nFeatures)
! STATISTICS
INTEGER, INTENT(OUT)                    :: nUpdates
INTEGER, INTENT(OUT)                    :: nCycles
INTEGER, INTENT(OUT)                    :: iError(5)

! -------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

! DIMENSIONS
INTEGER                                 :: nGroups
! GROUPS (RISK SETS)
INTEGER                                 :: obsBoundByTaskIn(nTasks)
INTEGER                                 :: obsBoundByTaskOut(nTasks)
INTEGER                                 :: groupBoundByTaskIn(nTasks)
INTEGER                                 :: groupBoundByTaskOut(nTasks)
INTEGER                                 :: obsGroupId(nObs)

INTEGER, ALLOCATABLE                    :: obsBoundByGroupIn(:)
INTEGER, ALLOCATABLE                    :: obsBoundByGroupOut(:)
DOUBLE PRECISION, ALLOCATABLE           :: tiesTotalWeight(:)
DOUBLE PRECISION, ALLOCATABLE           :: failureTimes(:)

DOUBLE PRECISION                        :: featureMean(nFeatures,nTasks)
DOUBLE PRECISION                        :: featureStdDev(nFeatures,nTasks)
! TEMPORARY STORAGE
INTEGER, ALLOCATABLE                    :: tmp_int(:)
DOUBLE PRECISION, ALLOCATABLE           :: tmp_dbl(:)

DOUBLE PRECISION                        :: lambda_current, tmp, lambda_decrease, lambda_previous

DOUBLE PRECISION                        :: logLik(nTasks)
DOUBLE PRECISION                        :: gradient(nTasks)
DOUBLE PRECISION                        :: hessian(nTasks)
DOUBLE PRECISION                        :: linearPredictor(nObs)
DOUBLE PRECISION                        :: beta(nTasks, nFeatures)
! ACTIVE SETS
INTEGER                                 :: iStrongRuleSet(nFeatures)

! COUNTERS
INTEGER                                 :: i,j,k,l,g
! FLAGS
INTEGER                                 :: allocationStatus


! -------------------------------------------------------------------------------------------------
! DIMENSIONS AND ALLOCATION
WRITE (*,*) 'allocate_all...'
call allocate_all
WRITE (*,*) 'done.'
! -------------------------------------------------------------------------------------------------
! PREPROCESSING (DIMENSIONS, GROUPS, STANDARDIZATION)
WRITE (*,*) 'preprocessing...'
call preprocessing()
WRITE (*,*) 'done.'
WRITE (*,*) 'nGroups = ', nGroups
! -------------------------------------------------------------------------------------------------
! REALLOCATE (size nObs to size nGroups)
WRITE (*,*) 'reallocate_all...'
call reallocate_all()
WRITE (*,*) 'done.'
! -------------------------------------------------------------------------------------------------
! IF REGULARIZATION PATH IS NOT PROVIDED BY USER
WRITE (*,*) 'initialize_lambda...'
linearPredictor = 0.0D0
gradient = 0.0D0
hessian = 0.0D0
logLik = 0.0D0
beta = 0.0D0
IF(lambdaFraction < 1.0D0) THEN
    call initialize_lambda()
ENDIF
WRITE (*,*) 'done.'
! -------------------------------------------------------------------------------------------------
! LAMBDA LOOP
WRITE (*,*) 'LAMBDA LOOP'
iStrongRuleSet = 1 ! start by excluding all variables
DO l=1,nLambda
    WRITE (*,*) '    ',l
    ! Identify lambda value
    call get_lambda()
    ! First lambda in log-decrease skipped because it must have beta = 0 by construction 
    IF(lambdaFraction >= 1.0D0 .OR. l>1) THEN
        ! Initialize strong set
        call initialize_strong_rule_set()
        ! Strong rule loop
        call strong_rule_loop()
    ENDIF
ENDDO
WRITE (*,*) 'END LAMBDA LOOP.'
! -------------------------------------------------------------------------------------------------






! -------------------------------------------------------------------------------------------------
CONTAINS

    ! ---------------------------------------------------------------------------------------------
    SUBROUTINE preprocessing()
    ! ---------------------------------------------------------------------------------------------
    IMPLICIT NONE
    ! counters
    INTEGER                                 :: newGroup
    DOUBLE PRECISION                        :: timeCurrentGroup
    INTEGER                                 :: currentGroupHasFailure
    ! for ranges 
    INTEGER, ALLOCATABLE                    :: ind(:)
    INTEGER                                 :: length
    ! temporary sums
    DOUBLE PRECISION                        :: weightSum
    DOUBLE PRECISION                        :: featureVariance
    ! ---------------------------------------------------------------------------------------------
    ! Task identification
    ! we suppose that iTask is a vector on integers 1..nTasks in order
    k = 0 !current task
    DO i=1,nObs
        IF (iTask(i)>k) THEN !change of task
            IF (k>0) THEN !end previous task
                obsBoundByTaskOut(k) = i - 1
            ENDIF
            !start new task
            k = k + 1
            obsBoundByTaskIn(k) = i
        ENDIF
    ENDDO
    ! last task was never closed
    obsBoundByTaskOut(nTasks) = nObs

    ! ---------------------------------------------------------------------------------------------
    ! Standardization
    featureMean = 0.0D0
    featureStdDev = 0.0D0 
    featureVariance = 0.0D0

    DO k=1,nTasks
        weightSum = 0.0D0
        length = obsBoundByTaskOut(k) - obsBoundByTaskIn(k) + 1
        ALLOCATE(ind(length), STAT = allocationStatus)
        IF(allocationStatus /= 0) THEN
            iError = (/1, 1, 0, k, 3/)
            RETURN
        ENDIF
        ind = (/ (i, i=obsBoundByTaskIn(k), obsBoundByTaskOut(k), 1 ) /)
        ! Weight standardization
        weightSum = sum(weight(ind))
        IF(weightSum == 0.0D0) THEN
            weight(ind) = 1.0D0 / length
            weightSum = 1.0D0
        ENDIF
        weight(ind) = weight(ind) / weightSum
        ! Feature standardization
        DO j=1,nFeatures
            featureMean(j,k) = sum(features(ind,j) * weight(ind))
            features(ind,j) = features(ind,j) - featureMean(j,k)
            featureVariance = sum(features(ind,j)*features(ind,j)*weight(ind))
            IF(abs(featureVariance)<1.0D-16) THEN
                ! fatal error number 1 at jth feature in task k, before lambda loop
                iError = (/1, 2, j, k, 0/)
                RETURN
            ENDIF
            featureStdDev(j,k) = sqrt(featureVariance)
            features(ind,j) = features(ind,j) / featureStdDev(j,k)
        ENDDO
        IF (allocated(ind)) DEALLOCATE(ind, STAT = allocationStatus)
        IF(allocationStatus /= 0 ) THEN
            iError = (/1, 1, 0, 0, 3/)
            RETURN
        ENDIF
    ENDDO

    ! ---------------------------------------------------------------------------------------------
    ! Groups identification
    k = 0 !current task
    g = 0 !current group
    currentGroupHasFailure = 0
    newGroup = 1 !first observation starts a new group in all cases
    timeCurrentGroup = time(1)
    failureTimes = 0.0D0
    DO i=1,nObs
        IF (iTask(i)==k) THEN !remain in same task
            ! current group already has failure time
            IF (currentGroupHasFailure == 1) THEN
                ! and the current observation is a failure time with lower time
                IF (timeCurrentGroup>time(i)) THEN
                newGroup = 1
                ! it is a failure time so we store the time
                currentGroupHasFailure = 1 - censored(i)
                IF(censored(i) == 0) timeCurrentGroup = time(i)
                ENDIF
                IF (censored(i) == 1) THEN
                    ! do nothing since same set
                ENDIF
            ENDIF
            IF (currentGroupHasFailure == 0) THEN
                ! check if we find one
                currentGroupHasFailure = 1 - censored(i)
                IF(censored(i) == 0) timeCurrentGroup = time(i)
            ENDIF
        ENDIF
        IF (iTask(i)>k) THEN !change of task, restart in any case
            groupBoundByTaskOut(k) = g
            k = k + 1
            newGroup = 1
            groupBoundByTaskIn(k) = g+1
            currentGroupHasFailure = 1 - censored(i)
            IF(censored(i) == 0) timeCurrentGroup = time(i)
        ENDIF


        ! we start a new group
        IF (newGroup == 1) THEN
            obsBoundByGroupOut(g) = i-1
            g = g + 1
            obsBoundByGroupIn(g) = i
            currentGroupHasFailure = 1 - censored(i)
            newGroup = 0
        ENDIF
        !store group id for debug
        obsGroupId(i) = g
        failureTimes(g) = timeCurrentGroup
    ENDDO
    nGroups = g
    ! last group was never closed
    obsBoundByGroupOut(nGroups) = nObs
    groupBoundByTaskOut(nTasks) = nGroups

    ! ---------------------------------------------------------------------------------------------
    ! Groups ties total weight
    tiesTotalWeight = 0.0D0
    DO g=1,nGroups
        length = obsBoundByGroupOut(g) - obsBoundByGroupIn(g) + 1
        ALLOCATE(ind(length), STAT = allocationStatus)
        IF(allocationStatus /= 0) THEN
            iError = (/1, 1, g, 0, 4/)
            RETURN
        ENDIF
        ind = (/ (i, i=obsBoundByGroupIn(g), obsBoundByGroupOut(g),1) /)
        tiesTotalWeight(g) = sum(weight(ind) * (1.0D0-censored(ind)))
        IF (allocated(ind)) DEALLOCATE(ind, STAT = allocationStatus)
        IF(allocationStatus /= 0 ) THEN
            iError = (/1, 1, 0, 0, 4/)
            RETURN
        ENDIF
    ENDDO

    ! ---------------------------------------------------------------------------------------------
    END SUBROUTINE preprocessing
    ! ---------------------------------------------------------------------------------------------



    ! ---------------------------------------------------------------------------------------------
    SUBROUTINE allocate_all
    ! ---------------------------------------------------------------------------------------------    
    ALLOCATE(obsBoundByGroupIn(nObs), STAT = allocationStatus)
    ALLOCATE(obsBoundByGroupOut(nObs), STAT = allocationStatus)
    ALLOCATE(tiesTotalWeight(nObs), STAT = allocationStatus)
    ALLOCATE(failureTimes(nObs), STAT = allocationStatus)
    
    IF(allocationStatus /= 0 ) THEN
        iError = (/1, 1, 0, 0, 1/)
        RETURN
    ENDIF
    ! ---------------------------------------------------------------------------------------------
    END SUBROUTINE allocate_all
    ! ---------------------------------------------------------------------------------------------




    ! ---------------------------------------------------------------------------------------------
    SUBROUTINE reallocate_all()
    ! given nGroups, we can reallocate some of the variables to save memory
    ! ---------------------------------------------------------------------------------------------    
    ALLOCATE(tmp_int(nGroups), STAT = allocationStatus)
    ALLOCATE(tmp_dbl(nGroups), STAT = allocationStatus)

    tmp_int = obsBoundByGroupIn(1:nGroups)
    DEALLOCATE(obsBoundByGroupIn)
    ALLOCATE(obsBoundByGroupIn(nGroups), STAT = allocationStatus)
    obsBoundByGroupIn = tmp_int

    tmp_int = obsBoundByGroupOut(1:nGroups)
    DEALLOCATE(obsBoundByGroupOut)
    ALLOCATE(obsBoundByGroupOut(nGroups), STAT = allocationStatus)
    obsBoundByGroupOut = tmp_int

    tmp_dbl = tiesTotalWeight(1:nGroups)
    DEALLOCATE(tiesTotalWeight)
    ALLOCATE(tiesTotalWeight(nGroups), STAT = allocationStatus)
    tiesTotalWeight = tmp_dbl

    tmp_dbl = failureTimes(1:nGroups)
    DEALLOCATE(failureTimes)
    ALLOCATE(failureTimes(nGroups), STAT = allocationStatus)
    failureTimes = tmp_dbl

    ALLOCATE(baselineHazard_path(nLambda, nGroups), STAT = allocationStatus)

    IF(allocationStatus /= 0 ) THEN
        iError = (/1, 1, 0, 0, 2/)
        RETURN
    ENDIF

    DEALLOCATE(tmp_int)
    DEALLOCATE(tmp_dbl)
    ! ---------------------------------------------------------------------------------------------
    END SUBROUTINE reallocate_all
    ! ---------------------------------------------------------------------------------------------




    ! ---------------------------------------------------------------------------------------------
    SUBROUTINE initialize_lambda()
    ! upon entering that subourtine, linearPredictor = 0 so that we are in the null model
    ! ---------------------------------------------------------------------------------------------
    lambda_current = 0.0D0
    DO j=1,nFeatures
        ! COMPUTE GRADIENT
        call PartialDerivatives(&
            nTasks, nObs, nGroups, obsBoundByGroupIn, obsBoundByGroupOut,&
            groupBoundByTaskIn, groupBoundByTaskOut, &
            weight, features(:,j), censored, tiesTotalWeight, linearPredictor, gradient, hessian&
        )
        ! UPDATE MAXIMUM LAMBDA VALUE
        SELECT CASE (regPower)
            CASE (0)
                tmp = maxval(abs(gradient))/penaltyFactor(j)
            CASE (2)
                tmp = sqrt(sum(gradient**2))/penaltyFactor(j)
        END SELECT
        lambda_current = max(lambda_current, tmp)
        WRITE (*,*) 'tmp            ', tmp
        WRITE (*,*) 'lambda_current ', lambda_current
    ENDDO
    ! PREPARE DECREASE FRACTION
    lambda_decrease = lambdaFraction ** (1.0D0 / (nLambda - 1.0D0))
    WRITE (*,*) 'lambda_decrease ', lambda_decrease
    ! INITIALIZE LOG-LIKELIHOOD
    call LogLikelihood(&
        nTasks, nObs, nGroups, obsBoundByGroupIn, obsBoundByGroupOut,&
        groupBoundByTaskIn, groupBoundByTaskOut, &
        weight, censored, tiesTotalWeight, linearPredictor, logLik&
    )
    ! ---------------------------------------------------------------------------------------------
    END SUBROUTINE initialize_lambda
    ! ---------------------------------------------------------------------------------------------



    
    ! ---------------------------------------------------------------------------------------------
    SUBROUTINE initialize_strong_rule_set()
    ! ---------------------------------------------------------------------------------------------
    WRITE (*,*) '        STRONG RULE SET CONDITION', 2.0D0*lambda_current - lambda_previous
    DO j=1,nFeatures
        IF (iExclude(j) == 1) CYCLE ! skip if excluded so it will always remain excluded
        ! Get gradient
        call PartialDerivatives(&
            nTasks, nObs, nGroups, obsBoundByGroupIn, obsBoundByGroupOut,&
            groupBoundByTaskIn, groupBoundByTaskOut, &
            weight, features(:,j), censored, tiesTotalWeight, linearPredictor, gradient, hessian&
        )
        ! Compute condition
        SELECT CASE (regPower)
            CASE (0)
                tmp = maxval(abs(gradient))/penaltyFactor(j)
            CASE (2)
                tmp = sqrt(sum(gradient**2))/penaltyFactor(j)
        END SELECT
        ! Check condition
        IF(tmp < 2.0D0*lambda_current - lambda_previous) iStrongRuleSet(j) = 0
        WRITE (*,*) '            ', j, tmp, iStrongRuleSet(j)
    ENDDO
    ! ---------------------------------------------------------------------------------------------
    END SUBROUTINE initialize_strong_rule_set
    ! ---------------------------------------------------------------------------------------------
    




    
    ! ---------------------------------------------------------------------------------------------
    SUBROUTINE strong_rule_loop()
    ! ---------------------------------------------------------------------------------------------
    INTEGER                                 :: nStrongRuleCycles, flStrongRule


    WRITE (*,*) '    STRONG RULE LOOP'
    nStrongRuleCycles = 0
    flStrongRule = 0 ! all success flag 
    DO 
        nStrongRuleCycles = nStrongRuleCycles + 1
        WRITE (*,*) '        ', nStrongRuleCycles
        ! Solved penalized MLE with GPG algorithm
        call gpg_descent(&
            nTasks, nFeatures, nObs, nGroups,&
            obsBoundByTaskIn, obsBoundByTaskOut,&
            groupBoundByTaskIn, groupBoundByTaskOut, &
            obsBoundByGroupIn, obsBoundByGroupOut, tiesTotalWeight,&
            weight, features, censored,&
            lambda_current, penaltyFactor, regPower, regMixing, iStrongRuleSet, &
            algorithm, threshold, backtrackingFraction, maxIteration,&
            beta, linearPredictor, logLik,&
            nCycles, nUpdates, iError&
        )
        ! WRITE (*,*) '    ', beta
        ! Catch errors
        IF(iError(1) == 1)THEN
            iError(3) = l
            RETURN
        ENDIF
        ! Check strong rule conditions
        DO j=1,nFeatures
            IF (iStrongRuleSet(j) == 1) CYCLE
            IF (iExclude(j) == 1) CYCLE ! skip if excluded so it will always remain excluded
            ! Get gradient
            call PartialDerivatives(&
                nTasks, nObs, nGroups, obsBoundByGroupIn, obsBoundByGroupOut,&
                groupBoundByTaskIn, groupBoundByTaskOut, &
                weight, features(:,j), censored, tiesTotalWeight, linearPredictor, gradient, hessian&
            )
            ! Compute condition
            SELECT CASE (regPower)
                CASE (0)
                    tmp = maxval(abs(gradient))/penaltyFactor(j)
                CASE (2)
                    tmp = sqrt(sum(gradient**2))/penaltyFactor(j)
            END SELECT
            ! Check condition
            IF(tmp > lambda_current) THEN
                ! we have found a violation, so we reintroduce that feature
                flStrongRule = 1
                iStrongRuleSet(j) = 1
            ENDIF
        ENDDO
        ! No violation found
        IF (flStrongRule == 0) EXIT
        ! Prevent infinite loops
        IF (nStrongRuleCycles >= 1) EXIT
    ENDDO   
    ! ---------------------------------------------------------------------------------------------
    END SUBROUTINE strong_rule_loop
    ! ---------------------------------------------------------------------------------------------
        
    




    
    ! ---------------------------------------------------------------------------------------------
    SUBROUTINE get_lambda()
    ! ---------------------------------------------------------------------------------------------
    IF (lambdaFraction >= 1.0D0) THEN
        ! User defined
        lambda_current = lambda(l)
    ELSE
        ! Logarithmic decrease if not the first one
        IF (l>1) lambda_current = lambda_current * lambda_decrease
    ENDIF
    WRITE (*,*) '        LAMBDA = ',lambda_current
    ! ---------------------------------------------------------------------------------------------
    END SUBROUTINE get_lambda
    ! ---------------------------------------------------------------------------------------------
        
            
    
            
    


! -------------------------------------------------------------------------------------------------
END SUBROUTINE mcox_solution_path
! -------------------------------------------------------------------------------------------------