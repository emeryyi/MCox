! -------------------------------------------------------------------------------------------------
SUBROUTINE preprocessing(&
    time, iTask, censored, weight, features, nObs, nFeatures, &
    nTasks, obsBoundByTaskIn, obsBoundByTaskOut,& 
    groupBoundByTaskIn, groupBoundByTaskOut, &
    nGroups, obsBoundByGroupIn, obsBoundByGroupOut, obsGroupId, failureTimes, &
    tiesTotalWeight, featureMean, featureStdDev, iError &
)
! Deprecated --- now in the CONTAINS clause of the main subroutine
! Preprosess the data set.
! This function standardize features and weight as well as creating the groups for risk set updates
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE

! INPUT
INTEGER, INTENT(IN)                     :: nObs
INTEGER, INTENT(IN)                     :: nFeatures

DOUBLE PRECISION, INTENT(IN)            :: time(nObs)
DOUBLE PRECISION, INTENT(INOUT)         :: weight(nObs)
INTEGER, INTENT(IN)                     :: iTask(nObs)
INTEGER, INTENT(IN)                     :: censored(nObs)
DOUBLE PRECISION, INTENT(INOUT)         :: features(nObs, nFeatures)

! OUTPUT
! observation ranges for each task
INTEGER                                 :: nTasks
INTEGER, INTENT(OUT)                    :: obsBoundByTaskIn(nTasks)
INTEGER, INTENT(OUT)                    :: obsBoundByTaskOut(nTasks)

! group ranges for each task
INTEGER, INTENT(OUT)                    :: groupBoundByTaskIn(nTasks)
INTEGER, INTENT(OUT)                    :: groupBoundByTaskOut(nTasks)

! observation range for each group
INTEGER                                 :: nGroups
INTEGER, INTENT(OUT)                    :: obsBoundByGroupIn(nGroups)
INTEGER, INTENT(OUT)                    :: obsBoundByGroupOut(nGroups)
INTEGER, INTENT(OUT)                    :: obsGroupId(nObs)
DOUBLE PRECISION, INTENT(OUT)           :: failureTimes(nGroups)

! total weight of observed ties in each group
DOUBLE PRECISION, INTENT(OUT)           :: tiesTotalWeight(nGroups)

! features distributions for destandardization
DOUBLE PRECISION, INTENT(OUT)           :: featureMean(nFeatures, nTasks)
DOUBLE PRECISION, INTENT(OUT)           :: featureStdDev(nFeatures, nTasks)

! error flags and info
INTEGER, INTENT(OUT)                    :: iError(5)

! LOCAL VARIABLES

! counters
INTEGER                                 :: i,g,j,k,newGroup
DOUBLE PRECISION                        :: timeCurrentGroup
INTEGER                                 :: currentGroupHasFailure

! for ranges 
INTEGER                                 :: allocationFlag
INTEGER, ALLOCATABLE                    :: ind(:)
INTEGER                                 :: length
! temporary sums
DOUBLE PRECISION                        :: weightSum
DOUBLE PRECISION                        :: featureVariance

iError = 0

! -------------------------------------------------------------------------------------------------
! Task identification
! we suppose that iTask is a vector on integers 1..nTasks in order
nTasks = maxval(iTask)
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

! -------------------------------------------------------------------------------------------------
! Standardization
featureMean = 0.0D0
featureStdDev = 0.0D0 
featureVariance = 0.0D0

DO k=1,nTasks
    weightSum = 0.0D0
    length = obsBoundByTaskOut(k) - obsBoundByTaskIn(k) + 1
    ALLOCATE(ind(length), STAT = allocationFlag)
    IF(allocationFlag /= 0) THEN
        iError = (/1,1,j,k,0/)
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
            iError = (/1,1,j,k,0/)
            RETURN
        ENDIF
        featureStdDev(j,k) = sqrt(featureVariance)
        features(ind,j) = features(ind,j) / featureStdDev(j,k)
    ENDDO
    IF (allocated(ind)) DEALLOCATE(ind, STAT = allocationFlag)
    IF(allocationFlag /= 0) WRITE (*,*) 'Error dealloating ind'
ENDDO

! -------------------------------------------------------------------------------------------------
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

! -------------------------------------------------------------------------------------------------
! Groups ties total weight
tiesTotalWeight = 0.0D0
DO g=1,nGroups
    length = obsBoundByGroupOut(g) - obsBoundByGroupIn(g) + 1
    ALLOCATE(ind(length), STAT = allocationFlag)
    IF(allocationFlag /= 0) THEN
        iError = (/1,1,g,0,0/)
        RETURN
    ENDIF
    ind = (/ (i, i=obsBoundByGroupIn(g), obsBoundByGroupOut(g),1) /)
    tiesTotalWeight(g) = sum(weight(ind) * (1.0D0-censored(ind)))
    IF (allocated(ind)) DEALLOCATE(ind, STAT = allocationFlag)
    IF(allocationFlag /= 0) WRITE (*,*) 'Error dealloating ind'
ENDDO

! -------------------------------------------------------------------------------------------------
END SUBROUTINE preprocessing
! -------------------------------------------------------------------------------------------------
