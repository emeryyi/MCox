! -------------------------------------------------------------------------------------------------
SUBROUTINE LogLikelihood(&
    nTasks, nObs, nGroups, obsBoundByGroupIn, obsBoundByGroupOut,&
    groupBoundByTaskIn, groupBoundByTaskOut, &
    weight, censored, tiesTotalWeight, linearPredictor, logLik&
)
! Computes the log-likelihood in each task
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE

! INPUT
INTEGER, INTENT(IN)                     :: nTasks
INTEGER, INTENT(IN)                     :: nObs
INTEGER, INTENT(IN)                     :: nGroups
INTEGER, INTENT(IN)                     :: obsBoundByGroupIn(nGroups)
INTEGER, INTENT(IN)                     :: obsBoundByGroupOut(nGroups)
INTEGER, INTENT(IN)                     :: groupBoundByTaskIn(nTasks)
INTEGER, INTENT(IN)                     :: groupBoundByTaskOut(nTasks)
DOUBLE PRECISION, INTENT(IN)            :: tiesTotalWeight(nGroups)
DOUBLE PRECISION, INTENT(IN)            :: weight(nObs)
INTEGER, INTENT(IN)                     :: censored(nObs)
DOUBLE PRECISION, INTENT(IN)            :: linearPredictor(nObs)

! OUTPUT
DOUBLE PRECISION, INTENT(OUT)           :: logLik(nTasks)

! LOCAL VARIABLES
! counters
INTEGER                                 :: i,g,k
! partial results
DOUBLE PRECISION                        :: partialSumWeightedExp
DOUBLE PRECISION                        :: partialSumWeightedFeaturesTies
DOUBLE PRECISION                        :: exponentialPredictor(nObs)
! for ranges 
INTEGER, ALLOCATABLE                    :: ind(:)
INTEGER                                 :: length

!--------------------------------------------------------------------------------------------------
! - INITIALIZATION
!--------------------------------------------------------------------------------------------------
logLik = 0.0D0
exponentialPredictor = exp(linearPredictor)
!--------------------------------------------------------------------------------------------------
! - COMPUTATION
!--------------------------------------------------------------------------------------------------
DO k=1,nTasks
    partialSumWeightedExp = 0.0D0
    DO g=groupBoundByTaskIn(k),groupBoundByTaskOut(k)
        length = obsBoundByGroupOut(g) - obsBoundByGroupIn(g) + 1
        ALLOCATE(ind(length))
        ind = (/ (i,i=obsBoundByGroupIn(g),obsBoundByGroupOut(g),1) /)
        partialSumWeightedExp = partialSumWeightedExp + &
            sum( weight(ind) * exponentialPredictor(ind) )
        partialSumWeightedFeaturesTies = sum(linearPredictor(ind) * weight(ind) * (1.0D0-censored(ind)) )
        DEALLOCATE(ind)
        logLik(k) = logLik(k) + partialSumWeightedFeaturesTies - &
            tiesTotalWeight(g) * log(partialSumWeightedExp)
    ENDDO
ENDDO
! -------------------------------------------------------------------------------------------------
END SUBROUTINE LogLikelihood
! -------------------------------------------------------------------------------------------------







! -------------------------------------------------------------------------------------------------
SUBROUTINE PartialDerivatives(&
    nTasks, nObs, nGroups, obsBoundByGroupIn, obsBoundByGroupOut,&
    groupBoundByTaskIn, groupBoundByTaskOut, &
    weight, features, censored, tiesTotalWeight, linearPredictor, gradient, hessian&
)
! Computes the gradient and Hessian for a feature
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE

! INPUT
INTEGER, INTENT(IN)                     :: nTasks
INTEGER, INTENT(IN)                     :: nObs
INTEGER, INTENT(IN)                     :: nGroups
INTEGER, INTENT(IN)                     :: obsBoundByGroupIn(nGroups)
INTEGER, INTENT(IN)                     :: obsBoundByGroupOut(nGroups)
INTEGER, INTENT(IN)                     :: groupBoundByTaskIn(nTasks)
INTEGER, INTENT(IN)                     :: groupBoundByTaskOut(nTasks)
DOUBLE PRECISION, INTENT(IN)            :: tiesTotalWeight(nGroups)
DOUBLE PRECISION, INTENT(IN)            :: weight(nObs)
DOUBLE PRECISION, INTENT(IN)            :: features(nObs)
INTEGER, INTENT(IN)                     :: censored(nObs)
DOUBLE PRECISION, INTENT(IN)            :: linearPredictor(nObs)

! OUTPUT
DOUBLE PRECISION, INTENT(OUT)           :: gradient(nTasks)
DOUBLE PRECISION, INTENT(OUT)           :: hessian(nTasks)

! LOCAL VARIABLES
! counters
INTEGER                                 :: i,g,k
! partial results
DOUBLE PRECISION                        :: partialSumWeightedExp
DOUBLE PRECISION                        :: partialSumWeightedExpFeatures
DOUBLE PRECISION                        :: partialSumWeightedExpFeaturesSquared
DOUBLE PRECISION                        :: partialSumWeightedFeaturesTies
DOUBLE PRECISION                        :: partialSumWeightedExpInverse
DOUBLE PRECISION                        :: exponentialPredictor(nObs)
! for ranges 
INTEGER, ALLOCATABLE                    :: ind(:)
INTEGER                                 :: length



!--------------------------------------------------------------------------------------------------
! - INITIALIZATION
!--------------------------------------------------------------------------------------------------
gradient = 0.0D0
hessian = 0.0D0
exponentialPredictor = exp(linearPredictor)
!--------------------------------------------------------------------------------------------------
! - ALGORITHM
!--------------------------------------------------------------------------------------------------
DO k=1,nTasks
    partialSumWeightedExp = 0.0D0
    partialSumWeightedExpFeatures = 0.0D0
    partialSumWeightedExpFeaturesSquared = 0.0D0
    partialSumWeightedFeaturesTies = 0.0D0
    partialSumWeightedExpInverse = 0.0D0
    DO g=groupBoundByTaskIn(k),groupBoundByTaskOut(k)
        length = obsBoundByGroupOut(g) - obsBoundByGroupIn(g) + 1
        ALLOCATE(ind(length))
        ind = (/ (i,i=obsBoundByGroupIn(g),obsBoundByGroupOut(g),1) /)
        partialSumWeightedExp = partialSumWeightedExp + sum( weight(ind)* exponentialPredictor(ind) )
        partialSumWeightedExpInverse = 1.0D0/partialSumWeightedExp
        partialSumWeightedExpFeatures = partialSumWeightedExpFeatures + &
            sum(weight(ind) * exponentialPredictor(ind) * features(ind) )
        partialSumWeightedExpFeaturesSquared = partialSumWeightedExpFeaturesSquared + &
            sum(weight(ind) * exponentialPredictor(ind) * features(ind) ** 2.0D0 )
        partialSumWeightedFeaturesTies = sum(features(ind) * weight(ind) * (1.0D0-censored(ind)) )
        gradient(k) = gradient(k) - partialSumWeightedFeaturesTies + &
            tiesTotalWeight(g) * partialSumWeightedExpInverse * partialSumWeightedExpFeatures
        hessian(k) = hessian(k) - &
            (tiesTotalWeight(g) * partialSumWeightedExpInverse * partialSumWeightedExpFeatures) ** 2.0D0 + &
            tiesTotalWeight(g) * partialSumWeightedExpInverse * partialSumWeightedExpFeaturesSquared
        DEALLOCATE(ind)
    ENDDO
ENDDO

! -------------------------------------------------------------------------------------------------
END SUBROUTINE PartialDerivatives
! -------------------------------------------------------------------------------------------------
    







! -------------------------------------------------------------------------------------------------
SUBROUTINE BaselineHazard(&
    nTasks, nObs, nGroups, obsBoundByGroupIn, obsBoundByGroupOut,&
    groupBoundByTaskIn, groupBoundByTaskOut, &
    weight, linearPredictor, baselineHazard_estimate&
)
! Computes the non-parametric estimates of the baseline hazard
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE

! INPUT
INTEGER, INTENT(IN)                     :: nTasks
INTEGER, INTENT(IN)                     :: nObs
INTEGER, INTENT(IN)                     :: nGroups
INTEGER, INTENT(IN)                     :: obsBoundByGroupIn(nGroups)
INTEGER, INTENT(IN)                     :: obsBoundByGroupOut(nGroups)
INTEGER, INTENT(IN)                     :: groupBoundByTaskIn(nTasks)
INTEGER, INTENT(IN)                     :: groupBoundByTaskOut(nTasks)
DOUBLE PRECISION, INTENT(IN)            :: weight(nObs)
DOUBLE PRECISION, INTENT(IN)            :: linearPredictor(nObs)

! OUTPUT
DOUBLE PRECISION, INTENT(OUT)           :: baselineHazard_estimate(nObs)

! LOCAL VARIABLES
! counters
INTEGER                                 :: i,g,k
! partial results
DOUBLE PRECISION                        :: partialSumWeightedExp
DOUBLE PRECISION                        :: exponentialPredictor(nObs)
! for ranges 
INTEGER, ALLOCATABLE                    :: ind(:)
INTEGER                                 :: length



!--------------------------------------------------------------------------------------------------
! - INITIALIZATION
!--------------------------------------------------------------------------------------------------
exponentialPredictor = exp(linearPredictor)
baselineHazard_estimate = 0.0D0
!--------------------------------------------------------------------------------------------------
! - ALGORITHM
!--------------------------------------------------------------------------------------------------
DO k=1,nTasks
    partialSumWeightedExp = 0.0D0
    DO g=groupBoundByTaskIn(k),groupBoundByTaskOut(k)
        length = obsBoundByGroupOut(g) - obsBoundByGroupIn(g) + 1
        ALLOCATE(ind(length))
        ind = (/ (i,i=obsBoundByGroupIn(g),obsBoundByGroupOut(g),1) /)
        partialSumWeightedExp = partialSumWeightedExp + sum( weight(ind)* exponentialPredictor(ind) )
        IF(partialSumWeightedExp > 1.0D-16) baselineHazard_estimate(g) = 1.0D0/partialSumWeightedExp
        DEALLOCATE(ind)
    ENDDO
ENDDO

! -------------------------------------------------------------------------------------------------
END SUBROUTINE BaselineHazard
! -------------------------------------------------------------------------------------------------