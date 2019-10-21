! -------------------------------------------------------------------------------------------------
SUBROUTINE Proximal(lambda, penaltyFactor, stepsize, regMixing, regPower, nTasks, beta, gradient)
! Performs the proximal operation in the following case
!   prox_{lambda * penaltyFactor * stepsize * P^{regMixing, regPower}} (beta - stepsize * gradient)
!   where P^{regMixing, regPower} = regMixing||.||_regPower + (1-regMixing)||.||_1 is the sparre group lasso penalty
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT
DOUBLE PRECISION, INTENT(IN)        :: lambda                   ! lambda
DOUBLE PRECISION, INTENT(IN)        :: penaltyFactor            ! penalty factor
DOUBLE PRECISION, INTENT(IN)        :: stepsize                 ! stepsize (related to inverse hessian)
DOUBLE PRECISION, INTENT(IN)        :: regMixing                ! mixing of regularization
INTEGER, INTENT(IN)                 :: regPower                 ! 0=infty of 2=2
INTEGER, INTENT(IN)                 :: nTasks                   ! dimension
DOUBLE PRECISION, INTENT(INOUT)     :: beta(nTasks)             ! previous estimate to update, contains update at exit
DOUBLE PRECISION, INTENT(IN)        :: gradient(nTasks)         ! gradientient
! LOCAL
DOUBLE PRECISION        :: s(nTasks)                            ! for soft-thresholding
DOUBLE PRECISION        :: sgn(nTasks)                          ! sign of s
INTEGER                 :: k                                    ! counter
DOUBLE PRECISION        :: proj(nTasks)                         ! for l1 projection
DOUBLE PRECISION        :: tmp, norm                            ! for l2 projection
! INITIALIZATION
s = beta - stepsize*gradient
sgn = 0.0D0
proj = 0.0D0
tmp = 0.0D0
norm = 0.0D0
! RETRIEVE SIGN
DO k=1,nTasks
    IF(s(k)>0.0D0) sgn(k) = 1.0D0
    IF(s(k)<0.0D0) sgn(k) = -1.0D0
ENDDO
! SOFT THRESHOLDING (l1 part)
s = abs(s) - regMixing*lambda*stepsize*penaltyFactor
DO k=1,nTasks
    IF(s(k)<0.0D0) s(k) = 0.0D0
ENDDO
s=s*sgn
! GROUP REGULARIZATION (lq part)
SELECT CASE (regPower)
    ! q = infty
    CASE (0)
        IF (sum(abs(s)) < lambda * penaltyFactor * stepsize) THEN
            beta = 0.0D0
        ELSE
            call ProjB1Mich(s, nTasks, lambda*penaltyFactor*(1.0D0 - regMixing)*stepsize, proj)
            beta = s - proj
        ENDIF
    ! q = 2
    CASE (2)
        norm = sqrt(dot_product(s,s))
        tmp = norm - penaltyFactor * lambda * (1.0D0 - regMixing) * stepsize
        IF (tmp > 0.0D0) THEN
            beta = s * tmp / norm
        ELSE
            beta = 0.0D0
        ENDIF
END SELECT
! -------------------------------------------------------------------------------------------------
END SUBROUTINE Proximal
! -------------------------------------------------------------------------------------------------



! -------------------------------------------------------------------------------------------------
SUBROUTINE ProjB1Mich(vector, dim, radius, projection)
! Performs projection onto L1 ball using Michelot's algorithm using active set
! See Condat, L. (2016) Fast projection onto the simplex and the l1 ball.
! -------------------------------------------------------------------------------------------------
    
! - VARIABLES DECLARATIONS -------------------------
    ! - INPUTS -
    INTEGER, INTENT(IN)             :: dim
    DOUBLE PRECISION, INTENT(IN)    :: vector(dim)
    DOUBLE PRECISION, INTENT(IN)    :: radius
    ! - OUTPUTS -
    DOUBLE PRECISION, INTENT(OUT)   :: projection(dim)
    ! - LOCAL VARS -
    DOUBLE PRECISION                :: vector_positive(dim)
    DOUBLE PRECISION                :: tau, rho
    INTEGER                         :: j,changed
    INTEGER                         :: sgn(dim)
    INTEGER                         :: active(dim)
! - END VARIABLE DECLARATIONS ----------------------

    ! - absolute components -
    vector_positive = abs(vector)
    ! - if already within the ball
    IF(sum(vector_positive) <= radius) THEN
        projection = vector
    ! - otherwise we have to project -
    ELSE
        ! - retrieve sgn -
        DO j=1,dim
            IF(vector(j) == 0.0D0) THEN
                sgn(j) = 0
            ELSEIF(vector(j) > 0.0D0)THEN
                sgn(j) = 1
            ELSE
                sgn(j) = -1
            ENDIF
        ENDDO
        ! - cycle -
        active = 1
        rho = (sum(vector_positive) - radius) / sum(active)
        DO
            changed = 0
            DO j=1,dim
                IF ( active(j) == 0 ) CYCLE
                IF ( vector_positive(j) <= rho ) THEN
                    active(j) = 0
                    changed = 1
                ENDIF
            ENDDO
            rho = ( sum(vector_positive * active) - radius ) / sum(active)
            IF ( changed == 0 ) EXIT
        ENDDO
        tau = rho
        ! - produce projection on simplex and adjust sign -
        DO j=1,dim
            projection(j) = max(vector_positive(j) - tau, 0.0D0)*sgn(j)
        ENDDO
    ENDIF
! -------------------------------------------------------------------------------------------------
END SUBROUTINE ProjB1Mich
! -------------------------------------------------------------------------------------------------

        