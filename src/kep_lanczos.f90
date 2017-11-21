    SUBROUTINE kep_lanczos(n,imax,its,ijob,v,bv,alpha,beta,theta)
    IMPLICIT NONE
    INTEGER :: n,imax,its,ijob
    DOUBLE PRECISION :: v,bv,alpha,beta,theta
!
    DIMENSION :: alpha(imax),theta(imax)
    DIMENSION :: beta(0:imax)
    DIMENSION :: v(n,imax+1)
    DIMENSION :: bv(n,0:imax+1)
!
!   Local arguments
    INTEGER :: i,j,k,info,seed_sz,seed
    DOUBLE PRECISION :: dtmp,sub,y,work
!
    DIMENSION :: sub(IMAX-1)
    DIMENSION :: work(2*IMAX-2)

    ALLOCATABLE :: seed
    DIMENSION :: seed(:)
!
!-----------------------------------------------------------------------
!   Initialize
    IF(its.EQ.1) THEN
        IF(ijob.EQ.1) THEN
            bv(:,0) = 0.0D0
!           v(1):=random vector
            CALL RANDOM_SEED(SIZE=seed_sz)
            ALLOCATE(seed(seed_sz))
            seed=1
            CALL RANDOM_SEED(PUT=seed)
            CALL RANDOM_NUMBER(v(:,1))
            v(:,1)=v(:,1)*2.0D0-1.0D0
            DEALLOCATE(seed)
!           Return to the caller for matvec bv(1):=B*v(1) >>>>>>>>>>>>>>
            ijob=2
            RETURN
        END IF
!>>>>>>>Return from the caller after matvec bv(1):=B*v(1)
        IF(ijob.EQ.2) THEN
!          beta(0):=|| v(1) ||_B
            CALL kep_innpro(n,v(:,1),bv(:,1),dtmp)
            IF(dtmp.LE.0.0D0) THEN
                WRITE(6,*) &
                & 'Subroutine kep_lanczos: Error in B-inner product'
                WRITE(6,*) 'beta(0) =', dtmp
                STOP
            END IF
            beta(0)=DSQRT(dtmp)
            ijob=10
        END IF
    END IF
!
!-----------------------------------------------------------------------
!   Main loop
    DO j=its,its
        IF(ijob.EQ.10) THEN
!           Normalize
            v(:,j)=v(:,j)/beta(j-1)
            bv(:,j)=bv(:,j)/beta(j-1)
!           Return to the caller for matvec bv(j+1):=A*v(j) >>>>>>>>>>>>
            ijob=20
            RETURN
        END IF
!>>>>>>>Return from the caller after matvec bv(j+1):=A*v(j)
        IF(ijob.EQ.20) THEN
!           Orthgonalize (Modified Gram-Schmidt)
!           bv(j+1):=bv(j+1)-B*v(j-1)*beta(j-1)
            bv(:,j+1)=bv(:,j+1)-bv(:,j-1)*beta(j-1)
!           alpha(j):=(v(j),bv(j+1))
            CALL kep_innpro(n,v(:,j),bv(:,j+1),dtmp)
            alpha(j)=dtmp
!           bv(j+1):=bv(j+1)-B*v(j)*alpha(j)
            bv(:,j+1)=bv(:,j+1)-bv(:,j)*alpha(j)
!           Return to the caller to solve linear equation >>>>>>>>>>>>>>
!               v(j+1):=B^(-1)*bv(j+1), B*v(j+1)=bv(j+1) for v(j+1)
            ijob=30
            RETURN
        END IF
!>>>>>>>Return from the caller after solving linear equation
        IF(ijob.EQ.30) THEN
!           beta(j):=|| v(j+1) ||_B
            CALL kep_innpro(n,v(:,j+1),bv(:,j+1),dtmp)
            IF(dtmp.LE.0.0D0) THEN
                WRITE(6,*) 'Subroutine kep_lanczos: Error in B-inner product'
                WRITE(6,*) 'j =', j, ',  beta(j) =', dtmp
                STOP
            END IF
            beta(j)=DSQRT(dtmp)
!           Eigenpair of tridiagonal matrix by LAPACK
            theta(1:j) = alpha(1:j)
            sub(1:j-1) = beta(1:j-1)
            CALL dstev('N',j,theta,sub,y,1,work,info)
            IF(info.NE.0) THEN
                WRITE(6,*) 'Subroutine kep_lanczos: Error in LAPACK dstev routine'
                WRITE(6,*) 'j =', j, ',  info =', info
                STOP
            END IF
            ijob=10
        END IF
    END DO
!
!-----------------------------------------------------------------------
    its=its+1
    RETURN
    END SUBROUTINE kep_lanczos
