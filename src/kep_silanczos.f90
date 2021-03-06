    SUBROUTINE kep_silanczos(n,imax,its,ijob,v,bv,alpha,beta,theta,y,m)
    IMPLICIT NONE
    INTEGER :: n,imax,its,ijob,m
    DOUBLE PRECISION :: v,bv,alpha,beta,theta,y
!
    DIMENSION :: alpha(imax),theta(imax)
    DIMENSION :: beta(0:imax)
    DIMENSION :: y(IMAX,IMAX)
    DIMENSION :: v(n,imax+1)
    DIMENSION :: bv(n,0:imax+1)
!
!   Local arguments
    INTEGER :: i,j,k,info,seed
    DOUBLE PRECISION :: dtmp,sub,work
!
    DIMENSION :: sub(IMAX-1)
    DIMENSION :: work(2*IMAX-2)
!
!-----------------------------------------------------------------------
!   Initialize
    IF(its.EQ.1) THEN
        IF(ijob.EQ.1) THEN
            bv(:,0) = 0.0D0
!           bv(1):=random vector
            seed=1
            CALL kep_randnum(seed,n,bv(:,1))
            bv(:,1)=bv(:,1)*2.0D0-1.0D0
!           Return to the caller for matvec v(1):=B*bv(1) >>>>>>>>>>>>>>
            ijob=2
            RETURN
        END IF
!>>>>>>>Return from the caller after matvec v(1):=B*bv(1)
        IF(ijob.EQ.2) THEN
!          beta(0):=|| v(1) ||_B^(-1)
            CALL kep_innpro(n,v(:,1),bv(:,1),dtmp)
            IF(dtmp.LE.0.0D0) THEN
                WRITE(6,*) &
                & 'Subroutine kep_silanczos: Error in B^(-1)-inner product'
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
!           Return to the caller to solve linear equation >>>>>>>>>>>>>>
!               bv(j+1):=(A-sigma*B)^(-1)*v(j+1),
!               (A-sigma*B)*bv(j+1)=v(j+1) for bv(j+1)
            ijob=20
            RETURN
        END IF
!>>>>>>>Return from the caller after solving linear equation
        IF(ijob.EQ.20) THEN
!           Orthgonalize (Modified Gram-Schmidt)
!           bv(j+1):=bv(j+1)-B*v(j-1)*beta(j-1)
            bv(:,j+1)=bv(:,j+1)-bv(:,j-1)*beta(j-1)
!           alpha(j):=(v(j),bv(j+1))
            CALL kep_innpro(n,v(:,j),bv(:,j+1),dtmp)
            alpha(j)=dtmp
!           bv(j+1):=bv(j+1)-B*v(j)*alpha(j)
            bv(:,j+1)=bv(:,j+1)-bv(:,j)*alpha(j)
!           Fully reorthogonalize (Modified Gram-Schmidt)
            DO k=1,j
!               dtmp:=(v(k),bv(j+1))
                CALL kep_innpro(n,v(:,k),bv(:,j+1),dtmp)
!               bv(j+1):=bv(j+1)-BV(k)*dtmp
                bv(:,j+1)=bv(:,j+1)-bv(:,k)*dtmp
            END DO
!           alpha(j):=alpha(j)+dtmp
            alpha(j)=alpha(j)+dtmp
!           Return to the caller for matvec bv(j+1):=B*v(j) >>>>>>>>>>>>
            ijob=30
            RETURN
        END IF
!>>>>>>>Return from the caller after matvec bv(j+1):=B*v(j)
        IF(ijob.EQ.30) THEN
!           beta(j):=|| v(j+1) ||_B^(-1)
            CALL kep_innpro(n,v(:,j+1),bv(:,j+1),dtmp)
            IF(dtmp.LE.0.0D0) THEN
                WRITE(6,*) 'Subroutine kep_silanczos: Error in B^(-1)-inner product'
                WRITE(6,*) 'j =', j, ',  beta(j) =', dtmp
                STOP
            END IF
            beta(j)=DSQRT(dtmp)
!           Eigenpair of tridiagonal matrix by LAPACK
!           Skip computation when subspace is smaller than m
            IF(j.GE.m) THEN
                theta(1:j) = alpha(1:j)
                sub(1:j-1) = beta(1:j-1)
                CALL dstev('V',j,theta,sub,y,imax,work,info)
                IF(info.NE.0) THEN
                    WRITE(6,*) 'Subroutine kep_silanczos: Error in LAPACK dstev routine'
                    WRITE(6,*) 'j =', j, ',  info =', info
                    STOP
                END IF
            END IF
            ijob=10
        END IF
    END DO
!
!-----------------------------------------------------------------------
    its=its+1
    RETURN
    END SUBROUTINE kep_silanczos
