    MODULE k2l_lanczos_mod
    USE k2l_inttype
    USE k2l_utility
    USE k2l_solver
    CONTAINS
!=======================================================================
    SUBROUTINE k2l_lanczos(n,imax,k2l_l)
    IMPLICIT NONE
    INTEGER :: n,imax
    TYPE(k2l_lanczos_type) :: k2l_l
!----------------------------------------------------------------------
    CALL k2l_lanczos2(n,imax,k2l_l%its,k2l_l%ijob,k2l_l%v,k2l_l%bv,k2l_l%alpha,k2l_l%beta,k2l_l%theta)
!    
    RETURN
    END SUBROUTINE k2l_lanczos
!=======================================================================
    SUBROUTINE k2l_lanczos2(n,imax,its,ijob,v,bv,alpha,beta,theta)
    IMPLICIT NONE
    INTEGER :: n,imax,its,ijob
    DOUBLE PRECISION :: v,bv,alpha,beta,theta
!
    DIMENSION :: alpha(:),theta(:)
    DIMENSION :: beta(0:)
    DIMENSION :: v(:,:)
    DIMENSION :: bv(:,0:)
!
!   Local arguments
    INTEGER :: i,j,k,info,seed
    DOUBLE PRECISION :: dtmp,y
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: sub, work
!
!----------------------------------------------------------------------
!   Initialize
    IF(its.EQ.1) THEN
        IF(ijob.EQ.1) THEN
            bv(:,0) = 0.0D0
!           v(1):=random vector
            seed=1
            CALL k2l_randnum(seed,n,v(:,1))
            v(:,1)=v(:,1)*2.0D0-1.0D0
!           Return to the caller for matvec bv(1):=B*v(1) >>>>>>>>>>>>>>
            ijob=2
            RETURN
        END IF
!>>>>>>>Return from the caller after matvec bv(1):=B*v(1)
        IF(ijob.EQ.2) THEN
!          beta(0):=|| v(1) ||_B
            CALL k2l_innpro(n,v(:,1),bv(:,1),dtmp)
            IF(dtmp.LE.0.0D0) THEN
                WRITE(6,*) '=====k2l: Error in subroutine k2l_lanczos (zero B-inner product)'
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
            CALL k2l_innpro(n,v(:,j),bv(:,j+1),dtmp)
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
            CALL k2l_innpro(n,v(:,j+1),bv(:,j+1),dtmp)
            IF(dtmp.LE.0.0D0) THEN
                WRITE(6,*) '=====k2l: Error in subroutine k2l_lanczos (zero B-inner product)'
                WRITE(6,*) 'j =', j, ',  beta(j) =', dtmp
                STOP
            END IF
            beta(j)=DSQRT(dtmp)
!           Eigenpair of tridiagonal matrix by LAPACK
            ALLOCATE(sub(imax-1),work(2*imax-2))
            theta(1:j) = alpha(1:j)
            sub(1:j-1) = beta(1:j-1)
            CALL dstev('N',j,theta,sub,y,1,work,info)
            DEALLOCATE(sub,work)
            IF(info.NE.0) THEN
                WRITE(6,*) '=====k2l: Error in subroutine k2l_lanczos (LAPACK dstev failed)'
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
    END SUBROUTINE k2l_lanczos2
!=======================================================================
    SUBROUTINE k2l_silanczos(n,imax,k2l_s)
    IMPLICIT NONE
    INTEGER :: n,imax
    TYPE(k2l_silanczos_type) :: k2l_s
!----------------------------------------------------------------------
    CALL k2l_silanczos2(n,imax,k2l_s%its,k2l_s%ijob,k2l_s%v,k2l_s%bv,&
    &   k2l_s%alpha,k2l_s%beta,k2l_s%theta,k2l_s%y,k2l_s%icnt)
!    
    RETURN
    END SUBROUTINE k2l_silanczos
!=======================================================================
    SUBROUTINE k2l_silanczos2(n,imax,its,ijob,v,bv,alpha,beta,theta,y,m)
    IMPLICIT NONE
    INTEGER :: n,imax,its,ijob,m
    DOUBLE PRECISION :: v,bv,alpha,beta,theta,y
!
    DIMENSION :: alpha(:),theta(:)
    DIMENSION :: beta(0:)
    DIMENSION :: y(:,:)
    DIMENSION :: v(:,:)
    DIMENSION :: bv(:,0:)
!
!   Local arguments
    INTEGER :: i,j,k,info,seed
    DOUBLE PRECISION :: dtmp
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: sub, work
!
!-----------------------------------------------------------------------
!   Initialize
    IF(its.EQ.1) THEN
        IF(ijob.EQ.1) THEN
            bv(:,0) = 0.0D0
!           bv(1):=random vector
            seed=1
            CALL k2l_randnum(seed,n,bv(:,1))
            bv(:,1)=bv(:,1)*2.0D0-1.0D0
!           Return to the caller for matvec v(1):=B*bv(1) >>>>>>>>>>>>>>
            ijob=2
            RETURN
        END IF
!>>>>>>>Return from the caller after matvec v(1):=B*bv(1)
        IF(ijob.EQ.2) THEN
!          beta(0):=|| v(1) ||_B^(-1)
            CALL k2l_innpro(n,v(:,1),bv(:,1),dtmp)
            IF(dtmp.LE.0.0D0) THEN
                WRITE(6,*) '=====k2l: Error in subroutine k2l_silanczos (zero B^(-1)-inner product)'
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
            CALL k2l_innpro(n,v(:,j),bv(:,j+1),dtmp)
            alpha(j)=dtmp
!           bv(j+1):=bv(j+1)-B*v(j)*alpha(j)
            bv(:,j+1)=bv(:,j+1)-bv(:,j)*alpha(j)
!           Fully reorthogonalize (Modified Gram-Schmidt)
            DO k=1,j
!               dtmp:=(v(k),bv(j+1))
                CALL k2l_innpro(n,v(:,k),bv(:,j+1),dtmp)
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
            CALL k2l_innpro(n,v(:,j+1),bv(:,j+1),dtmp)
            IF(dtmp.LE.0.0D0) THEN
                WRITE(6,*) '=====k2l: Error in subroutine k2l_silanczos (zero B^(-1)-inner product)'
                WRITE(6,*) 'j =', j, ',  beta(j) =', dtmp
                STOP
            END IF
            beta(j)=DSQRT(dtmp)
!           Eigenpair of tridiagonal matrix by LAPACK
!           Skip computation when subspace is smaller than m
            IF(j.GE.m) THEN
                ALLOCATE(sub(imax-1),work(2*imax-2))
                theta(1:j) = alpha(1:j)
                sub(1:j-1) = beta(1:j-1)
                CALL dstev('V',j,theta,sub,y,imax,work,info)
                DEALLOCATE(sub,work)
                IF(info.NE.0) THEN
                    WRITE(6,*) '=====k2l: Error in subroutine k2l_silanczos (LAPACK dstev failed)'
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
    END SUBROUTINE k2l_silanczos2
!=======================================================================    
    SUBROUTINE k2l_convtest(imax_s,k2l_s,k2l_c,shift,i,j)
    IMPLICIT NONE
    INTEGER :: imax_s,i,j
    DOUBLE PRECISION :: shift
    TYPE(k2l_silanczos_type) :: k2l_s
    TYPE(k2l_convtest_type) :: k2l_c
!-----------------------------------------------------------------------
    CALL k2l_convtest2(SIZE(k2l_c%egnvec,1),imax_s,k2l_c%mmax,k2l_c%tol_res,k2l_c%tol_dif,i,j, &
    &   k2l_c%ijob,k2l_s%v,k2l_s%bv,k2l_s%beta,k2l_s%theta, &
    &   k2l_s%y,k2l_s%icnt,k2l_c%int_length,shift, &
    &   k2l_c%egnval,k2l_c%egnvec,k2l_c%egnvec_old,k2l_c%resnorm,k2l_c%difnorm,k2l_c%info)
!
    RETURN
    END SUBROUTINE k2l_convtest
!=======================================================================    
    SUBROUTINE k2l_convtest2(n,imax,mmax,tol_res,tol_dif,its_si,its_test,&
    &   ijob,v,bv,beta,theta,&
    &   y,m,lngth,shift,&
    &   egnval,egnvec,egnvec_old,resnorm,difnorm,info)
    IMPLICIT NONE
    INTEGER :: n,imax,mmax,its_si,its_test,ijob,m,info
    DOUBLE PRECISION :: v,bv,beta,theta,y,lngth,shift,&
    &   egnval,egnvec,egnvec_old,resnorm,difnorm
!
    DIMENSION :: egnval(:),resnorm(:),difnorm(:)
    DIMENSION :: theta(:)
    DIMENSION :: beta(0:)
    DIMENSION :: y(:,:)
    DIMENSION :: v(:,:)
    DIMENSION :: bv(:,0:)
    DIMENSION :: egnvec(:,:),egnvec_old(:,:)
!
!   Local arguments
    INTEGER :: i,j,k,l,indx
    DOUBLE PRECISION :: tol_res,tol_dif,v_nrm,dtmp
    INTEGER, ALLOCATABLE, DIMENSION(:) :: theta_indx
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: theta_abs,x_nrm,x,bnd
!-----------------------------------------------------------------------
!
    IF(its_si.LT.m) THEN
        info=-1
        RETURN
    END IF
!
    IF(info.NE.0) THEN
        RETURN
    END IF
!
!-----------------------------------------------------------------------
    ALLOCATE(x(n))
    IF(ijob.EQ.1) THEN
!       Initialize
        ALLOCATE(theta_indx(imax),theta_abs(imax),x_nrm(mmax),bnd(0:2*mmax+1))
        DO i=1,its_si
            theta_indx(i)=i
        END DO
        theta_abs(1:its_si)=DABS(theta(1:its_si))
!
!-----------------------------------------------------------------------
!       Sort approximate eigenvalues in increasing order using modulus
        CALL k2l_qcksrt(its_si,theta_abs(1:its_si),theta_indx(1:its_si))
!
!       Check whether appoximate eigenvalues are within an interval
        IF(theta_abs(its_si-m+1)*lngth.LT.1.0D0) THEN
            info=1
            DEALLOCATE(theta_indx,theta_abs,x_nrm,bnd,x)
            RETURN
        END IF
!-----------------------------------------------------------------------
!       Sort approximate eigenvalues in the interval in increasing order
        DO i=1,m
            egnval(i)=(theta(theta_indx(its_si-m+i)))**(-1)
        END DO
        CALL k2l_qcksrt(m,egnval(1:m),theta_indx(its_si-m+1:its_si))
!
!       Compute an error bound
        bnd(0)=-lngth
        bnd(2*m+1)=lngth
        DO i=1,m
            l=theta_indx(its_si-m+i)
            dtmp=DABS(beta(its_si)*y(its_si,l)/theta(l))
            dtmp=DSIGN(1.0D0,theta(l))*dtmp/DSQRT(dtmp**2+1.0D0)
            bnd(2*i-1)=(-dtmp+1.0D0)/theta(l)
            bnd(2*i)=(dtmp+1.0D0)/theta(l)
        END DO
!
!       Check whether bounds are included in the interval and disjoint
        DO i=0,m
            IF(bnd(2*i).GE.bnd(2*i+1)) THEN
                info=2
                DEALLOCATE(theta_indx,theta_abs,x_nrm,bnd,x)
                RETURN
            END IF
        END DO
!-----------------------------------------------------------------------
!       Get approximate eigenvalue
        DO j=1,m
            egnval(j)=egnval(j)+shift
        END DO
!
!       Get approxmate eigenvector
        DO j=1,m
            l=theta_indx(its_si-m+j)
!
!           Approximate eigenvector x
            x=0.0D0
            DO k=1,its_si
                x(:)=x(:)+bv(:,k)*y(k,l)
            END DO
            x=x*theta(l)
            x(:)=x(:)+bv(:,its_si+1)*y(its_si,l)
!
!           Divisor of relative residual 2-norm
            CALL k2l_innpro(n,x,x,x_nrm(j))
            egnvec(:,j)=x(:)/DSQRT(x_nrm(j))
        END DO
!-----------------------------------------------------------------------
!       v_nrm:=( || v(its_si+1)*beta(its_si) ||_2 )^2
        CALL k2l_innpro(n,v(:,its_si+1),v(:,its_si+1),v_nrm)
!
!       Pre-checking relative residual 2-norm
        DO j=1,m
            l=theta_indx(its_si-m+j)
!
!           Dividend of relative residual 2-norm (squared)
            dtmp=DABS(y(its_si,l)/theta(l))
            dtmp=(dtmp**2)*v_nrm
!
!           Relative residual 2-norm (squared)
            dtmp=dtmp/x_nrm(j)
            IF(dtmp.GT.tol_res) THEN
                info=3
                DEALLOCATE(theta_indx,theta_abs,x_nrm,bnd,x)
                RETURN
            END IF
        END DO
!-----------------------------------------------------------------------
!
        ijob=2
        DEALLOCATE(theta_indx,theta_abs,x_nrm,bnd,x)
        RETURN
    END IF
!-----------------------------------------------------------------------
!   Checking relative residual 2-norm with matvec operation
    DO j=its_test,its_test
        IF(ijob.EQ.2) THEN
!           Return to the caller for matvec
!           egnvec(:,mmax+1):=A*egnvec(:,mmax+1) >>>>>>>>>>>>>>>>>>>>>>>
!           egnvec(:,mmax+1):=B*egnvec(:,mmax+1) >>>>>>>>>>>>>>>>>>>>>>>
            ijob=3
            IF(ALLOCATED(x)) DEALLOCATE(x)
            RETURN
        END IF
!>>>>>>>Return from the caller after matvec
        IF(ijob.EQ.3) THEN
!           Relative residual 2-norm (squared)
            x(:)=egnvec(:,mmax+1)-egnvec(:,mmax+2)*egnval(j)
            CALL k2l_innpro(n,x,x,dtmp)
!
            IF(dtmp.GT.tol_res) THEN
                info=4
                IF(ALLOCATED(x)) DEALLOCATE(x)
                RETURN
            END IF
!
            resnorm(j)=DSQRT(dtmp)
            ijob=2
            IF(ALLOCATED(x)) DEALLOCATE(x)
            RETURN
        END IF
    END DO
!-----------------------------------------------------------------------
!   Checking relative difference 2-norm
    IF(ijob.EQ.4) THEN
        DO j=1,m
            IF(DABS(MAXVAL(egnvec_old(:,j))) &
            &   .ge.DABS(MINVAL(egnvec_old(:,j)))) THEN
                indx=MAXLOC(egnvec_old(:,j),1)
            ELSE
                indx=MINLOC(egnvec_old(:,j),1)
            END IF
!
            IF((egnvec_old(indx,j)*egnvec(indx,j)).GT.0.0D0) THEN
                dtmp=1.0D0
            ELSE IF((egnvec_old(indx,j)*egnvec(indx,j)).LT.0.0D0) THEN
                dtmp=-1.0D0
            ELSE
                info=5
                IF(ALLOCATED(x)) DEALLOCATE(x)
                RETURN
            END IF
!
            x(:)=egnvec(:,j)*dtmp-egnvec_old(:,j)
!
            CALL k2l_innpro(n,x,x,dtmp)
            IF(dtmp.GT.tol_dif) THEN
                info=6
                IF(ALLOCATED(x)) DEALLOCATE(x)
                RETURN
            END IF
            difnorm(j)=DSQRT(dtmp)
        END DO
    END IF
!
!-----------------------------------------------------------------------
!
    IF(ALLOCATED(x)) DEALLOCATE(x)
    RETURN
    END SUBROUTINE k2l_convtest2
!=======================================================================
    RECURSIVE SUBROUTINE k2l_qcksrt(n,x,indx)
!   Quick sort in increasing order
!   Modification of the code on the following URL
!   (https://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Fortran)
    IMPLICIT NONE
    INTEGER :: n
    INTEGER, DIMENSION(:) :: indx
    DOUBLE PRECISION, DIMENSION(:) :: x
!
    INTENT(IN) :: n
    INTENT(INOUT) :: x,indx
!
!   Local arguments
    INTEGER :: left,right,itmp,marker
    DOUBLE PRECISION :: random,pivot,dtmp
!
!-----------------------------------------------------------------------
    IF((n.NE.SIZE(indx)).OR.(n.NE.SIZE(x))) THEN
        WRITE(*,*) 'Incompatible array size'
        STOP
    END IF
!
    IF(n>1) THEN
        ! Random pivot
        CALL RANDOM_NUMBER(random)
        pivot=x(1+INT(random*DBLE(n-1)))
        left=0
        right=n+1
!
        DO WHILE(left<right)
            right=right-1
            DO WHILE(x(right)>pivot)
                right=right-1
            END DO
            left=left+1
            DO WHILE(x(left)<pivot)
                left=left+1
            END DO
            IF(left<right) THEN
                dtmp=x(left)
                x(left)=x(right)
                x(right)=dtmp
                itmp=indx(left)
                indx(left)=indx(right)
                indx(right)=itmp
            END IF
        END DO
!
       IF(left.EQ.right) THEN
          marker=left+1
       ELSE
          marker=left
       END IF
!
       CALL k2l_qcksrt(marker-1,x(:marker-1),indx(:marker-1))
       CALL k2l_qcksrt(n-marker+1,x(marker:),indx(marker:))
    END IF
    END SUBROUTINE k2l_qcksrt
!=======================================================================    
    SUBROUTINE k2l_randnum(seed,n,randarray)
    USE xorshift128plus
    IMPLICIT NONE
    INTEGER :: seed,n
    DOUBLE PRECISION, DIMENSION(:) :: randarray
!
    INTEGER :: i
!-----------------------------------------------------------------------
    IF(n.NE.SIZE(randarray)) THEN
        WRITE(*,*) 'Incompatible array size'
        STOP
    END IF

    CALL rand_init(seed)
!
    DO i=1,n
        randarray(i)=rand_uniform()
    END DO
! 
    RETURN
    END SUBROUTINE k2l_randnum    
!=======================================================================
    END MODULE k2l_lanczos_mod
