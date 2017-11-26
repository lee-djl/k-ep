    SUBROUTINE kep_mtx2crs(n,nz,indx,jndx,rval,row_pntr,col_indx,a)
    IMPLICIT NONE
    INTEGER :: n,nz,indx,jndx,row_pntr,col_indx
    DOUBLE PRECISION :: rval,a
    DIMENSION :: row_pntr(n+1)
    DIMENSION :: indx(nz),jndx(nz),rval(nz),col_indx(nz),a(nz)
!
!   local arguments
    INTEGER :: temp_indx,temp_jndx
    DOUBLE PRECISION :: temp_rval
    DIMENSION :: temp_indx(nz),temp_jndx(nz),temp_rval(nz)
!
    INTEGER :: i,k,r,c,nzd
!
!-----------------------------------------------------------------------
!   Initialize
    temp_indx=indx
    temp_jndx=jndx
    temp_rval=rval
!
    row_pntr=0
    nzd=0
!
!   Extract the main diagonal
    DO i=1,nz
        r=temp_indx(i)
        c=temp_jndx(i)
        IF(r.NE.c) THEN
            k=i-nzd
            temp_indx(k)=r
            temp_jndx(k)=c
            temp_rval(k)=temp_rval(i)
            !Nonzero elements in each row
            row_pntr(r)=row_pntr(r)+1
        ELSE
            nzd=nzd+1
            col_indx(nzd)=c
            a(nzd)=temp_rval(i)
        END IF
    END DO
!
!   Strating pointer of each row
    k=nzd+1
    DO i=1,n+1
       r=row_pntr(i)
       row_pntr(i)=k
       k=k+r
    END DO
!
!   Substitute elements into A and col_indx
    DO i=1,nz-nzd
        r=temp_indx(i)
        k=row_pntr(r)
        col_indx(k)=temp_jndx(i)
        a(k)=temp_rval(i)
        row_pntr(r)=k+1
    END DO
!
!   Modifying pointer
    DO i=n,1,-1
        row_pntr(i+1)=row_pntr(i)
    END DO
    row_pntr(1) = nzd+1
!
    RETURN
    END SUBROUTINE kep_mtx2crs
!
!=======================================================================
    SUBROUTINE kep_matvec(n,nz,row_pntr,col_indx,a,x,y)
    IMPLICIT NONE
    INTEGER :: n,nz,row_pntr,col_indx
    DOUBLE PRECISION :: a,x,y
!
    DIMENSION :: x(n),y(n)
    DIMENSION :: row_pntr(n+1)
    DIMENSION :: col_indx(nz),a(nz)
!
    INTEGER :: i,j,k
!
!-----------------------------------------------------------------------
!   Initialize
    y=0.0D0
!
!   y:=A*x (strictly lower triangular)
    DO i=2,n
        DO k=row_pntr(i),row_pntr(i+1)-1
            y(i)=y(i)+a(k)*x(col_indx(k))
        END DO
    END DO
!
!   y:=A*x (strictly upper triangular)
    DO j=2,n
        DO k=row_pntr(j),row_pntr(j+1)-1
            y(col_indx(k))=y(col_indx(k))+a(k)*x(j)
        END DO
    END DO
!
!   y:=A*x (diagonal)
    DO k=1,row_pntr(1)-1
        i=col_indx(k)
        y(i)=y(i)+a(k)*x(i)
    END DO
!
    RETURN
    END SUBROUTINE kep_matvec
!
!=======================================================================
    SUBROUTINE kep_innpro(n,x,y,alpha)
    IMPLICIT NONE
    INTEGER :: n
    DOUBLE PRECISION :: x,y,alpha
!
    DIMENSION :: x(n),y(n)
!
    INTEGER :: i
!
!-----------------------------------------------------------------------
!   Initialize
    alpha=0.0D0
!
    DO i=1,n
        alpha=alpha+x(i)*y(i)
    END DO
!
    RETURN
    END SUBROUTINE kep_innpro
!
!=======================================================================
    SUBROUTINE kep_convtest(n,imax,mmax,tol_res,tol_dif,its_si,its_test,&
    &   ijob,v,bv,beta,theta,&
    &   y,m,lngth,shift,&
    &   egnval,egnvec,egnvec_old,resnorm,difnorm,info)
    IMPLICIT NONE
    INTEGER :: n,imax,mmax,its_si,its_test,ijob,m,info
    DOUBLE PRECISION :: v,bv,beta,theta,y,lngth,shift,&
    &   egnval,egnvec,egnvec_old,resnorm,difnorm
!
    DIMENSION :: egnval(mmax),resnorm(mmax),difnorm(mmax)
    DIMENSION :: theta(imax)
    DIMENSION :: beta(0:imax)
    DIMENSION :: y(IMAX,IMAX)
    DIMENSION :: v(n,imax+1)
    DIMENSION :: bv(n,0:imax+1)
    DIMENSION :: egnvec(n,mmax+2),egnvec_old(n,mmax)
!
!   Local arguments
    INTEGER :: i,j,k,l,theta_indx,indx
    DOUBLE PRECISION :: theta_abs,bnd,tol_res,tol_dif,v_nrm,x,x_nrm,dtmp
!
    DIMENSION :: theta_indx(imax),theta_abs(imax),x_nrm(mmax)
    DIMENSION :: bnd(0:2*mmax+1)
    DIMENSION :: x(n)
!
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
    IF(ijob.EQ.1) THEN
!       Initialize
        DO i=1,its_si
            theta_indx(i)=i
        END DO
        theta_abs(1:its_si)=DABS(theta(1:its_si))
!
!-----------------------------------------------------------------------
!       Sort approximate eigenvalues in increasing order using modulus
        CALL kep_qcksrt(its_si,theta_abs(1:its_si),theta_indx(1:its_si))
!
!       Check whether appoximate eigenvalues are within an interval
        IF(theta_abs(its_si-m+1)*lngth.LT.1.0D0) THEN
            info=1
            RETURN
        END IF
!-----------------------------------------------------------------------
!       Sort approximate eigenvalues in the interval in increasing order
        DO i=1,m
            egnval(i)=(theta(theta_indx(its_si-m+i)))**(-1)
        END DO
        CALL kep_qcksrt(m,egnval(1:m),theta_indx(its_si-m+1:its_si))
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
            CALL kep_innpro(n,x,x,x_nrm(j))
            egnvec(:,j)=x(:)/DSQRT(x_nrm(j))
        END DO
!-----------------------------------------------------------------------
!       v_nrm:=( || v(its_si+1)*beta(its_si) ||_2 )^2
        CALL kep_innpro(n,v(:,its_si+1),v(:,its_si+1),v_nrm)
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
                RETURN
            END IF
        END DO
!-----------------------------------------------------------------------
!
        ijob=2
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
            RETURN
        END IF
!>>>>>>>Return from the caller after matvec
        IF(ijob.EQ.3) THEN
!           Relative residual 2-norm (squared)
            x(:)=egnvec(:,mmax+1)-egnvec(:,mmax+2)*egnval(j)
            CALL kep_innpro(n,x,x,dtmp)
!
            IF(dtmp.GT.tol_res) THEN
                info=4
                RETURN
            END IF
!
            resnorm(j)=DSQRT(dtmp)
            ijob=2
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
                RETURN
            END IF
!
            x(:)=egnvec(:,j)*dtmp-egnvec_old(:,j)
!
            CALL kep_innpro(n,x,x,dtmp)
            IF(dtmp.GT.tol_dif) THEN
                info=6
                RETURN
            END IF
            difnorm(j)=DSQRT(dtmp)
        END DO
    END IF
!
!-----------------------------------------------------------------------
!
    RETURN
    END SUBROUTINE kep_convtest
!
!=======================================================================
!   Quick sort in increasing order
!
!   Modification of the code on the following URL
!   (https://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Fortran)
    RECURSIVE SUBROUTINE kep_qcksrt(n,x,indx)
    IMPLICIT NONE
    INTEGER :: n,indx
    DOUBLE PRECISION :: x
!
    DIMENSION :: indx(n),x(n)
!
    INTENT(IN) :: n
    INTENT(INOUT) :: x,indx
!
!   Local arguments
    INTEGER :: left,right,itmp,marker
    DOUBLE PRECISION :: random,pivot,dtmp
!
!-----------------------------------------------------------------------
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
       CALL kep_qcksrt(marker-1,x(:marker-1),indx(:marker-1))
       CALL kep_qcksrt(n-marker+1,x(marker:),indx(marker:))
    END IF
    END SUBROUTINE kep_qcksrt
!
!=======================================================================
    SUBROUTINE kep_summary(mmax,icnt,cmpt_time,&
    &   si_1,si_2_l,si_2_r,si_3,egnval,resnorm,difnorm,info)
    IMPLICIT NONE
    INTEGER :: ounit,mmax,icnt,info
    REAL :: cmpt_time
    DOUBLE PRECISION :: egnval,resnorm,difnorm
!
    DIMENSION :: icnt(10),cmpt_time(10)
    DIMENSION :: egnval(mmax),resnorm(mmax),difnorm(mmax)
!
    TYPE :: shift_inertia
        SEQUENCE
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: shift
        INTEGER, ALLOCATABLE, DIMENSION(:) :: inertia
    END TYPE shift_inertia
!
    TYPE(shift_inertia) :: si_1,si_2_l,si_2_r,si_3
!
!   Local arguments
    INTEGER :: i
    REAL :: rtmp
!
!-----------------------------------------------------------------------
!
    ounit=6
!
    WRITE(ounit,*)
    WRITE(ounit,*)'================================================================================'
    WRITE(ounit,*)
    WRITE(ounit,*)'-----Computation time (s)'
    rtmp=0.0E0
    DO i=1,10
        rtmp=rtmp+cmpt_time(i)
    END DO
    WRITE(ounit,"(A39,F15.7,/)") 'Overall                               :',rtmp
    WRITE(ounit,"(A39,F15.7)") 'Step 1 Task (i)   symbolic factor. B  :',cmpt_time(1)
    WRITE(ounit,"(A39,F15.7)") 'Step 1 Task (ii)  numer. factor. B    :',cmpt_time(2)
    WRITE(ounit,"(A39,F15.7)") 'Step 1 Task (iii) Lanczos             :',cmpt_time(3)
    WRITE(ounit,"(A39,F15.7)") 'Step 1 Task (iv)  inertia computation :',cmpt_time(4)
    WRITE(ounit,"(A39,F15.7)") 'Step 2 Task (v)   inertia computation :',cmpt_time(5)
    WRITE(ounit,"(A39,F15.7)") 'Step 3 Task (vi)  numer. factor. A-sB :',cmpt_time(6)
    WRITE(ounit,"(A39,F15.7,/)") 'Step 3 Task (vii) SI Lanczos          :',cmpt_time(7)
!
    WRITE(ounit,*)'-----Iteration count'
    WRITE(ounit,"(A8,I4)") 'Step 1 :',icnt(1)
    WRITE(ounit,"(A8,I4)") 'Step 2 :',icnt(2)
    WRITE(ounit,"(A8,I4,/)") 'Step 3 :',icnt(4)
!
    WRITE(ounit,*)'-----Step 1: set an initial interval'
    WRITE(ounit,"(A9,2X,A8,2X,A28)")'Iteration','Inertia','Shift                       '
    DO i=1,icnt(1)
        WRITE(ounit,"(I9,2X,I8,2X,E28.18)") &
        &   i,si_1%inertia(i),si_1%shift(i)
    END DO
!
    IF(info.EQ.1) THEN
        RETURN
    END IF
!
    WRITE(ounit,*)
    WRITE(ounit,*) '-----Step 2: narrow down the interval'
    WRITE(ounit,"(A9,2X,A12,2X,A12,2X,A28,2X,A28)")'Iteration','InertiaLower','InertiaUpper',&
    &   'ShiftLower                  ','ShiftUpper                  '
    DO i=0,icnt(2)
        WRITE(ounit,"(I9,2X,I12,2X,I12,2X,E28.18,2X,E28.18)") &
        &   i,si_2_l%inertia(i),si_2_r%inertia(i),&
        &   si_2_l%shift(i),si_2_r%shift(i)
    END DO
!
    IF(info.EQ.2) THEN
        RETURN
    END IF
!
    WRITE(ounit,*)
    WRITE(ounit,*) '-----Step 3: select the midpoint'
    WRITE(ounit,"(A8,2X,A28)")'Inertia','Shift                       '
    WRITE(ounit,"(I8,2X,E28.18)") &
    &   si_3%inertia(1),si_3%shift(1)
!
    IF(info.EQ.3) THEN
        RETURN
    END IF
!
    WRITE(ounit,*)
    WRITE(ounit,*) '-----Step 3: compute eipenpairs of the interval'
    WRITE(ounit,"(A8,2X,A28,2X,A28,2X,A28)")'Index   ','Eigenvalue                  ',&
    &   'Relative residual 2-norm    ','Relative difference 2-norm    '
    DO i=1,icnt(3)
        WRITE(ounit,"(I8,2X,E28.18,2X,E28.18,2X,E28.18)") &
        &   si_2_l%inertia(icnt(2))+i,egnval(i),resnorm(i),difnorm(i)
    END DO
!
    RETURN
    END SUBROUTINE kep_summary
!
!=======================================================================
    SUBROUTINE kep_checkprm(n,k,iprm,info)
    IMPLICIT NONE
    INTEGER :: n,k,iprm,info
    DIMENSION :: iprm(30)
!
!-----------------------------------------------------------------------
!
    IF((k.LE.1).OR.(k.GE.n)) THEN
        info=-11
    ELSE IF((iprm(1).LT.1).OR.(iprm(1).GE.n)) THEN
        info=-12
    ELSE IF((iprm(2).LT.1).OR.(iprm(2).GT.16)) THEN
        info=-13
    ELSE IF((iprm(3).LT.1).OR.(iprm(3).GT.16)) THEN
        info=-14
    ELSE IF((iprm(11).LT.2).OR.(iprm(11).GT.n)) THEN
        info=-15
    ELSE IF((iprm(12).LT.0).OR.(iprm(12).GT.64)) THEN
        info=-16
    ELSE IF((iprm(13).LT.iprm(1)).OR.(iprm(13).GT.n)) THEN
        info=-17
    END IF
!
!-----------------------------------------------------------------------
    RETURN
    END SUBROUTINE kep_checkprm
!
!=======================================================================
    SUBROUTINE kep_info(info)
    IMPLICIT NONE
    INTEGER :: info
!
!-----------------------------------------------------------------------
!
    WRITE(6,*)
    IF(info.EQ.-1) THEN
        WRITE(6,*)'=====kep info',info,': Allocate matrices'    
    ELSE IF(info.EQ.-11) THEN
        WRITE(6,*)'=====kep info',info,': k must satisfy 1 < k < n'
    ELSE IF(info.EQ.-12) THEN
        WRITE(6,*)'=====kep info',info,': iprm(1) must satisfy 1 <= iprm(1) < n'
    ELSE IF(info.EQ.-13) THEN
        WRITE(6,*)'=====kep info',info,': iprm(2) must satisfy 1 <= iprm(2) <= 16'
    ELSE IF(info.EQ.-14) THEN
        WRITE(6,*)'=====kep info',info,': iprm(3) must satisfy 1 <= iprm(3) <= 16'
    ELSE IF(info.EQ.-15) THEN
        WRITE(6,*)'=====kep info',info,': iprm(11) must satisfy 2 <= iprm(11) <= n'
    ELSE IF(info.EQ.-16) THEN
        WRITE(6,*)'=====kep info',info,': iprm(12) must satisfy 0 <= iprm(12) <= 64'
    ELSE IF(info.EQ.-17) THEN
        WRITE(6,*)'=====kep info',info,': iprm(13) must satisfy iprm(1) <= iprm(13) <= n'
    ELSE IF(info.EQ.1) THEN
        WRITE(6,*)'=====kep info',info,': Failed to set an initial interval'
    ELSE IF(info.EQ.2) THEN
        WRITE(6,*)'=====kep info',info,': Failed to narrow down the interval'
    ELSE IF(info.EQ.3) THEN
        WRITE(6,*)'=====kep info',info,': Failed to compute eipenpairs of the interval'
    ELSE IF((info.GE.11).AND.(info.LT.20)) THEN
        WRITE(6,*)'=====kep info',info,': Error in MUMPS'
    END IF
!
!-----------------------------------------------------------------------
    RETURN
    END SUBROUTINE kep_info

