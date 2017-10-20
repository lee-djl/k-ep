    SUBROUTINE kep(n,&
    &   nz_a,indx_a,jndx_a,rval_a,&
    &   nz_b,indx_b,jndx_b,rval_b,&
    &   k,iprm,keval,kevec,info)
!=======================================================================
!
!   n       : <input> matrix size
!
!   nz_a    : <input> number of non-zero elements of matrix A
!   indx_a  : <input> row index of non-zero elements of matrix A
!   jndx_a  : <input> column index of non-zero elements of matrix A
!   rval_a  : <input> value of non-zero elements of matrix A
!
!   nz_b    : <input> number of non-zero elements of matrix B
!   indx_b  : <input> row index of non-zero elements of matrix B
!   jndx_b  : <input> column index of non-zero elements of matrix B
!   rval_b  : <input> value of non-zero elements of matrix B
!
!   k       : <input> target index
!             Required to satisfy 1 < k < n
!
!   iprm    : <input> parameters for kep
!   iprm(1) : stopping criterion for bisection. See m_max in reference [1]
!             ( = maximum number of eigenpairs computed together)
!             Required to satisfy 1 <= iprm(1) < n
!   iprm(2) : 10**-iprm(2) is tolerance for relative residual 2-norm. See tau_res in [1]
!             Required to satisfy 1 <= iprm(2) <= 16
!   iprm(3) : 10**-iprm(3) is tolerance for relative difference 2-norm. See tau_diff in [1]
!             Required to satisfy 1 <= iprm(3) <= 16
!   iprm(12): maximum iteration count to set an initial interval. See Algorithm 2 in [1]
!             ( = maximum iteration count for Lanczos)
!             Required to satisfy 2 <= iprm(11) <= n
!   iprm(12): maximum iteration count to narrow down the interval. See Algorithm 1' in [1]
!             ( = maximum iteration count for Bisection)
!             Required to satisfy 0 <= iprm(12) <= 64
!   iprm(13): maximum iteration count to compute eigenpairs of the interval. See Algorithm 3 in [1]
!             ( = maximum iteration count for shift-and-invert Lanczos)
!             Required to satisfy iprm(1) <= iprm(13) <= n
!   iprm(21): Suppress detailed output to terminal if iprm(21) < 0
!
!   kevec   : <output> k-th eigenvalue
!   kevec   : <output> k-th eigenvector normalized with respect to 2-norm
!
!   info    : <output> diagonostic information
!   info = 0: successful exit
!   info < 0: illegal value in k, iprm
!   info > 0: failed to find the k-th eigenpair
!
!-----------------------------------------------------------------------
!
!   References
!   [1] D. Lee, T. Hoshi, T. Sogabe, Y. Miyatake, S.-L. Zhang,
!   Solution of the k-th eigenvalue problem in large-scale electronic structure calculations,
!   https://arxiv.org/abs/1710.05134.
!
!-----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER ::              &
    &   n,                  &
    &   nz_a,indx_a,jndx_a, &
    &   nz_b,indx_b,jndx_b, &
    &   k,iprm,info
    DOUBLE PRECISION ::     &
    &   rval_a,rval_b,      &
    &   keval,kevec
    DIMENSION ::            &
    &   indx_a(nz_a),jndx_a(nz_a),rval_a(nz_a),&
    &   indx_b(nz_b),jndx_b(nz_b),rval_b(nz_b),&
    &   iprm(30),kevec(n)
!
!=======================================================================
!
!   Convert data to Compressed Row Storage (CRS) format
    INTEGER :: row_pntr_a,col_indx_a,row_pntr_b,col_indx_b
    DOUBLE PRECISION :: a,b
    DIMENSION :: row_pntr_a(n+1),col_indx_a(nz_a),a(nz_a),  &
    &   row_pntr_b(n+1),col_indx_b(nz_b),b(nz_b)
!
!   MUMPS
    INCLUDE 'mpif.h'
    INCLUDE 'dmumps_struc.h'
    INTEGER :: ierr,m_ord
    TYPE (DMUMPS_STRUC) :: m_1,m_2,m_3
    ALLOCATABLE :: m_ord
    DIMENSION :: m_ord(:)
!
!   Lanczos
    INTEGER :: imax_lanc,its_lanc,ijob
    DOUBLE PRECISION :: v,bv,alpha,beta,theta
    ALLOCATABLE :: v,bv,alpha,beta,theta
    DIMENSION :: alpha(:),theta(:),beta(:)
    DIMENSION :: v(:,:),bv(:,:)
!
!   Bisection
    INTEGER :: imax_bi,mmax
!
!   Shift-and-invert Lanczos
    INTEGER :: imax_si,its_si,ijob_si,ijob_test,info_test
    DOUBLE PRECISION :: tol_res,tol_dif,                    &
    &   v_si,bv_si,alpha_si,beta_si,theta_si,               &
    &   y,lngth,egnval,egnvec,egnvec_old,res,errnorm
    ALLOCATABLE :: v_si,bv_si,alpha_si,beta_si,theta_si,y,  &
    &   egnval,egnvec,egnvec_old,res,errnorm
    DIMENSION :: alpha_si(:),theta_si(:),beta_si(:),        &
    &   egnval(:),res(:),errnorm(:)
    DIMENSION :: v_si(:,:),bv_si(:,:),y(:,:),egnvec(:,:),egnvec_old(:,:)
!
!   Time
    INTEGER :: clock_1,clock_2,clock_rate,clock_max
    REAL :: cmpt_time
    DIMENSION :: cmpt_time(10)
!
!   Results
    INTEGER :: icnt
    DIMENSION :: icnt(10)
    TYPE shift_inertia
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: shift
        INTEGER, ALLOCATABLE, DIMENSION(:) :: inertia
    END TYPE shift_inertia
    TYPE(shift_inertia) :: si_1,si_2_l,si_2_r,si_3
!
!   Local arguments
    INTEGER :: i,j
    DOUBLE PRECISION :: dtmp
!
!---Initialize----------------------------------------------------------
!
!   Check parameters
    info=0
    CALL kep_checkprm(n,k,iprm,info)
    IF(info.NE.0) THEN
        CALL kep_info(info)
        RETURN
    END IF
!
!   Initialize counts and info
    icnt=0
    info_test=0
!
!   Lanczos
    imax_lanc=iprm(11)
    its_lanc=1
    ijob=1
!
!   Bisection
    imax_bi=iprm(12)
    mmax=iprm(1)
!
!   Shift-and-invert Lanczos
    imax_si=iprm(13)
    its_si=1
    ijob_si=1
    ijob_test=1
    tol_res=1.0D1**(-REAL(2*iprm(2)))
    tol_dif=1.0D1**(-REAL(2*iprm(3)))
!
!   Time and results
    cmpt_time=0.0E0
    keval=0.0D0
    kevec=0.0D0
!
    ALLOCATE(si_1%shift(0:imax_lanc),si_1%inertia(0:imax_lanc))
    ALLOCATE(si_2_l%shift(0:imax_bi),si_2_l%inertia(0:imax_bi), &
    &   si_2_r%shift(0:imax_bi),si_2_r%inertia(0:imax_bi))
    ALLOCATE(si_3%shift(1),si_3%inertia(1))
    ALLOCATE(egnval(mmax),egnvec(n,mmax+2),                     &
    &   egnvec_old(n,mmax+2),res(mmax),errnorm(mmax))
!
!   Convert data to Compressed Row Storage (CRS) format
    CALL kep_mtx2crs(n,nz_a,indx_a,jndx_a,rval_a,row_pntr_a,col_indx_a,a)
    CALL kep_mtx2crs(n,nz_b,indx_b,jndx_b,rval_b,row_pntr_b,col_indx_b,b)
!
!   MUMPS
    CALL mpi_init(ierr)
    INCLUDE 'm_1'
    INCLUDE 'm_2'
!
!---Step 1 : set an initial interval------------------------------------
!
!   Task (i) : symbolic factorizaion of B
    CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
    m_1%job=1
    CALL dmumps(m_1)
    IF(m_1%infog(1).NE.0) THEN
        info=11
        CALL kep_info(info)
        RETURN
    END IF
    CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
    cmpt_time(1)=REAL(clock_2-clock_1)/REAL(clock_rate)
!
!   Recycle and reuse fill-reducing ordering
    m_ord=m_1%sym_perm
    m_2%perm_in=m_ord
!
!   Task (ii) : numerical factorizaion of B
    CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
    m_1%JOB=2
    CALL dmumps(m_1)
    IF(m_1%infog(1).NE.0) THEN
        info=12
        CALL kep_info(info)
        RETURN
    END IF
    CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
    cmpt_time(2)=REAL(clock_2-clock_1)/REAL(clock_rate)
!
!   Task (iii) : Lanczos method with reverse communication==============
!
    ALLOCATE(alpha(imax_lanc),theta(imax_lanc),beta(0:imax_lanc),&
    &   v(n,imax_lanc+1),bv(n,0:imax_lanc+1))
!
!   Initialize
    CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
    CALL kep_lanczos(n,imax_lanc,its_lanc,ijob,v,bv,alpha,beta,theta)
    CALL kep_matvec(n,nz_b,row_pntr_b,col_indx_b,b,v(:,1),bv(:,1))
    CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
    cmpt_time(3)=cmpt_time(3)+REAL(clock_2-clock_1)/REAL(clock_rate)
!
    DO i=1,imax_lanc
!       Expand subspace and solve tridiagonal eigenproblem
        CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
        CALL kep_lanczos(n,imax_lanc,its_lanc,ijob,v,bv,alpha,beta,theta)
        CALL kep_matvec(n,nz_a,row_pntr_a,col_indx_a,a,v(:,its_lanc),bv(:,its_lanc+1))
        CALL kep_lanczos(n,imax_lanc,its_lanc,ijob,v,bv,alpha,beta,theta)
        m_1%rhs=bv(:,its_lanc+1)
        m_1%job=3
        CALL dmumps(m_1)
        IF(m_1%infog(1).NE.0) THEN
            info=13
            CALL kep_info(info)
            RETURN
        END IF
        v(:,its_lanc+1)=m_1%rhs
        CALL kep_lanczos(n,imax_lanc,its_lanc,ijob,v,bv,alpha,beta,theta)
        CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
        cmpt_time(3)=cmpt_time(3)+REAL(clock_2-clock_1)/REAL(clock_rate)
!
!       Task (iv) : inertia computation
        CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
        IF(si_1%inertia(i-1).GE.k) THEN
            si_1%shift(i)=theta(1)
        ELSE
            si_1%shift(i)=theta(i)
        END IF
        m_2%a(1:nz_a)=a
        m_2%a(nz_a+1:nz_a+nz_b)=-si_1%shift(i)*b
        m_2%job=4
        CALL dmumps(m_2)
        IF(m_2%infog(1).NE.0) THEN
            info=14
            CALL kep_info(info)
            RETURN
        END IF
        si_1%inertia(i)=m_2%infog(12)
        CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
        cmpt_time(4)=cmpt_time(4)+REAL(clock_2-clock_1)/REAL(clock_rate)
!
        icnt(1)=i
        IF(i.NE.1) THEN
            IF((si_1%inertia(i-1).GE.k).AND.(si_1%inertia(i).LT.k)) THEN
                EXIT
            ELSE IF((si_1%inertia(i-1).LT.k).AND.(si_1%inertia(i).GE.k)) THEN
                EXIT
            ELSE IF(i.EQ.imax_lanc) THEN
!               Failed to find an initial interval
                info=1
                DEALLOCATE(alpha,theta,beta,v,bv)
                DEALLOCATE(m_1%irn,m_1%jcn,m_1%a,m_1%rhs)
                DEALLOCATE(m_ord)
                DEALLOCATE(m_2%irn,m_2%jcn,m_2%a,m_2%perm_in)
                m_1%job=-2
                m_2%job=-2
                CALL dmumps(m_1)
                CALL dmumps(m_2)
                CALL mpi_finalize(ierr)
                CALL kep_info(info)
                IF(iprm(21).GE.0) CALL kep_summary(mmax,icnt,cmpt_time,&
                &   si_1,si_2_l,si_2_r,si_3,egnval,res,errnorm,info)
                DEALLOCATE(si_1%shift,si_1%inertia,                         &
                &   si_2_l%shift,si_2_l%inertia,si_2_r%shift,si_2_r%inertia,&
                &   si_3%shift,si_3%inertia,egnval,res,errnorm,egnvec,egnvec_old)
                RETURN
            END IF
        END IF
    END DO
!
!   End of Lanczos method===============================================
!
    DEALLOCATE(alpha,theta,beta,v,bv)
    DEALLOCATE(m_1%irn,m_1%jcn,m_1%a,m_1%rhs)
    m_1%job=-2
    CALL dmumps(m_1)
!
!---Step 2 : narrow down the interval-----------------------------------
!
!   Initialize
    si_2_l%shift(0)=MIN(si_1%shift(icnt(1)-1),si_1%shift(icnt(1)))
    si_2_r%shift(0)=MAX(si_1%shift(icnt(1)-1),si_1%shift(icnt(1)))
    si_2_l%inertia(0)=MIN(si_1%inertia(icnt(1)-1),si_1%inertia(icnt(1)))
    si_2_r%inertia(0)=MAX(si_1%inertia(icnt(1)-1),si_1%inertia(icnt(1)))
!
!   Bisection===========================================================
!
    DO i=1,imax_bi+1
 !
        icnt(2)=i-1
        icnt(3)=si_2_r%inertia(i-1)-si_2_l%inertia(i-1)
        IF(icnt(3).LE.mmax) THEN
            EXIT
        ELSE IF(i-1.EQ.imax_bi) THEN
!           Failed to narrow down an interval
            info=2
            DEALLOCATE(m_ord)
            DEALLOCATE(m_2%irn,m_2%jcn,m_2%a,m_2%perm_in)
            m_2%job=-2
            CALL dmumps(m_2)
            CALL mpi_finalize(ierr)
            CALL kep_info(info)
            IF(iprm(21).GE.0) CALL kep_summary(mmax,icnt,cmpt_time,&
            &   si_1,si_2_l,si_2_r,si_3,egnval,res,errnorm,info)
            DEALLOCATE(si_1%shift,si_1%inertia,                         &
            &   si_2_l%shift,si_2_l%inertia,si_2_r%shift,si_2_r%inertia,&
            &   si_3%shift,si_3%inertia,egnval,res,errnorm,egnvec,egnvec_old)
            RETURN
        END IF
!
!       Task (v) : bisection and inertia computation
        CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
        dtmp=(si_2_l%shift(i-1)+si_2_r%shift(i-1))*0.5D0
        m_2%a(1:nz_a)=a
        m_2%a(nz_a+1:nz_a+nz_b)=-dtmp*b
        m_2%job=4
        CALL dmumps(m_2)
        IF(m_2%infog(1).NE.0) THEN
            info=15
            CALL kep_info(info)
            RETURN
        END IF
!
        si_2_l%shift(i) = si_2_l%shift(i-1)
        si_2_r%shift(i) = si_2_r%shift(i-1)
        si_2_l%inertia(i) = si_2_l%inertia(i-1)
        si_2_r%inertia(i) = si_2_r%inertia(i-1)
!
        IF(m_2%infog(12).GE.k) THEN
            si_2_r%shift(i)=dtmp
            si_2_r%inertia(i)=m_2%infog(12)
        ELSE
            si_2_l%shift(i)=dtmp
            si_2_l%inertia(i)=m_2%infog(12)
        END IF
        CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
        cmpt_time(5)=cmpt_time(5)+REAL(clock_2-clock_1)/REAL(clock_rate)
    END DO
!
!   End of bisection====================================================
!
    DEALLOCATE(m_2%irn,m_2%jcn,m_2%a,m_2%perm_in)
    m_2%job=-2
    CALL dmumps(m_2)
!
!---Step 3 : compute eigenpairs of the interval-------------------------
!
!   Initialize
    lngth=(si_2_r%shift(icnt(2))-si_2_l%shift(icnt(2)))*0.5D0
    egnvec_old=0.0D0
!
!   MUMPS
    INCLUDE 'm_3'
!
!   Reuse fill-reducing ordering
    m_3%perm_in=m_ord
    DEALLOCATE(m_ord)
!
!   Task (vi) : numerical factorization of A-sigma*B
    CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
    dtmp=(si_2_l%shift(icnt(2))+si_2_r%shift(icnt(2)))*0.5D0
    m_3%a(1:nz_a)=a
    m_3%a(nz_a+1:nz_a+nz_b)=-dtmp*b
    m_3%JOB=4
    CALL dmumps(m_3)
    IF(m_3%infog(1).NE.0) THEN
        info=16
        CALL kep_info(info)
        RETURN
    END IF
    si_3%shift(1)=dtmp
    si_3%inertia(1)=m_3%infog(12)
    CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
    cmpt_time(6)=REAL(clock_2-clock_1)/REAL(clock_rate)
!
!   Task (vii) : shift-and-invert Lanczos with reverse communication====
!
    ALLOCATE(alpha_si(imax_si),theta_si(imax_si),beta_si(0:imax_si),&
    &   v_si(n,imax_si+1),bv_si(n,0:imax_si+1),y(imax_si,imax_si))
!
!   Initialize shift-and-invert Lanczos
    CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
    CALL kep_silanczos(n,imax_si,its_si,ijob_si,&
    &   v_si,bv_si,alpha_si,beta_si,theta_si,y,icnt(3))
    CALL kep_matvec(n,nz_b,row_pntr_b,col_indx_b,b,bv_si(:,1),v_si(:,1))
    CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
    cmpt_time(7)=cmpt_time(7)+REAL(clock_2-clock_1)/REAL(clock_rate)
!
    DO i=1,imax_si
!       Expand subspace and solve tridiagonal eigenproblem
        CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
        CALL kep_silanczos(n,imax_si,its_si,ijob_si,&
        &   v_si,bv_si,alpha_si,beta_si,theta_si,y,icnt(3))
        m_3%rhs=v_si(:,its_si)
        m_3%job=3
        CALL dmumps(m_3)
        IF(m_3%infog(1).NE.0) THEN
            info=17
            CALL kep_info(info)
            RETURN
        END IF
        bv_si(:,its_si+1)=m_3%rhs
        CALL kep_silanczos(n,imax_si,its_si,ijob_si,&
        &   v_si,bv_si,alpha_si,beta_si,theta_si,y,icnt(3))
        CALL kep_matvec(n,nz_b,row_pntr_b,col_indx_b,b,&
        &   bv_si(:,its_si+1),v_si(:,its_si+1))
        CALL kep_silanczos(n,imax_si,its_si,ijob_si,&
        &   v_si,bv_si,alpha_si,beta_si,theta_si,y,icnt(3))
!
!       Convergence test with reverse communication
        CALL kep_convtest(n,imax_si,mmax,tol_res,tol_dif,i,j,&
        &   ijob_test,v_si,bv_si,beta_si,theta_si,&
        &   y,icnt(3),lngth,si_3%shift(1),&
        &   egnval,egnvec,egnvec_old,res,errnorm,info_test)
        DO j=1,icnt(3)
            IF(info_test.EQ.0) THEN
                CALL kep_convtest(n,imax_si,mmax,tol_res,tol_dif,i,j,&
                &   ijob_test,v_si,bv_si,beta_si,theta_si,&
                &   y,icnt(3),lngth,si_3%shift(1),&
                &   egnval,egnvec,egnvec_old,res,errnorm,info_test)
                CALL kep_matvec(n,nz_a,row_pntr_a,col_indx_a,a,&
                &   egnvec(:,j),egnvec(:,mmax+1))
                CALL kep_matvec(n,nz_b,row_pntr_b,col_indx_b,b,&
                &   egnvec(:,j),egnvec(:,mmax+2))
                CALL kep_convtest(n,imax_si,mmax,tol_res,tol_dif,i,j,&
                &   ijob_test,v_si,bv_si,beta_si,theta_si,&
                &   y,icnt(3),lngth,si_3%shift(1),&
                &   egnval,egnvec,egnvec_old,res,errnorm,info_test)
            END IF
        END DO
        IF(info_test.EQ.0) THEN
            ijob_test=4
            CALL kep_convtest(n,imax_si,mmax,tol_res,tol_dif,i,j,&
            &   ijob_test,v_si,bv_si,beta_si,theta_si,&
            &   y,icnt(3),lngth,si_3%shift(1),&
            &   egnval,egnvec,egnvec_old,res,errnorm,info_test)
        END IF
        CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
        cmpt_time(7)=cmpt_time(7)+REAL(clock_2-clock_1)/REAL(clock_rate)
!
        icnt(4)=i
        IF(info_test.EQ.0) THEN
            EXIT
        ELSE IF(i.EQ.imax_si) THEN
!           Failed to be converged
            info=3
            DEALLOCATE(alpha_si,theta_si,beta_si,v_si,bv_si,y)
            DEALLOCATE(m_3%irn,m_3%jcn,m_3%a,m_3%perm_in,m_3%rhs)
            m_3%job=-2
            CALL dmumps(m_3)
            CALL mpi_finalize(ierr)
            CALL kep_info(info)
            IF(iprm(21).GE.0) CALL kep_summary(mmax,icnt,cmpt_time,&
            &   si_1,si_2_l,si_2_r,si_3,egnval,res,errnorm,info)
            DEALLOCATE(si_1%shift,si_1%inertia,                         &
            &   si_2_l%shift,si_2_l%inertia,si_2_r%shift,si_2_r%inertia,&
            &   si_3%shift,si_3%inertia,egnval,res,errnorm,egnvec,egnvec_old)
            RETURN
        ELSE IF(info_test.GE.3) THEN
            egnvec_old=egnvec
        END IF
        info_test=0
        ijob_test=1
    END DO
!
    keval=egnval(k-si_2_l%inertia(icnt(2)))
    kevec=egnvec(:,k-si_2_l%inertia(icnt(2)))
!
!   End of Shift-and-invert Lanczos method==============================
!
    DEALLOCATE(alpha_si,theta_si,beta_si,v_si,bv_si,y)
    DEALLOCATE(m_3%irn,m_3%jcn,m_3%a,m_3%perm_in,m_3%rhs)
    m_3%job=-2
    CALL dmumps(m_3)
    CALL mpi_finalize(ierr)
!
!---Time and results----------------------------------------------------
!
    IF(iprm(21).GE.0) CALL kep_summary(mmax,icnt,cmpt_time,&
        &   si_1,si_2_l,si_2_r,si_3,egnval,res,errnorm,info)
    DEALLOCATE(si_1%shift,si_1%inertia,                         &
    &   si_2_l%shift,si_2_l%inertia,si_2_r%shift,si_2_r%inertia,&
    &   si_3%shift,si_3%inertia,egnval,res,errnorm,egnvec,egnvec_old)
!
!-----------------------------------------------------------------------
!
    RETURN
    END SUBROUTINE kep
