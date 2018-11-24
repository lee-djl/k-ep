    MODULE k2l_solver
    USE k2l_iotype
    USE k2l_inttype
    USE k2l_utility
    CONTAINS
!=======================================================================
    SUBROUTINE k2l_ldl_initialize(k2l_io,k2l_factor,shift,discardfactor)
    IMPLICIT NONE
    TYPE(k2l_io_type), INTENT(INOUT) :: k2l_io
    TYPE(k2l_factor_type), INTENT(OUT) :: k2l_factor
    DOUBLE PRECISION, INTENT(IN),OPTIONAL :: shift
    LOGICAL, INTENT(IN),OPTIONAL :: discardfactor
!
    SELECT CASE(TRIM(ADJUSTL(k2l_io%cprm(4))))
    CASE('dmumps')
        CALL k2l_ldl_initialize_dmumps(k2l_io,k2l_factor,shift,discardfactor)
    ! CASE('any_other_sparsesolver')
    !     CALL any_other_sparsesolver
    END SELECT
!
    RETURN
    END SUBROUTINE k2l_ldl_initialize
!=======================================================================
    SUBROUTINE k2l_ldl_initialize_dmumps(k2l_io,k2l_factor,shift,discardfactor)
    IMPLICIT NONE
    TYPE(k2l_io_type), INTENT(INOUT) :: k2l_io
    TYPE(k2l_factor_type), INTENT(OUT) :: k2l_factor
    DOUBLE PRECISION, INTENT(IN),OPTIONAL :: shift
    LOGICAL, INTENT(IN),OPTIONAL :: discardfactor
!
!   For dMUMPS
    INCLUDE 'mpif.h'
    INTEGER :: rank,ierr,i
!
!   Initialize dMUMPS
    k2l_factor%dmumps%comm=mpi_comm_world  ! MPI communicator: set by dMUMPS
    IF(PRESENT(shift)) THEN
        k2l_factor%dmumps%sym=2    ! matrix type: real symmetric indefinite
    ELSE
        k2l_factor%dmumps%sym=1    ! matrix type: real symmetric positive definite
    END IF
    k2l_factor%dmumps%par=1    ! parallelism for host processor: yes
    k2l_factor%dmumps%job=-1   ! job type: initialization
    CALL dmumps(k2l_factor%dmumps)
    IF(k2l_factor%dmumps%infog(1).LT.0) THEN
        CALL k2l_ldl_finalize_dmumps(k2l_factor)
        CALL k2l_dmumps_err(-1)
        k2l_io%info=15
        CALL k2l_info(k2l_io%info)
    END IF
!
!   Set parameters for dMUMPS
    IF(TRIM(ADJUSTL(k2l_io%cprm(1))).NE.'print') THEN
        k2l_factor%dmumps%icntl(3)=0    ! output stream for global information: suppressed
    END IF
    k2l_factor%dmumps%icntl(5)=0   ! input matrix format: assembled format,
                                   ! same as the Matrix Market format
    k2l_factor%dmumps%icntl(6)=0   ! no column permutation
    k2l_factor%dmumps%icntl(7)=5   ! fill-reducing ordering: METIS
    k2l_factor%dmumps%icntl(12)=1  ! usual ordering
    k2l_factor%dmumps%icntl(13)=1  ! parallelism for root node (MULTIFRONTAL node): NO,
                                   ! This parameter should NOT be changed.
    k2l_factor%dmumps%icntl(14)=20 ! extra working memory: 20% increase
    k2l_factor%dmumps%icntl(18)=0  ! distribution of input matrix: centralized on host processor
    k2l_factor%dmumps%icntl(22)=0  ! in-core factorization
    k2l_factor%dmumps%icntl(28)=1  ! parallelism for fill-reducing ordering: no,
                                   ! Currently using METIS (sequential),
                                   ! In case of using ParMETIS, set to icntl(28)=2
    k2l_factor%dmumps%icntl(29)=2  ! parallel fill-reducing ordering: ParMETIS
                                   ! This parameter is ignored when icntl(28)=1.
    IF(PRESENT(discardfactor).AND.(discardfactor.EQ..TRUE.)) THEN
        k2l_factor%dmumps%icntl(31)=1   ! discard factored matrix
                                        ! This parameter should NOT be changed.
    ELSE
        k2l_factor%dmumps%icntl(31)=0   ! keep factored matrix
                                        ! This parameter should NOT be changed.
    END IF
!
    CALL mpi_comm_rank(k2l_factor%dmumps%comm,rank,ierr)
    IF(rank.EQ.0) THEN
!   Allocate memory for dMUMPS & Load input matrices for dMUMPS
        IF(PRESENT(shift)) THEN
            k2l_factor%dmumps%n=k2l_io%n
            k2l_factor%dmumps%nz=k2l_io%nz_a
            IF(ALLOCATED(k2l_io%indx_b).AND.ALLOCATED(k2l_io%jndx_b).AND.ALLOCATED(k2l_io%rval_b)) THEN
                k2l_factor%dmumps%nz = k2l_factor%dmumps%nz + k2l_io%nz_b
            ELSE
                k2l_factor%dmumps%nz = k2l_factor%dmumps%nz + k2l_factor%dmumps%n
            END IF
            ALLOCATE(k2l_factor%dmumps%irn(k2l_factor%dmumps%nz),&
                &    k2l_factor%dmumps%jcn(k2l_factor%dmumps%nz))
            !
            k2l_factor%dmumps%irn(1:k2l_io%nz_a)=k2l_io%indx_a
            k2l_factor%dmumps%jcn(1:k2l_io%nz_a)=k2l_io%jndx_a
            IF(ALLOCATED(k2l_io%indx_b).AND.ALLOCATED(k2l_io%jndx_b).AND.ALLOCATED(k2l_io%rval_b)) THEN
                k2l_factor%dmumps%irn(k2l_io%nz_a+1:k2l_factor%dmumps%nz)=k2l_io%indx_b
                k2l_factor%dmumps%jcn(k2l_io%nz_a+1:k2l_factor%dmumps%nz)=k2l_io%jndx_b
            ELSE
                DO i=1,k2l_factor%dmumps%n,1
                    k2l_factor%dmumps%irn(k2l_io%nz_a+i)=i
                END DO
                DO i=1,k2l_factor%dmumps%n,1
                    k2l_factor%dmumps%jcn(k2l_io%nz_a+i)=i
                END DO
            END IF
        ELSE
            k2l_factor%dmumps%n=k2l_io%n
            k2l_factor%dmumps%nz=k2l_io%nz_b
            ALLOCATE(k2l_factor%dmumps%irn(k2l_factor%dmumps%nz),&
                &    k2l_factor%dmumps%jcn(k2l_factor%dmumps%nz))
            !
            k2l_factor%dmumps%irn(1:k2l_io%nz_a)=k2l_io%indx_b
            k2l_factor%dmumps%jcn(1:k2l_io%nz_a)=k2l_io%jndx_b
        END IF
    END IF
!
!   Compute fill-reducing ordering & Perform symbolic factorization
    k2l_factor%dmumps%job=1    ! job type: symbolic analysis
    CALL dmumps(k2l_factor%dmumps)
    IF(k2l_factor%dmumps%infog(1).LT.0) THEN
        CALL k2l_ldl_finalize_dmumps(k2l_factor)
        CALL k2l_dmumps_err(1)
        k2l_io%info=16
        CALL k2l_info(k2l_io%info)
    END IF
!
    IF(rank.EQ.0) ALLOCATE(k2l_factor%dmumps%a(k2l_factor%dmumps%nz))
!
    RETURN
    END SUBROUTINE k2l_ldl_initialize_dmumps
!=======================================================================
    SUBROUTINE k2l_ldl_factor(k2l_io,k2l_factor,nonzerofactor,shift,negativeinertia)
    IMPLICIT NONE
    TYPE(k2l_io_type), INTENT(INOUT) :: k2l_io
    TYPE(k2l_factor_type), INTENT(OUT) :: k2l_factor
    INTEGER(KIND=8), INTENT(OUT) :: nonzerofactor
    DOUBLE PRECISION, INTENT(IN),OPTIONAL :: shift
    INTEGER, INTENT(OUT),OPTIONAL :: negativeinertia
!
    SELECT CASE(TRIM(ADJUSTL(k2l_io%cprm(4))))
    CASE('dmumps')
        CALL k2l_ldl_factor_dmumps(k2l_io,k2l_factor,nonzerofactor,shift,negativeinertia)
    ! CASE('any_other_sparsesolver')
    !     CALL any_other_sparsesolver
    END SELECT
!
    RETURN
    END SUBROUTINE k2l_ldl_factor
!=======================================================================
    SUBROUTINE k2l_ldl_factor_dmumps(k2l_io,k2l_factor,nonzerofactor,shift,negativeinertia)
    IMPLICIT NONE
    TYPE(k2l_io_type), INTENT(INOUT) :: k2l_io
    TYPE(k2l_factor_type), INTENT(OUT) :: k2l_factor
    INTEGER(KIND=8), INTENT(OUT) :: nonzerofactor
    DOUBLE PRECISION, INTENT(IN),OPTIONAL :: shift
    INTEGER, INTENT(OUT),OPTIONAL :: negativeinertia
!
    INTEGER :: rank,ierr,i
!
!   Set parameters for dMUMPS
    k2l_factor%dmumps%icntl(31)=0  ! keep factored matrix
                                   ! This parameter should NOT be changed.
!
    CALL mpi_comm_rank(k2l_factor%dmumps%comm,rank,ierr)
    IF(rank.EQ.0) THEN
!       Load input matrices for dMUMPS
        IF(PRESENT(shift)) THEN
            k2l_factor%dmumps%a(1:k2l_io%nz_a)=k2l_io%rval_a
            IF(ALLOCATED(k2l_io%indx_b).AND.ALLOCATED(k2l_io%jndx_b).AND.ALLOCATED(k2l_io%rval_b)) THEN
                k2l_factor%dmumps%a(k2l_io%nz_a+1:k2l_factor%dmumps%nz)=-shift*k2l_io%rval_b
            ELSE
                k2l_factor%dmumps%a(k2l_io%nz_a+1:k2l_factor%dmumps%nz)=-shift
            END IF
        ELSE
            k2l_factor%dmumps%a(1:k2l_io%nz_a)=k2l_io%rval_b
        END IF
    END IF
!
!   Perform numerical factorization
    k2l_factor%dmumps%job=2    ! job type: numerical factorization
    CALL dmumps(k2l_factor%dmumps)
    IF(k2l_factor%dmumps%infog(1).LT.0) THEN
        CALL k2l_ldl_finalize_dmumps(k2l_factor)
        CALL k2l_dmumps_err(2)
        k2l_io%info=17
        CALL k2l_info(k2l_io%info)
    END IF
!
!   Get inertia
    IF(PRESENT(shift)) negativeinertia=k2l_factor%dmumps%infog(12)
!
!   Get the number of nonzero elements in the factor L of LDL factorization
    IF(k2l_factor%dmumps%infog(29).GE.0) THEN
        nonzerofactor=k2l_factor%dmumps%infog(29)
    ELSE
        nonzerofactor=INT(-k2l_factor%dmumps%infog(29) * (10**6),8)
    END IF
!
    RETURN
    END SUBROUTINE k2l_ldl_factor_dmumps
!=======================================================================
    SUBROUTINE k2l_ldl_inertia(k2l_io,k2l_factor,nonzerofactor,shift,negativeinertia)
    IMPLICIT NONE
    TYPE(k2l_io_type), INTENT(INOUT) :: k2l_io
    TYPE(k2l_factor_type), INTENT(OUT) :: k2l_factor
    INTEGER(KIND=8), INTENT(OUT) :: nonzerofactor
    DOUBLE PRECISION, INTENT(IN) :: shift
    INTEGER, INTENT(OUT) :: negativeinertia
!
    SELECT CASE(TRIM(ADJUSTL(k2l_io%cprm(4))))
    CASE('dmumps')
        CALL k2l_ldl_inertia_dmumps(k2l_io,k2l_factor,nonzerofactor,shift,negativeinertia)
    ! CASE('any_other_sparsesolver')
    !     CALL any_other_sparsesolver
    END SELECT
!
    RETURN
    END SUBROUTINE k2l_ldl_inertia
!=======================================================================
    SUBROUTINE k2l_ldl_inertia_dmumps(k2l_io,k2l_factor,nonzerofactor,shift,negativeinertia)
    IMPLICIT NONE
    TYPE(k2l_io_type), INTENT(INOUT) :: k2l_io
    TYPE(k2l_factor_type), INTENT(OUT) :: k2l_factor
    INTEGER(KIND=8), INTENT(OUT) :: nonzerofactor
    DOUBLE PRECISION, INTENT(IN) :: shift
    INTEGER, INTENT(OUT) :: negativeinertia
!
    INTEGER :: rank,ierr,i
!
!   Set parameters for dMUMPS
    k2l_factor%dmumps%icntl(31)=1  ! discard factored matrix
                                   ! This parameter should NOT be changed.
!
    CALL mpi_comm_rank(k2l_factor%dmumps%comm,rank,ierr)
    IF(rank.EQ.0) THEN
!       Load input matrices for dMUMPS
        k2l_factor%dmumps%a(1:k2l_io%nz_a)=k2l_io%rval_a
        IF(ALLOCATED(k2l_io%indx_b).AND.ALLOCATED(k2l_io%jndx_b).AND.ALLOCATED(k2l_io%rval_b)) THEN
            k2l_factor%dmumps%a(k2l_io%nz_a+1:k2l_factor%dmumps%nz)=-shift*k2l_io%rval_b
        ELSE
            k2l_factor%dmumps%a(k2l_io%nz_a+1:k2l_factor%dmumps%nz)=-shift
        END IF
    END IF
!
!   Perform numerical factorization
    k2l_factor%dmumps%job=2    ! job type: numerical factorization
    CALL dmumps(k2l_factor%dmumps)
    IF(k2l_factor%dmumps%infog(1).LT.0) THEN
        CALL k2l_ldl_finalize_dmumps(k2l_factor)
        CALL k2l_dmumps_err(2)
        k2l_io%info=17
        CALL k2l_info(k2l_io%info)
    END IF
!
!   Get inertia
    negativeinertia=k2l_factor%dmumps%infog(12)
!
!   Get the number of nonzero elements in the factor L of LDL factorization
    IF(k2l_factor%dmumps%infog(29).GE.0) THEN
        nonzerofactor=k2l_factor%dmumps%infog(29)
    ELSE
        nonzerofactor=INT(-k2l_factor%dmumps%infog(29) * (10**6),8)
    END IF
!
    RETURN
    END SUBROUTINE k2l_ldl_inertia_dmumps
!=======================================================================
    SUBROUTINE k2l_ldl_linsol(k2l_io,k2l_factor,x,y)
    IMPLICIT NONE
    TYPE(k2l_io_type), INTENT(INOUT) :: k2l_io
    TYPE(k2l_factor_type), INTENT(INOUT) :: k2l_factor
    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: x,y
    !
    SELECT CASE(TRIM(ADJUSTL(k2l_io%cprm(4))))
    CASE('dmumps')
        CALL k2l_ldl_linsol_dmumps(k2l_io,k2l_factor,x,y)
    ! CASE('any_other_sparsesolver')
    !     CALL any_other_sparsesolver
    END SELECT
    !
    RETURN
    END SUBROUTINE k2l_ldl_linsol
!=======================================================================
    SUBROUTINE k2l_ldl_linsol_dmumps(k2l_io,k2l_factor,x,y)
    IMPLICIT NONE
    TYPE(k2l_io_type), INTENT(INOUT) :: k2l_io
    TYPE(k2l_factor_type), INTENT(INOUT) :: k2l_factor
    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: x,y
!
    INTEGER :: rank,ierr,i
!
    IF((k2l_io%n.NE.SIZE(x)).OR.(k2l_io%n.NE.SIZE(y))) THEN
        WRITE(*,*) 'Incompatible array size'
        STOP
    END IF
!
!   Set parameters for dMUMPS
    k2l_factor%dmumps%icntl(10)=1   ! iterative refinement
!
    CALL mpi_comm_rank(k2l_factor%dmumps%comm,rank,ierr)
    IF(rank.EQ.0) THEN
!       Allocate memory for dMUMPS & Load input right-hand-side for dMUMPS
        ALLOCATE(k2l_factor%dmumps%rhs(k2l_io%n))
        k2l_factor%dmumps%rhs=x
    END IF
!
!   Compute fill-reducing ordering & Perform symbolic factorization
    k2l_factor%dmumps%job=3    ! job type: solution of linear systems
    CALL dmumps(k2l_factor%dmumps)
    IF(k2l_factor%dmumps%infog(1).LT.0) THEN
        CALL k2l_ldl_finalize_dmumps(k2l_factor)
        CALL k2l_dmumps_err(3)
        k2l_io%info=18
        CALL k2l_info(k2l_io%info)
    END IF
!
    IF(rank.EQ.0) y=k2l_factor%dmumps%rhs
!
    RETURN
    END SUBROUTINE k2l_ldl_linsol_dmumps
!=======================================================================
    SUBROUTINE k2l_ldl_finalize(k2l_io,k2l_factor)
    IMPLICIT NONE
    TYPE(k2l_io_type), INTENT(IN) :: k2l_io
    TYPE(k2l_factor_type), INTENT(INOUT) :: k2l_factor
!
    SELECT CASE(TRIM(ADJUSTL(k2l_io%cprm(4))))
    CASE('dmumps')
        CALL k2l_ldl_finalize_dmumps(k2l_factor)
    ! CASE('any_other_sparsesolver')
    !     CALL any_other_sparsesolver
    END SELECT
!
    RETURN
    END SUBROUTINE k2l_ldl_finalize
!=======================================================================
    SUBROUTINE k2l_ldl_finalize_dmumps(k2l_factor)
    IMPLICIT NONE
    TYPE(k2l_factor_type), INTENT(INOUT) :: k2l_factor
!
    INTEGER :: rank,ierr
!
    CALL mpi_comm_rank(k2l_factor%dmumps%comm,rank,ierr)
    ! IF(rank.EQ.0) THEN
!       Deallocate memory for dMUMPS
!       Note: intrinsic subroutine ALLOCATED fails to work when using gfortran compiler
!       Note: intrinsic subroutine ALLOCATED does work when using Intel compiler (ifort)
        ! IF(ALLOCATED(k2l_factor%dmumps%irn)) DEALLOCATE(k2l_factor%dmumps%irn)
        ! IF(ALLOCATED(k2l_factor%dmumps%jcn)) DEALLOCATE(k2l_factor%dmumps%jcn)
        ! IF(ALLOCATED(k2l_factor%dmumps%a)) DEALLOCATE(k2l_factor%dmumps%a)
        ! IF(ALLOCATED(k2l_factor%dmumps%rhs)) DEALLOCATE(k2l_factor%dmumps%rhs)
    ! END IF
!
    IF(rank.EQ.0) THEN
!       Deallocate memory for dMUMPS
        IF(k2l_factor%dmumps%job.NE.-1) THEN
            DEALLOCATE(k2l_factor%dmumps%irn)
            DEALLOCATE(k2l_factor%dmumps%jcn)
            DEALLOCATE(k2l_factor%dmumps%a)
        END IF
        IF(k2l_factor%dmumps%job.EQ.3) THEN
            DEALLOCATE(k2l_factor%dmumps%rhs)
        END IF
    END IF        
!
!   Finallize dMUMPS
    k2l_factor%dmumps%job=-2    ! job type: finalization
    CALL dmumps(k2l_factor%dmumps)
!
    RETURN
    END SUBROUTINE  k2l_ldl_finalize_dmumps
!=======================================================================
    SUBROUTINE k2l_dmumps_err(info)
    IMPLICIT NONE
    INTEGER :: info
!
!   For dMUMPS
    INCLUDE 'mpif.h'
    INTEGER :: rank,ierr
!
    CALL mpi_comm_rank(mpi_comm_world,rank,ierr)
    IF(rank.EQ.0) THEN
        SELECT CASE(info)
            CASE(-1)
                WRITE(6,*) '=====k2l: dMUMPS failed to be initialized'
            CASE(1)
                WRITE(6,*) '=====k2l: dMUMPS failed to perform symbolic factorization'
            CASE(2)
                WRITE(6,*) '=====k2l: dMUMPS failed to perform numerical factorization'
            CASE(3)
                WRITE(6,*) '=====k2l: dMUMPS failed to solve linear system'
            END SELECT
    END IF
!
    RETURN
    END SUBROUTINE k2l_dmumps_err
!=======================================================================
    SUBROUTINE k2l_matvec(selectmatrix,k2l_int,x,y)
    IMPLICIT NONE
    CHARACTER(*) :: selectmatrix
    TYPE(k2l_int_type), INTENT(IN) :: k2l_int
    DOUBLE PRECISION,DIMENSION(:) :: x,y
!-----------------------------------------------------------------------
    SELECT CASE(selectmatrix)
    CASE('A')
        CALL k2l_matvec2(SIZE(k2l_int%row_pntr_a)-1,SIZE(k2l_int%col_indx_a),&
            &   k2l_int%row_pntr_a,k2l_int%col_indx_a,k2l_int%a,x,y)
    CASE('B')
        CALL k2l_matvec2(SIZE(k2l_int%row_pntr_b)-1,SIZE(k2l_int%col_indx_b),&
            &   k2l_int%row_pntr_b,k2l_int%col_indx_b,k2l_int%b,x,y)
    CASE DEFAULT
        y=x
    END SELECT
!
    RETURN
    END SUBROUTINE k2l_matvec
!=======================================================================
    SUBROUTINE k2l_matvec2(n,nz,row_pntr,col_indx,a,x,y)
    IMPLICIT NONE
    INTEGER :: n,nz
    INTEGER, DIMENSION(:) :: row_pntr,col_indx
    DOUBLE PRECISION, DIMENSION(:) :: a,x,y
!
    INTEGER :: i,j,k
!-----------------------------------------------------------------------
    IF((n+1.NE.SIZE(row_pntr)).OR.(nz.NE.SIZE(col_indx)).OR.(nz.NE.SIZE(a))   &
    &   .OR.(n.NE.SIZE(x)).OR.(n.NE.SIZE(y))) THEN
        WRITE(*,*) 'Incompatible array size'
        STOP
    END IF

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
    END SUBROUTINE k2l_matvec2
!=======================================================================
    SUBROUTINE k2l_bracketing(k2l_io,n,lower,upper,iter,innerpoint)
    IMPLICIT NONE
    INTEGER :: n,iter
    TYPE(k2l_io_type) :: k2l_io
    DOUBLE PRECISION, DIMENSION(0:), INTENT(IN) :: lower,upper
    DOUBLE PRECISION, INTENT(out) :: innerpoint
    !
    SELECT CASE(TRIM(ADJUSTL(k2l_io%cprm(5))))
    CASE('bisection')
        CALL k2l_bracketing_bisection(lower(iter-1),upper(iter-1),innerpoint)
    ! CASE('any_other_bracketingalgorithm')
    !     CALL any_other_bracketingalgorithm
    END SELECT
    !
    RETURN
    END SUBROUTINE k2l_bracketing
!=======================================================================
    SUBROUTINE k2l_bracketing_bisection(lower,upper,innerpoint)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: lower,upper
    DOUBLE PRECISION, INTENT(INOUT) :: innerpoint
    !
    innerpoint=(lower+upper)*0.5D0
    !
    RETURN
    END SUBROUTINE k2l_bracketing_bisection
!=======================================================================
    SUBROUTINE k2l_innpro(n,x,y,alpha)
    IMPLICIT NONE
    INTEGER :: n
    DOUBLE PRECISION :: alpha
    DOUBLE PRECISION, DIMENSION(:) :: x,y
!
    INTEGER :: i
!-----------------------------------------------------------------------
    IF((n.NE.SIZE(x)).OR.(n.NE.SIZE(y))) THEN
        WRITE(*,*) 'Incompatible array size'
        STOP
    END IF

    alpha=0.0D0
!
    DO i=1,n
        alpha=alpha+x(i)*y(i)
    END DO
!
    RETURN
    END SUBROUTINE k2l_innpro
!=======================================================================
    SUBROUTINE k2l_ipr(k2l_io,k2l_int)
    IMPLICIT NONE
    TYPE(k2l_io_type) :: k2l_io
    TYPE(k2l_int_type) :: k2l_int
!
    INTEGER :: i
    DOUBLE PRECISION :: dtmp
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: vtmp
!-----------------------------------------------------------------------
!
    IF(TRIM(ADJUSTL(k2l_io%cprm(6))).EQ.'second') THEN
        RETURN
    END IF
!
    IF(TRIM(ADJUSTL(k2l_io%cprm(3))).NE.'ipr') THEN
        RETURN
    END IF
!
    ALLOCATE(vtmp(k2l_io%n))
    IF(ALLOCATED(k2l_io%kipr)) DEALLOCATE(k2l_io%kipr)
    ALLOCATE(k2l_io%kipr(k2l_io%k_upper-k2l_io%k_lower+1))
!
!   Compute Inverse Participation Ratio
    DO i=1,k2l_io%k_upper-k2l_io%k_lower+1
        CALL k2l_fournormquad(k2l_io%n,k2l_io%kvec(:,i),k2l_io%kipr(i))
        CALL k2l_matvec('B',k2l_int,k2l_io%kvec(:,i),vtmp)
        CALL k2l_innpro(k2l_io%n,k2l_io%kvec(:,i),vtmp,dtmp)
        dtmp=dtmp*dtmp
        k2l_io%kipr(i)=k2l_io%kipr(i)/dtmp
    END DO
    DEALLOCATE(vtmp)
!
    RETURN
    END SUBROUTINE k2l_ipr
!=======================================================================
    SUBROUTINE k2l_fournormquad(n,x,alpha)
    IMPLICIT NONE
    INTEGER :: n
    DOUBLE PRECISION :: alpha
    DOUBLE PRECISION, DIMENSION(:) :: x
!
    INTEGER :: i
    DOUBLE PRECISION :: temp
!-----------------------------------------------------------------------
    IF(n.NE.SIZE(x)) THEN
        WRITE(*,*) 'Incompatible array size'
        STOP
    END IF

    alpha=0.0D0
!
    DO i=1,n
        temp=x(i)*x(i)
        alpha=alpha+temp*temp
    END DO
!
    RETURN
    END SUBROUTINE k2l_fournormquad
!=======================================================================
    END MODULE k2l_solver
