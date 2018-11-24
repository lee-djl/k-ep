!-----------------------------------------------------------------------
!---Example program for demostrating how to use k2l---------------------
!-----------------------------------------------------------------------
    PROGRAM example
    USE example_utility     ! To read input matrices and write results
    USE k2l_iotype          ! To use k2l-related varialbes (derived type)
    USE k2l_mod             ! To use k2l
    IMPLICIT NONE
!---Temporary variables to read command line arguments
    CHARACTER(256), DIMENSION(10) :: cla
    INTEGER, INTRINSIC :: iargc
!
!---Input matrices A and B in coordinate format
    CHARACTER(256)  ::  &
    &   ipath_a,        &   ! file path for matrix A
    &   ipath_b             ! file path for matrix B
    INTEGER         ::  &
    &   n,              &   ! matrix size of A and B
    &   nz_a,           &   ! number of non-zero elements of A
    &   nz_b                ! number of non-zero elements of B
    INTEGER,ALLOCATABLE,DIMENSION(:) :: &
    &   indx_a,         &   ! row index of A
    &   jndx_a,         &   ! column index of A
    &   indx_b,         &   ! row index of B
    &   jndx_b              ! column index of B
    DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: &
    &   rval_a,         &   ! value of non-zero elements of A
    &   rval_b              ! value of non-zero elements of B
!
!---Input target index range [k_lower,k_upper]
    INTEGER ::              &
    &   targetindex_lower,  &   ! k_lower
    &   targetindex_upper       ! k_upper
!---Optional: User-defined initial interval [s_lower,s_upper) containing
!---the [k_lower,k_upper]-th eigenvalues
    DOUBLE PRECISION ::     &
    &   userinterval_lower, &   ! s_lower
    &   userinterval_upper      ! s_upper
!
!---MPI (for MUMPS)
    INTEGER :: ierr    
!
!---k2l
    TYPE(k2l_io_type) :: k2l_io ! k2l-related varialbes (derived type)
!-----------------------------------------------------------------------
!---Read command line arguments-----------------------------------------
!-----------------------------------------------------------------------
!---Provide the file path for A (coordinate format, Matrix Market MTX file)
!---as a command line argument
    CALL GETARG(1,cla(1))
    ipath_a=TRIM(ADJUSTL(cla(1)))
!---Provide the file path for B (coordinate format, Matrix Market MTX file)
!---as a command line argument
    CALL GETARG(2,cla(2))
    ipath_b=TRIM(ADJUSTL(cla(2)))
!---Provide the target index k_lower as a command line argument
    CALL GETARG(3,cla(3))
    READ(cla(3),*) targetindex_lower
!---Provide the target index k_upper as a command line argument
    CALL GETARG(4,cla(4)) 
    READ(cla(4),*) targetindex_upper
!
    IF(IARGC().GT.4) THEN
!-------Optional: Provide the lower endpoint of a user-defined initial
!-------interval as a command line argument
        CALL GETARG(5,cla(5))
        READ(cla(5),*) userinterval_lower
!-------Optional: Provide the upper endpoint of a user-defined initial
!-------interval as a command line argument
        CALL GETARG(6,cla(6))
        READ(cla(6),*) userinterval_upper
    END IF
!-----------------------------------------------------------------------
!---Read input matrices A and B from files------------------------------
!---(coordinate format, Matrix Market MTX file)-------------------------
!-----------------------------------------------------------------------
    CALL example_readmtx(ipath_a,n,nz_a,indx_a,jndx_a,rval_a)
    CALL example_readmtx(ipath_b,n,nz_b,indx_b,jndx_b,rval_b)
!-----------------------------------------------------------------------
!---Initialize k2l------------------------------------------------------
!-----------------------------------------------------------------------
    k2l_io%n=n              ! matrix size of A and B
    k2l_io%nz_a=nz_a        ! number of non-zero elements of A
    k2l_io%nz_b=nz_b        ! number of non-zero elements of B
!
    k2l_io%job=0            ! Initizalize k2l (job = 0)
    CALL k2l(k2l_io)
!-----------------------------------------------------------------------
!---Assign values of matrices A and B for k2l---------------------------
!-----------------------------------------------------------------------
    k2l_io%indx_a=indx_a    ! row index of A
    k2l_io%jndx_a=jndx_a    ! column index of A
    k2l_io%rval_a=rval_a    ! value of non-zero elements of A
!
    k2l_io%indx_b=indx_b    ! row index of B
    k2l_io%jndx_b=jndx_b    ! column index of B
    k2l_io%rval_b=rval_b    ! value of non-zero elements of B
!-----------------------------------------------------------------------
!---Set the target index range for k2l----------------------------------
!-----------------------------------------------------------------------
    k2l_io%k_lower=targetindex_lower
    k2l_io%k_upper=targetindex_upper
!-----------------------------------------------------------------------
!---Optional: Set input parameters for k2l------------------------------
!-----------------------------------------------------------------------
!---If cprm(1)='print', k2l prints details of its computation to terminal
    k2l_io%cprm(1)='print'
!
!---If cprm(2)='user', instead of setting an initial interval by itself,
!---k2l uses a user-defined initial interval [s_lower,s_upper) containing
!---the [k_lower,k_upper]-th eigenvalues
    ! k2l_io%cprm(2)='user'
    k2l_io%s_lower=userinterval_lower
    k2l_io%s_upper=userinterval_upper
!
!---If cprm(3)='ipr', k2l computes inverse participation ratio for
!---electronic structure calculations of materials
    k2l_io%cprm(3)='ipr'
!
!---If cprm(6)='second', k2l stops after the second stage (bisection) to compute
!---an interval [s_lower2,s_upper2) containing the [k_lower2,k_upper2]-th
!---eigenvalues (which include [k_lower,k_upper]-th eigenvalues).
    ! k2l_io%cprm(6)='second'
!
!---Stopping criterion for the second stage (bisection):
!---The initial interval is narrowed down until the number of eigenvalues
!---in the interval becomes smaller than equal to
!---MAX( k2l_io%iprm(10), CEILING( k2l_io%dprm(1) * (k_upper - k_lower + 1) ).
!---In particular, if k2l_io%iprm(10) = 1 and k2l_io%dprm(1) = 1.0,
!---the initial interval is narrowed down to contain only the
!---[k_lower,k_upper]-th eigenvalues.
!---It is recommended to change the following parameters only if users
!---know how the three-stage algorithm works.
    ! k2l_io%iprm(10)=1     ! Default value is 20
    ! k2l_io%dprm(1)=1.0D0  ! Default value is 1.0D0
!
!-----------------------------------------------------------------------
!---Compute the k_lower to k_upper eigenvalue pair(s)-------------------
!-----------------------------------------------------------------------
!---Initialize MPI (for using MUMPS)
    CALL mpi_init(ierr)
!        
    k2l_io%job=1            ! Solve an eigenproblem (job = 1)
    CALL k2l(k2l_io)
!    
!---Finalize MPI (for using MUMPS)
    CALL mpi_finalize(ierr)
!-----------------------------------------------------------------------
!---Write eigenpairs and inverse participation ratio (optional) to files
!-----------------------------------------------------------------------
    IF(k2l_io%info.EQ.0) THEN
        IF(TRIM(ADJUSTL(k2l_io%cprm(6))).NE.'second') THEN
!-----------Optional: Write the [k_lower,k_upper]-th eigenvalue(s) in one file
            CALL example_writekval(SIZE(k2l_io%kndx),k2l_io%kndx,k2l_io%kval)
!-----------Optional: Write the [k_lower,k_upper]-th eigenvector(s) in a separate file
            CALL example_writekvec(k2l_io%n,SIZE(k2l_io%kndx),k2l_io%kndx,k2l_io%kvec)
!
            IF(TRIM(ADJUSTL(k2l_io%cprm(3))).EQ.'ipr') THEN
!---------------Optional: Write the [k_lower,k_upper]-th inverse participation
!---------------ratios in one file        
                CALL example_writekipr(SIZE(k2l_io%kndx),k2l_io%kndx,k2l_io%kipr)
            END IF
        ELSE IF(TRIM(ADJUSTL(k2l_io%cprm(6))).EQ.'second') THEN
!-----------Optional: Write the interval [s_lower2,s_upper2) containing
!-----------the [k_lower2,k_upper2]-th eigenvalue(s) in one file
            CALL example_writekint(k2l_io%k_lower2,k2l_io%k_upper2,k2l_io%s_lower2,k2l_io%s_upper2)
        END IF
    END IF
!-----------------------------------------------------------------------
!---Finalize k2l to deallocate k2l-related memory-----------------------
!-----------------------------------------------------------------------
    k2l_io%job=-1            ! Finalize k2l (job = -1)
    CALL k2l(k2l_io)
!-----------------------------------------------------------------------
    DEALLOCATE(indx_a,jndx_a,rval_a,indx_b,jndx_b,rval_b)
    STOP
    END PROGRAM example
!-----------------------------------------------------------------------
!---End of the example program------------------------------------------
!-----------------------------------------------------------------------
