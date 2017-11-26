    PROGRAM example
    USE mod_example_readmtx
    USE mod_kep
    IMPLICIT NONE
!    
!   Input matrices A and B in Matrix Market exchange (MTX) format
    CHARACTER(256)  :: ipath_a,ipath_b
    INTEGER         :: n,nz_a,nz_b
    INTEGER,ALLOCATABLE,DIMENSION(:) :: indx_a,jndx_a,indx_b,jndx_b
    DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: rval_a,rval_b
!
!   For SUBROUTINE kep
    INCLUDE 'kep_struct.h'
    TYPE(kep_struct) :: kep_io
!
!---Read matrices A and B (MTX format)----------------------------------
!
!   Read A from file
    CALL getarg(1,ipath_a)
    CALL example_readmtx(ipath_a,n,nz_a,indx_a,jndx_a,rval_a)
!
!   Read B from file
    CALL getarg(2,ipath_b)
    CALL example_readmtx(ipath_b,n,nz_b,indx_b,jndx_b,rval_b)
!    
!===Assign values of A and B for kep====================================    
!    
!   Assign values of A
    kep_io%n=n              ! matrix size
    kep_io%nz_a=nz_a        ! number of non-zero elements of A
    ALLOCATE(kep_io%indx_a(kep_io%nz_a),&
        &    kep_io%jndx_a(kep_io%nz_a),&
        &    kep_io%rval_a(kep_io%nz_a))
    kep_io%indx_a=indx_a    ! row index of A
    kep_io%jndx_a=jndx_a    ! column index of A
    kep_io%rval_a=rval_a    ! value of A
!
!   Assign values of B
    kep_io%nz_b=nz_b        ! number of non-zero elements of B
    ALLOCATE(kep_io%indx_b(kep_io%nz_b),&
        &    kep_io%jndx_b(kep_io%nz_b),&
        &    kep_io%rval_b(kep_io%nz_b))
    kep_io%indx_b=indx_b    ! row index of B
    kep_io%jndx_b=jndx_b    ! columns index of B
    kep_io%rval_b=rval_b    ! value of B
!
!---Set parameters for kep----------------------------------------------
!
    kep_io%k       =2343! target index
!
    kep_io%iprm(1) =20  ! stopping criterion for bisection
                        ! (=maximum number of eigenpairs computed together)
    kep_io%iprm(2) =10  ! 10**-iprm(2)=tolerance for relative residual 2-norm
    kep_io%iprm(3) =10  ! 10**-iprm(3)=tolerance for relative difference 2-norm
    kep_io%iprm(11)=30  ! maximum iteration count for Lanczos
    kep_io%iprm(12)=30  ! maximum iteration count for Bisection
    kep_io%iprm(13)=300 ! maximum iteration count for shift-and-invert Lanczos
!
    kep_io%iprm(21)=1   ! If iprm(21) > 0, print details of computation to terminal
    kep_io%iprm(22)=0   ! If iprm(22) =< 0, output only the k-th eigenpair.
                        ! Otherwise, output all computed eigenpairs    
!
!---Compute the k-th eigenpair------------------------------------------
!
    CALL kep(kep_io)
!
!---Write k-th eigenvalue-----------------------------------------------
!
    IF(kep_io%info.EQ.0) THEN
        CALL example_writekval(SIZE(kep_io%kndx),kep_io%kndx,kep_io%kval)
    END IF
!
!---Write k-th eigenvector----------------------------------------------
!
    IF(kep_io%info.EQ.0) THEN
        CALL example_writekvec(&
            & kep_io%n,SIZE(kep_io%kndx),kep_io%kndx,kep_io%kvec)
    END IF
!
!=======================================================================
!
    DEALLOCATE(indx_a,jndx_a,rval_a,indx_b,jndx_b,rval_b)    
    DEALLOCATE(kep_io%indx_a,kep_io%jndx_a,kep_io%rval_a)
    DEALLOCATE(kep_io%indx_b,kep_io%jndx_b,kep_io%rval_b)
!
    STOP
    END PROGRAM example
