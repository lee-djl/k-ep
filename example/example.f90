    PROGRAM example
    USE mod_example_readmtx
    IMPLICIT NONE
!---Input: matrices A and B in Matrix Market exchange (MTX) format
!
    CHARACTER(256) :: ipath_a,ipath_b
    INTEGER :: iunit_a,iunit_b
!
    INTEGER :: n,nz_a,indx_a,jndx_a,nz_b,indx_b,jndx_b
    DOUBLE PRECISION :: rval_a,rval_b
    ALLOCATABLE :: indx_a,jndx_a,rval_a,indx_b,jndx_b,rval_b
    DIMENSION :: indx_a(:),jndx_a(:),rval_a(:),&
    &   indx_b(:),jndx_b(:),rval_b(:)
!
!---Input: parameters
!
    INTEGER :: k,iprm
    DIMENSION :: iprm(30)
!
!---Output: k-th eigenpair
!
    DOUBLE PRECISION :: keval,kevec
    ALLOCATABLE :: kevec
    DIMENSION :: kevec(:)
!
    INTEGER :: info
!
    CHARACTER(256) :: opath
    INTEGER :: ounit
!
!---Read matrices A and B (MTX format)----------------------------------
!
!   Read matrix A
    CALL getarg(1,ipath_a)
    iunit_a=999
    CALL example_readmtx(ipath_a,iunit_a,n,nz_a,indx_a,jndx_a,rval_a)
!
!   Read matrix B
    CALL getarg(2,ipath_b)
    iunit_b=998
    CALL example_readmtx(ipath_b,iunit_b,n,nz_b,indx_b,jndx_b,rval_b)
!
!---Set parameters------------------------------------------------------
!
    k=2343          ! target index
!
    iprm(1)=20      ! stopping criterion for bisection
                    ! (=maximum number of eigenpairs computed together)
    iprm(2)=10      ! 10**-iprm(2)=tolerance for relative residual 2-norm
    iprm(3)=10      ! 10**-iprm(3)=tolerance for relative difference 2-norm
    iprm(11)=30     ! maximum iteration count for Lanczos
    iprm(12)=30     ! maximum iteration count for Bisection
    iprm(13)=300    ! maximum iteration count for shift-and-invert Lanczos
!
    iprm(21)=1      ! suppress detailed output to terminal if iprm(21)<0
!
!---Compute the k-th eigenpair------------------------------------------
!
    ALLOCATE(kevec(n))
!
    CALL kep(n,&
    &   nz_a,indx_a,jndx_a,rval_a,&
    &   nz_b,indx_b,jndx_b,rval_b,&
    &   k,iprm,keval,kevec,info)
!
!---Write k-th eigenvector----------------------------------------------
!
    IF(info.EQ.0) THEN
        WRITE(opath,"(A6,I0.8,A4)") './evec',k,'.txt'
        ounit=899
        CALL example_writekevec(opath,ounit,n,kevec)
    END IF
!
!-----------------------------------------------------------------------
!
    DEALLOCATE(indx_a,jndx_a,rval_a,indx_b,jndx_b,rval_b)
    DEALLOCATE(kevec)
!
    STOP
    END PROGRAM example
