    MODULE k2l_iotype
!
    TYPE :: k2l_io_type
        SEQUENCE
!       Input
        INTEGER :: job,n,nz_a,nz_b,k_lower,k_upper,info,ipad
        DOUBLE PRECISION :: s_lower,s_upper
        CHARACTER(16),DIMENSION(10) :: cprm
        INTEGER,DIMENSION(30) :: iprm
        INTEGER,ALLOCATABLE,DIMENSION(:) :: indx_a,jndx_a
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: rval_a
        INTEGER,ALLOCATABLE,DIMENSION(:) :: indx_b,jndx_b
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: rval_b        
!       Output        
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: kval
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: kvec
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: kipr
        INTEGER,ALLOCATABLE,DIMENSION(:) :: kndx
    END TYPE k2l_io_type
!
    END MODULE k2l_iotype
