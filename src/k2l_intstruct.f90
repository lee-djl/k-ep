    MODULE k2l_inttype
!
    TYPE :: k2l_int_type
        SEQUENCE
!       Computation time
        REAL,DIMENSION(10) :: cmpt_time
!       Iteration count for internal job
        INTEGER(KIND=8),DIMENSION(10) :: icnt
!       Input matrix in CRS format        
        INTEGER,ALLOCATABLE,DIMENSION(:) :: row_pntr_a,row_pntr_b,  &
            &   col_indx_a,col_indx_b
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: a,b
!       Shift and inertia for eigenvalue interval
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: &
            &   shift_1_lower,shift_1_upper, &
            &   shift_2_lower,shift_2_upper, &
            &   shift_2_lower2,shift_2_upper2, &
            &   shift_3
        INTEGER,ALLOCATABLE,DIMENSION(:) :: &
            &   inertia_1_lower,inertia_1_upper, &
            &   inertia_2_lower,inertia_2_upper, &
            &   inertia_2_lower2,inertia_2_upper2, &
            &   inertia_3
    END TYPE k2l_int_type
!
    INCLUDE 'dmumps_struc.h'
!   INCLUDE 'any_other_sparsesolver_struc.h'
    TYPE :: k2l_factor_type
        TYPE(dmumps_struc) :: dmumps
!       TYPE(any_other_sparsesolver_struc) :: any_other_sparsesolver
    END TYPE k2l_factor_type
!    
    TYPE :: k2l_lanczos_type
        SEQUENCE
        INTEGER :: its,ijob
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: alpha,theta,beta
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: v,bv
    END TYPE k2l_lanczos_type
!
    TYPE :: k2l_silanczos_type
        SEQUENCE
        INTEGER :: its,ijob,icnt,ipad
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: alpha,theta,beta
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: v,bv,y
    END TYPE k2l_silanczos_type
!
    TYPE :: k2l_convtest_type
        SEQUENCE
        INTEGER :: ijob,info,mmax,ipad        
        DOUBLE PRECISION :: int_length,tol_res,tol_dif
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: egnval,resnorm,difnorm
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: egnvec,egnvec_old
    END TYPE k2l_convtest_type
!        
    END MODULE k2l_inttype
