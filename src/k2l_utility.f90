    MODULE k2l_utility
    USE k2l_iotype
    USE k2l_inttype
    CONTAINS
!=======================================================================
    SUBROUTINE k2l_initialize_external(k2l_io)
    IMPLICIT NONE
    TYPE(k2l_io_type) :: k2l_io
!-----------------------------------------------------------------------
!   Initialize parameters
    k2l_io%k_lower=0            ! target index (lower)
    k2l_io%k_upper=0            ! target index (upper)
    k2l_io%s_lower=0.0D0        ! lower endpoint of user-defined interval
    k2l_io%s_upper=0.0D0        ! upper endpoint of user-defined interval
!
!   Assign default values to parameters
    k2l_io%cprm(1)=''           ! If cprm(1) = 'print', print details of computation to terminal
    k2l_io%cprm(2)=''           ! If cprm(2) = 'user', use a user-defined initial interval
                                ! [s_lower,s_upper) containing the [k_lower to k_upper]-th eigenvalues
    k2l_io%cprm(3)=''           ! If cprm(3) = 'ipr', compute inverse participation ratio            
    k2l_io%cprm(4)='dmumps'     ! Sparse direct linear solver (Currently supports only double-precision MUMPS)
    k2l_io%cprm(5)='bisection'  ! Algorithm to narrowing down the initial internal (Currently only bisection is implemented) 
!
    k2l_io%iprm(1) =8   ! 10**-iprm(1)=tolerance for relative residual 2-norm
    k2l_io%iprm(2) =6   ! 10**-iprm(2)=tolerance for relative difference 2-norm
    k2l_io%iprm(10)=20  ! stopping criterion for bisection    
    k2l_io%iprm(11)=MIN(50,k2l_io%n)    ! maximum iteration count for Lanczos
    k2l_io%iprm(12)=50                  ! maximum iteration count for Bisection
    k2l_io%iprm(13)=MIN(1000,k2l_io%n)  ! maximum iteration count for shift-and-invert Lanczos
!
!   Allocate    
    ALLOCATE(k2l_io%indx_a(k2l_io%nz_a),k2l_io%jndx_a(k2l_io%nz_a),k2l_io%rval_a(k2l_io%nz_a))    
    ALLOCATE(k2l_io%indx_b(k2l_io%nz_b),k2l_io%jndx_b(k2l_io%nz_b),k2l_io%rval_b(k2l_io%nz_b))
!-----------------------------------------------------------------------    
    RETURN
    END SUBROUTINE k2l_initialize_external
!=======================================================================
    SUBROUTINE k2l_initialize_internal(k2l_io,k2l_int)
    IMPLICIT NONE
    TYPE(k2l_io_type) :: k2l_io
    TYPE(k2l_int_type) :: k2l_int
!-----------------------------------------------------------------------
!   Check input parameters
    CALL k2l_checkinput(k2l_io)
!    
!   Initialize internal variables
    k2l_int%icnt=0
    k2l_int%cmpt_time=0.0E0
!    
!   Convert coordinate format to Compressed Row Storage (CRS) format       
    CALL k2l_crd2crs(k2l_io,k2l_int)
!-----------------------------------------------------------------------    
    RETURN
    END SUBROUTINE k2l_initialize_internal    
!=======================================================================
    SUBROUTINE k2l_checkinput(k2l_io)
    IMPLICIT NONE
    TYPE(k2l_io_type) :: k2l_io
!-----------------------------------------------------------------------
    k2l_io%info=0
!    
    IF(.NOT.(ALLOCATED(k2l_io%indx_a).AND.ALLOCATED(k2l_io%jndx_a).AND.ALLOCATED(k2l_io%rval_a))   ) THEN
        k2l_io%info=-1    
    ELSE IF((k2l_io%k_lower.LE.1).OR.(k2l_io%k_lower.GE.k2l_io%n).OR.&
        &   (k2l_io%k_upper.LE.1).OR.(k2l_io%k_upper.GE.k2l_io%n).OR.&
        &   (k2l_io%k_lower.GT.k2l_io%k_upper)) THEN
        k2l_io%info=-11
    ELSE IF((k2l_io%iprm(10).LT.1).OR.(k2l_io%iprm(10).GE.k2l_io%n)) THEN
        k2l_io%info=-12
    ELSE IF((k2l_io%iprm(1).LT.1).OR.(k2l_io%iprm(1).GT.16)) THEN
        k2l_io%info=-13
    ELSE IF((k2l_io%iprm(2).LT.1).OR.(k2l_io%iprm(2).GT.16)) THEN
        k2l_io%info=-14
    ELSE IF((TRIM(ADJUSTL(k2l_io%cprm(2))).NE.'user').AND.((k2l_io%iprm(11).LT.2).OR.(k2l_io%iprm(11).GT.k2l_io%n))) THEN
        k2l_io%info=-15
    ELSE IF((k2l_io%iprm(12).LT.0).OR.(k2l_io%iprm(12).GT.64)) THEN
        k2l_io%info=-16
    ELSE IF((k2l_io%iprm(13).LT.k2l_io%iprm(10)).OR.(k2l_io%iprm(13).GT.k2l_io%n)) THEN
        k2l_io%info=-17
    ELSE IF((TRIM(ADJUSTL(k2l_io%cprm(2))).EQ.'user').AND.(k2l_io%s_lower.GE.k2l_io%s_upper)) THEN
        k2l_io%info=-18
    ELSE IF(TRIM(ADJUSTL(k2l_io%cprm(4))).NE.'dmumps') THEN
        k2l_io%info=-21
    ELSE IF(TRIM(ADJUSTL(k2l_io%cprm(5))).NE.'bisection') THEN
        k2l_io%info=-22
    END IF
!    
    CALL k2l_info(k2l_io%info)
!-----------------------------------------------------------------------
    RETURN
    END SUBROUTINE k2l_checkinput
!=======================================================================
    SUBROUTINE k2l_finalize_internal(k2l_int)
    IMPLICIT NONE
    TYPE(k2l_int_type) :: k2l_int
!-----------------------------------------------------------------------
    IF(ALLOCATED(k2l_int%row_pntr_a)) DEALLOCATE(k2l_int%row_pntr_a)
    IF(ALLOCATED(k2l_int%row_pntr_b)) DEALLOCATE(k2l_int%row_pntr_b)
    IF(ALLOCATED(k2l_int%col_indx_a)) DEALLOCATE(k2l_int%col_indx_a)
    IF(ALLOCATED(k2l_int%col_indx_b)) DEALLOCATE(k2l_int%col_indx_b)
    IF(ALLOCATED(k2l_int%a)) DEALLOCATE(k2l_int%a)
    IF(ALLOCATED(k2l_int%b)) DEALLOCATE(k2l_int%b)
!    
    IF(ALLOCATED(k2l_int%shift_1_lower)) DEALLOCATE(k2l_int%shift_1_lower)
    IF(ALLOCATED(k2l_int%shift_1_upper)) DEALLOCATE(k2l_int%shift_1_upper)
    IF(ALLOCATED(k2l_int%shift_2_lower)) DEALLOCATE(k2l_int%shift_2_lower)
    IF(ALLOCATED(k2l_int%shift_2_upper)) DEALLOCATE(k2l_int%shift_2_upper)
    IF(ALLOCATED(k2l_int%shift_2_lower2)) DEALLOCATE(k2l_int%shift_2_lower2)
    IF(ALLOCATED(k2l_int%shift_2_upper2)) DEALLOCATE(k2l_int%shift_2_upper2)
    IF(ALLOCATED(k2l_int%shift_3)) DEALLOCATE(k2l_int%shift_3)
!    
    IF(ALLOCATED(k2l_int%inertia_1_lower)) DEALLOCATE(k2l_int%inertia_1_lower)
    IF(ALLOCATED(k2l_int%inertia_1_upper)) DEALLOCATE(k2l_int%inertia_1_upper)
    IF(ALLOCATED(k2l_int%inertia_2_lower)) DEALLOCATE(k2l_int%inertia_2_lower)
    IF(ALLOCATED(k2l_int%inertia_2_upper)) DEALLOCATE(k2l_int%inertia_2_upper)
    IF(ALLOCATED(k2l_int%inertia_2_lower2)) DEALLOCATE(k2l_int%inertia_2_lower2)
    IF(ALLOCATED(k2l_int%inertia_2_upper2)) DEALLOCATE(k2l_int%inertia_2_upper2)
    IF(ALLOCATED(k2l_int%inertia_3)) DEALLOCATE(k2l_int%inertia_3)
!-----------------------------------------------------------------------    
    RETURN
    END SUBROUTINE k2l_finalize_internal
!=======================================================================
    SUBROUTINE k2l_finalize_external(k2l_io)
    IMPLICIT NONE
    TYPE(k2l_io_type) :: k2l_io
!-----------------------------------------------------------------------
    IF(ALLOCATED(k2l_io%indx_a)) DEALLOCATE(k2l_io%indx_a)
    IF(ALLOCATED(k2l_io%jndx_a)) DEALLOCATE(k2l_io%jndx_a)
    IF(ALLOCATED(k2l_io%rval_a)) DEALLOCATE(k2l_io%rval_a)
!
    IF(ALLOCATED(k2l_io%indx_b)) DEALLOCATE(k2l_io%indx_b)
    IF(ALLOCATED(k2l_io%jndx_b)) DEALLOCATE(k2l_io%jndx_b)
    IF(ALLOCATED(k2l_io%rval_b)) DEALLOCATE(k2l_io%rval_b)
!
    IF(ALLOCATED(k2l_io%kval)) DEALLOCATE(k2l_io%kval)
    IF(ALLOCATED(k2l_io%kvec)) DEALLOCATE(k2l_io%kvec)
    IF(ALLOCATED(k2l_io%kndx)) DEALLOCATE(k2l_io%kndx)
!
    IF(ALLOCATED(k2l_io%kipr)) DEALLOCATE(k2l_io%kipr)
!-----------------------------------------------------------------------    
    RETURN
    END SUBROUTINE k2l_finalize_external   
!=======================================================================
    SUBROUTINE k2l_info(info)
    IMPLICIT NONE
    INTEGER :: info
    INTEGER :: ierr        
!-----------------------------------------------------------------------
!
    IF(info.NE.0) THEN
        WRITE(6,*)
    END IF
!    
    IF(info.EQ.-1) THEN
        WRITE(6,*)'=====k2l info',info,': Allocate input matrices and assign values'    
    ELSE IF(info.EQ.-11) THEN
        WRITE(6,*)'=====k2l info',info,&
            &   ': k_lower and k_upper must satisfy 1 < k_lower <= k_upper < n'
    ELSE IF(info.EQ.-12) THEN
        WRITE(6,*)'=====k2l info',info,': iprm(10) must satisfy 1 <= iprm(10) < n'
    ELSE IF(info.EQ.-13) THEN
        WRITE(6,*)'=====k2l info',info,': iprm(1) must satisfy 1 <= iprm(1) <= 16'
    ELSE IF(info.EQ.-14) THEN
        WRITE(6,*)'=====k2l info',info,': iprm(2) must satisfy 1 <= iprm(2) <= 16'
    ELSE IF(info.EQ.-15) THEN
        WRITE(6,*)'=====k2l info',info,': iprm(11) must satisfy 2 <= iprm(11) <= n'
    ELSE IF(info.EQ.-16) THEN
        WRITE(6,*)'=====k2l info',info,': iprm(12) must satisfy 0 <= iprm(12) <= 64'
    ELSE IF(info.EQ.-17) THEN
        WRITE(6,*)'=====k2l info',info,': iprm(13) must satisfy iprm(10) <= iprm(13) <= n'
    ELSE IF(info.EQ.-18) THEN
        WRITE(6,*)'=====k2l info',info,': s_lower and s_upper must satisfy s_lower < s_upper'
    ELSE IF(info.EQ.-21) THEN
        WRITE(6,*)'=====k2l info',info,': Select a sparse direct linear solver'                 
    ELSE IF(info.EQ.-22) THEN
        WRITE(6,*)'=====k2l info',info,': Select a bracketing algorithm'         
    ELSE IF(info.EQ.1) THEN
        WRITE(6,*)'=====k2l info',info,': Failed to set an initial interval'
    ELSE IF(info.EQ.2) THEN
        WRITE(6,*)'=====k2l info',info,': Failed to narrow down the interval'
    ELSE IF(info.EQ.3) THEN
        WRITE(6,*)'=====k2l info',info,': Failed to compute eipenpairs of the interval'
    ELSE IF((info.GE.11).AND.(info.LT.20)) THEN
        WRITE(6,*)'=====k2l info',info,': Error in LDL factorization'
    ELSE IF(info.EQ.21) THEN
        WRITE(6,*)'=====k2l info',info,': Error in bracketing algorithm'        
    END IF
!-----------------------------------------------------------------------
    IF(info.NE.0) THEN
        STOP
    END IF
!    
    RETURN
    END SUBROUTINE k2l_info
!=======================================================================
    SUBROUTINE k2l_summary(k2l_io,k2l_int,k2l_c)
    IMPLICIT NONE
    TYPE(k2l_io_type) :: k2l_io
    TYPE(k2l_int_type) :: k2l_int
    TYPE(k2l_convtest_type), OPTIONAL :: k2l_c
!
    INTEGER :: ounit=6,i
    REAL :: rtmp
!-----------------------------------------------------------------------
!
    WRITE(ounit,*)
    WRITE(ounit,*)'========================================================================'
    WRITE(ounit,*)'====k2l summary========================================================='
    WRITE(ounit,*)'========================================================================'
    WRITE(ounit,*)
    WRITE(ounit,*)'-----Computation time (s)'
    rtmp=0.0E0
    DO i=1,10
        rtmp=rtmp+k2l_int%cmpt_time(i)
    END DO
    WRITE(ounit,"(A50,F15.7,/)") 'Overall                                          :',rtmp
    WRITE(ounit,"(A50,F15.7)") 'Step 1 Task (i)   symbolic factor. A-sB          :',k2l_int%cmpt_time(1)
    WRITE(ounit,"(A50,F15.7)") 'Step 1 Task (ii)  symbolic & numer. factor. B    :',k2l_int%cmpt_time(2)
    WRITE(ounit,"(A50,F15.7)") 'Step 1 Task (iii) Lanczos                        :',k2l_int%cmpt_time(3)
    WRITE(ounit,"(A50,F15.7)") 'Step 1 Task (iv)  inertia computation            :',k2l_int%cmpt_time(4)
    WRITE(ounit,"(A50,F15.7)") 'Step 2 Task (v)   inertia computation            :',k2l_int%cmpt_time(5)
    WRITE(ounit,"(A50,F15.7)") 'Step 3 Task (vi)  symbolic & numer. factor. A-sB :',k2l_int%cmpt_time(6)    
    WRITE(ounit,"(A50,F15.7,/)") 'Step 3 Task (vii) SI Lanczos                     :',k2l_int%cmpt_time(7)
!
    WRITE(ounit,*)'-----Iteration count'
    WRITE(ounit,"(A8,I4)") 'Step 1 :',k2l_int%icnt(1)
    WRITE(ounit,"(A8,I4)") 'Step 2 :',k2l_int%icnt(2)
    WRITE(ounit,"(A8,I4,/)") 'Step 3 :',k2l_int%icnt(4)
!
    WRITE(ounit,*)'-----Step 1: set an initial interval'
    WRITE(ounit,"(A9,2X,A12,2X,A12,2X,A28,2X,A28)")'Iteration','InertiaLower','InertiaUpper',&
    &   'ShiftLower                  ','ShiftUpper                        '
    DO i=1,k2l_int%icnt(1)
        WRITE(ounit,"(I9,2X,I12,2X,I12,2X,E28.18,2X,E28.18)") &
        &   i,k2l_int%inertia_1_lower(i),k2l_int%inertia_1_upper(i),&
        &   k2l_int%shift_1_lower(i),k2l_int%shift_1_upper(i)
    END DO
!
    IF(k2l_io%info.EQ.1) THEN
        RETURN
    END IF
!
    WRITE(ounit,*)
    WRITE(ounit,*) '-----Step 2: narrow down the interval'
    WRITE(ounit,"(A9,2X,A12,2X,A12,2X,A28,2X,A28)")'Iteration','InertiaLower','InertiaUpper',&
    &   'ShiftLower                  ','ShiftUpper                  '
    DO i=0,k2l_int%icnt(2)
        WRITE(ounit,"(I9,2X,I12,2X,I12,2X,E28.18,2X,E28.18)") &
        &   i,k2l_int%inertia_2_lower(i),k2l_int%inertia_2_upper(i),&
        &   k2l_int%shift_2_lower(i),k2l_int%shift_2_upper(i)
    END DO
!
    IF(k2l_io%info.EQ.2) THEN
        RETURN
    END IF
!
    WRITE(ounit,*)
    WRITE(ounit,*) '-----Step 3: select the midpoint'
    WRITE(ounit,"(A8,2X,A28)")'Inertia','Shift                       '
    WRITE(ounit,"(I8,2X,E28.18)") &
    &   k2l_int%inertia_3(1),k2l_int%shift_3(1)
!
    IF(k2l_io%info.EQ.3) THEN
        RETURN
    END IF
!
    IF(PRESENT(k2l_c)) THEN
        WRITE(ounit,*)
        WRITE(ounit,*) '-----Step 3: compute eipenpairs of the interval'
        WRITE(ounit,"(A8,2X,A28,2X,A28,2X,A28)")'Index   ','Eigenvalue                  ',&
        &   'Relative residual 2-norm    ','Relative difference 2-norm    '
        DO i=1,k2l_int%icnt(3)
            WRITE(ounit,"(I8,2X,E28.18,2X,E28.18,2X,E28.18)") &
            &   k2l_int%inertia_2_lower(k2l_int%icnt(2))+i, &
            &   k2l_c%egnval(i),k2l_c%resnorm(i),k2l_c%difnorm(i)
        END DO
    END IF
!    
    RETURN
    END SUBROUTINE k2l_summary
!=======================================================================    
    SUBROUTINE k2l_crd2crs(k2l_io,k2l_int)
    IMPLICIT NONE
    TYPE(k2l_io_type) :: k2l_io
    TYPE(k2l_int_type) :: k2l_int
!-----------------------------------------------------------------------
    ALLOCATE(k2l_int%row_pntr_a(k2l_io%n+1),k2l_int%col_indx_a(k2l_io%nz_a),k2l_int%a(k2l_io%nz_a)) 

    CALL k2l_crd2crs2(k2l_io%n,k2l_io%nz_a,k2l_io%indx_a,k2l_io%jndx_a,k2l_io%rval_a,&
        &   k2l_int%row_pntr_a,k2l_int%col_indx_a,k2l_int%a)    
!    
    IF(ALLOCATED(k2l_io%indx_b).AND.ALLOCATED(k2l_io%jndx_b).AND.ALLOCATED(k2l_io%rval_b)) THEN
        ALLOCATE(k2l_int%row_pntr_b(k2l_io%n+1),k2l_int%col_indx_b(k2l_io%nz_b),k2l_int%b(k2l_io%nz_b))   
        CALL k2l_crd2crs2(k2l_io%n,k2l_io%nz_b,k2l_io%indx_b,k2l_io%jndx_b,k2l_io%rval_b,&
            &   k2l_int%row_pntr_b,k2l_int%col_indx_b,k2l_int%b)
    END IF
!
    RETURN
    END SUBROUTINE k2l_crd2crs
!=======================================================================    
    SUBROUTINE k2l_crd2crs2(n,nz,indx,jndx,rval,row_pntr,col_indx,a)
    IMPLICIT NONE
    INTEGER :: n,nz
    INTEGER, DIMENSION(:) :: row_pntr
    INTEGER, DIMENSION(:) :: indx,jndx,col_indx
    DOUBLE PRECISION, DIMENSION(:) :: rval,a
!
!   local arguments
    INTEGER, ALLOCATABLE, DIMENSION(:) :: temp_indx,temp_jndx
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: temp_rval
!
    INTEGER :: i,k,r,c,nzd
!-----------------------------------------------------------------------
!   Initialize
    IF((n+1.NE.SIZE(row_pntr)).OR.(nz.NE.SIZE(indx)).OR.(nz.NE.SIZE(jndx))  &
    & .OR.(nz.NE.SIZE(col_indx)).OR.(nz.NE.SIZE(rval)).OR.(nz.NE.SIZE(a))) THEN
        WRITE(*,*) 'Incompatible array size'
        STOP
    END IF
!
    ALLOCATE(temp_indx(nz),temp_jndx(nz),temp_rval(nz))
!
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
!-----------------------------------------------------------------------        
    DEALLOCATE(temp_indx,temp_jndx,temp_rval)
    RETURN
    END SUBROUTINE k2l_crd2crs2
!=======================================================================    
    END MODULE k2l_utility
    