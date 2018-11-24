    MODULE k2l_stage
    USE k2l_iotype
    USE k2l_inttype
    USE k2l_utility
    USE k2l_solver
    CONTAINS
!=======================================================================
    SUBROUTINE k2l_setinterval(k2l_io,k2l_int,k2l_factor)
    USE k2l_lanczos_mod    
    IMPLICIT NONE
    TYPE(k2l_io_type) :: k2l_io
    TYPE(k2l_int_type) :: k2l_int
    TYPE(k2l_factor_type) :: k2l_factor
!
    TYPE(k2l_factor_type) :: k2l_factor_B    
    TYPE(k2l_lanczos_type) :: k2l_l
    LOGICAL :: userinterval=.FALSE.
    INTEGER :: imax_l
!
    INTEGER :: i,clock_1,clock_2,clock_rate,clock_max
    INTEGER(KIND=8) :: itmp8
!-----------------------------------------------------------------------        
!   Initialize
    IF(TRIM(ADJUSTL(k2l_io%cprm(2))).EQ.'user') userinterval=.TRUE.
!  
    IF(userinterval) THEN
        imax_l=1
        k2l_int%icnt(1)=1
    ELSE IF(.NOT.userinterval) THEN
        imax_l=k2l_io%iprm(11)
        k2l_l%its=1
        k2l_l%ijob=1
        ALLOCATE(k2l_l%alpha(imax_l),k2l_l%theta(imax_l),k2l_l%beta(0:imax_l),&
            &   k2l_l%v(k2l_io%n,imax_l+1),k2l_l%bv(k2l_io%n,0:imax_l+1))
    END IF
    ALLOCATE(k2l_int%shift_1_lower(imax_l),k2l_int%shift_1_upper(imax_l), &
        &   k2l_int%inertia_1_lower(imax_l),k2l_int%inertia_1_upper(imax_l))
    k2l_int%shift_1_lower=0.0D0
    k2l_int%shift_1_upper=0.0D0
    k2l_int%inertia_1_lower=0
    k2l_int%inertia_1_upper=0
!
!   Initialize sparse solver
    CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
    CALL k2l_ldl_initialize(k2l_io,k2l_factor,0.0D0,.TRUE.)
    CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
    k2l_int%cmpt_time(1)=k2l_int%cmpt_time(1)+REAL(clock_2-clock_1)/REAL(clock_rate)    
!
!-----------------------------------------------------------------------            
!   Check user-specificed initial interval contains k_lower to k_upper-th eigenvalues
    IF(userinterval) THEN
        k2l_int%shift_1_lower(1)=k2l_io%s_lower
        k2l_int%shift_1_upper(1)=k2l_io%s_upper
!
!       Compute inertia at lower endpoint
        CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
        CALL k2l_ldl_inertia(k2l_io,k2l_factor,itmp8,k2l_int%shift_1_lower(1),k2l_int%inertia_1_lower(1))
        CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
        k2l_int%cmpt_time(4)=k2l_int%cmpt_time(4)+REAL(clock_2-clock_1)/REAL(clock_rate)
!       Update the total number of LDL factorizations performed
        k2l_int%icnt(5)=k2l_int%icnt(5)+1
!       Update the number of nonzero elements in the factor L of LDL factorizations
        k2l_int%icnt(6)=k2l_int%icnt(6)+itmp8
!
        IF(k2l_int%inertia_1_lower(1).GE.k2l_io%k_lower) THEN
            CALL k2l_ldl_finalize(k2l_io,k2l_factor)
            k2l_io%info=1
            IF(TRIM(ADJUSTL(k2l_io%cprm(1))).EQ.'print') CALL k2l_summary(k2l_io,k2l_int)
            CALL k2l_info(k2l_io%info)            
        END IF        
!
!       Compute inertia at upper endpoint        
        CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
        CALL k2l_ldl_inertia(k2l_io,k2l_factor,itmp8,k2l_int%shift_1_upper(1),k2l_int%inertia_1_upper(1))
        CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
        k2l_int%cmpt_time(4)=k2l_int%cmpt_time(4)+REAL(clock_2-clock_1)/REAL(clock_rate)
!       Update the total number of LDL factorizations performed
        k2l_int%icnt(5)=k2l_int%icnt(5)+1
!       Update the number of nonzero elements in the factor L of LDL factorizations
        k2l_int%icnt(6)=k2l_int%icnt(6)+itmp8
!
        IF(k2l_int%inertia_1_upper(1).LT.k2l_io%k_upper) THEN
            CALL k2l_ldl_finalize(k2l_io,k2l_factor)
            k2l_io%info=1    
            IF(TRIM(ADJUSTL(k2l_io%cprm(1))).EQ.'print') CALL k2l_summary(k2l_io,k2l_int)
            CALL k2l_info(k2l_io%info)
        END IF
!        
        RETURN
    END IF
!-----------------------------------------------------------------------            
!   Set initial interval using Ritz values of Lanczos method
    IF(.NOT.userinterval) THEN
!
!       Factorize matrix B
        CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
        CALL k2l_ldl_initialize(k2l_io,k2l_factor_B)
        CALL k2l_ldl_factor(k2l_io,k2l_factor_B,itmp8)
        CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
        k2l_int%cmpt_time(2)=k2l_int%cmpt_time(2)+REAL(clock_2-clock_1)/REAL(clock_rate)
!       Update the total number of LDL factorizations performed
        k2l_int%icnt(5)=k2l_int%icnt(5)+1
!       Update the number of nonzero elements in the factor L of LDL factorizations
        k2l_int%icnt(6)=k2l_int%icnt(6)+itmp8
!
        CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
!       Initialize Lanczos
        CALL k2l_lanczos(k2l_io%n,imax_l,k2l_l)
        CALL k2l_matvec('B',k2l_int,k2l_l%v(:,k2l_l%its),k2l_l%bv(:,k2l_l%its))                        
!       Perform first iteration of Lanczos
        CALL k2l_lanczos(k2l_io%n,imax_l,k2l_l)
        CALL k2l_matvec('A',k2l_int,k2l_l%v(:,k2l_l%its),k2l_l%bv(:,k2l_l%its+1))
        CALL k2l_lanczos(k2l_io%n,imax_l,k2l_l)
        CALL k2l_ldl_linsol(k2l_io,k2l_factor_B,k2l_l%bv(:,k2l_l%its+1),k2l_l%v(:,k2l_l%its+1))
        CALL k2l_lanczos(k2l_io%n,imax_l,k2l_l)
        CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
        k2l_int%cmpt_time(3)=k2l_int%cmpt_time(3)+REAL(clock_2-clock_1)/REAL(clock_rate)        
!
!       Compute inertia
        CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
        k2l_int%shift_1_lower(1)=k2l_l%theta(1)
        k2l_int%shift_1_upper(1)=k2l_int%shift_1_lower(1)
        CALL k2l_ldl_inertia(k2l_io,k2l_factor,itmp8,k2l_int%shift_1_lower(1),k2l_int%inertia_1_lower(1))
        k2l_int%inertia_1_upper(1)=k2l_int%inertia_1_lower(1)
        CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
        k2l_int%cmpt_time(4)=k2l_int%cmpt_time(4)+REAL(clock_2-clock_1)/REAL(clock_rate)
!       Update the total number of LDL factorizations performed
        k2l_int%icnt(5)=k2l_int%icnt(5)+1
!       Update the number of nonzero elements in the factor L of LDL factorizations
        k2l_int%icnt(6)=k2l_int%icnt(6)+itmp8
!
        IF(k2l_int%inertia_1_upper(1).LT.k2l_io%k_lower) THEN
            DO i=2,imax_l
!               Perform one iteration of Lanczos
                CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
                CALL k2l_lanczos(k2l_io%n,imax_l,k2l_l)
                CALL k2l_matvec('A',k2l_int,k2l_l%v(:,k2l_l%its),k2l_l%bv(:,k2l_l%its+1))                
                CALL k2l_lanczos(k2l_io%n,imax_l,k2l_l)
                CALL k2l_ldl_linsol(k2l_io,k2l_factor_B,k2l_l%bv(:,k2l_l%its+1),k2l_l%v(:,k2l_l%its+1))
                CALL k2l_lanczos(k2l_io%n,imax_l,k2l_l)
                CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
                k2l_int%cmpt_time(3)=k2l_int%cmpt_time(3)+REAL(clock_2-clock_1)/REAL(clock_rate)                
!
!               Compute inertia
                CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
                IF(k2l_int%inertia_1_upper(i-1).LT.k2l_io%k_lower) THEN
                    k2l_int%shift_1_lower(i)=k2l_int%shift_1_upper(i-1)
                    k2l_int%inertia_1_lower(i)=k2l_int%inertia_1_upper(i-1)
                ELSE
                    k2l_int%shift_1_lower(i)=k2l_int%shift_1_lower(i-1)
                    k2l_int%inertia_1_lower(i)=k2l_int%inertia_1_lower(i-1)
                END IF
                k2l_int%shift_1_upper(i)=k2l_l%theta(i)
                CALL k2l_ldl_inertia(k2l_io,k2l_factor,itmp8,k2l_int%shift_1_upper(i),k2l_int%inertia_1_upper(i))
                CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
                k2l_int%cmpt_time(4)=k2l_int%cmpt_time(4)+REAL(clock_2-clock_1)/REAL(clock_rate)
!               Update the total number of LDL factorizations performed
                k2l_int%icnt(5)=k2l_int%icnt(5)+1
!               Update the number of nonzero elements in the factor L of LDL factorizations
                k2l_int%icnt(6)=k2l_int%icnt(6)+itmp8
!
                k2l_int%icnt(1)=i
                IF(i.NE.1) THEN
                    IF(k2l_int%inertia_1_upper(i).GE.k2l_io%k_upper) THEN
                        CALL k2l_ldl_finalize(k2l_io,k2l_factor_B)
                        EXIT
                    ELSE IF(i.EQ.imax_l) THEN
    !                   Failed to find an initial interval
                        CALL k2l_ldl_finalize(k2l_io,k2l_factor)
                        CALL k2l_ldl_finalize(k2l_io,k2l_factor_B)
                        k2l_io%info=1
                        IF(TRIM(ADJUSTL(k2l_io%cprm(1))).EQ.'print') CALL k2l_summary(k2l_io,k2l_int)
                        CALL k2l_info(k2l_io%info)
                    END IF 
                END IF
            END DO
        END IF
!               
        IF(k2l_int%inertia_1_lower(1).GE.k2l_io%k_upper) THEN
            DO i=2,imax_l
!               Perform one iteration of Lanczos
                CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
                CALL k2l_lanczos(k2l_io%n,imax_l,k2l_l)
                CALL k2l_matvec('A',k2l_int,k2l_l%v(:,k2l_l%its),k2l_l%bv(:,k2l_l%its+1))                
                CALL k2l_lanczos(k2l_io%n,imax_l,k2l_l)
                CALL k2l_ldl_linsol(k2l_io,k2l_factor_B,k2l_l%bv(:,k2l_l%its+1),k2l_l%v(:,k2l_l%its+1))                
                CALL k2l_lanczos(k2l_io%n,imax_l,k2l_l)
                CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
                k2l_int%cmpt_time(3)=k2l_int%cmpt_time(3)+REAL(clock_2-clock_1)/REAL(clock_rate)                
!
!               Compute inertia
                CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
                IF(k2l_int%inertia_1_lower(i-1).GE.k2l_io%k_upper) THEN
                    k2l_int%shift_1_upper(i)=k2l_int%shift_1_lower(i-1)
                    k2l_int%inertia_1_upper(i)=k2l_int%inertia_1_lower(i-1)
                ELSE
                    k2l_int%shift_1_upper(i)=k2l_int%shift_1_upper(i-1)
                    k2l_int%inertia_1_upper(i)=k2l_int%inertia_1_upper(i-1)
                END IF
                k2l_int%shift_1_lower(i)=k2l_l%theta(1)
                CALL k2l_ldl_inertia(k2l_io,k2l_factor,itmp8,k2l_int%shift_1_lower(i),k2l_int%inertia_1_lower(i))                
                CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
                k2l_int%cmpt_time(4)=k2l_int%cmpt_time(4)+REAL(clock_2-clock_1)/REAL(clock_rate)
!               Update the total number of LDL factorizations performed
                k2l_int%icnt(5)=k2l_int%icnt(5)+1
!               Update the number of nonzero elements in the factor L of LDL factorizations
                k2l_int%icnt(6)=k2l_int%icnt(6)+itmp8
!
                k2l_int%icnt(1)=i
                IF(i.NE.1) THEN
                    IF(k2l_int%inertia_1_lower(i).LT.k2l_io%k_lower) THEN
                        CALL k2l_ldl_finalize(k2l_io,k2l_factor_B)
                        EXIT
                    ELSE IF(i.EQ.imax_l) THEN
    !                   Failed to find an initial interval
                        CALL k2l_ldl_finalize(k2l_io,k2l_factor)
                        CALL k2l_ldl_finalize(k2l_io,k2l_factor_B)
                        k2l_io%info=1
                        IF(TRIM(ADJUSTL(k2l_io%cprm(1))).EQ.'print') CALL k2l_summary(k2l_io,k2l_int)
                        CALL k2l_info(k2l_io%info)
                    END IF 
                END IF
            END DO
        END IF
!
        IF((k2l_int%inertia_1_lower(1).GE.k2l_io%k_lower).AND.&
            &   (k2l_int%inertia_1_upper(1).LT.k2l_io%k_upper)) THEN
            DO i=2,imax_l
!               Perform one iteration of Lanczos
                CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
                CALL k2l_lanczos(k2l_io%n,imax_l,k2l_l)
                CALL k2l_matvec('A',k2l_int,k2l_l%v(:,k2l_l%its),k2l_l%bv(:,k2l_l%its+1))                
                CALL k2l_lanczos(k2l_io%n,imax_l,k2l_l)
                CALL k2l_ldl_linsol(k2l_io,k2l_factor_B,k2l_l%bv(:,k2l_l%its+1),k2l_l%v(:,k2l_l%its+1))                
                CALL k2l_lanczos(k2l_io%n,imax_l,k2l_l)
                CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
                k2l_int%cmpt_time(3)=k2l_int%cmpt_time(3)+REAL(clock_2-clock_1)/REAL(clock_rate)                
!
!               Compute inertia
                CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
                k2l_int%shift_1_lower(i)=k2l_int%shift_1_lower(i-1)
                k2l_int%shift_1_upper(i)=k2l_int%shift_1_upper(i-1)
                IF(k2l_int%inertia_1_lower(i-1).GE.k2l_io%k_lower) k2l_int%shift_1_lower(i)=k2l_l%theta(1)
                IF(k2l_int%inertia_1_upper(i-1).LT.k2l_io%k_upper) k2l_int%shift_1_upper(i)=k2l_l%theta(i)
!
                k2l_int%inertia_1_lower(i)=k2l_int%inertia_1_lower(i-1)
                k2l_int%inertia_1_upper(i)=k2l_int%inertia_1_upper(i-1)
                IF(k2l_int%inertia_1_lower(i-1).GE.k2l_io%k_lower) THEN
                    CALL k2l_ldl_inertia(k2l_io,k2l_factor,itmp8,k2l_int%shift_1_lower(i),k2l_int%inertia_1_lower(i))
                    k2l_int%icnt(5)=k2l_int%icnt(5)+1
                    k2l_int%icnt(6)=k2l_int%icnt(6)+itmp8
                END IF
                IF(k2l_int%inertia_1_upper(i-1).LT.k2l_io%k_upper) THEN
                    CALL k2l_ldl_inertia(k2l_io,k2l_factor,itmp8,k2l_int%shift_1_upper(i),k2l_int%inertia_1_upper(i))
                    k2l_int%icnt(5)=k2l_int%icnt(5)+1
                    k2l_int%icnt(6)=k2l_int%icnt(6)+itmp8
                END IF
                CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
                k2l_int%cmpt_time(4)=k2l_int%cmpt_time(4)+REAL(clock_2-clock_1)/REAL(clock_rate)
!
                k2l_int%icnt(1)=i
                IF(i.NE.1) THEN
                    IF((k2l_int%inertia_1_lower(1).LT.k2l_io%k_lower).AND.&
                    &   (k2l_int%inertia_1_upper(1).GE.k2l_io%k_upper)) THEN
                        CALL k2l_ldl_finalize(k2l_io,k2l_factor_B)
                        EXIT
                    ELSE IF(i.EQ.imax_l) THEN
    !                   Failed to find an initial interval                        
                        CALL k2l_ldl_finalize(k2l_io,k2l_factor)
                        CALL k2l_ldl_finalize(k2l_io,k2l_factor_B)
                        k2l_io%info=1
                        IF(TRIM(ADJUSTL(k2l_io%cprm(1))).EQ.'print') CALL k2l_summary(k2l_io,k2l_int)
                        CALL k2l_info(k2l_io%info)
                    END IF 
                END IF
            END DO            
        END IF        
!
        RETURN
    END IF
!-----------------------------------------------------------------------            
    END SUBROUTINE k2l_setinterval
!=======================================================================
    SUBROUTINE k2l_narrowinterval(k2l_io,k2l_int,k2l_factor)
    IMPLICIT NONE
    TYPE(k2l_io_type) :: k2l_io
    TYPE(k2l_int_type) :: k2l_int
    TYPE(k2l_factor_type) :: k2l_factor
!
    LOGICAL :: twointervals=.FALSE.
    INTEGER :: imax_b,mmax_b,mmax_b2
    DOUBLE PRECISION :: dtmp
    INTEGER :: i,itmp,clock_1,clock_2,clock_rate,clock_max
    INTEGER(KIND=8) :: itmp8
!-----------------------------------------------------------------------       
!   Initialize
    imax_b=k2l_io%iprm(12)
    itmp=k2l_io%k_upper-k2l_io%k_lower+1  ! Number of target eigenvalues
    mmax_b=MAX(CEILING(k2l_io%dprm(1)*itmp),k2l_io%iprm(10))  ! Stopping criterion
    mmax_b2=MAX(FLOOR((k2l_io%dprm(1)-1.0D0)*0.5D0*itmp),1)   ! Stopping criterion (internal)
    ALLOCATE(k2l_int%shift_2_lower(0:imax_b),k2l_int%shift_2_upper(0:imax_b), &
        &   k2l_int%inertia_2_lower(0:imax_b),k2l_int%inertia_2_upper(0:imax_b))
    ALLOCATE(k2l_int%shift_2_lower2(0:imax_b),k2l_int%shift_2_upper2(0:imax_b), &
        &   k2l_int%inertia_2_lower2(0:imax_b),k2l_int%inertia_2_upper2(0:imax_b))        
    k2l_int%shift_2_lower=0.0D0
    k2l_int%shift_2_upper=0.0D0
    k2l_int%inertia_2_lower=0
    k2l_int%inertia_2_upper=0
    k2l_int%shift_2_lower2=0.0D0
    k2l_int%shift_2_upper2=0.0D0
    k2l_int%inertia_2_lower2=0
    k2l_int%inertia_2_upper2=0
!
    k2l_int%shift_2_lower(0)=k2l_int%shift_1_lower(k2l_int%icnt(1))
    k2l_int%shift_2_upper(0)=k2l_int%shift_1_upper(k2l_int%icnt(1))
    k2l_int%inertia_2_lower(0)=k2l_int%inertia_1_lower(k2l_int%icnt(1))
    k2l_int%inertia_2_upper(0)=k2l_int%inertia_1_upper(k2l_int%icnt(1))
!
!   Bisection (or bisection-like bracketing)
    DO i=1,imax_b+1
        k2l_int%icnt(2)=i-1
        k2l_int%icnt(3)=k2l_int%inertia_2_upper(i-1)-k2l_int%inertia_2_lower(i-1)
        IF(k2l_int%icnt(3).LE.mmax_b) THEN
            CALL k2l_ldl_finalize(k2l_io,k2l_factor)
            EXIT
        ELSE IF(i-1.EQ.imax_b) THEN
!           Failed to narrow down an interval
            CALL k2l_ldl_finalize(k2l_io,k2l_factor)
            k2l_io%info=2    
            IF(TRIM(ADJUSTL(k2l_io%cprm(1))).EQ.'print') CALL k2l_summary(k2l_io,k2l_int)
            CALL k2l_info(k2l_io%info)
        END IF
!
        IF(.NOT.twointervals) THEN            
!           Bisection and inertia computation
            CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
            CALL k2l_bracketing(k2l_io,imax_b,k2l_int%shift_2_lower,k2l_int%shift_2_upper,i,dtmp)
            CALL k2l_ldl_inertia(k2l_io,k2l_factor,itmp8,dtmp,itmp)
!           Update the total number of LDL factorizations performed
            k2l_int%icnt(5)=k2l_int%icnt(5)+1
!           Update the number of nonzero elements in the factor L of LDL factorizations
            k2l_int%icnt(6)=k2l_int%icnt(6)+itmp8
!
            k2l_int%shift_2_lower(i)=k2l_int%shift_2_lower(i-1)
            k2l_int%shift_2_upper(i)=k2l_int%shift_2_upper(i-1)
            k2l_int%inertia_2_lower(i)=k2l_int%inertia_2_lower(i-1)
            k2l_int%inertia_2_upper(i)=k2l_int%inertia_2_upper(i-1)
!
            IF(itmp.GE.k2l_io%k_upper) THEN
                k2l_int%shift_2_upper(i)=dtmp
                k2l_int%inertia_2_upper(i)=itmp
            ELSE IF(itmp.LT.k2l_io%k_lower) THEN
                k2l_int%shift_2_lower(i)=dtmp
                k2l_int%inertia_2_lower(i)=itmp
            ELSE
                twointervals=.TRUE.
                k2l_int%shift_2_lower2(i)=dtmp
                k2l_int%shift_2_upper2(i)=dtmp
                k2l_int%inertia_2_lower2(i)=itmp
                k2l_int%inertia_2_upper2(i)=itmp
            END IF
            CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
            k2l_int%cmpt_time(5)=k2l_int%cmpt_time(5)+REAL(clock_2-clock_1)/REAL(clock_rate)
!
        ELSE IF(twointervals) THEN
!           Bisection and inertia computation
            CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
            k2l_int%shift_2_lower(i)=k2l_int%shift_2_lower(i-1)
            k2l_int%shift_2_upper2(i)=k2l_int%shift_2_upper2(i-1)
            k2l_int%inertia_2_lower(i)=k2l_int%inertia_2_lower(i-1)
            k2l_int%inertia_2_upper2(i)=k2l_int%inertia_2_upper2(i-1)
!            
            IF(k2l_int%inertia_2_upper2(i-1)-k2l_int%inertia_2_lower(i-1).GT.mmax_b2) THEN
                CALL k2l_bracketing(k2l_io,imax_b,k2l_int%shift_2_lower,k2l_int%shift_2_upper2,i,dtmp)
                CALL k2l_ldl_inertia(k2l_io,k2l_factor,itmp8,dtmp,itmp)
!               Update the total number of LDL factorizations performed
                k2l_int%icnt(5)=k2l_int%icnt(5)+1
!               Update the number of nonzero elements in the factor L of LDL factorizations
                k2l_int%icnt(6)=k2l_int%icnt(6)+itmp8
!                
                IF(itmp.GE.k2l_io%k_lower) THEN
                    k2l_int%shift_2_upper2(i)=dtmp
                    k2l_int%inertia_2_upper2(i)=itmp
                ELSE IF(itmp.LT.k2l_io%k_lower) THEN
                    k2l_int%shift_2_lower(i)=dtmp
                    k2l_int%inertia_2_lower(i)=itmp
                END IF
            END IF
!        
            k2l_int%shift_2_lower2(i)=k2l_int%shift_2_lower2(i-1)
            k2l_int%shift_2_upper(i)=k2l_int%shift_2_upper(i-1)
            k2l_int%inertia_2_lower2(i)=k2l_int%inertia_2_lower2(i-1)
            k2l_int%inertia_2_upper(i)=k2l_int%inertia_2_upper(i-1)
!            
            IF(k2l_int%inertia_2_upper(i-1)-k2l_int%inertia_2_lower2(i-1).GT.mmax_b2) THEN
                CALL k2l_bracketing(k2l_io,imax_b,k2l_int%shift_2_lower2,k2l_int%shift_2_upper,i,dtmp)
                CALL k2l_ldl_inertia(k2l_io,k2l_factor,itmp8,dtmp,itmp)
!               Update the total number of LDL factorizations performed
                k2l_int%icnt(5)=k2l_int%icnt(5)+1
!               Update the number of nonzero elements in the factor L of LDL factorizations
                k2l_int%icnt(6)=k2l_int%icnt(6)+itmp8
!                
                IF(itmp.GE.k2l_io%k_upper) THEN
                    k2l_int%shift_2_upper(i)=dtmp
                    k2l_int%inertia_2_upper(i)=itmp
                ELSE IF(itmp.LT.k2l_io%k_upper) THEN
                    k2l_int%shift_2_lower2(i)=dtmp
                    k2l_int%inertia_2_lower2(i)=itmp
                END IF
            END IF
!
            CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
            k2l_int%cmpt_time(5)=k2l_int%cmpt_time(5)+REAL(clock_2-clock_1)/REAL(clock_rate)
        END IF
    END DO
!
    IF(TRIM(ADJUSTL(k2l_io%cprm(6))).EQ.'second') THEN
        k2l_io%k_lower2=k2l_int%inertia_2_lower(k2l_int%icnt(2))+1
        k2l_io%k_upper2=k2l_int%inertia_2_upper(k2l_int%icnt(2))
        k2l_io%s_lower2=k2l_int%shift_2_lower(k2l_int%icnt(2))
        k2l_io%s_upper2=k2l_int%shift_2_upper(k2l_int%icnt(2))
        k2l_io%info=2    
        IF(TRIM(ADJUSTL(k2l_io%cprm(1))).EQ.'print') CALL k2l_summary(k2l_io,k2l_int)
        k2l_io%info=0
        CALL k2l_info(k2l_io%info)
    END IF
!
    RETURN
    END SUBROUTINE k2l_narrowinterval
!=======================================================================
    SUBROUTINE k2l_pair(k2l_io,k2l_int,k2l_factor)
    USE k2l_lanczos_mod
    IMPLICIT NONE
    TYPE(k2l_io_type) :: k2l_io
    TYPE(k2l_int_type) :: k2l_int
    TYPE(k2l_factor_type) :: k2l_factor
!    
    TYPE(k2l_silanczos_type) :: k2l_s
    TYPE(k2l_convtest_type) :: k2l_c
    INTEGER :: imax_s
    INTEGER :: i,j,itmp,clock_1,clock_2,clock_rate,clock_max
    INTEGER(KIND=8) :: itmp8
!-----------------------------------------------------------------------       
    IF(TRIM(ADJUSTL(k2l_io%cprm(6))).EQ.'second') THEN
        RETURN
    END IF
!
!   Initialize
    imax_s=k2l_io%iprm(13)
    k2l_s%its=1
    k2l_s%ijob=1
    k2l_s%icnt=INT(k2l_int%icnt(3))
    k2l_c%ijob=1
    k2l_c%info=0
    ALLOCATE(k2l_s%alpha(imax_s),k2l_s%theta(imax_s),k2l_s%beta(0:imax_s),&
        &   k2l_s%v(k2l_io%n,imax_s+1),k2l_s%bv(k2l_io%n,0:imax_s+1),k2l_s%y(imax_s,imax_s))
    ALLOCATE(k2l_int%shift_3(1),k2l_int%inertia_3(1))
    k2l_int%shift_3=0.0D0
    k2l_int%inertia_3=0
!
    k2l_c%mmax=k2l_io%k_upper-k2l_io%k_lower+1
    k2l_c%mmax=MAX(2*k2l_c%mmax,k2l_io%iprm(10))
    k2l_c%int_length=(k2l_int%shift_2_upper(k2l_int%icnt(2))-k2l_int%shift_2_lower(k2l_int%icnt(2)))*0.5D0
    k2l_c%tol_res=1.0D1**(-DBLE(2*k2l_io%iprm(1)))
    k2l_c%tol_dif=1.0D1**(-DBLE(2*k2l_io%iprm(2)))    
    ALLOCATE(k2l_c%egnval(k2l_c%mmax),k2l_c%resnorm(k2l_c%mmax),k2l_c%difnorm(k2l_c%mmax), &
    &   k2l_c%egnvec(k2l_io%n,k2l_c%mmax+2),k2l_c%egnvec_old(k2l_io%n,k2l_c%mmax))
    k2l_c%egnvec_old=0.0D0
!
!   Initialize sparse solver
    CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
    CALL k2l_ldl_initialize(k2l_io,k2l_factor,0.0D0,.FALSE.)
    CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
    k2l_int%cmpt_time(6)=k2l_int%cmpt_time(6)+REAL(clock_2-clock_1)/REAL(clock_rate)
!----------------------------------------------------------------------- 
!   Factorize matrix A-sigma*B
    CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
    CALL k2l_bracketing_bisection(k2l_int%shift_2_lower(k2l_int%icnt(2)),k2l_int%shift_2_upper(k2l_int%icnt(2)),k2l_int%shift_3(1))
    CALL k2l_ldl_factor(k2l_io,k2l_factor,itmp8,k2l_int%shift_3(1),k2l_int%inertia_3(1))    
    CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
    k2l_int%cmpt_time(6)=k2l_int%cmpt_time(6)+REAL(clock_2-clock_1)/REAL(clock_rate)
!   Update the total number of LDL factorizations performed
    k2l_int%icnt(5)=k2l_int%icnt(5)+1
!   Update the number of nonzero elements in the factor L of LDL factorizations
    k2l_int%icnt(6)=k2l_int%icnt(6)+itmp8
!
    CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
!   Initialize Lanczos
    CALL k2l_silanczos(k2l_io%n,imax_s,k2l_s)
    CALL k2l_matvec('B',k2l_int,k2l_s%bv(:,k2l_s%its),k2l_s%v(:,k2l_s%its))                        
    CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
    k2l_int%cmpt_time(7)=k2l_int%cmpt_time(7)+REAL(clock_2-clock_1)/REAL(clock_rate)                
!
    DO i=1,imax_s
!       Perform one iteration of Lanczos
        CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
        CALL k2l_silanczos(k2l_io%n,imax_s,k2l_s)
        CALL k2l_ldl_linsol(k2l_io,k2l_factor,k2l_s%v(:,k2l_s%its),k2l_s%bv(:,k2l_s%its+1))                        
        CALL k2l_silanczos(k2l_io%n,imax_s,k2l_s)
        CALL k2l_matvec('B',k2l_int,k2l_s%bv(:,k2l_s%its+1),k2l_s%v(:,k2l_s%its+1))                
        CALL k2l_silanczos(k2l_io%n,imax_s,k2l_s)
        CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
        k2l_int%cmpt_time(7)=k2l_int%cmpt_time(7)+REAL(clock_2-clock_1)/REAL(clock_rate)                
!
!       Convergence test
        CALL SYSTEM_CLOCK(clock_1,clock_rate,clock_max)
        CALL k2l_convtest(imax_s,k2l_s,k2l_c,k2l_int%shift_3(1),i,j)
        DO j=1,k2l_s%icnt
            IF(k2l_c%info.EQ.0) THEN
                CALL k2l_convtest(imax_s,k2l_s,k2l_c,k2l_int%shift_3(1),i,j)
                CALL k2l_matvec('A',k2l_int,k2l_c%egnvec(:,j),k2l_c%egnvec(:,k2l_c%mmax+1))                
                CALL k2l_matvec('B',k2l_int,k2l_c%egnvec(:,j),k2l_c%egnvec(:,k2l_c%mmax+2))                
                CALL k2l_convtest(imax_s,k2l_s,k2l_c,k2l_int%shift_3(1),i,j)
            END IF
        END DO
        IF(k2l_c%info.EQ.0) THEN
            k2l_c%ijob=4
            CALL k2l_convtest(imax_s,k2l_s,k2l_c,k2l_int%shift_3(1),i,j)
        END IF
        CALL SYSTEM_CLOCK(clock_2,clock_rate,clock_max)
        k2l_int%cmpt_time(7)=k2l_int%cmpt_time(7)+REAL(clock_2-clock_1)/REAL(clock_rate)                
!
        k2l_int%icnt(4)=i
        IF(k2l_c%info.EQ.0) THEN
            CALL k2l_ldl_finalize(k2l_io,k2l_factor)
            EXIT
        ELSE IF(i.EQ.imax_s) THEN
!           Failed to be converged
            CALL k2l_ldl_finalize(k2l_io,k2l_factor)
            k2l_io%info=3
            IF(TRIM(ADJUSTL(k2l_io%cprm(1))).EQ.'print') CALL k2l_summary(k2l_io,k2l_int,k2l_c)
            CALL k2l_info(k2l_io%info)
        ELSE IF(k2l_c%info.GE.3) THEN
            k2l_c%egnvec_old(:,1:k2l_s%icnt)=k2l_c%egnvec(:,1:k2l_s%icnt)
        END IF
        k2l_c%ijob=1
        k2l_c%info=0
    END DO
!
!   Get the k_lower-th to k_upper-th eigenpairs
    IF(ALLOCATED(k2l_io%kndx).AND.ALLOCATED(k2l_io%kval).AND.ALLOCATED(k2l_io%kvec)) THEN
        DEALLOCATE(k2l_io%kndx,k2l_io%kval,k2l_io%kvec)
    END IF
    itmp=k2l_io%k_upper-k2l_io%k_lower+1
    ALLOCATE(k2l_io%kndx(itmp),k2l_io%kval(itmp),k2l_io%kvec(k2l_io%n,itmp))
    k2l_io%kndx=k2l_io%k_lower-1
    DO i=1,itmp
        k2l_io%kndx(i)=k2l_io%kndx(i)+i
    END DO        
    k2l_io%kval=k2l_c%egnval(k2l_io%k_lower-k2l_int%inertia_2_lower(k2l_int%icnt(2)): &
    &   k2l_io%k_upper-k2l_int%inertia_2_lower(k2l_int%icnt(2)))
    k2l_io%kvec=k2l_c%egnvec(:,k2l_io%k_lower-k2l_int%inertia_2_lower(k2l_int%icnt(2)): &
    &   k2l_io%k_upper-k2l_int%inertia_2_lower(k2l_int%icnt(2)))
!   
    IF(TRIM(ADJUSTL(k2l_io%cprm(1))).EQ.'print') CALL k2l_summary(k2l_io,k2l_int,k2l_c)
!
    RETURN
    END SUBROUTINE k2l_pair
!=======================================================================
    END MODULE k2l_stage