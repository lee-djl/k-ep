    MODULE k2l_mod
    USE k2l_iotype
    USE k2l_inttype
    USE k2l_utility
    USE k2l_solver
    USE k2l_stage
    CONTAINS
!=======================================================================    
    SUBROUTINE k2l(k2l_io)
    IMPLICIT NONE    
    TYPE(k2l_io_type) :: k2l_io         ! Input/Output variables
    TYPE(k2l_int_type) :: k2l_int       ! Internal variables
    TYPE(k2l_factor_type) :: &
    &   k2l_factor_12, k2l_factor_3     ! Internal LDL factors
!-----------------------------------------------------------------------
    IF(k2l_io%job.EQ.0) THEN
!       Assign default values to parameters
!       & Allocate IO variables
        CALL k2l_initialize_external(k2l_io)
    END IF
!
    IF(k2l_io%job.EQ.1) THEN
!       Check input parameters
!       & Initialize internal variables
!       & Convert coordinate format to Compressed Row Storage (CRS) format     
        CALL k2l_initialize_internal(k2l_io,k2l_int)
!-----------------------------------------------------------------------    
!       1st Stage: Set an intial interval    
        CALL k2l_setinterval(k2l_io,k2l_int,k2l_factor_12)
!    
!       2nd Stage: Narrow down the interval    
        CALL k2l_narrowinterval(k2l_io,k2l_int,k2l_factor_12)
!    
!       3rd Stage: Compute eigenpairs in the interval    
        CALL k2l_pair(k2l_io,k2l_int,k2l_factor_3)
!
!       Optional: Compute inverse participation ratio
!       for electronic structure calculations
        CALL k2l_ipr(k2l_io,k2l_int)
!-----------------------------------------------------------------------
!       Deallocate internal variables
        CALL k2l_finalize_internal(k2l_int)
    END IF
!
    IF(k2l_io%job.EQ.-1) THEN
!       Deallocate IO variables    
        CALL k2l_finalize_external(k2l_io)
    END IF
!-----------------------------------------------------------------------    
    RETURN
    END SUBROUTINE k2l
!=======================================================================      
    END MODULE k2l_mod
