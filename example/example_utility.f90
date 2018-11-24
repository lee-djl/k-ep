    MODULE example_utility
    CONTAINS
!=======================================================================    
    SUBROUTINE example_readmtx(ipath,n,nz,indx,jndx,rval)
    IMPLICIT NONE
    CHARACTER(256) :: ipath
    INTEGER :: n,nz
    INTEGER, ALLOCATABLE, DIMENSION(:) :: indx,jndx
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: rval
!
    CHARACTER(7) :: field
    CHARACTER(10) :: rep
    CHARACTER(19) :: symm
!
    INTEGER :: iunit=999,nz_mx,itmp
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ival
    COMPLEX, ALLOCATABLE, DIMENSION(:) :: cval
!
    PARAMETER (nz_mx=2147483647)
!
!-----------------------------------------------------------------------
    OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(ipath)),STATUS='OLD',&
    &   ACCESS='SEQUENTIAL',ACTION='READ')
    CALL mminfo(iunit,rep,field,symm,itmp,n,nz)
!
    ALLOCATE(indx(nz),jndx(nz),rval(nz),ival(1),cval(1))
    CALL mmread(iunit,rep,field,symm,n,n,nz,nz_mx,&
    &   indx,jndx,ival,rval,cval)
    DEALLOCATE(ival,cval)
!
    CLOSE(iunit)
    RETURN
    END SUBROUTINE example_readmtx
!=======================================================================
    SUBROUTINE example_writekval(n,kndx,kval)
    IMPLICIT NONE
    INTEGER :: n
    INTEGER,DIMENSION(:) :: kndx
    DOUBLE PRECISION,DIMENSION(:) :: kval
!   
    CHARACTER(256) :: opath
    INTEGER :: ounit=888,i
!
!-----------------------------------------------------------------------
!
    IF((n.NE.SIZE(kndx)).OR.(n.NE.SIZE(kval))) THEN
        WRITE(*,*) 'Incompatible array size'
        RETURN
    END IF
!    
    IF(n.EQ.1) THEN
        WRITE(opath,"(A6,I0.8,A4)") './eval',kndx(1),'.txt'
    ELSE
        WRITE(opath,"(A6,I0.8,A2,I0.8,A4)") &
            & './eval',kndx(1),'to', kndx(n),'.txt'        
    END IF
!   
    OPEN(UNIT=ounit,FILE=TRIM(ADJUSTL(opath)),STATUS='REPLACE',&
    &   ACCESS='SEQUENTIAL',ACTION='WRITE')
!
    DO i=1,n
        WRITE(ounit,"(I8,2X,E28.18)") kndx(i),kval(i)
    END DO
!
    CLOSE(ounit)
    RETURN
    END SUBROUTINE example_writekval    
!  
!=======================================================================
    SUBROUTINE example_writekvec(m,n,kndx,kvec)
    IMPLICIT NONE
    INTEGER :: m,n
    INTEGER,DIMENSION(:) :: kndx
    DOUBLE PRECISION,DIMENSION(:,:) :: kvec
!   
    CHARACTER(256) :: opath
    INTEGER :: ounit=777,i,j
!
!-----------------------------------------------------------------------
!
    IF((n.NE.SIZE(kndx)).OR.(m.NE.SIZE(kvec,1)).OR.(n.NE.SIZE(kvec,2))) THEN
        WRITE(*,*) 'Incompatible array size'
        RETURN
    END IF
!
    DO j=1,n
        WRITE(opath,"(A6,I0.8,A4)") './evec',kndx(j),'.txt'
        OPEN(UNIT=ounit,FILE=TRIM(ADJUSTL(opath)),STATUS='REPLACE',&
            &   ACCESS='SEQUENTIAL',ACTION='WRITE')
!        
        DO i=1,m
            WRITE(ounit,"(E28.18)") kvec(i,j)
        END DO
!
        CLOSE(ounit)    
    END DO        
!
    RETURN
    END SUBROUTINE example_writekvec

!=======================================================================
    SUBROUTINE example_writekipr(n,kndx,kipr)
    IMPLICIT NONE
    INTEGER :: n
    INTEGER,DIMENSION(:) :: kndx
    DOUBLE PRECISION,DIMENSION(:) :: kipr
!   
    CHARACTER(256) :: opath
    INTEGER :: ounit=666,i
!
!-----------------------------------------------------------------------
!
    IF((n.NE.SIZE(kndx)).OR.(n.NE.SIZE(kipr))) THEN
        WRITE(*,*) 'Incompatible array size'
        RETURN
    END IF
!    
    IF(n.EQ.1) THEN
        WRITE(opath,"(A5,I0.8,A4)") './ipr',kndx(1),'.txt'
    ELSE
        WRITE(opath,"(A5,I0.8,A2,I0.8,A4)") &
            & './ipr',kndx(1),'to', kndx(n),'.txt'        
    END IF
!   
    OPEN(UNIT=ounit,FILE=TRIM(ADJUSTL(opath)),STATUS='REPLACE',&
    &   ACCESS='SEQUENTIAL',ACTION='WRITE')
!
    DO i=1,n
        WRITE(ounit,"(I8,2X,E28.18)") kndx(i),kipr(i)
    END DO
!
    CLOSE(ounit)
    RETURN
    END SUBROUTINE example_writekipr    
!=======================================================================
    SUBROUTINE example_writekint(k_lower,k_upper,s_lower,s_upper)
    IMPLICIT NONE
    INTEGER :: k_lower,k_upper
    DOUBLE PRECISION :: s_lower,s_upper
!   
    CHARACTER(256) :: opath
    INTEGER :: ounit=666,i
!
!-----------------------------------------------------------------------
!
    IF(k_upper-k_lower.EQ.0) THEN
        WRITE(opath,"(A5,I0.8,A4)") './int',k_lower,'.txt'
    ELSE
        WRITE(opath,"(A5,I0.8,A2,I0.8,A4)") &
            & './int',k_lower,'to', k_upper,'.txt'        
    END IF
!   
    OPEN(UNIT=ounit,FILE=TRIM(ADJUSTL(opath)),STATUS='REPLACE',&
    &   ACCESS='SEQUENTIAL',ACTION='WRITE')
!
    WRITE(ounit,"(2E28.18)") s_lower,s_upper
    WRITE(ounit,"(2I8)") k_lower,k_upper
!
    CLOSE(ounit)
    RETURN
    END SUBROUTINE example_writekint    
!=======================================================================
    END MODULE example_utility
    