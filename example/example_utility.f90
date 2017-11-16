    MODULE mod_example_readmtx
    CONTAINS
    SUBROUTINE example_readmtx(ipath,iunit,n,nz,indx,jndx,rval)
    IMPLICIT NONE
    CHARACTER(256) :: ipath
    INTEGER :: iunit,n,nz
    INTEGER, ALLOCATABLE, DIMENSION(:) :: indx,jndx
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: rval
!
    CHARACTER(7) :: field
    CHARACTER(10) :: rep
    CHARACTER(19) :: symm
!
    INTEGER :: nz_mx,itmp,ival
    COMPLEX :: cval
!
    ALLOCATABLE :: ival,cval
    DIMENSION :: ival(:),cval(:)
!
    PARAMETER (nz_mx=2147483647)
!
!-----------------------------------------------------------------------
    OPEN(UNIT=iunit,FILE=TRIM(ipath),STATUS='OLD',&
    &   ACCESS='SEQUENTIAL',ACTION='READ')
    CALL mminfo(iunit,rep,field,symm,itmp,n,nz)
!
    ALLOCATE(indx(nz),jndx(nz),rval(nz))
    CALL mmread(iunit,rep,field,symm,n,n,nz,nz_mx,&
    &   indx,jndx,ival,rval,cval)
!
    CLOSE(iunit)
    RETURN
    END SUBROUTINE example_readmtx
    END MODULE mod_example_readmtx
!
!=======================================================================
    SUBROUTINE example_writekevec(opath,ounit,n,kevec)
    IMPLICIT NONE
    CHARACTER(256) :: opath
    INTEGER :: ounit,n
    DOUBLE PRECISION :: kevec
    DIMENSION :: kevec(n)
!
    INTEGER :: i
!
!-----------------------------------------------------------------------
!
    OPEN(UNIT=ounit,FILE=TRIM(opath),STATUS='REPLACE',&
    &   ACCESS='SEQUENTIAL',ACTION='WRITE')
!
    DO i=1,n
        WRITE(ounit,"(E28.18)") kevec(i)
    END DO
!
    CLOSE(ounit)
    RETURN
    END SUBROUTINE example_writekevec
