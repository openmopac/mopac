      MODULE prtdrc_I   
      INTERFACE
      SUBROUTINE prtdrc (DELTT, XPARAM, REF, EKIN, GTOT, ETOT, VELO0&
        , mcoprt, ncoprt, parmax) 
      use molkst_C, only : numat
      DOUBLE PRECISION, INTENT(IN) :: DELTT
      DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: XPARAM 
      DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: REF 
      DOUBLE PRECISION, INTENT(IN) :: EKIN 
      DOUBLE PRECISION, INTENT(INOUT) :: GTOT 
      DOUBLE PRECISION, INTENT(INOUT) :: ETOT 
      DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: VELO0 
      logical, intent(in) :: parmax
      INTEGER, INTENT(IN) :: ncoprt
      integer, dimension(3,numat), intent(in) :: mcoprt
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
