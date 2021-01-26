      MODULE initsn_I   
      INTERFACE
      SUBROUTINE initsn (INDEPS, DIRSM, DIRSMH) 
      use cosmo_C, only : nppa 
      integer :: INDEPS 
      DOUBLE PRECISION, DIMENSION(3,NPPA) :: DIRSM 
      DOUBLE PRECISION, DIMENSION(3,NPPA/3) :: DIRSMH 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
