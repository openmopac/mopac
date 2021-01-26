      MODULE delri_I   
      INTERFACE
      SUBROUTINE delri (DG, NI, NJ, RR, DEL1) 
      DOUBLE PRECISION, DIMENSION(22), INTENT(OUT) :: DG 
      integer, INTENT(IN) :: NI 
      integer, INTENT(IN) :: NJ 
      DOUBLE PRECISION, INTENT(IN) :: RR 
      DOUBLE PRECISION, INTENT(IN) :: DEL1 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
