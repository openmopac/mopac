      MODULE nuchar_I   
      INTERFACE
      SUBROUTINE nuchar (LINE, L_LINE, VALUE, NVALUE) 
      CHARACTER, INTENT(INOUT) :: LINE*(*)
      integer, intent(in) :: L_Line
      DOUBLE PRECISION, DIMENSION(40), INTENT(OUT) :: VALUE 
      INTEGER, INTENT(OUT) :: NVALUE 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
