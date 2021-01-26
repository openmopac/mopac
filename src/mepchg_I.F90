      MODULE mepchg_I   
      INTERFACE
      SUBROUTINE mepchg (CO, POTPT, NMEP) 
      DOUBLE PRECISION, DIMENSION(3,*), INTENT(IN) :: CO 
      DOUBLE PRECISION, DIMENSION(3,*), INTENT(IN) :: POTPT 
      integer, INTENT(IN) :: NMEP      
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
