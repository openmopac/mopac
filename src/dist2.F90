      real(kind(0.0d0)) function dist2 (a, b) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!
!     DETERMINE DISTANCES BETWEEN NEIGHBORING ATOMS
!
!...Translated by Pacific-Sierra Research 77to90  4.4G  09:34:29  03/16/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double) , intent(in) :: a(3) 
      real(double) , intent(in) :: b(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      dist2 = (a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2 
      return  
      end function dist2 


      real(kind(0.0d0)) function dot1 (a, b) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  09:34:29  03/16/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double) , intent(in) :: a(3) 
      real(double) , intent(in) :: b(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      dot1 = a(1)*b(1) + a(2)*b(2) + a(3)*b(3) 
      return  
      end function dot1 
