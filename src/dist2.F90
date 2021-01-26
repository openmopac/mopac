      double precision function dist2 (a, b) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!
!     DETERMINE DISTANCES BETWEEN NEIGHBORING ATOMS
!
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision , intent(in) :: a(3) 
      double precision , intent(in) :: b(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      dist2 = (a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2 
      return  
      end function dist2 


      double precision function dot1 (a, b) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision , intent(in) :: a(3) 
      double precision , intent(in) :: b(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      dot1 = a(1)*b(1) + a(2)*b(2) + a(3)*b(3) 
      return  
      end function dot1 
