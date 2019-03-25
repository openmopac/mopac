      subroutine dang(a1, a2, b1, b2, rcos) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:33:27  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double) , intent(inout) :: a1 
      real(double) , intent(inout) :: a2 
      real(double) , intent(inout) :: b1 
      real(double) , intent(inout) :: b2 
      real(double) , intent(out) :: rcos 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double) :: zero, anorm, bnorm, sinth, costh 
!-----------------------------------------------
!*********************************************************************
!
!    DANG  DETERMINES THE ANGLE BETWEEN THE POINTS (A1,A2), (0,0),
!          AND (B1,B2).  THE RESULT IS PUT IN RCOS.
!
!*********************************************************************
      zero = 1.0D-6 
      if (abs(a1) >= zero .or. abs(a2) >= zero) then 
        if (abs(b1) >= zero .or. abs(b2) >= zero) then  
          anorm = 1.0D0/sqrt(a1**2 + a2**2) 
          bnorm = 1.0D0/sqrt(b1**2 + b2**2) 
          a1 = a1*anorm 
          a2 = a2*anorm 
          b1 = b1*bnorm 
          b2 = b2*bnorm 
          sinth = a1*b2 - a2*b1 
          costh = a1*b1 + a2*b2 
          costh = min(1.0D0,costh) 
          costh = dmax1(-1.0D0,costh) 
          rcos = acos(costh) 
          if (abs(rcos) >= 4.0D-5) then 
            if (sinth > 0.D0) rcos = 6.28318530717959D0 - rcos 
            rcos = -rcos 
            return  
          endif 
        endif 
      endif 
      rcos = 0.0D0 
      return  
      end subroutine dang 
