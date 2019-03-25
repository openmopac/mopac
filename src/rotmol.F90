      subroutine rotmol(numat, coord, sina, cosa, i, j, r) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double  
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:15:45  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use symopr_I 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: numat 
      integer , intent(in) :: i 
      integer , intent(in) :: j 
      real(double) , intent(in) :: sina 
      real(double) , intent(in) :: cosa 
      real(double)  :: coord(3,numat) 
      real(double)  :: r(3,3) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k 
      real(double) :: buff 
!-----------------------------------------------
      call symopr (numat, coord, -1, r) 
      do k = 1, 3 
        buff = (-sina*r(k,i)) + cosa*r(k,j) 
        r(k,i) = cosa*r(k,i) + sina*r(k,j) 
        r(k,j) = buff 
      end do 
      call symopr (numat, coord, 1, r) 
      return  
      end subroutine rotmol 
