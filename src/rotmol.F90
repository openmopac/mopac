      subroutine rotmol(numat, coord, sina, cosa, i, j, r) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!***********************************************************************
!DECK MOPAC
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
      double precision , intent(in) :: sina 
      double precision , intent(in) :: cosa 
      double precision  :: coord(3,numat) 
      double precision  :: r(3,3) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k 
      double precision :: buff 
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
