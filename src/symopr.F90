      subroutine symopr(numat, coord, jump, r) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:18:55  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: numat 
      integer , intent(in) :: jump 
      real(double) , intent(inout) :: coord(3,numat) 
      real(double) , intent(in) :: r(3,3) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j
      real(double), dimension(3) :: help 
!-----------------------------------------------
!***********************************************************************
!
!   SYMOPR performs the symmetry operation stored in the 3 by 3 array
!          R on the coordinates.   If JUMP is +1, the operation is
!          X.R, if JUMP is -1, the operation is R.X
!
!   On output, COORD = COORD.R or R.COORD
!
!***********************************************************************
      if (jump >= 0) then 
        do i = 1, numat 
          help = coord(:,i) 
          do j = 1, 3 
            coord(j,i) = 0.D0 
            coord(j,i) = coord(j,i) + sum(r(:,j)*help) 
          end do 
        end do 
        return  
      endif 
      do i = 1, numat 
        help = coord(:,i) 
        do j = 1, 3 
          coord(j,i) = 0.D0 
          coord(j,i) = coord(j,i) + sum(r(j,:)*help) 
        end do 
      end do 
      return  
      end subroutine symopr 
