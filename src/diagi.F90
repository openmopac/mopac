      real(kind(0.0d0)) function diagi (ialpha, ibeta, eiga, xy, nmos) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  10:47:08  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nmos 
      integer , intent(in) :: ialpha(nmos) 
      integer , intent(in) :: ibeta(nmos) 
      real(double) , intent(in) :: eiga(*) 
      real(double) , intent(in) :: xy(nmos,nmos,nmos,nmos) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j 
      real(double) :: x 
!-----------------------------------------------
!***********************************************************************
!
!  CALCULATES THE ENERGY OF A MICROSTATE DEFINED BY IALPHA AND IBETA
!
!***********************************************************************
      x = 0.0D0 
      do i = 1, nmos 
        if (ialpha(i) == 0) cycle  
        x = x + eiga(i) 
        do j = 1, nmos 
          x = x + ((xy(i,i,j,j)-xy(i,j,i,j))*ialpha(j)*0.5D0+xy(i,i,j,j)*ibeta(&
            j)) 
        end do 
      end do 
      do i = 1, nmos 
        if (ibeta(i) == 0) cycle  
        x = x + eiga(i) 
        do j = 1, i - 1 
          x = x + (xy(i,i,j,j)-xy(i,j,i,j))*ibeta(j) 
        end do 
      end do 
      diagi = x 
      return  
      end function diagi 
