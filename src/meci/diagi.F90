      double precision function diagi (ialpha, ibeta, eiga, xy, nmos) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!***********************************************************************
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
      double precision , intent(in) :: eiga(*) 
      double precision , intent(in) :: xy(nmos,nmos,nmos,nmos) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j 
      double precision :: x 
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
