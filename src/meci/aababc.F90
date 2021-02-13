      double precision  function aababc (iocca1, ioccb1, iocca2, nmos, xy)  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE meci_C, only : occa 
!***********************************************************************
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nmos 
      integer , intent(in) :: iocca1(nmos) 
      integer , intent(in) :: ioccb1(nmos) 
      integer , intent(in) :: iocca2(nmos) 
      double precision , intent(in) :: xy(nmos,nmos,nmos,nmos) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ij, j, k 
      double precision :: sum 
!-----------------------------------------------
!**********************************************************************
!
! AABABC EVALUATES THE C.I. MATRIX ELEMENT FOR TWO MICROSTATES DIFFERING
!       BY ONE ALPHA ELECTRON. THAT IS, ONE MICROSTATE HAS AN ALPHA
!       ELECTRON IN PSI(I) WHICH, IN THE OTHER MICROSTATE, IS IN PSI(J)
!
!**********************************************************************
      do i = 1, nmos 
        if (iocca1(i) == iocca2(i)) cycle  
        exit  
      end do 
      ij = ioccb1(i) 
      do j = i + 1, nmos 
        if (iocca1(j) /= iocca2(j)) exit  
        ij = ij + iocca1(j) + ioccb1(j) 
      end do 
      sum = 0.D0 
      do k = 1, nmos 
        sum = sum + (xy(i,j,k,k)-xy(i,k,j,k))*(iocca1(k)-occa(k)) + &
        xy(i,j,k,k)*(ioccb1(k)-occa(k)) 
      end do 
      if (mod(ij,2) == 1) sum = -sum 
      aababc = sum 
      return  
      end function aababc 
