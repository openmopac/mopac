      real(kind(0.0d0)) function capcor (nat, nfirst, nlast, p, h) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use molkst_C, only : numat
!...Translated by Pacific-Sierra Research 77to90  4.4G  10:47:02  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nat(numat) 
      integer , intent(in) :: nfirst(numat) 
      integer , intent(in) :: nlast(numat) 
      real(double) , intent(in) :: p(*) 
      real(double) , intent(in) :: h(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ni, il, iu, j, ii, k, jl, kk 
      real(double) :: sum 
!-----------------------------------------------
!*****************************************************************
!
!    CORRECTION TO ELECTRONIC ENERGY DUE TO CAPPED BONDS
!
!*****************************************************************
      sum = 0.D0 
      do i = 1, numat 
        ni = nat(i) 
        il = nfirst(i) 
        iu = nlast(i) 
        if (ni == 102) then 
!
!   DO ENTIRE ROW - NO NEED TO CHECK FURTHER.
!
          j = (nlast(i)*(nlast(i)+1))/2 
          ii = iu - 1 
          do k = 1, ii 
            j = j - 1 
            sum = sum + p(j)*h(j) 
          end do 
        else 
          do j = 1, i 
            jl = nfirst(j) 
            if (nat(j) /= 102) cycle  
            do k = il, iu 
              kk = (k*(k - 1))/2 + jl 
              sum = sum + p(kk)*h(kk) 
            end do 
          end do 
        endif 
      end do 
!
!   DOUBLE SUM SINCE WE ONLY CALCULATED LOWER HALF, AND CAPCOR
!   WILL APPEAR IN 1/2*P(H+F).  ONLY H PART OF F WILL BE USED.
      capcor = -sum*2.D0 
      return  
      end function capcor 
