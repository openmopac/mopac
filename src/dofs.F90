      subroutine dofs(eref, mono3, n, dd, m, bottom, top) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use chanel_C, only : iw
!...Translated by Pacific-Sierra Research 77to90  4.4G  10:47:10  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: mono3 
      integer , intent(in) :: n 
      integer , intent(in) :: m 
      real(double) , intent(in) :: bottom 
      real(double) , intent(in) :: top 
      real(double) , intent(inout) :: eref(mono3,n) 
      real(double) , intent(inout) :: dd(m) 

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, ii, k
      real(double) :: range, x, b, a, spread, partj, partk 
!-----------------------------------------------
      dd(:500) = 0.D0 
!
!   SPREAD OUT THE ENERGIES OVER THE ENERGY SPECTRUM, TOP TO BOTTOM
!
      range = m/(top - bottom) 
      do j = 1, mono3 
        do i = 1, n 
          x = eref(j,i) 
          if (x<bottom .or. x>top) x = -1.D7 
          eref(j,i) = (x - bottom)*range 
        end do 
      end do 
      do ii = 1, mono3 
        do i = 2, n 
          b = eref(ii,i-1) 
          if (b < 1) cycle  
          a = eref(ii,i) 
          if (a < 1) cycle  
          if (b > a) then 
            x = b 
            b = a 
            a = x 
          endif 
          j = int(b) 
          k = int(a) 
!
! IF J EQUALS K THE INTERVAL FALLS WITHIN ONE BIN
!
          if (j == k) then 
            dd(k) = dd(k) + 1.D0 
          else 
            spread = 1.D0/(a - b + 1.D-12) 
            partj = (j + 1 - b)*spread 
            partk = (a - k)*spread 
            dd(j) = dd(j) + partj 
            dd(k) = dd(k) + partk 
!
! IF K EQUALS J+1 THE INTERVAL STRADDLES TWO BINS
!
            if (k /= j + 1) then 
!
! IF K IS GREATER THAN J+1 THE INTERVAL COVERS MORE THAN TWO BINS
!
              j = j + 1 
              k = k - 1 
              dd(j:k) = dd(j:k) + spread 
            endif 
          endif 
        end do 
      end do 
      x = m/((n - 1)*(top - bottom)) 
      dd = dd*x 
      write (iw, '(A)') ' NORMALIZED DENSITY OF STATES' 
!
!  THE FIRST 'BIN' HAS LOWER BOUND AT BOTTOM AND UPPER BOUND
!  AT BOTTOM+RANGE, THEREFORE THE FIRST 'BIN' IS FOR BOTTOM+0.5*RANGE
!  THE LAST 'BIN' HAS BOUNDS TOP-RANGE AND TOP,
!  THEREFOR THE LAST 'BIN' IS FOR TOP-0.5*RANGE
      range = m/(top - bottom) 
      do i = 1, m 
        write (iw, '(F9.2,F12.6)') bottom + (i - 0.5D0)/range, dd(i) 
      end do 
      return  
      end subroutine dofs 
