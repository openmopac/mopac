      subroutine minv(a, n, d) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  10:47:25  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      real(double) , intent(out) :: d 
      real(double) , intent(inout) :: a(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(n) :: m, l 
      integer :: nk, k, kk, j, iz, i, ij, ki, ji, jp, jk, ik, kj, jq, jr 
      real(double) :: biga, hold 
!-----------------------------------------------
!*********************************************************************
!
!     INVERT A MATRIX USING GAUSS-JORDAN METHOD.  PART OF DIIS
!     A - INPUT MATRIX (MUST BE A GENERAL MATRIX), DESTROYED IN
!        COMPUTATION AND REPLACED BY RESULTANT INVERSE.
!     N - ORDER OF MATRIX A
!     D - RESULTANT DETERMINANT
!
!*********************************************************************
!
!     SEARCH FOR LARGEST ELEMENT
!
      d = 1.0D0 
      nk = -n 
      do k = 1, n 
        nk = nk + n 
        l(k) = k 
        m(k) = k 
        kk = nk + k 
        biga = a(kk) 
        do j = k, n 
          iz = n*(j - 1) 
          do i = k, n 
            ij = iz + i 
            if (abs(biga) >= abs(a(ij))) cycle  
            biga = a(ij) 
            l(k) = i 
            m(k) = j 
          end do 
        end do 
!
!     INTERCHANGE ROWS
!
        j = l(k) 
        if (j - k > 0) then 
          ki = k - n 
          do i = 1, n 
            ki = ki + n 
            hold = -a(ki) 
            ji = ki - k + j 
            a(ki) = a(ji) 
            a(ji) = hold 
          end do 
        endif 
!
!     INTERCHANGE COLUMNS
!
        i = m(k) 
        if (i - k > 0) then 
          jp = n*(i - 1) 
          do j = 1, n 
            jk = nk + j 
            ji = jp + j 
            hold = -a(jk) 
            a(jk) = a(ji) 
            a(ji) = hold 
          end do 
        endif 
!
!     DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS
!     CONTAINED IN BIGA)
!
        if (biga == 0.D0) then 
          d = 0.D0 
          return  
        endif 
        do i = 1, n 
          if (i - k == 0) cycle  
          ik = nk + i 
          a(ik) = a(ik)/(-biga) 
        end do 
!  REDUCE MATRIX
        do i = 1, n 
          ik = nk + i 
          hold = a(ik) 
          ij = i - n 
          if (i - k == 0) then 
          else 
            do j = 1, n 
              ij = ij + n 
              if (j - k == 0) cycle  
              kj = ij - i + k 
              a(ij) = hold*a(kj) + a(ij) 
            end do 
          endif 
        end do 
!
!     DIVIDE ROW BY PIVOT
!
        kj = k - n 
        do j = 1, n 
          kj = kj + n 
          if (j - k == 0) cycle  
          a(kj) = a(kj)/biga 
        end do 
!
!     PRODUCT OF PIVOTS
!
        d = max(-1.D25,min(1.D25,d)) 
        d = d*biga 
!
!     REPLACE PIVOT BY RECIPROCAL
!
        a(kk) = 1.D0/biga 
      end do 
!
!     FINAL ROW AND COLUMN INTERCHANGE
!
      k = n 
  190 continue 
      k = k - 1 
      if (k <= 0) go to 260 
      i = l(k) 
      if (i - k > 0) then 
        jq = n*(k - 1) 
        jr = n*(i - 1) 
        do j = 1, n 
          jk = jq + j 
          hold = a(jk) 
          ji = jr + j 
          a(jk) = -a(ji) 
          a(ji) = hold 
        end do 
      endif 
      j = m(k) 
      if (j - k <= 0) go to 190 
      ki = k - n 
      do i = 1, n 
        ki = ki + n 
        hold = a(ki) 
        ji = ki - k + j 
        a(ki) = -a(ji) 
        a(ji) = hold 
      end do 
      go to 190 
  260 continue 
      return  
      end subroutine minv 
