      subroutine osinv(a, n, d) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  10:47:31  03/09/06  
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
      integer , dimension(n) :: l, m 
      integer :: nk, k, kk, j, iz, i, ij, ki, ji, jp, jk, ik, kj, jq, jr 
      real(double) :: tol, biga, holo 
!-----------------------------------------------
!***********************************************************************
!
!    OSINV INVERTS A GENERAL SQUARE MATRIX OF ORDER UP TO MAXORB. SEE
!          DIMENSION STATEMENTS BELOW.
!
!   ON INPUT       A = GENERAL SQUARE MATRIX STORED LINEARLY.
!                  N = DIMENSION OF MATRIX A.
!                  D = VARIABLE, NOT DEFINED ON INPUT.
!
!   ON OUTPUT      A = INVERSE OF ORIGINAL A.
!                  D = DETERMINANT OF ORIGINAL A, UNLESS A WAS SINGULAR,
!                      IN WHICH CASE D = 0.0
!
!***********************************************************************
!***********************************************************************
!
!    IF THE VALUE OF TOL GIVEN HERE IS UNSUITABLE, IT CAN BE CHANGED.
      tol = 1.D-8 
!
!
!***********************************************************************
      d = 1.D0 
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
!
!     10 FOLLOWS
!
            if (abs(biga) - abs(a(ij)) >= 0.D0) cycle  
            biga = a(ij) 
            l(k) = i 
            m(k) = j 
          end do 
        end do 
        j = l(k) 
        if (j - k > 0) then 
          ki = k - n 
          do i = 1, n 
            ki = ki + n 
            holo = -a(ki) 
            ji = ki - k + j 
            a(ki) = a(ji) 
            a(ji) = holo 
          end do 
        endif 
        i = m(k) 
        if (i - k > 0) then 
          jp = n*(i - 1) 
          do j = 1, n 
            jk = nk + j 
            ji = jp + j 
            holo = -a(jk) 
            a(jk) = a(ji) 
            a(ji) = holo 
          end do 
        endif 
        if (abs(biga) - tol < 0.D0) then 
          d = 0.D0 
          return  
        endif 
        do i = 1, n 
          if (i - k == 0) cycle  
          ik = nk + i 
          a(ik) = a(ik)/(-biga) 
        end do 
        do i = 1, n 
          ik = nk + i 
          ij = i - n 
          if (i - k == 0) then 
          else 
            do j = 1, n 
              ij = ij + n 
              if (j - k == 0) cycle  
              kj = ij - i + k 
              a(ij) = a(ik)*a(kj) + a(ij) 
            end do 
          endif 
        end do 
        kj = k - n 
        do j = 1, n 
          kj = kj + n 
          if (j - k == 0) cycle  
          a(kj) = a(kj)/biga 
        end do 
        d = min(d*biga,1.D10) 
        a(kk) = 1.D0/biga 
      end do 
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
          holo = a(jk) 
          ji = jr + j 
          a(jk) = -a(ji) 
          a(ji) = holo 
        end do 
      endif 
      j = m(k) 
      if (j - k <= 0) go to 190 
      ki = k - n 
      do i = 1, n 
        ki = ki + n 
        holo = a(ki) 
        ji = ki + j - k 
        a(ki) = -a(ji) 
        a(ji) = holo 
      end do 
      go to 190 
  260 continue 
      return  
!
      end subroutine osinv 
