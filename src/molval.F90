      subroutine molval(c, p, rhfuhf) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE molkst_C, only : numat, norbs, mpack
      USE chanel_C, only : iw
      use common_arrays_C, only : nfirst, nlast
!***********************************************************************
!DECK MOPAC
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision, dimension (norbs) :: val
      double precision , intent(in) :: rhfuhf 
      double precision , intent(in) :: c(norbs,norbs) 
      double precision , intent(in) :: p(mpack) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer ::  i, jj, jl, ju, j, kk, kl, ku, k, l1, l2, l 
      double precision :: sum 
!----------------------------------------------- 
      do i = 1, norbs 
        sum = 0.D0 
        do jj = 1, numat 
          jl = nfirst(jj) 
          ju = nlast(jj) 
          do j = jl, ju 
            do kk = 1, numat 
              if (kk == jj) cycle  
              kl = nfirst(kk) 
              ku = nlast(kk) 
              do k = kl, ku 
                l1 = max(j,k) 
                l2 = j + k - l1 
                l = (l1*(l1 - 1))/2 + l2 
                sum = sum + c(j,i)*c(k,i)*p(l) 
              end do 
            end do 
          end do 
        end do 
        val(i) = sum*rhfuhf 
      end do 
      write (iw, '(10F8.4)') (val(i),i=1,norbs) 
      return  
      end subroutine molval 
