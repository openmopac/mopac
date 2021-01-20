      module overlaps_C 
      USE vast_kind_param, ONLY:  double 
!...Created by Pacific-Sierra Research 77to90  4.4G  10:47:07  03/09/06  
      real(double) :: &
      cutof1,   &  !  For distances beyond cutof1 overlaps are set to zero
      cutof2       !  For distances beyond cutof2, NDDO is replaced by point-charges 
      real(double), dimension(60,6) :: ccc, zzz 
      real(double), dimension(6,6,2) :: allc, allz 
      integer :: isp, ips 
      real(double), dimension(7) :: a, b 
      real(double) :: sa, sb 
      real(double), dimension(0:17) :: fact    ! Factorials:  fact(n) = n!
      data fact/ 1.d0, 1.D0, 2.D0, 6.D0, 24.D0, 120.D0, 720.D0, 5040.D0, 40320.D0, &
        362880.D0, 3628800.D0, 39916800.D0, 479001600.D0, 6227020800.D0, &
        8.71782912D10, 1.307674368D12, 2.092278989D13, 3.556874281D14/  
      end module overlaps_C 
