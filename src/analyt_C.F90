      module analyt_C 
!
!  Variables used in STO-6G expansion of Slater orbitals
      USE vast_kind_param, ONLY:  double 
      implicit none
      integer, dimension(107) :: nztype 
      integer, dimension(30) :: mtype 
      integer :: ltype 
      real(double), dimension(16) :: ds 
      real(double), dimension(22) :: dg 
      real(double), dimension(100) :: dr 
      real(double), dimension(3) :: tdx, tdy, tdz 
      real(double), dimension(22) :: g 
      real(double), dimension(3) :: tx, ty, tz 
      end module analyt_C 
