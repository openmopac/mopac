      module funcon_C 
      USE vast_kind_param, only:  double 
!...Created by Pacific-Sierra Research 77to90  4.4G  08:33:33  03/09/06  
      real(double) :: &
      &  fpc_1,     & !
      &  fpc_2,     & !
      &  a0,        & !
      &  ev,        & !
      &  fpc_5,     & !
      &  fpc_6,     & !
      &  fpc_7,     & !
      &  fpc_8,     & !
      &  fpc_9,     & !
      &  fpc_10   
      double precision, dimension (10) :: fpc
      double precision, parameter :: pi = 3.14159265358979323846d0
      double precision, parameter :: twopi = 2.0d0 * pi
      equivalence (fpc(1), fpc_1), (fpc(2), fpc_2), (fpc(3), a0), (fpc(4), ev), &
      & (fpc(5), fpc_5), (fpc(6), fpc_6), (fpc(7), fpc_7), (fpc(8), fpc_8), &
      & (fpc(9), fpc_9), (fpc(10), fpc_10)
      end module funcon_C 
