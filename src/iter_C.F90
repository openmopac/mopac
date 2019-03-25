      module iter_C 
      USE vast_kind_param, ONLY:  double 
      real(double), dimension(:), allocatable :: pold, pold2, pold3, pbold, &
      pbold2, pbold3,  pgasa, pgasb, psona, psonb
!       
! Arrays used in interp.F90
!
      real(double), dimension(:), allocatable :: h_ai, vecl_ai, h_bi, vecl_bi
      real(double), dimension(:,:), allocatable :: vec_ai, fock_ai, &
       & p_ai, vec_bi, fock_bi, p_bi

      end module iter_C 
