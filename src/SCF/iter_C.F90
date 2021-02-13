      module iter_C 
      double precision, dimension(:), allocatable :: pold, pold2, pold3, pbold, &
      pbold2, pbold3,  pgasa, pgasb, psona, psonb
!       
! Arrays used in interp.F90
!
      double precision, dimension(:), allocatable :: h_ai, vecl_ai, h_bi, vecl_bi
      double precision, dimension(:,:), allocatable :: vec_ai, fock_ai, &
       & p_ai, vec_bi, fock_bi, p_bi

      end module iter_C 
