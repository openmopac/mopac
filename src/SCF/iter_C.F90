! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

      module iter_C
      double precision, dimension(:), allocatable :: pold, pold2, pold3, pbold, &
      pbold2, pbold3,  pgasa, pgasb, psona, psonb
!
! Arrays used in interp.F90
!
      double precision, dimension(:), allocatable :: h_ai, vecl_ai, h_bi, vecl_bi
      double precision, dimension(:,:), allocatable :: vec_ai, fock_ai, &
       & p_ai, vec_bi, fock_bi, p_bi, pulay_work1, pulay_work2, pulay_work3

      end module iter_C
