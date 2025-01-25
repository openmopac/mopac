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

      module rotate_C
!
!   Set up an equivalence between the core-core array ccore and the individual
!   terms, for use in ROTATE
!
      double precision :: css1, csp1, cpps1, cppp1, css2, csp2, cpps2, cppp2
      double precision, dimension(4,2) :: ccore = 0.d0
      equivalence (ccore(1,1), css1), (ccore(2,1), csp1), (ccore(3,1), cpps1), &
      & (ccore(4,1), cppp1), (ccore(1,2), css2), (ccore(2,2), csp2), &
      & (ccore(3,2), cpps2), (ccore(4,2), cppp2)
      end module rotate_C
