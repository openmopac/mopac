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

module mod_vars_cuda

  integer, parameter :: nthreads_gpu = 256, nblocks_gpu = 256
  logical :: lgpu = .false.
  real, parameter :: real_cuda = selected_real_kind(8)
  integer, parameter :: prec = 8
  integer :: ngpus,gpu_id
  logical, parameter :: exe_gpu_kepler = .true.
end module mod_vars_cuda

