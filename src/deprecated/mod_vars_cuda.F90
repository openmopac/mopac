! Molecular Orbital PACkage (MOPAC)
! Copyright (C) 2021, Virginia Polytechnic Institute and State University
!
! MOPAC is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! MOPAC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module mod_vars_cuda

  integer, parameter :: nthreads_gpu = 256, nblocks_gpu = 256
  logical :: lgpu = .false.
  real, parameter :: real_cuda = selected_real_kind(8)
  integer, parameter :: prec = 8
  integer :: ngpus,gpu_id
  logical, parameter :: exe_gpu_kepler = .true.
end module mod_vars_cuda

