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

subroutine buildf (f, partf, mode)
    use molkst_C, only: numat, mpack, id
    use common_arrays_C, only : w, h, wk
    implicit none
    integer, intent (in) :: mode
    double precision, dimension (mpack), intent (in) :: partf
    double precision, dimension (mpack), intent (inout) :: f
    integer :: alloc_stat
    double precision, dimension(:), allocatable :: q, qe
    double precision, dimension(:,:), allocatable :: ptot2
   !
   !  MODE= 1:  Add on terms to a partial Fock matrix to make it complete.
   !  MODE= 0:  Build a Fock matrix, starting with the one-electron matrix.
   !  MODE=-1:  Remove terms from a complete Fock matrix to make a partial.
   !
    allocate (q(numat), qe(numat), ptot2(numat,81), stat=alloc_stat)
    if (alloc_stat /= 0) then
       call memory_error ("buildf")
       goto 100
    end if
    if (mode ==-1) f(:mpack) = partf(:mpack) - h(:mpack)
    if (mode == 0) f(:mpack) =                 h(:mpack)
    if (mode == 1) f(:mpack) = partf(:mpack) + h(:mpack)
    if (id == 0) then
      !
      !  Discrete system
      !
      call fock2z (f, q, qe, w, w, ptot2, mode, 1)
    else
      !
      !  A solid (polymer, layer, or 3-D solid)
      !
      call fock2z (f, q, qe, w, wk, ptot2, mode, 0)
    end if
100 continue
    deallocate (q, qe, ptot2)
end subroutine buildf
