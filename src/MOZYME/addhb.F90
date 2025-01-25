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

subroutine addhb (nocc1, nvir1,  idiagg, nij, nhb)
   !***********************************************************************
   !
   !   CHECK FOCK MATRIX AND DENSITY MATRIX TO ENSURE THAT ALL
   !   FOCK MATRIX ELEMENTS ARE BEING CORRECTLY HANDLED BY THE
   !   DIAGONALISER.
   !
   !**********************************************************************
   !
   !  nocc1    Number of occupied M.O.s (might be a subset of noccupied)
   !  nvir1    Number of virtual M.O.s (might be a suset of nvirtual)
   !
   !**********************************************************************
    use molkst_C, only : numat, norbs
    use common_arrays_C, only: eigs, f
    implicit none
    integer, intent (in) :: idiagg, nhb, nocc1, nvir1
    integer, intent (out) :: nij
    integer :: alloc_stat
    logical, dimension (:), allocatable :: latom_loc
    double precision, dimension(:), allocatable :: storei, storej
    integer, allocatable, dimension(:) :: iused
    double precision, dimension (4) :: hblims
    data hblims / 1.d0, 0.1d0, 0.01d0, 0.001d0 /
    allocate (latom_loc(numat), iused(Max (norbs, numat)), &
         & storei(norbs), storej(norbs), &
         & stat=alloc_stat)
    if (alloc_stat /= 0) then
      call memory_error ("addhb")
      return
    end if
    call hbonds (f, nocc1, nvir1, iused, nij, hblims(nhb))
    if (nij /= 0) then
      call diagg2 (nocc1, nvir1, eigs(nocc1+1:), iused, latom_loc, &
      nij, idiagg, storei, storej)
    end if
    deallocate (latom_loc, iused, storei, storej)
    if (alloc_stat /= 0) then
      call memory_error ("addhb")
    end if
end subroutine addhb
