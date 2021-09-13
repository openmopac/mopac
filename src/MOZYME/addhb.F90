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
