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

subroutine diagg (fao, nocc, nvir, idiagg, partp, indi)
!**********************************************************************
!
!  This is the most important operation of MOZYME - the matrix elements
!  involving occupied and virtual LMO's are annihilated by a single pass
!  Euler rotation.
!
!**********************************************************************
    use MOZYME_C, only: icocc_dim, fmo, ifmo, fmo_dim
    use common_arrays_C, only: eigs, p
    use molkst_C, only : mpack, norbs, numat
!
    implicit none
    integer, intent (in) :: idiagg, nocc, nvir, indi
    double precision, dimension (mpack), intent (inout) :: partp
    double precision, dimension (mpack), intent (in) :: fao
!
    integer :: nij, alloc_stat
    logical, dimension (:), allocatable :: latoms
    integer, dimension (:), allocatable :: iused
    double precision, dimension(:), allocatable :: storei, storej, ws, aov, aocc, avir
    allocate (storei(norbs), storej(norbs), ws(norbs), aov(numat), &
         & avir(norbs), aocc(Max(1, icocc_dim)), iused(numat), latoms(numat), &
         & stat=alloc_stat)
    if (alloc_stat /= 0) then
      call mopend("Insufficient memory to run DIAGG")
      return
    end if
    nij = fmo_dim
!**********************************************************************
!
!  In diagg1, all significant matrix elements are constructed.  These
!  nij elements are stored in fmo, and their indices are stored in ifmo.
!
!**********************************************************************

!
    call diagg1 (fao, nocc, nvir,  eigs(nocc+1:), ws, latoms,  ifmo, fmo, &
     & fmo_dim, nij, idiagg, avir, aocc, aov)
!
!**********************************************************************
!
!  In diagg2, the significant matrix elements are annihilated by a two by
!  two rotation.
!
!**********************************************************************

    call diagg2 (nocc, nvir, eigs(nocc+1:), iused, latoms, nij, idiagg, storei, storej)
!
    call density_for_MOZYME (p, indi, nocc, partp)
!
    deallocate (storei, storej, ws, aov, avir, aocc, iused, latoms)
    return
end subroutine diagg
