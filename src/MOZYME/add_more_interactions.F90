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

subroutine add_more_interactions()
!
!   In a MOZYME calculation, during a geometry optimization, the distance between a
!  pair of atoms might become less, so that instead of an interaction being considered as a
!  point charge it might need to be considered as a point charge plus dipole, or instead of being a
!  point charge plus dipole, it might need to be considered using the full NDDO approximation.
!
!  Here, the matrix of interactions is re-evaluated and if an iteraction needs to be promoted,
!  it would be promoted.
!
!  The array nijbo is the only array whose contents are changed.
!  other arrays might need to be increased in size.  If so, then they would be increased by
!  20% more than needed, to leave room for future expansion.
!
  use molkst_C, only: mpack, numcal
!
  use common_arrays_C, only: p, f, h
!
  use MOZYME_C, only : lijbo, direct, semidr, partp, parth, partf, &
    rapid
!
  use iter_C, only : pold
!
  implicit none
  double precision, dimension (:), allocatable :: temp_store
!
!  Write out the contents of ijbo at the entry to add_more_interactions
!
  integer :: i, j, alloc_stat, imol = 0
  save imol
!
!  Method is limited to direct SCF
!
  if (.not. direct .or. .not. semidr .or. .not. lijbo) return
  if (imol /= numcal) then
    imol = numcal
    return  !  First time add_more_interactions is called, do nothing
  end if
  i = mpack
  call fillij (.false.)
  if (mpack > i) then
!
! Some arrays are now too small - increase their size, but also preserve their contents
!
    allocate (temp_store(i), stat=alloc_stat)
    if (alloc_stat /= 0) then
        call memory_error ("add_more_interactions")
        return
    end if
    j = Nint(mpack*1.2)
!
    temp_store = pold
    deallocate(pold)
    allocate(pold(j), stat=alloc_stat)
    if (alloc_stat /= 0) then
        call memory_error ("add_more_interactions")
        return
    end if
    pold(1:i) = temp_store(1:i)
    pold(i+1:) = 0.d0
!
    temp_store = p
    deallocate(p)
    allocate(p(j), stat=alloc_stat)
    if (alloc_stat /= 0) then
        call memory_error ("add_more_interactions")
        return
    end if
    p(1:i) = temp_store(1:i)
    p(i+1:) = 0.d0
!
    temp_store = h
    deallocate(h)
    allocate(h(j), stat=alloc_stat)
    if (alloc_stat /= 0) then
        call memory_error ("add_more_interactions")
        return
    end if
    h(1:i) = temp_store(1:i)
    h(i+1:) = 0.d0
!
    temp_store = f
    deallocate(f)
    allocate(f(j), stat=alloc_stat)
    if (alloc_stat /= 0) then
        call memory_error ("add_more_interactions")
        return
    end if
    f(1:i) = temp_store(1:i)
    f(i+1:) = 0.d0
!
!
!
    if (rapid) then
      temp_store = partp
      deallocate(partp)
      allocate(partp(j), stat=alloc_stat)
      if (alloc_stat /= 0) then
          call memory_error ("add_more_interactions")
          return
      end if
      partp(1:i) = temp_store(1:i)
      partp(i+1:) = 0.d0
  !
  !
      temp_store = parth
      deallocate(parth)
      allocate(parth(j), stat=alloc_stat)
      if (alloc_stat /= 0) then
          call memory_error ("add_more_interactions")
          return
      end if
      parth(1:i) = temp_store(1:i)
      parth(i+1:) = 0.d0
  !
      temp_store = partf
      deallocate(partf)
      allocate(partf(j), stat=alloc_stat)
      if (alloc_stat /= 0) then
          call memory_error ("add_more_interactions")
          return
      end if
      partf(1:i) = temp_store(1:i)
      partf(i+1:) = 0.d0
    end if
  end if
  return
end subroutine add_more_interactions
