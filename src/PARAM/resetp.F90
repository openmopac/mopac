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

subroutine resetp (mode, loop)
    use molkst_C, only : mpack, keywrd
    use param_global_C, only : pas, pbs
    use common_arrays_C, only : p, pa, pb
!--------------------------------------------------------------------
    implicit none
    integer, intent (in) :: loop, mode
!--------------------------------------------------------------------
    integer :: i
    integer, save :: ilin
    if (loop == 1) then
      ilin = 0
    end if
    if (mode == 0) then
    !
    !  Restore density matrix from store (PAS and PBS) to current array (P)
    !
      do i = 1, mpack
        pa(i) = pas(i+ilin)
        pb(i) = pbs(i+ilin)
        p(i) = pa(i) + pb(i)
      end do
    else
    !
    !  Store density matrix
    !
      if(index(keywrd," UHF") /=0)then
        do i = 1, mpack
          pas(i+ilin) = pa(i)
          pbs(i+ilin) = pb(i)
        end do
      else
       do i = 1, mpack
          pas(i+ilin) = pa(i)
          pbs(i+ilin) = pa(i)
        end do
      end if
    end if
    ilin = ilin + mpack
end subroutine resetp
