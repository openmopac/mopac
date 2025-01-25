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
