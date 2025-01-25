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

subroutine filusp (nat, nfirst, nlast, uspd)
    use molkst_C, only :  norbs, numat
    use parameters_C, only : uss, upp, udd
  !
  !.. Implicit Declarations ..
    implicit none
  !
  !.. Formal Arguments ..
    integer, dimension (numat), intent (in) :: nat, nfirst, nlast
    double precision, dimension (norbs), intent (out) :: uspd
  !
  !.. Local Scalars ..
    integer :: i, ia, ib, j, ni_loc
  !
  ! ... Executable Statements ...
  !
    do i = 1, numat
      ni_loc = nat(i)
      ia = nfirst(i)
      ib = nlast(i)
      if (ia <= ib) then
      !
      !  Atom has an "s" orbital.
      !
        uspd(ia) = uss(ni_loc)
        if (ia /= ib) then
        !
        !  Atom has a "p" shell.
        !
          do j = ia + 1, ia + 3
            uspd(j) = upp(ni_loc)
          end do
        !
        !  Set the "d" shell, if present.
        !
          do j = ia + 4, ib
            uspd(j) = udd(ni_loc)
          end do
        end if
      end if
    end do
end subroutine filusp
