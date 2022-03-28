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

  subroutine txtype (jj, jtype, letter)
!
!  Construct the Greek letter for each atom using
!  the PDB convention.  The letter will go in character 15
!  of the entries in array "txtatm"  If there is also a number,
!  it will go in character 16.
!
    use common_arrays_C, only : txtatm, nat, nbonds, ibonds
    implicit none
    character, intent (in) :: letter
    integer, intent (inout) :: jj
    integer, dimension (jj), intent (inout) :: jtype
!
    integer :: i, j, k, loop, m
   !
   !   Remove any duplicates from JTYPE
   !
    j = 1
    k = 1
    outer_loop: do i = 2, jj
      do k = 1, j
        if (jtype(k) == jtype(i)) cycle outer_loop
      end do
      j = j + 1
      jtype(j) = jtype(i)
    end do outer_loop
    jj = j
    do loop = 1, 4
      m = 0
      do i = 1, jj
        if (nat(jtype(i)) /= 1) then
          m = m + 1
          k = i
        end if
      end do
      !
      !   Make symbol for atom.  Greek letter goes in 15:15
      !                            and number goes in 16:16
      !
      if (m == 1) then
        txtatm (jtype(k)) (15:15) = letter
        j = jtype(k)
        m = 0
        do i = 1, nbonds(j)
          if (txtatm(ibonds(i, j)) (16:16) /= " ") then
            txtatm(j) (16:16) = txtatm(ibonds(i, j)) (16:16)
            m = m + 1
          end if
        end do
        if (m /= 1) then
          txtatm (j) (16:16) = " "
        end if
        if (txtatm(j)(18:20) == "TRP" .and. letter == "H") then
          txtatm (j) (16:16) = "2"
        end if
      else
        m = 0
!
!  Special case: Tryptophan zeta starts with 2, not 1
!
        if (txtatm(j)(18:20) == "TRP" .and. letter == "Z") m = 1
        do i = 1, jj
          j = jtype(i)
          if (nat(jtype(i)) /= 1) then
            m = m + 1
            txtatm(j) (15:15) = letter
            if (txtatm(j)(18:20) /= "UNK" .or. txtatm (j) (16:16) == " " &
              .or. txtatm (j) (16:16) > Char (min(9,m)+Ichar ("0"))) then
              txtatm (j) (16:16) = Char (min(9,m)+Ichar ("0"))
            end if
          end if
        end do
      end if
    end do
end subroutine txtype
