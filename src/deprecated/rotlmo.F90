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

subroutine rotlmo (rotvec)
   !***********************************************************************
   !                                                                      *
   !   ROTLMO rotates the atomic orbitals in all LMOs.  The rotation
   !   matrix used is ROTVEC.
   !                                                                      *
   !***********************************************************************
    use MOZYME_C, only: iorbs, ncf, nce, ncocc, ncvir, nncf, nnce, icocc, &
     icvir, cocc, cvir
    use molkst_C, only: nelecs, norbs
    implicit none
    double precision, dimension (3, 3), intent (in) :: rotvec
    integer :: i, j, jj, k, ka, loop, nj, nocc, nvir
    double precision :: sum
    double precision, dimension (3) :: vec
   !
   !   Rotate the occupied LMOs
   !
    nocc = nelecs / 2
    do i = 1, nocc
      loop = ncocc(i) + 1
      do jj = nncf(i) + 1, nncf(i) + ncf(i)
        nj = iorbs(icocc(jj))
        if (nj == 4) then
          ka = loop
          do j = 1, 3
            sum = 0.0d00
            do k = 1, 3
              sum = sum + cocc(ka+k) * rotvec(k, j)
            end do
            vec(j) = sum
          end do
          do j = 1, 3
            cocc(ka+j) = vec(j)
          end do
        end if
        loop = loop + nj
      end do
    end do
   !
   !   Rotate the virtual LMOs
   !
    nvir = norbs - nocc
    do i = 1, nvir
      loop = ncvir(i) + 1
      do jj = nnce(i) + 1, nnce(i) + nce(i)
        nj = iorbs(icvir(jj))
        if (nj == 4) then
          ka = loop
          do j = 1, 3
            sum = 0.0d00
            do k = 1, 3
              sum = sum + cvir(ka+k) * rotvec(k, j)
            end do
            vec(j) = sum
          end do
          do j = 1, 3
            cvir(ka+j) = vec(j)
          end do
        end if
        loop = loop + nj
      end do
    end do
end subroutine rotlmo
