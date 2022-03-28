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

subroutine convert_mat_packed_to_triangle (matrix_packed, matrix_triangle)
   !***********************************************************************
   !                                                                      *
   !        Converts matrix from the storage format used                  *
   !        in MOZYME into lower triangular form.
   !                                                                      *
   !***********************************************************************
    use molkst_C, only: numat, norbs, mpack
    use common_arrays_C, only: nfirst, nlast
    implicit none
    double precision, dimension (norbs*(norbs+1)/2), intent (out) :: &
         & matrix_triangle
    double precision, dimension (mpack), intent (in) :: matrix_packed
    integer :: i, ii, ij, il, iu, j, jj, jl, ju, linear
    integer, external :: ijbo
   !
    linear = (norbs*(norbs+1)) / 2
    do i = 1, linear
      matrix_triangle(i) = 0.d0
    end do
    do i = 1, numat
      do j = 1, i
        if (ijbo(i, j) >= 0) then
          ij = ijbo (i, j)
          il = nfirst(i)
          iu = nlast(i)
          jl = nfirst(j)
          ju = nlast(j)
          do ii = il, iu
            do jj = jl, min(ju, ii)
              ij = ij + 1
              matrix_triangle((ii*(ii-1))/2+jj) = matrix_packed(ij)
            end do
          end do
        end if
      end do
    end do
end subroutine convert_mat_packed_to_triangle

subroutine convert_lmo_packed_to_square (c_square)
   !***********************************************************************
   !                                                                      *
   !        Converts LMOs from the storage format used in MOZYME to that  *
   !        used in MOPAC.
   !                                                                      *
   !***********************************************************************
    use molkst_C, only: norbs, nelecs
    use MOZYME_C, only: isort, nce, ncf, ncocc, ncvir, nnce, nncf, &
          icocc, icvir, cocc, cvir
    use common_arrays_C, only : nfirst, nlast
    implicit none
    double precision, dimension (norbs, norbs), intent (out) :: c_square
    integer :: i, ii, iunsrt, j, jj, k, ka, nocc, nvir, alloc_stat
    nocc = nelecs / 2
    nvir = norbs - nocc
    if (.not. Allocated (isort)) then
!
!  isort has not been assigned, therefore assume that eigenvectors are not wanted.
!  in which case, use the LMOs in the order in which they occur
!
      allocate (isort(norbs), stat=alloc_stat)
      if (alloc_stat /= 0) then
        call memory_error ("convert_lmo_packed_to_square")
        call mopend ("Error in converting storage format of orbitals")
      end if
      do i = 1, nocc
        isort(i) = i
      end do
      do i = 1,nvir
        isort(nocc+i) = i
      end do
    end if
      !
      !   Put occupied LMOs into conventional form
      !
      do iunsrt = 1, nocc
        i = isort(iunsrt)
        do j = 1, norbs
          c_square(j, iunsrt) = 0.d0
        end do
        ka = ncocc(i)
        do jj = nncf(i) + 1, nncf(i) + ncf(i)
          j = icocc(jj)
          do k = nfirst(j), nlast(j)
            ka = ka + 1
            c_square(k, iunsrt) = cocc(ka)
          end do
        end do
      end do
      !
      !   Put virtual LMOs into conventional form
      !
      do iunsrt = 1, nvir
        i = isort(iunsrt+nocc)
        ii = iunsrt + nocc
        do j = 1, norbs
          c_square(j, ii) = 0.d0
        end do
        ka = ncvir(i)
        do jj = nnce(i) + 1, nnce(i) + nce(i)
          j = icvir(jj)
          do k = nfirst(j), nlast(j)
            ka = ka + 1
            c_square(k, ii) = cvir(ka)
          end do
        end do
      end do
end subroutine convert_lmo_packed_to_square
