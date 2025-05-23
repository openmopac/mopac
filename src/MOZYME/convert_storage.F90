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
