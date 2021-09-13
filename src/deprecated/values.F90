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

subroutine values (type)
    use molkst_C, only: norbs, mpack
    use MOZYME_C, only : icocc_dim, cocc_dim, icvir_dim, cvir_dim, &
    isort, nncf, ncf, icocc, ncocc, iorbs, cocc, nnce, nce, icvir, ncvir, &
    cvir, noccupied, nvirtual
    use chanel_C, only: iw
    use common_arrays_C, only: f, eigs
    implicit none
    character (len=*), intent (in) :: type
    integer :: alloc_stat
    double precision, dimension(:), allocatable :: fdiat
    if (.not. Allocated (isort)) then
      allocate (isort(norbs), stat=alloc_stat)
      if (alloc_stat /= 0) then
        call memory_error ("values")
        go to 1100
      end if
    end if
    allocate (fdiat(norbs), stat=alloc_stat)
    if (alloc_stat /= 0) then
      call memory_error ("values")
      go to 1100
    end if
!
!
!
    if (type == "OCCUPIED") then
      call valuen (f, mpack, noccupied, nncf, ncf, norbs, icocc, icocc_dim, &
        ncocc, iorbs, cocc, cocc_dim, fdiat, eigs, isort)
!
!
!
    else if (type == "VIRTUAL") then
      !
      call valuen (f, mpack, nvirtual, nnce, nce, norbs, icvir, icvir_dim, ncvir, &
        iorbs, cvir, cvir_dim, fdiat, eigs(noccupied+1:), isort(noccupied+1:))
    else
      write (iw,*) " Error"
      call mopend ("Error")
    end if
    deallocate (fdiat)
1100 continue
end subroutine values
subroutine valuen (fao, nfao, nocc, nncf, ncf, nnn, icocc, nico, ncocc, iorbs, &
   cocc, nco, fdiat, eigf, isort)
    implicit none
    integer, intent (in) :: nco, nfao, nico, nnn, nocc
    integer, dimension (nico), intent (in) :: icocc
    integer, dimension (*), intent (in) :: ncf, ncocc, nncf
    integer, dimension (nocc), intent (out) :: isort
    integer, dimension (*), intent (in) :: iorbs
    double precision, dimension (nco), intent (in) :: cocc
    double precision, dimension (nfao), intent (in) :: fao
    double precision, dimension (nnn), intent (out) :: fdiat
    double precision, dimension (nocc), intent (inout) :: eigf
    integer :: i, i4, ii, j, j1, j4, jl, jx, k, k1, kl, l, loopi
    double precision :: sum, sum1
    integer, external :: ijbo
   !**********************************************************************
   !
   !   Calculates the energy levels of the LMOs, and sorts the LMOs into
   !   increasing energy order.
   !
   !**********************************************************************
    do i = 1, nocc
      loopi = ncocc(i)
      l = 0
      do j = nncf(i) + 1, nncf(i) + ncf(i)
        j1 = icocc(j)
        sum = 0.d0
        do k = l + 1, l + iorbs(j1)
          sum = sum + cocc(k+loopi) ** 2
        end do
        l = l + iorbs(j1)
      end do
      sum = 0.d0
      jl = loopi
      do j = nncf(i) + 1, nncf(i) + ncf(i)
        j1 = icocc(j)
        do jx = 1, iorbs(j1)
          jl = jl + 1
          kl = loopi
          do k = nncf(i) + 1, nncf(i) + ncf(i)
            k1 = icocc(k)
            if (ijbo(k1, j1) >= 0) then
                  !
                  !  Extract the atom-atom intersection of FAO
                  !
              if (k1 > j1) then
                     !
                     !   LOWER TRIANGLE
                     !
                ii = ijbo (k1, j1) + jx - iorbs(j1)
                do i4 = 1, iorbs(k1)
                  ii = ii + iorbs(j1)
                  fdiat(i4) = fao(ii)
                end do
              else if (k1 < j1) then
                     !
                     !   UPPER TRIANGLE
                     !
                ii = ijbo (k1, j1) + iorbs(k1) * (jx-1)
                do i4 = 1, iorbs(k1)
                  ii = ii + 1
                  fdiat(i4) = fao(ii)
                end do
              else
                     !
                     !   DIAGONAL TERM
                     !
                ii = ijbo (k1, j1) + (jx*(jx+1)) / 2
                do i4 = jx + 1, iorbs(k1)
                  ii = ii + jx
                  fdiat(i4) = fao(ii)
                  ii = ii + i4 - jx
                end do
                ii = ijbo (k1, j1) + (jx*(jx-1)) / 2
                do j4 = 1, jx
                  ii = ii + 1
                  fdiat(j4) = fao(ii)
                end do
              end if
              sum1 = 0.d0
              do i4 = 1, iorbs(k1)
                kl = kl + 1
                sum1 = sum1 + fdiat(i4) * cocc(kl)
              end do
              sum = sum + cocc(jl) * sum1
            end if
          end do
        end do
      end do
      eigf(i) = sum
    end do
   !
   !  Now sort eigenvalues.
   !
    do i = 1, nocc
      fdiat(i) = eigf(i)
    end do
    do i = 1, nocc
      k = 0
      sum = 1.d9
      do j = 1, nocc
        if (fdiat(j) < sum) then
          k = j
          sum = fdiat(j)
        end if
      end do
      fdiat(k) = 1.d10
      isort(i) = k
    end do
end subroutine valuen
