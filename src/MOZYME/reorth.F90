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

subroutine reorth (ws)
    use molkst_C, only: numat, norbs
    use MOZYME_C, only : nvirtual, noccupied, icocc_dim, &
       & icvir_dim, cocc_dim, cvir_dim, ncf, nce, nncf, nnce, ncocc, ncvir, iorbs, &
       & icocc, icvir, cocc, cvir
    use common_arrays_C, only : nfirst
    implicit none
    double precision, dimension (norbs) :: ws
!
    integer :: i, ii, j, j1, jj, jx, loopi, loopii, alloc_stat
    double precision :: sum, sumtot
    logical, dimension (:), allocatable :: latom_loc
    integer, dimension (:), allocatable :: iused
!
    allocate (latom_loc(numat), iused(numat), stat=alloc_stat)
    if (alloc_stat /= 0) then
      call memory_error ("reorth")
      go to 1100
    end if
    sumtot = 0.d0
    do i = 1, nvirtual
      loopi = ncvir(i)
      do j = 1, numat
        latom_loc(j) = .false.
      end do
      do jj = nnce(i) + 1, nnce(i) + nce(i)
        j1 = icvir(jj)
        latom_loc(j1) = .true.
        j = nfirst(j1) - 1
        do jx = 1, iorbs(j1)
          loopi = loopi + 1
          j = j + 1
          ws(j) = cvir(loopi)
        end do
      end do
      !
      !    WS holds the LMO coefficients, in atom number order.
      !
      do ii = i + 1, nvirtual
        j1 = icvir(nnce(ii)+1)
        sum = 0.d0
        loopii = ncvir(ii)
        do jj = nnce(ii) + 1, nnce(ii) + nce(ii)
          j1 = icvir(jj)
          if (latom_loc(j1)) then
            j = nfirst(j1) - 1
            do jx = 1, iorbs(j1)
              loopii = loopii + 1
              j = j + 1
              sum = sum + ws(j) * cvir(loopii)
            end do
          else
            loopii = loopii + iorbs(j1)
          end if
        end do
        call adjvec (cvir, cvir_dim, icvir, icvir_dim, nnce, nce, nvirtual, &
             & ncvir, i, iorbs, cvir, cvir_dim, icvir, icvir_dim, nnce, &
             & nce, nvirtual, ncvir, ii, sum, iused, sumtot)
      end do
      !
      !   EVALUATE THE VIRTUAL-OCCUPIED LMO OVERLAPS
      !
      do ii = 1, noccupied
        j1 = icocc(nncf(ii)+1)
        sum = 0.d0
        loopii = ncocc(ii)
        do jj = nncf(ii) + 1, nncf(ii) + ncf(ii)
          j1 = icocc(jj)
          if (latom_loc(j1)) then
            j = nfirst(j1) - 1
            do jx = 1, iorbs(j1)
              loopii = loopii + 1
              j = j + 1
              sum = sum + ws(j) * cocc(loopii)
            end do
          else
            loopii = loopii + iorbs(j1)
          end if
        end do
        call adjvec (cvir, cvir_dim, icvir, icvir_dim, nnce, nce, nvirtual, &
             & ncvir, i, iorbs, cocc, cocc_dim, icocc, icocc_dim, nncf, &
             & ncf, noccupied, ncocc, ii, sum, iused, sumtot)
      end do
    end do
    do i = 1, noccupied
      loopi = ncocc(i)
      do j = 1, numat
        latom_loc(j) = .false.
      end do
      do jj = nncf(i) + 1, nncf(i) + ncf(i)
        j1 = icocc(jj)
        latom_loc(j1) = .true.
        j = nfirst(j1) - 1
        do jx = 1, iorbs(j1)
          loopi = loopi + 1
          j = j + 1
          ws(j) = cocc(loopi)
        end do
      end do
      !
      !    WS holds the LMO coefficients, in atom number order.
      !
      !  EVALUATE THE OCCUPIED-OCCUPIED LMO OVERLAPS
      !
      do ii = i + 1, noccupied
        j1 = icocc(nncf(ii)+1)
        sum = 0.d0
        loopii = ncocc(ii)
        do jj = nncf(ii) + 1, nncf(ii) + ncf(ii)
          j1 = icocc(jj)
          if (latom_loc(j1)) then
            j = nfirst(j1) - 1
            do jx = 1, iorbs(j1)
              loopii = loopii + 1
              j = j + 1
              sum = sum + ws(j) * cocc(loopii)
            end do
          else
            loopii = loopii + iorbs(j1)
          end if
        end do
        call adjvec (cocc, cocc_dim, icocc, icocc_dim, nncf, ncf, &
             & noccupied, ncocc, i, iorbs, cocc, cocc_dim, icocc, &
             & icocc_dim, nncf, ncf, noccupied, ncocc, ii, sum, &
			 & iused, sumtot)
      end do
    end do
    deallocate (latom_loc, iused)
1100 continue
end subroutine reorth
