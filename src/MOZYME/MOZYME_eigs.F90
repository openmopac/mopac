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

subroutine MOZYME_eigs (nocc)
   !**********************************************************************
   !
   !  f:      FOCK MATRIX OVER ATOMIC ORBITALS
   !  VECTOR:   EIGENVECTORS OF M.O.S
   !  FILLED:   OCCUPIED M.O.S
   !  eigs:     EIGENVALUES OF OCCUPIED M.O.S
   !  NOCC:     NUMBER OF OCCUPIED M.O.S
   !  EIG:      ALL THE EIGENVALUES
   !  N:        NUMBER OF M.O.S = NUMBER OF A.O.S*
   !
   !**********************************************************************
   !
   !  ALL OCCUPIED LMO ENERGY LEVELS ARE CALCULATED.
   !
   !**********************************************************************
   !
    use molkst_C, only: norbs
    use MOZYME_C, only : icocc_dim, &
       lijbo, nijbo, &
       ijc, ncf, nncf, ncocc, &
     & iorbs, icocc, cocc
    use common_arrays_C, only : eigs, nfirst, nlast, p, f
    implicit none
    integer, intent (in) :: nocc
    integer :: i, i4, ii, j, j1, j4, &
         & jl, jx, k, k1, kk, kl, l, loopi, loopj, k1j1
    double precision, allocatable :: avir(:), aocc(:)
    double precision :: cutoff, sum, sum1
    integer, external :: ijbo
    allocate (avir(norbs), aocc(icocc_dim))
    !
    !
    aocc(:) = 0.d0

    do j = 1, nocc
        loopj = ncocc(j)
        kl = 0
        do kk = nncf(j) + 1, nncf(j) + ncf(j)
          k1 = icocc(kk)
          sum = 0.d0
          do k = nfirst(k1), nlast(k1)
            kl = kl + 1
            sum = sum + cocc(kl+loopj) ** 2
          end do
          !
          !   AOCC(KK) HOLDS THE SQUARE OF THE CONTRIBUTION OF THE KK'TH ATOM
          !   IN THE OCCUPIED SET.  NOTE:  THIS IS NOT ATOM KK.
          !
          aocc(kk) = sum
          !
        end do
      !
    end do
    cutoff = 1.d-8
    ijc = 0
!
!  EVALUATE THE OCCUPIED ENERGY LEVELS
!
!
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
            !
            !   AOCC(KK) HOLDS THE SQUARE OF THE CONTRIBUTION OF THE KK'TH ATOM
            !   IN THE OCCUPIED SET.  NOTE:  THIS IS NOT ATOM KK.
            !
            !   AVIR(J1) HOLDS THE SQUARE OF THE CONTRIBUTION OF ATOM J1
            !   IN THE OCCUPIED SET.  NOTE:  THIS IS NOT THE J1'th ATOM.
            !
            avir(j1) = sum
            !
            !   AOCC(KK) HOLDS THE SQUARE OF THE CONTRIBUTION OF THE KK'TH ATOM
            !   IN THE OCCUPIED SET.  NOTE:  THIS IS NOT ATOM KK.
            !
            !   AVIR(J1) HOLDS THE SQUARE OF THE CONTRIBUTION OF ATOM J1
            !   IN THE OCCUPIED SET.  NOTE:  THIS IS NOT THE J1'th ATOM.
            !
          end do
          sum = 0.d0
          jl = loopi
          !
          if (lijbo) then
            do j = nncf(i) + 1, nncf(i) + ncf(i)
              j1 = icocc(j)
              kl = loopi
              do k = nncf(i) + 1, nncf(i) + ncf(i)
                k1 = icocc(k)
                k1j1 = nijbo(k1, j1)
                if (k1j1 >= 0) then
                  if (avir(k1)*p(k1j1+1)*aocc(k) >= cutoff) then
                    !
                    !  EXTRACT THE ATOM-ATOM INTERSECTION OF f
                    !
                    if (k1 > j1) then
                      !
                      !   LOWER TRIANGLE
                      !
                      do jx = 1, iorbs(j1)
                        sum1 = 0.d0
                        do i4 = 1, iorbs(k1)
                          ii = k1j1 + (i4-1) * iorbs(j1) + jx
                          sum1 = sum1 + f(ii) * cocc(kl+i4)
                        end do
                        sum = sum + cocc(jl+jx) * sum1
                      end do
                    else if (k1 < j1) then
                      !
                      !   UPPER TRIANGLE
                      !
                      do jx = 1, iorbs(j1)
                        sum1 = 0.d0
                        do i4 = 1, iorbs(k1)
                          ii = k1j1 + (jx-1) * iorbs(k1) + i4
                          sum1 = sum1 + f(ii) * cocc(kl+i4)
                        end do
                        sum = sum + cocc(jl+jx) * sum1
                      end do
                    else
                      !
                      !   DIAGONAL TERM
                      !
                      do jx = 1, iorbs(j1)
                        sum1 = 0.d0
                        do j4 = 1, jx
                          ii = k1j1 + (jx*(jx-1)) / 2 + j4
                          sum1 = sum1 + f(ii) * cocc(kl+j4)
                        end do
                        ii = k1j1 + (jx*(jx+1)) / 2
                        do i4 = jx + 1, iorbs(k1)
                          ii = k1j1 + (i4*(i4-1)) / 2 + jx
                          sum1 = sum1 + f(ii) * cocc(kl+i4)
                        end do
                        sum = sum + cocc(jl+jx) * sum1
                      end do
                    end if
                  end if
                end if
                kl = kl + iorbs(k1)
              end do
              !
              jl = jl + iorbs(j1)
            end do
          else
            do j = nncf(i) + 1, nncf(i) + ncf(i)
              j1 = icocc(j)
              kl = loopi
              do k = nncf(i) + 1, nncf(i) + ncf(i)
                k1 = icocc(k)
                k1j1 = ijbo (k1, j1)
                if (k1j1 >= 0) then
                  if (avir(k1)*p(k1j1+1)*aocc(k) >= cutoff) then
                    !
                    !  EXTRACT THE ATOM-ATOM INTERSECTION OF f
                    !
                    if (k1 > j1) then
                      !
                      !   LOWER TRIANGLE
                      !
                      do jx = 1, iorbs(j1)
                        sum1 = 0.d0
                        do i4 = 1, iorbs(k1)
                          ii = k1j1 + (i4-1) * iorbs(j1) + jx
                          sum1 = sum1 + f(ii) * cocc(kl+i4)
                        end do
                        sum = sum + cocc(jl+jx) * sum1
                      end do
                    else if (k1 < j1) then
                      !
                      !   UPPER TRIANGLE
                      !
                      do jx = 1, iorbs(j1)
                        sum1 = 0.d0
                        do i4 = 1, iorbs(k1)
                          ii = k1j1 + (jx-1) * iorbs(k1) + i4
                          sum1 = sum1 + f(ii) * cocc(kl+i4)
                        end do
                        sum = sum + cocc(jl+jx) * sum1
                      end do
                    else
                      !
                      !   DIAGONAL TERM
                      !
                      do jx = 1, iorbs(j1)
                        sum1 = 0.d0
                        do j4 = 1, jx
                          ii = k1j1 + (jx*(jx-1)) / 2 + j4
                          sum1 = sum1 + f(ii) * cocc(kl+j4)
                        end do
                        ii = k1j1 + (jx*(jx+1)) / 2
                        do i4 = jx + 1, iorbs(k1)
                          ii = k1j1 + (i4*(i4-1)) / 2 + jx
                          sum1 = sum1 + f(ii) * cocc(kl+i4)
                        end do
                        sum = sum + cocc(jl+jx) * sum1
                      end do
                    end if
                  end if
                end if
                kl = kl + iorbs(k1)
              end do
              !
              jl = jl + iorbs(j1)
            end do
          end if
          !
          eigs(i) = sum
          !
      end do

end subroutine MOZYME_eigs
