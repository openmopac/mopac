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

subroutine hbonds (fao, nocc, nvir, iused, nij_loc, cutoff)
   !**********************************************************************
   !                                                                     *
   !   HBONDS IDENTIFIES HYDROGEN BONDS WHICH HAVE NOT BEEN CREATED IN   *
   !          EARLIER STEPS                                              *
   !                                                                     *
   !**********************************************************************
   !                                                                     *
   !  FAO:      Fock matrix over atomic orbitals                         *
   !  FMO:      Fock matrix over molecular orbitals                      *
   !  VECTOR:   Eigenvectors of M.O.s                                    *
   !  FILLED:   Occupied M.O.s                                           *
   !  EMPTY:    Virtual M.O.s                                            *
   !  NOCC:     Number of occupied M.O.s                                 *
   !  EIG:      All the eigenvalues                                      *
   !  N:        Number of M.O.s = number of A.O.s*                       *
   !                                                                     *
   !  On exit:  NIJ = NUMBER OF HYDROGEN BONDS WHICH NEED TO BE CREATED  *
   !                                                                     *
   !**********************************************************************
    use molkst_C, only: numat, norbs, mpack
    use MOZYME_C, only : nncf, ncf, nnce, nce, iorbs, icocc, icvir, &
    fmo, ifmo
    use common_arrays_C, only:  p
    implicit none
    integer, intent (out) :: nij_loc
    integer, intent (in) :: nocc, nvir
    double precision, intent (in) :: cutoff
    integer, dimension (norbs), intent (out) :: iused
    double precision, dimension (mpack), intent (in) :: fao
    integer :: i
    integer :: j, k, kk, l, ll, m
    double precision :: sumf, sump
    integer, external :: ijbo
    intrinsic Min
   !
   !   CHECK FOCK MATRIX AND DENSITY MATRIX TO ENSURE THAT ALL
   !   FOCK MATRIX ELEMENTS ARE BEING CORRECTLY HANDLED BY THE
   !   DIAGONALISER.
   !
    nij_loc = 0
    do i = 1, numat
      do j = 1, i - 1
        if (ijbo(i, j) >= 0) then
          ll = ijbo (i, j)
          sumf = 0.d0
          sump = 0.d0
          do k = 1, iorbs(i)
            do l = 1, iorbs(j)
              ll = ll + 1
              sumf = sumf + fao(ll) ** 2
              sump = sump + p(ll) ** 2
            end do
          end do
          if (sump < 1.d-10 .and. sumf > cutoff) then
            m = 0
            do k = 1, nocc
              !
              !  If the LMO only contains one atom, do NOT check the
              !  second atom.
              !
              kk = Min (ncf(k)-1, 1)
              l = nncf(k) + 1
              if (icocc(l) == i .or. icocc(l+kk) == i .or. icocc(l) == j .or. &
             & icocc(l+kk) == j) then
                m = m + 1
                iused(m) = k
              end if
            end do
            do k = 1, nvir
              !
              !  If the LMO only contains one atom, do NOT check the
              !  second atom.
              !
              kk = Min (nce(k)-1, 1)
              l = nnce(k) + 1
              if (icvir(l) == i .or. icvir(l+kk) == i .or. icvir(l) == j .or. &
             & icvir(l+kk) == j) then
                do l = 1, m
                  nij_loc = nij_loc + 1
                  fmo(nij_loc) = 0.1d0
                  ifmo(2, nij_loc) = iused(l)
                  ifmo(1, nij_loc) = k
                end do
              end if
            end do
          end if
        end if
      end do
    end do
end subroutine hbonds
