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

subroutine outer2 (ni, nj, xi, xj, w, kr, e1b, e2a, enuc, mode, direct)
    use molkst_C, only: numcal, l1u, l2u, l3u, clower
    use parameters_C, only: natorb, tore
    use common_arrays_C, only: tvec
    implicit none
    logical, intent (in) :: direct
    integer, intent (in) :: mode, ni, nj
    integer, intent (inout) :: kr
    double precision, intent (out) :: enuc
    double precision, dimension (3), intent (in) :: xi, xj
    double precision, dimension (7), intent (inout) :: w
    double precision, dimension (45), intent (out) :: e1b, e2a
!
    logical :: si, sj
    integer :: icalcn = 0, i, j, k, ki, l
    double precision :: a, rij, rijx, w1, w2, w3, w4, w5, w6, w7, cutof2, enubit
    double precision, save :: gab
    double precision, dimension (3) :: x
    double precision, dimension (7) :: wb
    double precision, dimension (22) :: ri
    double precision, dimension(45) :: e1bits, e2bits
    save icalcn, cutof2
    if (mode == 0) then
!
!   Simple di-atomic
!
      x(1) = xi(1) - xj(1)
      x(2) = xi(2) - xj(2)
      x(3) = xi(3) - xj(3)
      rij = x(1) * x(1) + x(2) * x(2) + x(3) * x(3)
      rijx = Sqrt (rij)
      rij = rijx
!
! *** COMPUTE INTEGRALS IN DIATOMIC FRAME
!
!   The last argument should be QS if the Tomasi model is
!   supported
!
      call reppd (ni, nj, rij, ri, gab)
      a = 1.d0 / rijx
      x(1) = x(1) * a
      x(2) = x(2) * a
      x(3) = x(3) * a
      if (Abs (x(3)) > 0.99999999d0) then
        x(3) = Sign (1.d0, x(3))
      end if
      si = (natorb(ni) > 1)
      sj = (natorb(nj) > 1)
      w1 = ri(1)
      ki = 1
        e1b = 0.d0
        e2a = 0.d0
      e1b(1) = -w1 * tore(nj)
      e2a(1) = -w1 * tore(ni)
!
      if ( .not. direct) then
        w(1) = w1
      end if
      if (sj) then
        w2 = -ri(5) * x(1)
        w3 = -ri(5) * x(2)
        w4 = -ri(5) * x(3)
        e2a(2) = -w2 * tore(ni)
        e2a(3) = -w1 * tore(ni)
        e2a(4) = -w3 * tore(ni)
        e2a(6) = -w1 * tore(ni)
        e2a(7) = -w4 * tore(ni)
        e2a(10) = -w1 * tore(ni)
        if (natorb(nj) > 4) then
          e2a(15) = e2a(1)
          e2a(21) = e2a(1)
          e2a(28) = e2a(1)
          e2a(36) = e2a(1)
          e2a(45) = e2a(1)
        end if
        ki = 4
        if ( .not. direct) then
          w(2) = w2
          w(3) = w3
          w(4) = w4
        end if
      end if
!
      if (si) then
        if (sj) then
          w5 = -ri(2) * x(1)
          w6 = -ri(2) * x(2)
          w7 = -ri(2) * x(3)
          ki = 7
          e1b(2) = -w5 * tore(nj)
          e1b(3) = -w1 * tore(nj)
          e1b(4) = -w6 * tore(nj)
          e1b(6) = -w1 * tore(nj)
          e1b(7) = -w7 * tore(nj)
          e1b(10) = -w1 * tore(nj)
          if (natorb(ni) > 4) then
            e1b(15) = e1b(1)
            e1b(21) = e1b(1)
            e1b(28) = e1b(1)
            e1b(36) = e1b(1)
            e1b(45) = e1b(1)
          end if
!
          if ( .not. direct) then
            w(5) = w5
            w(6) = w6
            w(7) = w7
          end if
        else
          w2 = -ri(2) * x(1)
          w3 = -ri(2) * x(2)
          w4 = -ri(2) * x(3)
          ki = 4
          e1b(2) = -w2 * tore(nj)
          e1b(3) = -w1 * tore(nj)
          e1b(4) = -w3 * tore(nj)
          e1b(6) = -w1 * tore(nj)
          e1b(7) = -w4 * tore(nj)
          e1b(10) = -w1 * tore(nj)
          if (natorb(ni) > 4) then
            e1b(15) = e1b(1)
            e1b(21) = e1b(1)
            e1b(28) = e1b(1)
            e1b(36) = e1b(1)
            e1b(45) = e1b(1)
          end if
!
          if ( .not. direct) then
            w(2) = w2
            w(3) = w3
            w(4) = w4
          end if
        end if
      end if
      enuc = tore(ni) * tore(nj) * w1
    else
!
!   SOLID SYSTEM
!
      if (icalcn /= numcal) then
        icalcn = numcal
        cutof2 = clower**2
      end if
      do i = 1, 10
        e1b(i) = 0.d0
        e2a(i) = 0.d0
      end do
      si = (natorb(ni) > 1)
      sj = (natorb(nj) > 1)
      ki = 1
      if (si .or. sj) then
        ki = 4
      end if
      if (si .and. sj) then
        ki = 7
      end if
      do i = 1, ki
        w(i) = 0.d0
      end do
      do i = -l1u, l1u
        do j = -l2u, l2u
          do k = -l3u, l3u
            do l = 1, 3
              x(l) = xi(l) + tvec(l, 1) * i + tvec(l, 2) * j + tvec(l, 3) &
             & * k - xj(l)
            end do
            rij = x(1) * x(1) + x(2) * x(2) + x(3) * x(3)
            if (rij >  cutof2) then
!
!  Interaction distance is greater than cutoff, so use point-charge
!
              rij = sqrt(rij)
              call point(rij, ni, nj, ri, l, e1bits, e2bits, enubit)
              ri(4) = ri(1)
              ri(6) = 0.d0
            else
              rij = Sqrt (rij)
!
! *** COMPUTE INTEGRALS IN DIATOMIC FRAME
!
              call reppd (ni, nj, rij, ri, gab)
            end if
!
              a = 1.d0 / rij
              x(1) = x(1) * a
              x(2) = x(2) * a
              x(3) = x(3) * a
              if (Abs (x(3)) > 0.99999999d0) x(3) = Sign (1.d0, x(3))
              wb(1) = ri(1)
              e1b(1) = e1b(1) - wb(1) * tore(nj)
              e2a(1) = e2a(1) - wb(1) * tore(ni)
              if (sj) then
                wb(2) = -ri(5) * x(1)
                wb(3) = -ri(5) * x(2)
                wb(4) = -ri(5) * x(3)
                e2a(2) = e2a(2) - wb(2) * tore(ni)
                e2a(3) = e2a(3) - wb(1) * tore(ni)
                e2a(4) = e2a(4) - wb(3) * tore(ni)
                e2a(6) = e2a(6) - wb(1) * tore(ni)
                e2a(7) = e2a(7) - wb(4) * tore(ni)
                e2a(10) = e2a(10) - wb(1) * tore(ni)
              end if
!
              if (si) then
                if (sj) then
                  wb(5) = -ri(2) * x(1)
                  wb(6) = -ri(2) * x(2)
                  wb(7) = -ri(2) * x(3)
                  e1b(2) = e1b(2) - wb(5) * tore(nj)
                  e1b(3) = e1b(3) - wb(1) * tore(nj)
                  e1b(4) = e1b(4) - wb(6) * tore(nj)
                  e1b(6) = e1b(6) - wb(1) * tore(nj)
                  e1b(7) = e1b(7) - wb(7) * tore(nj)
                  e1b(10) = e1b(10) - wb(1) * tore(nj)
                else
                  wb(2) = -ri(2) * x(1)
                  wb(3) = -ri(2) * x(2)
                  wb(4) = -ri(2) * x(3)
                  e1b(2) = e1b(2) - wb(2) * tore(nj)
                  e1b(3) = e1b(3) - wb(1) * tore(nj)
                  e1b(4) = e1b(4) - wb(3) * tore(nj)
                  e1b(6) = e1b(6) - wb(1) * tore(nj)
                  e1b(7) = e1b(7) - wb(4) * tore(nj)
                  e1b(10) = e1b(10) - wb(1) * tore(nj)
              end if
            end if
            do l = 1, ki
              w(l) = w(l) + wb(l)
            end do
          end do
        end do
      end do
      enuc = tore(ni) * tore(nj) * w(1)
    end if
!
    if ( .not. direct) then
      if (natorb(ni)*natorb(nj) == 0) then
        ki = 0
      end if
      kr = kr + ki
    end if
end subroutine outer2
