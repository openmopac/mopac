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

subroutine outer1 (ni_loc, nj, c1, c2, w, kr, e1b, e2a, enuc, mode, direct)
    use molkst_C, only: l1u, l2u, l3u, l_feather
    use parameters_C, only: natorb, tore, am
    use funcon_C, only: a0, ev
    use common_arrays_C, only: tvec
    implicit none
    logical, intent (in) :: direct
    integer, intent (in) :: mode, ni_loc, nj
    integer, intent (inout) :: kr
    double precision, intent (out) :: enuc
    double precision, dimension (*), intent (inout) :: w
    double precision, dimension (3), intent (in) :: c1, c2
    double precision, dimension (45), intent (out) :: e1b, e2a
!
    integer :: i, j, k, l
    double precision :: aee, r, rij, ww, point, const
    double precision, dimension (3) :: c2s
    double precision, external :: trunk
    if (mode == 0) then
      !
      !                                            MOLECULAR TERM
      !
     e1b = 0
     e2a = 0
      !
      rij = Sqrt ((c1(1)-c2(1))**2+ (c1(2)-c2(2))**2+ (c1(3)-c2(3))**2)
      r = rij / a0
      !
      aee = 0.5d0 / am(ni_loc) + 0.5d0 / am(nj)
      aee = aee * aee
      ww = ev / Sqrt (r*r+aee)
      if (l_feather) then
        call to_point(rij, point, const)
        ww = ww*const + (1.d0 - const)*point
      end if
      e1b(1) = -ww * tore(nj)
      e1b(3) = -ww * tore(nj)
      e1b(6) = -ww * tore(nj)
      e1b(10) = -ww * tore(nj)
      e1b(15) = -ww * tore(nj)
      e1b(21) = -ww * tore(nj)
      e1b(28) = -ww * tore(nj)
      e1b(36) = -ww * tore(nj)
      e1b(45) = -ww * tore(nj)
      e2a(1) = -ww * tore(ni_loc)
      e2a(3) = -ww * tore(ni_loc)
      e2a(6) = -ww * tore(ni_loc)
      e2a(10) = -ww * tore(ni_loc)
      e2a(15) = -ww * tore(ni_loc)
      e2a(21) = -ww * tore(ni_loc)
      e2a(28) = -ww * tore(ni_loc)
      e2a(36) = -ww * tore(ni_loc)
      e2a(45) = -ww * tore(ni_loc)
      enuc = tore(ni_loc) * tore(nj) * ww
      !
      if ( .not. direct) then
        w(1) = ww
        if (natorb(ni_loc)*natorb(nj) /= 0) then
          kr = kr + 1
        end if
      end if
    else
      !
      !                                            SOLID-STATE TERM
      !
      w(1) = 0.d0
      do i = 1, 10
        e1b(i) = 0.d0
        e2a(i) = 0.d0
      end do
      do i = -l1u, l1u
        do j = -l2u, l2u
          do k = -l3u, l3u
            do l = 1, 3
              c2s(l) = c2(l) + tvec(l, 1) * i + tvec(l, 2) * j + tvec(l, &
             & 3) * k - c1(l)
            end do
            r = Sqrt (c2s(1)**2+c2s(2)**2+c2s(3)**2)
            r = trunk (r)
               !
            w(1) = w(1) + a0 * ev / r
          end do
        end do
      end do
      e1b(1) = -w(1) * tore(nj)
      e1b(3) = -w(1) * tore(nj)
      e1b(6) = -w(1) * tore(nj)
      e1b(10) = -w(1) * tore(nj)
      e1b(15) = -w(1) * tore(nj)
      e1b(21) = -w(1) * tore(nj)
      e1b(26) = -w(1) * tore(nj)
      e1b(36) = -w(1) * tore(nj)
      e1b(45) = -w(1) * tore(nj)
      e2a(1) = -w(1) * tore(ni_loc)
      e2a(3) = -w(1) * tore(ni_loc)
      e2a(6) = -w(1) * tore(ni_loc)
      e2a(10) = -w(1) * tore(ni_loc)
      e2a(15) = -w(1) * tore(ni_loc)
      e2a(21) = -w(1) * tore(ni_loc)
      e2a(26) = -w(1) * tore(ni_loc)
      e2a(36) = -w(1) * tore(ni_loc)
      e2a(45) = -w(1) * tore(ni_loc)
      enuc = tore(ni_loc) * tore(nj) * w(1)
      if (natorb(ni_loc)*natorb(nj) /= 0) then
        kr = kr + 1
      end if
    end if
!
end subroutine outer1
