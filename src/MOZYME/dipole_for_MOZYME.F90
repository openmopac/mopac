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

double precision function dipole_for_MOZYME (dipvec, mode)
    use molkst_C, only: numcal, numat, keywrd
    use chanel_C, only: iw
    use MOZYME_C, only : iorbs
    use common_arrays_C, only : nat, q, p, coord
    use parameters_C, only: dd, ams, ddp
    use funcon_C, only: fpc_8, fpc_1, a0
    implicit none
    double precision, intent(out) :: dipvec(3)
    integer, intent (in) :: mode
!.
    logical, save :: chargd, first, force
    integer, save :: icalcn = 0
    integer :: i, j, k, l, ni
    double precision :: const, sum, dx, dy, dz, hyfpd, xt
    double precision, save :: wtmol_loc
    double precision, dimension (3) :: center
    double precision, dimension (107), save :: hyf
    double precision, dimension (4, 3), save :: dip
    integer, external :: ijbo
    data hyf(1) / 0.0d00 /
    first = (icalcn /= numcal)
    icalcn = numcal
   !
   !   CONST = c.C.a0.2
   !
    const = fpc_8 * fpc_1 * a0 * 2.d0 * 1.d-10
    if (first) then
      hyf(2:107) = const * dd(2:107)
      wtmol_loc = 0.d0
      sum = 0.d0
      do i = 1, numat
        wtmol_loc = wtmol_loc + ams(nat(i))
        sum = sum + q(i)
      end do
      chargd = (Abs (sum) > 0.5d0)
      force = index(keywrd,'FORCE') + index(keywrd," THERMO") + index(keywrd,'IRC') /= 0
    end if
    if (chargd) then
      !
      !   NEED TO RESET ION'S POSITION SO THAT THE CENTER OF MASS IS AT THE
      !   ORIGIN.
      !
      center = 0.d0
      do i = 1, 3
        do j = 1, numat
          center(i) = center(i) + ams(nat(j)) * coord(i, j)
        end do
      end do
      center = center / wtmol_loc
      do i = 1, 3
        do j = 1, numat
          coord(i, j) = coord(i, j) - center(i)
        end do
      end do
    end if
    dip = 0.0d00
    do i = 1, numat
      ni = nat(i)
      l = iorbs(i) - 1
      if (l >= 3) then
        do j = 1, 3
          k = ijbo (i, i) + 1 + (j*(j+1)) / 2
          dip(j, 2) = dip(j, 2) - hyf(ni) * p(k)
        end do
        if (l == 8) then
          xt = 1.d0 / Sqrt (3.d00)
          hyfpd = 2.0d0 * ddp(5, ni) * a0 * fpc_8 * fpc_1 * 1.d-10
          k = ijbo (i, i)
          !
          ! x:   <x|x|x2-y2> + <z|x|xz> - 1/root(3)<x|x|z2> + <y|x|xy>
          ! y: - <y|y|x2-y2> - 1/root(3)<y|y|z2> + <z|y|yz> + <x|y|xy>
          ! z:   <x|z|xz> + 2/root(3)<z|z|z2> + <y|z|yz>
          !
          dx = p(k+12) + p(k+19) - p(k+23) * xt + p(k+39)
          dy = -p(k+13) - p(k+24) * xt + p(k+32) + p(k+38)
          dz = p(k+17) + 2.d0 * p(k+25) * xt + p(k+31)
          dip(1, 2) = dip(1, 2) - dx * hyfpd
          dip(2, 2) = dip(2, 2) - dy * hyfpd
          dip(3, 2) = dip(3, 2) - dz * hyfpd
        end if
      end if
      do j = 1, 3
         !
         !  FPC(8)=SPEED OF LIGHT,  FPC(1)=CHARGE ON ELECTRON.
         !
        dip(j, 1) = dip(j, 1) + fpc_8 * fpc_1 * 1.d-10 * q(i) * coord(j, i)
      !
      !  FPC(8)=SPEED OF LIGHT,  FPC(1)=CHARGE ON ELECTRON.
      !
      end do
    end do
    do j = 1, 3
      dip(j, 3) = dip(j, 2) + dip(j, 1)
    end do
    dip(4, :) = Sqrt (dip(1, :)**2+dip(2, :)**2+dip(3, :)**2)
    if (force) then
      dipvec(1) = dip(1,3)
      dipvec(2) = dip(2,3)
      dipvec(3) = dip(3,3)
    end if
    dipole_for_MOZYME = dip(4, 3)
    if ( chargd) then
      do i = 1, 3
        coord(i, :numat) = coord(i, :numat) + center(i)
      end do
    end if
    if (mode == 0) return
   !
10000 format (" DIPOLE", 11x, "X         Y         Z       TOTAL", / &
   & , " POINT-CHG.", 4f10.3/, " HYBRID", 4x, 4f10.3/, " SUM", 7 x, 4 f10.3)
    write (iw, 10000) ((dip(i, j), i=1, 4), j=1, 3)
!
end function dipole_for_MOZYME
