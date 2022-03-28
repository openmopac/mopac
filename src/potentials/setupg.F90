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

      subroutine setupg
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE common_arrays_C, only : nat, nfirst, nlast
      use molkst_C, only : numat
      USE overlaps_C, only : ccc, zzz, allc, allz
      use parameters_C, only : zs, zp
!***********************************************************************
      implicit none

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ni, ia, nqn
      double precision :: xi
!-----------------------------------------------
!     SET-UP THE STEWART'S STO-6G EXPANSIONS
!                                     1S
      allz(1,1,1) = 2.310303149D01
      allz(2,1,1) = 4.235915534D00
      allz(3,1,1) = 1.185056519D00
      allz(4,1,1) = 4.070988982D-01
      allz(5,1,1) = 1.580884151D-01
      allz(6,1,1) = 6.510953954D-02
!
      allc(1,1,1) = 9.163596280D-03
      allc(2,1,1) = 4.936149294D-02
      allc(3,1,1) = 1.685383049D-01
      allc(4,1,1) = 3.705627997D-01
      allc(5,1,1) = 4.164915298D-01
      allc(6,1,1) = 1.303340841D-01
!                                     2S
      allz(1,2,1) = 2.768496241D01
      allz(2,2,1) = 5.077140627D00
      allz(3,2,1) = 1.426786050D00
      allz(4,2,1) = 2.040335729D-01
      allz(5,2,1) = 9.260298399D-02
      allz(6,2,1) = 4.416183978D-02
!
      allc(1,2,1) = -4.151277819D-03
      allc(2,2,1) = -2.067024148D-02
      allc(3,2,1) = -5.150303337D-02
      allc(4,2,1) = 3.346271174D-01
      allc(5,2,1) = 5.621061301D-01
      allc(6,2,1) = 1.712994697D-01
!                                     2P
      allz(1,2,2) = 5.868285913D00
      allz(2,2,2) = 1.530329631D00
      allz(3,2,2) = 5.475665231D-01
      allz(4,2,2) = 2.288932733D-01
      allz(5,2,2) = 1.046655969D-01
      allz(6,2,2) = 4.948220127D-02
!
      allc(1,2,2) = 7.924233646D-03
      allc(2,2,2) = 5.144104825D-02
      allc(3,2,2) = 1.898400060D-01
      allc(4,2,2) = 4.049863191D-01
      allc(5,2,2) = 4.012362861D-01
      allc(6,2,2) = 1.051855189D-01
!                                     3S
      allz(1,3,1) = 3.273031938D00
      allz(2,3,1) = 9.200611311D-01
      allz(3,3,1) = 3.593349765D-01
      allz(4,3,1) = 8.636686991D-02
      allz(5,3,1) = 4.797373812D-02
      allz(6,3,1) = 2.724741144D-02
      allc(1,3,1) = -6.775596947D-03
      allc(2,3,1) = -5.639325779D-02
      allc(3,3,1) = -1.587856086D-01
      allc(4,3,1) = 5.534527651D-01
      allc(5,3,1) = 5.015351020D-01
      allc(6,3,1) = 7.223633674D-02
!                                     3P
      allz(1,3,2) = 5.077973607D00
      allz(2,3,2) = 1.340786940D00
      allz(3,3,2) = 2.248434849D-01
      allz(4,3,2) = 1.131741848D-01
      allz(5,3,2) = 6.076408893D-02
      allz(6,3,2) = 3.315424265D-02
      allc(1,3,2) = -3.329929840D-03
      allc(2,3,2) = -1.419488340D-02
      allc(3,3,2) = 1.639395770D-01
      allc(4,3,2) = 4.485358256D-01
      allc(5,3,2) = 3.908813050D-01
      allc(6,3,2) = 7.411456232D-02
!                                     4S
      allz(1,4,1) = 1.365346D+00
      allz(2,4,1) = 4.393213D-01
      allz(3,4,1) = 1.877069D-01
      allz(4,4,1) = 9.360270D-02
      allz(5,4,1) = 5.052263D-02
      allz(6,4,1) = 2.809354D-02
      allc(1,4,1) = 3.775056D-03
      allc(2,4,1) = -5.585965D-02
      allc(3,4,1) = -3.192946D-01
      allc(4,4,1) = -2.764780D-02
      allc(5,4,1) = 9.049199D-01
      allc(6,4,1) = 3.406258D-01
!                                    4P
      allc(1,4,2) = -7.052075D-03
      allc(2,4,2) = -5.259505D-02
      allc(3,4,2) = -3.773450D-02
      allc(4,4,2) = 3.874773D-01
      allc(5,4,2) = 5.791672D-01
      allc(6,4,2) = 1.221817D-01
      allz(1,4,2) = 1.365346D+00
      allz(2,4,2) = 4.393213D-01
      allz(3,4,2) = 1.877069D-01
      allz(4,4,2) = 9.360270D-02
      allz(5,4,2) = 5.052263D-02
      allz(6,4,2) = 2.809354D-02
!                                   5S
      allz(1,5,1) = 7.701420258D-01
      allz(2,5,1) = 2.756268915D-01
      allz(3,5,1) = 1.301847480D-01
      allz(4,5,1) = 6.953441940D-02
      allz(5,5,1) = 4.002545502D-02
      allz(6,5,1) = 2.348388309D-02
      allc(1,5,1) = 1.267447151D-02
      allc(2,5,1) = 3.266734789D-03
      allc(3,5,1) = -4.307553999D-01
      allc(4,5,1) = -3.231998963D-01
      allc(5,5,1) = 1.104322879D+00
      allc(6,5,1) = 4.368498703D-01
!                                   5P
      allz(1,5,2) = 7.701420258D-01
      allz(2,5,2) = 2.756268915D-01
      allz(3,5,2) = 1.301847480D-01
      allz(4,5,2) = 6.953441940D-02
      allz(5,5,2) = 4.002545502D-02
      allz(6,5,2) = 2.348388309D-02
      allc(1,5,2) = -1.105673292D-03
      allc(2,5,2) = -6.243132446D-02
      allc(3,5,2) = -1.628476766D-01
      allc(4,5,2) = 3.210328714D-01
      allc(5,5,2) = 6.964579592D-01
      allc(6,5,2) = 1.493146125D-01
!                                   6S
      allz(1,6,1) = 5.800292686D-01
      allz(2,6,1) = 2.718262251D-01
      allz(3,6,1) = 7.938523262D-02
      allz(4,6,1) = 4.975088254D-02
      allz(5,6,1) = 2.983643556D-02
      allz(6,6,1) = 1.886067216D-02
      allc(1,6,1) = 4.554359511D-03
      allc(2,6,1) = 5.286443143D-02
      allc(3,6,1) = -7.561016358D-01
      allc(4,6,1) = -2.269803820D-01
      allc(5,6,1) = 1.332494651D+00
      allc(6,6,1) = 3.622518293D-01
!                                   6P
      allz(1,6,2) = 6.696537714D-01
      allz(2,6,2) = 1.395089793D-01
      allz(3,6,2) = 8.163894960D-02
      allz(4,6,2) = 4.586329272D-02
      allz(5,6,2) = 2.961305556D-02
      allz(6,6,2) = 1.882221321D-02
      allc(1,6,2) = 2.782723680D-03
      allc(2,6,2) = -1.282887780D-01
      allc(3,6,2) = -2.266255943D-01
      allc(4,6,2) = 4.682259383D-01
      allc(5,6,2) = 6.752048848D-01
      allc(6,6,2) = 1.091534212D-01
      ia = 0
      do i = 1, numat
        if (nat(i) == 0) cycle
        ni = nat(i)
        xi = zs(ni)
        ia = ia + 1
        if (ni < 2) then
          nqn = 1
        else if (ni < 10) then
          nqn = 2
        else if (ni < 18) then
          nqn = 3
        else if (ni < 36) then
          nqn = 4
        else if (ni < 54) then
          nqn = 5
        else
          nqn = 6
        end if
        ccc(ia,:) = allc(:,nqn,1)
        zzz(ia,:) = allz(:,nqn,1)*xi**2
        if (nlast(i) > nfirst(i)) then
          ia = ia + 1
          xi = zp(ni)
          ccc(ia,:) = allc(:,nqn,2)
          zzz(ia,:) = allz(:,nqn,2)*xi**2
        end if
      end do
      return
      end subroutine setupg
