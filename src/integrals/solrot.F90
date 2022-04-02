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

      subroutine solrot(ni, nj, xi, xj, wj, wk, kr, e1b, e2a, enuc)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use common_arrays_C, only : tvec
      use molkst_C, only : numcal, l1u, l2u, l3u, clower
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ni
      integer  :: nj
      integer , intent(inout) :: kr
      double precision , intent(out) :: enuc
      double precision , intent(in) :: xi(3)
      double precision , intent(in) :: xj(3)
      double precision , intent(out) :: wj(2025)
      double precision , intent(out) :: wk(2025)
      double precision , intent(inout) :: e1b(45)
      double precision , intent(inout) :: e2a(45)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icalcn, i, j, k, kb
      double precision, dimension(:), allocatable :: wsum, wmax, wbits
      double precision, dimension(3) :: xjuc
      double precision, dimension(45) :: e1bits = 0.d0, e2bits = 0.d0
      double precision, dimension(3) :: xdumy
      double precision :: one, cutof2, r, enubit

      save xdumy, icalcn, cutof2
!-----------------------------------------------
!***********************************************************************
!
!   SOLROT FORMS THE TWO-ELECTRON TWO-ATOM J AND K INTEGRAL STRINGS.
!          ON EXIT WJ = "J"-TYPE INTEGRALS
!                  WK = "K"-TYPE INTEGRALS
!
!      FOR MOLECULES, WJ = WK.
!***********************************************************************
      data icalcn/ 0/
      data xdumy/ 3*0.D0/
      if (icalcn /= numcal) then
        icalcn = numcal
        cutof2 = clower**2
      end if
      one = 1.D0
      if (Abs(xi(1) - xj(1)) <1.d-20 .and. Abs(xi(2) - xj(2)) <1.d-20 .and. &
      Abs(xi(3) - xj(3)) <1.d-20) one = 0.5D0
      allocate (wsum(2025), wmax(2025), wbits(2025))
      wmax = 0.D0
      wsum = 0.D0
      wbits = 0.D0
      e1b = 0.D0
      e2a = 0.D0
      enuc = 0.D0
      do i = -l1u, l1u
        do j = -l2u, l2u
          do k = -l3u, l3u
            xjuc = xj + tvec(:,1)*i + tvec(:,2)*j + tvec(:,3)*k - xi
            r = xjuc(1)**2 + xjuc(2)**2 + xjuc(3)**2
            if (r > cutof2) then
!
!  Interaction distance is greater than cutoff, so use point-charge
!
              r = sqrt(r)
              call point(r, ni, nj, wbits, kb, e1bits, e2bits, enubit)
            else
!
!  Interaction distance is less than cutoff
!
              kb = 0
              call rotate (ni, nj, xdumy, xjuc, wbits, kb, e1bits, e2bits, enubit)
            end if
            wsum(:kb) = wsum(:kb) + wbits(:kb)
            if (wmax(1) < wbits(1)) wmax(:kb) = wbits(:kb) ! "K" integrals apply only to nearest pair
            e1b = e1b + e1bits
            e2a = e2a + e2bits
            enuc = enuc + enubit*one
          end do
        end do
      end do
      if (one < 0.9D0) wmax(:kb) = 0.D0
      wk(:kb) = wmax(:kb)
      wj(:kb) = wsum(:kb)
      kr = kb + kr
      deallocate(wmax, wbits, wsum)
      return
      end subroutine solrot
      subroutine nddo_to_point(wbits, e1bits, e2bits, enubit, r, ni, nj)
!
! NDDO_to_point smoothly transitions the two center integrals from NDDO to point-charge
! as the distance increases.  This is done by steadily mixing more and more point-charge
! contribution as the distance increases.
!
!  Up to 3 Angstroms, use "exact" NDDO, beyond 3 Angstroms, use a Gauusian mixture.
!
      use molkst_C, only : l_feather
        implicit none
        double precision, dimension(45) :: e1bits, e2bits, e1bits_p = 0.d0, e2bits_p = 0.d0
        double precision, dimension(2025) :: wbits
        double precision, intent (inout):: r
        double precision, intent (inout):: enubit
        integer, intent (in) :: ni, nj
        integer :: kb
        double precision :: const, enubit_p, dummy
        double precision, dimension(2025) :: wbits_p

        if (l_feather) then
          call to_point(r, dummy, const)
        else
          if (r < 3.0d0) return
          const = Exp (-0.025d0*(r - 3.d0)**2)
        end if
        call point (r, ni, nj, wbits_p, kb, e1bits_p, e2bits_p, enubit_p)
        wbits(1:kb) = const*wbits(1:kb) +(1.d0-const)*wbits_p(1:kb)
        e1bits = const*e1bits + (1.d0 - const)*e1bits_p
        e2bits = const*e2bits + (1.d0 - const)*e2bits_p
        enubit = const*enubit + (1.d0 - const)*enubit_p
      end subroutine nddo_to_point

  subroutine point (r, ni, nj, w, kr, e1b, e2a, enuc)
  use funcon_C, only : a0, ev
  use parameters_C, only : tore, natorb
  implicit none
  integer, intent (in) :: ni, nj
  integer, intent (out) :: kr
  double precision, intent (inout) :: r
  double precision, intent (out) :: enuc
  double precision, dimension (45), intent (out) :: e1b, e2a
  double precision, dimension (*), intent (out) :: w
  integer :: i, j, ii, jj, nii, njj
  double precision :: ee, ee1, ee2
  double precision, external :: trunk
!
!   In solid-state work, if an interatomic distance is larger than
!   clower, then the NDDO approximation is replaced by a point-charge
!   approximation, in which the point-charge is located at a distance
!   that depends on the interatomic distance (see trunk)
!
    r = trunk (r)
    ee = ev * a0 / r
    ee1 = -ee * tore(nj)
    ee2 = -ee * tore(ni)
    i = natorb(ni)
    j = natorb(nj)
    nii = (i*(i+1)) / 2
    njj = (j*(j+1)) / 2
    kr = nii * njj
    !
    !  Fill "w" array
    !
    w(1:kr) = 0.d0
    i_loop: do ii = 1, i
      j_loop: do jj = 1, j
        w(((ii*(ii+1))/2-1)*njj+ (jj*(jj+1))/2) = ee
      end do j_loop
    end do i_loop
    !
    ! Fill e1b and e2a
    !
    e1b(1:nii) = 0
    e2a(1:njj) = 0
    do ii = 1, i
      e1b((ii*(ii+1))/2) = ee1
    end do
    do jj = 1, j
      e2a((jj*(jj+1))/2) = ee2
    end do
    !
    !  Compute nuclear contribution
    !
    enuc = -ee1 * tore(ni)
end subroutine point
double precision function trunk (r)
!
!  In solids, change the apparent distance so that the Madelung sum can be
!  solved.
!
!   r: Distance in Angstroms
!
    use molkst_C, only: numcal, clower, cutofp, cupper
    implicit none
    double precision, intent (in) :: r
    integer, save :: icalcn = 0
    double precision, save :: bound1, bound2, c, clim, cr, cr2, range
    if (icalcn /= numcal) then
      bound1 = clower / cutofp
      !
      !   CUPPER = upper bound of truncation function, as a function of CUTOFP
      !
      cupper = cutofp + 0
      bound2 = cupper / cutofp
      range = bound2 - bound1
      !
      !   Truncation function = C+CR*R+CR2*R**2
      !
      c = -0.5d0 * bound1 ** 2 * cutofp / range
      cr = 1.0d0 + bound1 / range
      cr2 = -1 / (cutofp*2*range)
      !
      !   Above CUPPER function = constant
      !
      clim = c + cr * cupper + cr2 * cupper ** 2
      icalcn = numcal
    end if
   !
   !  Need a smooth function in the region of CUTOFP
   !
   !  Function has form:
   !    Up to CLOWER  R=R
   !    At CLOWER,  slope = 1.0
   !    Between CLOWER and CUPPER = Monotomic increase,
   !    with rate of increase dropping from 1.0 to 0.0
   !    At CUPPER, slope = 0
   !    Above CUTOFP   R=CUTOFP
    if (r > clower) then
      if (r > cupper) then
        trunk = clim
      else
        trunk = c + cr * r + cr2 * r ** 2
      end if
    else
      trunk = r
    end if
end function trunk
