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

subroutine delsta (nat, iorbs, p, cdi, dstat, ii, jj)
   !***********************************************************************
   !                                                                      *
   !   DELSTA calculates the derivative of the energy WRT Cartesian       *
   !   coordinates.  Called by DCART only.                                *
   ! On input: CDI:    The coordinates of the two atoms                   *
   !           P:      The density matrix (only the diagonal will be used *
   !           II,JJ:  The atom numbers                                   *
   ! On exit:  DSTAT:  Derivatives (3 numbers)                            *
   !                                                                      *
   !***********************************************************************
    use molkst_C, only: numat, mpack, cutofp
    use funcon_C, only: fpc_9, ev
    use parameters_C, only: tore
   !
   !.. Implicit Declarations ..
    implicit none
   !
   !.. Formal Arguments ..
    integer, intent (in) :: ii, jj
    integer, dimension (numat), intent (in) :: iorbs, nat
    double precision, dimension (3), intent (out) :: dstat
    double precision, dimension (mpack), intent (in) :: p
    double precision, dimension (3, 2), intent (in) :: cdi
   !
   !.. Local Scalars ..
    integer :: k, l
    double precision :: qii, qjj, rij, sum
   !
   !.. Local Arrays ..
    double precision, dimension (3) :: vect
    integer, external :: ijbo
   !
   !.. Intrinsic Functions ..
    intrinsic Sqrt
   !
   ! ... Executable Statements ...
   !
    qii = tore(nat(ii))
    l = ijbo (ii, ii)
    do k = 1, iorbs(ii)
      l = l + k
      qii = qii - p(l)
    end do
    qjj = tore(nat(jj))
    l = ijbo (jj, jj)
    do k = 1, iorbs(jj)
      l = l + k
      qjj = qjj - p(l)
    end do
    rij = Sqrt ((cdi(1, 1)-cdi(1, 2))**2+ (cdi(2, 1)-cdi(2, 2))**2+ (cdi(3, &
   & 1)-cdi(3, 2))**2)
    if (rij > cutofp) then
      do k = 1, 3
        dstat(k) = 0.d0
      end do
    else
      do k = 1, 3
        vect(k) = (cdi(k, 1)-cdi(k, 2)) / rij
      end do
      !
      !   This component of the derivative has an electrostatic term only
      !
      sum = fpc_9 * ev / rij ** 2
      do k = 1, 3
        dstat(k) = -0.5d0 * qjj * qii * sum * vect(k)
      end do
    end if
end subroutine delsta
