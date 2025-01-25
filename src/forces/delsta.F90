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
