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

  subroutine rotate(ni, nj, xi, xj, w, kr, e1b, e2a, enuc)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
  use molkst_C, only : id, method_pm7
  use parameters_C, only: natorb
!***********************************************************************
!
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
  double precision , intent(out) :: w(2025)
  double precision , intent(out) :: e1b(45)
  double precision , intent(out) :: e2a(45)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: i, li, lj, ik, k
  double precision ::  rijx, rij
  double precision, dimension(3) :: x
  double precision, dimension(171) :: en
!-----------------------------------------------
!***********************************************************************
!
!..IMPROVED SCALAR VERSION
!..WRITTEN BY ERNEST R. DAVIDSON, INDIANA UNIVERSITY.
!
!
!   ROTATE CALCULATES THE TWO-PARTICLE INTERACTIONS.
!
!   ON INPUT  NI     = ATOMIC NUMBER OF FIRST ATOM.
!             NJ     = ATOMIC NUMBER OF SECOND ATOM.
!             XI     = COORDINATE OF FIRST ATOM.
!             XJ     = COORDINATE OF SECOND ATOM.
!
! ON OUTPUT W      = ARRAY OF TWO-ELECTRON REPULSION INTEGRALS.
!           E1B,E2A= ARRAY OF ELECTRON-NUCLEAR ATTRACTION INTEGRALS,
!                    E1B = ELECTRON ON ATOM NI ATTRACTING NUCLEUS OF NJ.
!           ENUC   = NUCLEAR-NUCLEAR REPULSION TERM.
!
!
! *** THIS ROUTINE COMPUTES THE REPULSION AND NUCLEAR ATTRACTION
!     INTEGRALS OVER MOLECULAR-FRAME COORDINATES.  THE INTEGRALS OVER
!     LOCAL FRAME COORDINATES ARE EVALUATED BY SUBROUTINE REPP AND
!     STORED AS FOLLOWS (WHERE P-SIGMA = O,   AND P-PI = P AND P* )
!     IN RI
!     (SS/SS)=1,   (SO/SS)=2,   (OO/SS)=3,   (PP/SS)=4,   (SS/OS)=5,
!     (SO/SO)=6,   (SP/SP)=7,   (OO/SO)=8,   (PP/SO)=9,   (PO/SP)=10,
!     (SS/OO)=11,  (SS/PP)=12,  (SO/OO)=13,  (SO/PP)=14,  (SP/OP)=15,
!     (OO/OO)=16,  (PP/OO)=17,  (OO/PP)=18,  (PP/PP)=19,  (PO/PO)=20,
!     (PP/P*P*)=21,   (P*P/P*P)=22.
!
!***********************************************************************
  x(1) = xi(1) - xj(1)
  x(2) = xi(2) - xj(2)
  x(3) = xi(3) - xj(3)
  rij = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
  if (rij < 0.00002D0) then
!
!     SMALL RIJ CASE
!
    e1b = 0.D0
    e2a = 0.D0
    w = 0.D0
    enuc = 0.D0
    return
  end if
!
!
! *** THE REPULSION INTEGRALS OVER MOLECULAR FRAME (W) ARE STORED IN THE
!     ORDER IN WHICH THEY WILL LATER BE USED.  IE.  (I,J/K,L) WHERE
!     J.LE.I  AND  L.LE.K     AND L VARIES MOST RAPIDLY AND I LEAST
!     RAPIDLY.  (ANTI-NORMAL COMPUTER STORAGE)
!
!
  rijx = sqrt(rij)
  rij = rijx
!
! *** COMPUTE INTEGRALS IN DIATOMIC FRAME
!
  li = natorb(ni)
  lj = natorb(nj)
  call rotatd(ni, nj, xi, xj, w, kr, enuc)
  en(:) = 0.d0
  call elenuc (1, li, li + 1, li + lj, en)
  ik = 0
  do i = 1, li
    do k = 1, i
      ik = ik + 1
      e1b(ik) = en((i*(i - 1))/2 + k)
    end do
  end do
  ik = 0
  do i = li + 1, li + lj
    do k = li + 1, i
      ik = ik + 1
      e2a(ik) = en((i*(i - 1))/2 + k)
    end do
  end do
  if (id == 3 .and. .not. method_PM7) call nddo_to_point(w, e1b, e2a, enuc, rij, ni, nj) ! In solids, transition to point-charge as Rij increases
  return
  end subroutine rotate
