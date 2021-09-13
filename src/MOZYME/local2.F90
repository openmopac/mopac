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

subroutine local2 (c, mdim, nmos, nfirst, nlast, numat)
   !
   !.. Implicit Declarations ..
    implicit none
   !
   !.. Formal Arguments ..
    integer, intent (in) :: mdim, nmos, numat
    integer, dimension (*), intent (in) :: nfirst, nlast
    double precision, dimension (mdim, mdim), intent (inout) :: c
   !
   !.. Local Scalars ..
    integer :: i, j, k, k1, loop, n, niter
    double precision :: aij, bij, ca, dii, dij, djj, sa, sum, xiiii, xiijj, &
   & xijij, xijjj, xjiii, xjjjj
   !
   !.. Local Arrays ..
    double precision, dimension (20) :: psi1, psi2
   !
   !.. Intrinsic Functions ..
    intrinsic Abs, Sqrt
   !
   ! ... Executable Statements ...
   !
    n = nlast(numat)
    niter = 1
    if (nmos == 4) then
      niter = 19
    end if
    do loop = 1, niter
      sum = 0.d0
      do i = 1, nmos
        do j = 1, nmos
          if (j /= i) then
            xijjj = 0.0d0
            xjiii = 0.0d0
            xiiii = 0.0d0
            xjjjj = 0.0d0
            xijij = 0.0d0
            xiijj = 0.0d0
            do k = 1, n
              psi1(k) = c(k, i)
              psi2(k) = c(k, j)
            end do
            do k1 = 1, numat
              dij = 0.d0
              dii = 0.d0
              djj = 0.d0
              do k = nfirst(k1), nlast(k1)
                dij = dij + psi1(k) * psi2(k)
                dii = dii + psi1(k) * psi1(k)
                djj = djj + psi2(k) * psi2(k)
              end do
              xijjj = xijjj + dij * djj
              xjiii = xjiii + dij * dii
              xiiii = xiiii + dii * dii
              xjjjj = xjjjj + djj * djj
              xijij = xijij + dij * dij
              xiijj = xiijj + dii * djj
            end do
            aij = xijij - (xiiii+xjjjj-2.0d0*xiijj) / 4.0d0
            bij = xjiii - xijjj
            ca = Sqrt (aij*aij+bij*bij)
            sa = aij + ca
            if (sa > 1.0d-14) then
              ca = (1.0d0+Sqrt((1.0d0-aij/ca)/2.0d0)) / 2.0d0
              sa = Sqrt (1.0d0-ca)
              sum = sum + Abs (sa)
              ca = Sqrt (ca)
              do k = 1, n
                c(k, i) = ca * psi1(k) + sa * psi2(k)
                c(k, j) = -sa * psi1(k) + ca * psi2(k)
              end do
            end if
          end if
        end do
      end do
      if (sum < 1.d-5) exit
    end do
end subroutine local2
