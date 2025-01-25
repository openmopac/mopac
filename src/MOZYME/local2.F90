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
