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

      subroutine tred2e(nm, n, a, d, e, z)

      implicit none
      integer :: i, j, k, l, n, ii, nm, jp1
      double precision :: a(nm, n), d(n), e(n), z(nm, n)
      double precision :: f, g, h, hh, scale

!
!     this subroutine is a translation of the algol procedure tred2,
!     num. math. 11, 181 - 195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii - linear algebra, 212 - 226(1971).
!
!     this subroutine reduces a real symmetric matrix to a
!     symmetric tridiagonal matrix using and accumulating
!     orthogonal similarity transformations.
!
!     on input -
!
!        nm must be set to the row dimension of two - dimensional
!          array parameters as declared in the calling program
!          dimension statement,
!
!        n is the order of the matrix,
!
!        a contains the real symmetric input matrix.  only the
!          lower triangle of the matrix need be supplied.
!
!     on output -
!
!        d contains the diagonal elements of the tridiagonal matrix,
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n - 1 positions.  e(1) is set to zero,
!
!        z contains the orthogonal transformation matrix
!          produced in the reduction,
!
!        a and z may coincide.  if distinct, a is unaltered.
!
!     questions and comments should be directed to b. s. garbow,
!     applied mathematics division, argonne national laboratory
!
!    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      do i = 1, n
         do j = 1, i
            z(i, j) = a(i, j)
         end do
      end do
!
      if (n  ==  1) go to 320
!     ********** for i = n step - 1 until 2 do - - **********
      do ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0
         scale = 0.0
         if (l  <  2) go to 130
!     ********** scale row (algol tol then not needed) **********
         do k = 1, l
            scale = scale + abs(z(i, k))
         end do
!
         if (scale /= 0.0) go to 140
  130    e(i) = z(i, l)
         go to 290
!
  140    do k = 1, l
            z(i, k) = z(i, k) / scale
            h = h + z(i, k) * z(i, k)
         end do
!
         f = z(i, l)
         g = -sign(sqrt(h), f)
         e(i) = scale * g
         h = h - f * g
         z(i, l) = f - g
         f = 0.0
!
         do j = 1, l
            z(j, i) = z(i, j) / (scale * h)
            g = 0.0
!     ********** form element of a*u **********
            do k = 1, j
               g = g + z(j, k) * z(i, k)
            end do
!
            jp1 = j + 1
            if (l  <  jp1) go to 220
!
            do k = jp1, l
               g = g + z(k, j) * z(i, k)
            end do
!     ********** form element of p **********
  220       e(j) = g / h
            f = f + e(j) * z(i, j)
         end do
!
         hh = f / (h + h)
!     ********** form reduced a **********
         do j = 1, l
            f = z(i, j)
            g = e(j) - hh * f
            e(j) = g
            do k = 1, j
               z(j, k) = z(j, k) - f * e(k) - g * z(i, k)
            end do
         end do
!
         do k = 1, l
            z(i, k) = scale * z(i, k)
         end do
!
  290    d(i) = h
      end do
!
  320 d(1) = 0.0
      e(1) = 0.0
!     ********** accumulation of transformation matrices **********
      do i = 1, n
         l = i - 1
         if (d(i)  ==  0.0) go to 380
!
         do j = 1, l
            g = 0.0
!
            do k = 1, l
               g = g + z(i, k) * z(k, j)
            end do
!
            do k = 1, l
               z(k, j) = z(k, j) - g * z(k, i)
            end do
         end do
!
  380    d(i) = z(i, i)
         z(i, i) = 1.d0
         if (l  >=  1) then
            do j = 1, l
               z(i, j) = 0.0
               z(j, i) = 0.0
            end do
         endif
!
      end do
!
      return
!     ********** last card of tred2 **********
      end subroutine tred2e
!
!    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      subroutine tql2e(nm, n, d, e, z, ierr)

      implicit none
      integer :: i, j, k, l, m, n, ii, nm, mml, ierr
      double precision :: d(n), e(n), z(nm, n)
      double precision :: b, c, f, g, h, p, r, s, machep
!
!     this subroutine is a translation of the algol procedure tql2,
!     num. math. 11, 293 - 306(1968) by bowdler, martin, reinsch, and
!     wilkinson.
!     handbook for auto. comp., vol.ii - linear algebra, 227 - 240(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     on input -
!
!        nm must be set to the row dimension of two - dimensional
!          array parameters as declared in the calling program
!          dimension statement,
!
!        n is the order of the matrix,
!
!        d contains the diagonal elements of the input matrix,
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n - 1 positions.  e(1) is arbitrary,
!
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output -
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1, 2, ..., ierr - 1,
!
!        e has been destroyed,
!
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues,
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j - th eigenvalue has not been
!                     determined after 30 iterations.
!
!     questions and comments should be directed to b. s. garbow,
!     applied mathematics division, argonne national laboratory
!
!    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!     ********** machep is a machine dependent parameter specifying
!                the relative precision of floating point arithmetic.
!
!                **********
      machep = 2.d0**(- 52)

      ierr = 0
      if (n  ==  1) go to 1001

      do i = 2, n
         e(i - 1) = e(i)
      end do

      f = 0.0
      b = 0.0
      e(n) = 0.0

      do 240 l = 1, n
         j = 0
         h = machep * (abs(d(l)) + abs(e(l)))
         if (b  <  h) b = h
!     ********** look for small sub - diagonal element **********
         do 110 m = l, n
            if (abs(e(m))  <=  b) go to 120
!     ********** e(n) is always zero, so there is no exit
!                through the bottom of the loop **********
  110    continue

  120    if (m  ==  l) go to 220
  130    if (j  ==  30) go to 1000
         j = j + 1
!     ********** form shift **********
         p = (d(l + 1) - d(l)) / (2.0 * e(l))
         r = sqrt(p*p + 1.d0)
         h = d(l) - e(l) / (p + sign(r, p))

         do i = l, n
            d(i) = d(i) - h
         end do

         f = f + h
!     ********** ql transformation **********
         p = d(m)
         c = 1.d0
         s = 0.0
         mml = m - l
!     ********** for i = m - 1 step - 1 until l do - - **********
         do 200 ii = 1, mml
            i = m - ii
            g = c * e(i)
            h = c * p
            if (abs(p)  <  abs(e(i))) go to 150
            c = e(i) / p
            r = sqrt(c*c + 1.d0)
            e(i + 1) = s * p * r
            s = c / r
            c = 1.0 / r
            go to 160
  150       c = p / e(i)
            r = sqrt(c*c + 1.d0)
            e(i + 1) = s * e(i) * r
            s = 1.0 / r
            c = c * s
  160       p = c * d(i) - s * g
            d(i + 1) = h + s * (c * g + s * d(i))
!     ********** form vector **********
            do 180 k = 1, n
               h = z(k, i + 1)
               z(k, i + 1) = s * z(k, i) + c * h
               z(k, i) = c * z(k, i) - s * h
  180       continue

  200    continue

         e(l) = s * p
         d(l) = c * p
         if (abs(e(l))  >  b) go to 130
  220    d(l) = d(l) + f
  240 continue
!     ********** order eigenvalues and eigenvectors **********
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)

         do 260 j = ii, n
            if (d(j)  >=  p) go to 260
            k = j
            p = d(j)
  260    continue

         if (k  ==  i) go to 300
         d(k) = d(i)
         d(i) = p

         do 280 j = 1, n
            p = z(j, i)
            z(j, i) = z(j, k)
            z(j, k) = p
  280    continue

  300 continue

      go to 1001
!     ********** set error - - no convergence to an
!                eigenvalue after 30 iterations **********
 1000 ierr = l
 1001 return
!     ********** last card of tql2 **********
      end subroutine tql2e
