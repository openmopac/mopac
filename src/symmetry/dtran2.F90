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

      subroutine dtran2(r, t, ioper)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use funcon_C, only : pi
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ioper
      double precision , intent(inout) :: r(3,3)
      double precision , intent(inout) :: t(5,5,12)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      double precision, dimension(2,4) :: f
      double precision :: tol, s12, s3, one, r1, r2, r3, check, arg, b, sina, g&
        , a, e1, x1, e2, e3, e4, x2, x3, x4, ta, tg
      logical :: right

      save tol, s12, s3, one
!-----------------------------------------------
!
!  3.464... =SQRT(12)
!
      data tol, s12/ 0.001D0, 3.46410161513D0/
      data s3, one/ 1.73205080756D0, 1.D0/
      r1 = r(2,1)*r(3,2) - r(3,1)*r(2,2)
      r2 = r(3,1)*r(1,2) - r(1,1)*r(3,2)
      r3 = r(1,1)*r(2,2) - r(2,1)*r(1,2)
      check = r1*r(1,3) + r2*r(2,3) + r3*r(3,3)
      right = check > 0.D0
      r(1,3) = r1
      r(2,3) = r2
      r(3,3) = r3
      arg = r3
      if (abs(arg) > one) arg = sign(one,arg)
      b = acos(arg)
      sina = sqrt(1.D0 - arg*arg)
      if (sina >= tol) then
        arg = r(3,2)/sina
        if (abs(arg) > one) arg = sign(one,arg)
        g = asin(arg)
        arg = r(2,3)/sina
        if (abs(arg) > one) arg = sign(one,arg)
        a = asin(arg)
      else
        arg = r(1,2)
        if (abs(arg) > one) arg = sign(one,arg)
        g = asin(arg)
        a = 0.D0
      end if
      f(1,1) = a
      f(1,2) = a
      f(1,3) = pi - a
      f(1,4) = pi - a
      f(2,1) = g
      f(2,3) = g
      f(2,2) = pi - g
      f(2,4) = pi - g
      do i = 1, 4
        a = f(1,i)
        g = f(2,i)
        check = abs(sin(b)*cos(a) + r(1,3))
        if (check > tol) cycle
        check = (-sin(g)*cos(b)*sin(a)) + cos(g)*cos(a)
        if (abs(check - r(2,2)) > tol) cycle
        check = sin(a)*cos(g) + cos(a)*cos(b)*sin(g)
        if (abs(check - r(1,2)) > tol) cycle
        exit
      end do
      g = -g
      a = -a
      b = -b
      e1 = cos(b*0.5D0)
      x1 = -sin(b*0.5D0)
      e2 = e1*e1
      e3 = e1*e2
      e4 = e2*e2
      x2 = x1*x1
      x3 = x1*x2
      x4 = x2*x2
      ta = 2.D0*a
      tg = 2.D0*g
      t(1,1,ioper) = e4*cos(ta + tg) + x4*cos(ta - tg)
      t(1,2,ioper) = 2.D0*e3*x1*cos(a + tg) - 2.D0*e1*x3*cos(a - tg)
      t(1,3,ioper) = 2.D0*s3*e2*x2*cos(tg)
      t(1,4,ioper) = 2.D0*e3*x1*sin(a + tg) - 2.D0*e1*x3*sin(a - tg)
      t(1,5,ioper) = e4*sin(ta + tg) + x4*sin(ta - tg)
      t(2,1,ioper) = 2.D0*e1*x3*cos(ta - g) - 2.D0*e3*x1*cos(ta + g)
      t(2,2,ioper) = (e4 - 3.D0*e2*x2)*cos(a + g) - (3.D0*e2*x2 - x4)*cos(a - g)
      t(2,3,ioper) = 2.D0*s3*(e3*x1 - e1*x3)*cos(g)
      t(2,4,ioper) = (e4 - 3.D0*e2*x2)*sin(a + g) - (3.D0*e2*x2 - x4)*sin(a - g)
      t(2,5,ioper) = (-2.D0*e3*x1*sin(ta + g)) + 2.D0*e1*x3*sin(ta - g)
      t(3,1,ioper) = s12*e2*x2*cos(ta)
      t(3,2,ioper) = -s12*(e3*x1 - e1*x3)*cos(a)
      t(3,3,ioper) = e4 - 4.D0*e2*x2 + x4
      t(3,4,ioper) = -s12*(e3*x1 - e1*x3)*sin(a)
      t(3,5,ioper) = s12*e2*x2*sin(ta)
      t(4,1,ioper) = 2.D0*e1*x3*sin(ta - g) + 2.D0*e3*x1*sin(ta + g)
      t(4,2,ioper) = (-(e4 - 3.D0*e2*x2)*sin(a + g)) - (3.D0*e2*x2 - x4)*sin(a - g)
      t(4,3,ioper) = -2.D0*s3*(e3*x1 - e1*x3)*sin(g)
      t(4,4,ioper) = (e4 - 3.D0*e2*x2)*cos(a + g) + (3.D0*e2*x2 - x4)*cos(a - g)
      t(4,5,ioper) = (-2.D0*e3*x1*cos(ta + g)) - 2.D0*e1*x3*cos(ta - g)
      t(5,1,ioper) = (-e4*sin(ta + tg)) + x4*sin(ta - tg)
      t(5,2,ioper) = (-2.D0*e3*x1*sin(a + tg)) - 2.D0*e1*x3*sin(a - tg)
      t(5,3,ioper) = -2.D0*s3*e2*x2*sin(tg)
      t(5,4,ioper) = 2.D0*e3*x1*cos(a + tg) + 2.D0*e1*x3*cos(a - tg)
      t(5,5,ioper) = e4*cos(ta + tg) - x4*cos(ta - tg)
      if (right) return
      t(2,:,ioper) = -t(2,:,ioper)
      t(4,:,ioper) = -t(4,:,ioper)
      return
      end subroutine dtran2
