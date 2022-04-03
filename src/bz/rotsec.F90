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

subroutine rotsec (t1, transform, new_coord, loop)
!
!  Using the symmetry operation stored in xop(:,loop), operate on the coordinates
!  in new_coord to generate a new set of coordinates.
!
!  If transform is true:
!
!  Using the symmetry operation stored in xop(:,loop), rotate the unit matrix t1
!  to form a unitary matrix representing that operation.
!
!  For electronics, t1 represents atomic orbitals, 
!  for phonons, t1 represents the Cartesian directions, i.e., the vibrational coordinates
!
  use common_common, only : numat, id, tvec, xop, phonon, title, trans, keywrd, &
    iw_new, nijk, ixyz, nijk_CUC
  implicit none
  double precision, dimension(9, 9), intent(out) :: t1
  logical, intent(in) :: transform
  integer, intent(in) :: loop
  double precision, dimension(3, numat), intent(inout) :: new_coord
! 
  logical :: first = .true.
  logical :: prt
  integer :: i, j, k, l, phon
  double precision :: c1, c2, c3, ca, cb, one, r, sa, sb, six, sum, theta, &
       & tnpx, tnpy, tnpz, two, x, xnew, xx, xy, y, ynew, yx, z, znew, zx, &
       & a, b, c, d, e, f, g, h, o, x1, y1, z1, e1(9,9), e2(9,9), sum1
  double precision, dimension(3, 3) :: atrans, atxyz, bb, ctox, xtoc
  double precision, dimension(3, 3) :: at
  double precision, dimension(5, 5) :: aa, cc
  save
! 
!.. Equivalences .. 
  equivalence (at(1, 1), a)
  equivalence (at(2, 1), b)
  equivalence (at(3, 1), c)
  equivalence (at(1, 2), d)
  equivalence (at(2, 2), e)
  equivalence (at(3, 2), f)
  equivalence (at(1, 3), g)
  equivalence (at(2, 3), h)
  equivalence (at(3, 3), o)
  if (first) then
    prt = (Index (keywrd, " ROTFOK") /= 0)
    do j = 1, id
      if (tvec(j, j) < 0) then
        do k = 1, 3
          tvec(j, k) = -tvec(j, k)
        end do
      end if
    end do
    do j = id + 1, 3
      tvec(j, j) = 1.d0
    end do
    do j = 1, 3
      do k = 1, 3
        sum = 0.d0
        do l = 1, 3
          sum = sum + trans(l, j)*tvec(l, k)
        end do
        xtoc(k, j) = sum
        ctox(k, j) = sum
      end do
    end do
    if (prt) then
      write(iw_new, '(/,a)') " Transform CRYSTAL to XYZ Coordinates"
    end if
    if (prt) then
      write(iw_new, "(3f13.6)") ctox
    end if
    call minv (xtoc, 3, sum)
    if (prt) then
      write(iw_new, '(/,a)') " Transform XYZ to CRYSTAL Coordinates (The inverse of the previous transform)"
    end if
    if (prt) then
      write(iw_new, "(3f13.6)") xtoc
    end if
    if (prt) then
      write(iw_new, '(/,a)') " Unitary Transform: XYZ to CRYSTAL"
    end if
    if (prt) then
      write(iw_new, "(3f13.6)") trans
    end if
    first = .false.
  end if
  if (transform) then
    t1 = 0.d0
!
!     Specification of the symmetry operations
!
!  xop(1,  loop) = 1 if inversion, 0 otherwise
!  xop(2,  loop) = motion in fractional unit cell coordinates in direction "a"
!  xop(3,  loop) = motion in fractional unit cell coordinates in direction "b"
!  xop(4,  loop) = motion in fractional unit cell coordinates in direction "c"
!  xop(5,  loop) = Rotation angle 
!  xop(6,  loop) = "a" component of axis
!  xop(7,  loop) = "b" component of axis
!  xop(8,  loop) = "c" component of axis
!  xop(9,  loop) = "x" component of Cartesian location of the operation
!  xop(10, loop) = "y" component of Cartesian location of the operation
!  xop(11, loop) = "z" component of Cartesian location of the operation

    i = Nint (xop(1, loop))
    tnpx = xop(2, loop)
    tnpy = xop(3, loop)
    tnpz = xop(4, loop)
    theta = xop(5, loop)
    c1 = xop(7, loop)
    c2 = xop(6, loop)
    c3 = xop(8, loop)
    one = 1.d0
    if (i == 1) one = -1.d0
    theta = 2.d0 * 3.141592653589d0 * theta
    bb(1, 3) = 0.d0
    bb(2, 3) = 0.d0
    bb(3, 1) = 0.d0
    bb(3, 2) = 0.d0
    bb(1, 1) = Cos(theta) * one
    bb(2, 2) = bb(1, 1)
    bb(1, 2) = Sin(theta) * one
    bb(2, 1) = -bb(1, 2)
    bb(3, 3) = one
    if (i == 2) bb(3, 3) = -1.d0
!
!  bb now holds the symmetry operation, at this point the operation
!  is about the "Z" axis.
!
!  Now rotate the operation to be about the axis defined by xop(6:8,loop)
!  in crystallographic coordinates
!
    xy = c1*c1 + c2*c2
    if (xy == 0.d0) then
      ca = 1.d0
      cb = 1.d0
      sa = 0.d0
      sb = 0.d0
    else
      r = Sqrt (xy + c3*c3)
      xy = Sqrt (xy)
      ca = c1 / xy
      cb = c3 / r
      sa = c2 / xy
      sb = xy / r
    end if
    a = ca
    b = sa * cb
    c = sa * sb
    d = -sa
    e = ca * cb
    f = ca * sb
    g = 0.d0
    h = -sb
    o = cb
    do i = 1, 3
      do j = 1, 3
        x = 0.d0
        do k = 1, 3
          x = x + at(k, i)*bb(j, k)
        end do
        cc(i, j) = x
      end do
    end do
    do i = 1, 3
      do j = 1, 3
        x = 0.d0
        do k = 1, 3
          x = x + cc(j, k)*at(k, i)
        end do
        aa(i, j) = x
      end do
    end do
    do i = 1, 3
      do j = 1, 3
        at(i, j) = aa(i, j)
      end do
    end do
  else
    if(phonon) then
      phon = 0
    else
      phon = 1
    end if
    do i = 1, 3
      l = i + phon
      do j = 1, 3
        k = j + phon
        at(i, j) = t1(k, l)
      end do
    end do
  end if
!
!  Now to convert from crystal coordinates to Cartesian coordinates,
!  in order to operate on the position of the atoms.
!
!   The operation is TRANS*AT*TRANS
!
  do i = 1, 3
    do j = 1, 3
      sum = 0.d0
      do k = 1, 3
        sum = sum + at(k, i)*trans(k, j)
      end do
      atrans(j, i) = sum
    end do
  end do
  do i = 1, 3
    do j = 1, 3
      sum = 0.d0
      do k = 1, 3
        sum = sum + atrans(i, k)*trans(k, j)
      end do
      atxyz(i, j) = at(i, j)
      at(i, j) = sum
    end do
  end do
  if (prt .and. transform) then
    write(iw_new, '(/,a,i2,a)') " Unitary Transform representing Operation ", loop, ", "//trim(title(loop))
    write(iw_new, '(a)') "            Cartesian Coordinates"    
    if (abs(xop(2,loop)) + abs(xop(3,loop)) + abs(xop(4,loop)) > 0.1d0) then
      write(iw_new, '(a)')"     From:  x         y         z            Fractional crystal translation"
      write(iw_new, "(a,3f10.6, 10x, f13.4)") "     x ", at(1,:), xop(2,loop)
      write(iw_new, "(a,3f10.6, 10x, f13.4)") "To:  y ", at(2,:), xop(3,loop)
      write(iw_new, "(a,3f10.6, 10x, f13.4)") "     z ", at(3,:), xop(4,loop)
    else
      write(iw_new, '(a)')"     From:  x         y         z"
      write(iw_new, "(a,3f10.6)") "     x ", at(1,:)
      write(iw_new, "(a,3f10.6)") "To:  y ", at(2,:)
      write(iw_new, "(a,3f10.6)") "     z ", at(3,:)
    end if
  end if
  do i = 1, numat
!
!  GET CRYSTAL COORDINATES
!
    xx = new_coord(1, i)
    yx = new_coord(2, i)
    zx = new_coord(3, i)
!
!  CONVERT FROM CRYSTALLOGRAPHIC TO CARTESIAN COORDINATES
!
    x1 = xx*ctox(1, 1) + yx*ctox(2, 1) + zx*ctox(3, 1) - xop(9, loop)
    y1 = xx*ctox(1, 2) + yx*ctox(2, 2) + zx*ctox(3, 2) - xop(10, loop)
    z1 = xx*ctox(1, 3) + yx*ctox(2, 3) + zx*ctox(3, 3) - xop(11, loop)
!
!  PERFORM THE OPERATION IN CARTESIAN COORDINATES
!
    xnew = x1*a + y1*d + z1*g + xop(9, loop)
    ynew = x1*b + y1*e + z1*h + xop(10, loop)
    znew = x1*c + y1*f + z1*o + xop(11, loop)
!
!  CONVERT BACK TO CRYSTAL
!
    x = xnew*xtoc(1, 1) + ynew*xtoc(2, 1) + znew*xtoc(3, 1) + tnpx
    y = xnew*xtoc(1, 2) + ynew*xtoc(2, 2) + znew*xtoc(3, 2) + tnpy
    z = xnew*xtoc(1, 3) + ynew*xtoc(2, 3) + znew*xtoc(3, 3) + tnpz
!
!  Store the new coordinates in new_coord
!
    if (prt .and. .not. transform) then
      write(iw_new, "(i4, i8,2i3, i4,a,i4,4(f14.4, 2f10.4))") ixyz, nijk(1,ixyz) + nijk_CUC(1,i), &
      nijk(2,ixyz) + nijk_CUC(2,i), nijk(3,ixyz) + nijk_CUC(3,i),  i, "  =", &
      (ixyz - 1)*numat + i, xx, yx, zx, x1, y1, z1, xnew, ynew, znew, x,y,z
    end if
    new_coord(1, i) = x
    new_coord(2, i) = y
    new_coord(3, i) = z
  end do
  if (.not. transform) return
  if( phonon ) then
!
!  A phonon spectrum - use "p" orbital rotation matrix only
!
    do i = 1,3
      do j = 1,3
         t1(i,j) = atxyz(j,i)
      end do
    end do
    return
  end if
!
!  "d"-orbital transform
!
  two = Sqrt (0.5d0)
  six = Sqrt (1/6.d0)
  at = atxyz
!
!   d(x^2 - y^2)
!
    aa(1, 1) = (a*a - b*b - d*d + e*e)/2.d0
    aa(1, 2) = a*c - d*f  
    aa(1, 3) = (2.d0*c*c - b*b - a*a - 2.d0*f*f + e*e + d*d) * two * six  
    aa(1, 4) = b*c - e*f  
    aa(1, 5) = a*b - d*e 
!
!  d(xz)
! 
    aa(2, 1) = a*g - b*h  
    aa(2, 2) = a*o + c*g  
    aa(2, 3) = (2.d0*c*o - b*h - a*g) * six/two  
    aa(2, 4) = b*o + c*h  
    aa(2, 5) = a*h + b*g  
!
!  d(2z^2 - x^2 - y^2)
!
    aa(3, 1) = (2.d0*g*g - d*d - a*a - 2.d0*h*h + e*e + b*b) * two * six  
    aa(3, 2) = (2.d0*g*o - d*f - a*c) * six/two
    aa(3, 3) = (4.d0*o*o - 2.d0*(h*h + g*g + f*f + c*c) + e*e + d*d + b*b + a*a)/6.d0 
    aa(3, 4) = (2.d0*h*o - e*f - c*b) * six/two  
    aa(3, 5) = (2.d0*g*h - e*d - a*b) * six/two 
!
!  d(yz)
!
    aa(4, 1) = d*g - e*h  
    aa(4, 2) = d*o + f*g  
    aa(4, 3) = (2.d0*o*f - e*h - d*g) * six/two  
    aa(4, 4) = e*o + f*h  
    aa(4, 5) = d*h + e*g  
!
!  d(xy)
!
    aa(5, 1) = a*d - b*e  
    aa(5, 2) = a*f + c*d  
    aa(5, 3) = (2.d0*c*f - e*b - a*d) * six/two  
    aa(5, 4) = b*f + c*e  
    aa(5, 5) = a*e + d*b
! 
!  Use the following block for diagnostics only
!
    if (.false. .and. loop > 1) then
      write(*,*)"  atxyz"
      write(*,'(3f16.10)')((atxyz(i,j), i = 1,3),j = 1,3)
      write(*,*)
      do i = 1, 3
        do j = 1, 3
          sum = 0.d0
          sum1 = 0.d0
          do k = 1, 3
            sum = sum + atxyz(i,k)*atxyz(j,k)
            sum1 = sum1 + atxyz(k,i)*atxyz(k,j)
          end do
          e1(i,j) = sum
          e2(i,j) = sum1
        end do
      end do
      write(*,*) " atxyz squared"
      write(*,'(3f16.10)')((e1(i,j), i = 1,3),j = 1,3)
      write(*,*)
      write(*,'(3f16.10)')((e2(i,j), i = 1,3),j = 1,3)
      write(*,*)  
      do i = 1, 5
        do j = 1, 5
          sum = 0.d0
          sum1 = 0.d0
          do k = 1, 5
            sum = sum + aa(i,k)*aa(j,k)
            sum1 = sum1 + aa(k,i)*aa(k,j)
          end do
          e1(i,j) = sum
          e2(i,j) = sum1
        end do
      end do
      write(*,*) " d_rot "
      write(*,'(5f16.10)')((aa(i,j), i = 1, 5),j = 1, 5)
      write(*,*)
      write(*,*) " d_rot squared"
      write(*,'(5f16.10)')((e1(i,j), i = 1, 5),j = 1, 5)
      write(*,*)
      write(*,'(5f16.10)')((e2(i,j), i = 1, 5),j = 1, 5)
      i = i
    end if
!
! Put the "d" transform into t1
!
  l = 0
  do i = 5, 9
    l = l + 1
    k = 0
    do j = 5, 9
      k = k + 1
      t1(i, j) = aa(l,k)
    end do
  end do
!
! Put the "p" transform into t1
!
  do i = 1, 3
    l = i + 1
    do j = 1, 3
      k = j + 1
      t1(k, l) = atxyz(i, j)
    end do
  end do
  t1(1, 1) = 1.d0
end subroutine rotsec
