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

subroutine prtgra ()
    use molkst_C, only: natoms, nvar, maxtxt, na1, keywrd
    use chanel_C, only: iw
    use elemts_C, only: elemnt
    use common_arrays_C, only : grad, xparam, txtatm, labels, loc, na
    implicit none
    character (len=30) :: blank
    character (len=100) :: fmat
    logical :: lint, lxyz
    integer :: i, j, k, l, m, n
    double precision :: glim, gmax, sum
    double precision, dimension (3, 2) :: tmp
    j = 0
    do i = 2, natoms
      if (na(i) /= 0) j = j + 1
    end do
    lxyz = (j == 0)
    lint = (j == natoms - 1)
    glim = 0.1d0 ! Only print large gradients
    i = Index (keywrd, " GRAD=") + Index (keywrd,"DERIV") + Index (keywrd," GRAD")
    if (i /= 0) glim = -0.1d0
    blank = " "
    if (maxtxt == 0) then
      i = 1
    else
      i = maxtxt + 3
    end if
    if (lxyz) then   !  All coordinates are Cartesian
      write (iw, '("       ", a, "     Cartesian Gradients         Cartesian Coordinates          |Gradient|")') &
        blank(:i)
      write (iw, '("  Atom ", a, "     X        Y        Z        X          Y          Z")') blank(:i)
    else if (lint) then  !  All coordinates are internal
      write (iw, '("  Atom ", a, "         Gradients                      Geometry    ")') blank(:i)
      write (iw, '("  Atom ", a, "    Bond    Angle  Dihedral    Bond       Angle    Dihedral")') blank(:i)
    else !  Coordinates are mixed - Cartesian and internal
      write (iw, '("  Atom ", a, "         Gradients                      Geometry                  |Grad|")')  &
        blank(:i)
    end if
    write (iw,*)
    l = 1
    do i = 1, natoms
      !
      !                       1         2         3         4
      !              12345678901234567890123456789012345678901234567
      !
      fmat = "(I5,A,"
      j = 0
      gmax = 0.d0
      if (l <= nvar) then
        k = 7
        n = 38
        do m = 1, 3
          if (loc(1, l) == i .and. loc(2, l) == m) then
            j = j + 1
            tmp(j, 1) = grad(l)
            tmp(j, 2) = xparam(l)
            if (na1 /= 99 .and. m > 1 .and. na(i) /= 0) then
              tmp(j, 2) = tmp(j, 2) * 180.d0 / 3.14159265358979d0
            end if
            gmax = Max (gmax, Abs (grad(l)))
            fmat (k:k+4) = "F9.3,"
            fmat (n:n+5) = "F11.5,"
            l = l + 1
            k = k + 5
            n = n + 6
          else
            fmat (k:k+12) = "'     -   ',"
            fmat (n:n+14) = "'     -     ',"
            k = k + 13
            n = n + 15
          end if
        end do
        fmat (n:n+8) = "f12.3,1X)"
      end if
      if (j /= 0 .and. glim < gmax) then
        if (maxtxt == 0) then
          if (j == 3 .and. na(i) == 0) then
            sum = sqrt(tmp(1, 1)**2 + tmp(2, 1)**2 + tmp(3, 1)**2)
            write (iw, fmat) i, " " // elemnt (labels(i)), ((tmp(m, n), m=1, j), n=1, 2), sum
          else
            write (iw, fmat) i, " " // elemnt (labels(i)), ((tmp(m, n), m=1, j), n=1, 2)
          end if
        else
         if (j == 3 .and. na(i) == 0) then
            sum = sqrt(tmp(1, 1)**2 + tmp(2, 1)**2 + tmp(3, 1)**2)
            write (iw, fmat) i, " " // elemnt (labels(i)) // "(" // txtatm (i) &
         & (:maxtxt) // ")" // blank (maxtxt+3:17), ((tmp(m, n), m=1, j), n=1, 2), sum
          else
            write (iw, fmat) i, " " // elemnt (labels(i)) // "(" // txtatm (i) &
         & (:maxtxt) // ")" // blank (maxtxt+3:17), ((tmp(m, n), m=1, j), n=1, 2)
          end if
        end if
      end if
    end do
    write(iw,*)
  end subroutine prtgra
  subroutine prt_sorted_gradients
    use molkst_C, only: natoms, nvar, maxtxt, pdb_label
    use chanel_C, only: iw
    use elemts_C, only: elemnt
    use common_arrays_C, only : grad, txtatm, loc, labels, coord
    implicit none
    character (len=30) :: blank
    integer :: i, j, k, l, m
    double precision :: glim, gmax, sum
    double precision, dimension (3) :: tmp
    double precision, allocatable :: mod_grad(:), grads(:,:)
    if (nvar < 3) return
    allocate (mod_grad(natoms), grads(3,natoms))
    mod_grad = 0.d0
    l = 1
    do i = 1, natoms
      j = 0
      if (l <= nvar) then
        tmp = 0.d0
        do m = 1, 3
          if (loc(1, l) == i .and. loc(2, l) == m) then
            j = j + 1
            tmp(m) = grad(l)
            l = l + 1
          end if
        end do
        mod_grad(i) = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
        grads(:,i) = tmp(:)
      end if
    end do
    if (pdb_label) then
      write(iw,'(/10x, a, //1x, a, 16x, a, 19x,a, 11x, a, 7x, a)')"LARGEST ATOMIC GRADIENTS","Atom","Label", &
      "Cartesian Gradients", "Cartesian Coordinates", " |Gradient|"
    else
      write(iw,'(/10x, a, //10x, a, 5x, a, 21x,a, 11x, a, 6x, a)')"LARGEST ATOMIC GRADIENTS","Atom","Label", &
      "Cartesian Gradients", "Cartesian Coordinates", " |Gradient|"
    end if
    write(iw,'(41x, a, 8x, a, /)')"  X          Y          Z", "X          Y          Z"
    gmax = 0.d0
    glim = 1.d6
    blank = " "
    do i = 1, min(20, natoms)
      sum = gmax
      do j = 1, natoms
        if (mod_grad(j) > sum .and. mod_grad(j) < glim) then
          sum = mod_grad(j)
          k = j
        end if
      end do
      if (sum < 1.d-3) exit
      if (pdb_label) then
      write (iw, '(i5, a, 3f11.3, 3f11.5, f11.3)') k, " " // elemnt (labels(k)) // "(" // txtatm (k) &
         & (:maxtxt) // ")" // blank (maxtxt+3:17), grads(:,k), coord(:,k), sum
      else
         write (iw, '(i14, 5x, a, 13x, 3f11.3, 3f11.5, f11.3)') k, " " // elemnt (labels(k)),grads(:,k), coord(:,k), sum
      end if
      glim = sum
      if (i > 9 .and. sum < 50) exit
    end do
    return
  end subroutine prt_sorted_gradients
