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

      subroutine dimens(coord, iw)
!
!  Work out the three dimensions of a molecule.
!
!  Dimension 1: The largest distance between any two atoms.
!  Dimension 2: The largest distance in the plane perpendicular to Dimension 1, between any two atoms.
!  Dimension 3: The largest distance in the direction orthogonal to dimensions 1 and 2, between any two atoms.
!
      use molkst_C, only : numat
      use common_arrays_C, only : nat
      use elemts_C, only : elemnt
!***********************************************************************
      implicit none
      integer , intent(in) :: iw
      double precision  :: coord(3,numat)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(3,2) :: ij
      integer :: l, loop, j, k, kmax, lmax, i, kk, ll
      double precision, dimension(3,3) :: c
      double precision, dimension(3) :: dim
      double precision :: x1, y1, z1, rabmax, rmax, r, xy, ca, cb, sa, sb, ymin, &
        ymax, store_coords(3,numat), sum
!-----------------------------------------------
    if (numat == 1) return
    store_coords(:,:numat) = coord(:,:numat)
    rmax = 0.d0
    kmax = 0
    lmax = 0
    do k = 1, numat
      x1 = coord(1, k)
      y1 = coord(2, k)
      z1 = coord(3, k)
      do l = 1, k - 1
        r = (x1-coord(1, l)) ** 2 + (y1-coord(2, l)) ** 2 + (z1-coord(3, l)) ** 2
        if (r > rmax) then
          rmax = r
          kmax = k
          lmax = l
        end if
      end do
    end do
    k = kmax
    l = lmax
   !
   !   Determine vector joining most distant atoms.
   !
    x1 = coord(1, k) - coord(1, l)
    y1 = coord(2, k) - coord(2, l)
    z1 = coord(3, k) - coord(3, l)
   !
   !  Rotate the system so that the most distant atoms have the same
   !  "Z" coordinates.
   !
    xy = x1 ** 2 + y1 ** 2
    r = Sqrt (xy+z1**2)
    xy = Sqrt (xy)
    if (xy < 1.d-10) then
      if (z1 < 0.0d0) then
        ca = -1.d0
        cb = -1.d0
        sa = 0.d0
        sb = 0.d0
      else if (z1 > 0.0d0) then
        ca = 1.d0
        cb = 1.d0
        sa = 0.d0
        sb = 0.d0
      else
        ca = 0.d0
        cb = 0.d0
        sa = 0.d0
        sb = 0.d0
      end if
    else
      ca = x1 / xy
      cb = z1 / r
      sa = y1 / xy
      sb = xy / r
    end if
    c(1, 3) = ca * cb
    c(1, 2) = -sa
    c(1, 1) = ca * sb
    c(2, 3) = sa * cb
    c(2, 2) = ca
    c(2, 1) = sa * sb
    c(3, 3) = -sb
    c(3, 2) = 0.d0
    c(3, 1) = cb
    call symopr (numat, coord, 1, c)

    dim(1) = r
    ij(1, 1) = k
    ij(1, 2) = l
    dim_loop: do
     !
     !   The longest dimension is now "X"
     !   Find the most distant atom in "Y-Z" plane'
     !
      y1 = coord(2, l)
      z1 = coord(3, l)
      rabmax = 0.d0
      kk = 0
      ll = 0
      loop = 0
      do
        loop = loop + 1
        if (loop > 10) then
          k = kk
          l = ll
          exit
        else
         !
         !    Find the atom most distant from atom "L"
         !
          rmax = 0.d0
          do j = 1, numat
            r = (y1-coord(2, j)) ** 2 + (z1-coord(3, j)) ** 2
            if (r > rmax) then
              rmax = r
              k = j
            end if
          end do
         !
         !   Atom K is most distant from L in the Y-Z plane.
         !
          if (Abs (rmax-rabmax) < 1.d-5) exit
          if (rmax > rabmax) then
            rabmax = rmax
            kk = k
            ll = l
          end if
         !
         !   We don't know if RMAX is the maximum, therefore go to the
         !   mid-point, and repeat the test.
         !
         !   Now find the mid-point between the two atoms
         !
          y1 = 0.5d0 * (y1+coord(2, k))
          z1 = 0.5d0 * (z1+coord(3, k))
         !
         !    Find the atom most distant
         !
          rmax = 0.d0
          do j = 1, numat
            r = (y1-coord(2, j)) ** 2 + (z1-coord(3, j)) ** 2
            if (r > rmax) then
              rmax = r
              l = j
            end if
          end do
          y1 = coord(2, l)
          z1 = coord(3, l)
        end if
      end do
     !
     !   Determine vector joining most distant atoms (K and L).
     !
      y1 = coord(2, k) - coord(2, l)
      z1 = coord(3, k) - coord(3, l)
      r = Sqrt (y1**2 + z1**2 + 1.d-20)
      ca = y1 / r
      sa = z1 / r
     !
     !  Rotate system so that atoms K and L have the same Y coordinate.
     !
      do i = 1, numat
        sum         =  ca * coord(2, i) + sa * coord(3, i)
        coord(3, i) = -sa * coord(2, i) + ca * coord(3, i)
        coord(2, i) =  sum
      end do
      if (r > dim(1) + 1.d-10) then
        dim(1) = r
        ij(1, 1) = k
        ij(1, 2) = l
        cycle dim_loop
      else
        dim(2) = r
        ij(2, 1) = k
        ij(2, 2) = l
        exit dim_loop
      end if
    end do dim_loop
   !
   !   The longest dimension is now "X", the second longest dimension
   !   is "Y".
   !   Find the largest dimension in the "Z" direction
   !
    ymin = 1.d16
    ymax = -1.d16
    do i = 1, numat
      if (coord(3, i) > ymax) then
        k = i
        ymax = coord(3, i)
      end if
      if (coord(3, i) < ymin) then
        l = i
        ymin = coord(3, i)
      end if
    end do
    dim(3) = ymax - ymin
    ij(3, 1) = k
    if (k == l) then
      if (l == 1) then
        l = 2
      else
        l = 1
      end if
    end if
    ij(3, 2) = l
    write (iw, "(/,9X,A,/)") " MOLECULAR DIMENSIONS (Angstroms)"
    write (iw, "(9X,A)") "   Atom       Atom       Distance"
    write (iw, "(11X,A2,I6,3X,A2,I6,F12.5)") ((elemnt(nat(ij(j, i))), &
             & ij(j, i), i=1, 2), dim(j), j=1, 3)
    coord(:,:numat) =  store_coords(:,:numat)
end subroutine dimens
