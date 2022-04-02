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

      subroutine molsym(coord, ierror, r)
      use common_arrays_C, only : nat, nbonds, ibonds
      use symmetry_C, only : ielem, cub, name
      use chanel_C, only : iw
      use molkst_C, only : numat, keywrd, moperr
      implicit none
      integer , intent(out) :: ierror
      double precision  :: coord(3,numat)
      double precision  :: r(3,3)
      integer , dimension(6) :: icyc
      integer :: i, j, k, l, ij, iturn, iqual, kndex, icheck, naxes, iz, ix, iy
      double precision, dimension(6) :: f
      double precision, dimension(3) :: ew, help
      double precision, dimension(3,3) :: rhelp
      double precision, dimension(3) :: shift
      double precision :: toler, wmol, sum, rxy, tole, distxy, rmin, sina, &
        cosa, theta, total
      logical :: linear, cubic, axis, sphere, debug, reorie, iscube, used(10)
!***************************************************************
!                                                              *
!     Credits not yet written                                  *
!                                                              *
!***************************************************************
!
!    MOLSYM CALCULATES THE MOLECULAR SYMMETRY AS A POINT-GROUP SYMBOL
!
      debug = index(keywrd,' MOLSYM') /= 0
!
!   REORIE IS NORMALLY .TRUE.   IF SET .FALSE., THEN THE ORIENTATION
!   SUPPLIED WILL BE USED IN DETERMINING SYMMETRY
!
      ew = 0.d0
      reorie = index(keywrd,' NOREOR') == 0
      toler = 0.1D0
      ierror = 0
      kndex = 0
      name = '????'
      do i = 1, 3
        cub(i,:) = 0.D0
        cub(i,i) = 1.D0
      end do
!
!  Now to put all the symmetry operations into the array ELEM
!   In order, the operations in ELEM are
!
!   1 C2(X)       6 Sigma(YZ)   11 C6    16 S8
!   2 C2(Y)       7 inversion   12 C7    17 S10
!   3 C2(Z)       8 C3          13 C8    18 S12
!   4 Sigma(XY)   9 C4          14 S4    19 1 if cubic, 0 otherwise.
!   5 Sigma(XZ)  10 C5          15 S6    20 1 if infinite, 0 otherwise.
!
!
      do i = 1, 18
        call bldsym (i, i)
        ielem(i) = 0
      end do
      ielem(19) = 0
      ielem(20) = 0
!
!   Calculate Moments of Inertia and Axes of Inertia.
!   First, center the molecule.
!
      shift = 0.D0
      wmol = 0.D0
      do i = 1, numat
        wmol = wmol + nat(i)
        shift = shift + nat(i)*coord(:,i)
      end do
      ij = 0
      do i = 1, 3
        shift(i) = shift(i)/wmol
        coord(i,:numat) = coord(i,:numat) - shift(i)
        do j = 1, i
          ij = ij + 1
          f(ij) = ij*1.D-8
          do k = 1, numat
            f(ij) = f(ij) + nat(k)*coord(i,k)*coord(j,k)
          end do
        end do
      end do
      iscube = .true.
      do
      if (.not.reorie) then
!
!     USE ORIENTATION AS SUPPLIED
!
!     FOR SYMMETRY, ORIENTATION IS: Z=PRINCIPAL AXIS, X,Y=SECONDARY
!
!      WRITE(IW,*)'H IN MOLSYM'
        sum = 0.D0
        do i = 1, 5
          sum = sum + abs(f(i))
        end do
        linear = sum < 0.01D0
        sphere = linear .and. abs(f(6))<0.01D0
        cubic = abs(f(1)-f(3))<0.01D0 .and. &
                abs(f(1)-f(6))<0.01D0 .and. &
                abs(f(2))+abs(f(4))+abs(f(5))<0.01D0
        do i = 1, 3
          r(i,:) = 0.D0
          r(i,i) = 1.D0
        end do
      else
        call rsp (f, 3, ew, r)
        sum = r(1,1)*(r(2,2)*r(3,3)-r(3,2)*r(2,3)) + &
              r(1,2)*(r(2,3)*r(3,1)-r(2,1)*r(3,3)) + &
              r(1,3)*(r(2,1)*r(3,2)-r(2,2)*r(3,1))
        if (sum > 1.d0) then
          sum = 1.01D0
          do j = 1, 3
            if (r(j,j) >= sum) cycle
            sum = r(j,j)
            i = j
          end do
          r(:,i) = -r(:,i)
        end if
        r(1,3) = r(2,1)*r(3,2) - r(3,1)*r(2,2)
        r(2,3) = r(3,1)*r(1,2) - r(1,1)*r(3,2)
        r(3,3) = r(1,1)*r(2,2) - r(2,1)*r(1,2)
!
!   Determine if molecule belongs to a special group
!
        linear = ew(2) < 1.D-2
        sphere = ew(3) < 1.D-2
        cubic = ew(3) - ew(1) < 5.D-3*max(ew(3),40.D0)
      end if
      if (sphere) then
!
!   Set flags 8 and 12 simultaneously to 1 - impossible for
!   non-spherical systems.
!
        ielem(7) = 1
        ielem(8) = 1
        ielem(10) = 1
        ielem(20) = 1
!
!  Make sure that C5 axis is different from C3 - this eliminates
!  accidental degeneracy.
!
        cub = 0.D0
        cub(1,2) = 1.D0
        cub(2,3) = 1.D0
        cub(3,1) = 1.D0
        go to 270
      else if (linear) then
!
!   Molecule is C-infinity-v or D-infinity-h
!
        call symopr (numat, coord, 1, r)
        ielem(20) = 1
        go to 250
      end if
      if (reorie) then
        if (.not.cubic .and. ew(3) - ew(2) < 1.D-2*ew(3)) then
!
!   Molecule has two-fold degeneracy ( Cn, Cnv, Dn, Dnh, Sn, n>2)
!
          do i = 1, 3
            rxy = -r(i,1)
            r(i,1) = r(i,3)
            r(i,3) = rxy
          end do
          rxy = ew(1)
          ew(1) = ew(3)
          ew(3) = rxy
        end if
        axis = abs(ew(1) - ew(2)) < 0.01D0*ew(2)
      else
        axis = abs(f(1) - f(3)) < 0.01D0*f(6)
      end if
!
!   Is there a plane of symmetry perpendicular to the Z axis?
!
      call symopr (numat, coord, 1, r)
      if ( .not. cubic) go to 1100
        ielem(19) = 1
        call plato (coord, r, cubic)
        if (moperr) return
      if ( .not. iscube .or. cubic) go to 1000
         !
         !  Three moments of inertia were equal, but the system does
         !  not have cubic symmetry.  Use an additional test to check
         !  for cubic nature.
         !
        do i = 1, 6
          f(i) = 0.d0
        end do
        ij = 0
        do i = 1, 3
          do j = 1, i
            ij = ij + 1
            f(ij) = ij * 1.d-8
            do k = 1, numat
              sum = 0.d0
              do l = 1, numat
                sum = sum + Sqrt ((coord(1, k)-coord(1, l))**2 + &
                                  (coord(2, k)-coord(2, l))**2 + &
                                  (coord(3, k)-coord(3, l))**2)
              end do
              f(ij) = f(ij) + sum * coord(i, k) * coord(j, k)
            end do
          end do
        end do
         !
         !  All species of this type (pseudo-cubic) have similar
         !  moments of inertia.  Increase the relative difference of
         !  the moments.
        f(1) = f(1) - 0.8d0 * f(1)
        f(3) = f(3) - 0.8d0 * f(1)
        f(6) = f(6) - 0.8d0 * f(1)
         !
         !  Reset IELEM(19) - system is NOT cubic!
         !
        ielem(19) = 0
        iscube = .false.
    end do
   !
   !   Set flags 8 and 12 simultaneously to 1 - impossible for
   !   non-spherical systems.
   !
    ielem(7) = 1
    ielem(8) = 1
    ielem(10) = 1
    ielem(20) = 1
   !
   !  Make sure that C5 axis is different from C3 - this eliminates
   !  accidental degeneracy.
   !
    do i = 1, 3
      do j = 1, 3
        cub(i, j) = 0.d0
      end do
    end do
    cub(1, 2) = 1.d0
    cub(2, 3) = 1.d0
    cub(3, 1) = 1.d0
    go to 270
1000 if (moperr) return
1100 if (axis) then
!
!  Molecule has degeneracy.  At this point, the Z-axis is
!  defined.  The molecule is Cn Cnv Cnh Dn, etc, or a special group,
!  e.g. Td, Oh, Ih.
!
  130   continue
        iturn = 7
        j = 0
!
!   Check for the existance of a Cn(2<n<9) or Sn(1<0.5*n<7)
!
        do i = 8, 18
!
!   Make certain that TOLE gets tighter as the n in Cn increases.
!
          if (i < 14) then
            tole = toler*9.D0/(i - 5)**2
          else
            tole = toler*16.D0/(2*i - 24)**2
          end if
          call chi (tole, coord, i, iqual)
!
!   If Cn, then set ITURN = n  (Remember ELEM(8) is C3, therefore offset
!                               by 5)
!
          if (ielem(i)/=1 .or. i>=14) cycle
          if (iturn > 9) then
            toler = toler*0.5D0
            go to 130
          end if
          iturn = i
        end do
        if (ielem(14) + ielem(15) + ielem(17) > 1 .or. &
            ielem(15) + ielem(16) + ielem(17) > 1 .or. &
            ielem(16) + ielem(17) + ielem(18) > 1) then
          toler = toler*0.5D0
          go to 130
        end if
        iturn = iturn - 5
        if (debug) then
          write (iw, '(A)') ' after checking 8-18'
          write (iw, '(20I3)') ielem
        end if
!
!  Now use two adjacent equivalent atoms, not on the
!  Z axis, to define the X-axis.
!
        do i = 1, numat
          distxy = coord(1,i)**2 + coord(2,i)**2
          if (distxy < toler) cycle
!
!   Atom I is the first atom.
!
          rmin = 1000.D0
          kndex = 0
          do j = i + 1, numat
            if (abs(abs(coord(3,i)) - abs(coord(3,j))) > 0.2D0) cycle
            rxy = coord(1,j)**2 + coord(2,j)**2
            if (abs(rxy - distxy) > toler .or. nat(i) /= nat(j)) cycle
!
!   Atoms i and j are equidistant from the Z axis, and are of the same type.
!   Now check that the atoms connected to them are also identical.
!
            if (nbonds(i) /= nbonds(j)) cycle
            used(:nbonds(i)) = .false.
            do k = 1, nbonds(i)
              do l = 1, nbonds(j)
                if (nat(ibonds(k,i)) == nat(ibonds(l,j)) .and. .not. used(l)) exit
              end do
              if (l <= nbonds(i)) used(l) = .true.
            end do
            l = 0
            do k = 1, nbonds(i)
              if (used(k)) l = l + 1
            end do
            if (l < nbonds(i)) cycle
            rxy = (coord(1,i) - coord(1,j))**2 + (coord(2,i) - coord(2,j))**2
            if (rxy > rmin) cycle
            kndex = j
            rmin = rxy
          end do
!
!   Atom KNDEX is the second, adjacent, atom, equivalent to I.
!
          exit
        end do
        if (kndex < 1) then
!
!  System does not have a Cn axis!  Go back and treat it as an Abelian
!  system.
!
          axis = .FALSE.
          go to 190
        end if
        help(1) = coord(1,i) + coord(1,kndex)
        help(2) = coord(2,i) + coord(2,kndex)
        distxy = sqrt(help(1)**2+help(2)**2)
        sina = help(2)/distxy
        cosa = help(1)/distxy
        call rotmol (numat, coord, sina, cosa, 1, 2, r)
!
!   Is there a Sigma(XZ) plane of symmetry?
!
        call chi (toler, coord, 5, iqual)
        if (ielem(5) /= 1) then
!
!    Is there a C2(X) axis of rotation?
!
          call chi (toler, coord, 1, iqual)
          if (ielem(1) /= 0) then
!
!    Check for an improper axis of rotation.
!
            theta = 1.5707963268D0/dble(iturn)
            sina = sin(theta)
            cosa = cos(theta)
            icheck = 0
  180       continue
            call rotmol (numat, coord, sina, cosa, 1, 2, r)
            if (icheck > 0) then
              axis = (iturn /= 2)
              go to 190
            end if
            call chi (toler, coord, 5, iqual)
            if (ielem(5) > 0) go to 190
            icheck = 1
            sina = -sina
            go to 180
          end if
        end if
      end if
  190 continue
      if (cubic) call orient (numat, coord, r)
      if (debug) then
        write (iw, '(A)') ' C2 and sigma-v     '
        write (iw, '(20I3)') ielem
      end if
      if (.not.axis) then
        toler = 0.2D0
!
!   Molecule belongs to one of the 8 Abelian groups
!   (C1, C2, Ci, Cs, C2v, D2, C2h, or D2h)
!
        do i = 1, 6
          call chi (toler, coord, i, iqual)
          icyc(i) = (1 + iqual)*ielem(i)
        end do
        if (reorie) then
          naxes = ielem(1) + ielem(2) + ielem(3)
!
!    Special handling for the Abelian Groups:  Determine
!    the principal axis by atom counts.
!
          if (naxes <= 1) then
            iz = 1
            if (ielem(1) == 1) go to 220
            iz = 2
            if (ielem(2) == 1) go to 220
            iz = 3
            if (ielem(3) == 1) go to 220
            if (icyc(5) > icyc(4)) iz = 2
            if (icyc(6) > icyc(7-iz)) iz = 1
            go to 220
          end if
          iz = 1
          if (icyc(2) > icyc(1)) iz = 2
          if (icyc(3) > icyc(iz)) iz = 3
  220     continue
          icyc(7-iz) = -1
          ix = 1
          if (icyc(5) > icyc(6)) ix = 2
          if (icyc(4) > icyc(7-ix)) ix = 3
          iy = 6 - ix - iz
!
!  Whew!   Now to re-orient the molecule so the the principal
!  axis is 'Z'
!
          rhelp(:,1) = r(:,ix)
          rhelp(:,2) = r(:,iy)
          rhelp(1,3) = r(2,ix)*r(3,iy) - r(3,ix)*r(2,iy)
          rhelp(2,3) = r(3,ix)*r(1,iy) - r(1,ix)*r(3,iy)
          rhelp(3,3) = r(1,ix)*r(2,iy) - r(2,ix)*r(1,iy)
          call symopr (numat, coord, -1, r)
          r = rhelp
          call symopr (numat, coord, 1, r)
        end if
!
!   And re-calculate the first 7 Characters.
!   (C2(X), C2(Y), C2(Z), Sigma(XY), Sigma(XZ), Sigma(YZ), i)
!
      end if
  250 continue
      do i = 1, 7
        call chi (toler, coord, i, iqual)
      end do
      if (debug) then
        write (iw, '(A)') ' After re-doing 1-7'
        write (iw, '(20I3)') ielem
      end if
  270 continue
      call symopr (numat, coord, -1, r)
      total = ew(1) + ew(2) + ew(3)
      ew = total - ew
      call cartab
      return
      end subroutine molsym
      subroutine chi(toler, coord, ioper, iqual)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use symmetry_C, only : ielem, elem, jelem
      use molkst_C, only : numat
      use common_arrays_C, only : nat
!***********************************************************************
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ioper
      integer , intent(out) :: iqual
      double precision , intent(in) :: toler
      double precision , intent(in) :: coord(3,numat)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iresul, i, j
      double precision, dimension(3) :: help
!-----------------------------------------------
!***********************************************************************
!
!   CHI RETURNS A '1' IN IELEM(IOPER) IF THE SYMMETRY OPERATION ELEM(IOP
!       LEAVES THE SYSTEM UNCHANGED.  OTHERWISE CHI RETURNS '0' IN IELEM
!
!   ON INPUT: COORD        = CARTESIAN COORDINATES
!             IOPER        = SYMMETRY OPERATION TO BE PERFORMED
!             ELEM         = SYMMETRY OPERATORS AS 3*3 MATRICES
!   ON OUTPUT IELEM(IOPER) = 1 OR 0.
!
!***********************************************************************
      iresul = 1
      iqual = 0
      l20: do i = 1, numat
        help(1) = coord(1,i)*elem(1,1,ioper) + &
                  coord(2,i)*elem(1,2,ioper) + &
                  coord(3,i)*elem(1,3,ioper)
        help(2) = coord(1,i)*elem(2,1,ioper) + &
                  coord(2,i)*elem(2,2,ioper) + &
                  coord(3,i)*elem(2,3,ioper)
        help(3) = coord(1,i)*elem(3,1,ioper) + &
                  coord(2,i)*elem(3,2,ioper) + &
                  coord(3,i)*elem(3,3,ioper)
        do j = 1, numat
          if (nat(i) /= nat(j)) cycle
          if (abs(coord(1,j)-help(1)) > toler) cycle
          if (abs(coord(2,j)-help(2)) > toler) cycle
          if (abs(coord(3,j)-help(3)) > toler) cycle
          jelem(ioper,i) = j
          if (i == j) iqual = iqual + 1
          cycle  l20
        end do
        iresul = 0
      end do l20
      ielem(ioper) = iresul
      return
      end subroutine chi
      subroutine makopr(numat, coord, ierror, r)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use symmetry_C, only : ielem, jy, nclass
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: numat
      integer , intent(out) :: ierror
      double precision  :: coord(3,numat)
      double precision  :: r(3,3)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, iqual
      double precision :: toler
!-----------------------------------------------
!*********************************************************************
!
!   MAKOPR builds the operations based on the point group
!          of the system.  A check is made to verify that the
!          operations are valid.
!
!*********************************************************************
      call symopr (numat, coord, 1, r)
      if (nclass < 2) return
!
!   NCLASS is the number of Classes in the Group.
!   Construct the Operations corresponding to the Classes
!   These are stored in ELEM.
!
      do i = 2, nclass
        call bldsym (jy(i), i)
      end do
!
!   Use a more tolerant criterion for recognizing operations because
!   the point-group has already been identified.
!
      toler = 0.2D0
      do i = 2, nclass
        call chi (toler, coord, i, iqual)
        if (ielem(i) >= 1) cycle
        ierror = 5
      end do
      call symopr (numat, coord, -1, r)
      if (iqual < 0) return ! dummy use of iqual
      return
      end subroutine makopr
      subroutine orient(numat, coord, r)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use symmetry_C, only : cub, ielem
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: numat
      double precision  :: coord(3,numat)
      double precision  :: r(3,3)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, jota, iqual, j
      double precision, dimension(2) :: wink
      double precision :: toler, wink2, sina, cosa, sinb, cosb

      save wink, toler
!-----------------------------------------------
!***********************************************************************
!
!    ORIENT is part of the SYMMETRY package
!           In cubic systems the three moments of inertia are identical.
!           Therefore, use the nearest atom to the center to determine
!           the principal axis (Z)
!
!**********************************************************************
      data toler/ 0.1D0/
      data wink(1), wink(2)/ 0.955316618125D0, 0.65235813978437D0/
!
!   WINK(1) rotates a cubic molecule through half the tetrahedral angle.
!   WINK(2) does something similar for an icosahedral system.
!
!   WINK(1)=ACOS(SQRT(1/3))
!   WINK(2)=ACOS((47+21*SQRT(5))/(75+33*SQRT(5)))
!           = angle from a vertex of an icosahedron to the center and th
!             to the middle of the triangular face
!           = 37.377366 degrees
!
!
      wink2 = 0.D0
      if (ielem(8) >= 1) then
!
!   Check for a S4 or C5 axis
!
        do i = 1, 2
          jota = 18 - 4*i
          wink2 = wink(i)
          sina = sin(wink2)
          cosa = cos(wink2)
          call rotmol (numat, coord, sina, cosa, 1, 3, r)
          call chi (toler, coord, jota, iqual)
          if (ielem(jota) > 0) exit
          if (i == 1) then
            call chi (toler, coord, 3, iqual)
            if (ielem(3) == 1) exit
          end if
          wink2 = -wink2
          sinb = sin(2.D0*wink2)
          cosb = cos(2.D0*wink2)
          call rotmol (numat, coord, sinb, cosb, 1, 3, r)
          call chi (toler, coord, jota, iqual)
          if (ielem(jota) > 0) exit
          if (i == 1) then
            call chi (toler, coord, 3, iqual)
            if (ielem(3) == 1) exit
          end if
          call rotmol (numat, coord, sina, cosa, 1, 3, r)
        end do
        call chi (toler, coord, 9, iqual)
!
!   Check on all IELEM registers
        if (ielem(10) > 0) call chi (toler, coord, 17, iqual)
      else
!
!   No C3 axis, therefore not T, Td, Th, O, or Oh.
!
        wink2 = -wink(1)
        if (ielem(10) > 0) wink2 = -wink(2)
        sina = -sin(wink2)
        cosa = cos(wink2)
        call rotmol (numat, coord, sina, cosa, 1, 3, r)
        call chi (toler, coord, 8, iqual)
        call rotmol (numat, coord, (-sina), cosa, 1, 3, r)
        if (ielem(8) <= 0) then
          if (ielem(9) <= 0) then
            wink2 = -wink2
          else
            call rotmol (numat, coord, 0.707106781186D0, 0.707106781186D0, 1, 2, r)
          end if
        end if
      end if
      j = sum(ielem(:17))
      if (j==2 .and. ielem(1)+ielem(8)==2) return
      cub(1,1) = cos(wink2)
      cub(3,3) = cub(1,1)
      cub(1,3) = sin(wink2)
      cub(3,1) = -cub(1,3)
      call mult33 (cub, 8)
      call mult33 (cub, 15)
      call chi (toler, coord, 8, iqual)
      j = iqual ! dummy use of iqual
      call chi (toler, coord, 15, iqual)
      if (j + iqual == 0) return ! dummy use of iqual
      return
      end subroutine orient
      subroutine plato(coord, r, cubic)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use chanel_C, only : iw
      use molkst_C, only : numat
      use common_arrays_C, only : nat
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, intent(out) :: cubic
      double precision  :: coord(3,numat)
      double precision  :: r(3,3)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(numat) :: near
      integer , dimension(3) :: ipoly
      integer ::  i, l, j, ii, k, i1, j1, m, i3, i2, ligand_type, &
      store_ligand, j2, k1, k2
      double precision, dimension(3) :: allr
      double precision, dimension(3,3) :: xyz
      double precision :: toler = 0.1d0, xmin, dist, r2j, angle, sum, buff, buff1, &
      rmin, ymin, dist1, dist2

      save toler
!-----------------------------------------------
!***********************************************************************
!
!    PLATO  is part of the SYMMETRY package.  It generates the
!           unitary transform which will orientate the system
!           so that the 'Z' axis is an axis of rotation.
!
!    PLATO is called when the system belongs to a cubic point group.
!          Atoms nearest to the center of symmetry are identified.
!          If there are 4, 6, 8, 12, or 20, then they outline one
!          of the Platonic solids (tetrahedron, octahedron, cube,
!          icosahedron, or pentagonal dodecahedron).  In that case,
!          each atom lies on a 3, 4, or 5 fold symmetry axis.
!
!          In other cases the axis of symmetry lies goes through a
!          regular polygon, either an equilateral triangle, a square,
!          or a regular pentagon.
!
!          The unitary transform is returned through R.
!
!**********************************************************************
      call symopr (numat, coord,-1, r)
      rmin = 0.1d0
      xmin = 0.d0
      ymin = 0.d0
      store_ligand = 0
      ligand_type = 0
      k = 100
      do
  !
  !   Find the smallest number of atoms that are at the same distance
  !   from the center.  This is done by steadily moving out, rmin, from the center,
  !   and counting the number of atoms in the shell.  The minimum number of atoms, k,
  !   is at a distance of ymin.
  !
      do i = 1, numat
        xmin = coord(1, i) ** 2 + coord(2, i) ** 2 + coord(3, i) ** 2
        if (xmin >= rmin) then
          ligand_type = nat(i)
          exit
        end if
      end do
      if (xmin < rmin) then
        xmin = ymin
        exit
      end if
  !
  !  How many atoms are at the same distance from the center?
  !
      l = 0
      do i = 1, numat
        dist = coord(1, i) ** 2 + coord(2, i) ** 2 + coord(3, i) ** 2
        if (Abs(dist - xmin) < toler .and. nat(i) == ligand_type) l = l + 1
      end do
  !
  !  The current "shell" has less than an earlier shell, therefore use it.
  !
      if (l < k) then
        k = l
        ymin = xmin
        store_ligand = ligand_type
      end if
  !
  !  l cannot be less than 4 (a tetrahedron)
  !
      if (l < 5) exit
  !
  ! Increase the value of rmin to start the check of the next shell
  !
      rmin = xmin + 2.d0*toler
      end do
  !
  !  How many atoms are at the same distance from the center?
  !
      ligand_type = store_ligand
      l = 0
      do i = 1, numat
        dist = coord(1, i) ** 2 + coord(2, i) ** 2 + coord(3, i) ** 2
        if (dist >= toler .and. dist <= xmin+toler) then
          if (Abs (dist-xmin) < toler .and. nat(i) == ligand_type) then
            l = l + 1
            near(l) = i
          end if
        end if
      end do
      if (l /= 0) then
        if (l == 1) then
  !
  !  System was identified as cubic because of moments of inertia only.
  !  This was a coincidence.  Correct error and continue
  !
          cubic = .false.
          return
        end if
        if (l == 4 .or. l == 6 .or. l == 8 .or. l == 12 .or. l == 20) then
           !
           !   How many near neighbors has atom NEAR(1) got?
           !
          xmin = 10000.d0
          j = near(1)
          do ii = 2, l
            i = near(ii)
            dist = (coord(1, j)-coord(1, i)) ** 2 + &
                   (coord(2, j)-coord(2, i)) ** 2 + &
                   (coord(3, j)-coord(3, i)) ** 2
            if (dist < xmin) then
              xmin = dist
            end if
          end do
          j1 = near(1)
          j2 = near(2)
          k1 = 0
          k2 = 0
          do ii = 1, l
            i = near(ii)
            dist1 = (coord(1, j1)-coord(1, i)) ** 2 + &
                    (coord(2, j1)-coord(2, i)) ** 2 + &
                    (coord(3, j1)-coord(3, i)) ** 2
            dist2 = (coord(1, j2)-coord(1, i)) ** 2 + &
                    (coord(2, j2)-coord(2, i)) ** 2 + &
                    (coord(3, j2)-coord(3, i)) ** 2
            if (Abs (dist2-xmin) < toler) then
              k2 = k2 + 1
            end if
            if (Abs (dist1-xmin) < toler) then
              k1 = k1 + 1
            end if
          end do
          if (k1 == k2) then
            k = k1
          else if (k1 > 2 .or. k2 > 2) then
              !
              !  There are no cubic solids that (a) have three or more equidistant
              !  vertices and (b) are not platonic
              !
            cubic = .false.
            return
          else
            k = 0
          end if   !   # points    #neighbors   Polyhedron
  !
          if (          l == 4  .and. k == 3 &  !   Tetrahedron
          .or.          l == 6  .and. k == 4 &  !   Cube
          .or.          l == 8  .and. k == 3 &  !   Octahedron
          .or.          l == 12 .and. k == 5 &  !   Icosahedron
          .or.          l == 20 .and. k == 3 &  !   Pentagonal dodecahedron
          ) then
              !
              !     The system is a Platonic solid.  This is the simplest case.
              !     An atom lies on a high-symmetry axis.
              !
            i = near(l)
            r(1, 3) = coord(1, i) / dist
            r(2, 3) = coord(2, i) / dist
            r(3, 3) = coord(3, i) / dist
            go to 1100
          end if
        end if
        !
        !  Pick one atom (here NEAR(1)), and find the three nearest
        !  neighbors in the same set.  Put these in IPOLY(1:3).
        !
        i1 = near(1)
        do i = 1, 3
          xmin = 10000.d0
          ipoly(i) = 0
          loop: do j1 = 2, l
            j = near(j1)
            do m = 1, i - 1
              if (ipoly(m) == j) cycle loop
            end do
            r2j = (coord(1, i1)-coord(1, j)) ** 2 + (coord(2, i1)-coord(2, j)) &
           & ** 2 + (coord(3, i1)-coord(3, j)) ** 2
            if (xmin > r2j) then
              xmin = r2j
              ipoly(i) = j
              allr(i) = r2j
              if (xmin < 0.01d0) then
                write (iw,*) " Geometry is apparently cubic, but the distance"
                write (iw, "(A,I5,A,I5,A,F12.8,A)") " between atoms", i1, " and",&
               &  j, " is", Sqrt (xmin), " Angstroms"
                write (iw, "(/,A,/)") " Here is the current geometry"
                call geout (iw)
                call mopend("Geometry is severely faulty")
                return
              end if
            end if
          end do loop
        end do
        do i = 1, 3
          allr(i) = Sqrt (allr(i))
        end do
        if (ipoly(1)*ipoly(2)*ipoly(3) == 0) then
  !
  !  System is obviously not cubic.
  !
          cubic = .false.
          return
        end if
        !
        !   Identify the two atoms in the regular polygon.
        !
        do i = 1, 2
          do j = i + 1, 3
            if (Abs (allr(i)-allr(j)) < 0.01d0) go to 1000
          end do
        end do
        !
        !  There are not many cubic solids in which a vertex nearest to
        !  the center has three neighbors all at different distances.
        !  Since this is very rare, relative to the number of truly non-
        !  cubic systems, ignore all such cases.
        !
        cubic = .false.
        return
1000  do k = 1, 3
          xyz(k, 3) = coord(k, ipoly(i))
          xyz(k, 2) = coord(k, ipoly(j))
          xyz(k, 1) = coord(k, i1)
        end do
        i3 = ipoly(i)
        i2 = ipoly(j)
        ipoly(1) = i3
        ipoly(2) = i2
        call bangle (xyz, 3, 1, 2, angle)
        if (Abs (angle-1.0472d0) < 0.1d0) then
           !
           !   It's a triangle!
           !
          do i = 1, 3
            r(i, 3) = coord(i, i1) + coord(i, i2) + coord(i, i3)
          end do
        else if (Abs (angle-1.5707963d0) < 0.1d0) then
           !
           !   It's a square!
           !
          do i = 1, 3
            r(i, 3) = coord(i, i2) + coord(i, i3)
          end do
        else if (Abs (angle-1.885d0) < 0.1d0) then
           !
           !  It's a pentagon!
           !
          do i = 1, 3
            r(i, 3) = 1.6180341d0 * (coord(i, i2)+coord(i, i3)) - coord(i, i1)
          end do
        else
           !
           !  There are no cubic solids that (a) have two or more equidistant
           !  vertices and (b) these vertices are not part of a regular
           !  polygon.
           !
          cubic = .false.
          return
        end if
1100    sum = Sqrt (r(1, 3)**2+r(2, 3)**2+r(3, 3)**2)
        do i = 1, 3
          r(i, 3) = r(i, 3) / sum
        end do
        buff = Sqrt (r(1, 3)**2+r(2, 3)**2)
        buff1 = Sqrt (r(1, 3)**2+r(3, 3)**2)
        if (buff <= buff1) then
          r(1, 1) = r(3, 3) / buff1
          r(2, 1) = 0.d0
          r(3, 1) = -r(1, 3) / buff1
        else
          r(1, 1) = r(2, 3) / buff
          r(2, 1) = -r(1, 3) / buff
          r(3, 1) = 0.d0
        end if
        r(1, 2) = r(2, 3) * r(3, 1) - r(2, 1) * r(3, 3)
        r(2, 2) = r(3, 3) * r(1, 1) - r(3, 1) * r(1, 3)
        r(3, 2) = r(1, 3) * r(2, 1) - r(1, 1) * r(2, 3)
        !
        !   The molecule is now orientated along one of the high-symmetry axes
        !
        call symopr (numat, coord, 1, r)
        return
      end if
      return
      end subroutine plato
      subroutine cartab
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use symmetry_C, only : nclass, igroup, group, allrep, &
      & nallop,ntbs, ntab, nallg, name, nirred, ielem, jx, jy
      use molkst_C, only : numcal, keywrd
      use chanel_C, only : iw
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: ngps = 57
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(ngps) :: nope, nrep
      integer :: icalcn, i, j, k, l, nopers, nreprs, nstabl, kl, ku, istart, &
        igp, nzz, nz, prnt
      double precision :: buff, fz, fn
      logical :: debug, first, large
      character, dimension(20) :: class*9
      logical, external :: symdec

      save nrep, class, debug, first, icalcn, large, nope
!-----------------------------------------------
!***********************************************************************
!
! CARTAB constructs Character Tables for Point Groups.  Each Group is
!        defined by an array having the same name as the Group, e.g. C2v
!        The format of each group is explained in BLOCK TABLES
!
!
!   The Magic Number is the decimal representation of a 20-digit binary
!   number, each digit of which is '1' if the associated operation is
!   present in the system, '0' otherwise.  The 20 operations are, in
!   order
!
!  1 C2(X)         6 Sigma(YZ)    11 C6     16 S8
!  2 C2(Y)         7 inversion    12 C7     17 S10
!  3 C2(Z)         8 C3           13 C8     18 S12
!  4 Sigma(XY)     9 C4           14 S4     19 Cubic group
!  5 Sigma(XZ)    10 C5           15 S6     20 Infinite group.
!
!***********************************************************************
      data first/ .TRUE./
      data prnt/ 1/
      data class/ 'C2(x)', 'C2(y)', 'C2(z)', 'Sigma(XY)', 'Sigma(v)', &
        'Sigma(d)', 'Inversion', ' C3', ' C4', ' C5', ' C6', ' C7', ' C8', &
        ' S4', ' S6', ' S8', ' S10', ' S12', ' ?? ', 'C(Inf)'/
      data icalcn/ 0/
!
      if (numcal /= icalcn) then
        icalcn = numcal
        debug = index(keywrd,'CARTAB') /= 0
        large = index(keywrd,'LARGE')/=0 .and. debug
      end if
      large = large .and. prnt==1
      if (large) prnt = 0
      if (first) then
        first = .FALSE.
!
!   Generate array to hold addresses of groups
!
!  NOPE will apply to operations
!  NREP will apply to group name and irreducible representations
!
        nope(1) = 1
        nrep(1) = 1
        do i = 2, ngps
          nope(i) = nope(i-1) + nallop(nope(i-1)) + 4
          nrep(i) = nrep(i-1) + nallop(nope(i-1)+1) + 1
        end do
!
!   Generate array to hold addresses of character tables
!
        j = 1
        do i = 1, ntbs
          k = ntab(i)
          ntab(i) = j
          j = j + k
        end do
!
!  Use the following lines to check that the common blocks
!  have the correct sizes.
!
!#      WRITE(*,*)' Size of Groups:',NOPE(NGPS)+NALLOP(NOPE(NGPS))+3
!#      WRITE(*,*)' Size of Reps  :',NREP(NGPS)+NALLOP(NOPE(NGPS)+1)
!#      WRITE(*,*)' Size of TABLES:',J-1
      end if
!
!  Debug code:  Print all the point groups, in order
!
      if (large) then
        do i = 1, ngps
          k = nope(i)
          l = nrep(i)
          nopers = nallop(k)
          nreprs = nallop(k+1)
          nstabl = ntab(nallop(k+2))
          write (iw, *) '         Point Group: ', allrep(l)
          write (iw, "(a,i4,a,i4,a,i4,a,i9)") ' No. Ops:', nopers + 1, ' No. Reps:', &
          nreprs, ' Table:', nallop(k+2), ' "Magic" No.:', nallop(k+3)
          write (iw, '(/10X,8(A,1X))') '   Identity   ', &
          (class(nallop(j)),j=k + 4,k + 3 + nallop(k))
          write (iw, '(A,10I10)') '    '//allrep(l+1), (1,j=1,nallop(k) + 1)
          do j = 2, nreprs
            kl = nstabl + (j - 2)*(nopers + 1)
            ku = kl + nopers
            write (iw, '(A,10I10)') '    '//allrep(l+j), (nallg(k),k=kl,ku)
          end do
        end do
      end if
!
!    Identify the Point Group of the molecule.
!
      do igroup = ngps, 1, -1
        if (.not.symdec(nallop(nope(igroup)+3),ielem)) cycle
        exit
      end do
      if (index(keywrd,' NOSYM') /= 0) igroup = 1
      istart = nope(igroup)
      igp = nrep(igroup)
      name = allrep(igp)
      nclass = nallop(istart) + 1
      nirred = nallop(istart+1)
!
!  IGROUP:  The number of the point-group of the system.
!  ISTART:  Starting address of the group.
!  IGP   :  Starting address of the names used in the group.
!  NAME  :  The name of the group.
!  NCLASS:  Number of operations used to represent the group, plus 1.
!  NIRRED:  Number of Irreducible Representations.
!  GROUP :  The full point-group character table for group IGROUP.
!
      group(1,:nclass) = 1.D0
      jy(2:nclass) = nallop(istart+4:nclass+2+istart)
      jx(:nirred) = allrep(igp+1:nirred+igp)
      istart = ntab(nallop(istart+2)) - 1
      do i = 2, nirred
        do k = 1, nclass
          istart = istart + 1
          buff = nallg(istart)
          if (buff>=10.D0 .and. &
          (name/='R3' .or. name=='R3' .and. k/=1 .and. k /= nclass)) then
            nzz = nallg(istart)
            nz = nzz/10
            fz = nz
            fn = nzz - 10*nz
            buff = 2.D0*cos(6.283185307179D0*fn/fz)
          end if
          group(i,k) = buff
        end do
      end do
      jy(1) = 0
      if (debug) then
        write (iw, '(2A)') ' Character Table for Group  ', name
        write (iw, '(7X,10A10)') '   E    ', (class(jy(k)),k=2,nclass)
        do i = 1, nirred
          write (iw, '(A4,10F10.4)') jx(i), (group(i,k),k=1,nclass)
        end do
      end if
      return
      end subroutine cartab
      logical function symdec (n1, ielem)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n1
      integer , intent(in) :: ielem(20)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: in1, i, ibin
!-----------------------------------------------
!
!*********************************************************************
!
!   SYMDEC matches up a set of symmetry operations for the system
!          with a point-group.  The set of symmetry operations are
!          stored in IELEM.  Point-groups are represented by a
!          number, N1, which is expanded into binary.
!
!          A system has the symmetry of a specific point-group if
!          every operation of the point-group is present in the
!          system.  Extra operations may also be present, but are
!          ignored.  This allows a controlled descent in symmetry.
!          The symmetry groups in the calling routine, CARTAB, are
!          stored in order of symmetry - C1 to R3.
!  SYMDEC  returns a .TRUE. if the system matches point-group N1.
!
!*********************************************************************
      in1 = n1
      do i = 1, 20
        ibin = mod(in1,2)
        if (ielem(i)/=1 .and. ibin==1) go to 20
        in1 = in1/2
      end do
      symdec = .TRUE.
      return
   20 continue
      symdec = .FALSE.
      return
      end function symdec
