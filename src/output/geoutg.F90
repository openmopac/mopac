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

      subroutine geoutg(iprt)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use common_arrays_C, only : na, coord, geo, loc, labels, &
      & simbol, nb, nc, txtatm
      use molkst_C, only : natoms, nvar, ndep, maxtxt
      USE symmetry_C, ONLY: locpar, idepfn, locdep
      use chanel_C, only : iscr
      use elemts_C, only : elemnt
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
      integer , intent(in) :: iprt
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(3,natoms) :: igeo
      integer :: i, j, nopt, nbi, nci, l, nchars
      double precision :: degree, w, x
      character , dimension(3,natoms) :: line*14
      character , dimension(3) :: type
      character , dimension(3*natoms) :: optdat*14
      character :: blank*80
      save type
!-----------------------------------------------
!***********************************************************************
!
!   GEOUTG WRITES OUT THE GEOMETRY IN GAUSSIAN-8X STYLE
!
!***********************************************************************
      data type/ 'r', 'a', 'd'/
      degree = 57.29577951308232D0
      do i = 2, natoms
        if (na(i) == 0) exit
      end do
      if (i <= natoms) then
!
!  At least one atom is not in Cartesian coordinates
!  Convert everything to Cartesian coordinates
!
        coord(:,:natoms) = geo(:,:natoms)
        call xyzint (coord, natoms, na, nb, nc, 1.D0, geo)
        nvar = 0
        do i = 1, natoms
          do j = 1, min(3,i - 1)
            nvar = nvar + 1
            loc(1,nvar) = i
            loc(2,nvar) = j
          end do
        end do
      else
!
!   CONSTRAIN ANGLES TO THE RANGE 0 to 180 DEGREES
!
        do i = 4, natoms
          w = geo(2,i)*degree
          x = geo(3,i)*degree
          w = w - aint(w/360.D0)*360.D0
          if (w < 0) w = w + 360.D0
          if (w > 180.D0) then
            x = x + 180.D0
            w = 360.D0 - w
          end if
!
!  CONSTRAIN DIHEDRAL TO DOMAIN -180 - 180 DEGREES
!
          x = x - aint(x/360.D0 + sign(0.5D0 - 1.D-9,x) - 1.D-9)*360.D0
          geo(2,i) = w/degree
          geo(3,i) = x/degree
        end do
      end if
      igeo(:,:natoms) = -1
      do i = 1, nvar
        igeo(loc(2,i),loc(1,i)) = -2
      end do
      do i = 1, ndep
        if (idepfn(i) == 14) then
          igeo(3,locdep(i)) = -locpar(i)
        else
          if (idepfn(i) > 3) cycle
          igeo(idepfn(i),locdep(i)) = locpar(i)
        end if
      end do
      open(unit=iscr, status='SCRATCH', position='asis')
      nopt = 0
      do i = 1, natoms
        do j = 1, 3
          line(j,i) = ' '
          if (igeo(j,i) == (-1)) then
            rewind iscr
            if (j /= 1) then
              write (iscr, '(F12.6)') geo(j,i)*degree
            else
              write (iscr, '(F12.6)') geo(j,i)
            end if
            rewind iscr
            read (iscr, '(A)') line(j,i)
          else if (igeo(j,i) == (-2)) then
            nopt = nopt + 1
            if (simbol(nopt) /= '---------') then
              if (simbol(nopt)(1:1) == '-') then
                line(j,i)(4:) = simbol(nopt)(2:)
              else
                line(j,i)(4:) = simbol(nopt)
              end if
            else
              nbi = nb(i)
              nci = nc(i)
              if (j /= 3) nci = 0
              if (j == 1) nbi = 0
              call xxx (type(j), i, na(i), nbi, nci, line(j,i)(4:))
            end if
            optdat(nopt) = line(j,i)
          else if (igeo(j,i) < 0) then
            line(3,i) = line(3,(-igeo(j,i)))
            line(3,i)(3:3) = '-'
          else
            line(j,i) = line(j,igeo(j,i))
          end if
        end do
        if (maxtxt /= 0 .and. maxtxt /= 27) then
          blank = elemnt(labels(i))//txtatm(i)//'  '
          nchars = max(4,maxtxt + 2)
        else
          blank = elemnt(labels(i))
          nchars = 3
        end if
        if (labels(i) == 99) blank(1:1) = ' '
        select case (i)
        case (1)
          write (iprt, '(1X,A,I4,A,I4,A,I4,A,I4)') blank(:nchars)
        case (2)
          write (iprt, '(1X,A,I4,A,I4,A,I4,A,I4)') blank(:nchars), na(i), line(1,i)
        case (3)
          write (iprt, '(1X,A,I4,A,I4,A,I4,A,I4)') &
          blank(:nchars), na(i), line(1,i), nb(i), line(2,i)
        case default
          write (iprt, '(1X,A,I4,A,I4,A,I4,A,I4)') &
          blank(:nchars), na(i), line(1,i), nb(i), line(2,i), nc(i), line(3,i), 0
        end select
      end do
      write (iprt, *)
      do l = 1, 3
        do i = 1, nopt
          if (loc(2,i) /= l) cycle
          if (loc(2,i) /= 1) then
            write (iprt, '(A,F12.6)') optdat(i), geo(loc(2,i),loc(1,i))*degree
          else
            write (iprt, '(A,F12.6)') optdat(i), geo(loc(2,i),loc(1,i))
          end if
        end do
      end do
      return
      end subroutine geoutg
