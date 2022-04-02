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

      subroutine gmetry(geo, coord)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : natoms, numcal, step, id, nvar, keywrd
      use common_arrays_C, only : labels, geoa, na, nb, nc, tvec, loc, &
      xparam, txtatm
      use elemts_C, only : elemnt
      use chanel_C, only : iw
      use funcon_C, only : pi
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision :: geo(3,natoms)
      double precision , intent(out) :: coord(3,natoms)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icalcn, j, i, mb, mc, ma, k, l, counter
      double precision :: sum, error, ccos, cosa, xb, yb, zb, rbc, xa, ya, za, xyb&
        , xpa, xpb, costh, sinth, ypa, sinph, cosph, zqa, yza, coskh, sinkh, &
        sina, sind, cosd, xd, yd, zd, ypd, zpd, xpd, zqd, xqd, yqd, xrd
      character, dimension(4) :: ndimen*16
      double precision, dimension(3, natoms) :: geovec
      save ndimen, icalcn, counter
!-----------------------------------------------
      data icalcn/ 0/, counter/0/
      data ndimen/ ' MOLECULE     ', ' POLYMER       ', ' LAYER STRUCTURE', &
        ' SOLID         '/
!***********************************************************************
!
!    GMETRY  COMPUTES COORDINATES FROM BOND-ANGLES AND LENGTHS.
! *** IT IS ADAPTED FROM THE PROGRAM WRITTEN BY M.J.S. DEWAR.
!
!    (A) IF STEP IS NON-ZERO (THIS IS THE CASE WHEN "SADDLE" IS USED)
!        THEN GEO IS FIRST MODIFIED BY SHIFTING THE INTERNAL COORDINATES
!        ALONG A RADIUS FROM GEOA TO PLACE GEO AT A DISTANCE STEP FROM GEOA.
!    (B) NORMAL CONVERSION FROM INTERNAL TO CARTESIAN COORDINATES IS DONE.
!
!  ON INPUT:
!         GEO    = ARRAY OF INTERNAL COORDINATES.
!         NATOMS = NUMBER OF ATOMS, INCLUDING DUMMIES.
!         NA     = ARRAY OF ATOM LABELS FOR BOND LENGTHS.
!
!  ON OUTPUT:
!         COORD  = ARRAY OF CARTESIAN COORDINATES
!
!***********************************************************************
!                                     OPTION (A)
      if (abs(step) > 1.D-4) then
        sum = 0.D0
        do j = 1, 3
          do i = 1, natoms
            geovec(j,i) = geo(j,i) - geoa(j,i)
            sum = sum + geovec(j,i)**2
          end do
        end do
        sum = sqrt(sum)
        error = (sum - step)/sum
      else
        error = 0.D0
        geovec = 0.d0
      end if
      geo = geo - error*geovec(:,:natoms)
!                                     OPTION (B)
      coord(:,1) = geo(:,1)
      if (natoms == 1) return
      if (na(2) == 1) then
        coord(1,2)   = coord(1,1) + geo(1,2)
        coord(2:3,2) = coord(2:3,1)
      else
        coord(:,2) = geo(:,2)
      end if
      if (natoms /= 2) then
        if (na(3) == 0) then
          coord(:,3) = geo(:,3)
        else
          ccos = cos(geo(2,3))
          if (na(3) == 1) then
            coord(1,3) = coord(1,1) + geo(1,3)*ccos
          else
            coord(1,3) = coord(1,2) - geo(1,3)*ccos
          end if
          coord(2,3) = coord(2,2) + geo(1,3)*sin(geo(2,3))
          coord(3,3) = coord(3,2)
        end if

        if(size(na) == 4) call exit(1)
        do i = 4, natoms
          if (na(i) == 0) then
            coord(:,i) = geo(:,i)  ! Coordinate is already Cartesian
          end if
        end do
        do i = 4, natoms
          if (na(i) /= 0) then
            cosa = cos(geo(2,i))
            mb = nb(i)
            mc = na(i)
            xb = coord(1,mb) - coord(1,mc)
            yb = coord(2,mb) - coord(2,mc)
            zb = coord(3,mb) - coord(3,mc)
            rbc = xb*xb + yb*yb + zb*zb
            if (rbc < 1.D-16) then
  !
  !     TWO ATOMS ARE COINCIDENT.  A FATAL ERROR.
  !
              write (iw, '(/10x, A,I4,A,I4,A, 18x, a)') 'ATOMS', mb, ' AND', mc, ' ARE COINCIDENT', 'CARTESIAN COORDINATES'
              write(iw,'(58x, a, 2(11x,a))')"X", "Y", "Z"
              write(iw, '(10x, a, i5, a, a, 3f12.5)')"Atom", mb, ": ", elemnt(labels(mb))//"("//trim(txtatm(mb))//")",coord(:,mb)
              write(iw, '(10x, a, i5, a, a, 3f12.5,/)')"Atom", mc, ": ", elemnt(labels(mc))//"("//trim(txtatm(mc))//")",coord(:,mc)
              write (iw,'(a)')" Geometry at point of failure:"
              call geout(1)
              call mopend ('TWO ATOMS ARE COINCIDENT.  A FATAL ERROR.')
              return
            else
              rbc = 1.0D00/sqrt(rbc)
            end if
            ma = nc(i)
            xa = coord(1,ma) - coord(1,mc)
            ya = coord(2,ma) - coord(2,mc)
            za = coord(3,ma) - coord(3,mc)
  !
  !     ROTATE ABOUT THE Z-AXIS TO MAKE YB=0, AND XB POSITIVE.  IF XYB IS
  !     TOO SMALL, FIRST ROTATE THE Y-AXIS BY 90 DEGREES.
  !
            xyb = sqrt(xb*xb + yb*yb)
            k = -1
            if (xyb <= 0.009d0) then
              xpa = za
              za = -xa
              xa = xpa
              xpb = zb
              zb = -xb
              xb = xpb
              xyb = sqrt(xb*xb + yb*yb)
              if (xyb < 0.009d0) then
                write (iw, '(/10x, A,I4,A,I4,A, 18x, a)') 'ATOMS', ma, ' AND', mc, ' ARE COINCIDENT', 'CARTESIAN COORDINATES'
                write(iw,'(58x, a, 2(11x,a))')"X", "Y", "Z"
                write(iw, '(10x, a, i5, a, a, 3f12.5)')"Atom", ma, ": ", elemnt(labels(ma))//"("//trim(txtatm(ma))//")",coord(:,ma)
                write(iw, '(10x, a, i5, a, a, 3f12.5,/)')"Atom", mc, ": ", &
                  elemnt(labels(mc))//"("//trim(txtatm(mc))//")",coord(:,mc)
                write (iw,'(a)')" Geometry at point of failure:"
                call geout(1)
                call mopend ('TWO ATOMS ARE COINCIDENT.  A FATAL ERROR.')
                return
              end if
              k = 1
            end if
  !
  !     ROTATE ABOUT THE Y-AXIS TO MAKE ZB VANISH
  !
            costh = xb/xyb
            sinth = yb/xyb
            xpa = xa*costh + ya*sinth
            ypa = ya*costh - xa*sinth
            sinph = zb*rbc
            cosph = sqrt(abs(1.D00 - sinph*sinph))
            zqa = za*cosph - xpa*sinph
  !
  !     ROTATE ABOUT THE X-AXIS TO MAKE ZA=0, AND YA POSITIVE.
  !
            yza = sqrt(ypa**2 + zqa**2)
            if (yza >= 1.D-4) then
              coskh = ypa/yza
              sinkh = zqa/yza
            else
  !
  !   ANGLE TOO SMALL TO BE IMPORTANT
  !
              coskh = 1.D0
              sinkh = 0.D0
            end if
  !
  !     COORDINATES :-   A=(???,YZA,0),   B=(RBC,0,0),  C=(0,0,0)
  !     NONE ARE NEGATIVE.
  !     THE COORDINATES OF I ARE EVALUATED IN THE NEW FRAME.
  !
            sina = sin(geo(2,i))
            sind = -sin(geo(3,i))
            cosd = cos(geo(3,i))
            xd = geo(1,i)*cosa
            yd = geo(1,i)*sina*cosd
            zd = geo(1,i)*sina*sind
  !
  !     TRANSFORM THE COORDINATES BACK TO THE ORIGINAL SYSTEM.
  !
            ypd = yd*coskh - zd*sinkh
            zpd = zd*coskh + yd*sinkh
            xpd = xd*cosph - zpd*sinph
            zqd = zpd*cosph + xd*sinph
            xqd = xpd*costh - ypd*sinth
            yqd = ypd*costh + xpd*sinth
            if (k >= 1) then
              xrd = -zqd
              zqd = xqd
              xqd = xrd
            end if
            coord(1,i) = xqd + coord(1,mc)
            coord(2,i) = yqd + coord(2,mc)
            coord(3,i) = zqd + coord(3,mc)
!
!  If the bond-angle cosine is near zero or 180 degrees, then re-define the torsion or dihedral connectivity
!  unless it is a solid.  In solids, the re-definition would need to be restricted to the same unit cell.
!
            if (Abs (cosa) > 0.9998d0 .and. i > 4 .and. id == 0) then
              call bangle (coord, na(i), nb(i), nc(i), sum)
                if (sum > 3.05d0 .or. sum < 0.09d0) then
                  j = nc(i)
                  call renum (coord, na, nb, nc, i, natoms)

                  if (nc(i) == 0) nc(i) = j
                  !
                  !  Correct GEO to show the new dihedral
                  !
                  if (nc(i) /= j) then
                    call dihed (coord, i, na(i), nb(i), nc(i), geo(3, i))
                    if (sum >= pi/2) sum = pi - sum
!
!  If angle is not in the domain 0 - 180 degrees, increase dihedral by 180 degrees,
!  because the dihedral is for the wrong angle.
!
                    if (Mod (Int(geo(2,i)/pi + 100.d0), 2) == 1) &
                 &  geo(3,i) = geo(3,i) + pi
!
! Update xparam, if necessary
!
                  do j = 1, nvar
                    if(loc(2,j) == 3 .and. loc(1,j) == i) then
                      xparam(j) = geo(3,i)
                      exit
                    end if
                  end do
                end if
              end if
            end if
          end if
        end do
!
! *** NOW REMOVE THE TRANSLATION VECTORS, IF ANY, FROM THE ARRAY COORD
!
      end if
      if (natoms == 0) then
        call mopend("SYSTEMS HAVE NO ATOMS IN COMMON!")
        return
      end if
      k = natoms
      do while(labels(k) == 107)
        k = k - 1
        if (k == 0) then
           call mopend ('SOLIDS WITHOUT ATOMS ARE NOT ALLOWED.')
           return
        end if
      end do
      k = k + 1
      if (icalcn /= numcal) id = 0
      if (k <= natoms) then
!
!   SYSTEM IS A SOLID, OF DIMENSION NATOMS+1-K
!
        l = 0
        do i = k, natoms
          l = l + 1
          mc = na(i)
          tvec(1,l) = coord(1,i)
          tvec(2,l) = coord(2,i)
          tvec(3,l) = coord(3,i)
        end do
        id = l
        if (icalcn /= numcal) then
          if (counter < 1) then
          counter = counter + 1
          go to 170
          else
          icalcn = numcal
          counter = 0
          end if
          if (index(keywrd, " SILENT") == 0) then
            write (iw, 140) trim(ndimen(id+1))
  140       format(/,19x,'THE SYSTEM IS A',a,/)
            if (id == 0) go to 170
            write (iw, 150)
            write (iw, 160) (i,(tvec(j,i),j=1,3),i=1,id)
  150       format(/,'                UNIT CELL TRANSLATION VECTORS',/,/,&
            '              X                Y              Z')
160         format('    T',i1,' =',f14.7,' ',f14.7,' ',f14.7)
          end if
        end if
      end if
  170 continue
      j = 0
      do i = 1, natoms
        if (labels(i)==99) cycle
        j = j + 1
        coord(:,j) = coord(:,i)
      end do
      return
      end subroutine gmetry
      subroutine renum (coord, na, nb, nc, ii, natoms)
      implicit none
      integer, intent (in) :: ii, natoms
      integer, dimension (natoms), intent (in) :: na, nb
      integer, dimension (natoms), intent (inout) :: nc
      double precision, dimension (3, natoms), intent (in) :: coord
      integer :: i, jj, nai, nbi
      double precision :: angle, rab, rmin, theta
      intrinsic Asin
 !***********************************************************************
 !
 !  Renumber the NC of atom II.  On input, the angle NA(II)-NB(II)-NC(II)
 !  is too near to 0 or 180 degrees.  Find a new atom for NC(II), so that
 !  the angle will be acceptable (as large as possible)
 !***********************************************************************
      nai = na(ii)
      nbi = nb(ii)
 !
 !   Theta = 45 degrees
 !
      theta = 0.7853d0
      jj = 0
      rmin = 1.d10
      do
        do i = 1, ii - 1
          if (i /= nai .and. i /= nbi) then
            call bangle (coord, nai, nbi, i, angle)
            if (angle > 1.5707963d0) then
              angle = 2.d0 * Asin (1.d0) - angle
            end if
            if (angle >= theta) then
             !
             !   Angle is OK.  Now find atom of lowest distance
             !
              rab = (coord(1, nbi)-coord(1, i)) ** 2 + (coord(2, &
             & nbi)-coord(2, i)) ** 2 + (coord(3, nbi)-coord(3, i)) ** 2
              if (rab < rmin) then
                jj = i
                rmin = rab
              end if
            end if
          end if
        end do
        if (jj /= 0) then
       !
       !  Best NC is JJ; best angle is THMIN
       !
          nc(ii) = jj
          exit
        end if
    !
    !   No atom inside the allowed angle - reduce the angle
    !
        theta = theta * 0.5d0
        if (theta < 0.0174533d0) theta = 0.d0
      end do
   end subroutine renum
