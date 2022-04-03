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

subroutine cell(coord, tvec)
!*********************************************************************
!
!   GEOUT PRINTS THE CURRENT GEOMETRY.  IT CAN BE CALLED ANY TIME,
!         FROM ANY POINT IN THE PROGRAM AND DOES NOT AFFECT ANYTHING.
!         IF MODE1 .EQ.1 THEN GEOMETRY IS PRINTED IN USUAL .OUT FORMAT
!                  .GT.1 THEN GEOMETRY IS PRINTED IN .DAT FORMAT to
!                  CHANNEL MODE1
!
!*********************************************************************
    use common_systm, only : natoms, iw, maxtxt, ndep, numat, line
    use common_symult, only : depmul
    use common_keywrd, only : keywrd
    use common_geosym, only : locpar, idepfn, locdep
    use common_elemts, only : elemnt
    implicit none
    integer, parameter :: n_spacer = 4
    double precision :: coord(3,natoms + 3), tvec(3,3)
    integer :: i, j, k, l, m, ic, iend
    double precision :: cell_paras(6), a, b, c, alpha, beta, gamma
    character :: spacer(n_spacer)*1
    double precision, external :: reada
    equivalence (cell_paras(1), a), (cell_paras(2), b), (cell_paras(3), c), (cell_paras(4), alpha), &
      (cell_paras(5), beta), (cell_paras(6), gamma)
    data spacer /" ", "=", ";", ","/
      i = index(keywrd, " CELL") + 4
      keywrd(i + 1:i + 1) = " "
      iend = index(keywrd(i:), ")")
      if (iend == 0) then
        write(iw,'(/10x, a)')" Closing parethesis for keyword ""CELL"" not found"
        stop
      else
        iend = iend + i
      end if
      do j = 1, 6
!
!  Find a spacer
!
        m = iend
        do k = 1, n_spacer
          l = index(keywrd(i:iend), spacer(k))
          if (l /= 0) m = min(m, l)
        end do
        i = i + m
        if (i > iend) exit
        cell_paras(j) = reada(keywrd, i)
!
!  Find the end of the spacer
!
        do l = 1, 60
          i = i + 1
          if (i > iend) exit
          do k = 1, n_spacer
            if (keywrd(i:i) == spacer(k)) exit
          end do
          if (k > n_spacer) exit
        end do
      end do
      if (j < 7) then
        write(iw,'(/10x, a, i1, a)')"Only ", j, " cell parameters found, 6 were expected"
        i = index(keywrd, " CELL") + 6
        write(iw,'(/10x, a)')"Keyword: ""CELL("//keywrd(i:iend-1)//""""
        stop
      end if
      alpha = alpha*0.01745329252d0
      beta  = beta*0.01745329252d0
      gamma = gamma*0.01745329252d0
      tvec = 0.d0
      tvec(1,1) = a
      tvec(1,2) = b*cos(gamma)
      tvec(2,2) = b*sin(gamma)
      tvec(1,3) = c*cos(beta)
      tvec(2,3) = c*(cos(alpha) - cos(beta)*cos(gamma))/sin(gamma)
      tvec(3,3) = sqrt(c**2 - tvec(1,3)**2 - tvec(2,3)**2)
      return
    end subroutine cell
