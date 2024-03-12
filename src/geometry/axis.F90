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

      subroutine axis(a, b, c, evec)
      use molkst_C, only : numcal, keywrd, mol_weight, numat
      use common_arrays_C, only : atmass, coord
      USE chanel_C, only : iw
      use to_screen_C, only : rot, xyzmom
      USE funcon_C, only : fpc_6, fpc_8, fpc_10, pi
      implicit none
      double precision , intent(out) :: a
      double precision , intent(out) :: b
      double precision , intent(out) :: c
      double precision  :: evec(3,3)
      integer :: icalcn, i, j
      double precision, dimension(6) :: t
      double precision, dimension(numat) :: x, y, z
      double precision, dimension(3) ::  eig
      double precision :: const1, const2, sumwx, sumwy, sumwz, sum
      logical :: first
      save t, eig, first, icalcn
!***********************************************************************
!
!  AXIS CALCULATES THE THREE MOMENTS OF INERTIA AND THE MOLECULAR
!       WEIGHT.  THE MOMENTS OF INERTIA ARE RETURNED IN A, B, AND C.
!       THE MOLECULAR WEIGHT IN mol_weight.
!       THE UNITS OF INERTIA ARE 10**(-40)GRAM-CM**2,
!       AND MOL.WEIGHT IN ATOMIC-MASS-UNITS. (AMU'S)
!***********************************************************************
      data t/ 6*0.D0/
      data icalcn/ 0/
      if (icalcn /= numcal) then
        icalcn = numcal
        first = .TRUE.
      end if
!***********************************************************************
!     CONST1 =  10**40/(N*A*A)
!               N = AVERGADRO'S NUMBER
!               A = CM IN AN ANGSTROM
!               10**40 IS TO ALLOW UNITS TO BE 10**(-40)GRAM-CM**2
!
!***********************************************************************
      const1 = 10.d0**24/fpc_10
!***********************************************************************
!
!     CONST2 = CONVERSION FACTOR FROM ANGSTROM-AMU TO CM**(-1)
!
!            = (PLANCK'S CONSTANT*N*10**16)/(8*PI*PI*C)
!            = H[ERG-SEC]*N*10**16/
!              (8*PI**2*C[CM/SEC])
!
!***********************************************************************
      const2 = fpc_6*fpc_10*1.D16/(8.D0*pi**2*fpc_8)
!    FIRST WE CENTRE THE MOLECULE ABOUT THE CENTRE OF GRAVITY,
!    THIS DEPENDS ON THE ISOTOPIC MASSES, AND THE CARTESIAN GEOMETRY.
!
      sumwx = 0.D0
      sumwy = 0.D0
      sumwz = 0.D0
!
      if (mol_weight > 0) then
        do i = 1, numat
          sumwx = sumwx + atmass(i)*coord(1,i)
          sumwy = sumwy + atmass(i)*coord(2,i)
          sumwz = sumwz + atmass(i)*coord(3,i)
        end do
      else
        mol_weight = mol_weight + dble(numat)
        do i = 1, numat
          sumwx = sumwx + coord(1,i)
          sumwy = sumwy + coord(2,i)
          sumwz = sumwz + coord(3,i)
        end do
      end if
!
      if (mol_weight>0 .and. first) &
      write (iw, '(/10X,''MOLECULAR WEIGHT ='',F8.2,/)') min(99999.99D0,mol_weight)
      sumwx = sumwx/mol_weight
      sumwy = sumwy/mol_weight
      sumwz = sumwz/mol_weight
      x(:numat) = coord(1,:numat) - sumwx
      y(:numat) = coord(2,:numat) - sumwy
      z(:numat) = coord(3,:numat) - sumwz
!***********************************************************************
!
!    MATRIX FOR MOMENTS OF INERTIA IS OF FORM
!
!           |   Y**2+Z**2                         |
!           |    -Y*X       Z**2+X**2             | -I =0
!           |    -Z*X        -Z*Y       X**2+Y**2 |
!
!***********************************************************************
      do i = 1, 6
        t(i) = dble(i)*1.0D-10
      end do
!
      if (mol_weight > 0) then
        do i = 1, numat
          t(1) = t(1) + atmass(i)*(y(i)**2+z(i)**2)
          t(2) = t(2) - atmass(i)*x(i)*y(i)
          t(3) = t(3) + atmass(i)*(z(i)**2+x(i)**2)
          t(4) = t(4) - atmass(i)*z(i)*x(i)
          t(5) = t(5) - atmass(i)*y(i)*z(i)
          t(6) = t(6) + atmass(i)*(x(i)**2+y(i)**2)
        end do
      else
        do i = 1, numat
          t(1) = t(1) + (y(i)**2+z(i)**2)
          t(2) = t(2) - x(i)*y(i)
          t(3) = t(3) + (z(i)**2+x(i)**2)
          t(4) = t(4) - z(i)*x(i)
          t(5) = t(5) - y(i)*z(i)
          t(6) = t(6) + (x(i)**2+y(i)**2)
        end do
      end if
!
      call rsp (t, 3, eig, evec)
      if (mol_weight>0 .and. first .and. index(keywrd,'RC=')==0) then
        write (iw,'(2/9X,'' ROTATIONAL CONSTANTS IN CM(-1)'',/)')
        where (eig < 3.D-4)
          eig = 0.D0
          rot = 0.D0
        elsewhere
          rot = const2/eig
        end where
        xyzmom = eig*const1
        write (iw, &
      '(10X,''A ='',F14.8,''   B ='',F14.8,''   C ='',F14.8,/)') (rot(i),i=1,3)
        if (index(keywrd,'RC=') == 0) write (iw, &
      '(2/10X,'' PRINCIPAL MOMENTS OF INERTIA IN UNITS OF 10**(-40)*GRAM-CM**2'',/)')
        write (iw, &
      '(10X,''A ='',F14.4,''   B ='',F14.4,''   C ='',F14.4,/)') (xyzmom(i),i=1,3)
        c = rot(1)
        b = rot(2)
        a = rot(3)
      end if
!
!     MAKE DIAGONAL TERMS OBLIGATE POSITIVE
!
      do i = 1, 3
        if (evec(i,i) >= 0.D0) cycle
        evec(:,i) = -evec(:,i)
      end do
!
!   NOW TO ORIENT THE MOLECULE SO THE CHIRALITY IS PRESERVED
!   CHIRALITY CAN ONLY BE LOST IF ONE OR MORE EVEC(I,I) ARE ZERO
!
      sum = evec(1,1)*(evec(2,2)*evec(3,3)-evec(3,2)*evec(2,3)) + &
            evec(1,2)*(evec(2,3)*evec(3,1)-evec(2,1)*evec(3,3)) + &
            evec(1,3)*(evec(2,1)*evec(3,2)-evec(2,2)*evec(3,1))
      if (sum < 0) then
        sum = 1.D0
        do j = 1, 3
          if (evec(j,j) >= sum) cycle
          sum = evec(j,j)
          i = j
        end do
        evec(:,i) = -evec(:,i)
      end if
      if (index(keywrd,' NOREOR') + index(keywrd,' FORCETS') == 0) then
        coord(1,:numat) = x(:numat)
        coord(2,:numat) = y(:numat)
        coord(3,:numat) = z(:numat)
      end if
      if (mol_weight > 0) first = .FALSE.
      return
      end subroutine axis
