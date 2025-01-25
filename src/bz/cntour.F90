! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

subroutine cntour (a, mdim, imax, jmax, pz, value, i1, i2, i3, r11)
  use common_common, only : im, r1, jm
  implicit none
  integer, intent(in) :: i1, i2, i3, imax, jmax, mdim
  double precision, intent(in) :: pz, value
  double precision, intent(out) :: r11
  double precision, dimension(mdim, *), intent(inout) :: a
  ! 
  !.. Local Scalars .. 
  integer :: i, imax1, j, jmax1, k
  double precision :: aa, ab
  logical, dimension(40000) :: b
  logical, dimension(20000, 2) :: c
  ! 
  !***********************************************************************
  !
  ! GENERAL CONTOUR-DRAWING ROUTINE.
  !
  !  ON INPUT:
  !      A      = 2-D ARRAY OF POINTS DEFINING HEIGHT IN MAP-SPACE.
  !      MDIM   = SIZE OF FIRST DIMENSION OF A.
  !      IMAX   = LENGTH OF ONE SIDE OF ARRAY WITHIN A.
  !      JMAX   = LENGTH OF THE OTHER SIDE OF THE ARRAY WITHIN A.
  !      PZ     = VALUE OF Z DIRECTION IF 3-D PLOTS ARE BEING DRAWN.
  !      VALUE  = VALUE OF CONTOUR.
  !      I1,I2,I3 = 1,2,3 IF AN ORTHOSCOPIC PLOT IS WANTED, I.E. IF THE
  !               VIEW IS ALONG Z, AND THE PLOT IS IN X AND Y, THIS IS THE
  !               'NORMAL' PLOT. THIS OPTION HAS NEVER BEEN EXCERCISED
  !               FULLY.
  !  ON EXIT R11= LENGTH OF CONTOUR. THIS IS THE TOTAL LENGTH OF ALL LINES
  !               MAKING UP THE CONTOUR, OPEN AND CLOSED.
  !***********************************************************************
  im = imax
  r1 = 0.d0
  jm = jmax
  !
  !  SET UP STRING FOR LATER ANNOTATION OF CONTOURS.
  !
  !
  !  THE FOLLOWING CODE SUBTRACTS THE VALUE OF THE VALUE BEING SOUGHT
  !  FROM ALL ELEMENTS OF A, SO THAT THE ZERO VALUE WILL BE SOUGHT
  !  HEREAFTER.
  !
  do i = 1, imax
    do j = 1, jmax
      a(i, j) = a(i, j) - value
    end do
  end do
  r11 = r1
  jmax1 = jmax - 1
  imax1 = imax - 1
  do i = 1, imax
    do j = 1, jmax1
      aa = a(i, j)
      ab = a(i, j+1)
      !
      !  IF A CONTOUR  PASSES CLOSE TO AN ELEMENT OF A, THE CORRESPONDING
      !  ELEMENT OF LOGICAL ARRAY B IS SET.
      !
      b(imax*(j-1)+i) = ((aa>=0..and.ab<0.).or. (aa<0..and.ab>=0.))
    end do
  end do
  k = 1
  do j = 1, jmax, jmax1
    ab = a(1, j)
    do i = 1, imax1
      aa = ab
      ab = a(i+1, j)
      !
      !  THE ARRAY C HAS THE SAME FUNCTION AS B, BUT ON THE BOUNDARY OF THE
      !  AREA.
      !
      c(i, k) = ((aa>=0..and.ab<0.).or. (aa<0..and.ab>=0.))
    !
    !  THE ARRAY C HAS THE SAME FUNCTION AS B, BUT ON THE BOUNDARY OF THE
    !  AREA.
    !
    end do
    k = 2
  end do
  !
  !   DRAW UNCLOSED CONTOURS (IE STARTING ON A BOUNDARY).
  !
  do i = 1, imax1
    if (c(i, 1)) then
      call draw (i, 1, 1, a, b, c, pz, i1, i2, i3, mdim)
    end if
    if (c(i, 2)) then
      call draw (i+1, jmax, 3, a, b, c, pz, i1, i2, i3, mdim)
    end if
  end do
  !
  !  DRAW INTERIOR CONTOURS PASSING CLOSE TO BOUNDARY POINTS.
  !
  do j = 1, jmax1
    if (b((j-1)*imax+1)) then
      call draw (1, j+1, 4, a, b, c, pz, i1, i2, i3, mdim)
    end if
    if (b(j*imax)) then
      call draw (imax, j, 2, a, b, c, pz, i1, i2, i3, mdim)
    end if
  end do
  !
  !  DRAW OTHER INTERIOR CONTOURS.
  !
  do i = 2, imax1
    do j = 1, jmax1
      if (b((j-1)*imax+i)) then
        call draw (i, j, 2, a, b, c, pz, i1, i2, i3, mdim)
      end if
    end do
  end do
  !
  !  RESTORE "A" MATRIX TO ITS ORIGINAL VALUE.
  !
  do i = 1, imax
    do j = 1, jmax
      a(i, j) = a(i, j) + value
    end do
  end do
end subroutine cntour
