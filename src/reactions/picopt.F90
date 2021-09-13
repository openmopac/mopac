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

subroutine picopt (loop)
    use molkst_C, only: natoms, numat, nvar, ndep, numcal
    use common_arrays_C, only: loc, nbonds, ibonds, labels
    use symmetry_C, only : locdep
    use MOZYME_C, only : jopt, numred
    implicit none
    integer, intent (in) :: loop
    integer, allocatable :: iopt_loc(:)
!
    integer :: i, j, imol = 0
   !***********************************************************************
   !                                                                      *
   !   PICOPT DETERMINES WHICH ATOMS ARE TO BE USED IN THE SCF.           *
   !                                                                      *
   !   THIS DETERMINATION IS BASED ON THE FOLLOWING:                      *
   !                                                                      *
   !     LOOP=-1:     USE ALL ATOMS                                       *
   !                  ELSE USE ATOMS MARKED FOR OPTIMIZATION IN LOC.      *
   !                  PLUS ANY ATOMS DEFINED BY SYMMETRY                  *
   !                                                                      *
   !    ON EXIT, JOPT   = ATOM NUMBERS OF ATOMS TO BE USED IN SCF         *
   !             NUMRED = NUMBER OF ATOMS IN JOPT                         *
   !                                                                      *
   !   iopt is scratch,                                                   *
   !   LOC, NBONDS, IBONDS, LABELS, and LOCDEP are reference arrays       *
   !                                                                      *
   !***********************************************************************
    allocate (iopt_loc(natoms))
    if (loop ==-1) then
      iopt_loc = 1
    else
      iopt_loc = 0
      !
      !   USE THE GEOMETRY SUPPLIED
      !
      iopt_loc = 0
      do i = 1, nvar
        iopt_loc(loc(1, i)) = 2
      end do
      do i = 1, ndep
        iopt_loc(locdep(i)) = 2
      end do
      !
      !   AT THIS POINT, iopt_loc REFERS TO ALL ATOMS, INCLUDING DUMMIES.
      !   IN THE FOLLOWING SECTION, WE NEED TO WORK WITH REAL ATOMS
      !   ONLY
      !
      j = 0
      do i = 1, natoms
        if (labels(i) /= 99) then
          j = j + 1
          iopt_loc(j) = iopt_loc(i)
        end if
      end do
      !
      !   IOPT IS NOW OVER REAL ATOMS ONLY
      !
      if (imol == numcal) then
        do i = 1, numat
          if (iopt_loc(i) == 2) then
            do j = 1, nbonds(i)
              if (iopt_loc(ibonds(j, i)) == 0) then
                iopt_loc(ibonds(j, i)) = 1
              end if
            end do
          end if
        end do
      end if
    end if
    imol = numcal
    numred = 0
    jopt = 0
    do i = 1, numat
      if (iopt_loc(i) /= 0) then
        numred = numred + 1
        jopt(numred) = i
      end if
    end do
    return
end subroutine picopt
