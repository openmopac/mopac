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

      subroutine perm(iperm, nels, nmos, nperms, limci)
      use meci_C, only : maxci
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nels
      integer , intent(in) :: nmos
      integer , intent(inout) :: nperms
      integer , intent(in) :: limci
      integer , intent(inout) :: iperm(nmos,maxci*4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: Upper, Lower, Prm, El, i
      integer, dimension (:), allocatable :: Occ
!-----------------------------------------------
!***********************************************************************
!
!  PERM PERMUTES NELS ENTITIES AMONG NMOS LOCATIONS. THE ENTITIES AND
!       LOCATIONS ARE EACH INDISTINGUISHABLE. THE PAULI EXCLUSION
!       PRINCIPLE IS FOLLOWED. THE NUMBER OF STATES PRODUCED IS GIVEN
!       BY NMOS!/(NELS!*(NMOS-NELS)!).
! ON INPUT: NELS  = NUMBER OF INDISTINGUISHABLE ENTITIES
!           NMOS  = NUMBER OF INDISTINGUISHABLE LOCATIONS
!
! ON OUTPUT IPERM = ARRAY OF PERMUTATIONS, A 0 INDICATES NO ENTITY,
!                   A 1 INDICATES AN ENTITY.
!           NPERM = NUMBER OF PERMUTATIONS.
!
!***********************************************************************
    Upper = nmos - nels
    Lower = 1
    Prm = 1
    El = nels
    allocate (Occ(-1:nels))
!
!  Occ initially contains the upper bounds that the various electrons
!  can go to.  Thus electron 1 can end up in M.O. "nmos",
!  and electron "nels-1" will end up in M.O. "nmos-nels".
!
    Occ(-1:nels-1) = (/(nmos+1-i, i=-1, nels-1) /)
!
!  For electron "El", the starting position is "1".
!
    call rperm (iperm, nperms, nels, Occ, Lower, Upper, El, limci)
    nperms = Prm-1
    deallocate (Occ)
  contains
    recursive subroutine rperm (iperm, nperms, nels, Occ, Lower, Upper, &
         & El, limci)
      use meci_C, only: nmos
    implicit none
      integer, intent (in) :: nperms, nels, limci
      integer, intent (in) :: Lower, Upper
      integer  :: El
      integer, dimension (-1:nels), intent(inout) :: Occ
      integer, dimension (nmos, nperms), intent(inout) :: iperm
      integer :: i, j, k
      if (El /= 0) then
        do j = Lower, Upper
          Occ(El) = j
          call rperm (iperm, nperms, nels, Occ, Occ(El)+1, Occ(El-2)-2, &
               & El-1, limci)
        end do
!
!  "do" loop is complete.  Now add in the final permutation.
!  (This is the reason for the second "-2" in "Occ(El-2)-2" in the call above.)
!
        Occ(El) = Upper + 1
!
!  Simplest test for exit - is the electron in an unacceptable M.O.,
!  or have all the permutations been calculated.
!
        if (Upper+1 > nmos .or. Prm > nperms) return
      end if
!
!  Generate the permutation.
!
      iperm(1:nmos, Prm) = 0
      do i = 1, nels
        iperm(Occ(i), Prm) = 1
      end do
      ! This statement only has an effect during the first call to rperm
      El = nels + 1
!
!   Eliminate microstates which are not allowed by "limci".
!
      if (limci /= 0 .and. Prm > 1) then
        k = 0
        do j = 1, nmos
          k = k + Abs (iperm(j, Prm)-iperm(j, 1))
        end do
        if (k > limci) then
          Prm = Prm - 1
        end if
      end if
      Prm = Prm + 1
      return
    end subroutine rperm
end subroutine perm
