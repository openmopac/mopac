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

subroutine greek (n1)
    use MOZYME_C, only: nbackb
    use common_arrays_C, only : txtatm, nat, nbonds, ibonds
    implicit none
    integer, intent (in) :: n1
!
    integer :: i, ico, ihcr, ii, j, jatom1, jjj, k1, l, ll, m, ncurr, nnext
    integer :: nold4 = 0, nprev = 0
    character, dimension (22) :: types
    integer, dimension (150) :: icurr, inext, iprev
    logical :: Trp, Ile
    data types / "G", "D", "E", "Z", "H", "T", "I", "K", "L", "M", &
   & "N", "X", "O", "P", "R", "S", "T", "U", "F", "C", "Y", "W" /
   !
   !             NBACKB(1) =  Carbon CHR  of N-CHR-CO-N'
   !             NBACKB(2) =  Carbon CO   of N-CHR-CO-N'
   !             NBACKB(4) =  Nitrogen N' of N-CHR-CO-N'
    ihcr = nbackb(1)
    ico  = nbackb(2)
    jatom1 = nold4
    nold4 = nbackb(4)
   !
   !  Backbone atom 'IHCR' is ALPHA.
   !
   !   Special case:  avoid last residue being re-labeled.
   !
    do i = 1, nbonds(ihcr)
      if (nat(ibonds(i, ihcr)) == 8) go to 1000
    end do
    txtatm (ihcr) (15:15) = "A"
   !
   !  Find the BETA carbon, if present
   !
1000 do i = 1, nbonds(ihcr)
      j = ibonds(i, ihcr)
      if (nat(j) == 6) then
        if (j /= ico) go to 1010
      end if
    end do
    return
1010 txtatm (j) (15:15) = "B"
   !
   !  Now slowly go through the rest of the residue, assigning
   !  letters and numbers.
   !
    Trp = txtatm(j)(18:20) == "TRP"
    Ile = txtatm(j)(18:20) == "ILE"
    nnext = 1
    icurr(1) = j
    ncurr = 1
    iprev(1) = ihcr
    do jjj = 1, 22
      !
      !   NNEXT  =           Number of atoms of next type
      !   INEXT =     Atom numbers of atoms of next type
      !
      !   NCURR  =        Number of atoms of current type
      !   ICURR =  Atom numbers of atoms of current type
      !
      !   NPREV  =       Number of atoms of previous type
      !   IPREV = Atom numbers of atoms of previous type
      !
      nprev = ncurr
      ncurr = nnext
      nnext = 0
      do ii = 1, ncurr
        j = icurr(ii)
        if (j == n1) cycle  !  Unconditionally, do not label the N of the peptide.
        k1 = nbonds(j)
        loop1: do l = 1, k1
          m = ibonds(l, j)
            !
            !  Check that atom 'M' is not (a) a previous atom,
            !  (b) a current atom, or (c) the backbone nitrogen of a proline.
            !
          do ll = 1, nprev
            if (m == iprev(ll)) cycle loop1
          end do
          do ll = 1, ncurr
            if (m == icurr(ll)) cycle loop1
          end do
          if (m /= jatom1) then
            nnext = nnext + 1
            inext(nnext) = m
          end if
          if (nat(m) == 7 .and. txtatm(ihcr)(18:20) == "PRO") then
            nnext = max(0, nnext - 1)
          end if
        end do loop1
      end do
      if (nnext == 0) exit
      if (Trp) then
        if (jjj == 2) then
!
!  Make sure that the carbon attached to a nitrogen is first
!
          j = inext(1)
          do l = 1, nbonds(j)
            if (nat(ibonds(l,j)) == 7) exit
          end do
          if (l > nbonds(j)) then
            l = inext(1)
            inext(1) = inext(2)
            inext(2) = j
          end if
        end if
      else if (Ile) then
        if (jjj == 1) then
!
!  Make sure that CG1 is the carbon of the ethyl group in isoleucine
!
         j = inext(1)
          m = 0
          do l = 1, nbonds(j)
            if (nat(ibonds(l,j)) == 6) m = m + 1
          end do
          if (m == 1) then
            l = inext(1)
            inext(1) = inext(2)
            inext(2) = j
          end if

        end if
      end if
      call txtype (nnext, inext, types(jjj))
      !
      !  Copy current atoms into previous atoms
      !
       iprev(:ncurr) = icurr(:ncurr)
      !
      !  Copy next atoms into current atoms
      !
      icurr(:nnext) = inext(:nnext)
    end do
    return  ! go to 1020
end subroutine greek
