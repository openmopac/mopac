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

subroutine names (ioptl, lused, n1, ires, nfrag, io, uni_res, mres)
    use molkst_C, only: natoms, numat, keywrd, id
    use MOZYME_C, only : nres, allres, maxres, nbackb, nxeno, mxeno, k, iatom, jatom, &
       & loop, bbone, angles, txeno, afn, allr, lstart_res
    use common_arrays_C, only : nat, labels, txtatm, txtatm1, ibonds, coord, breaks, nbonds
    use funcon_C, only : pi
    use chanel_C, only: iw
    implicit none
    integer, intent (in) :: n1, nfrag, io
    integer, intent (inout) :: ires, uni_res, mres
    logical, dimension (numat), intent (inout) :: ioptl
    integer, dimension (natoms), intent (inout) :: lused
!
    logical :: lreseq, first_res, iatoms(numat)
    integer :: i, iold, irold, j, l, ires_start, ires_loop, &
      nbreaks
    intrinsic Index, Max
!
    lreseq = (Index (keywrd, " RESEQ") /= 0)
    do i = 1, 4
      nbackb(i) = 0
      nxeno(i, 1) = 0
    end do
    mxeno = 1
!
!   WORK OUT THE BACKBONE
!
    iatom = n1
    iold = -999
    if (nfrag == 1) txtatm(:numat) = " "
    irold = ires + 1
    bbone(1,ires + 1) = n1
    ires = ires + 1
    ires_start = ires
    first_res = .true.
    iatoms = .false.
    jatom = 0
    do ires_loop = ires_start, maxres
      if (iatom == 0) goto 999
      if (iatoms(iatom)) goto 1000
      iatoms(iatom) = .true.
!
!  NXTMER identifies the backbone atoms to the next residue.
!  The nitrogen of the next residue is in NBACKB(4)
!
      ioptl(iatom) = .true.
      call nxtmer (iatom, nbackb)
      jatom = nbackb(4)
!
!  The residue is attached to IATOM and does not involve JATOM
!
!  NOW TO WORK OUT THE ATOMS IN THE RESIDUE
!
      uni_res = uni_res + 1
      i = io
      if (.not. first_res) then
        do i = 1, numat
          if (nat(i) /= 8) exit
        end do
      end if
      call atomrs (lused, ioptl, ires, n1, i, uni_res, first_res)
      lstart_res(ires) = .true.
      first_res = .false.
      bbone(1,ires + 1) = jatom!  N of residue
      bbone(2,ires) = nbackb(1)!  C(alpha) of residue
      bbone(3,ires) = nbackb(2)!  C1 of residue (Carboxylic)
      if (loop /= 1 .and. k > 0 .and. k /= 23 .and. i /= 1) then
        i = Max (Index (txeno(loop), " "), 2) - 1
        if (i /= 1 .and. index(keywrd, " ADD-H") == 0) &
          write (iw,"(a,i5, 4a,i5,a)") "         Residue:", ires, &
            " ("//allres(ires)(:3)//") contains a xeno group."
        if (allres(ires) (4:4) == "+") then
          write (iw, "(A,I5,A)") "         Residue:", ires, &
               & " plus xeno group has a positive charge"
        else if (allres(ires) (4:4) == "-") then
          write (iw, "(A,I5,A)") "         Residue:", ires, &
               & " plus xeno group has a negative charge"
        end if
        allres (ires) = allres (ires) (1:3) // "*"
      end if
      if (iatom == jatom) go to 1000
      if (iold == jatom) then
        write (iw,*) " STRUCTURE UNRECOGNIZABLE"
        call mopend("Structure Unrecognizable")
        exit
      end if
      iold = iatom
      iatom = jatom
      if (iatom == 0 .and. ires_loop == ires_start) go to 1010
      ires = ires + 1
      if (ires == 0) then
        if (txtatm1(1) == " " .or. index(keywrd, " CHECKZERO") /= 0) then
          call l_control("CHECKZERO", len("CHECKZERO"), 1)
            if (index(keywrd, " NOZERO") /= 0) then
              ires = 1
            else if (index(keywrd, " ZERO") /= 0) then
              ires = 0
            else
              write(iw,'(/10x, a)') "A residue that can have either the number 1 or 0 has been detected"
              write(iw,'(10x, a)') "Add keyword ""NOZERO"" if the number 1 is to be used for this residue"
              write(iw,'(10x, a)') "Add keyword ""ZERO"" if the number 0 is to be used for this residue"
              call mopend("DURING ASSIGNMENT OF RESIDUES AN AMBIGUITY IN RESIDUE NUMBER 0 WAS FOUND.")
              return
            end if
          end if
        end if
    end do
    return
 999 ires = ires - 1
1000 continue
!
!  Check to see if there is a stray N attached to the end of the chain
!
    if (jatom > 0) then
      if (nat(jatom) == 7 .and. txtatm(jatom) == " ") then
        l = 0
        do i = 1, nbonds(jatom)
          k = ibonds(i,jatom)
          if ( .not. ioptl(k) .and. nat(k) > 1) l = l + 1
        end do
        if (l == 0) then
          write(txtatm(jatom),'(a,i5,a,i5)')  "ATOM  ", jatom, "  N   UNK ", ires + 1
        end if
      end if
    end if
    if (index(keywrd, " ADD-H") /= 0) return
    if (ires > 0) then
      first_res = .true.
      do i = irold, ires
        mres = mres + 1
       if (.not. first_res) then
!
! Work out phi
!
          if (bbone(3,i - 1)*bbone(1,i) /= 0 .and. bbone(2,i)*bbone(3,i) /= 0) then
            call dihed(coord, bbone(3,i - 1), bbone(1,i), bbone(2,i), bbone(3,i), angles(1,mres))
          else
            angles(1,mres) = 0.d0
          end if
        end if
        if (i /= ires) then
!
! Work out psi
!
          if (bbone(1,i)*bbone(2,i) /= 0 .and. bbone(3,i)*bbone(1,i + 1) /= 0) then
            call dihed(coord, bbone(1,i), bbone(2,i), bbone(3,i), bbone(1,i + 1), angles(2,mres))
          else
            angles(2,mres) = 0.d0
          end if
!
! Work out omega
!
          if (bbone(2,i)*bbone(3,i) /= 0 .and. bbone(2,i + 1)*bbone(1,i + 1) /= 0) then
            call dihed(coord, bbone(2,i), bbone(3,i), bbone(1,i + 1), bbone(2,i + 1), angles(3,mres))
          else
            angles(3,mres) = 0.d0
          end if
        end if
        do j = 1,3
          if (angles(j,mres) >  pi) angles(j,mres) = angles(j,mres) - 2.d0*pi
          if (angles(j,mres) < -pi) angles(j,mres) =0.d0
          angles(j,mres) = angles(j,mres)*180.d0/pi
        end do
        first_res = .false.
      end do
      do i = 1, ires
        do j = 1, 20
          if(allres(i)(:3) == afn(j)(:3)) exit
        end do
        if ( j > 20) allr(i) = "X"
      end do
    end if
!
!  LABEL THE HYDROGEN ATOMS
!
    nbreaks = 1  !  In this subroutine, "nbreaks" is local
    do j = 1, numat
      if (nat(j) == 1) then
        if (j /= breaks(nbreaks)) then
          if (ibonds(1, j) > 0) lused(j) = lused(ibonds(1, j))
        else
          nbreaks = nbreaks + 1
        end if
      end if
    end do
    if (numat /= natoms + id.and. .not. lreseq) then
      j = natoms - numat
!
!  MODIFY TXTATM TO ALLOW FOR DUMMY ATOMS
!
      l = numat
      do i = natoms - id, 1, -1
        if (labels(i) /= 99) then
          txtatm(i) = txtatm(l)
          lused(i) = lused(l)
          l = l - 1
        else
          write (txtatm(i), "(I6)") j
          j = j - 1
        end if
      end do
    end if
    nres = ires
    goto 1020
1010 nres = 1
    if (allres (ires) (:3) /= "   ") then
      if (index(keywrd, " RESIDUE") /= 0) &
        write (iw,'(10x,a,i4,a,a)') " AMINO ACID (residue:",ires,") = ", allres (ires) (:3)
    end if
!
!   LABEL EVERY NON-BACKBONE ATOM AS A SIDE-CHAIN ATOM
!
    do i = 1, natoms
      if (lused(i) /=-1) then
        lused(i) = 1
      end if
    end do
    do j = 1, numat
      if (nat(j) == 1) then
        lused(j) = lused(ibonds(1, j))
      end if
    end do
    if (numat == natoms) return
    j = natoms - numat
!
!  MODIFY TXTATM TO ALLOW FOR DUMMY ATOMS
!
    l = numat
    do i = natoms - id, 1, -1
      if (labels(i) /= 99) then
        txtatm(i) = txtatm(l)
        lused(i) = lused(l)
        l = l - 1
      else
        write (txtatm(i), "(I6)") j
        j = j - 1
      end if
    end do
1020 continue
     return
end subroutine names
