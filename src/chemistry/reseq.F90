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

subroutine reseq (iopt, lused, n1, new, io)
!
!  on exit, lused: New sequence of atoms.  Atom n is original lused(n)
!
    use molkst_C, only: numat, natoms, keywrd
    use chanel_C, only: iw
    use elemts_C, only: elemnt, atom_names
    use common_arrays_C, only : nat, ibonds, nbonds, txtatm, coord
    implicit none
    integer, intent (in) :: n1
    integer, intent (inout) :: new
    integer, intent (out) :: io
    logical, dimension (numat), intent (inout) :: iopt
    integer, dimension (natoms), intent (inout) :: lused
!
    integer, parameter :: natomr = 10000, max_n_live = 200
    integer :: i, i2, i3, iatom, ico, ihcr, ires_loc, j, j2, jatom, k, l, &
   & ninres, nlive
    character :: num*1
    integer, dimension (4) :: nbackb
    integer, dimension (max_n_live) :: live
    integer, allocatable :: inres(:)
    logical :: iopt_ico
    logical, external :: peptide_n
!
    allocate (inres(natomr))
    nbackb = 0
    inres = 0    
    call lewis(.false.)
    i = Index (keywrd, " CVB")
    if (i /= 0) then
      call check_cvs(.true.)
    end if
!
!  Break all intra-chain bonds, so that the residues can easily be
!  identified.
!
    call lyse! Break all S-S bonds
!
!   WORK OUT THE BACKBONE
!
    iatom = n1
    do ires_loc = 1, 100000
      call nxtmer (iatom, nbackb)
      ihcr = nbackb(1)
      ico = nbackb(2)
      io = nbackb(3)
      jatom = nbackb(4)
!
!    STORE INDICES FOR  C (OF -HCR-) AND C (OF -CO-) in 'N-C(H,R)-C(=O)'
!   
      if (ico > 0) then
        iopt_ico = iopt(ico)
      else
        iopt_ico = .false.
      end if
      if (iopt(iatom) .or. iopt(ihcr) .or. iopt_ico) then
        l = 0
        j = 0
        io = 0
        if (iatom /= 0 .and. iatom == jatom) then
!
!   An amine or ammine group
!
          jatom = 0
          io = -1
        end if
      else if (io == 0) then
        lused(new+1) = iatom
        lused(new+2) = ihcr
        lused(new+3) = ico
        iopt(iatom) = .true.
        iopt(ihcr) = .true.
        iopt(ico) = .true.
        new = new + 3
        l = ihcr
        j = ico
      else if (iopt(io)) then
        l = 0
        j = 0
        io = 0
        if (iatom /= 0 .and. iatom == jatom) then
!
!  Not a residue, but there is something attached to the N at the END
!  of the protein (at the COO end), so turn off everything to do with
!  the backbone, and treat it as a xeno-group.
!
          jatom = 0
          ico = 0
          io = 0
          ihcr = 0
          if ( .not. iopt(iatom)) then
            iopt(iatom) = .true.
            new = new + 1
            lused(new) = iatom
          end if
        end if
      else
        lused(new+1) = iatom
        lused(new+2) = ihcr
        lused(new+3) = ico
        lused(new+4) = io
        iopt(iatom) = .true.
        iopt(ihcr) = .true.
        iopt(ico) = .true.
        iopt(io) = .true.
        new = new + 4
        l = ihcr
        j = ico
      end if
!
!  The residue is attached to IATOM and does not involve JATOM
!
!  NOW TO WORK OUT THE ATOMS IN THE RESIDUE
!
      ninres = 0
      if (io < 1) then
        nlive = 1
        live(1) = iatom
      else
        nlive = 0
        ninres = 0
        do i2 = 1, nbonds(iatom)
          if (nat(ibonds(i2, iatom)) == 1 .or. ibonds(i2, iatom) == ihcr) then
            nlive = nlive + 1
            live(nlive) = ibonds(i2, iatom)
          end if
        end do
      end if
      if (l /= 0) then
        do i2 = 1, nbonds(l)
          live(nlive+i2) = ibonds(i2, l)
        end do
        nlive = nlive + nbonds(l)
        do i2 = 1, nbonds(j)
          live(nlive+i2) = ibonds(i2, j)
        end do
        nlive = nlive + nbonds(j)
      end if
      do
        l = live(1)
        if (l == jatom .or. iopt(l) .or. (l == iatom .or. peptide_n(l)) .and. io > -1) then
          if (nlive == 0) exit
!
!  Pick the lowest numbered atom
!
          i3 = 0
          i2 = 1000000
          do i = 1, nlive
            if (i2 > live(i)) then
              i3 = i
              i2 = live(i)
            end if
          end do
          live(1) = live(i3)
          nlive = nlive - 1
          do i = i3, nlive
            live(i) = live(i + 1)
          end do          
        else
          iopt(l) = .true.
!
! Complicated way of saying "If there is a nitrogen atom on its own, add it to the set"
!
          if (io < 0 .and. (l == iatom .or. l == ihcr)) then
            new = new + 1
            lused(new) = l
          end if
          ninres = ninres + 1
          if (ninres > natomr) then
            write (iw, "(A,I4,A,I4)") " There are more than", natomr, &
           & " atoms in residue", ires_loc
            write (iw,*) " Atoms in residue", ires_loc
            write (iw, "(10(1X,A2,I5))") (elemnt(nat(inres(l))), inres(l), &
           & l=1, natomr)
            go to 1000
          else
            inres(ninres) = l
            if (nbonds(l) /= 0) then
!
!  THERE IS AT LEAST ONE ATOM ATTACHED TO THE 'LIVE' ATOM
!
              loop: do i2 = 1, nbonds(l)
                j = ibonds(i2, l)
                do i3 = 1, nlive
                  if (live(i3) == j) cycle loop
                end do
                nlive = nlive + 1
                if (nlive > max_n_live) then
                  go to 1000
                else
                  live(nlive) = j
                end if
              end do loop
              live(1) = ibonds(1, l)
            else
              if (nlive == 0) exit
              live(1) = live(nlive)
              nlive = nlive - 1
            end if
          end if
        end if
      end do
!
!   Check that atoms are not counted twice
!
      do j = 1, ninres
        do k = j + 1, ninres
          if (inres(j) == inres(k)) then
            inres(j) = 0
          end if
        end do
      end do
!
!  Put hydrogen atoms at the end of the list
!
      k = ninres
      do j = 1, ninres
        l = inres(j)
        if (l /= 0) then
          if (nat(l) == 1) then
            k = k + 1
            if (k > natomr) then
              goto 1000
            end if
            inres(k) = l
            inres(j) = 0
          end if
        end if
      end do
      ninres = k
      if (io == 0 .and. iatom /= 0 .and. .not. iopt(iatom)) then
        new = new + 1
        lused(new) = iatom
        iopt(iatom) = .true.        
      end if
      do i2 = 1, ninres
        j2 = inres(i2)
        if (j2 /= 0 .and. j2 /= iatom .and. j2 /= ihcr .and. j2 /= ico) then
          new = new + 1
          lused(new) = j2
        end if
      end do
      if (ico > 0) then
!
!  Check to see if a OXT atom exists.  If it does, then move it to the end of the heavy atoms
!
        j = 0
        do i = 1, nbonds(ico)
          if (nat(ibonds(i,ico)) == 8) then
            j = j + 1
            if (ibonds(i,ico) /= io) k = ibonds(i,ico)
          end if
        end do
        if (j == 2) then
  !
  !   Found a terminal -COO structure.  Now move the OXT to the end of the heavy atoms
  !
          do i = new, 1, -1
            if (lused(i) == k) exit
          end do
          do i = i + 1, numat
            if (lused(i) == -1000) exit
            if (nat(lused(i)) == 1) exit
            lused(i - 1) = lused(i)
          end do
          lused(i - 1) = k        
        end if
      end if
      j = new
      if (iatom == jatom .or. jatom == 0) go to 1100
      iatom = jatom
    end do
    return
1000 i = 0
    if (ihcr > 0 .and. ihcr <= numat) then
      i = ihcr
    else if (iatom > 0 .and. iatom <= numat) then
      i = iatom
    else if (ico > 0 .and. ico <= numat) then
      i = ico
    end if    
    if (index(keywrd, " ADD-H") > 0) then
      write(iw,'(/10x,a)') "Attempt to resequence the atoms after adding hydrogen atoms (ADD-H) failed"
    else if (index(keywrd, " SITE=") > 0) then
      write(iw,'(/10x,a)') "Attempt to resequence the atoms after adding or deleting hydrogen atoms (SITE) failed"
    else
      write(iw,'(/10x,a)') "Attempt to resequence the atoms by using keyword RESEQ failed"      
    end if      
    if (i > 0) then
      do j = 1, 10
        if (atom_names(nat(i))(j:j) /= " ") exit
      end do
      num = char(ichar("1") + int(log10(i*1.0001)))
      write (iw, "(/10X,A,I"//num//",a,/10x,a, /10x,a, 3F8.3)") &
      "The fault is likely to be near to "//atom_names(nat(i))(j:)//" atom, atom number: ", i, ",", &
      "PDB line label: """//txtatm(i)//""",", &
      "at coordinates:", coord(:,i)
    else
      write(iw,'(/10x,a)')" RESEQ procedure is completely confused!"
    end if
    if (index(keywrd, " ADD-H") + index(keywrd, " SITE=") > 0)  &
      write(iw,'(/10x,a)')"(To avoid this error, add keyword ""NORESEQ"")"
    call mopend("ATTEMPT TO RESEQUENCE ATOMS FAILED")
    call web_message(iw,"residue_too_big.html")
1100 continue
     if (nbackb(4) > 0) then
!
!  SPECIAL CASE - IS THERE A HYDROGEN ATTACHED TO THE FINAL NITROGEN?
       jatom = nbackb(4)
       do i = 1, nbonds(jatom)
         j = ibonds(i, jatom)
         if (.not. iopt(j) .and. nat(j) == 1) then
           iopt(j) = .true.
           new = new + 1
           lused(new) = j
         end if
       end do
     end if       
     if (io < 1) return
!
!  SPECIAL CASE - IS THERE A HYDROGEN ATTACHED TO THE FINAL OXYGEN?
!
    do i = 1, nbonds(io)
      j = ibonds(i, io)
      if ( .not. iopt(j) .and. nat(j) == 1) then
        iopt(j) = .true.
        new = new + 1
        lused(new) = j
      end if
    end do
    return
end subroutine reseq
