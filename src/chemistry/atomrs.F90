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

subroutine atomrs (lused, ioptl, ires, n1, io, uni_res, first_res)
    use molkst_C, only: natoms, numat, line, keywrd, ncomments, numat_old
    use chanel_C, only: iw, log, ilog
    use MOZYME_C, only: nbackb, jatom, iatom, k, loop, mxeno, txeno, &
       & nxeno,  allres, maxres,  ions, res_start, &
       tyres, size_mres, odd_h, atomname, afn, n_add
    use common_arrays_C, only : txtatm, nat, ibonds, nbonds, coord, all_comments, coorda, &
      txtatm1
    use elemts_C, only: elemnt, atom_names, cap_elemnt
    use funcon_C, only : pi
    implicit none
    integer, intent (in) :: io, n1, uni_res
    integer, intent (inout) :: ires
    logical, dimension (numat), intent (inout) :: ioptl
    integer, dimension (natoms), intent (inout) :: lused
    logical, intent (in) :: first_res
    integer, parameter :: natomr = 10000
!
    integer :: i2
    integer :: i, i3, ico, ihcr, j, j2, j3, j4, j5, j6, j7, OXT, live(500), &
      j8, j9, j10, jofco, l, ninres, nlive, nh, jj3, extra_atoms(2), once(5), n_once = 0
      double precision :: sum
    character (len=3) :: loc_tyres, tmp
    integer, dimension (natomr) :: inres
    integer, dimension (natomr) :: npack
    integer, dimension (5, size_mres) :: mres
    logical :: okay, found(14), UNK
    intrinsic Abs, Max
    logical, external :: peptide_n
    save :: once, n_once
    data inres / natomr * 0 /
! No. of        C  N  O  S  H    C  N  O  S  H     C  N  O  S  H
    data mres / 2, 1, 1, 0, 3,   3, 1, 1, 0, 5,    5, 1, 1, 0, 9, & !  GLY ALA VAL
                6, 1, 1, 0,11,   6, 1, 1, 0,11,    3, 1, 2, 0, 5, & !  LEU ILE SER
                4, 1, 2, 0, 7,   4, 1, 3, 0, 5,    4, 2, 2, 0, 6, & !  THR ASP ASN
                6, 2, 1, 0,13,   5, 1, 3, 0, 7,    5, 2, 2, 0, 8, & !  LYS GLU GLN
                6, 4, 1, 0,13,   6, 3, 1, 0, 7,    9, 1, 1, 0, 9, & !  ARG HIS PHE
                3, 1, 1, 1, 5,  11, 2, 1, 0,10,    9, 1, 2, 0, 9, & !  CYS TRP TYR
                5, 1, 1, 1, 9,   6, 2, 1, 0, 7,    5, 1, 2, 0, 9, & !  MET PRO PRO
                5, 1, 2, 0, 9,   0, 0, 0, 0, 0/                     !  PRO
    ninres = 0
   !             NBACKB(1) =  Carbon CHR  of N-CHR-CO-N'
   !             NBACKB(2) =  Carbon CO   of N-CHR-CO-N'
   !             NBACKB(3) =  Oxygen CO   of N-CHR-CO-N'
   !             NBACKB(4) =  Nitrogen N' of N-CHR-CO-N'
    ihcr = nbackb(1)
    ico = nbackb(2)
    jofco = nbackb(3)
    jatom = nbackb(4)
   !
   !   Starting with the nitrogen, work out the atoms in the residue.
   !
    nlive = nbonds(iatom)
    j3 = 0
    j4 = 0
    extra_atoms = 0
    do i2 = 1, nlive
      j2 = ibonds(i2, iatom)
      if (nat(j2) == 1 .or. j2 == ihcr) then
 !
 ! Select atoms that are definitely in the residue
 !
        j3 = j3 + 1
        live(j3) = j2
      else
 !
 ! Store list of atoms that are probably not in the residue. (Think of proline)
 !
        sum = 0.d0
        do i = 1, nlive
          if (nat(ibonds(i,iatom)) == 1) then
            do j = 1, nbonds(j2)
              if (nat(ibonds(j,j2)) == 8) then
!
!  Found a H-N-C-O sequence.  Check torsion angle
!
                call dihed (coord, ibonds(i,iatom), iatom, j2, ibonds(j,j2), sum)
              end if
            end do
          end if
        end do
        if (sum < pi*0.5d0 .and. sum > -pi*0.5d0) then
          j4 = min(2, j4 + 1)
          extra_atoms(j4) = j2
        end if
      end if
    end do
    nlive = j3
    outer_loop: do
      l = live(1)
      if (l == 0) go to 1100
      if (l == jatom .or. ioptl(l) .or. l == iatom .or. peptide_n(l)) then
         !
         !   Atom L is either already known to be in the residue, or is known
         !   to not be in the residue.
         !
        if (nlive == 0) go to 1100
        live(1) = live(nlive)
        nlive = nlive - 1
      else
        ioptl(l) = .true.
        ninres = ninres + 1
        if (ninres > natomr) then
          write (iw, "(A,I4,A,I4)") " There are more than", natomr, &
         & " atoms in residue", ires
          write (iw,*) " Atoms in residue", ires
          write (iw, "(10(1X,A2,I5))") (elemnt(nat(inres(l))), &
               & inres(l), l=1, natomr)
          if (Index (keywrd, " ADD-H") /= 0)  &
            write(iw,'(/10x,a)')"(To avoid this error, add keyword ""NORESEQ"" and re-run.)"
          call mopend ("Too many atoms in residue")
          go to 1000
        else
            !
            !  Assign atom L to the residue.
            !
          inres(ninres) = l
          if (nbonds(l) /= 0) then
               !
               !  THERE IS AT LEAST ONE ATOM ATTACHED TO THE 'LIVE' ATOM
               !
            loop1: do i2 = 2, nbonds(l)
              j = ibonds(i2, l)
              do i3 = 1, nlive
                if (live(i3) == j) cycle loop1
              end do
              nlive = nlive + 1
                  !
                  !   Panic if number of live atoms is greater than 500
                  !
              if (nlive > 500) exit outer_loop
              live(nlive) = j
            end do loop1
            live(1) = ibonds(1, l)
          else
            if (nlive == 0) go to 1100
            live(1) = live(nlive)
            nlive = nlive - 1
          end if
        end if
      end if
    end do outer_loop
    write (iw, "(10x,A,I4,A)") " Number of live atoms in residue", ires, &
                         & " greater than 500"
    do i = 1, numat
      if (index(txtatm(i), "UNK") > 0) exit
    end do
    if (i > 1 .and. i <= numat) &
      write(iw,'(/10x,a)')"The fault is probably after atom:  """//trim(txtatm(i - 1))//""""
    if (txtatm1(jatom) /= " ") then
      write(iw,'(/10x,a,a)')"The fault is probably before atom: """//trim(txtatm1(jatom))//""""
    else
      write(iw,'(/10x,a, i5)')"The fault is probably before atom: "//elemnt(nat(jatom)), jatom
    end if
    call mopend ("More than 500 atoms in residue")
1000 if (nat(io) == 8) then
      write (iw, "(10X,A,I5,A,3(/10X,A))") "Oxygen atom", io, &
     & " was incorrectly assumed to be in the backbone ", &
     & "of a residue.  To correct the fault, please move this atom", &
     & "to the end of the data set."
     end if
     if (Index (keywrd, " ADD-H") /= 0)  &
       write(iw,'(/10x,a)')"(To avoid this error, add keyword ""NORESEQ"" and re-run.)"
    return
   !
   !   Check that atoms are not counted twice
   !
1100 do j = 1, ninres
      do k = j + 1, ninres
        if (inres(j) == inres(k)) then
          inres(j) = 0
        end if
      end do
    end do
    npack(:) = 0
!
!  If the starting atom of the residue was a nitrogen, then count it
!
    if (nat(iatom) == 7) npack(7) = 1
    jj3 = 0
    do i2 = 1, ninres
      if (inres(i2) /= 0) then
        j2 = inres(i2)
!
!  Test for strange atoms
!
        do j4 = 1, nbonds(j2)
          j5 = nat(ibonds(j4,j2))
          select case (j5)
          case (1,6:8,16)
          case default
            jj3 = j5
          end select
        end do
        npack(nat(j2)) = npack(nat(j2)) + 1
      end if
      lused(inres(i2)) = ires
    end do
    lused(iatom) = -ires
    lused(ihcr) = -ires
    if (ico >0) lused(ico) = -ires
    if (jofco /= 0) then
      lused(jofco) = -ires
    end if
   !
   !   Determine the residue name (it must be one of the 21 amino acids,
   !   or an amino acid plus something attached.
   !
   OXT = 0
   if (ico > 0) then
!
!  Check for terminal -COOH  If present, delete one oxygen from formula count.
!
      k = 0
      i2 = 0
      j3 = 0
      j4 = 0
      do j = 1, nbonds(ico)
        j2 = ibonds(j,ico)
        if (nat(j2) == 8) then
          k = k + 1
          if (nbonds(j2) == 2) i2 = j2
          if (nbonds(j2) == 1) j3 = j3 + 1
          j4 = j2
        end if
      end do
      if (j3 == 2) i2 = j4
    end if
    if (k == 2) then
      OXT = i2
      npack(8) = npack(8) - 1
    else
      OXT = 0
    end if
    nh = npack(1)  !  Store number of hydrogen atoms in the residue
    k = 0
    loc_tyres = " "
    if (jj3 == 34) then
!
!  Convert Se into selenomethionine
!
      npack(16) = npack(16) + 1
    end if
    if (jj3 == 0 .or. jj3 == 34) then
      do loop = 1, mxeno
        npack(1) = npack(6) - nxeno(1, loop)
        npack(2) = npack(7) - nxeno(2, loop)
        npack(3) = npack(8) - nxeno(3, loop)
        npack(4) = npack(16) - nxeno(4, loop)
        if (iatom == jatom .and. ico > 0) then
!
!  Terminal group.  Check for -COOH
!
          j = 0
          do i = 1, nbonds(ico)
            if (nat(ibonds(i,ico)) == 8) j = j + 1
          end do
          if (j == 2) npack(3) = Max (0, npack(3)-1)
        end if
        loop2: do j = 1, size_mres - 1
          do k = 1, 4
            if (npack(k) /= mres(k, j)) cycle loop2
          end do
          if (j == 21 .or. j == 22) then
!
!  Check for hydroxyproline
!
            j2 = 0
            do j3 = 1, nbonds(ico)
              if (nat(ibonds(j3,ico)) == 8) j2 = j2 + 1
            end do
            if (j2 /= 2) cycle loop2
          end if
          if (j == 3) then
  !
  !  Check that the side chain involves isopropyl.  If so, it's valine
  !
  !
            j4 = 0
            do k = 1, nbonds(ihcr)
              l = 0
              j2 = ibonds(k,ihcr)
              if (nat(j2) == 6) then
                do j3 = 1, nbonds(j2)
                  if (nat(ibonds(j3,j2)) == 6) l = l + 1
                end do
              end if
              j4 = max(j4,l)
            end do
            if (j4 == 2) then
  !
  ! No, it's not valine.  Check for terminal proline
  !
              okay = .false.
              proline_check: do k = 1, nbonds(ihcr)
                if (nat(ibonds(k,ihcr)) == 7) then
  !
  ! Found the nitrogen
  !
                  j2 = ibonds(k,ihcr)
                  do j3 = 1, nbonds(j2)
                    if (nat(ibonds(j3,j2)) == 6 .and. ibonds(j3,j2) /= ihcr ) then
  !
  ! Found a C attached to N, and the C is not the C of -CHR-
  !
                      j4 = ibonds(j3,j2)
                      do j5 = 1,nbonds(j4)
                        if (nat(ibonds(j5,j4)) == 6) then
                          j6 = ibonds(j5,j4)
                          do j7 = 1, nbonds(j6)
                            if (nat(ibonds(j7,j6)) == 6) then
  !
  ! Found the C of N-C-C
  !
                              j8 = ibonds(j7,j6)
                              do j9 = 1, nbonds(j8)
                                do j10 = 1, nbonds(ihcr)
                                  if (ibonds(j9,j8) == ibonds(j10,ihcr)) then
                                    okay = .true.
                                    exit proline_check
                                  end if
                                end do
                              end do
                            end if
                          end do
                        end if
                      end do
                    end if
                  end do
                end if
              end do proline_check
              if (okay) then
                k = 20
              else
                cycle loop2
              end if
            end if
          end if
          k = j
          go to 1200
        end do loop2
        j = 0
        if (iatom == jatom) then
           !
           !   SPECIAL CASE:  THE -COOH END OF THE PROTEIN IS ONLY -C=O
           !
          if (npack(3) /= 0) then
            k = 1
            npack(3) = npack(3) + 1
          end if
        else if (jatom == 0) then
           !
           !   SPECIAL CASE:  SINGLE AMINO ACID, THEREFORE MAKE END -COO
           !
          k = 1
          npack(3) = npack(3) - 1
        end if
        if (k /= 0) then
           !
           !   CHECK ALL SPECIAL CASES
           !
          loop3: do k = 1, size_mres - 1
            do j = 1, 4
              if (npack(j) /= mres(j, k)) cycle loop3
            end do
            if (k == 21 .or. k == 22) then
!
!  Check for hydroxyproline
!
              j2 = 0
              do j3 = 1, nbonds(ico)
                if (nat(ibonds(j3,ico)) == 8) j2 = j2 + 1
              end do
              if (j2 /= 2) then
                loc_tyres = "HYP"
                cycle loop3
              end if
            end if
            if (k == 3) then
              j4 = 0
              do j = 1, nbonds(ihcr)
                l = 0
                j2 = ibonds(j,ihcr)
                do j3 = 1, nbonds(j2)
                  if (nat(ibonds(j3,j2)) == 6) l = l + 1
                end do
                j4 = max(j4,l)
              end do
              if (j4 /= 3) cycle loop3
            end if
            go to 1200
          end do loop3
        end if
      end do
    else
      do j = 1, 10
        if (atom_names(jj3)(j:j) /= " ") exit
      end do
      do j6 = 1, numat
        if (abs(coord(1,ihcr) - coorda(1,j6)) < 0.0005d0) then
          if (abs(coord(2,ihcr) - coorda(2,j6)) < 0.0005d0) then
            if (abs(coord(3,ihcr) - coorda(3,j6)) < 0.0005d0) exit
          end if
        end if
      end do
      do j7 = 1, n_once
        if (once(j7) == j6) exit
      end do
      if (j7 > n_once) then
        n_once = j7
        once(n_once) = j6
        write(iw,"(a,a,a)")"         Residue: '",txtatm1(j6)(18:), &
          "' contains "//atom_names(jj3)(j:)//" so identification not done"
      end if
      loop = 1
      goto 1200
    end if
    if (npack (6) == 0 .and. npack (7) == 1 .and. npack (8) == 0 .and. npack (16) == 0) then
      ioptl(jatom) = .true.
      ires = ires - 1
      return
    end if
1200 if (k == 0) then
       k = size_mres
     end if
     if (k == 4) then
      !
      !   RESIDUE IS EITHER LEU OR ILE.  TO DISTINGUISH BETWEEN THEM,
      !   CHECK THE NUMBER OF CARBON ATOMS ATTACHED TO THE CARBON
      !   OF THE SIDE-CHAIN.  IF IT IS 3 THEN ILE, ELSE LEU.
      !
      do i2 = 1, nbonds(iatom)
        j = ibonds(i2, iatom)
         !
         !  IS ATOM ATTACHED TO NITROGEN A CARBON?
         !
        if (nat(j) == 6) then
          loop4: do j2 = 1, nbonds(j)
            j3 = ibonds(j2, j)
               !
               !  IS THE ATOM ATTACHED TO THE CARBON THAT IS ATTACHED
               ! TO THE NITROGEN A CARBON?
               !
            if (nat(j3) == 6) then
              l = 0
              do j4 = 1, nbonds(j3)
                j5 = ibonds(j4, j3)
                     !
                     !  ARE ALL THE ATTACHED ATOMS, CARBON ATOMS?
                     !
                if (nat(j5) == 6) then
                  l = l + 1
                end if
                if (nat(j5) /= 1 .and. nat(j5) /= 6) cycle loop4
              end do
              if (l == 3) then
                k = 5
              end if
            end if
          end do loop4
        end if
      end do
    else if (k == 3) then
      !
      !    RESIDUE IS EITHER VAL OR PRO
      !
      !
      !  Check for both start of chain and (middle or end) of chain
      !
      l = 0
      do j = 1, nbonds(iatom)
        if (nat(ibonds(j, iatom)) /= 1) then
          l = l + 1
        end if
      end do
      if (l == 3 .and. iatom /= n1 .or. l == 2 .and. iatom == n1) then
        k = 20
      end if
    end if
    if (nat(n1) == 1) then  !  If the starting nitrogen atom is in fact a hydrogen,
      ninres = ninres + 1   !  then add it to the atoms in the residue.
      inres(ninres) = n1
    end if
    if (loc_tyres == " ") loc_tyres = tyres(k)
    if (k == 19 .and. jj3 == 34) then
      loc_tyres = "MSE"
    end if
    if (k == 23) then
      okay = .false.
      if (npack(6) == 2 .and. npack(7) == 0 .and. npack(8) == 1 .and. npack(16) == 0) then
        loc_tyres = "ACE"
        line = "Acetyl group"
        okay = .true.
      end if
      if (npack(6) == 1 .and. npack(7) == 1 .and. npack(8) == 0 .and. npack(16) == 0) then
        loc_tyres = "CH3"
        line = "Methyl group"
        okay = .true.
      end if
      if (okay .and. index(keywrd," RESID") /= 0) then
        ncomments = ncomments + 1
        write(all_comments(ncomments),'(a,i5)') &
               "*REMARK   2   "//loc_tyres//" = "//line(:47)//"res: ", ires
      end if
    end if
    if (len_trim(txeno(loop)) == 3) loc_tyres = txeno(loop)(:3)
    do j = 1, ninres
      l = inres(j)
      if (l /= 0) then
        write (txtatm(l), "(a,i6,1x,a3,a5,i6)")"ATOM ", l, cap_elemnt(nat(l)), loc_tyres, ires
      end if
    end do
    res_start(uni_res) = iatom
    call greek (jatom)
    if (k == 20) then            !  Check that N of proline is present
      do i = 1, ninres            !  If it is, and the nitrogen is in the
        if (inres(i) == jatom) exit    !  residue, remove the label
      end do
      if (i == ninres + 1 .and. jatom > 0) then
        txtatm(jatom) = " "
      end if
    end if
    if (txtatm(iatom)(:15) /= " ") then
      allres (ires) = loc_tyres // " "
      return
    end if
    if (ires == 1 .and.  txtatm(ihcr)(18:20) == "PRO") then
      write(txtatm(n1), "(I6,1X,A3, i4)") n1, "PRO", 1  ! Label N of proline in
    end if
    if (OXT /= 0) txtatm(OXT)(15:16) = "XT"
!
!  Run a check to make sure that all atoms are properly labeled.
!  If they are not properly labeled, re-label them as "unknown"
!
    do j2 = 1, 20
      if (afn(j2) == txtatm(inres(1))(18:20)) exit
    end do
    found = .false.
    UNK = .false.
    if (j2 < 21) then
      do i = 1, ninres
        l = inres(i)
        line = elemnt(nat(l))//txtatm(l)(15:15)//txtatm(l)(16:16)
        do j = 1, n_add(j2)
          if (line(:4) == atomname(j2,j)) found(j) = .true.
        end do
      end do
      do i = 5, n_add(j2)
        if (.not. found(i)) exit
      end do
!
!  Check iatom to make sure it's a backbone nitrogen
!
      if (j2 /= 8) then
        l = 0
        do j = 1, nbonds(iatom)
          if (nat(ibonds(j, iatom)) > 1) l = l + 1
        end do
        if (l > 2) i = 0
      end if
      if (i .le. n_add(j2)) then
        loc_tyres = "UNK"
        do i = 1, ninres
          txtatm(inres(i))(18:20) = loc_tyres
        end do
        UNK = .true.
      end if
    end if
    if (extra_atoms(1) > 0) then
!
! Assign labels to atoms that are associated with this residue, but call them "UNK"
!
      live(:2) = extra_atoms(:2)
      nlive = 1
      if (live(2) /= 0) nlive = 2
      ninres = 0
      outer_loop_extra: do
        l = live(1)
        if (l == jatom .or. ioptl(l) .or. l == iatom) then
           !
           !   Atom L is either already known to be in the residue, or is known
           !   to not be in the residue.
           !
          if (nlive == 0) go to 1101
          live(1) = live(nlive)
          nlive = nlive - 1
        else
          ioptl(l) = .true.
          ninres = ninres + 1
          if (ninres > natomr) then
            write (iw, "(A,I4,A,I4)") " There are more than", natomr, &
           & " atoms in residue", ires
            write (iw,*) " Atoms in residue", ires
            write (iw, "(10(1X,A2,I5))") (elemnt(nat(inres(l))), &
                 & inres(l), l=1, natomr)
            write(iw,'(/10x,a,/)')" Atoms recognized"
            do i = 1, numat
              if (txtatm(i) /= " ") write(iw,'(a)')txtatm(i)
            end do
            write(iw,'(/10x,a)')"(The fault is probably near the end of the list of recognized atoms)"
            call mopend ("Too many atoms in residue")
            go to 1000
          else
              !
              !  Assign atom L to the residue.
              !
            inres(ninres) = l
            if (nbonds(l) /= 0) then
                 !
                 !  THERE IS AT LEAST ONE ATOM ATTACHED TO THE 'LIVE' ATOM
                 !
              loop1_extra: do i2 = 2, nbonds(l)
                j = ibonds(i2, l)
                do i3 = 1, nlive
                  if (live(i3) == j) cycle loop1_extra
                end do
                nlive = nlive + 1
                    !
                    !   Panic if number of live atoms is greater than 500
                    !
                if (nlive > 500) exit outer_loop_extra
                live(nlive) = j
              end do loop1_extra
              live(1) = ibonds(1, l)
            else
              if (nlive == 0) go to 1101
              live(1) = live(nlive)
              nlive = nlive - 1
            end if
          end if
        end if
      end do outer_loop_extra
1101  continue
      tmp = "UNK"
      npack = 0
      do j = 1, ninres
        npack(nat(inres(j))) = npack(nat(inres(j))) + 1
      end do
      okay = .false.
      select case (npack(6))
      case(1)
        if (npack(7) == 1 .and. npack(8) == 0 .and. npack(16) == 0) then
          tmp = "CH3"
          line = "Methyl group"
          okay = .true.
        end if
      case(2)
        if (npack(7) == 0 .and. npack(8) == 1 .and. npack(16) == 0) then
          tmp = "ACE"
          line = "Acetyl group"
          okay = .true.
        end if
      case(4)
        if (npack(7) == 0 .and. npack(8) == 0 .and. npack(16) == 0) then
          tmp = "TBU"
          line = "Tertiary butyl group"
          okay = .true.
        end if
      end select
      if (okay .and. index(keywrd," RESID") /= 0) then
        ncomments = ncomments + 1
        write(all_comments(ncomments),'(a,i5)')  &
               "*REMARK   2   "//tmp//" = "//line(:47)//"res: ", ires
        do j = 1, ninres
          l = inres(j)
          if (l /= 0) then
            write (txtatm(l), "(a,i6,1x,a3,a5,i6)")"ATOM ", l, cap_elemnt(nat(l)), tmp, ires
          end if
        end do
      else
        do j = 1, ninres
          l = inres(j)
          if (l /= 0) then
            write (txtatm(l), "(a,i6,1x,a3,a5,i6)")"ATOM ", l, cap_elemnt(nat(l)), tmp, ires
          end if
        end do
      end if
    end if
    write (txtatm(iatom), "(a,i6,1x,a3,a5,i6)")"ATOM ", iatom, cap_elemnt(nat(iatom)), loc_tyres, ires
    if (index(keywrd, " ADD-H") /= 0) return
!
! residue 1 "by hand"
!
    i = Abs(mres(5,k) - nh)
    if (ico /= 0 .and. .not. UNK .and. i > 1 .and. .not. first_res .and. loop == 1) then
      if (nbonds(ico) /= 4 .or. nbonds(jofco) /= 2 .or. i /= 3) then  ! Exclude end group being -CH2OH instead of -COOH
        if(odd_h .and. k /= size_mres) then
          odd_h = .false.
          write(line,"(a, a)")"         Residues where the number of hydrogen", &
                       " atoms is not equal to that expected"
          write(iw,'(/,a,/)')trim(line)
          if (log) write(ilog,'(/,a,/)')trim(line)
          write(line,"(15x, a)")"Atom label              Coordinates of atom ""N""   No. of extra hydrogen atoms"
          write(iw,'(a,/)')trim(line)
          if (log) write(ilog,'(a)')trim(line)
        end if
        if (k /= size_mres) then
          do j = 1, numat_old
            if (abs(coord(1,iatom) - coorda(1,j)) > 0.1d0) cycle
            if (abs(coord(2,iatom) - coorda(2,j)) > 0.1d0) cycle
            if (abs(coord(3,iatom) - coorda(3,j)) > 0.1d0) cycle
            exit
          end do
          line = "          "//loc_tyres
          if (txtatm1(j) /= " ") line = txtatm1(j)
          write(line,"(6x, a, 4x, 3F8.3, i16)")'"'//line(:27)//'"', coord(:,iatom), nh - mres(5,k)
          write(iw,'(a)')trim(line)
          if (log) write(ilog,'(a)')trim(line)
        end if
      end if
    end if
    if (ires > maxres) then
      write (iw,*) " MODIFY SUBROUTINE NAMES"
      write (iw,*) " MAXIMUM NUMBER OF RESIDUES ALLOWED:", maxres
      call mopend ("Too many residues")
    end if
    allres (ires) = loc_tyres // " "
    do i = 1, numat
      if (ions(i) == 1 .and. Abs (lused(i)) == ires) then
        allres (ires) (4:4) = "+"
      end if
    end do
    j = 0
    l = 0
    do i = 1, numat
      if (ions(i) > 0 .and. Abs (lused(i)) == ires) j = j + ions(i)
      if (ions(i) < 0 .and. Abs (lused(i)) == ires) l = l + ions(i)
    end do
    if (j /= 0 .or. l /= 0) then
      if (j + l == 0) then
        if (j == 1) write (iw, "(8X,2A,I4,A)") " Residue:  '", allres(ires)(1:3), &
               & ires, "' is Zwitterionic"
      else if (j > 0) then
        allres (ires) (4:4) = "+"
      else
        allres (ires) (4:4) = "-"
      end if
    end if
    return
  end subroutine atomrs
  logical function peptide_n(l)
    use common_arrays_C, only : nat, ibonds, nbonds
    implicit none
    integer :: l
    integer :: i, j, k, nc, nh, nco
    if (nat(l) /= 7) then
      peptide_n = .false.
      return
    else
!
!  "l" is a nitrogen.  Check that it has three bonds, two to carbon and one to hydrogen,
!  and that one carbon is attached to an oxygen
!
      if (nbonds(l) /= 3) then
        peptide_n = .false.
        return
      end if
      nc = 0
      nh = 0
      nco = 0
      do i = 1, 3
        if (nat(ibonds(i,l)) == 6)then
          nc = nc + 1
          k = ibonds(i,l)
          if (nbonds(k) == 3) then
            do j = 1, 3
              if (nat(ibonds(j,k)) == 8) then
                if (nbonds(ibonds(j,k)) == 1) then
                  nco = nco + 1
                else
                  peptide_n = .false.         ! Structure N-C-O-R detected
                  return
                end if
              end if
            end do
          else if (nbonds(k) == 4) then
            do j = 1, 4
              if (nat(ibonds(j,k)) == 8) then ! Structure N-C-O-R detected
                peptide_n = .false.
                return
              end if
            end do
          end if
        end if
        if (nat(ibonds(i,l)) == 1) nh = nh + 1
      end do
      if (nc /= 2 .or. nh /= 1 .or. nco /= 1) then
        peptide_n = .false.   ! Structure C-NH-C=O not detected
        return
      end if
!
!  "l" is a nitrogen.  Check that it has three bonds, two to carbon and one to hydrogen
!
      peptide_n = .true.
      return
    end if
  end function peptide_n
