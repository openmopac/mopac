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

subroutine ligand (ires, start_res, nfrag)
    use elemts_C, only: cap_elemnt
    use common_arrays_C, only : txtatm, txtatm1, nat, labels, nbonds, ibonds, &
      all_comments
    use MOZYME_C, only : allres
    use molkst_C, only : natoms, id, ncomments, line, keywrd, numat
    use chanel_C, only: log, ilog, iw
    implicit none
    integer, intent (inout) :: ires, nfrag
    integer, intent (in) :: start_res(*)
!
    integer, parameter :: maxlive = 15000
    integer :: i, j, k, l, m, n, ii, jj, kk, mm, nn, res, change_no, max_comments, &
      nlive, n_C, n_O, n_N, n_S, n_P, n_ele(107), ninres, i_panic, ele_order(107)
    logical, allocatable :: l_used(:)
    integer, allocatable :: inres(:)
    integer, allocatable :: live(:)
    logical :: attached, l_chain, first, l_write = .true.
    integer, external :: nheavy
    double precision, external :: distance, reada
    character :: het*3, het_group*120, num, num1, el*2, het2*3, line1*120
    save :: l_write
    allocate (l_used(natoms), inres(natoms), live(maxlive))
    ele_order(1) = 6
    ele_order(2) = 1
    ele_order(3) = 7
    ele_order(4) = 8
    k = 4
    do i = 2, 107
       select case (i)
       case (6, 7, 8)
       case default
         k = k + 1
         ele_order(k) = i
       end select
    end do
    j = 0
    do i = 1, ncomments
      if (index(all_comments(i), "REMARK   3") == 0) then
        j = j + 1
        all_comments(j) = all_comments(i)
      end if
    end do
    ncomments = j
    max_comments = size(all_comments)
    first = .true.
    change_no = 0
!
    i_loop: do i = 1, natoms - id
      if (txtatm(i) (26:26) == " ") then
!
!  Residue number is missing, therefore not a residue.
!
!
!    Test for specific molecules
!
        if (labels(i) == 15) then
          if (nbonds(i) == 4) then
            k = 0
            do j = 1, 4
              if (labels(ibonds(j, i)) /= 8) go to 1000
              if (nbonds(ibonds(j,i)) > 1) k = k + 1
            end do
            if (k > 1) go to 1000
!
!                                 Phosphate
!
            call inc_res(ires, start_res, nfrag)
            write (txtatm(i), "(a6,i5,1x,a3,a5,i6)")"HETATM", i, " P ", "PO4", ires
            do j = 1, 4
              write (txtatm(ibonds(j, i)), "(a6,i5,1x,a3,a5,i6)")"HETATM", ibonds(j, i), " O ", "PO4", ires
            end do
          end if
        else if (labels(i) == 16) then
          if (nbonds(i) == 4) then
            do j = 1, 4
              if (labels(ibonds(j, i)) /= 8) go to 1000
            end do
!
!                                 Sulfate
!
            call inc_res(ires, start_res, nfrag)
            write (txtatm(i), "(a6,i5,1x,a3,a5,i6)")"HETATM", i, " S ", "SO4", ires
            do j = 1, 4
              write (txtatm(ibonds(j, i)), "(a6,i5,1x,a3,a5,i6)")"HETATM", ibonds(j, i), " O ", "SO4", ires
            end do
          end if
        else if (labels(i) == 8 .and. nheavy(i) == 0) then
!
!                                 Water
!
          call inc_res(ires, start_res, nfrag)
          write (txtatm(i), "(a6,i5,1x,a3,a5,i6)")"HETATM", i, " O ", "HOH", ires
            do j = 1, nbonds(i)
              write (txtatm(ibonds(j, i)), "(a6,i5,1x,a3,a5,i6)")"HETATM", ibonds(j, i), " H ", "HOH", ires
            end do
        else if (labels(i) == 8) then
!
!                                     Ethylene glycol and glycerol
!                                     H-O-CH2-CH2-O-H and HOCH2-HCOH-H2COH
!
          if (nheavy(i) == 1) then
            do ii = 1, nbonds(i)
              j = ibonds(ii,i)
              if (nat(j) /= 6) cycle
              if (nheavy(j) /= 2) cycle! Oxygen joined to carbon(1)
                do jj = 1, nbonds(j)
                k = ibonds(jj,j)
                if (nat(k) /= 6) cycle
                if (nheavy(k) == 1) then
!
!  Found methanol
!
                else if (nheavy(k) == 2) then! Carbon(1) joined to oxygen and carbon(2)

                  do kk = 1,nbonds(k)
                    l = ibonds(kk,k)
                    if (nat(l) == 8 .and. l /= j) then
                      l = l
                    end if
                    if (nat(l) /= 8) cycle
                    if (nheavy(l) == 1) then! Oxygen joined to carbon(2)
!
! Found O-C-C-O
!
                      call inc_res(ires, start_res, nfrag)
                      write (txtatm(i), "(a6,i5,1x,a3,a5,i6)")"HETATM", i, " O ", "EDO", ires
                      write (txtatm(j), "(a6,i5,1x,a3,a5,i6)")"HETATM", i, " C ", "EDO", ires
                      write (txtatm(k), "(a6,i5,1x,a3,a5,i6)")"HETATM", i, " C ", "EDO", ires
                      write (txtatm(l), "(a6,i5,1x,a3,a5,i6)")"HETATM", i, " O ", "EDO", ires
                    end if
                  end do
                else if (nheavy(k) == 3) then

                  do kk = 1,nbonds(k)
                    l = ibonds(kk,k)
                    if (nat(l) /= 8) cycle
                    if (nheavy(l) /= 1) cycle! Oxygen joined to carbon(2)
                    do mm = 1,nbonds(k)
                      m = ibonds(mm,k)
                      if (nat(m) /= 6) cycle
                      if (m == j) cycle
                      if (nheavy(m) /= 2) cycle! Carbon(3) joined to carbon(2)
                      do nn = 1, nbonds(m)
                        n = ibonds(nn,m)
                        if (nat(n) /= 8) cycle
!
! Found O-C-CO-C-O
!
                        call inc_res(ires, start_res, nfrag)
                        write (txtatm(i), "(a6,i5,1x,a3,a5,i6)")"HETATM", i, " O ", "GOL", ires
                        write (txtatm(j), "(a6,i5,1x,a3,a5,i6)")"HETATM", i, " C ", "GOL", ires
                        write (txtatm(k), "(a6,i5,1x,a3,a5,i6)")"HETATM", i, " C ", "GOL", ires
                        write (txtatm(l), "(a6,i5,1x,a3,a5,i6)")"HETATM", i, " O ", "GOL", ires
                        write (txtatm(m), "(a6,i5,1x,a3,a5,i6)")"HETATM", i, " C ", "GOL", ires
                        write (txtatm(n), "(a6,i5,1x,a3,a5,i6)")"HETATM", i, " O ", "GOL", ires
                      end do
                    end do
                  end do
                end if
              end do
            end do
          end if
        end if
1000    continue
      end if
      n_ele = 0
      ninres = 0
      i_panic = 0
      l_used = .false.
      if (txtatm(i) (26:26) == " " .and. nat(i) /= 1) then
        ninres = 1
        inres(1) = i
        l_used(i) = .true.
        nlive = nbonds(i)
        live(:nlive) = ibonds(:nlive,i)
        select case (nat(i))
        case (6:9, 15:17, 34:35, 53)
        case default
          nlive = 0
        end select
        call inc_res(ires, start_res, nfrag)
        outer_loop: do
          l = live(1)
          if (nlive == 0) exit
          if (l == 0) exit
          if (l_used(l) .or. txtatm(l)(11:14) /= " ") then
!
!   Atom L is already labeled
!
            live(1) = live(nlive)
            nlive = nlive - 1
          else
!
!  Assign atom L to the hetero group.
!
            l_used(l) = .true.
            do ii = 1, nbonds(l)
              jj = ibonds(ii,l)
              if (txtatm(jj)(11:14) /= " ") then
                select case (nat(i))
                case (6:9, 15:17, 34:35, 53)
                  i_panic = jj
                end select
              end if
            end do
            select case (nat(i))
            case (6:9, 15:17, 34:35, 53)
              ninres = ninres + 1
              inres(ninres) = l
            end select
            if (nbonds(l) /= 0) then
!
!  THERE IS AT LEAST ONE ATOM ATTACHED TO THE 'LIVE' ATOM
!
              loop1: do ii = 2, nbonds(l)
                j = ibonds(ii, l)
                do jj = 1, nlive
                  if (live(jj) == j) cycle loop1
                end do
                if (.not. l_used(j) .and. txtatm(j)(11:14) == " ") then
                  nlive = nlive + 1
                  if (nlive > maxlive) return
                  live(nlive) = j
                end if
              end do loop1
              live(1) = ibonds(1, l)
            else
              if (nlive == 0) exit
              live(1) = live(nlive)
              nlive = nlive - 1
            end if
          end if
        end do outer_loop
        if (i_panic == 0) then
          n_ele = 0
          l = 0
          do j = 1, ninres
            k = inres(j)
            if (k <= numat) then
              if (nat(k) < 99 .or. (nat(k) > 99 .and. nat(k) < 107)) then
                l = l + 1
                n_ele(nat(k)) = n_ele(nat(k)) + 1
              end if
            end if
          end do
          if (l == 0) return
          ninres = l
          n_C = n_ele(6)
          n_N = n_ele(7)
          n_O = n_ele(8)
          n_S = n_ele(16)
          n_P = n_ele(15)
          attached = .false.
          res = -10000
          n = 0
          loop_k: do k = 1, ninres
            m = inres(k)
            do l = 1, nbonds(m)
              n = ibonds(l,m)
              if (txtatm(n)(26:26) /= " ") then
                attached = .true.
                res = nint(reada(txtatm(n), 23))
                exit loop_k
              end if
            end do
          end do loop_k
          if (n > 0 .and. res /= -10000) then
            if (txtatm(n)(13:16) == " N  " .and. nat(m) == 6) then
!
!  Inserting a fragment at the start of a chain -
!  so move all residue numbers in the rest of the chain up by 1
!
              k = res
              do
                j = nint(reada(txtatm(n), 23))
                if (j == k) then
                  write(txtatm(n)(23:26),'(i4)') k + 1
                else if (j == k + 1) then
                  k = k + 1
                  write(txtatm(n)(23:26),'(i4)') k + 1
                else
                  exit
                end if
                n = n + 1
              end do
              k = res
              het  = allres(k)(:3)
              allres(k) = "UNK"
              do
                het2 = allres(k + 1)(:3)
                if (het2 == " ") exit
                allres(k + 1) = het
                het = het2
                k = k + 1
              end do
            end if
          end if
!
!  Formula of hetero group is C(n_C) N(n_N) O(n_O)
!  and atoms are in inres(:ninres)
!
!  Now assign "obvious" hetero groups
!
          het = "HET"
          line = "Not identified;"
          do m = 1, 84
            k = ele_order(m)
            if (n_ele(k) > 0) then
              l = n_ele(k)
              num = char(ichar("1") +int(log10(l + 0.05)))
              if (cap_elemnt(k)(1:1) == " ") then
                num1 = "2"
              else
                num1 = "3"
              end if
              if (l == 1) then
                write(line(len_trim(line) + 1:),'(a'//num1//')')cap_elemnt(k)
              else
                write(line(len_trim(line) + 1:),'(a'//num1//',i'//num//')')cap_elemnt(k), l
              end if
            end if
          end do
          write(line(35:),'(a, i6, a)')" =", ninres, " atoms"
          het_group = trim(line)
!
!  Make sure all oxygen atoms are counted
!
          kk = ninres
          do ii = 1, kk
            k = inres(ii)
            if (nat(k) /= 6) cycle
            jj = ibonds(nbonds(k) + 1,k)
            if (jj > 0) then
              if (nat(jj) == 8) then
                if (distance(k,jj) < 1.9d0) then
                  do k = 1, ninres
                    if (inres(k) == jj) exit
                  end do
                  if (k > ninres) then
                    n_O = n_O + 1
                    ninres = ninres + 1
                    inres(ninres) = jj
                  end if
                end if
              end if
            end if
          end do
          select case (n_C)
            case (0)
              j = 0
              do k = 1, ninres
                if (nat(inres(k)) /= 1) then
                  j = j + 1
                  l = inres(k)
                end if
              end do
              if (ninres == 1) then
                j = 0
                select case (nat(i))
                  case(3:4, 11:12, 19:30, 37:48, 55:80)
                    j = 1
                end select
                if (j == 1 .or. nbonds(i) == 0) then
                  het = " "//cap_elemnt(nat(inres(1)))
                  if (het(3:3) >= "a" .and. het(3:3) <= "z") het(3:3) = char(ichar(het(3:3)) + ichar("A") - ichar("a"))
                  het_group = "Isolated element"
                else if (nbonds(i) > 0) then
                  het = " "//cap_elemnt(nat(inres(1)))
                  if (het(3:3) >= "a" .and. het(3:3) <= "z") het(3:3) = char(ichar(het(3:3)) + ichar("A") - ichar("a"))
                  het_group = "Covalently bound element"
                end if
                k = ibonds(1,i)
                if (nbonds(i) == 1) then
                  if (txtatm(k) /= " ") then
                    el = cap_elemnt(nat(i))
                    if (el(2:2) == " ") el = " "//el(1:1)
                    if (nat(i) == 1) then
                      txtatm(i) = txtatm(k)(:12)//" H"//txtatm(k)(15:)
                    else
!
!  A single element is bonded to something.
!  If possible, use the original label (in txtatm1), if not, then use the label of the atom it's attached to.
!
                      if (txtatm1(i) == " ") then
                        j = nint(reada(txtatm(k), 21)) + 1
                        write(line,'(8x,i4)')j
                      else
                        line =txtatm1(i)(15:)
                      end if
                      if(cap_elemnt(nat(i))(2:2) == " ") then
                        txtatm(i) = txtatm(k)(:12)//" "//cap_elemnt(nat(i))(1:1)//trim(line)
                      else
                        txtatm(i) = txtatm(k)(:12)//cap_elemnt(nat(i))//txtatm(k)(15:)
                      endif
                    end if
                    cycle
                  end if
                end if
              end if
              if (ninres == 2 .and. n_O == 1 .and. n_ele(1) == 1) then
                het = " OH"
                het_group = "Hydroxide ion"
              end if
              if (ninres == 3 .and. n_O == 1 .and. n_ele(1) == 2) then
                het = "H2O"
                het_group = "Complexed water"
              end if
              if (nat(l) == 7) then
                if (nheavy(l) > 1) then
                  do k = 1, nbonds(l)
                    j = ibonds(k,l)
                    if (nat(j) /= 6) cycle
!
!  Found carbon attached to nitrogen
!
                      do m = 1, nbonds(j)
                        if (nat(ibonds(m,j)) == 8) exit
                      end do
                      if (m <= nbonds(j)) then
                        write(txtatm(l),"(a6,i5,a3,a12)")txtatm(j)(:6), l, "  N", txtatm(j)(15:)
                        cycle i_loop
                      end if
                  end do
                else
                  het = "NH3"
                  het_group = "Ammonia"
                end if
              end if
              if  (n_N == 0 .and. n_O == 2 .and. n_S == 0 .and. n_P == 0) then
                het = "O-O"
                het_group = "Hydrogen peroxide?"
              end if
             if  (n_N == 1 .and. n_O == 3 .and. n_S == 0 .and. n_P == 0) then
                het = "NO3"
                het_group = "Nitrate"
              end if
             if  (n_P == 1 .and. n_O == 4 .and. n_S == 0) then
                het = "PO4"
                het_group = "Phosphate"
              end if
            case (1)
              if  (n_N == 0 .and. n_O == 1 .and. n_S == 0 .and. n_P == 0) then
                het = "MOH"
                het_group = "Methanol"
              end if
              if  (n_N == 0 .and. n_O == 2 .and. n_S == 0 .and. n_P == 0) then
                het = "FMT"
                het_group = "Formic acid"
              end if
              if  (n_N == 0 .and. n_O == 3 .and. n_S == 0 .and. n_P == 0) then
                het = "CO3"
                het_group = "Carbonate ion"
              end if
              if  (n_N == 1 .and. n_O == 0 .and. n_S == 0 .and. n_P == 0) then
                het = "CYN"
                het_group = "Cyanide ion"
              end if
              if  (n_N == 0 .and. n_O == 0 .and. n_ele(1) == 3 .and. n_S == 0 .and. n_P == 0) then
                het = "CH3"
                het_group = "Methyl group"
              end if
              if  (n_N == 2 .and. n_O == 1 .and. n_ele(1) == 4 .and. n_S == 0 .and. n_P == 0) then
                het = "URE"
                het_group = "UREA"
              end if
            case (2)
              if  (n_N == 0 .and. n_O == 1 .and. n_S == 0 .and. n_P == 0) then
                het = "EOH"
                het_group = "Ethanol"
              else if (n_N == 0 .and. n_O == 2 .and. n_S == 0 .and. n_P == 0) then
                do ii = 1, ninres
                  if (nat(inres(ii)) == 6) exit
                end do
                if (nheavy(inres(ii)) == 2) then
                  het = "EDO"
                  het_group = "1,2-Ethanediol "
                else
                  het = "ACY"
                  het_group = "Acetic acid "
                end if
              end if
            case (3)
              if  (n_N == 0 .and. n_O == 2 .and. n_S == 0 .and. n_P == 0) then
                het = "PGO"
                het_group = "S-1,2-Propanediol"
              end if
              if  (n_N == 0 .and. n_O == 3 .and. n_S == 0 .and. n_P == 0) then
                het = "GOL"
                het_group = "Glycerol"
              end if
              if  (n_N == 0 .and. n_O == 4 .and. n_S == 0 .and. n_P == 0) then
                het = "MLI"
                het_group = "Malonate ion (-)"
              end if
              if  (n_N == 2 .and. n_O == 0 .and. n_S == 0 .and. n_P == 0) then
                het = "IMD"
                het_group = "Imidazole"
              end if
            case (4)
              if  (n_N == 0 .and. n_O == 3 .and. n_S == 0 .and. n_P == 0) then
                het = "PEG"
                het_group = "Di(hydroxyethyl) ether"
              end if
              if  (n_N == 0 .and. n_O == 6 .and. n_S == 0 .and. n_P == 0) then
                het = "TLA"
                het_group = "L(+)-Tartaric acid"
              end if
              if  (n_N == 1 .and. n_O == 0 .and. n_S == 0 .and. n_P == 0) then
                het = "NTB"
                het_group = "N-Tertiary butyl"
              end if
              if  (n_N == 1 .and. n_O == 3 .and. n_S == 0 .and. n_P == 0) then
                het = "TRS"
                het_group = "2-Amino-2-hydroxymethyl-propane-1,3-diol"
              end if
            case (5)
              if  (n_N == 1 .and. n_O == 3 .and. n_S == 0 .and. n_P == 0) then
                het = "PCA"
                het_group = "Pyroglutamic acid"
              end if
              if  (n_N == 1 .and. n_O == 7 .and. n_S == 0 .and. n_P == 0) then
                het = "10E"
                het_group = "4-Amino-3-methylbut-2-en-1-yl diphosphate"
              end if
            case (6)
               if  (n_N == 0 .and. n_O == 1 .and. n_P == 0) then
                het = "IPH"
                het_group = "Phenol"
              end if
               if  (n_N == 0 .and. n_O == 1 .and. n_P == 1) then
                het = "HXP"
                het_group = "Hexyl bonded to phosphorus"
              end if
              if  (n_N == 0 .and. n_O == 2 .and. n_S == 0 .and. n_P == 0) then
                het = "RCO"
                het_group = "Resorcinol"
              end if
              if  (n_N == 0 .and. n_O == 9 .and. n_P == 0) then
                het = "SGA"
                het_group = "O3-Sulfonylgalactose"
              end if
              if  (n_N == 0 .and. n_O == 5 .and. n_S == 0 .and. n_P == 0) then
                het = "FUC"
                het_group = "Alpha-L-fucose"
              end if
              if  (n_N == 0 .and. n_O == 6 .and. n_S == 0 .and. n_P == 0) then
                het = "HEX"
                het_group = "A hexose, e.g. glucose or mannose"
                call identify_hexose(ninres, inres, het, het_group)
              end if
              if  (n_N == 0 .and. n_O == 2 .and. n_S == 0 .and. n_P == 0) then
                het = "MPD"
                het_group = "(4S)-2-Methyl-2,4-pentanediol"
              end if
              if  (n_N == 0 .and. n_O == 4 .and. n_S == 0 .and. n_P == 0) then
                het = "PGE"
                het_group = "Triethylene glycol"
              end if
              if  (n_N == 1 .and. n_O == 4 .and. n_P == 0) then
                het = "MES"
                het_group = "2-(N-Morpholino)-ethanesulfonic acid"
              end if
            case (7)
              if  (n_N == 1 .and. n_O == 1 .and. n_S == 0 .and. n_P == 0) then
                het = "DBZ"
                het_group = "Benzoylamino (H6C7NO)"
              end if
              if  (n_N == 1 .and. n_O == 3 .and. n_S == 0 .and. n_P == 0) then
                het = "TAM"
                het_group = "Tris(hydroxyethyl) aminomethane"
              end if
            case (8)
              if  (n_N == 0 .and. n_O == 0 .and. n_S == 0 .and. n_P == 0) then
                het = "OCT"
                het_group = "N-Octane"
              end if
              if  (n_N == 0 .and. n_O == 5 .and. n_S == 0 .and. n_P == 0) then
                het = "PG4"
                het_group = "Tetraethylene glycol"
              end if
              if  (n_N == 0 .and. n_O == 1 .and. n_S == 0 .and. n_P == 0) then
                het = "0VT"
                het_group = "6-Methyl-5-hepten-2-one"
              end if
              if  (n_N == 1 .and. (n_O > 3 .and. n_O < 7) .and. n_S == 0 .and. n_P == 0) then
                het = "NAG"
                het_group = "N-Acetyl-D-glucosamine"
              end if
              if  (n_N == 1 .and. n_O == 0 .and. n_S == 0 .and. n_P == 0) then
                het = "TBA"
                het_group = "Tetrabutylammonium ion"
              end if
              if  (n_N == 1 .and. n_O == 2 .and. n_P == 0) then
                het = "AES"
                het_group = "4-(2-Aminoethyl)benzenesulfonyl fluoride"
              end if
              if  (n_N == 2 .and. n_O == 4 .and. n_P == 0) then
                het = "EPE"
                het_group = "Hepes (C8 H18 N2 O4 S)"
              end if
            case (9)
              if  (n_N == 1 .and. n_O == 2 .and. n_ele(17) == 1.and. n_S == 0 .and. n_P == 0) then
                het = "GM5"
                het_group = "4-Chlorocinnamylhydroxamate"
              end if
              if  (n_N == 3 .and. n_O == 13 .and. n_S == 0 .and. n_P == 3) then
                het = "DCP"
                het_group = "2'-deoxycytidine-5'-triphosphate"
              end if
              if (n_N == 3 .and. n_O == 14 .and. n_S == 0 .and. n_P == 3) then
                het = "CTP"
                het_group = "Cytidine-5'-triphosphate"
              end if
              if (n_N == 2 .and. n_O == 15 .and. n_S == 0 .and. n_P == 3) then
                het = "UTP"
                het_group = "Uridine-5'-triphosphate"
              end if
            case (10)
              if  (n_N == 0 .and. n_O == 6 .and. n_S == 0 .and. n_P == 0) then
                het = "TSA"
                het_group = "Endo-oxabicyclic transition state analogue"
              end if
              if  (n_N == 1 .and. n_O == 1 .and. n_S == 0 .and. n_P == 0) then
                het = "T55"
                het_group = "8-Methylnonanoic acid"
              end if
              if (n_N == 1 .and. n_O == 3 .and. n_S == 0) then
                het = "VPF"
                het_group = "C22 H32 N3 O8 P (incomplete)"
              end if
              if (n_N == 2 .and. n_O == 14 .and. n_S == 0 .and. n_P == 3) then
                het = "TTP"
                het_group = "Thymidine-5'-triphosphate"
              end if
              if (n_N == 5 .and. n_O == 10 .and. n_S == 0 .and. n_P == 2) then
                het = "ADP"
                het_group = "Adenosine-5'-diphosphate"
              end if
              if (n_N == 5 .and. n_O == 13 .and. n_S == 0 .and. n_P == 3) then
                het = "ATP"
                het_group = "Adenosine-5'-triphosphate"
              end if
              if (n_N == 5 .and. n_O == 14 .and. n_S == 0 .and. n_P == 3) then
                het = "GTP"
                het_group = "Guanosine-5'-triphosphate"
              end if
              if  (n_N == 6 .and. n_O == 12 .and. n_S == 0) then
                het = "ANP"
                het_group = "Phosphoaminophosphonic acid-adenylate ester"
              end if
            case (13)
              if  (n_N == 2 .and. n_O == 0 .and. n_S == 0 .and. n_P == 0) then
                het = "THA"
                het_group = "Tacrine"
              end if
            case (14)
              if  (n_N == 0 .and. n_O == 4 .and. n_S == 0 .and. n_P == 0) then
                het = "HKA"
                het_group = "3-Methoxy-4-phenoxybenzoic acid"
              end if
            case (15)
              if  (n_N == 1 .and. n_O == 5 .and. n_S == 0) then
                het = "1CT"
                het_group = "Oxo-trifluoro methyl phenyl ethoxy "// &
                  &"phenyl phosphonic acid"
              end if
              if  (n_N == 2 .and. n_O == 12 .and. n_S == 0 .and. n_P == 2) then
                het = "5GW"
                het_group = "5-Phenyluridine 5'-(trihydrogen diphosphate)"
              end if
            case (16)
              if  (n_N == 2 .and. n_O == 11 .and. n_S == 0 .and. n_P == 0) then
                het = "NA2"
                het_group = "N-Acetyl-D-glucosamine, dimer"
              end if
              if  (n_N == 0 .and. n_O == 1 .and. n_S == 0 .and. n_P == 0) then
                het = "BOM"
                het_group = "Hexadecyl-10,12-dien-1-ol"
              end if
            case (17)
              if (n_N == 4 .and. n_O == 9 .and. n_S == 0 .and. n_P == 1) then
                het = "FMN"
                het_group = "Riboflavin monophosphate "
              end if
              if (n_N == 2 .and. n_O == 2) then
                het = "641"
                het_group = "Oxopyrrolidine-3-carboxamide"
              end if
            case (18)
              if (n_N == 2 .and. n_O == 0) then
                het = "HUX"
                het_group = "(-)-Huprine X (C18 H19 N2 Cl)"
              end if
            case (19)
              if (n_N == 0 .and. n_O == 2) then
                het = "ASD"
                het_group = "4-Androstene-3-17-dione"
              end if
              if (n_N == 4 .and. n_O == 2) then
                het = "PNT"
                het_group = "1,5-Bis(4-amidinophenoxy)pentane"
              end if
            case (20)
              if (n_N == 0 .and. n_O == 1) then
                het = "ARC"
                het_group = "3,7,11,15-Tetramethyl-hexadecan-1-ol"
              end if
              if (n_N == 0 .and. n_O == 10) then
                het = "BHE"
                het_group = "Galactopyranoside"
              end if
              if (n_N == 4 .and. n_O == 2) then
                het = "DID"
                het_group = "4,4'[1,6-Hexanediylbis(oxy)]"// &
                &"bisbenzenecarboximidamide"
              end if
              if (n_N == 7 .and. n_O == 9) then
                het = "BT5"
                het_group = "Biotinyl-5-amp"
              end if
            case (21)
              if (n_N == 2 .and. n_O == 17 .and. n_P == 2) then
                het = "2GW"
                het_group = "5-Phenyl-uridine-5'-alpha-D-galactosyl-diphosphate"
              end if
              if (n_N == 7 .and. n_O == 14) then
                het = "NAD"
                het_group = "Nicotinamide-adenine-dinucleotide"
              end if
              if (n_N == 7 .and. n_O == 17 .and. n_P > 0) then
                het = "NAP"
                het_group = "Nicotinamide-adenine-dinucleotide phosphate"
              end if
              if (n_N == 10 .and. n_O == 1) then
                het = "32G"
                het_group = "C21 H30 N10 O S"
              end if
            case (22)
              if (n_N == 2 .and. n_O == 2) then
                het = "CZM"
                het_group = "3,3'-Me2-salophen"
              end if
              if (n_N == 3 .and. n_O == 8) then
                het = "VPF"
                het_group = "C22 H32 N3 O8 P"
              end if
            case (23)
              if (n_N == 0 .and. n_O == 11) then
                het = "CM5"
                het_group = "5-Cyclohexyy-1-pentyl-beta-d-maltoside"
              end if
              if (n_N == 3 .and. n_O == 3) then
                het = "4A2"
                het_group = "C23 H17 F4 N3 O3"
              end if
            case (26)
              if (n_N == 0 .and. n_O == 8) then
                het = "0DV"
                het_group = "Fusicoccin H"
              end if
            case (34)
              if (n_N == 4 .and. n_O == 4) then
                het = "HEM"
                het_group = "Heme ring"
              end if
          end select
          if (n_ele(15) > 3) then
            het = "NUC"
            het_group = "Nucleic acid (DNA or RNA type)"
          end if
          l = index(keywrd, " XENO")
          if (l /= 0) then
            l_chain = (index(keywrd, " CHAINS=(") == 0)
            j = index(keywrd(l:), ") ")
            if (j /= 0) then
              line = " "//keywrd(l + 5:j + l - 1)
              do
                do l = 1, 10
                  line = trim(line(2:))
                  if (line(1:1) == "(" .or. line(1:1) == "," .or. line(1:1) == ";" .or. line(1:1) == ")") exit
                end do
                if (line == " ") exit
                num = line(2:2)
                k = 0
                if (l_chain) then
                  if (num >= "A" .and. num <= "Z") then
                    if (txtatm1(i)(22:22) >= "A" .and. txtatm1(i)(22:22) <= "Z") then
                      if (num /= txtatm1(i)(22:22)) k = -1000
                    end if
                  end if
                end if
                k = nint(reada(line,1)) + k
                if (k == ires) then
                  do l = 1, 10
                    if (line(1:1) == "=") exit
                    line = trim(line(2:))
                  end do
                  if (het_group(1:3) /= "Not") then
                    line1 = "Defined using keyword XENO as "//trim(het_group)
                    het_group = trim(line1)
                  end if
                  if (l_write) then
                    if (first) then
                      first = .false.
                      write(iw,'(/,a)') "      Ligand names that have been changed by XENO keyword"
                      write(iw,'(/,a)') "      Residue No.  Calculated name   XENO name"
                    end if
                      change_no = change_no + 1
                      write(iw,'(i3,i9,3x,a1,8x,a3,13x,a3)')change_no, ires, num, het, line(2:4)
                      het = line(2:4)
                      exit
                    end if
                  end if
              end do
            end if
          end if
          if (nat(i) /= 1 .or. nbonds(i) /= 1) then
            if (ncomments < max_comments .and. index(keywrd, " RESID") /= 0) then
              ncomments = ncomments + 1
              if (index(het_group, "XENO") == 0 ) then
                write(line,'(a,i4)') "*REMARK   3 "//het//" = "//trim(het_group)//" res: ", ires
              else
                write(line,'(a,i4)') "*REMARK   3 "//het//" = "//trim(het_group)//"   res: ", ires
              end if
              j = len_trim(line)
              if (j > 80) then
                j = index(line, "XENO")
                write(all_comments(ncomments),'(a)') line(1:j +7)
                ncomments = ncomments + 1
                write(all_comments(ncomments),'(a)') "*REMARK   3"//trim(line(j + 7:))
              else
                write(all_comments(ncomments),'(a)') trim(line)
              end if
              if (log) then
                 write(ilog,'(a)')  trim(line)
                 if (het == "HET")  write(ilog,'(20x,a,i5,a)')" (Includes atom number:", inres(1),")"
              end if
            end if
            do j = 1, ninres
              k = inres(j)
              el = cap_elemnt(nat(k))
              if (el(2:2) == " ") el = " "//el(1:1)
              if (attached) then
                write (txtatm(k), "(a6,i5,a3,a6,i6)")"ATOM  ", k, " "//el, het, res
              else
                write (txtatm(k), "(a6,i5,a3,a6,i6)")"HETATM", k, " "//el, het, ires
              end if
            end do
          end if
        else
          do j = 1, ninres
            k = inres(j)
            el = cap_elemnt(nat(k))
            if (el(2:2) == " ") el = " "//el(1:1)
            write(txtatm(k),"(a6,i5,a3,a6,i6)")"HETATM", k, " "//el, "UNK", ires
          end do
          call inc_res(ires, start_res, nfrag)
        end if
      end if
    end do i_loop
!
!  Label any remaining atoms
!
    do i = 1, natoms - id
      if (txtatm(i) (26:26) == " ") then
        el = cap_elemnt(nat(i))
        if (el(2:2) == " ") el = " "//el(1:1)
        if (el == " H" .and. i > 1 .and. nbonds(i) == 1) then
          j = ibonds(1,i)
          write(txtatm(i),"(a6,i5,a3,a6)")txtatm(j)(:6), j, " "//el, txtatm(j)(18:20)
        else
         write(txtatm(i),"(a6,i5,a3,a6)")"HETATM", i, " "//el, "UNK"
        end if
      end if
    end do
    l_write = .false.
    return
end subroutine ligand
subroutine moiety (iopt, lused, istart, new)
    use molkst_C, only: numat, natoms
    use common_arrays_C, only : nat, nbonds, ibonds
    use chanel_C, only: iw
    use elemts_C, only: elemnt
!
!.. Implicit Declarations ..
    implicit none
!
!.. Formal Arguments ..
    integer, intent (in) :: istart
    integer, intent (inout) :: new
    logical, dimension (numat), intent (inout) :: iopt
    integer, dimension (natoms), intent (inout) :: lused
!
!.. Local Scalars ..
    integer :: i2, i3, iatom, j, j2, k, l, ninbit, nlive
!
!.. Local Arrays ..
    integer, allocatable :: live(:), inres(:)
!
! ... Executable Statements ...
!
   i2 = max(2000, natoms*2)
    allocate (live(i2), inres(i2))
    inres = 0
    iatom = istart
!
!   The first atom identified in the moiety is atom IATOM.
!
    iopt(iatom) = .true.
!
!  NOW TO WORK OUT THE ATOMS IN THE MOIETY
!
    nlive = nbonds(iatom)
    if (nlive == 0) then
      ninbit = 1
      inres(1) = iatom
    else
      ninbit = 0
      do i2 = 1, nlive
        live(i2) = ibonds(i2, iatom)
      end do
      do
        l = live(1)
        if (iopt(l) .or. l == iatom) then
          if (nlive < 1) exit
          live(1) = live(nlive)
          nlive = nlive - 1
        else
          iopt(l) = .true.
          ninbit = ninbit + 1
          if (ninbit > natoms) then
            write (iw, "(A,I4,A)") " There are more than", natoms, &
                 & " atoms in moiety "
            write (iw,*) " Atoms in moiety"
            write (iw, "(10(1X,A2,I5))") (elemnt(nat(inres(l))), inres(l), &
                 & l=1, natoms)
            call mopend ("Too many atoms in moiety")
            return
          end if
          inres(ninbit) = l
          if (nbonds(l) /= 0) then
!
!  THERE IS AT LEAST ONE ATOM ATTACHED TO THE 'LIVE' ATOM
!
            loop: do i2 = 2, nbonds(l)
              j = ibonds(i2, l)
              do i3 = 1, nlive
                if (live(i3) == j) cycle loop
              end do
              nlive = nlive + 1
              live(nlive) = j
            end do loop
            live(1) = ibonds(1, l)
          else
            if (nlive == 0) exit
            live(1) = live(nlive)
            nlive = nlive - 1
          end if
        end if
      end do
    end if
!
!   Check that atoms are not counted twice
!
    do j = 1, ninbit
      do k = j + 1, ninbit
        if (inres(j) == inres(k)) then
          inres(j) = 0
        end if
      end do
    end do
!
!  Put hydrogen atoms at the end of the list
!
    k = ninbit
    do j = 1, ninbit
      l = inres(j)
      if (l /= 0) then
        if (nat(l) == 1) then
          k = k + 1
          inres(k) = l
          inres(j) = 0
        end if
      end if
    end do
    ninbit = k
    if (iatom /= 0) then
      new = new + 1
      lused(new) = iatom
      iopt(iatom) = .true.
    end if
    do i2 = 1, ninbit
      j2 = inres(i2)
      if (j2 /= 0 .and. j2 /= iatom) then
        new = new + 1
        lused(new) = j2
      end if
    end do
!
!  If more than one hydrogen is attached to an atom, then the order of the hydrogen atoms is incorrect.
!  To correct this, the sequence of hydrogen atoms is reversed.  This is necessary in order to prevent
!  hydrogen atoms being swapped around by repeated calls to RESEQ.
!
    do j = 1, new
      if (nat(lused(j)) == 1) exit
    end do
    if (j > -new) return
!
! Start of hydrogen atoms = j
!
    k = j
    do
      do
        k = k + 1
        if (k == new) exit
        if (ibonds(1,lused(k)) /= ibonds(1,lused(j))) then
          k = k - 1
          exit
        end if
      end do
      if (k - j > 0) then
!
!  Two or three hydrogen atoms attached to the same atom
!  Swapping the atoms around is unusually simple!
!
        i2 = lused(j)
        lused(j) = lused(k)
        lused(k) = i2
      else
        k = k + 1
      end if
      if (k == new) exit
      j = k
    end do
    return
end subroutine moiety
integer function nheavy(icc)
  use common_arrays_C, only : nat, nbonds, ibonds
  implicit none
  integer, intent (in) :: icc
  integer :: i, n
  n = 0
  do i = 1, nbonds(icc)
    if (nat(ibonds(i,icc)) > 1) n = n + 1
  end do
  nheavy = n
  return
end function nheavy
subroutine identify_hexose(ninres, inres, nam, name)
  use common_arrays_C, only : nat, nbonds, ibonds, coord
  use funcon_C, only : pi
  implicit none
  integer, intent (in) :: ninres, inres(ninres)
  character, intent (inout) :: nam*3, name*60
!
!  Identify hexose by looking at the five chiral sites, C1 - C5
!  Assume that hexose exists as hemi-acetal
!
  double precision :: torsion
  integer :: i, j, k, l, m, backbone(6), C1, Cn, Cm, Cp, On, Hn, chiral(6)
  integer, external :: nheavy
  character :: hexose_aldose(8)*9, hexose_ketose(4)*9
  logical :: aldose
  data hexose_aldose /"Allose   ", "Altrose  ", "Glucose  ", "Mannose  ", &
    "Gulose   ", "Idose    ", "Galactose", "Talose   "/
   data hexose_ketose /"Psicose  ", "Fructose ", "Sorbose ", "Tagatose "/
   C1 = 0
   l = 0
   m = 0
!
! First, check that all carbon atoms have four ligands
!
    do i = 1, ninres
      if (nat(inres(i)) == 6) then
        if (nbonds(inres(i)) /= 4) then
          if (nbonds(inres(i)) < 3) return
          j = inres(i)
          k = ibonds(4,j)
          if (k == 0) return
          if (nat(k) /= 8) return
        end if
      end if
    end do
    aldose = .true.
!
!  Locate C1
!
    do i = 1, ninres
      if (nat(inres(i)) == 6) then
        C1 = inres(i)
        k = 0
        do j = 1, 4
          if (nat(ibonds(j,C1)) == 8) k = k + 1
        end do
        if (k == 2) exit
      end if
    end do
    do i = 1, nbonds(C1)
      if (nat(ibonds(i,C1)) == 6) then
        j = ibonds(i,C1)
        if (nheavy(j) == 2) then
          l = 0
          do k = 1, nbonds(j)
            if (nat(ibonds(k,j)) == 8) l = l + 1
          end do
          if (l == 1) then
            C1 = j
            aldose = .false.
            exit
          end if
        end if
      end if
    end do
    backbone(1) = C1
    Cn = C1
    Cm = Cn
!
!  Now locate C2 - C6
!
    do i = 2, 6
      do j = 1, 4
        l = ibonds(j,Cn)
        if (l == 0) return
        if (nat(l) == 6 .and. l /= Cm) exit
      end do
      backbone(i) = l
      Cm = backbone(i - 1)
      Cn = l
    end do
    chiral = 0
!
!  Work out chirality of C1
!
    On = 0
    do i = 1, 4
      j = ibonds(i,C1)
       if (j == 0) return
      if (nat(j) == 8 .and. nheavy(j) == 2) then
!
!   Make sure that the atom the oxygen is attached to is not in the ring
!
        do l = 1, nbonds(j)
          m = ibonds(l,j)
           if (m == 0) return
           if (nat(m) > 1 .and. m /= C1) exit
        end do
        do l = 1, ninres
          if (inres(l) == m) exit
        end do
        if (l > ninres) On = j
      end if
      if (nat(j) == 8 .and. nheavy(j) == 1) On = j
      if (nat(j) == 1                     ) Hn = j
    end do
    if (On == 0) return!  Failed to find oxygen
    Cp = backbone(2)
    call dihed (coord, Hn, Cp, C1, On, torsion)
    if (torsion > pi) torsion = torsion - 2*pi
    if (torsion < 0) chiral(1) = 1
!
!  Work out chirality of C2 - C5
!
    do k = 2, 5
      Cn = backbone(k)
      Cp = backbone(k + 1)
      l = 0
      do i = 1, 4
        if (nat(ibonds(i,Cn)) == 8) l = l + 1
      end do
      do i = 1, 4
        j = ibonds(i,Cn)
        if (l == 1) then
          if (nat(j) == 8) On = j
          if (nat(j) == 1) Hn = j
        else
          if (nat(j) == 8 .and. nheavy(j) == 2) Hn = j
          if (nat(j) == 8 .and. nheavy(j) == 1) On = j
        end if
      end do
      call dihed (coord, Hn, Cn, Cp, On, torsion)
      if (torsion > pi) torsion = torsion - 2*pi
      if (torsion > 0) chiral(k) = 1
    end do
!
!  Determine the hexose
!
    if (aldose) then
      i = chiral(2) + 2*chiral(3) + 4*chiral(4) + 1
      name = hexose_aldose(i)
    else
      i = chiral(3) + 2*chiral(4) + 1
      name = hexose_ketose(i)
    end if
    if (chiral(5) == 0) then
      name = "D-"//trim(name)
    else
      name = "L-"//trim(name)
    end if
    if (aldose) then
      if (chiral(1) == 0) then
        name = "alpha-"//trim(name)
      else
        name = "beta-"//trim(name)
      end if
    else
      if (chiral(2) == 0) then
        name = "alpha-"//trim(name)
      else
        name = "beta-"//trim(name)
      end if
    end if
    select case(name)
      case("alpha-D-Allose")
        nam = "ALO"
      case("alpha-D-Altrose")
        nam = "ALT"!
      case("alpha-D-Glucose")
        nam = "GLC"
      case("alpha-D-Mannose")
        nam = "MAN"!
      case("alpha-D-Gulose")
        nam = "GUL"!
      case("alpha-D-Idose")
        nam = "IDO"!
      case("alpha-D-Galactose")
        nam = "GAL"!
      case("alpha-D-Talose")
        nam = "TAL"!
      case("beta-D-Allose")
        nam = "BAL"
      case("beta-D-Altrose")
        nam = "BAT"!
      case("beta-D-Glucose")
        nam = "BGC"
      case("beta-D-Mannose")
        nam = "BMA"!
      case("beta-D-Gulose")
        nam = "BGU"!
      case("beta-D-Idose")
        nam = "BID"!
      case("beta-D-Galactose")
        nam = "BGA"!
      case("beta-D-Talose")
        nam = "BTA"!
      case("alpha-D-Psicose")
        nam = "ADP"!
      case("alpha-D-Fructose")
        nam = "ADF"!
      case("alpha-D-Sorbose")
        nam = "ADS"!
      case("alpha-D-Tagatose")
        nam = "ADT"!
    end select
  end subroutine identify_hexose
  subroutine inc_res(ires, start_res, nfrag)
    implicit none
    integer, intent (inout) :: ires, nfrag
    integer, intent (in) :: start_res(*)
!
    if (start_res(max(1,nfrag)) /= -200) then
      ires = start_res(nfrag)
      nfrag = nfrag + 1
    end if
    ires = ires + 1
    return
  end subroutine inc_res
