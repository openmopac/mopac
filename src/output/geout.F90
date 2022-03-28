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

      subroutine geout(mode1)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use common_arrays_C, only : labels, na, nb, nc, geo, nat, loc, txtatm, &
      & p, atmass, all_comments, pibonds_txt, l_atom
      USE parameters_C, only : tore, ams
      USE molkst_C, ONLY: numat, natoms, ndep, keywrd, maxtxt, line, &
      & ncomments, moperr, nl_atoms, gui
      USE symmetry_C, ONLY: depmul, locpar, idepfn, locdep
      USE elemts_C, only : elemnt
      USE chanel_C, only : iw
      use maps_C, only : lpara1, latom1, lpara2, latom2, latom, lparam
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: mode1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mode, iprt, i, j, n, ia, ii, k, igui, store_maxtxt
      double precision, dimension(natoms) :: q2
      double precision :: degree, w, x, y, z
      logical :: cart, lxyz, isotopes, charge, store_moperr, nabc_pdb, sym_pdb
      character , dimension(3) :: q*2
      character :: flag1*2, flag0*2, flagn*2, blank*80, fmt1*4, fmt2*4, fmt3*4, num*1
!*********************************************************************
!
!   GEOUT PRINTS THE CURRENT GEOMETRY.  IT CAN BE CALLED ANY TIME,
!         FROM ANY POINT IN THE PROGRAM AND DOES NOT AFFECT ANYTHING.
!
!   mode1:   1 write geometry in normal MOPAC *.out format
!           -n write geometry in normal MOPAC *.arc format, but do not print keywords, title
!              symmetry data, etc.
!            n write geometry in normal MOPAC *.arc format
!
!*********************************************************************
      if (nl_atoms == 0 .and. (index(keywrd,' 0SCF') /= 0)) l_atom(:natoms) = .true.
      nabc_pdb = (index(keywrd, "CONTROL_NABC_in_PDB") /= 0)
      sym_pdb = (index(keywrd, "CONTROL_SYM_in_PDB") /= 0)
      mode = mode1
      igui = -10 ! Set to impossible value
      store_maxtxt = maxtxt
      if (index(keywrd, " NOTXT") /= 0) maxtxt = 0
      store_moperr = moperr
      moperr = .false.
      charge = (index(keywrd, " PRTCHAR") /= 0)
      lxyz = index(keywrd,' COORD') + index(keywrd,'VELO') /= 0
      if (index(keywrd,' 0SCF') /= 0 .and. lxyz) then
!
!  If 0SCF and coord and a polymer, then get rid of TV.
!
        natoms = numat
        lxyz = .FALSE.
      end if
      if (mode == 1) then
        flag1 = ' *'
        flag0 = '  '
        flagn = ' +'
        iprt = iw
      else
        if (gui) then
          flag1 = ' 1'
          flag0 = ' 0'
        else
          flag1 = '+1'
          flag0 = '+0'
        end if
        flagn = '-1'
        iprt = abs(mode)
      end if
      degree = 57.29577951308232D0
      if (lxyz) degree = 1.D0
      blank = ' '
      if (mode /= 1 .and. index(keywrd, " RESEQ") /= 0) then
!
! Atoms may have moved, so re-set atmass
!
         do i = 1, numat
           atmass(i) = ams(nat(i))
         end do
      end if
      if (mode /= 1 .and. nat(1) /= 0) then
        do i = 1, numat
          if (Abs(ams(nat(i)) - atmass(i)) > 1.d-6) exit
        end do
        isotopes = (i <= numat)
      else
        isotopes = .false.
      end if
      cart = .true.
      do i = 1, natoms
        if (na(i) > 0) cart = .false.
      end do
      if (mode > 1 .and. ncomments > 0) &
          write(iprt,"(a)")(all_comments(i)(:len_trim(all_comments(i))), i = 1, ncomments)
      fmt1  = "13.8"
      fmt2  = "13.7"
      fmt3  = "13.7"
      if (maxtxt /= 0)  maxtxt = maxtxt + 2
      if (cart) then
        x = 1.d-8
        y = 1.d-8
        z = 1.d-8
        do i = 1, natoms
          x = max(x, abs(geo(1,i)))
          y = max(y, abs(geo(2,i)))
          z = max(z, abs(geo(3,i)))
        end do
        fmt1 = "1"//char(ichar("2") + max(1,int(log10(x))))//".8"
        fmt2 = "1"//char(ichar("2") + max(1,int(log10(y))))//".8"
        fmt3 = "1"//char(ichar("2") + max(1,int(log10(z))))//".8"
        if (mode == 1) then
          if (maxtxt == 0) then
            i = 6
          else
            i = maxtxt/2 + 1
          end if
          write (iprt,'(4a)') &
             & "   ATOM "//blank(:maxtxt/2 + 1)//" CHEMICAL  "//blank(:i)//"  X               Y               Z"
          write (iprt,'(4a)') "  NUMBER "//blank(:maxtxt/2 + 1)//" SYMBOL" &
          &//blank(:i)//"(ANGSTROMS)     (ANGSTROMS)     (ANGSTROMS)"
          write (iprt,*)
        else if (mode > 0) then
          call wrttxt (iprt)
          if (moperr) return
        end if
      else if (mode == 1) then
        j = max(9,maxtxt + 2)
        if (maxtxt == 0) j = 8
        write (iprt,'(4a)') &
           & "  ATOM"//blank(:j/2)//"CHEMICAL "//blank(:(j + 1)/2)//" BOND LENGTH      BOND ANGLE     TWIST ANGLE "
        write (iprt,'(4a)') &
           & " NUMBER"//blank(:j/2)//"SYMBOL   "//blank(:(j + 1)/2)//"(ANGSTROMS)      (DEGREES)       (DEGREES) "
        write (iprt, "(A)") &
           & "   (I)       "//blank(:j)//"      NA:I           NB:NA:I    " // &
           & "   NC:NB:NA:I " // "      NA    NB    NC "
      else if (mode > 0) then
        call wrttxt (iprt)
      end if
      if (mode /= 1 .and. allocated(p) .and. index(keywrd, " 0SCF") == 0) then
        call chrge (p, q2)
        q2(:numat) = tore(nat(:numat)) - q2(:numat)
      else
        q2(:numat) = 0.D0
      end if
      n = 1
      ia = loc(1,1)
      ii = 0
      blank = " "
      num = char(ichar("2") + int(log10(natoms*1.01)))
      do i = 1, natoms
        do j = 1, 3
          q(j) = flag0
          if (ia /= i) cycle
          if (j /= loc(2,n)) cycle
          q(j) = flag1
          n = n + 1
          ia = loc(1,n)
        end do
        if (na(i) > 0) then
          w = geo(2,i)*degree
          x = geo(3,i)*degree
          if (i > 3) then
!
!  CONSTRAIN ANGLE TO DOMAIN 0 - 180 DEGREES
!
            w = w - aint(w/360.D0)*360.D0
            if (w < -1.d-6) w = w + 360.D0
            if (w > 180.000001D0) then
              x = x + 180.D0
              w = 360.D0 - w
            end if
!
!  CONSTRAIN DIHEDRAL TO DOMAIN -180 - 180 DEGREES
!
            x = x - aint(x/360.D0 + sign(0.5D0 - 1.D-9,x) - 1.D-9)*360.D0
          end if
        else
          w = geo(2,i)
          x = geo(3,i)
        end if
        if (latom == i) q(lparam) = flagn
        if ( .not. gui .and. latom1 == i) q(lpara1) = flagn
        if (latom2 == i) q(lpara2) = flagn
        blank = elemnt(labels(i))
        if (labels(i) /= 99 .and. labels(i) /= 107) ii = ii + 1
        if (labels(i) == 1) then  ! Check for Deuterium and tritium
          if (Abs(atmass(ii) - 2.014d0) < 1.d-3) blank(1:2) = " D"
          if (Abs(atmass(ii) - 3.016d0) < 1.d-3) blank(1:2) = " T"
        end if
        k = 4
        if (maxtxt > 0) then
          if (txtatm(i)(15:16) == "**") txtatm(i)(15:16) = "99"
          if (txtatm(i)(23:26) == "****") txtatm(i)(23:26) = "9999"
          line = txtatm(i)
          if (gui .and. txtatm(i)(3:) /= " ") then
            blank = trim(blank)//"("//txtatm(i)(3:)//'  '
            igui = -10
          else
            if (line /= " ") blank = trim(blank)//"("//line(:store_maxtxt)//')  '
            igui = 0
          end if
        end if
        k = 0
        if (mode /= 1 .and. ii > 0 .and. nat(1) /= 0) then
          if (Abs(ams(nat(ii)) - atmass(ii)) > 1.d-6) then
            j = int(log10(atmass(ii)))
            line(40:) = "(f8."//char(6 - j + ichar("0"))//")"
            if (blank(1:2) /= " D" .and. blank(1:2) /= " T") then
              write(line,line(40:))atmass(ii)
              do j = len_trim(line), 1, -1
                if (line(j:j) /= "0") exit
                line(j:j) = " "
              end do
              j = len_trim(line) + 1
              line(j:j) = "0"
              j = index(blank,"(")
              if (j > 0) then
                blank = blank(:j - 1)//trim(line)//blank(j:)
              else
                blank = blank(:2)//trim(line)
              end if
            end if
          end if
          j = max(4,maxtxt + 2)
          if (isotopes) j = j + 8
          k = max(0,8 - j)
        else
          if (mode /= 1) then
            j = max(4,maxtxt + 2)
          else
            j = max(9,maxtxt + 3)
          end if
        end if
        if (index(blank(:j),"(") /= 0) then
          if (index(blank(:j),")") == 0) then
            if (index(blank(j + 1:j + 1), ")") /= 0) j = j + 1
          end if
        end if
        if (labels(i) == 0) cycle
        if ( .not. l_atom(i)) cycle
        if (na(i) == igui .or. cart) then
          if (mode /= 1) then !  Print suitable for reading as a data-set
            if (labels(i)/=99 .and. labels(i)/=107) then

              if (maxtxt == 0 .and. .not. isotopes) j = 4
              if (charge) then
                write(line,'(a,f8.4)')blank(41:59 + k), q2(ii)
              else
                write(line,'(a,f8.4)')blank(41:59 + k)
              end if
              write (iprt, '(1X,A,F'//fmt1//',1X,A2,F'//fmt2//',1X,A2,F'//fmt3//',1X, A2, A)') &
                 blank(:j), geo(1,i), q(1), w, q(2), x, q(3),trim(line)
            else
              write (iprt, '(1X,A,F'//fmt1//',1X,A2,F'//fmt2//',1X,A2,F'//fmt3//',1X,A2,a)') &
              blank(:j), geo(1,i), q(1), w, q(2), x, q(3), " "
            end if
          else !  Print in output style
            if (maxtxt == 0) j = 9
            write (iprt, '(i6,6x,A,F'//fmt1//',1X,A2,F'//fmt2//',1X,A2,F'//fmt3//',1X,A2,a)')i  &
              , blank(:j), geo(1,i), q(1), w, q(2), x, q(3)
          end if
        else
          if (mode /= 1) then  !  Print suitable for reading as a data-set
            if (maxtxt == 0 .and. .not. isotopes) j = 4
            if (labels(i) /= 107) then
              if (nabc_pdb) then
                line = " "
                call atom_no_to_txt(na(i), line(2:))
                call atom_no_to_txt(nb(i), line(20:))
                call atom_no_to_txt(nc(i), line(38:))
                write(line(56:),'(a,3i'//num//')')"=", na(i), nb(i), nc(i)
              else
                write(line,'(3i6)')na(i), nb(i), nc(i)
              end if
              write(line(len_trim(line) + 1:),'(a,f8.4)') blank(41:41 + k)
              if (charge) write(line(len_trim(line) + 1:),'(f8.4)')q2(ii)
              write (iprt, &
                '(1X,A,F'//fmt1//',1X,A2,F'//fmt2//',1X,A2,F'//fmt3//',1X,A2,A)') &
                & blank(:j), geo(1,i), q(1), w, q(2), x, q(3), trim(line)
            else
              write (iprt, '(1X,A,F'//fmt1//',1X,A2,F'//fmt2//',1X,A2,F'//fmt3//',1X,A2,3I6)') &
                blank(:j), geo(1,i), q(1), w, q(2), x, q(3), na(i), nb(i), nc(i)
            end if
          else !  Print in output style
            if (maxtxt == 0) j = 9
            write (iprt, &
              '(I6,6X,A,F'//fmt1//',1X,A2,F'//fmt2//',1X,A2,F'//fmt3//',1X,A2,I6,2I6)') i, &
              blank(:j), geo(1,i), q(1), w, q(2), x, q(3), na(i), nb(i), nc(i)
          end if
        end if
      end do
      moperr = (moperr .or. store_moperr)
      maxtxt = store_maxtxt
      if (mode == 1) return
      write (iprt, *)
      if (ndep /= 0) then
!
!   OUTPUT SYMMETRY DATA.
!
        n = 1
        i = 1
        outer_loop: do
          j = i
          do
            if (j == ndep) exit outer_loop
             !
             !  Group together symmetry functions of the same type
             !  (same reference atom, same reference function, same multiplier,
             !   if function 18 or 19)
             !  (Maximum number of dependent atoms on a line: 9)
             !
            if (locpar(j) /= locpar(j+1) .or. idepfn(j) /= idepfn(j+1) .or. &
                 & j-i >= 9) exit
            if (idepfn(i) == 18 .or. idepfn(i) == 19) then
              if (Abs(depmul(n) - depmul(n+1)) > 1.d-10) exit
              n = n + 1
            end if
            j = j + 1
          end do
          if (idepfn(i) == 18 .or. idepfn(i) == 19) then
            write (iprt, "(I4,I3,F13.9,10I5)") locpar (i), idepfn (i), &
                 & depmul(n), (locdep(k), k=i, j)
            n = n + 1
          else
            if (sym_pdb) then
              line = " "
              call atom_no_to_txt(locpar(i), line(2:))
              write(line(len_trim(line) + 1:),'(i3)')idepfn(i)
              call atom_no_to_txt(locdep(i), line(len_trim(line) + 3:))
              do k = i + 1, j
                call atom_no_to_txt(locdep(k), line(len_trim(line) + 2:))
              end do
              write(iprt, '(a)')trim(line)
            else
              write (iprt, "(I4,I3,10I5)") locpar (i), idepfn (i), (locdep(k), k=i, j)
            end if
          end if
          i = j + 1
        end do outer_loop
        if (na(i) > 0 .and. (idepfn(i) == 19 .or. idepfn(i) == 18)) then
          write (iprt, "(I4,I3,F13.9,10I5)") locpar (i), idepfn (i), &
               & depmul(n), (locdep(k), k=i, j)
        else
          if (sym_pdb) then
            line = " "
            call atom_no_to_txt(locpar(i), line(2:))
            write(line(len_trim(line) + 1:),'(i3)')idepfn(i)
            call atom_no_to_txt(locdep(i), line(len_trim(line) + 3:))
            do k = i + 1, j
              call atom_no_to_txt(locdep(k), line(len_trim(line) + 2:))
            end do
            write(iprt, '(a)')trim(line)
          else
            write (iprt, "(I4,I3,10I5)") locpar (i), idepfn (i), (locdep(k), k=i, j)
          end if
        end if
        write (iprt,*)
      end if
      if (index(keywrd, " SETPI") /= 0 .and. allocated(pibonds_txt)) then
!
!  The user has supplied some pi bonds.  Write them out.
!
        do n = 1, 10
          if (pibonds_txt(n) == " ") exit
          if (n == 1) then
            write(iprt,"(a)")trim(pibonds_txt(n))//"  User-supplied pi bonds"
          else
            write(iprt,"(a)")trim(pibonds_txt(n))
          end if
        end do
      end if
      moperr = (moperr .or. store_moperr)
      end subroutine geout
