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

subroutine lewis (use_cvs)
   !***********************************************************************
   !
   !
   !  LEWIS works out what atom is bonded to what atom for the purpose of
   !  generating a Lewis structure.  It starts with the topography generated
   !  by "set_up_dentate" (this is the raw topography), then tidies up the
   !  topography using a series of filters - keywords such as CVB, and criteria
   !  such as a hydrogen atom must bond to exactly one other atom, and that must
   !  be a non-hydrogen atom.
   !
   ! On exit, NBONDS holds the number of bonds to each atom, and
   !          IBONDS contains the atom numbers, in order of increasing
   !                 atomic number.
   !
   ! The maximum number of bonded atoms per atom is 4.  More are allowed
   ! for checking purposes, but for a SCF calculation, 4 is the maximum.
   !
   !***********************************************************************
    use molkst_C, only: natoms, numat, maxtxt, keywrd, line, moperr, prt_topo
    use chanel_C, only: iw
    use elemts_C, only: elemnt
    use common_arrays_C, only: txtatm, nat, nbonds, ibonds
    implicit none
    logical, intent (in) :: use_cvs
!
    logical :: debug, let
    integer :: i, ibad, j, jbad, k, l, m, n
    integer, dimension (:), allocatable :: iz
    l = 0
!
! Most of the time, LET should be false, the commonest exceptions being 0SCF and RESEQ
!
    allocate (iz(natoms))
    let = (Index (keywrd, " 0SCF") + Index (keywrd, " RESEQ") + Index (keywrd, " SITE=") +  &
      Index (keywrd, " GEO-OK") + Index (keywrd, " ADD-H") + Index (keywrd, " LET") /= 0)
    debug = (Index (keywrd, " LEWIS") /= 0)
!
!  Work out raw connectivity.  This depends only on the topology.
!  Arrays nbonds and ibonds are filled here.
!
    call set_up_dentate()
    if (use_cvs) call check_cvs(.true.)
!
! Sanity check - Are all TXTATM filled
!
     do i = 1, numat
       do j = 1, nbonds(i)
         if (ibonds(j,i) == 0) then
            call mopend("A SEVERE ERROR OCCURRED WHILE ATTEMPTING TO RESEQUENCE THE ATOMS.")
            if (index(keywrd," RESEQ") == 0) &
              write(iw,'(10x, a)')"(Try adding keyword ""NORESEQ"" and re-submit.)"
            write(iw,'(10x, a)')"(Ignore any error messages after this line.)"
            return
         end if
       end do
     end do
!
!  Put limit on maximum number of bonds to an atom.
!
!  This is a practical issue.  In theory, it should be possible to make Lewis structures
!  for all systems, but in practice some are just too difficult.  For example, [Fe(ii)(NH3)6](2+)
!  and [Ni(ii)(NH3)6](2+).  In these systems, consider the metal as bonding to a nitrogen atom.
!  The metal obviously has a formal charge of -4 and the nitrogen atoms have a formal charge of +1,
!  to give a net charge on the complex of +2.  However, writing software for this has been too difficult,
!  so it is easier to re-word such complexes as [Fe(ii)(NH3)5](2+).(NH3)  This does not affect the results,
!  but is a useful expedient.
!
    do i = 1, numat
      select case (nat(i))
      case (26, 28, 46, 78)  ! Fe, Ni, Pd, Pt
        do
          if (nbonds(i) > 5) then
            call remove_bond(i)
          else
            exit
          end if
        end do
      case default
      end select
    end do
!
!  Check bonding of hydrogen atoms - check for bridges
!
    call check_H(ibad)
    if (ibad == 1) then
      write (iw,*)
    end if
    if (ibad /= 0 .and. .not. let) then
      write(iw,'(//10x,a)') "Add keyword ""LET"" to allow job to continue"
      call mopend ("A hydrogen atom is badly positioned")
    end if
    jbad = ibad
    ibad = 0
!
!  Special check for nitro groups
!
    do i = 1, numat
      if (nat(i) == 7) then
        k = 0
        l = 0
        do j = 1, nbonds(i)
          if (nat(ibonds(j,i)) == 8 .and. nbonds(ibonds(j,i)) == 1) then
            k = ibonds(j,i)
            exit
          end if
        end do
        do j = j + 1, nbonds(i)
          if (nat(ibonds(j,i)) == 8 .and. nbonds(ibonds(j,i)) == 1) then
            l = ibonds(j,i)
            exit
          end if
        end do
!
!  Found a system with N bonded to two oxygen atoms (k and l),
!  which are not bonded to anything else.
!
        if (k /= 0 .and. l /= 0) then
!
!  Put a bond between the two oxygen atoms
!
          nbonds(k) = 2
          nbonds(l) = 2
          ibonds(2,k) = l
          ibonds(2,l) = k
        end if
      end if
    end do
!
!  End of special check for nitro groups
!
!
!  Put atoms into order of atomic number
!
    do i = 1, numat
      if (nbonds(i) > 1) then
        n = nbonds(i)
        iz(:n) = ibonds(:n,i)
        do j = 1, n
          k = 0
          do m = 1, n
            if (iz(m) /= 0) then
              if (nat(iz(m)) > k) then
                l = m
                 k = nat(iz(m))
              end if
            end if
          end do
          ibonds(j,i) = iz(l)
          iz(l) = 0
        end do
      end if
    end do
    if (prt_topo .and. debug .and. jbad /= 0 .and. .not. let) then
      if (maxtxt == 0) then
        j = 1
      else
        j = maxtxt/2 + 2
      end if
      write (iw, "(/,A,/)") "   TOPOGRAPHY OF SYSTEM"
      line = " "
      write (iw,*) "  ATOM No. "//line(:j)//"  LABEL  "//line(:j)//"Atoms connected to this atom"
      if (j == 0) then
        do i = 1, numat
          write (iw, "(I7,9X,A,9I7)") i, elemnt (nat(i)) // "  ", (ibonds(j, i), j=1, nbonds(i))
        end do
      else
        if (maxtxt > 2) then
          do i = 1, numat
            write (iw, "(I7,9X,A,9I7)") i, elemnt (nat(i)) // " (" // txtatm(i) (:maxtxt) // ") ", (ibonds(j, i), j=1, nbonds(i))
          end do
        else
          do i = 1, numat
            write (iw, "(I7,9X,A,9I7)") i, elemnt (nat(i)), (ibonds(j, i), j=1, nbonds(i))
          end do
        end if
      end if
    end if
    if (jbad /= 0 .and. .not. let) call mopend ("Geometry is faulty")
    if (let) moperr = .false.
end subroutine lewis
!
!
!
subroutine remove_bond(i)
  use common_arrays_C, only: coord, nbonds, ibonds
!
!  Remove the longest bond from atom i
!
  implicit none
  integer, intent (in) :: i
  double precision :: rmax, sum
  integer :: j, k, l
  rmax = 0.d0
  l = 0
  do j = 1, nbonds(i)
    k = ibonds(j,i)
    sum = (coord(1, i)-coord(1, k)) ** 2 + &
          (coord(2, i)-coord(2, k)) ** 2 + &
          (coord(3, i)-coord(3, k)) ** 2
    if (sum > rmax) then
      rmax = sum
      l = k
    end if
  end do
  k = 0
  do j = 1, nbonds(i)
    if (ibonds(j,i) /= l) then
      k = k + 1
      ibonds(k,i) = ibonds(j,i)
    end if
  end do
  nbonds(i) = nbonds(i) - 1
  k = 0
  do j = 1, nbonds(l)
    if (ibonds(j,l) /= i) then
    k = k + 1
      ibonds(k,l) = ibonds(j,l)
    end if
  end do
  nbonds(l) = nbonds(l) - 1
  end subroutine remove_bond
  subroutine check_h(ibad)
    use molkst_C, only: natoms, numat, line, mozyme, keywrd, numcal, maxtxt, pdb_label
    use chanel_C, only: iw, ilog, log, log_fn
    use elemts_C, only: elemnt
    use atomradii_C, only: radius
    use common_arrays_C, only: nat, labels, nbonds, ibonds, txtatm
    implicit none
    integer, intent (out) :: ibad
    integer :: i, j, m, n, k, l, ii, heavy(10), icalcn = -1
    integer, allocatable :: im(:)
    double precision :: rmin, r
    logical :: opend, prt
    double precision, external :: distance
    character :: num1*1, num2*1, num3*1
    save :: icalcn
    if ( .not. mozyme) then
      ibad = 0
      return
    end if
    prt = (icalcn /= numcal)
    if (prt) icalcn = numcal
    allocate (im(natoms))
    inquire(unit=ilog, opened=opend)
    if (.not. opend) open(unit=ilog, form='FORMATTED', status='UNKNOWN', file=log_fn, position='asis')
    j = 0
    k = 0
    do i = 1, natoms
      if (labels(i) /= 99) then
        j = j + 1
        im(j) = i
      end if
    end do
!
!  Check bonding of hydrogen atoms - check for bridges
!
    ibad = 0
    do i = 1, numat
      if (nat(i) == 1) then
        if (nbonds(i) == 0) then
!
!  Find nearest heavy atom to hydrogen atom i
!
          rmin = 1.d12
          do j = 1, numat
            if (j == i .or. nat(j) == 1) cycle
            r = distance(i, j)
            if (r < rmin) then
              rmin = r
              k = j
            end if
          end do
          nbonds(i) = 1
          ibonds(1,i) = k
          nbonds(k) = nbonds(k) + 1
          ibonds(nbonds(k),k) = i
          if (rmin < max(1.3d0, 2.d0*radius(k))) cycle
          if (ibad < 20 .and. &
            index(keywrd, " 0SCF") + index(keywrd, " LET") + index(keywrd, " ADD-H") + index(keywrd, " RESEQ") + &
            index(keywrd, " SITE") == 0) then
            num1 = char(Int(log10(i     + 1.0)) + ichar("1"))
            num2 = char(Int(log10(rmin  + 0.05d0)) + ichar("5"))
            num3 = char(Int(log10(im(k) + 1.0)) + ichar("1"))
            j = 1
            if (elemnt (nat(k))(1:1) == " ") j = 2
            write (line, "(a,i"//num1//",a,f"//num2//".2,a, a,i"//num3//",a)") &
              "Hydrogen atom ", i, " is", rmin, " Angstroms from the nearest non-hydrogen atom (", &
              elemnt (nat(k))(j:2), im(k),")"
            call mopend(trim(line))
            write(iw,"(10x,a,i"//num1//",a,/)")"(Label for atom ", i, ": """//trim(txtatm(i))//""")"
            inquire(unit=ilog, opened=opend)
            if (log) write (ilog, "(a)")trim(line)
            ibad = ibad + 1
          end if
          cycle
        else if (nbonds(i) == 1) then
!
!  The normal case
!
          r = distance(i, ibonds(1,i))
        else
!
!  Find the heavy atom, k, nearest to hydrogen atom i
!
          rmin = 1.d12
          m = 0
          k = 0
          do j = 1, nbonds(i)
            if (nat(ibonds(j,i)) == 1) cycle
            m = m + 1
            heavy(m) = ibonds(j,i)
            r = distance(i, ibonds(j,i))
            if (r < rmin .and. r > 0.95d0) then
              rmin = r
              k = ibonds(j,i)
            end if
          end do
          if (k == 0) then
            do j = 1, nbonds(i)
              if (nat(ibonds(j,i)) == 1) cycle
              m = m + 1
              heavy(m) = ibonds(j,i)
              r = distance(i, ibonds(j,i))
              if (r < rmin) then
                rmin = r
                k = ibonds(j,i)
              end if
            end do
          end if
!
! Write out error message
!
          if (m == 2 .and. prt) then
            j = 1
            m = heavy(1)
            n = heavy(2)
            num1 = char(Int(log10(m     + 1.0)) + ichar("1"))
            num2 = char(Int(log10(n     + 1.0)) + ichar("1"))
            num3 = char(Int(log10(im(i) + 1.0)) + ichar("1"))
            if (pdb_label) then
              write (iw, "(/,a,i"//num3//",6a)") " Hydrogen atom ", im(i), &
              & " bonds to """, txtatm(m), """ and """, txtatm(n),""""
            if (log) write (ilog, "(a,i"//num3//",6a)") " Hydrogen atom ", im(i), &
              & " bonds to """, txtatm(m), """ and """, txtatm(n),""""
            else
              write (iw, "(/,a,i"//num3//",a,a2,i"//num1//",a,a2,i"//num2//")") " Hydrogen atom ", im(i), &
                & " has covalent bonds to ", elemnt (nat(m)), im(m), " and ", elemnt (nat(n)), im(n)
              if (log) write (ilog, "(a,i"//num3//",a,a2,i"//num1//",a,a2,i"//num2//")") " Hydrogen atom ", im(i), &
                & " has covalent bonds to ", elemnt (nat(m)), im(m), " and ", elemnt (nat(n)), im(n)
            end if
          else if (m > 2 .and. prt) then
            if (ibad < 20) then
              num3 = char(Int(log10(im(i) + 1.0)) + ichar("1"))
              write (iw, "(/,a,i"//num3//",a,a2,I5,a,a2,i5)") " Hydrogen atom ", im(i), " has more than two covalent bonds"
              if (txtatm(im(i)) /= " ") then
                j = im(i)
                write(iw, '(5x, a)')"Atom label: """//txtatm(j)(:maxtxt)//""""
                do n = 1, nbonds(j)
                  l = ibonds(n,j)
                  write(iw,'(5x, a, f8.3, a)')"The bond-length to """//txtatm(l)(:maxtxt)//""" is ",distance(l, im(i)), " Angstroms"
                end do
              end if
              if (log) write (ilog, "(a,i"//num3//",a,a2,I5,a,a2,i5)") &
                 " Hydrogen atom ", im(i), " has more than two covalent bonds"
            end if
          end if
          if (m > 1) then
!
!  Remove all bonds to hydrogen atom i except the shortest bond, to atom "k"
!
            do j = 1, nbonds(i)
              ii = ibonds(j,i)
              if (ii == k .or. nat(ii) == 1) cycle
              m = 0
              do l = 1, nbonds(ii)
                if (ibonds(l,ii) /= i) then
                  m = m + 1
                  ibonds(m,ii) = ibonds(l,ii)
                end if
              end do
              nbonds(ii) = m
            end do
            nbonds(i) = 1
            ibonds(1,i) = k
            num1 = char(Int(log10(k     + 1.0)) + ichar("1"))
            if (prt) then
              if (pdb_label) then
                write (iw, "(12x,3a)") "(Keeping bond to """, txtatm(k), """)"
              else
                write (iw, "(25x,a,a2,i"//num1//",a)") "(Keeping bond to ", elemnt (nat(k)), im(k), ")"
              end if
            end if
          end if
        end if
      end if
      if (ibad == 21) write (iw,*) " Remaining errors suppressed"
    end do
    return
  end subroutine check_h


  subroutine check_CVS(let)
  use molkst_C, only: numat, keywrd, line, numcal, moperr, nscf, maxtxt, natoms
  use common_arrays_C, only: txtatm, nat, coord, &
      nbonds, ibonds
  use elemts_C, only: elemnt
  use atomradii_C, only: radius, is_metal
  use chanel_C, only: iw
  implicit none
  logical, intent (in) :: let
  integer :: i, j, k, l, n, m, ii, ncvb, icalcn = -1, ik, cvb_atoms(60)
  logical :: lbond, error, PDB_input, OO_first = .true., made, first = .true.
  double precision :: r
  double precision, external :: reada
  double precision, allocatable :: store_coord(:,:)
  double precision, external :: distance
  character :: num*1, num_1*1, store*2000
  save :: store_coord, OO_first, first
  error = .false.
  PDB_input = .false.
!
! Most of the time, LET should be false, the commonest exceptions being 0SCF and RESEQ
!
    if (allocated (store_coord) .and. icalcn /= numcal) then
      deallocate (store_coord)
      icalcn = numcal
    end if
    if (.not. allocated (store_coord)) then
      i = size(coord)/3
      allocate (store_coord(3,i))
      store_coord = coord
    end if
!
!  Parse CVB(2:3,45:65)
!
     store = trim(keywrd)
     i = index(keywrd, " CVB")
     if (i > 0) then
       j = index(keywrd(i:), ") ") + i
       line = keywrd(i + 1: j - 1)
       ii = len_trim(line)
!
!    If everything is in quotation marks, use the CVB keyword
!
       k = 0
       do
         k = k + 1
         if (k >= ii) exit
         if (line(k:k) == '"') then
           do
             k = k + 1
             if (line(k:k) == '"') exit
           end do
         end if
         if (line(k:k) >= "0" .and. line(k:k) <= "9") exit
       end do
       PDB_input = (k > ii - 3)
    end if
    i = Index (keywrd, " CVB")
    k = 0
    if (i /= 0) then
      k = Index (keywrd(i:), ")") + i
      if (keywrd(k:k) /= " ") then
        write(iw,"(//10x,a)")" There must be a space after the closing parenthesis of the CVB keyword"
        do
          k = k + 1
          if (keywrd(k:k) == " ") exit
        end do
        write(iw,'(/10x,a)')" CVB keyword: """//keywrd(i + 1:k-1)//'"'
        if (index(keywrd, " 0SCF") /= 0) return
        call mopend ("ERROR IN CVB KEYWORD")
        return
      end if
    end if
    if (k == i .and. i /= 0) then
      k = Index (keywrd(i + 3:), " ") + i + 3
      write(iw,"(//10x,a)")" No closing parenthesis for CVB keyword"
      write(iw,'(10x,a)')" CVB keyword ="//'"'//keywrd(i:k)//'"'
      if (index(keywrd, " 0SCF") /= 0) return
      call mopend ("ERROR IN CVB KEYWORD")
      return
    end if
    if (i /= 0) then
      do n = 1, 50
        k = Index (keywrd(i:), ") ") + i
        j = index(keywrd(i:k),"""") + i
        if (j == i) exit
        if (keywrd(j: j + 3) == 'O:-O') then
!
!  Run O-O check: If two oxygen atoms are bonded together, and both oxygen atoms are
!  bonded to other atoms, then delete the O-O bond.
!
!
! Jump over the text "O:-O"
!
          i = j + 5
          lbond = .false.
          do j = 1, numat
            if (nat(j) == 8) then
              do k = 1, nbonds(j)
                l = ibonds(k,j)
                if (nat(l) == 8 .and. nbonds(l) + nbonds(j) > 3) then
!
!  Delete O-O bond j-l
!
                  do m = k + 1, nbonds(j)
                    ibonds(m - 1,j) = ibonds(m,j)
                  end do
                  nbonds(j) = nbonds(j) - 1
                  do m = 1, nbonds(l)
                    if (ibonds(m,l) == j) then
                      do ik = m + 1, nbonds(l)
                        ibonds(ik - 1,l) = ibonds(ik,l)
                      end do
                      nbonds(l) = nbonds(l) - 1
                    end if
                  end do
                  lbond = .true.
                  if (OO_first) &
                    write(iw,'(10x,a)')"O - O Bond between """//txtatm(j)//""" and """//txtatm(l)//""" deleted"
                end if
              end do
            end if
          end do
          if (.not. lbond) then
             if (OO_first .and. index(keywrd, " 0SCF") == 0) &
               write(iw,'(/10x,"Option ""O:-O"" present in CVB, but no O-O bonds found.  This option ignored, job continuing")')
          end if
          OO_first = .false.
          cycle
        end if
        call txt_to_atom_no(keywrd, j - 1, let)
        if (moperr) return
      end do
      k = Index (keywrd(i:), ")") + i
      ncvb = 0
      if (first) then
        write(iw,'(20x,a,/)')" Bonds made and broken by keyword CVB"
        write(iw,'(9x,a)')"Atom #          Label                  Atom #          Label"
      end if
      do
        made = .true.
        j = Nint (reada (keywrd(i:k), 1))
        ii = Index (keywrd(i:k), ":") + i
        if (ii == i) then
          if (index(keywrd, '("O:-O")') /= 0) exit
          write (iw,*)
          write (iw,*) " ERROR IN CVB KEYWORD - ':' EXPECTED BUT NOT FOUND"
          j = Index (keywrd, " CVB")
          write(iw,'(10x,a)')" CVB keyword ="//'"'//keywrd(j:k)//'"'
          call mopend ("ERROR IN CVB KEYWORD")
          return
        end if
        i = ii
        if (j > natoms .or. j == 0 .or. j < -natoms) then
          write (iw,*) " AT LEAST ONE ATOM DEFINED BY CVB IS FAULTY"
          write (iw,*) " THE FAULTY ATOM NUMBER IS", j
          j = Index (keywrd, " CVB")
          write(iw,'(10x,a)')" CVB keyword ="//'"'//keywrd(j:k)//'"'
          error = .true.
        end if
        l = Nint (reada (keywrd(i:k), 1))
        if (l > natoms .or. l == 0 .or. l < -natoms) then
          write (iw,*) " AT LEAST ONE ATOM DEFINED BY CVB IS FAULTY"
          write (iw,*) " THE FAULTY ATOM NUMBER IS", l
          m = Index (keywrd, " CVB")
          write(iw,'(10x,a)')" CVB keyword ="//'"'//keywrd(m:k)//'"'
          error = .true.
        end if
        if (error) goto 99
        n = abs(l)
        m = abs(j)
        if (allocated(store_coord) .and. .not. PDB_input) then
!
!   If an atom has moved, find its new location
!
          do ii = 1, numat
            if (abs(coord(1,ii) - store_coord(1,n)) < 1.d-1) then
               if (abs(coord(2,ii) - store_coord(2,n)) < 1.d-1) then
                 if (abs(coord(3,ii) - store_coord(3,n)) < 1.d-1) exit
               end if
            end if
          end do
          l = sign(ii,l)
          do ii = 1, numat
            if (abs(coord(1,ii) - store_coord(1,m)) < 1.d-1) then
               if (abs(coord(2,ii) - store_coord(2,m)) < 1.d-1) then
                 if (abs(coord(3,ii) - store_coord(3,m)) < 1.d-1) exit
               end if
            end if
          end do
          j = sign(ii,j)
        else
          j = sign(m,j)
          l = sign(n,l)
        end if
        do ii = 1, ncvb
          if (cvb_atoms(ii) == m) exit
        end do
        if (ii > ncvb) then
          ncvb = ncvb + 1
          cvb_atoms(ncvb) = m
        end if
        m = abs(l)
        do ii = 1, ncvb
          if (cvb_atoms(ii) == m) exit
        end do
        if (ii > ncvb) then
          ncvb = ncvb + 1
          cvb_atoms(ncvb) = m
        end if
        if (j > 0 .and. l > 0) then
          do ii = 1, nbonds(j)
            if (ibonds(ii, j) == l .and. .not. let .and. index(keywrd, " GEO-OK") == 0 .and. &
              index(keywrd, " 0SCF") == 0) then
              r = distance(j,l)
              write(iw,"(/,a,i5,a,i5,a, f6.3, a)")" The bond to be made by CVB between atoms",j," and", &
              l, " already exists, bond length =", r," Angstrom."
              write(iw,'(/30x,a,/)')"PDB label              Coordinates of the two atoms"
              write(iw,"(a,i5,a,3x,a, 3f12.6)")" Atom:", j, ": "//elemnt(nat(j)), "("//txtatm(j)//")", coord(:,j)
              write(iw,"(a,i5,a,3x,a, 3f12.6)")" Atom:", l, ": "//elemnt(nat(l)), "("//txtatm(l)//")", coord(:,l)
              write(iw,"(/10x,a,/)") "(If an extra bond needs to be added, use keyword 'SETPI')"
              error = .true.
            end if
          end do
          r = distance(j,l)
          if (r > 3.d0*(radius(j) + radius(l)) .and. .not. let .and. index(keywrd, " GEO-OK") == 0 .and. &
              index(keywrd, " 0SCF") == 0) then
            if (.not. is_metal(nat(j)) .and. .not. is_metal(nat(l))) then
              num = char(ichar("5") + int(log10(r)))
              write(iw,'(/,a,i5,a,i5,a,f'//num//'.2,a)')" The bond to be made by CVB between atoms", j, &
              " and", l, " would be", r, " Angstroms long.  This is unrealistic."
              write(iw,'(/30x,a,/)')"PDB label              Coordinates of the two atoms"
              write(iw,"(a,i5,a,3x,a, 3f12.6)")" Atom:", j, ": "//elemnt(nat(j)), "("//txtatm(j)//")", coord(:,j)
              write(iw,"(a,i5,a,3x,a, 3f12.6)")" Atom:", l, ": "//elemnt(nat(l)), "("//txtatm(l)//")", coord(:,l)
              error = .true.
            end if
          end if
!
!  Do NOT add a bond if that bond already exists
!
          do ii = 1, nbonds(j)
            if (ibonds(ii,j) == l) exit
          end do
          if (ii > nbonds(j)) then
            nbonds(j) = nbonds(j) + 1
            nbonds(l) = nbonds(l) + 1
            ibonds(nbonds(j), j) = l
            ibonds(nbonds(l), l) = j
          end if
        else
!
!  make sure that the faulty bond actually exists.
!
          j = Abs(j)
          l = Abs(l)
          if (j == 0 .or. l == 0) then
            write(iw,"(10x,a)")" Error detected while reading CVB"
            j = Index (keywrd, " CVB")
            write(iw,'(10x,a)')" CVB keyword ="//'"'//keywrd(j:k)//'"'
            error = .true.
          end if
          lbond = .false.
          do n = 1, nbonds(j)
            if (ibonds(n,j) == l) then
!
!  Faulty bond exists, therefore delete it.
!
              ik = l
              do m = n + 1,nbonds(j)
                ibonds(m - 1,j) = ibonds(m,j)
              end do
              ibonds(nbonds(j),j) = ik
              nbonds(j) = nbonds(j) - 1
              lbond = .true.
              exit
            end if
          end do
          do n = 1, nbonds(l)
            if (ibonds(n,l) == j) then
!
!  Faulty bond exists, therefore delete it.
!
              ik = j
              do m = n + 1,nbonds(l)
                ibonds(m - 1,l) = ibonds(m,l)
              end do
              ibonds(nbonds(l),l) = ik
              nbonds(l) = nbonds(l) - 1
              lbond = .true.
              exit
            end if
          end do
          made = .false.
          if (.not. lbond .and. .not. let) then
            if (radius(l) > 0.1d0 .and. radius(j) > 0.1d0) then
            if (index(keywrd, " GEO-OK") == 0 .and. index(keywrd, " 0SCF") == 0) then
              num = char(ichar("2") + int(log10(j*1.0001)))
              num_1 = char(ichar("2") + int(log10(l*1.0001)))

              write(iw,"(/,a,i"//num//",a,i"//num_1//",a)")" The bond to be deleted by CVB between atoms",j," and", &
                l, " did not exist.  Correct error and re-submit"
              write(iw,'(/30x,a,/)')"PDB label              Coordinates of the two atoms"
              write(iw,"(a,i5,a,3x,a, 3f12.6)")" Atom:", j, ": "//elemnt(nat(j)), "("//txtatm(j)//")", coord(:,j)
              write(iw,"(a,i5,a,3x,a, 3f12.6)")" Atom:", l, ": "//elemnt(nat(l)), "("//txtatm(l)//")", coord(:,l)
            else
              if (index(keywrd, " LOCATE-TS") == 0 .and. index(keywrd, " 0SCF") == 0) then
                write(iw,"(/,a,i5,a,i5,a,/,a)")" The bond to be deleted by CVB between atoms",j," and", &
                l, " did not exist,", " but, because GEO-OK is present, the job will continue."
              end if
            end if
            error = .true.
            end if
          end if
        end if
99      continue
        if (.not. error .and. first) then
          if (made) then
            write(iw,'(a8, i6,a,i8,a)')" Made",   j, " "//elemnt(nat(j))//"("//txtatm(j)(:maxtxt)//")", &
              l, " "//elemnt(nat(l))//"("//txtatm(l)(:maxtxt)//")"
          else
            write(iw,'(a8, i6,a,i8,a)')" Broken", j, " "//elemnt(nat(j))//"("//txtatm(j)(:maxtxt)//")", &
              l, " "//elemnt(nat(l))//"("//txtatm(l)(:maxtxt)//")"
          end if
        end if
        do j = i, k
          if (keywrd(j:j) == "," .or. keywrd(j:j) == ";") exit
        end do
        if (j <= k) then
          i = j
        else
          exit
        end if
      end do
      do i = 1, ncvb
        j = cvb_atoms(i)
        if (nbonds(j) == 0 .and. radius(j) > 0.1d0) then
          if (nat(j) == 1) then
!
!  A bond to hydrogen has been broken, so see if there is a different atom within range.
!
            do_k: do k = 1, numat
              do l = 1, ncvb
                if (cvb_atoms(l) == k) cycle do_k
              end do
              if (distance(j,k) < 1.2d0) then
                nbonds(j) = 1
                ibonds(1,j) = k
                nbonds(k) = nbonds(k) + 1
                ibonds(nbonds(k),k) = j
                exit
              end if
            end do do_k
            if (nbonds(j) == 1) cycle
          end if
          if (nscf < 3) then
            if (txtatm(j) /= " ") then
            num = char(ichar("1") + int(log10(j*1.01)))
            write(line,"(a, i"//num//",a)")"All bonds to atom ", j, ", label: """ &
            //txtatm(j)//""" have been deleted"
            else
              write(line,"(a,i5,a,3f10.5)") &
              " All bonds to atom "//elemnt(nat(j)), j, " have been deleted. Coords:",coord(:,j)
            end if
            call mopend(trim(line))
            error = .true.
          end if
        end if
      end do
    end if
    if (error .and. .not. let) then
      if (index(keywrd, " GEO-OK") + index(keywrd, " 0SCF") == 0  ) then
        write(iw,"(/10x,a)")"To continue add ""GEO-OK"" or correct reported faults"
        call mopend("Fault in CVB keyword")
        return
      end if
    end if
    if (let) moperr = .false.
    keywrd = store
    first = .false.
    return
  end subroutine check_CVS
  subroutine txt_to_atom_no(text, j_in, let)
!
!  Convert atom number from text to a number.  Text can be PDB or Jmol
!
!  On input, "text": Line of text containing PDB or Jmol
!            j_in: Start of text (location of the character after the leading '"')
!            let:  TRUE if the job is to continue even if a fault is found.
!
! On exit, "text": PDB or Jmol text has been replaced by numbers
!
  use molkst_C, only: natoms, pdb_label
  use common_arrays_C, only: txtatm
  use chanel_C, only: iw
  implicit none
  character :: text*(*)
  integer, intent (in) :: j_in
  logical, intent (in) :: let
  integer :: i, j, l, n, k, m, p, q
  character :: line*80, txt*30, num*1, txt2*30
  if (.not. pdb_label) then
    call mopend("Labeled atoms can only be used when atoms have labels")
    m = 0
    return
  end if
  j = j_in + 1
  l = index(text(j:),"""") + j
  call upcase(text, len_trim(text))
  line = text(j:l - 2)
  txt2 = trim(line)
  if (line(1:1) == "[") then
    k = len_trim(line)
    if (k > 16) then
      write(line,'(a)')"Atom defined by '"//trim(line)//"' is not in JSmol format"
      call mopend(trim(line))
      return
    end if
!
! Atom defined using Jmol format
!
    n = index(text(j:l - 2),".") + j
    k = index(text(j:l - 2),":") + j
    m = index(text(j:l - 2),"]") + j
    if (n == j .or. m == j) then
      write(line,'(a)')"Atom defined by '"//trim(line)//"' is not in JSmol format"
      call mopend(trim(line))
      if (n == j) call mopend("(The dot separator, ""."", is missing.)")
      if (m == j) call mopend("(The close square bracket, ""]"", is missing.)")
      return
    end if
    if (k == j) then
      line = text(n:l - 2)//text(j + 1:m - 2)//text(m:n - 2)
    else
      line = text(n:l - 2)//text(j + 1:m - 2)//text(k:k)//text(m:k - 2)
    end if
  end if
  m = 0
  do k = 1, len_trim(line)
    if (line(k:k) /= " ") then
      m = m + 1
      line(m:m) = line(k:k)
    end if
  end do
  line(m + 1:) = " "
  if (m > 12) then
     write(line,'(a)')"Atom defined by '"//trim(line)//"' is not in PDB format"
     call mopend(trim(line))
     return
  end if
  line(m + 1:m + 5) = " "
  q = 0
  do i = 1, natoms
    txt = txtatm(i)
    call upcase(txt, len_trim(txt))
    n = 0
    do k = 13, 27
      if (txt(k:k) /= " ") then
        n = n + 1
        txt(n:n) = txt(k:k)
      end if
    end do
    txt(n + 1:m ) = "?"
    n = max(n, m)
    do k = 1, n
      if (txt(k:k) /= line(k:k) .and. line(k:k) /= "*") exit
    end do
    if (k > n) then
      q = q + 1
      if (q == 2) then
          write(iw,'(/1x,a)')"An atom has been defined using a PDB label, "// &
            "but more than one atom in the data-set has this label:"
          write(iw,'(/10x,a)')"Atom number      PDB label"
          write(iw,'(10x, i5, 4x, a)')p, '"'//txtatm(p)//'"'
          write(iw,'(10x, i5, 4x, a)')i, '"'//txtatm(i)//'"'        
          call mopend("TWO ATOMS HAVE THE SAME PDB LABEL")
          return
      end if
      p = i
    end if
  end do
  if (q == 1) then
    i = p
    num = char(ichar("1") + int(log10(i*1.0001)))
    write(line,'(i'//num//')')i
    text(j - 1:) = trim(line)//text(l:)          
  end if
  if (q == 0) then
    k = len_trim(line)
    if (index(txt2,"[") /= 0) then
      write(line,'(a)')"Atom defined by JSmol label '"//trim(txt2)//"' was not found in the data set"
    else
      write(line,'(a)')"Atom defined by PDB label '"//trim(txt2)//"' was not found in the data set"
    end if
    if(.not. let) then
      call mopend(trim(line))
      write(iw,'(10x,a)') &
      "(Hint: Check that the atom-label in the data set matches the atom-label in the PDB file.)"
      write(iw,'(10x,a)')"(Hint: Atom-labels used in connectivity can only use atoms that have already been defined.)"
      return
    end if
  end if
  return
  end subroutine txt_to_atom_no
  subroutine atom_no_to_txt(atom_no, text)
!
! Convert an atom number to its PDB label
!
    use common_arrays_C, only: txtatm
    implicit none
    integer, intent (in) :: atom_no
    character, intent (out) :: text*(*)
!
!  Local
!
    character :: txt*18
    if (atom_no == 0) then
      text = " 0 "
      return
    end if
    txt = txtatm(atom_no)(12:)
    do
      if (txt(1:1) /= " ") exit
      txt = trim(txt(2:))
    end do
    text = '"'//trim(txt)//'"'
    return
  end subroutine atom_no_to_txt
