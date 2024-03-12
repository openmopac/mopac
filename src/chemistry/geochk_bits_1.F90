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

subroutine extvdw_for_MOZYME (radius, refvdw)
    use molkst_C, only: numat, keywrd, line
    use chanel_C, only: iw
    use common_arrays_C, only: nat
    use atomradii_C, only: is_metal
    use elemts_C, only : cap_elemnt
    implicit none
    double precision, dimension (numat), intent (out) :: radius
    double precision, dimension (107), intent (in) :: refvdw
    integer :: i, j, k
    double precision, dimension (107) :: vdw
    character (len=80) :: txt_rad
    logical :: paren
    double precision, external :: reada
!
!  Modify Van der Waals' radii of various atoms.  The format of
! the keyword is: VDW(:chemical symbol=n.nn:chemical symbol=n.nn ...)
! e.g. VDWM(:H=1.0:Cl=1.7)
!
    txt_rad = " "
    j = -1
    i = index(keywrd," METAL")
    if (i /= 0) j = index(keywrd(i:),") ") + i
    if (i /= 0 .and. keywrd(i + 6:i + 6) == " " .or. index(keywrd, " ADD-H") /= 0) then
      is_metal = .false.
      do i = 1, 102
        select case (i)
        case (3, 4, 11, 12, 19:30, 37:48, 55:80)
          is_metal(i) = .true.
        end select
      end do
    else if (i > 0) then
      i = index(keywrd(i:j), "(") + i
      txt_rad = keywrd(i:j - 1)
      do k = 1, 83
        if (.not. is_metal(k)) is_metal(k) = (index(txt_rad,cap_elemnt(k)) /= 0)
      end do
    end if
    i = index(keywrd," VDWM(")
    if (i == 0) then
      line = " "
    else
      if (keywrd(i + 6:i + 6) /= ";" .and. keywrd(i + 6:i + 6) /= ":") keywrd(i + 6:) = ";"//keywrd(i + 6:)
      j = index(keywrd(i + 6:),")") + i + 6
      do k = i + 5, j
         if (keywrd(k:k) == ":") keywrd(k:k) = ";"
         if (keywrd(k:k) == ",") keywrd(k:k) = ";"
      end do
      line = keywrd(i + 5:j - 1)
    end if
    vdw(:107) = refvdw(:107)
    if (line /= " ") then
      do i = 1, 107
        j = 2
        if (cap_elemnt(i)(2:2) == " ") j = 1
        k = index(line, ";"//cap_elemnt(i)(:j)//"=")
        if (k > 0) then
          vdw(i) = reada(line, k)
          is_metal(i) = .false.
        end if
      end do
    end if
!
!  Verify that all radii that will be used are, in fact, set correctly
!
    do i = 1, numat
      if (nat(i) > 98) cycle
      if (vdw(nat(i)) > 900.d0) then
        write (line,*) "MISSING VAN DER WAALS RADIUS FOR " // cap_elemnt (nat(i))
        call mopend(trim(line))
        j = 2
        if (cap_elemnt(nat(i)) (2:2) == " ") j = 1
        write (iw, "(2x,3a)") "To correct this, add keyword 'VDWM(", &
       & cap_elemnt (nat(i)) (1:j), "=n.nn)'"
        return
      end if
    end do
!
! Flag which elements are to be treated as metals by giving them negative
! atomic radii
    do i = 1, 102
      if (is_metal(i)) vdw(i) = -1.0d0
    end do
    do i = 1, numat
      radius(i) = vdw(nat(i))
    end do
!
!  Set VDW radius for individual atoms to -1
!
    if (txt_rad /= " ") then
      paren = .false.
      i = 0
      do k = 1, len_trim(txt_rad)
        i = i + 1
        if (paren) then
          paren = (txt_rad(i:i) /= ")")
        else
           paren = (txt_rad(i:i) == "(")
          if (ichar(txt_rad(i:i)) - ichar("0") <= 9 .and. ichar(txt_rad(i:i)) - ichar("0") > 0) then
            j = nint(reada(txt_rad(i:), 1))
            i = i + int(log10(j*1.0001))
            radius(j) = -1.d0
            do j = i, len_trim(txt_rad)
              if (ichar(txt_rad(i:i)) - ichar("0") > 9 .or. &
              ichar(txt_rad(i:i)) - ichar("0") < 0) exit
              txt_rad(i:i) = " "
            end do
          end if
        end if
      end do
    end if
    return
  end subroutine extvdw_for_MOZYME
    subroutine rectify_sequence()
    use common_arrays_C, only : geo, coord, txtatm, txtatm1, nat, labels, coorda
    use molkst_C, only: numat, pdb_label
    implicit none
!
!  Search atom sequence for the following structure:
!  A break in a chain: same chain letter, but the residue sequence increases suddenly.
!  After the increase, atoms of a residue in the same chain that are out-of-sequence are found.
!
! If this situation occurs, move the out-of-sequence atoms to their correct place.
!
! This should only be done if RESEQ or similar keyword is present.
!
    integer :: i_lower, res2, res1, res3, res, ii, jj, num_min, num_max
    character :: chain1*1, chain2*1, chain*1, used(numat)*1
    double precision, external :: reada
    if (pdb_label) return
    txtatm1(:numat) = txtatm(:numat)
    coorda(:,:numat) = coord(:,:numat)
    used(:numat) = txtatm1(:numat)(22:22)
    used(1) = "a"
    labels(:numat) = nat(:numat)
    geo(:,:numat) = coord(:,:numat)
    res1 = nint(reada(txtatm(1), 23))
    chain1 = "a"
!
!  ii keeps track of the new numbering sequence
!
    ii = 1
    do i_lower = 2, numat
      chain2 = used(i_lower)
      if (chain2 == "a") cycle           ! Atom already used
!
!   res1: Residue number of previous residue
!   res2: Residue number of current residue
!
      res2 = nint(reada(txtatm1(i_lower), 23))
!
! Look for a fault at the end of the first bit of chain
! This will occur in atoms in the domain res1 to res1 + 1
!
      do res3 = res1, res1 + 1
        num_min = max(i_lower - 50, 1)
        num_max = min(i_lower + 50, numat)
        do
          do jj = num_min, num_max
            res = nint(reada(txtatm1(jj), 23))
            chain = used(jj)
            if (res == res3 .and. chain == chain1 .and. used(jj) /= "a") exit
          end do
          if (jj > num_max) exit
          ii = ii + 1
          coord(:,ii) = geo(:,jj)
          txtatm(ii) = txtatm1(jj)
          nat(ii) = labels(jj)
          used(jj) = "a"
        end do
      end do
!
! Look for a fault at the start of the next bit of chain
! This will occur in atoms in the domain res2 - 1 to res2
!
      do res3 = res2 - 1, res2
        num_min = max(i_lower - 50, 1)
        num_max = min(i_lower + 50, numat)
        do
          do jj = num_min, num_max
            res = nint(reada(txtatm1(jj), 23))
            chain = used(jj)
            if (res == res3 .and. chain == chain1 .and. used(jj) /= "a") exit
          end do
          if (jj > num_max) exit
          ii = ii + 1
          coord(:,ii) = geo(:,jj)
          txtatm(ii) = txtatm1(jj)
          nat(ii) = labels(jj)
          used(jj) = "a"
        end do
      end do
      chain1 = chain2
      res1 = res2
    end do
    if (ii < numat) then
      do jj = 1, numat
        if (used(jj) /= "a") then
          ii = ii + 1
          coord(:,ii) = geo(:,jj)
          txtatm(ii) = txtatm1(jj)
          nat(ii) = labels(jj)
          used(jj) = "a"
        end if
      end do
    end if
    return
  end subroutine rectify_sequence
  subroutine write_sequence
    use common_arrays_C, only : txtatm, nat
!
    use MOZYME_C, only : ions, allres, tyr, allr, tyres, maxres
!
    use chanel_C, only: iw
!
    use molkst_C, only: numat, line, job_no, maxtxt
    implicit none
    integer :: i, nfrag, jj, ii, ires, kl, ku, irold, l, j, k, charge, i_job_no = -50, &
      ifrag, i_alt
    integer, allocatable :: res_no(:)
    character, allocatable :: res_charge(:)*1
    logical, save :: prt, use_alt, use_alt_2
    character :: chain, ch*2, ilet*1, jj_let*1
    double precision, external :: reada
    save :: i_job_no
    prt = (i_job_no /= job_no)
    if (prt) then
      i_job_no = job_no
    else
     return
    end if
!
!  Write out all residue names, use information in txtatm to define sequence
!
      ii = 1
      ifrag = 0
      if (maxtxt == 27) allocate(res_no(numat), res_charge(numat))
      use_alt_2 = .false.
      do nfrag = 1, 100
        if (ii > numat) exit
        if (txtatm(ii) == " ") exit
        chain = txtatm(ii)(22:22)
        ires = nint(reada(txtatm(ii), 23))
        allres = " "
        allres(ires) = txtatm(ii)(18:20)
        ilet = txtatm(ii)(27:27)
        charge = ions(ii)
        irold = ires
        j = 0
        i_alt = 0
        use_alt = .false.
        do ii = ii + 1, numat
          if (nat(ii) /= 1 .and. txtatm(ii)(22:22) /= chain) exit !  chain letter has changed.
          jj = nint(reada(txtatm(ii),23))
          jj_let = txtatm(ii)(27:27)
          if (nat(ii - 1) /= 1) j = ii - 1
          charge = charge + ions(ii)
          if (ires == jj .and. ilet == jj_let) cycle              !  Residue is same
!
! At this point, atom "ii" marks the start of the next residue
! and the charge on the current residue is "charge"
!
          charge = charge - ions(ii)
          if (txtatm(ii)(14:14) == "H") cycle!  Ignore hydrogen atoms
          if (ires + 1 /= jj .and. (ires /= jj .and. ilet == jj_let)) exit           !  Residue is not contiguous
!
!  Take residue name from the previous atom (the atom label has just changed)
!
!  Charge is the net charge on the previous residue, so:
!
          if (charge == 1) then
            allres(ires)(4:4) = "+"
          else if (charge == -1) then
            allres(ires)(4:4) = "-"
          end if
          ires = ires + 1
          allres(ires)(1:3) = txtatm(ii)(18:20)
          if (maxtxt == 27) then
            i_alt = i_alt + 1
            res_no(i_alt) = ii - 1
            res_charge(i_alt) = allres(ires)(4:4)
            if (.not. use_alt) use_alt = (ilet /= jj_let)
          end if
          charge = 0
          if (ires /= jj) ires = ires + 1
          ilet = jj_let
          if (ires > maxres) then
            write(line,'(a,i5)')" Maximum residue number allowed:", maxres
            call mopend(trim(line))
            return
          end if
        end do
        if (ii > numat) j = numat
        if (j /= 0) then
          if (txtatm(j)(18:20) /= "   " .and. allres(ires) == " ") allres(ires) = txtatm(j)(18:20)
        end if
        if (charge == 1) then
          allres(ires)(4:4) = "+"
        else if (charge == -1) then
          allres(ires)(4:4) = "-"
        end if
        do j = irold, ires
          do k = 1, 20
            if (allres(j) == tyres(k)) exit
          end do
          if (k < 21) exit
        end do
        if (j > ires .or. (j == ires .and. j == 0)) cycle
        do i = ires, 1, -1
          if (allres(i) /= "HOH" .and. allres(i) /= "SO4") exit
        end do
        ires = i
        ifrag = ifrag + 1
        if (ifrag == 1) then
          write (iw, "(/,16X,A,/)") "RESIDUE SEQUENCE IN PROTEIN Chain: "//chain
        else
          write (iw, "(/,16X,A,I2, a/)") "RESIDUE SEQUENCE IN PROTEIN FRAGMENT:", ifrag, " Chain: "//chain
        end if
        if (use_alt) then
          do i = 1, i_alt, 10
            write(iw,'(3x,10a)')(txtatm(res_no(j))(23:27)//": "//txtatm(res_no(j))(18:20)//res_charge(j)//" ", &
            j = i, min(i + 9, i_alt))
          end do
          use_alt_2 = .true.
          cycle
        end if
        if (irold < 1) then
          jj = 100
          i = mod(irold + jj, 10)
          if (i == 0) i = 10
          kl = irold
          ku = min(ires, kl - i + 10)
          line = " "
          l = 6*(i - 1) + 1
          j = ((kl + jj)/10)*10 - jj
          if (i == 10) j = j - 10
          write(iw,"(8x,10i6)")((k - 10), k = 1, 10)
          write(iw,"(5x,i4,2x,a,10(a4,'  '))")j, line(:l), allres(kl:ku)
          do
            kl = ku + 1
            if (kl > min(0, ires)) exit
            j = j + 10
            ku = min(ires, ku + 10)
            write(iw,"(5x,i4,3x,10(a4,'  '))")j, allres(kl:ku)
         end do
         write(iw,'(a)')" "
        end if
        kl = max(irold, 1)
        i = mod(kl, 10)
        if (i == 0) i = 10
        ku = min(ires, kl - i + 10)
        line = " "
        l = 6*(i - 1) + 1
        j = (kl/10)*10
        if (i == 10) j = j - 10
        write(iw,"(8x,10i6)")(k, k = 1, 10)
        write(iw,"(5x,i4,2x,a,10(a4,'  '))")j, line(:l), allres(kl:ku)
        do
          kl = ku + 1
          if (kl > ires) exit
          j = j + 10
          ku = min(ires, ku + 10)
          write(iw,"(5x,i4,3x,10(a4,'  '))")j, allres(kl:ku)
        end do
      end do
      if (use_alt_2) return
      ii = 1
      ifrag = 0
      do nfrag = 1, 100
        j = 0
        if (ii > numat) exit
        if (txtatm(ii) == " ") exit
        chain = txtatm(ii)(22:22)
        ires = nint(reada(txtatm(ii), 23))
        irold = ires
        do ii = ii + 1, numat
          if (nat(ii) /= 1 .and. txtatm(ii)(22:22) /= chain) exit !  chain letter has changed.
          jj = nint(reada(txtatm(ii),23))
          if (nat(ii - 1) /= 1) j = ii - 1
          if (ires == jj) cycle              !  Residue is same
          if (txtatm(ii)(14:14) == "H") cycle!  Ignore hydrogen atoms
          if (ires + 1 /= jj) exit           !  Residue is not contiguous
          allres(ires) = txtatm(ii - 1)(18:20)
          ires = ires + 1
        end do
        if (j /= 0) then
          if (txtatm(j)(18:20) /= "   ") allres(ires) = txtatm(j)(18:20)
        end if
        do j = irold, ires
          do k = 1, 20
            if (allres(j) == tyres(k)) exit
          end do
          if (k < 21) exit
        end do
        if (j > ires .or. (j == ires .and. j == 0)) cycle
        do i = ires, 1, -1
          if (allres(i) /= "HOH" .and. allres(i) /= "SO4") exit
        end do
        ires = i
        do i = irold, ires
          do k = 1, 20
            if (tyres(k) == allres(i)) exit
          end do
          if (k > 0 .and. k < 21) then
            allr(i) = tyr(k)
          else if (allres(i) /= " ") then
            allr(i) = "X"
          else
            allr(i) = " "
          end if
        end do
        ifrag = ifrag + 1
        if (ifrag == 1) then
          write (iw, "(/,16X,A,/)") "RESIDUE SEQUENCE IN PROTEIN Chain: "//chain
        else
          write (iw, "(/,16X,A,I2, a/)") "RESIDUE SEQUENCE IN PROTEIN FRAGMENT:", ifrag, " Chain: "//chain
        end if
        jj = 100
        i = mod(irold + jj, 10)
        if (irold < 0) i = i - 10
        if (i == 0) i = 10
        kl = irold
        if (i < 0) then
          ku = min(ires, kl - i + 40)
        else
          ku = min(ires, kl - i + 50)
        end if
        line = " "
        j = ((kl + jj)/10)*10 - jj
        if (i == 10) j = j - 10
        if (i == 1) then
          ch = "10"
        else
          ch(1:1) = " "
          if (i < 0) then
             ch(2:2) = char(9 + i + ichar("0"))
             if (ch(2:2) == "0") then
               ch(1:1) = "1"
               i = 12
             else
               i = 2 - i
             end if
          else
             ch(2:2) = char(11 - i + ichar("0"))
          end if
        end if
        if (i /= 0) write (iw, "(5x,i4,a,"//ch//"a,x,10(10a1,x))")j, line(:i + 1), (allr(i), i=kl, ku)
        do
          kl = ku + 1
          if (kl > ires) exit
          j = j + 50
          ku = min(ires, kl + 49)
          write (iw, "(5X,I4,2X,5(10A1,1X))") j, (allr(i), i=kl, ku)
        end do
      end do
      return
  end subroutine write_sequence
  subroutine fix_charges(ichrge)
  use molkst_C, only: refkey, keywrd, line
  implicit none
  integer, intent (in) :: ichrge
  integer :: i, j
!
!  Modify charges in keyword line so that the system will run using MOZYME.
!  Steps:
!     If "old" CHARGE=n keyword exists in the reference keyword line, delete it.
!     Write new CHARGE keyword, if charge is non-zero
!     Rewind the data set file.
!     Write out the new data set.
!     Rewind it again, so it is ready for reading.
!     Rewind the output file, so the user does not see the faulty output.
!
      line = trim(refkey(1))
      call upcase(line,len_trim(line))
      i = Index(line," CHARGE=")
      if (i /= 0) then
        j = Index(refkey(1)(i+2:)," ")
        refkey(1)(i + 1:) = refkey(1)(i + j + 2:)
      end if
!                 12345678901
      i = Index(refkey(1),"            ")
      if (ichrge /= 0) then
        if (ichrge > 99) then
          write(refkey(1)(i:i+11),'(" CHARGE=",i3)')ichrge
        else if (ichrge > 9) then
          write(refkey(1)(i:i+11),'(" CHARGE=",i2)')ichrge
        else if (ichrge > 0) then
          write(refkey(1)(i:i+11),'(" CHARGE=",i1)')ichrge
        else if (ichrge > -10) then
          write(refkey(1)(i:i+11),'(" CHARGE=",i2)')ichrge
        else
          write(refkey(1)(i:i+11),'(" CHARGE=",i3)')ichrge
        end if
      end if
       i = Index(keywrd," CHARGE=")
      if (i /= 0) then
        j = Index(keywrd(i+2:)," ")
        keywrd(i + 1:) = keywrd(i + j + 2:)
      end if
      i = Index(keywrd,"            ")
      if (ichrge /= 0) then
        if (ichrge > 99) then
          write(keywrd(i:i+11),'(" CHARGE=",i3)')ichrge
        else if (ichrge > 9) then
          write(keywrd(i:i+11),'(" CHARGE=",i2)')ichrge
        else if (ichrge > 0) then
          write(keywrd(i:i+11),'(" CHARGE=",i1)')ichrge
        else if (ichrge > -10) then
          write(keywrd(i:i+11),'(" CHARGE=",i2)')ichrge
        else
          write(keywrd(i:i+11),'(" CHARGE=",i3)')ichrge
        end if
      end if
end subroutine fix_charges
subroutine add_sp_H(i1, i, i2)
!
! Given three atoms, i1, i, and i2, put a hydrogen atom at the point
! of a triangle
!
  use molkst_C, only: natoms, maxatoms
  use common_arrays_C, only : geo, nat, txtatm
  use chanel_C, only: iw
  implicit none
  integer :: i1, i, i2
  logical :: first = .true.
  double precision :: scale
  natoms = natoms + 1
  if (natoms > maxatoms) then
    if (first) then
      write(iw,*)" Too many changes. Re-run using the data set generated by this job"
      first = .false.
    end if
    natoms = natoms - 1
    return
  end if
  geo(:,natoms) = 2.d0*geo(:,i1) - 2.d0*geo(:,i) + geo(:,i2)
  scale = 1.1d0/sqrt( (geo(1,i1) - geo(1,natoms))**2 + &
                      (geo(2,i1) - geo(2,natoms))**2 + &
                      (geo(3,i1) - geo(3,natoms))**2 )
  geo(:,natoms) = geo(:,i1) + scale*(geo(:,natoms) - geo(:,i1))
  nat(natoms) = 1
  txtatm(natoms) = " "
  return
end subroutine add_sp_H
subroutine add_sp2_H(i1, i, i2)
!
! Given three atoms, i1, i, and i2, put a hydrogen atom at the point
! of a triangle
!
  use molkst_C, only: natoms, maxatoms
  use common_arrays_C, only : geo, nat, txtatm
  use chanel_C, only: iw
  implicit none
  integer :: i1, i, i2
  logical :: first = .true.
  double precision :: scale
  natoms = natoms + 1
  if (natoms > maxatoms) then
    if (first) then
      write(iw,*)" Too many changes. Re-run using the data set generated by this job"
      first = .false.
    end if
    natoms = natoms - 1
    return
  end if
  geo(:,natoms) = 3.d0*geo(:,i) - geo(:,i1) - geo(:,i2)
  scale = 1.1d0/sqrt( (geo(1,i) - geo(1,natoms))**2 + &
                      (geo(2,i) - geo(2,natoms))**2 + &
                      (geo(3,i) - geo(3,natoms))**2 )
  geo(:,natoms) = geo(:,i) + scale*(geo(:,natoms) - geo(:,i))
  nat(natoms) = 1
  txtatm(natoms) = " "
end subroutine add_sp2_H
subroutine add_sp3_H(i1, i, i2, i3)
!
! Given four atoms, i1, i, and i2, put a hydrogen atom at the point
! of a tetrahedron
!
  use molkst_C, only: natoms, maxatoms
  use common_arrays_C, only : geo, nat, txtatm
  use chanel_C, only: iw
  implicit none
  integer :: i1, i, i2, i3
  logical :: first = .true.
  double precision :: scale
  natoms = natoms + 1
  if (natoms > maxatoms) then
    if (first) then
      write(iw,*)" Too many changes. Re-run using the data set generated by this job"
      first = .false.
    end if
    natoms = natoms - 1
    return
  end if
  geo(:,natoms) = 4.d0*geo(:,i) - geo(:,i1) - geo(:,i2) - geo(:,i3)
  scale = 1.1d0/sqrt( (geo(1,i) - geo(1,natoms))**2 + &
                      (geo(2,i) - geo(2,natoms))**2 + &
                      (geo(3,i) - geo(3,natoms))**2 )
  geo(:,natoms) = geo(:,i) + scale*(geo(:,natoms) - geo(:,i))
  nat(natoms) = 1
  txtatm(natoms) = " "
end subroutine add_sp3_H

subroutine compare_sequence(n_new)
!
!  Compare calculated residue with residue name from data-set
!  and print any differences
!
  use common_arrays_C, only : txtatm, txtatm1, nat, chains, breaks
  use molkst_C, only: numat, line, maxtxt, keywrd
  use MOZYME_C, only : tyr, tyres
  use chanel_C, only: iw, log, ilog
  implicit none
  integer :: n_new
!
! Local
!
  integer :: i_atom, i, j, i_delta = 0, new_res, old_res, previous = -200, &
    mbreaks, loop
  character :: old*3, new*3
  double precision, external :: reada
  logical :: first = .true.
  character :: num
  first = .true.
  mbreaks = 1
  i_atom = 0
  do loop = 1, numat
!
!  Find first non-hydrogen atom
!
    do
      i_atom = i_atom + 1
      if (i_atom > numat) exit
       if (i_atom == breaks(mbreaks)) mbreaks = mbreaks + 1
      if (nat(i_atom) /= 1) exit
    end do
    if (i_atom > numat) exit
    new_res = nint(reada(txtatm(i_atom),23))
    new = txtatm (i_atom)(18:20)
    if (maxtxt == 14) then
      old_res = nint(reada(txtatm1(i_atom),12))
      old = txtatm1(i_atom)(8:10)
    else
      old_res = nint(reada(txtatm1(i_atom),23))
      old = txtatm1(i_atom)(18:20)
    end if
    if (new_res - i_delta /= old_res) then
      i_delta = new_res - old_res
    end if
    if (new_res /= previous) then
      if (old /= new .and. index(keywrd, "CHANGED ALREADY PRINTED") == 0) then
        if (old /= "HET" .or. new /= "HOH") then
          if (first) then
            write(iw,'(/16x,a)')"Residue names that have changed"
            write(iw,'(7x,a,/)')"Original residue name   Calculated residue name"
            if (log) then
              write(ilog,'(/16x,a)')"Residue names that have changed"
              write(ilog,'(7x,a,/)')"Original residue name   Calculated residue name"
            end if
            first = .false.
            line="XENO=("
          end if
          do i = 1, 20
            if (old == tyres(i)) exit
          end do
          num = tyr(i)
          if (i == 21) num = " "
          write(iw,'(i14,a, 2x,a3,3x, 11x,i5,a,2x,a3)')old_res, " "//chains(mbreaks), old, &
            new_res, " "//chains(mbreaks), new
          if (log) write(ilog,'(i14,a, 2x,a3,3x, a1, 11x,i5,a,2x,a3)')old_res, " "//chains(mbreaks), &
            old, num, new_res, " "//chains(mbreaks), new
          if (i < 21) then
            j = len_trim(line) + 1
            if (j > 200) exit
            if (old_res > -1) then
              num = char(ichar("1") + int(log10(old_res + 0.05)))
              write(line(j:),'(a1,i'//num//',a1,a1,a1)')chains(mbreaks), old_res, "=", tyr(i),","
            else
              num = char(ichar("1") + int(log10(-old_res + 0.05)) + 1)
              write(line(j:),'(a1,i'//num//',a1,a1,a1)')chains(mbreaks), old_res, "=", tyr(i),","
            end if
          end if
        end if
      end if
      previous = new_res
    end if
  end do
  call l_control("CHANGED ALREADY PRINTED", len("CHANGED ALREADY PRINTED"), 1)
  if (first .and. n_new == 0) then
  else
    j = len_trim(line)
    if (j > 6 .and. n_new == 0) then
      line(j:j) = ")"
      write(line,'(a)')"(Use the XENO keyword to re-define unrecognized residues.)"
      write(iw,'(/2x,a)')trim(line)
      if (log) write(ilog,'(/a)')trim(line)
    end if
  end if
  return
end subroutine compare_sequence
  subroutine find_salt_bridges(in_cat, in_ani, n_cat, n_ani)
!
!  Find Salt Bridges:  Identify all ionizable sites that can be used for making salt bridges
!
    use common_arrays_C, only : nat, txtatm, nbonds, ibonds
    use molkst_C, only : numat, line, keywrd, moperr, maxtxt
    use chanel_C, only: iw
    use MOZYME_C, only : tyres
    use atomradii_C, only: is_metal
    implicit none
    integer, intent(in) :: n_cat, n_ani, in_cat(n_cat), in_ani(n_ani)
    integer :: i, j, k, kk, l, m, n, n_cations, n_anions, n_C, n_H, n_O, C, pairs(2,200), n_pairs, &
      salt_bridges(2,200), n_salt, i_length, ii, jj, nh_cat, nh_ani, nmetals
    double precision :: r, cutoff, Rab(200), R_min, R_sorted(200)
    character :: bits*10, num
    logical :: het
    double precision, external :: distance, reada
    logical, external :: near_a_metal
    integer, allocatable :: cations(:), anions(:), metals(:)
    allocate (cations(numat), anions(numat), metals(numat))
    het = .false.
    C = 0
    ii = 0
    jj = 0
    nmetals = 0
     do i = 1, numat
      do j = 1, 102
        if (is_metal(nat(i))) exit
      end do
      if (j < 103) then
        nmetals = nmetals + 1
        metals(nmetals) = i
      end if
    end do
    if (n_cat == 0) then
      i = index(keywrd,"SALT=")
      if (i /= 0) then
        cutoff = reada(keywrd, i + 4)
      else
        cutoff = 4.d0
      end if
      nh_cat = 0
      nh_ani = 0
    else
      cutoff = 8.d0
      nh_cat = 1
      nh_ani = -1
    end if
!
!  Locate all Arg guanidine "C" atoms and all -C(R)H-NH2 "N" atoms
!
    n_cations = 0
    do i = 1, numat
      if (nat(i) == 6) then
        if (nbonds(i) == 3) then
          if (nat(ibonds(1,i)) == 7 .and. nat(ibonds(2,i)) == 7 .and. nat(ibonds(3,i)) == 7) then
!
!  Make sure it has the correct charge
!
            m = 0
            do j = 1, 3
              m = m + nbonds(ibonds(j,i))
            end do
            if (m == 8 + nh_cat) then
              if (n_cat > 0) then
                do n = 1, n_cat
                  if (txtatm(in_cat(n))(18:maxtxt) == txtatm(i)(18:maxtxt)) then
                    exit
                  end if
                end do
              else
                n = 0
              end if
              if (n <= n_cat .and. txtatm(i)(18:20) /= "UNK" &
                .and. txtatm(i)(22:22) >= "A" .and. txtatm(i)(22:22) <= "Z") then
                n_cations = n_cations + 1
                cations(n_cations) = i
              end if
            end if
          end if
        end if
      end if
    end do
    do i = 1, numat
      if (nat(i) == 7) then
        if (nbonds(i) == 3 + nh_cat) then
!
!  Test for -NH2, e.g. N-terminus and lysine
!
          n_C = 0
          n_H = 0
          do j = 1, 3 + nh_cat
            if (nat(ibonds(j,i)) == 6  .and. txtatm(ibonds(j,i))(18:20) /= "UNK" &
                .and. txtatm(ibonds(j,i))(22:22) >= "A" .and. txtatm(ibonds(j,i))(22:22) <= "Z") then
              n_C = n_C + 1
              C = ibonds(j,i)
            end if
            if (nat(ibonds(j,i)) == 1) n_H = n_H + 1
          end do
!
!  Make sure it has the correct charge
!
          if (n_C == 1 .and. n_H == 2 + nh_cat) then
            do j = 1, nbonds(C)
              if (nat(ibonds(j,C)) == 8) exit
              if (nat(ibonds(j,C)) == 7 .and. ibonds(j,C) /= i) exit
            end do
            if (j <= nbonds(C) .and. nbonds(C) == 3) then
!
!  See if it's a -CO-NH3
!
               j = 1
               do k = 1, nbonds(i)
                 if (nat(ibonds(k,i)) == 1) j = j + 1
               end do
            end if
            if (j > nbonds(C)) then
              if (n_cat > 0) then
                do n = 1, n_cat
                  if (txtatm(in_cat(n))(18:maxtxt) == txtatm(i)(18:maxtxt)) then
                    exit
                  end if
                end do
              else
                n = 0
              end if
              if (n <= n_cat .and. txtatm(C)(18:20) /= "UNK" &
                .and. txtatm(C)(22:22) >= "A" .and. txtatm(C)(22:22) <= "Z") then
                n_cations = n_cations + 1
                cations(n_cations) = C
              end if
            end if
          end if
        else if (nbonds(i) == 2 + nh_cat .and. txtatm(i)(18:20) == "HIS") then
!
!  Check for Histidine
!
          do j = 1, nbonds(i)
            k = ibonds(j,i)
            if (nbonds(k) < 3) cycle
            do l = 1, nbonds(k)
              if (nat(ibonds(l,k)) == 7 .and. ibonds(l,k) /= i) exit
            end do
            if (l <= nbonds(k)) exit
          end do
          if (j > nbonds(i)) cycle
          if (n_cat > 0) then
            do n = 1, n_cat
              if (txtatm(in_cat(n))(18:maxtxt) == txtatm(i)(18:maxtxt)) then
                exit
              end if
            end do
          else
            n = 0
          end if
!
!  Avoid a "double count"
!
          do k = 1, n_cations
            if (txtatm(cations(k))(18:maxtxt) == txtatm(ibonds(j,i))(18:maxtxt)) exit
          end do
          if (n <= n_cat .and. k > n_cations) then
            if ( .not. near_a_metal(i, i, metals, nmetals)) then
              n_cations = n_cations + 1
              cations(n_cations) = ibonds(j,i)
            end if
          end if
        end if
      end if
    end do
!
!  Locate all -COOH "C" atoms
!
    n_anions = 0
    do i = 1, numat
      if (nat(i) == 6) then
      if (nbonds(i) == 3) then
        n_C = 0
        n_O = 0
        do j = 1,3
          if (nat(ibonds(j,i)) == 6) n_C = n_C + 1
          if (nat(ibonds(j,i)) == 8) n_O = n_O + 1
        end do
        if (n_C == 1 .and. n_O == 2) then
          m = 0
          ii = 0
          do j = 1, 3
            if (nat(ibonds(j,i)) == 8) then
              l = ibonds(j,i)
              ii = ii + nbonds(l)
              do k = 1, nbonds(l)
                if (nat(ibonds(k,l)) == 1) m = m + 1
              end do
            end if
          end do
          if (m == 1 + nh_ani .and. ii < 4) then
            if (n_ani > 0) then
              do n = 1, n_ani
                if (txtatm(in_ani(n))(18:maxtxt) == txtatm(i)(18:maxtxt)) then
                  exit
                end if
              end do
            else
              n = 0
            end if
            if (n <= n_ani .and. txtatm(i)(18:20) /= "UNK" &
                .and. txtatm(i)(22:22) >= "A" .and. txtatm(i)(22:22) <= "Z") then
                n_anions = n_anions + 1
                anions(n_anions) = i
              end if
            end if
          end if
        end if
      end if
    end do
!
!  Find interatomic distances and weight them by type
!
    n_pairs = 0
    do i = 1, n_cations
      do j = 1, n_anions
        r = distance(cations(i), anions(j))
        if (r < cutoff + 2.d0) then
          m = cations(i)
          n = anions(j)
          R_min = 100.d0
          do k = 1, nbonds(m)
            if (nat(ibonds(k,m)) /= 7) cycle
            do l = 1, nbonds(n)
              if (nat(ibonds(l,n)) /= 8) cycle
              r = distance(ibonds(k,m), ibonds(l,n))
              if (r < R_min) then
                R_min = r
                ii = ibonds(k,m)
                jj = ibonds(l,n)
!
!  Check that hetero-groups are excluded, i.e., by default, do not form salt-bridges between the protein and any substrate molecules.
!  This is because substrates might have multiple -COO groups.
!
                do kk = 1, 20
                  if (txtatm(ii)(18:20) == tyres(kk)) exit
                end do
                het = (kk > 20)
                if (.not. het) then
                  do kk = 1, 20
                    if (txtatm(jj)(18:20) == tyres(kk)) exit
                  end do
                  het = (kk > 20)
                end if
              end if
            end do
          end do
          if (R_min < cutoff .and. .not. het) then
            n_pairs = n_pairs + 1
            pairs(1,n_pairs) = ii
            pairs(2,n_pairs) = jj
            Rab(n_pairs) = R_min
          end if
        end if
      end do
    end do
!
!  Sort interatomic distances
!
    n_salt = 0
    do
      if (n_pairs == 0) exit
      r = 1.d10
      do i = 1, n_pairs
        if (r > Rab(i)) then
          r = Rab(i)
          j = i
        end if
      end do
      if (r > 40.d0) exit
      n_salt = n_salt + 1
      salt_bridges(:, n_salt) = pairs(:, j)
      R_sorted(n_salt) = r
      n_pairs = n_pairs - 1
      do k = j, n_pairs
        Rab(k) = Rab(k + 1)
        pairs(:, k) = pairs(:, k + 1)
      end do
    end do
    if (n_cat == 0) then
!
!  Eliminate duplicates
!
      do i = 1, n_salt
        if (salt_bridges(1,i) == 0) cycle
        do j = n_salt, i + 1, -1
          if (salt_bridges(1,j) == 0) cycle
          m = salt_bridges(1,i)
          n = salt_bridges(2,i)
          ii = salt_bridges(1,j)
          jj = salt_bridges(2,j)
          if (salt_bridges(1,j) == salt_bridges(1,i) .or. salt_bridges(2,j) == salt_bridges(2,i)) then
            salt_bridges(:,j) = 0
          else if (txtatm(m)(18:) == txtatm(jj)(18:)) then
            salt_bridges(:,j) = 0
          else
            do ii = 1,2
              m = salt_bridges(ii,i)
              n = salt_bridges(ii,j)
              do k = 1, nbonds(m)
                if (nat(ibonds(k,m)) /= 6) cycle
                do l = 1, nbonds(n)
                  if (nat(ibonds(l,n)) /= 6) cycle
                  if (ibonds(k,m) == ibonds(l,n))  salt_bridges(:,j) = 0
                end do
              end do
              if (salt_bridges(1,j) == 0) exit
            end do
          end if
        end do
      end do
    end if
!
!  Build new SITE keyword
!
    if (n_salt > 0) then
      num = char(ichar("4") + int(log10(cutoff + 0.04999d0)))
      write(iw,'(//19x, a, f'//num//'.1, a)')"Salt Bridges Found (Up to", cutoff," Angstroms)"
      write(iw,'(/4x,a,10x,a,28x,a,16x,a,/)')" No.","Cationic site","Anionic site","Dist. (Angstroms)"
     if (n_cat == 0) call update_txtatm(.true., .true.)
    else
       if (n_cat == 0) then
         call mopend("SALT option used in SITE command, but no Salt Bridges found")
         moperr = .false.
         return
       else
         write(iw,'(//19x,a)')"No Salt Bridges Found"
       end if
    end if
    line = " "
    j = 0
    do i = 1, n_salt
      if (salt_bridges(1,i) == 0) cycle
      k = salt_bridges(1,i)
      do m = 26, 23, -1
        if (txtatm(k)(m:m) == " ") exit
      end do
      m = max(23,m + 1)
      l = salt_bridges(2,i)
      do n = 26, 23, -1
        if (txtatm(l)(n:n) == " ") exit
      end do
      n = max(23,n + 1)
      i_length = len_trim(line)
      if (i_length > 1900 - len_trim(keywrd)) then
        call mopend("TOO MANY SALT BRIDGES FOUND")
        call mopend("(TO ADD MORE SALT BRIDGES, USE THE ""ARC"" OR ""PDB"" FILE FROM THIS JOB)")
        moperr = .false.
        exit
      end if
      write(line(i_length + 1:),'(2a)')","//txtatm(k)(22:22)//txtatm(k)(m:maxtxt)//"(+),", &
                                        txtatm(l)(22:22)//txtatm(l)(n:maxtxt)//"(-)"
      j = j + 1
      write(iw,'(i7,2x,a,a9,4x,a,a9, f10.2)') &
      j, "("//txtatm(k)(:maxtxt)//")", txtatm(k)(22:22)//txtatm(k)(m:maxtxt)//"(+)", &
      "("//txtatm(l)(:maxtxt)//")", txtatm(l)(22:22)//txtatm(l)(n:maxtxt)//"(-)", R_sorted(i)
      if (txtatm(k)(22:22) < "A" .or. txtatm(k)(22:22) > "Z") then
        write(line,'(a)')"Chain letter for the cation is not in the range 'A' to 'Z'"
        call mopend(trim(line))
        write(iw,'(10x,a)')"Chain letter = """//txtatm(k)(22:22)//""""
        return
      end if
      if (txtatm(l)(22:22) < "A" .or. txtatm(l)(22:22) > "Z") then
        write(line,'(a)')"Chain letter for the anion is not in the range 'A' to 'Z'"
        call mopend(trim(line))
        write(iw,'(10x,a)')"Chain letter = """//txtatm(l)(22:22)//""""
        return
      end if
    end do
    if (n_cat > 0) return
    line = line(2:)
!
!  Delete the word "SALT"
!
    k = index(keywrd, "(SALT") + index(keywrd, ",SALT") + 1
    i = index(keywrd,"(SALT=") + index(keywrd, ",SALT=")
    if (i > 0) then
      do j = i, len_trim(keywrd)
        if (keywrd(j:j) == "," .or. keywrd(j:j) == ")") exit
      end do
      if (keywrd(j:j) == ",") j = j + 1
      keywrd(i + 1:) = keywrd(j:)
      if (keywrd(i:i) == ",") keywrd(i:) = trim(keywrd(i + 1:))
    else
      i = index(keywrd, "(SALT") + index(keywrd, ",SALT") + 1
      keywrd(i:) = keywrd(i + 4:)
      if (keywrd(i:i) == ",") keywrd(i:) = trim(keywrd(i + 1:))
    end if
!
!   Check for duplicate sites in keywrd
!
    m = index(keywrd," SITE")
    m = index(keywrd(m:), "(") + m
    j = index(keywrd(m:),") ") + m
    if (keywrd(m:m) == ",") m = m + 1
    do
      k = index(keywrd(m:j), ",")
      if (k == 0) k = index(keywrd(m:j + 1), ") ")
      if (k == 0) exit
      k = k + m - 2
      bits = keywrd(m:k)
      if (k - m < 3) exit
      l = index(line, trim(bits))
      if (l > 0) then
!
!  Remove the comma
!
        if (line(l - 1:l - 1) == ",") then
          line = line(:l - 2)//line(l + k - m + 1:)
        else
          line = line(:l - 1)//line(l + k - m:)
        end if
      end if
      m = k + 2
      if (keywrd(m:m) == ",") m = m + 1
    end do
!
!  Insert the new keyword
!
    i = index(keywrd," SITE")
    if (keywrd(i + 7:i + 7) /= ")") then
!
! Keyword is NOT SITE=(SALT)
! it is either SITE=(text,SALT,text), SITE=(SALT,text), or SITE=(text,SALT),
! whichever it is, add a comma
!
      keywrd(i + 7:) = trim(line)//","//trim(keywrd(i + 7:))
    else
      keywrd(i + 7:) = trim(line)//trim(keywrd(i + 7:))
    end if
    return
  end subroutine find_salt_bridges
subroutine update_txtatm(output, sort)
!
!  If "output" is true. and keyword RESIDUES is not present, then use input text in TXTATM1
!
!  if "sort" is true, then use the geometry in coorda to work out which atom is which.
!
  use common_arrays_C, only : nat, txtatm, txtatm1, coorda, nbonds, ibonds, coord, &
    breaks
  USE molkst_C, ONLY: numat, numat_old, keywrd, pdb_label, line
  use MOZYME_C, only: tyres
  implicit none
  logical :: output, sort, pdb, update_chain
!
! Local
!
  logical :: L_all
  integer :: i, j, k, l, H_Z
  integer, allocatable :: n_H(:), nn_H(:)
  if (.not. pdb_label) return
!
!   Do the hydrogen atom labels need to be re-done?
!  "Yes," if atoms are added or removed or shuffled.
!
  if (index(keywrd," RESID") + index(keywrd," ADD-H") + &
    index(keywrd," SITE=") + index(keywrd," RESEQ") /= 0 .or. output) then
    H_Z = 1
  else
    H_Z = 0
  end if
  pdb = (index(keywrd, " PDBOUT") /= 0)
  L_all = (output .and. index(keywrd," RESID") == 0)
  if (L_all .and. .not. sort .and. numat == numat_old) then
!
!  Do nothing!
!
    do i = 1, numat
      if (txtatm1(i) == " ") then
        txtatm1(i) = txtatm(i)
      else
        txtatm(i) = txtatm1(i)
      end if
    end do
    goto 99
  end if
!
!  Add text to TXTATM to label hydrogen atoms and to add chain letter
!
  if (index(keywrd," ADD-H") + index(keywrd," SITE=") + index(keywrd," RESEQ") /= 0) &
    call set_up_dentate()
  if (index(keywrd, "0SCF") == 0 .and. index(keywrd, "SITE") == 0) call check_cvs(.true.)
  call check_H(i)
  update_chain = (index(keywrd, " CHAINS=(") == 0)
  do i = 1, numat
    if (sort) then
      if (nat(i) > H_Z) then
        if (L_all .or. update_chain) then
          do j = 1, numat_old
            if (abs(coord(1,i) - coorda(1,j)) > 0.01d0) cycle
            if (abs(coord(2,i) - coorda(2,j)) > 0.01d0) cycle
            if (abs(coord(3,i) - coorda(3,j)) > 0.01d0) cycle
            if (L_all) then
              if (txtatm1(j) /= " ") txtatm(i) = txtatm1(j)
            else
              if (txtatm1(j)(22:22) /= " ") txtatm(i)(22:22) = txtatm1(j)(22:22)
            end if
            exit
          end do
        end if
      else
        if (nbonds(i) > 0) txtatm(i) = txtatm(ibonds(1,i))
      end if
    else
      if (L_all) then
        if (nat(i) > 1) then
          if (txtatm1(i) /= " ") txtatm(i) = txtatm1(i)
        else if (nbonds(i) > 0) then
          txtatm(i) = txtatm(ibonds(1,i))
        end if
      else
        if (update_chain .and. txtatm1(i)(22:22) /= " ") txtatm(i)(22:22) = txtatm1(i)(22:22)
      end if
    end if
  end do
  99 continue
!
!  Number hydrogen atoms, if more than one on a heavy atom
!
  allocate(n_H(numat), nn_H(numat))
  n_H = 0
  nn_H = 0
  do i = 1, numat
    if (nat(i) > H_Z) then
      j = 0
      do k = 1, nbonds(i)
        if (nat(ibonds(k,i)) == 1) j = j + 1
      end do
      n_H(i) = j
    end if
  end do
  do i = 1, numat
    if (nat(i) <= H_Z .and. nbonds(i) > 0) then
      k = ibonds(1,i)
      txtatm(i) = txtatm(k)(:12)//" H"//txtatm(k)(15:)
!
!  If a hydrogen atom is attached to an unlabeled carbon atom, make the hydrogen atom a terminal hydrogen
!
      if (nat(k) == 6 .and. txtatm(i)(15:16) == " ") then
        do j = 1, 20
          if (txtatm(i)(18:20) == tyres(j)) exit
        end do
        if (j < 21) txtatm(i)(15:16) = "XT"
      end if
      if (n_H(k) == 1) then
        txtatm(i)(13:14) = " H"
      else
        nn_H(k) = nn_H(k) + 1
        txtatm(i)(13:14) = char(nn_H(k) + ichar("0"))//"H"
      end if
    end if
    if (numat == numat_old) then
      if (txtatm1(i) == " ") txtatm1(i) = txtatm(i)
    end if
  end do
!
!  Check for duplicate hydrogen atom labels
!
  do i = 1, numat
    if (nat(i) == 1) then
      l = 1
      do j = i + 1, numat
        if (nat(j) == 1) then
          if (txtatm(i)(13:) == txtatm(j)(13:)) then
            do k = 13, 16
              if (txtatm(j)(k:k) == " ") then
                l = l + 1
                txtatm(j)(k:k) = char(l + ichar("0"))
                if (numat == numat_old) then
                  if (txtatm1(j) == " ") txtatm1(j) = txtatm(j)
                end if
                exit
              end if
            end do
          end if
        end if
        if (l == 9) exit
      end do
      if (l > 1) then
        do k = 13, 16
            if (txtatm(i)(k:k) == " ") then
              txtatm(i)(k:k) = "1"
              if (numat == numat_old) then
                if (txtatm1(i) == " ") txtatm1(i) = txtatm(i)
              end if
              exit
            end if
        end do
      end if
    end if
  end do

!
!  Add atom numbering, using PDB format, i.e., a TER group has its own number
!  (TER groups are represented by BREAKS)
!
  j = 1
  do i = 1, numat
    write(line,'(a6,i5,a)')txtatm(i)(:6),i + j - 1,txtatm(i)(12:)
    txtatm(i) = trim(line)
    if (PDB .and. i == breaks(j)) j = j + 1
  end do
  return
  end subroutine update_txtatm
