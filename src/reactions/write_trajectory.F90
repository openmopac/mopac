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

subroutine write_trajectory(xyz, escf, ekin, rc_dipo, time, xtot, l_dipole)
!
! Write a single point in a trajectory to a file of type <file>.xyz or <file>.pdb
!
!  If the type is <file>.xyz then the single point consists of:
!
!  A line containing the number of atoms in the system.
!  A line containing the text "Profile." followed by the model number,  some more text ending in the character "="
!    and a real datum, either the heat of formation or the dipole moment.
!    The geometry in Cartesian coordinates
!
!  Example of a single point in such a file (6 lines):
!     4
!Profile.   2 HEAT OF FORMATION =    -0.467 KCAL =    -1.955 KJ
!    N        0.00000        0.00000       -0.00092
!    H        0.97231        0.00000        0.00427
!    H       -0.48616        0.84205        0.00427
!    H       -0.48616       -0.84205        0.00427
!
!  If the type if <file>.pdb then the single point consists of:
!
!  A line containing the text "MODEL" followed by the model number.
!  A PDB file giving a header line, the geometry, and an "END" line
!  A A line containing the text "ENDMDL"
!
!  Example of a single point in such a file (8 lines)
!MODEL      2
!HEADER Heat of Formation =       2.677 Kcal/mol
!HETATM    1  N   NH3 A   1       0.000   0.000   0.000  1.00 -0.83      PROT N
!HETATM    2 1H   NH3 A   1       0.920   0.000   0.000  1.00  0.28      PROT H
!HETATM    3 2H   NH3 A   1      -0.460   0.797   0.000  1.00  0.28      PROT H
!HETATM    4 3H   NH3 A   1      -0.460  -0.797  -0.000  1.00  0.28      PROT H
!END
!ENDMDL
!
  use molkst_C, only : step_num, numat, jloop => itemp_1, line, keywrd, backslash, nvar, &
  ncomments
  use chanel_C, only : ixyz, xyz_fn
  use drc_C, only : georef
  USE elemts_C, only : elemnt
  use common_arrays_C, only : nat, l_atom
  implicit none
  double precision, intent (in) :: escf, ekin, time, xtot, rc_dipo
  double precision, intent (in) :: xyz(3,numat)
  logical, intent (in) :: l_dipole
  integer :: icalcn = 0, i, k, j, ipdb = 14, imodel = 0, npt = 0, ixyz1
  character :: num1*1, num2*1, num3*1
  double precision :: factor
  save :: icalcn, ipdb, imodel, npt, ixyz1
  if (icalcn /= step_num) then
    ncomments = -1
    if (nvar >= numat*3 .and. index(keywrd, " MINI ") == 0 ) l_atom(:numat) = .true.
    do i = len_trim(xyz_fn), 1, -1
      if (xyz_fn(i:i) == "/" .or. xyz_fn(i:i) == backslash) exit
    end do
    open(unit=ixyz, file=xyz_fn)
    if (l_dipole) then
      ixyz1 = ixyz + 1
      open(unit=ixyz1, file=xyz_fn(:len_trim(xyz_fn) - 4)//" for dipole.xyz")
    end if
    if (index(keywrd, " PDBOUT") /= 0) &
      open(unit = ipdb, file = xyz_fn(:len_trim(xyz_fn) - 3)//"pdb")
    icalcn = step_num
  end if
!
!  Write out "xyz" file
!
  write(ixyz,"(i6,a)") numat," "
  if (index(keywrd, " REVERSE") /= 0) then
    npt = npt - 1
    imodel = imodel -1
    call l_control("REVERSE", len("REVERSE"), -1)
  end if
  num1 = char(ichar("1") + int(log10(jloop*1.01)))
  factor = abs(escf)
  num2 = char(max(ichar("0"), ichar("0") + min(9, int(log10(factor)) - 1)))
  num3 = char(max(ichar("0"), ichar("0") + min(9, int(log10(4.184d0*factor)) - 1)))
  npt = npt + 1
  if (l_dipole) then
    write(ixyz1,"(i6,a)") numat," "
    write(line,'(a, i'//num1//', a, f1'//num2//'.4, a)')"Profile.", jloop, &
    " DIPOLE =", rc_dipo, " DEBYE"
     k = index(line, "DIPO")
     write(ixyz1,'(a,i4,a)')line(:8), npt," "//trim(line(k:))
  end if
  write(line,'(a, i'//num1//', a, f1'//num2//'.3, a, f1'//num3//'.3, a)')"Profile.", jloop, &
    " HEAT OF FORMATION =", escf, " KCAL "
  k = index(line, "HEAT")
  write(ixyz,'(a,i4,a)')line(:8), npt," "//trim(line(k:))
  k = 0
  do i = 1, numat
    if (l_atom(i)) then
      k = k + 1
      if (l_dipole) write(ixyz1,"(3x,a2,3f15.5)")elemnt(nat(i)), (xyz(j,k), j = 1, 3)
      write(ixyz,"(3x,a2,3f15.5)")elemnt(nat(i)), (xyz(j,k), j = 1,3)
    else
      if (l_dipole) write(ixyz1,"(3x,a2,3f15.5)")elemnt(nat(i)), (georef(j,i), j = 1, 3)
      write(ixyz,"(3x,a2,3f15.5)")elemnt(nat(i)), (georef(j,i), j = 1, 3)
    end if
  end do
  if (mod(jloop,10) == 0) then! write out every 10'th iteration.
    if (time > 1.d-6) then
      write (line, "(a7,i6,a8,f17.5,a8,f9.4,a7,f9.4,a7,f9.3)") " CYCLE:", jloop, &
        & "  Pot.E:", escf, "  Kin.E:", ekin, "  Move:", xtot,"  Time:", Min(time,99999.999d0)
      else
        write (line, "(a7,i6,a19,f17.5,a8,f9.4,a7,f9.4,a7,f9.3)") " CYCLE:", jloop, &
    & "  Potential energy:", escf, " Diff.:", ekin, "  Move:", xtot
    end if
    call to_screen(line)
  end if
  if (index(keywrd, " PDBOUT") /= 0) then
    imodel = imodel + 1
    write(line,'(F13.5)') escf
    do i = 1, 12
      if (line(i:i) /= " ") exit
    end do
    write(ipdb,'(a, i9, 2x, a)')"MODEL",imodel, trim(line(i:))
    call pdbout(ipdb)
    write(ipdb,'(a)')"ENDMDL"
  end if
end subroutine write_trajectory
subroutine reverse_trajectory(mode)
!
! Read in a trajectory that was created earlier in this run, and
! write it back out, in reverse order of points.
!
  use molkst_C, only : step_num, jloop => itemp_1, line, nl_atoms, keywrd, backslash
  use chanel_C, only : ixyz, xyz_fn
  implicit none
  integer, intent (in) :: mode
  integer :: icalcn = 0, i, k, io_stat, ipdb = 14, imodel = 0, npdb, npt, ixyzn
  character :: dummy_char*1, xyz_fnn*200
  character, allocatable :: store_path(:,:)*200, store_hofs(:)*200
  save :: icalcn, ipdb, imodel, npt
  ixyzn = 0
  if (icalcn /= step_num) then
    do i = len_trim(xyz_fn), 1, -1
      if (xyz_fn(i:i) == "/" .or. xyz_fn(i:i) == backslash) exit
    end do
    if (mode == 1) then
      ixyzn = ixyz
      xyz_fnn = trim(xyz_fn)
    else
      ixyzn = ixyz + 1
      xyz_fnn = xyz_fn(:len_trim(xyz_fn) - 4)//" for dipole.xyz"
    end if
    open(unit=ixyzn, file=xyz_fnn)
    if (index(keywrd, " PDBOUT") /= 0) &
      open(unit = ipdb, file = xyz_fnn(:len_trim(xyz_fnn) - 3)//"pdb")
    if (mode == 1) icalcn = step_num
  end if

  close (ixyzn)
  open(ixyzn, file = xyz_fnn)
  rewind(ixyzn)
  allocate (store_path(nl_atoms, jloop), store_hofs(jloop))
!
!  Read in the path already generated
!
  do i = 1, 100000
    read(ixyzn,*, iostat=io_stat)dummy_char
    if (dummy_char =="0" .and. i < 0) return  ! dummy use of dummy_char
    read(ixyzn,'(a)', iostat=io_stat)line
    if (io_stat /= 0) exit
    store_hofs(i) = trim(line)
    read(ixyzn,"(a)", iostat=io_stat) &
      (store_path(k,i), k = 1, nl_atoms)
  end do
  close (ixyzn)
  open(ixyzn, file = xyz_fnn)
  rewind(ixyzn)
!
! Write out the path already generated, in reverse
!
  jloop = i - 1
  npt = 0
  do i = jloop, 2, -1  ! Exclude the first point - it's common to both paths
    write(ixyzn,"(i6,a)") nl_atoms," "
    k = index(store_hofs(i), "HEAT") + index(line, "DIPO")
    npt = npt + 1
    write(ixyzn,'(a,i4,1x,a)')store_hofs(i)(:8), npt,trim(store_hofs(i)(k:))
    do k = 1, nl_atoms
      write(ixyzn,"(a)") trim(store_path(k,i))
    end do
  end do
!
!  Do the same with the PDB file
!
  if (index(keywrd, " PDBOUT") /= 0) then
    deallocate (store_path)
    open(unit = ipdb, file = xyz_fnn(:len_trim(xyz_fnn) - 3)//"pdb")
    rewind (ipdb)
!
! How many lines are there in each PDB file?
!
    do npdb = 0, 10000
      read(ipdb,'(a)') line
      if (line(:6) == "ENDMDL") exit
    end do
    allocate (store_path(npdb, jloop))
    rewind (ipdb)
    do i = 1, 100000
      read(ipdb,*, iostat=io_stat)dummy_char
      if (dummy_char =="0" .and. i < 0) return  ! dummy use of dummy_char
      read(ipdb,"(a)", iostat=io_stat) &
        (store_path(k,i), k = 1, npdb)
      if (io_stat /= 0) exit
    end do
    close (ipdb)
    open(unit = ipdb, file = xyz_fnn(:len_trim(xyz_fnn) - 3)//"pdb")
    rewind(ipdb)
!
! Write out the path already generated, in reverse
!
    imodel = 0
    jloop = i - 1
    do i = jloop, 2, -1! Exclude the first point - it's common to both paths
      imodel = imodel + 1
      line = store_hofs(i)(33:43)
      do k = 1, 12
        if (line(k:k) /= " ") exit
      end do
      write(ipdb,'(a, i9, 2x, a)')"MODEL",imodel, trim(line(k:))
      do k = 1, npdb
        write(ipdb,"(a)") trim(store_path(k,i))
      end do
    end do
  end if
  deallocate (store_path)
end subroutine reverse_trajectory
subroutine reverse_aux
!
!  Read in the current AUX file, then write out everything down to the start of the IRC,
!  then write out the IRC in reverse.
!
  use molkst_C, only : npoints => itemp_1, line, keywrd
  use chanel_C, only : output_fn
  implicit none
  integer :: hook = 50, i, j, n_head, n_lines
  character, allocatable :: start(:)*100, irc(:,:)*100
  if (index(keywrd, " AUX") == 0) return
  close (hook)
  open(unit=hook,file=output_fn(:len_trim(output_fn) - 4)//".aux")
  do i = 1, 10000
    read(hook,'(a)') line
    if (index(line, "#               IRC                #") /= 0) exit
  end do
  n_head = i + 2
  read(hook,'(a)') line
  if (line /= "dummy") read(hook,'(a)') line! a clumsy way to "fool" forcheck
  if (line /= "dummy") read(hook,'(a)') line
  if (line == "dummy") return
  do i = n_head + 1, n_head + 10000
    read(hook,'(a)') line
    if (index(line, "REF_POINT") /= 0) exit
  end do
  n_lines = i - n_head
!
!  The sizes of the data blocks are now known,
!  so rewind, and store everything.
!
  rewind (hook)
  allocate (start(n_head), irc(n_lines, npoints))
  read(hook,'(a)')start
  read(hook,'(a)')irc
!
! Now write out the AUX file again, in reverse
!
  rewind (hook)
  do i = 1, n_head
    write(hook,'(a)')trim(start(i))
  end do
  do i = npoints, 2, -1
    write(hook,'(a,i5.5)')" REF_POINT=",npoints - i + 1
    irc(2,i)(21:21) = "-"
    do j = 2, n_lines
      write(hook,'(a)')trim(irc(j,i))
    end do
  end do
  write(hook,'(a)')" ######################################"
  write(hook,'(a)')" #                                    #"
  write(hook,'(a)')" # Next point is the transition state #"
  write(hook,'(a)')" #                                    #"
  write(hook,'(a)')" ######################################"
  deallocate(start, irc)
end subroutine reverse_aux
