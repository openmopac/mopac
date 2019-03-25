subroutine write_trajectory(xyz, mode, charge, escf, ekin, time, xtot)
  use molkst_C, only : step_num, numat, jloop => itemp_1, line, nl_atoms, keywrd
  use chanel_C, only : ixyz, xyz_fn
  USE elemts_C, only : elemnt 
  use common_arrays_C, only : nat, l_atom
  implicit none
  double precision, optional :: escf, ekin, time, xtot
  double precision, intent (in) :: xyz(3,numat)
  double precision, optional :: charge(numat)
  integer, intent (in) :: mode
  integer :: icalcn = 0, i, k, j, io_stat, ipdb = 14, imodel = 0, npdb, npt
  character :: dummy_char*1, num1*1, num2*1, num3*1
  character, allocatable :: store_path(:,:)*100, store_hofs(:)*100  
  double precision :: factor
  save :: icalcn, ipdb, imodel, npt
  if (icalcn /= step_num) then
    do i = len_trim(xyz_fn), 1, -1
      if (xyz_fn(i:i) == "/" .or. xyz_fn(i:i) == "\") exit
    end do
    open(unit=ixyz, file=xyz_fn)
    if (index(keywrd, " PDBOUT") /= 0) &
      open(unit = ipdb, file = xyz_fn(:len_trim(xyz_fn) - 3)//"pdb")
    icalcn = step_num
  end if
  select case (mode)  
  case (1)!  write out a trajectory
!
!  Write out "xyz" file
!
    write(ixyz,"(i6,a)") nl_atoms," "
    num1 = char(ichar("1") + int(log10(jloop*1.01)))
    factor = abs(escf)
    num2 = char(max(ichar("0"), ichar("0") + min(9, int(log10(factor)) - 1)))
    num3 = char(max(ichar("0"), ichar("0") + min(9, int(log10(4.184d0*factor)) - 1)))     
    write(line,'(a, i'//num1//', a, f1'//num2//'.3, a, f1'//num3//'.3, a)')"Profile.", jloop, &
      " HEAT OF FORMATION =", escf, " KCAL =", escf*4.184d0, " KJ"
    k = index(line, "HEAT")
    npt = npt + 1
    write(ixyz,'(a,i4,a)')line(:8), npt," "//trim(line(k:))
    do i = 1, numat
      if (l_atom(i)) write(ixyz,"(3x,a2,3f15.5)")elemnt(nat(i)), (xyz(j,i),j=1,3)
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
    endif
    if (index(keywrd, " PDBOUT") /= 0) then
      imodel = imodel + 1
      write(ipdb,'(a,i7)')"MODEL",imodel 
      call pdbout(ipdb)
      write(ipdb,'(a)')"ENDMDL"
    end if
    if (charge(1) > -100.d0) return  
    case (2)!  Reverse the path
    close (ixyz)
    open(ixyz, file = xyz_fn)
    rewind(ixyz)
    allocate (store_path(nl_atoms, jloop), store_hofs(jloop))
!
!  Read in the path already generated
!
    do i = 1, 100000
      read(ixyz,*, iostat=io_stat)dummy_char
      if (dummy_char =="0" .and. i < 0) return  ! dummy use of dummy_char
      read(ixyz,'(a)', iostat=io_stat)line
      if (io_stat /= 0) exit
      store_hofs(i) = trim(line)
      read(ixyz,"(a)", iostat=io_stat) &
        (store_path(k,i), k = 1, nl_atoms)
    end do
    close (ixyz)
    open(ixyz, file = xyz_fn)
    rewind(ixyz)
!
! Write out the path already generated, in reverse
!
    jloop = i - 1
    npt = 0
    do i = jloop, 2, -1  ! Exclude the first point - it's common to both paths
      write(ixyz,"(i6,a)") nl_atoms," "
      k = index(store_hofs(i), "HEAT")
      npt = npt + 1
      write(ixyz,'(a,i4,1x,a)')store_hofs(i)(:8), npt,trim(store_hofs(i)(k:))
      do k = 1, nl_atoms
        write(ixyz,"(a)") trim(store_path(k,i))
      end do
    end do  
!
!  Do the same with the PDB file
!
    if (index(keywrd, " PDBOUT") /= 0) then
      deallocate (store_path)
      open(unit = ipdb, file = xyz_fn(:len_trim(xyz_fn) - 3)//"pdb")
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
      open(unit = ipdb, file = xyz_fn(:len_trim(xyz_fn) - 3)//"pdb")
      rewind(ipdb)
!
! Write out the path already generated, in reverse
!
      imodel = 0
      jloop = i - 1
      do i = jloop, 2, -1! Exclude the first point - it's common to both paths
        imodel = imodel + 1
        write(ipdb,"(a,i6)") "MODEL ", imodel
        do k = 1, npdb
          write(ipdb,"(a)") trim(store_path(k,i))
        end do
      end do          
    end if
    deallocate (store_path)
  end select   
end subroutine write_trajectory
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

