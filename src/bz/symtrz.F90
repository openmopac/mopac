! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

subroutine symtrz (sec_det, nvecs, coord)
!
!  symtrz uses symmetry operations supplied by the user to symmetrize the
!  secular determinant.  This is a complicated operation, because the effect
!  of a symmetry operation usually involves moving an atom out of its
!  primitive unti cell and into a different unit cell.
!
  use common_common, only : numat, id, iop, nfirst, nlast, per_atom, ncell, nijk, &
    title, allt1, mop, mr1, mr2, mr3, bcc, keywrd, iw_new, ixyz, nijk_CUC, line, &
    data_set_name, jobnam
  implicit none
  integer, intent(in) :: nvecs
  double precision, dimension(3, numat), intent(in) :: coord
  double precision, dimension(nvecs, nvecs*(mr1*mr2*mr3 + 1)), intent(inout) :: sec_det
! 
  logical :: debugs, debugf, opend
  integer :: i, i1, ii, il, iu, j, j1, jj, ix, &
       & jl, jop, ju, jx, k, k1, kk, sym_op, m1, m2, m3, mm, newcel, &
       & nx, ny, nz, alloc_stat
  integer :: new(numat), newcuc(numat), incell(3,numat), itrans(3,numat), &
    nxyz(-11:11, -11:11, -11:11)
  double precision :: sum, new_coord(3,numat), t1(9,9), t2(9,9), new_cuc(3,numat)
  double precision, dimension(:,:), allocatable :: new_sec_det
!
!    Put the unit cell indices into the 3 dimensional integer array NXYZ
!    and also store the indices in a linear integer array NIJK. 
!
!    Initially, unit cells are in the range (1:mr1), (1:mr2), and (1:mr3), but
!    for this operation, they must be as near to the CUC as possible, so
!    they are re-numbered as (0:mr1/2, -mr1/2:-1), etc. I.e., if mr1=9, 
!    then the order of unit cells would be (0,1,2,3,4,-4,-3,-2,-1).
!
! 
  incell = 0
  itrans = 0
  sum = 0.d0
!
! Make a list, nijk, of unit cell indices, nx, ny, nz, and a map
! to that list, nxyz
!
  ncell = 0
  do m3 = 1, mr3
    if (m3 <= (mr3+1)/2) then
      nz = m3 - 1
    else
      nz = m3 - mr3 - 1
    end if
    do m2 = 1, mr2
      if (m2 <= (mr2+1)/2) then
        ny = m2 - 1
      else
        ny = m2 - mr2 - 1
      end if
      do m1 = 1, mr1
        if (m1 <= (mr1+1)/2) then
          nx = m1 - 1
        else
          nx = m1 - mr1 - 1
        end if
        if (.not.bcc .or. Mod(nx + ny + nz, 2) == 0) then
          ncell = ncell + 1        
          if (nx <= -mr1/2) nx = nx + mr1
          if (ny <= -mr2/2) ny = ny + mr2
          if (nz <= -mr3/2) nz = nz + mr3
          if (nx >  mr1/2) nx = nx - mr1
          if (ny >  mr2/2) ny = ny - mr2 
          if (nz >  mr3/2) nz = nz - mr3
          nxyz(nx, ny, nz) = ncell
          nijk(1, ncell) = nx
          nijk(2, ncell) = ny
          nijk(3, ncell) = nz
        end if
      end do
    end do
  end do
  if (index(keywrd," NOSYM") /= 0) then
    do i = 1, numat
      mop(1, i, 1) = i
      mop(2:4, i, 1) = 0 
    end do
    iop = 1
    return
  end if 
!
!  Make a list of the unit cell that each atom in the Central Unit Cell is in
!  (Some atoms in the CUC can be in the adjacent unit cell.)
!
!  Coordinates are in crystallographic units.
!
  do i = 1, numat
    do j = 1, id
      k = Int (coord(j, i) + 40.1d0) - 40
      nijk_CUC(j,i) = k
    end do
  end do
  debugs = (Index (keywrd, " SYMTRZ") /= 0)
  debugf = (Index (keywrd, " ROTFOK") /= 0)
  if (debugs .or. debugf) write(*, "(/,2a)") "  (For debug information, please look in the output file created by this run)"   
  allocate(new_sec_det(nvecs, ncell*nvecs), stat = alloc_stat)
  if (alloc_stat /= 0) then
    write(iw_new,'(a)')"Failed to allocate array new_sec_det in symtrz"
    stop
  end if
  new_sec_det = 0.d0
!
!   Outer-most sym_op:  Over all operations
!
  jop = iop
  big_sym_op: do sym_op = 1, iop
    if (debugf) write(iw_new, "(/,2A)") " Operation: ", title(sym_op)
!
!   Work out where the atoms of the CUC go
!
    new_coord = coord
    call rotsec (t1, .true., new_coord, sym_op)
!
!  On exit from rotsec, t1 holds the unitary rotation matrix.
!
    do i = 1, 9
      do j = 1, 9
        allt1(i, j, sym_op) = t1(i, j)
      end do
    end do
!
!   Work out the behavior of the atoms in the CUC.
!
!   Atoms in the CUC might, in fact, not be in the CUC, i.e., they might not have indices = (0,0,0)
!   so keep track of the unit cell they are in.  For atom "i" this is nijk_CUC(:,i)

!   Clear newcuc so that if a mistake occurs, it will be obvious.
!
    do j = 1, numat
      newcuc(j) = 0
    end do
    do i = 1, numat
      do j = 1, id
        k = Int (new_coord(j, i) + 40.1d0) - 40 - nijk_CUC(j,i)
        new_coord(j, i) = new_coord(j, i) - k 
!
!  itrans holds the indices of the atoms of the CUC after they have been moved by operation sym_op
!  These indices are "relative" in the sense that the fractional coordinates are relative to those
!  of the atoms they came from in the CUC.  I.e., if an atom in the CUC has coordinates (0.1 -0.1 0.1)
!  then after the operation the indices of that atom will have (0,1,0) added to it.
!
        itrans(j, i) = k
      end do
      new_cuc = new_coord
      if (Mod (itrans(1, i) + itrans(2, i) + itrans(3, i), 2) /= 0 .and. bcc) then
        inquire(unit=iw_new, opened=opend) 
        if (.not. opend) open (unit = iw_new, file = trim(data_set_name)//".txt", &
        form = "FORMATTED", status = "UNKNOWN")
        write(iw_new, '(a,i3,a)') " Operation,", sym_op, " is not an allowed operation."
        write(iw_new, "(A,I3,3F12.6,3I3)") &
          & " New coordinates of CUC atom", i, (new_coord(j, i), j = 1, 3), &
          & (itrans(j, i), j = 1, 3)
        jop = jop - 15
        cycle big_sym_op
      end if
      sym_op1: do j = 1, numat
        do k = 1, 3
          if (Abs (new_coord(k, i)  - nijk_CUC(k,i) - coord(k, j) + nijk_CUC(k,j)) > 1d-4) cycle sym_op1
        end do
        new_coord(:, i) = new_coord(:, i) + nijk_CUC(:,j) - nijk_CUC(:,i)
        itrans(:, i) = itrans(:, i) - nijk_CUC(:,j) + nijk_CUC(:,i)
        newcuc(i) = j
        exit 
      end do sym_op1
    end do 
!
!  Next sym_op:  Over all (supplied) unit cells
!
   if (debugf) then
      write(iw_new, "(/,2a)") " Unit Cell Indices    Atom               Crystal coordinates             ", &
        "Cartesian coordinates            New Cartesian coordinates         New Crystal coordinates"
      write(iw_new, "(2a,/)") "           i  j  k                    a         b         c             ", &
        "x         y         z             x         y         z             a         b         c"                             
    end if
    do ixyz = 1, ncell
!
!  Put into new_coord the positions of all atoms in unit cell ixyz
!  after they have been moved by operation sym_op
!
!  The steps involved are:
!  1.  Move the atom into unit cell nijk(1:3,ixyz)
!  2.  Operate on the atom's position using operation "sym_op".
!
      do i = 1, numat
        new_coord(:, i) = coord(:, i) + nijk(:,ixyz)
      end do
      call rotsec (t1, .false., new_coord, sym_op)
!
!   Find out which unit cell each atom has been moved into.
!
      do ii = 1, numat
        do jj = 1, id
          k1 = Int (new_coord(jj, ii) + 40.1d0) - 40 - nijk_CUC(jj,ii)
          new_coord(jj, ii) = new_coord(jj, ii) - k1
!
!  incell holds the indices of the atoms of the  unit cell ixyz after they have been moved by operation sym_op
!  These indices are "relative" in the sense that the fractional coordinates are relative to those
!  of the atoms they came from in the CUC.  I.e., if an atom in the CUC has coordinates (0.1 -0.1 0.1)
!  then after the operation the indices of that atom will have (0,1,0) added to it.
!
          incell(jj, ii) = k1
        end do
        if (Mod (incell(1, ii) + incell(2, ii) + incell(3, ii), 2) /= 0 .and. bcc) then
          inquire(unit=iw_new, opened=opend) 
          if (.not. opend) open (unit = iw_new, file = trim(data_set_name)//".txt", &
          form = "FORMATTED", status = "UNKNOWN")  
          write(iw_new, '(/,a, i3, a,i3,a)') " For atom",ii, " Operation,", sym_op, &
          ", "//trim(title(sym_op))//", is not an allowed operation."
          write(iw_new, '(a)') &
          " Original position of atom     Transformed position of atom (in fractional unit cell coordinates)"
          write(iw_new,'(f12.4,2f8.4, f12.4, 2f8.4)') coord(:,ii) + nijk(:,ixyz), &
          (new_coord(jj,ii) + incell(jj,ii), jj = 1,id)
          jop = jop - 1
          cycle
        end if
        sym_op2: do jj = 1, numat
          do kk = 1, 3
            sum = new_coord(kk, ii) - nijk_CUC(kk,ii) - coord(kk, jj) + nijk_CUC(kk,jj)
            if (Abs (sum) > 1d-4) cycle sym_op2
          end do
!
!  Atom II has been moved to become atom new(ii) in unit cell incell(1-3,ii)
!
          mop(1, ii, sym_op) = jj
          mop(2, ii, sym_op) = itrans(1, ii)
          mop(3, ii, sym_op) = itrans(2, ii)
          mop(4, ii, sym_op) = itrans(3, ii) 
          new(ii) = jj
          new_coord(:,ii) = new_coord(:,ii) + nijk_CUC(:,jj) - nijk_CUC(:,ii)
          incell(:,ii) = incell(:,ii) - nijk_CUC(:,jj) + nijk_CUC(:,ii)
          exit
        end do sym_op2
        if (jj > numat) then
!
!  Error!
!
          inquire(unit=iw_new, opened=opend) 
          if (.not. opend) open (unit = iw_new, file = trim(data_set_name)//".txt", &
          form = "FORMATTED", status = "UNKNOWN")      
          write(line,'(a)')"Symmetry data were supplied from file """//trim(jobnam)//".ops"""
          write(*,'(/10x,a)')trim(line)
          write(iw_new,'(/10x,a)')trim(line)
          write(line,'(a,i3,a)') &
          "Operation", sym_op, " ("//trim(title(sym_op))//") is not an allowed operation."
          write(*,'(/10x,a)')trim(line)
          write(iw_new,'(/10x,a)')trim(line)
          write(line,'(a,i3,a)')"(It moved atom", ii, " into a position that is not part of the lattice)" 
          write(*,'(10x,a,/)')trim(line)
          write(iw_new,'(10x,a,/)')trim(line)
          write(line,'(a)')"(Edit """//trim(jobnam)//".mop"" or  """//trim(jobnam)//".ops"" or "
          write(*,'(10x,a)')trim(line)
          write(iw_new,'(10x,a)')trim(line)
          write(line,'(a)')"delete or rename file """//trim(jobnam)//".ops"" to avoid this error)"
          write(*,'(10x,a,/)')trim(line)
          write(iw_new,'(10x,a,/)')trim(line)
          stop 
        end if
      end do
      do ii = 1, numat
        if (new(ii) == 0) then
          write(*,*)ii
          stop "Error in SYMTRZ"
        end if
      end do
!
!   Now comes the difficult part:  Take each pair of atoms in unit
!   cell IXYZ interacting with the central unit cell, and work
!   out the new unit cell indices.
      do i = 1, numat
        do j = 1, numat
          il = nfirst(i)
          iu = nlast(i)
          jl = nfirst(j) + (ixyz - 1)*nvecs
          ju = nlast(j)  + (ixyz - 1)*nvecs
!
!  Extract that part of the secular determinant which represents atom I in 
!  the central unit cell interacting with atom J in unit cell IXYZ.
!
          i1 = 0
          do ii = il, iu
            i1 = i1 + 1
            j1 = 0
            do jj = jl, ju
              j1 = j1 + 1
              t2(i1, j1) = sec_det(ii, jj)
            end do
          end do
          nx = incell(1, j) - itrans(1, i)
          ny = incell(2, j) - itrans(2, i)
          nz = incell(3, j) - itrans(3, i)
!
!  Make sure that ix, iy, and iz are in the correct range
!
          if (nx <= -mr1/2) nx = nx + mr1
          if (ny <= -mr2/2) ny = ny + mr2
          if (nz <= -mr3/2) nz = nz + mr3
          if (nx >  mr1/2) nx = nx - mr1
          if (ny >  mr2/2) ny = ny - mr2
          if (nz >  mr3/2) nz = nz - mr3
          if (nxyz(nx, ny, nz) == 0) then
            if (bcc) then
              inquire(unit=iw_new, opened=opend) 
              if (.not. opend) open (unit = iw_new, file = trim(data_set_name)//".txt", &
                form = "FORMATTED", status = "UNKNOWN")
              write(iw_new,'(a)') "There is a known bug in some BCC systems.  Symmetry will not be used"
              do ii = 1, numat
                mop(1, ii, 1) = ii
                mop(2:4, ii, 1) = 0 
              end do
              iop = 1
              return             
            end if
            write(*,'(//10x,a,3i2,a)')" Unit cell (",nx,ny,nz,") does not exist!" 
            write(*,'(a,//)')"     (For more details, see '"//trim(data_set_name)//".txt')"
            write(iw_new,'(//10x,a,3i2,a,//)')" Unit cell (",nx,ny,nz,") does not exist!" 
            write(iw_new,'(a,i3,a,a,i3,a,i3,a,i3,a,3i3,a)') &
            " Faulty operation: Symmetry operation", sym_op, ", ("//trim(title(sym_op)), &
              "), Atom:",i," in CUC with atom:",j, " in unit cell",ixyz, " indices: (",nijk(:,ixyz),")"
            write(iw_new,'(a,i3,a,3i3,a)')" Indices of atom", i," in CUC:          (", itrans(:,i),")"
            write(iw_new,'(a,i3,a,i3,a,3i3,a)') &
            " Indices of atom", j," in unit cell", ixyz,": (", incell(:, j),")"
            stop
          end if
          newcel = nxyz(nx, ny, nz)
          call rotcel (t2,per_atom(i), per_atom(j), t1)
!
! Now put the rotated interaction matrix into the correct unit cell in new_sec_det.
!
          il = nfirst(newcuc(i))
          iu = nlast(newcuc(i))
          jl = nfirst(new(j)) + (newcel - 1)*nvecs 
          ju = nlast(new(j))  + (newcel - 1)*nvecs 
          i1 = 0
          do ii = il, iu
            i1 = i1 + 1
            j1 = 0
            do jj = jl, ju
              j1 = j1 + 1
              sum = t2(i1, j1)-sec_det(ii, jj)
              if (Abs (sum) > 0.5d0) then
                write(*,*)" Error detected. see file '"//trim(data_set_name)//".txt' for more information"
                call sleep(2)
                inquire(unit=iw_new, opened=opend) 
                if (.not. opend) open (unit = iw_new, file = trim(data_set_name)//".txt", &
                form = "FORMATTED", status = "UNKNOWN")
                write(iw_new, '(/,a)') " Secular determinant generated by symmetry operations may be faulty."
                write(iw_new, '(/,a,i3, a,i3,a,i3,a,i3,a,i3)') &
                  " Operation", sym_op," ("//trim(title(sym_op))//"), atom",i, &
                " in the CUC, and atom",j," in unit cell",ixyz
                write(iw_new,'(/,a,i3,a,3i3,a,3f8.3)')" New indices of atom",i," in the CUC:", itrans(:,i), &
                "   Coordinates:",new_cuc(:,i)
                write(iw_new,'(a,i3,a,3i3,a,3f8.3)')" New indices of atom",j,"            ", incell(:,j), &
                "   Coordinates:", new_coord(:,j)
                write(iw_new,'(a,i3,a,i3,a,i3,a,3i3)')" Atom",j," maps to atom",new(j), " in cell", newcel, &
                " at indices", incell(:,j)
                write(iw_new,'(/,a,i3,a,i3,a, 3i3)') " F(s,s) terms for unit cells involving atom", i, &
                " in the CUC and", new(j), " in all unit cells"
                write(iw_new,'(a,i3,a,i3,a)')" Unit cell  Indices   Atom ",new(j), "   F(ss)     Atom", &
                new(j), "  fractional coords"
                do mm = 1, ncell
                  write(iw_new, '(i6,i7,2i3,i6,f14.4,f12.4,2f8.4)')mm, nijk(:,mm),(mm - 1)*numat + new(j), & 
                    sec_det(nfirst(newcuc(i)), (mm - 1)*nvecs + nfirst(new(j))), &
                    coord(1,new(j)) + nijk(1,mm), coord(2,new(j)) + nijk(2,mm), coord(3,new(j)) + nijk(3,mm)
                end do
                if (.not. debugs .or. .not. debugf) write(iw_new,'(a)') &
                " Use keywords SYMTRZ and ROTFOK for more information."
                stop
              else
                new_sec_det(ii, jj) = new_sec_det(ii, jj) + t2(i1, j1)
              end if 
            end do
          end do
        end do
      end do
    end do
  end do big_sym_op
  iop = jop
  if (jop == 0) then
    write(iw_new, '(a)') "System has NO symmetry operations"
    stop  
  end if
!
!  Re-normalize new_sec_det, then copy it back into sec_det
!  Identify the largest error
!
  sum = 0.d0
  do i = 1, nvecs
    do j = 1, ncell*nvecs
      new_sec_det(i, j) = new_sec_det(i, j)/iop
      if (Abs (new_sec_det(i, j) - sec_det(i, j)) > sum) then 
        sum = Abs (new_sec_det(i, j) - sec_det(i, j))
        ix = i
        jx = j
      end if
      sec_det(i, j) = new_sec_det(i, j)
    end do
  end do
  if (sum > 1.d-1) then
    inquire(unit=iw_new, opened=opend) 
    if (.not. opend) open (unit = iw_new, file = trim(data_set_name)//".txt", form = "FORMATTED", status = "UNKNOWN")
    write(iw_new, '(a,f10.2,a,2i4)') " Largest Difference in symmetrized secular determinant:", sum, " at location", ix, jx
  end if
  return
end subroutine symtrz
