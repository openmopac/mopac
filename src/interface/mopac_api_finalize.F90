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

submodule (mopac_api:mopac_api_operations) mopac_api_finalize
  use iso_c_binding, only: c_int, c_double, c_char, c_ptr, c_loc, c_null_char, c_null_ptr

  use chanel_C, only : iw ! file handle for main output file
  use Common_arrays_C, only : xparam, & ! values of coordinates undergoing optimization
    geo, & ! raw coordinates of atoms
    grad, & ! gradients of heat
    p, & ! total density matrix
    q, & ! partial charges
    nat, & ! atomic number of each atom
    loc, & ! indices of atoms and coordinates marked for optimization
    bondab ! bond order matrix in packed triangular format
  use molkst_C, only : escf, & ! heat of formation
    moperr, & ! error status
    numat, & ! number of real atoms
    nvar, & ! number of coordinates to be optimized
    id, & ! number of lattice vectors
    keywrd, & ! keyword string to adjust MOPAC behavior
    voigt, & ! Voigt stress tensor
    mozyme, & ! logical flag for MOZYME calculations
    dummy, & ! dummy integer
    errtxt, & ! most recent error message
    use_disk ! logical flag to enable disk access
  use MOZYME_C, only : iorbs ! number of atomic orbitals for each atom
  use parameters_C, only : tore ! number of valence electrons per element
  use to_screen_C, only : dip, & ! dipole moments
    freq, cnorml ! vibrational properties

  implicit none

contains

  ! save properties and clean up after a MOPAC/MOZYME calculation
  module subroutine mopac_finalize(properties)
    type(mopac_properties), intent(out) :: properties
    integer :: status, i, j, size
    character(kind=c_int), pointer :: cptr(:)
    type(c_ptr), pointer :: pptr(:)

    ! close dummy output file to free up /dev/null
    close(iw)

    ! record properties
    if (.not. moperr) call mopac_record(properties)

    ! collect error messages
    if (moperr) then
      call summary("",0)
      properties%nerror = dummy
      allocate(pptr(properties%nerror), stat=status)
      properties%error_msg = c_loc(pptr)
      if (status /= 0) then
        properties%nerror = -properties%nerror
        return
      else
        do i=1, properties%nerror
          call summary("",-i)
          size = len_trim(errtxt)
          allocate(cptr(size+1), stat=status)
          pptr(i) = c_loc(cptr)
          if (status /= 0) then
            properties%nerror = -properties%nerror
            return    
          end if
          do j=1, size
            cptr(j) = errtxt(j)
          end do
          cptr(size+1) = c_null_char
        end do
      end if
      call summary("",-abs(properties%nerror)-1)
    else
      properties%nerror = 0
    end if

    ! deallocate memory
    call setup_mopac_arrays(0,0)
    if (mozyme) call delete_MOZYME_arrays()
    ! turn use_disk back on
    use_disk = .true.
  end subroutine mopac_finalize

  subroutine mopac_record(properties)
    type(mopac_properties), intent(out) :: properties
    integer, external :: ijbo
    double precision, external :: dipole, dipole_for_MOZYME
    integer :: status, i, j, k, kk, kl, ku, io, jo, natom_move, nlattice_move
    double precision :: valenc, sum, dumy(3)
    integer(c_int), pointer :: bond_index(:), bond_atom(:)
    real(c_double), pointer :: rptr(:), bond_order(:)

    ! trigger charge & dipole calculation
    call chrge (p, q)
    q(:numat) = tore(nat(:numat)) - q(:numat)
    if (id == 0) then
      if (mozyme) then
        sum = dipole_for_MOZYME(dumy, 2)
        properties%dipole = dumy
      else
        sum = dipole(p, xparam, dumy, 1)
        properties%dipole = dip(:3,3)
      end if
    else
      properties%dipole = 0.d0
    end if
    ! save basic properties
    properties%heat = escf
    allocate(rptr(numat), stat=status)
    if (status /= 0) then
      call mopend("Failed to allocate memory in MOPAC_FINALIZE")
      return
    end if
    rptr = q(:numat)
    properties%charge = c_loc(rptr)
    properties%stress = voigt
    natom_move = nvar/3
    nlattice_move = 0
    do i=nvar/3-id, nvar/3
      if (nat(loc(1,3*i)) == 107) then
        natom_move = natom_move - 1
        nlattice_move = nlattice_move + 1
      end if
    end do
    ! save properties of moveable coordinates
    allocate(rptr(3*natom_move), stat=status)
    if (status /= 0) then
      call mopend("Failed to allocate memory in MOPAC_FINALIZE")
      return
    end if
    rptr = xparam(:3*natom_move)
    properties%coord_update = c_loc(rptr)
    allocate(rptr(3*natom_move), stat=status)
    if (status /= 0) then
      call mopend("Failed to allocate memory in MOPAC_FINALIZE")
      return
    end if
    rptr = grad(:3*natom_move)
    properties%coord_deriv = c_loc(rptr)
    if (nlattice_move > 0) then
      allocate(rptr(3*nlattice_move), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOPAC_FINALIZE")
        return
      end if
      rptr = xparam(3*natom_move+1:)
      properties%lattice_update = c_loc(rptr)
      allocate(rptr(3*nlattice_move), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOPAC_FINALIZE")
        return
      end if
      rptr = grad(3*natom_move+1:)
      properties%lattice_deriv = c_loc(rptr)
    end if
    ! save vibrational properties if available
    if (index(keywrd, " FORCETS") /= 0) then
      allocate(rptr(nvar), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOPAC_FINALIZE")
        return
      end if
      rptr = freq
      properties%freq = c_loc(rptr)
      allocate(rptr(nvar*nvar), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOPAC_FINALIZE")
        return
      end if
      rptr = cnorml
      properties%disp = c_loc(rptr)
    else
      properties%freq = c_null_ptr
      properties%disp = c_null_ptr
    end if
    ! prune bond order matrix
    allocate(bond_index(numat+1), stat=status)
    properties%bond_index = c_loc(bond_index)
    if (status /= 0) then
      call mopend("Failed to allocate memory in MOPAC_FINALIZE")
      return
    end if
    if (mozyme) then
      ! 1st pass to populate bond_index
      bond_index(1) = 1
      do i = 1, numat
        io = iorbs(i)
        bond_index(i+1) = bond_index(i)
        valenc = 0.d0
        if (io == 1) then
          kk = ijbo (i, i) + 1
          valenc = 2.d0 * p(kk) - p(kk) ** 2
        else
          kk = ijbo (i, i)
          do j = 1, io
            do k = 1, j
              kk = kk + 1
              valenc = valenc - 2.d0*p(kk) * p(kk)
            end do
            valenc = valenc + p(kk) * p(kk)
            valenc = valenc + 2.d0 * p(kk)
          end do
        end if
        if (valenc > 0.01d0) then
          bond_index(i+1) = bond_index(i+1) + 1
        end if
        do j = 1, numat
          jo = iorbs(j)
          if (i /= j .and. ijbo (i, j) >= 0) then
            kl = ijbo (i, j) + 1
            ku = kl + io * jo - 1
            sum = 0.d0
            do k = kl, ku
              sum = sum + p(k) ** 2
            end do
            if (sum > 0.01d0) then
              bond_index(i+1) = bond_index(i+1) + 1
            end if
          end if
        end do
      end do
      ! 2nd pass to populate bond_atom and bond_order
      allocate(properties%bond_atom(properties%bond_index(numat+1)), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOPAC_FINALIZE")
        return
      end if
      properties%bond_atom = c_loc(bond_atom)
      allocate(properties%bond_order(properties%bond_index(numat+1)), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOPAC_FINALIZE")
        return
      end if
      properties%bond_order = c_loc(bond_order)
      do i = 1, numat
        io = iorbs(i)
        valenc = 0.d0
        if (io == 1) then
          kk = ijbo (i, i) + 1
          valenc = 2.d0 * p(kk) - p(kk) ** 2
        else
          kk = ijbo (i, i)
          do j = 1, io
            do k = 1, j
              kk = kk + 1
              valenc = valenc - 2.d0*p(kk) * p(kk)
            end do
            valenc = valenc + p(kk) * p(kk)
            valenc = valenc + 2.d0 * p(kk)
          end do
        end if
        kk = bond_index(i)
        do j = 1, numat
          jo = iorbs(j)
          if (i /= j .and. ijbo (i, j) >= 0) then
            kl = ijbo (i, j) + 1
            ku = kl + io * jo - 1
            sum = 0.d0
            do k = kl, ku
              sum = sum + p(k) ** 2
            end do
            if (sum > 0.01d0) then
              bond_atom(kk) = j
              bond_order(kk) = sum
              kk = kk + 1
            end if
          else if (valenc > 0.01d0) then
            bond_atom(kk) = j
            bond_order(kk) = valenc
            kk = kk + 1
          end if
        end do
      end do
    else
      call bonds()
      ! 1st pass to populate bond_index
      bond_index(1) = 1
      do i = 1, numat
        bond_index(i+1) = bond_index(i)
        ku = i*(i-1)/2 + 1
        kl = (i+1)*(i+2)/2 - 1
        do j = 1, i
          if (bondab(ku) > 0.01d0) then
            bond_index(i+1) = bond_index(i+1) + 1
          end if
          ku = ku + 1
        end do
        do j = i+1, numat
          if (bondab(kl) > 0.01d0) then
            bond_index(i+1) = bond_index(i+1) + 1
          end if
          kl = kl + j
        end do
      end do
      ! 2nd pass to populate bond_atom and bond_order
      allocate(bond_atom(properties%bond_index(numat+1)), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOPAC_FINALIZE")
        return
      end if
      properties%bond_atom = c_loc(bond_atom)
      allocate(bond_order(properties%bond_index(numat+1)), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOPAC_FINALIZE")
        return
      end if
      properties%bond_order = c_loc(bond_order)
      do i = 1, numat
        ku = i*(i-1)/2 + 1
        kl = (i+1)*(i+2)/2 - 1
        kk = bond_index(i)
        do j = 1, i
          if (bondab(ku) > 0.01d0) then
            bond_atom(kk) = j
            bond_order(kk) = bondab(ku)
            kk = kk + 1
          end if
          ku = ku + 1
        end do
        do j = i+1, numat
          if (bondab(kl) > 0.01d0) then
            bond_atom(kk) = j
            bond_order(kk) = bondab(kl)
            kk = kk + 1
          end if
          kl = kl + j
        end do
      end do
    end if
  end subroutine mopac_record

end submodule mopac_api_finalize
