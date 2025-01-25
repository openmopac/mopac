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

submodule (mopac_api:mopac_api_operations) mopac_api_finalize
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
    type(c_ptr), allocatable :: pptr(:)

    ! record properties
    if (.not. moperr) call mopac_record(properties)

    ! collect error messages & assign NULL pointers for memory safety
    if (moperr) then
      properties%charge = c_null_ptr
      properties%coord_update = c_null_ptr
      properties%coord_deriv = c_null_ptr
      properties%lattice_update = c_null_ptr
      properties%lattice_deriv = c_null_ptr
      properties%freq = c_null_ptr
      properties%disp = c_null_ptr
      properties%bond_index = c_null_ptr
      properties%bond_atom = c_null_ptr
      properties%bond_order = c_null_ptr
      call summary("",0)
      properties%nerror = dummy
      allocate(pptr(properties%nerror), stat=status)
      if (status /= 0) then
        properties%nerror = -properties%nerror
        properties%error_msg = c_null_ptr
        return
      else
        do i=1, properties%nerror
          call summary("",-i)
          size = len_trim(errtxt) + 1
          pptr(i) = create_copy(errtxt, [size])
        end do
      end if
      properties%error_msg = create_copy(pptr, [properties%nerror])
      call summary("",-abs(properties%nerror)-1)
    else
      properties%nerror = 0
    end if

    ! deallocate memory
    call setup_mopac_arrays(0,0)
    if (mozyme) call delete_MOZYME_arrays()

    ! turn use_disk back on
    use_disk = .true.

    ! close dummy output file to free up /dev/null
    close(iw)
  end subroutine mopac_finalize

  subroutine mopac_record(properties)
    type(mopac_properties), intent(out) :: properties
    integer, external :: ijbo
    double precision, external :: dipole, dipole_for_MOZYME
    integer :: status, i, j, k, kk, kl, ku, io, jo, natom_move, nlattice_move
    double precision :: valenc, sum, dumy(3)
    integer(c_int), pointer :: bond_index(:), bond_atom(:)
    real(c_double), pointer :: bond_order(:)

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
    properties%charge = create_copy(q, [numat])
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
    properties%coord_update = create_copy(xparam, [3*natom_move])
    properties%coord_deriv = create_copy(grad, [3*natom_move])
    if (nlattice_move > 0) then
      properties%lattice_update = create_copy(xparam(3*natom_move+1:), [3*nlattice_move])
      properties%lattice_deriv = create_copy(grad(3*natom_move+1:), [3*nlattice_move])
    else
      properties%lattice_update = c_null_ptr
      properties%lattice_deriv = c_null_ptr
    end if
    ! save vibrational properties if available
    if (index(keywrd, " FORCETS") /= 0) then
      properties%freq = create_copy(freq, [nvar])
      properties%disp = create_copy(cnorml, [nvar*nvar])
    else
      properties%freq = c_null_ptr
      properties%disp = c_null_ptr
    end if
    ! prune bond order matrix
    properties%bond_index = create_int(numat+1)
    call c_f_pointer(properties%bond_index, bond_index, [numat+1])
    if (mozyme) then
      ! 1st pass to populate bond_index
      bond_index(1) = 0
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
      properties%bond_atom = create_int(bond_index(numat+1))
      call c_f_pointer(properties%bond_atom, bond_atom, [bond_index(numat+1)])
      properties%bond_order = create_real(bond_index(numat+1))
      call c_f_pointer(properties%bond_order, bond_order, [bond_index(numat+1)])
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
        kk = bond_index(i) + 1
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
      bond_index(1) = 0
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
      properties%bond_atom = create_int(bond_index(numat+1))
      call c_f_pointer(properties%bond_atom, bond_atom, [bond_index(numat+1)])
      properties%bond_order = create_real(bond_index(numat+1))
      call c_f_pointer(properties%bond_order, bond_order, [bond_index(numat+1)])
      do i = 1, numat
        ku = i*(i-1)/2 + 1
        kl = (i+1)*(i+2)/2 - 1
        kk = bond_index(i) + 1
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
