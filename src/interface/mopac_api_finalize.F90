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
    numat, & ! number of real atoms
    nvar, & ! number of coordinates to be optimized
    id, & ! number of lattice vectors
    keywrd, & ! keyword string to adjust MOPAC behavior
    voigt, & ! Voigt stress tensor
    mozyme, & ! logical flag for MOZYME calculations
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
    integer, external :: ijbo
    double precision, external :: dipole, dipole_for_MOZYME
    integer :: i, j, k, kk, kl, ku, io, jo, natom_move, nlattice_move
    double precision :: valenc, sum, dumy(3)

    ! close dummy output file to free up /dev/null
    close(iw)
    ! deallocate any prior arrays
    if (allocated(properties%coord_update)) deallocate(properties%coord_update)
    if (allocated(properties%coord_deriv)) deallocate(properties%coord_deriv)
    if (allocated(properties%freq)) deallocate(properties%freq)
    if (allocated(properties%disp)) deallocate(properties%disp)
    if (allocated(properties%charge)) deallocate(properties%charge)
    if (allocated(properties%bond_index)) deallocate(properties%bond_index)
    if (allocated(properties%bond_atom)) deallocate(properties%bond_atom)
    if (allocated(properties%bond_order)) deallocate(properties%bond_order)
    if (allocated(properties%lattice_update)) deallocate(properties%lattice_update)
    if (allocated(properties%lattice_deriv)) deallocate(properties%lattice_deriv)
    ! trigger charge & dipole calculation
    call chrge (p, q)
    q(:numat) = tore(nat(:numat)) - q(:numat)
    if (mozyme) then
      sum = dipole_for_MOZYME(dumy, 2)
      properties%dipole = dumy
    else
      sum = dipole(p, xparam, dumy, 1)
      properties%dipole = dip(:3,3)
    end if
    ! save basic properties
    properties%heat = escf
    allocate(properties%charge(numat))
    properties%charge = q(:numat)
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
    allocate(properties%coord_update(3*natom_move))
    properties%coord_update = xparam(:3*natom_move)
    allocate(properties%coord_deriv(3*natom_move))
    properties%coord_deriv = grad(:3*natom_move)
    if (nlattice_move > 0) then
      allocate(properties%lattice_update(3*nlattice_move))
      properties%lattice_update = xparam(3*natom_move+1:)
      allocate(properties%lattice_deriv(3*nlattice_move))
      properties%lattice_deriv = grad(3*natom_move+1:)
    end if
    ! save vibrational properties if available
    if (index(keywrd, " FORCE") /= 0) then
      properties%calc_vibe = .true.
      allocate(properties%freq(nvar))
      allocate(properties%disp(nvar,nvar))
      properties%freq = freq
      properties%disp = reshape(cnorml,[nvar, nvar])
    else
      properties%calc_vibe = .false.
    end if
    ! prune bond order matrix
    allocate(properties%bond_index(numat+1))
    if (mozyme) then
      ! 1st pass to populate bond_index
      properties%bond_index(1) = 1
      do i = 1, numat
        io = iorbs(i)
        properties%bond_index(i+1) = properties%bond_index(i)
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
        if (valenc > 0.001d0) then
          properties%bond_index(i+1) = properties%bond_index(i+1) + 1
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
            if (sum > 0.001d0) then
              properties%bond_index(i+1) = properties%bond_index(i+1) + 1
            end if
          end if
        end do
      end do
      ! 2nd pass to populate bond_atom and bond_order
      allocate(properties%bond_atom(properties%bond_index(numat+1)))
      allocate(properties%bond_order(properties%bond_index(numat+1)))
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
        kk = properties%bond_index(i)
        do j = 1, numat
          jo = iorbs(j)
          if (i /= j .and. ijbo (i, j) >= 0) then
            kl = ijbo (i, j) + 1
            ku = kl + io * jo - 1
            sum = 0.d0
            do k = kl, ku
              sum = sum + p(k) ** 2
            end do
            if (sum > 0.001d0) then
              properties%bond_atom(kk) = j
              properties%bond_order(kk) = sum
              kk = kk + 1
            end if
          else if (valenc > 0.001d0) then
            properties%bond_atom(kk) = j
            properties%bond_order(kk) = valenc
            kk = kk + 1
          end if
        end do
      end do
    else
      call bonds()
      ! 1st pass to populate bond_index
      properties%bond_index(1) = 1
      do i = 1, numat
        properties%bond_index(i+1) = properties%bond_index(i)
        ku = i*(i-1)/2 + 1
        kl = (i+1)*(i+2)/2 - 1
        do j = 1, i
          if (bondab(ku) > 0.001d0) then
            properties%bond_index(i+1) = properties%bond_index(i+1) + 1
          end if
          ku = ku + 1
        end do
        do j = i+1, numat
          if (bondab(kl) > 0.001d0) then
            properties%bond_index(i+1) = properties%bond_index(i+1) + 1
          end if
          kl = kl + j
        end do
      end do
      ! 2nd pass to populate bond_atom and bond_order
      allocate(properties%bond_atom(properties%bond_index(numat+1)))
      allocate(properties%bond_order(properties%bond_index(numat+1)))
      do i = 1, numat
        ku = i*(i-1)/2 + 1
        kl = (i+1)*(i+2)/2 - 1
        kk = properties%bond_index(i)
        do j = 1, i
          if (bondab(ku) > 0.001d0) then
            properties%bond_atom(kk) = j
            properties%bond_order(kk) = bondab(ku)
            kk = kk + 1
          end if
          ku = ku + 1
        end do
        do j = i+1, numat
          if (bondab(kl) > 0.001d0) then
            properties%bond_atom(kk) = j
            properties%bond_order(kk) = bondab(kl)
            kk = kk + 1
          end if
          kl = kl + j
        end do
      end do
    end if
    ! deallocate memory
    call setup_mopac_arrays(0,0)
    if (mozyme) call delete_MOZYME_arrays()
    ! turn use_disk back on
    use_disk = .true.
    ! mark the job as successful
    properties%status = 0
  end subroutine mopac_finalize

end submodule mopac_api_finalize
