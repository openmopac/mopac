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
  use Common_arrays_C, only : grad, & ! gradients of heat
    p, & ! total density matrix
    q, & ! partial charges
    nat, & ! atomic number of each atom
    grad, & ! heat gradients
    bondab ! bond order matrix in packed triangular format
  use molkst_C, only : escf, & ! heat of formation
    numat, & ! number of real atoms
    id, & ! number of lattice vectors
    voigt, & ! Voigt stress tensor
    mozyme, & ! logical flag for MOZYME calculations
    use_disk ! logical flag to enable disk access
  use MOZYME_C, only : iorbs ! number of atomic orbitals for each atom
  use parameters_C, only : tore ! number of valence electrons per element
  use to_screen_C, only : dip ! dipole moments
  implicit none

contains

  ! save properties and clean up after a MOPAC/MOZYME calculation
  module subroutine mopac_finalize(properties)
    type(mopac_properties), intent(out) :: properties
    integer, external :: ijbo
    integer :: i, j, k, kk, kl, ku, io, jo
    double precision :: valenc, sum

    ! deallocate any prior arrays
    if (allocated(properties%coord_deriv)) deallocate(properties%coord_deriv)
    if (allocated(properties%lattice_deriv)) deallocate(properties%lattice_deriv)
    if (allocated(properties%charge)) deallocate(properties%charge)
    if (allocated(properties%bond_index)) deallocate(properties%bond_index)
    if (allocated(properties%bond_atom)) deallocate(properties%bond_atom)
    if (allocated(properties%bond_order)) deallocate(properties%bond_order)
    ! save properties
    properties%heat = escf
    allocate(properties%coord_deriv(3*numat))
    properties%coord_deriv = grad(:3*numat)
    if (id > 0) then
      allocate(properties%lattice_deriv(3*id))
      properties%lattice_deriv = grad(3*numat+1:3*(numat+id))
    end if
    call chrge (p, q)
    q(:numat) = tore(nat(:numat)) - q(:numat)
    allocate(properties%charge(numat))
    properties%charge = q(:numat)
    properties%dipole = dip(:3,3)
    properties%stress = voigt
    ! prune bond order matrix
    allocate(properties%bond_index(numat+1))
    if (mozyme) then
      ! 1st pass to populate bond_index
      properties%bond_index(i) = 1
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
      properties%bond_index(i) = 1
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