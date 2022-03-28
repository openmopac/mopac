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

subroutine eigen (write_gpt, write_out)
    use molkst_C, only: norbs, numat, nelecs, keywrd
    use chanel_C, only: igpt, iw
    use elemts_C, only: elemnt
    use common_arrays_C, only: eigs, nfirst, nlast, nat, c
    use MOZYME_C, only : icocc, icvir, ncocc, ncvir, nvirtual, noccupied, &
      & nnce, nncf, cocc, cvir, ncf, nce,  cocc_dim, &
      & cvir_dim, icocc_dim, icvir_dim
    implicit none
    logical, intent (in) :: write_gpt, write_out
    integer :: alloc_stat, vec_ao_size
    integer :: i, j, k, l, m
    integer :: print_set, print_first_orbital
    integer :: print_nocc, print_nvir, print_ntot, print_shift
    integer :: print_npanels, print_ncolumns, print_columns_left
    integer, parameter :: print_columns_per_panel = 8
    integer, parameter :: do_occ = 1, do_vir = 2
    integer, parameter :: default_print_nocc = 8
    integer, parameter :: default_print_nvir = 8
    character (len=248) :: tmpkey
    character (len=2), dimension (9), save :: print_orbital_name
    double precision, external :: reada
    data print_orbital_name / "S ", "Px", "Py", "Pz", "x2", "xz", "z2", "yz", "xy" /
    noccupied = nelecs/2
    nvirtual = norbs - noccupied
   !
   ! ==============================
   ! Write eigenvectors to gpt file
   ! ==============================
    if (write_gpt) then
      !
      ! Allocate local arrays
      allocate (c(norbs, norbs), stat=alloc_stat)
      if (alloc_stat /= 0) then
        call memory_error ("eigen")
        goto 1100
      end if
  !
  ! Convert LMOs to Eigenvectors
  !
      call lmo_to_eigenvectors(noccupied, ncf, nncf, ncocc, noccupied, &
           & icocc, icocc_dim, cocc, cocc_dim, eigs, c)
      call lmo_to_eigenvectors(nvirtual, nce, nnce, ncvir, nvirtual, icvir, &
           & icvir_dim, cvir, cvir_dim, eigs(noccupied + 1), c(1, noccupied + 1))
      !
      ! Write out results
      write (igpt) c
    end if
   !
   ! =================================================
   ! Write eigenvalues and eigenvectors to output file
   ! =================================================
    if (write_out) then
      !
      ! Determine how much needs to be written out
      if (Index (keywrd, " ALLVEC") /= 0) then
        !
        ! Everything requested
        print_nocc = noccupied
        print_nvir = nvirtual
      else if ((Index (keywrd, " VECTORS(") /= 0) .or. &
             & (Index (keywrd, " VECTORS=(") /= 0)) then
        !
        ! User has requested limits
        !
        ! Extract user limits from keyword string
        !
        ! Build temporary keyword string, and erase all text from it
        ! except VECTORS data
        tmpkey = trim(keywrd)
        i = Index (tmpkey, "VECTORS=(") + Index (tmpkey, "VECTORS(")
        tmpkey (:i) = " "
        j = Index (tmpkey, ")")
        tmpkey (j:) = " "
        ! Extract number of occupied orbitals to be written out
        print_nocc = Nint (reada (tmpkey, i))
        ! Extract number of virtual orbitals to be written out
        i = Index (tmpkey, ",") + Index (tmpkey, ":")
        if (i /= 0) then
          print_nvir = Nint (reada (tmpkey, i + 1))
        else
          print_nvir = default_print_nvir
        end if
      else
        !
        ! Use defaults
        print_nocc = default_print_nocc
        print_nvir = default_print_nvir
      end if
     !
     ! Make sure number of orbitals requested is no larger than
     ! the number available
      print_nocc = Max (0, Min (print_nocc, noccupied))
      print_nvir = Max (0, Min (print_nvir, nvirtual))
      print_ntot = print_nocc + print_nvir
      if (print_ntot > 0) then
        if (.not. write_gpt) then
          !
          ! Allocate local arrays
          vec_ao_size = Max (noccupied, nvirtual)
          allocate (c(norbs, vec_ao_size), stat=alloc_stat)
          if (alloc_stat /= 0) then
            call memory_error ("eigen")
            goto 1100
          end if
        end if
        c = 0.d0
        do print_set = do_occ, do_vir
          if (.not. write_gpt) then
            !
            ! Convert LMOs to Eigenvectors
            !
            if (print_set == do_occ) then
              call lmo_to_eigenvectors(noccupied, ncf, nncf, ncocc, noccupied, &
           & icocc, icocc_dim, cocc, cocc_dim, eigs, c)
              print_shift = noccupied - print_nocc
            else
              call lmo_to_eigenvectors(nvirtual, nce, nnce, ncvir, nvirtual, &
           & icvir, icvir_dim, cvir, cvir_dim, eigs, c)
              print_shift = 0
            end if
          else
            if (print_set == do_occ) then
              print_shift = noccupied - print_nocc
            else
              print_shift = noccupied
            end if
          end if
          !
          ! Write out results
          !
          ! Find number of print panels
          if (print_set == do_occ) then
            write (iw, "(/'  Occupied Orbitals'/)")
            print_first_orbital = noccupied - print_nocc
            print_npanels = 1 + (print_nocc - 1) / print_columns_per_panel
            print_columns_left = print_nocc
          else
            write (iw, "(/'  Virtual Orbitals'/)")
            print_first_orbital = noccupied
            print_npanels = 1 + (print_nvir - 1) / print_columns_per_panel
            print_columns_left = print_nvir
          end if
          !
          ! For each panel ...
          do i = 1, print_npanels
            !
            ! Evaluate how many columns to print in this panel
            if (print_columns_left >= print_columns_per_panel) then
              print_ncolumns = print_columns_per_panel
            else
              print_ncolumns = print_columns_left
            end if
            print_columns_left = print_columns_left - print_ncolumns
            !
            ! Print header for this panel
            write (iw, "(//'    Root No.')", advance="no")
            do j = 1, print_ncolumns
              write (iw, "(i6,'  ')", advance="no") print_first_orbital + &
                  & (i - 1) * print_columns_per_panel + j
            end do
            write (iw, "(/)")
            !
            ! Print eigenvalues
            write (iw, "('            ')", advance="no")
            do j = 1, print_ncolumns
              write (iw, "(f8.3)", advance="no") eigs(print_shift + &
                  & (i - 1) * print_columns_per_panel + j)
            end do
            write (iw, "(/)")
            !
            ! Print Eigenvectors
            k = 1
            do l = 1, numat
              do m = 1, 1 + nlast(l) - nfirst(l)
                write (iw, "(' ',a2,' ',a2,' ',i5)", advance="no") &
                    & print_orbital_name(m), elemnt(nat(l)), l
                do j = 1, print_ncolumns
                  write (iw, "(f8.4)", advance="no") c &
                 & (k, print_shift + (i - 1) * print_columns_per_panel + j)
                end do
                write (iw, "(' ')")
                k = k + 1
              end do
            end do
          end do
        end do
      end if
    end if
1100 continue
end subroutine eigen
subroutine FlmoFromFao (nlmo, ncx, nncx, ncxxx, n_dim, &
     & latoms, ws, icxxx, icxxx_dim, cxxx, cxxx_dim, &
     & ac, fmo_expanded)
  use molkst_C, only: numat, norbs
  use common_arrays_C, only : f, nfirst, nlast
  use MOZYME_C, only : lijbo, nijbo, iorbs
  !
  !.. Implicit Declarations ..
  implicit none
  !
  !.. Formal Arguments ..
  integer, intent (in) :: nlmo, icxxx_dim, cxxx_dim, n_dim
  logical, dimension (numat), intent (out) :: latoms
  integer, dimension (n_dim), intent (in) :: ncx, ncxxx, nncx
  integer, dimension (icxxx_dim), intent (in) :: icxxx
  double precision, dimension (nlmo * nlmo), intent (out) :: fmo_expanded
  double precision, dimension (norbs), intent (out) :: ws
  double precision, dimension (cxxx_dim), intent (in) :: cxxx
  double precision, dimension (numat), intent (out) :: ac
  integer :: i, i1, i2, i4, ii, j, j1, j2, &
       & jj, jx, k, k1, kk, kl, l, loopi, loopj, kj
  integer :: i1j1, i1j2, i2j1, i2j2
  integer :: k_fmo
  double precision :: flim, sum, cutoff
  integer, external :: ijbo

  ! Tolerance parameters
  flim = 1.0d-8
  cutoff = 1.0d-8
  !
  ! For each lmo do ...
  fmo_expanded = 0.0d0
  k_fmo = 1
  do i = 1, nlmo
    loopi = ncxxx(i)
    !
    ! Evaluate weight on each atom from lmo and build list of flags
    ! indicating if atoms exist in the lmo
    do j = 1, numat
      latoms(j) = .false.
    end do
    l = 0
    do j = nncx(i) + 1, nncx(i) + ncx(i)
      j1 = icxxx(j)
      sum = 0.d0
      do k = l + 1, l + iorbs(j1)
        sum = sum + cxxx(k+loopi) ** 2
      end do
      ac(j1) = sum
      latoms(j1) = .true.
      l = l + iorbs(j1)
    end do
    !
    ! Evaluate F * C
    !
    ! For each atom in the lmo do ...
    do jj = nncx(i) + 1, nncx(i) + ncx(i)
      j1 = icxxx(jj)
      do jx = 1, iorbs(j1)
        ws(nfirst(j1)+jx-1) = 0.0d00
      end do
      kl = loopi
      !
      ! For each neighboring atom in the lmo do ...
      do kk = nncx(i) + 1, nncx(i) + ncx(i)
        k1 = icxxx(kk)
        if (lijbo) then
          kj = nijbo (k1, j1)
        else
          kj = ijbo (k1, j1)
        end if
        !
        ! If the neighbor contributes to the lmo then evaluate
        ! contribution to F * C
        if (kj >= 0) then
          if (ac(k1) > cutoff) then
            if (iorbs(k1) .eq. 1 .and. iorbs(j1) .eq. 1) then
              ws(nfirst(j1)) = ws(nfirst(j1)) + f(kj+1) * cxxx(kl+1)
            else
              if (k1 > j1) then
                ii = kj
                do i4 = 1, iorbs(k1)
                  do jx = 1, iorbs(j1)
                    ii = ii + 1
                    ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + f(ii) &
                         & * cxxx(kl+i4)
                  end do
                end do
              else if (k1 < j1) then
                ii = kj
                do jx = 1, iorbs(j1)
                  do i4 = 1, iorbs(k1)
                    ii = ii + 1
                    ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + f(ii) &
                         & * cxxx(kl+i4)
                  end do
                end do
              else
                do jx = 1, iorbs(j1)
                  do i4 = 1, iorbs(j1)
                    ii = kj + (jx*(jx-1)) / 2 + i4
                    if (i4 .gt. jx) ii = kj + (i4*(i4-1)) / 2 + jx
                    ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + f(ii) &
                         & * cxxx(kl+i4)
                  end do
                end do
              end if
            end if
          end if
        end if
        kl = kl + iorbs(k1)
      end do
    end do
    !
    ! Evaluate C * (F * C)
    i1 = icxxx(nncx(i)+1)
    if (ncx(i) > 1) then
      i2 = icxxx(nncx(i)+2)
    else
      i2 = i1
    end if
    !
    ! For each second orbital do ...
    do j = 1, i
        !
        ! FAST TEST TO SEE IF THE INTEGRAL IS WORTH EVALUATING
        j1 = icxxx(nncx(j)+1)
        if (ncx(j) > 1) then
          j2 = icxxx(nncx(j)+2)
        else
          j2 = j1
        end if
        if (lijbo) then
          i1j1 = nijbo (i1, j1)
          i1j2 = nijbo (i1, j2)
          i2j1 = nijbo (i2, j1)
          i2j2 = nijbo (i2, j2)
        else
          i1j1 = ijbo (i1, j1)
          i1j2 = ijbo (i1, j2)
          i2j1 = ijbo (i2, j1)
          i2j2 = ijbo (i2, j2)
        end if
        if (i1j1 >= 0 .or. i1j2 >= 0 .or. i2j1 >= 0 .or. i2j2 >= 0) &
             & then
          sum = 0.d0
          if (i1j1 >= 0) then
            sum = Abs (f(i1j1+1))
          end if
          if (i1j2 >= 0) then
            sum = sum + Abs (f(i1j2+1))
          end if
          if (i2j1 >= 0) then
            sum = sum + Abs (f(i2j1+1))
          end if
          if (i2j2 >= 0) then
            sum = sum + Abs (f(i2j2+1))
          end if
          !
          ! If test passed, build C * (F * C)
          if (sum >= flim) then
            loopj = ncxxx(j)
            sum = 0.d0
            kl = 0
            !
            ! For each atom in second lmo do ...
            do kk = nncx(j) + 1, nncx(j) + ncx(j)
              k1 = icxxx(kk)
              !
              ! If atom in first lmo, add contribution to C * (F * C)
              if (.not. latoms(k1)) then
                kl = kl + iorbs(k1)
              else
                do k = nfirst(k1), nlast(k1)
                  kl = kl + 1
                  sum = sum + ws(k) * cxxx(kl+loopj)
                end do
              end if
            end do
            !
            ! Store final integral
            fmo_expanded (k_fmo) = sum
          end if
        end if
        k_fmo = k_fmo + 1
    end do
  end do
end subroutine FlmoFromFao

!-------------------------------------------------------------------------------

subroutine eigen_limits (print_nocc, print_nvir)
   !
   ! This routine evaluates the number of occupied and virtual orbitals
   ! to be written out to the output file.
    use molkst_C, only: norbs, nelecs, keywrd
    implicit none
    integer, intent (inout) :: print_nocc, print_nvir
    integer :: i, j
    integer :: noccupied, nvirtual
    integer, parameter :: default_print_nocc = 8
    integer, parameter :: default_print_nvir = 8
    character (len=248) :: tmpkey
    double precision, external :: reada
    noccupied = nelecs/2
    nvirtual = norbs - noccupied
    if (Index (keywrd, " ALLVEC") /= 0) then
      !
      ! Everything requested
      print_nocc = noccupied
      print_nvir = nvirtual
    else if ((Index (keywrd, " VECTORS(") /= 0) .or. &
           & (Index (keywrd, " VECTORS=(") /= 0)) then
      !
      ! User has requested limits
      !
      ! Extract user limits from keyword string
      !
      ! Build temporary keyword string, and erase all text from it
      ! except VECTORS data
      tmpkey = trim(keywrd)
      i = Index (tmpkey, "VECTORS=(") + Index (tmpkey, "VECTORS(")
      tmpkey (:i) = " "
      j = Index (tmpkey, ")")
      tmpkey (j:) = " "
      ! Extract number of occupied orbitals to be written out
      print_nocc = Nint (reada (tmpkey, i))
      ! Extract number of virtual orbitals to be written out
      i = Index (tmpkey, ",") + Index (tmpkey, ":")
      if (i /= 0) then
        print_nvir = Nint (reada (tmpkey, i + 1))
      else
        print_nvir = default_print_nvir
      end if
    else
      !
      ! Use defaults
      print_nocc = default_print_nocc
      print_nvir = default_print_nvir
    end if
   !
   ! Make sure number of orbitals requested is no larger than
   ! the number available
    print_nocc = Max (0, Min (print_nocc, noccupied))
    print_nvir = Max (0, Min (print_nvir, nvirtual))
end subroutine eigen_limits
subroutine lmo_to_eigenvectors(nlmo, ncx, nncx, ncxxx, n_dim, &
     & icxxx, icxxx_dim, cxxx, cxxx_dim, eigs, c)
  use molkst_C, only : numat, norbs
  use common_arrays_C, only : nfirst, nlast
  implicit none
  integer, intent (in) :: nlmo, icxxx_dim, cxxx_dim, n_dim
  integer, dimension (n_dim), intent (in) :: ncx, ncxxx, nncx
  integer, dimension (icxxx_dim), intent (in) :: icxxx
  double precision, dimension (norbs), intent (out) :: eigs
  double precision, dimension (norbs, nlmo), intent (out) :: c
  double precision, dimension (cxxx_dim), intent (in) :: cxxx
!
  integer :: i, ka, jj, j, k, l
  double precision, dimension(:), allocatable :: fmo_expanded, &
  vec_lmo, c_lmo, ws, ac
  logical, allocatable, dimension (:) :: latoms
  allocate (fmo_expanded (nlmo*nlmo), vec_lmo(nlmo*nlmo), c_lmo(norbs), &
         & latoms(numat), ws(norbs), &
         & ac(numat), stat=i)

  if (i /= 0) then
    call memory_error ("DiagonaliseFockLMOScheme")
    return
  end if
  call FlmoFromFao (nlmo, ncx, nncx, ncxxx, nlmo, &
           & latoms, ws, icxxx, icxxx_dim, cxxx, cxxx_dim, &
           & ac, fmo_expanded)
!
! Diagonalise Fock matrix
   call rsp(fmo_expanded, nlmo, eigs, vec_lmo)
!
! Reconstruct eigenvectors in AO basis
!
   c(:norbs,:nlmo) = 0.d0
  do i = 1, nlmo
    !
    ! Expand i'th LMO into uncompressed format
    c_lmo = 0.d0
    ka = ncxxx(i)
    do jj = nncx(i) + 1, nncx(i) + ncx(i)
      j = icxxx(jj)
      do k = nfirst(j), nlast(j)
        ka = ka + 1
        c_lmo(k) = cxxx(ka)
      end do
    end do
    !
    ! Add contribution from i'th LMO to all eigenstates
    !
    do k = 1, norbs
      if (Abs (c_lmo(k)) > 1.0d-10) then
        l = i
        do j = 1, nlmo
          c(k,j) = c(k,j) + vec_lmo(l) * c_lmo(k)
          l = l + nlmo
        end do
      end if
    end do
  end do
  deallocate (fmo_expanded, vec_lmo, c_lmo, latoms, ws, ac)
  end subroutine lmo_to_eigenvectors
