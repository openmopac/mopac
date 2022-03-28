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

subroutine makvec ()
!
!   Construct the starting occupied and virtual localized molecular orbitals.
!   These are on one or two atoms.  The LMOs are stored in
!     Occupied   Virtual
!      cocc        cvir    Coefficients (this is a dilute array)
!      ncocc       ncvir   Starting addresses of LMO coefficients (dense)
!      icocc       icvir   Atoms numbers (this is a dilute array)
!      nncf        nnce    Starting addresses of atom numbers (dense)
!      ncf         nce     Number of atoms in LMOs (dense)
!
    use MOZYME_C, only: iorbs, icocc, icvir, nce, ncf, morb, Lewis_elem, &
         & ncocc, ncvir, nnce, nncf, Lewis_tot, &
         & cocc, cvir, partf, cocc_dim, cvir_dim
!
    use common_arrays_C, only : nfirst, nlast, pdiag, f, p
    use molkst_C, only: numat, norbs, keywrd, moperr, mpack
    implicit none
    logical :: lok, times
    integer :: i, ii, j, jj, k, locc, lvir, m, m1, ne, nf_loc, ni, nj, &
         & nocc, nvir, alloc_stat, nLewis
    double precision, dimension(:,:), allocatable :: catom
    integer, dimension(:), allocatable :: iz, ib
    logical, dimension (:), allocatable :: u
    integer, external :: ijbo
    times = (Index (keywrd, " TIMES") /= 0)
    !***************************************************************
    !
    !  Set up starting eigenvector matrix, based on the Fock matrix
    !  and the type of atom
    !
    !***************************************************************
    p = 0.d0
    m = 0
    do i = 1, numat
      j = ijbo (i, i)
      do k = 1, iorbs(i)
        m = m + 1
        j = j + k
        p(j) = pdiag(m)
      end do
    end do
    if (times) call timer (" After entry to Makvec")
    call buildf (f, partf, 0)
    if (times) call timer (" After BUILDF in MAKV")
   !
    allocate (iz(numat), ib(numat), catom(morb, norbs), u(norbs), stat=alloc_stat)
    if (alloc_stat /= 0) then
      call memory_error ("makvec")
      go to 1200
    end if
   !
    moperr = .false.
!
!   Construct hybrid atomic orbitals for each atom. Put these in catom
!
    call hybrid (catom)
    !
    !     Now to make the starting localized molecular orbitals
    !
    u(:norbs) = .false.
    nocc = 0
    nvir = 0
    nf_loc = 0
    ne = 0
    locc = 0
    lvir = 0
    do nLewis = 1, Lewis_tot
      ii = Lewis_elem(1, nLewis)
      jj = Lewis_elem(2, nLewis)
      if (ii > 0 .and. jj > 0) then
!
! Di-atomic Lewis element
!
        ni = nfirst(ii)
        nj = nfirst(jj)
        call mbonds (locc, lvir, f, catom, nfirst, ii, jj, u(ni), u(nj), &
               & lok, iorbs, cocc, cvir, cocc_dim, cvir_dim, numat, norbs, &
               & morb, mpack)
        if (lok) then
          nncf(nocc+1) = nf_loc
          nnce(nvir+1) = ne
          call mlmo (locc, lvir, ii, jj, nf_loc, ne, nocc, nvir, iz, ib, &
                 & nce, ncf, ncocc, ncvir, iorbs, icocc, icvir, cocc, cvir)
        else
         ii = 0 ! Use for error trap
        end if
      else if (ii > 0) then
!
! Lone pair
!
        do k = nlast(ii), nfirst(ii), -1
          if ( .not. u(k)) then
            u(k) = .true.
            m = locc
            do m1 = 1, iorbs(ii)
              m = m + 1
              cocc(m) = catom(m1, k)
            end do
            nncf(nocc+1) = nf_loc
            call mlmo (locc, lvir, ii, 0, nf_loc, ne, nocc, nvir, iz, ib, &
               & nce, ncf, ncocc, ncvir, iorbs, icocc, icvir, cocc, cvir)
            exit
          end if
        end do
      else
!
! Virtual Lone pair
!
        do k = nlast(jj), nfirst(jj), -1
          if ( .not. u(k)) then
            u(k) = .true.
            m = lvir
            do m1 = 1, iorbs(jj)
              m = m + 1
              cvir(m) = catom(m1, k)
            end do
            nnce(nvir+1) = ne
            call mlmo (locc, lvir, 0, jj, nf_loc, ne, nocc, nvir, iz, ib, &
                 & nce, ncf, ncocc, ncvir, iorbs, icocc, icvir, cocc, cvir)
            exit
          end if
        end do
      end if
    end do
 1200 continue
    if (Allocated (iz))   deallocate ( iz, ib, catom, u)
end subroutine makvec
