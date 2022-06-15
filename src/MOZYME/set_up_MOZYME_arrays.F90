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

  subroutine set_up_MOZYME_arrays()
    use MOZYME_C, only : jopt, ncocc, ncvir, nncf, nnce, icocc, icvir, kopt, &
    noccupied, nvirtual, icocc_dim, icvir_dim, ncf, nce, cocc, cvir, &
    ipad2,  ipad4, cocc_dim, cvir_dim, numred, norred, nelred, &
    fmo_dim, part_dxyz, p1, p2, p3, ws, partf, parth, &
    partp, idiag, nfmo, fmo, ifmo, iorbs, mode, rapid
!
    use molkst_C, only : n2elec, mpack, numat, nelecs, norbs, l123, &
    keywrd, id
!
    use common_arrays_C, only : h, w, wk, p, f, ifact, q, &
    dxyz, eigs, nfirst, nlast, errfn
!
    use iter_C, only : pold, pold2, pold3, pbold, pbold2
!
    implicit none
    integer :: i, j
!
! Deallocate arrays before use (not essential, but prevents possible errors if a complicated
! job involving several calculations is run)
!
  call delete_MOZYME_arrays()
  allocate (iorbs(numat), stat = i)
    if (i /= 0) then

    end if
    iorbs = nlast - nfirst + 1
!
! The sizes of the arrays used by MOZYME must be evaluated by a call to fillij.
! Once the sizes are known, the arrays can be created.
!



  call fillij(.true.)
  noccupied = nelecs/2
  nvirtual = norbs - noccupied
  ipad2 = min(numat + 10, nint((real(numat))**0.25*50))
  ipad4 = nint(real(ipad2*norbs)/numat)
  icocc_dim = ipad2*noccupied + 1000
  icvir_dim = ipad2*nvirtual + 1000
  cocc_dim = ipad4*noccupied + 1000 !  100 is a safety factor, in case no Lewis structure is possible
  cvir_dim = ipad4*nvirtual + 1000
!
!  The two main sizes are:
!     mpack:  the number of elements in the density - one electron - Fock matrices,
!     n2elec: the number of two-electron integrals that are stored
!


!
!  Allocate the arrays.  If there is not enough memory, it's better to find out before
!  any real work is done.
!
  allocate(h(mpack), w(n2elec + 2025), stat = i)
  j = i
  if (id > 0) then
    allocate( wk(n2elec + 2025), stat = i)
    j = j + i
  end if
  allocate (pold(mpack), pold2(1), pold3(1), pbold(1), pbold2(1), stat = i)
  j = j + i
  allocate (ncocc(noccupied+1), ncvir(nvirtual+1), nncf(noccupied+1), nnce(nvirtual+1), &
               & icocc(icocc_dim), icvir(icvir_dim), jopt(numat), &
               & kopt(numat), ncf(noccupied + 1), nce(nvirtual + 1), p(mpack), &
               f(mpack), q(numat), stat = i)
  j = j + i
  allocate(ifact(norbs), stat = i)
  j = j + i
  allocate(cocc(cocc_dim), cvir(cvir_dim), stat = i)
  j = j + i
  if ( index(keywrd, " RESTART") /= 0 .or. &
       index(keywrd," 1SCF") == 0     .or. &
       index(keywrd, " GRAD") /= 0) then
    if (rapid) then
      allocate(part_dxyz(3,numat*l123), stat = i)
      j = j + i
    end if
    allocate(dxyz(3*numat*l123), errfn(3*numat*l123), stat = i)
    j = j + i
    dxyz = 0.d0
  end if
  allocate(p1(norbs), p2(norbs), p3(norbs), nfmo(norbs), idiag(norbs), stat = i)
  j = j + i
  allocate( ws(norbs), eigs(norbs), stat = i)
  j = j + i
  if (rapid) then
    allocate(  partf(mpack), partp(mpack), parth(mpack), stat = i)
    j = j + i
  else
    allocate(  partf(1), partp(1), parth(1), stat = i)
    mode = 0
  end if
  fmo_dim = Min(300, norbs)*norbs
  allocate (fmo(fmo_dim), ifmo(2,fmo_dim), stat = i)
  if (j /= 0) then
    call mopend("Failed to assign memory in MOZYME")
    return
  end if
  icocc = 0
  icvir = 0
  jopt = 0
  kopt = 0
! Initialize LMOs to zero size
  ncf = 0
  nce = 0
  nncf = 0
  nnce = 0
  ncocc = 0
  ncvir = 0
!
!   Create the ijbo array, if it exists?
!

  call fillij(.false.)

  ifact(1) = 0
  do i = 1, norbs - 1
    ifact(i + 1) = i + ifact(i)
  end do


  numred = numat
  norred = norbs
  nelred = nelecs

  end subroutine set_up_MOZYME_arrays
  subroutine delete_MOZYME_arrays()
    use molkst_C, only : numat
    use MOZYME_C, only : jopt, ncocc, ncvir, nncf, nnce, icocc, icvir, kopt, &
    ncf, nce, cocc, cvir, ions, iopt, nijbo, &
    part_dxyz, p1, p2, p3, ws, partf, parth, &
    partp, idiag, nfmo, fmo, ifmo, iorbs
!
    use common_arrays_C, only : h, w, wk, p, f, ifact, fb, q, c, dxyz, eigs, errfn
!
    use iter_C, only : pold, pold2, pold3, pbold, pbold2
    implicit none
    if (numat == 0) then
      if (allocated(ions))     deallocate(ions)
      if (allocated(iopt))     deallocate(iopt)
    end if

    if (allocated(c))          deallocate(c)
    if (allocated(h))          deallocate(h)
    if (allocated(jopt))       deallocate(jopt)
    if (allocated(wk))         deallocate(wk)
    if (allocated(w))          deallocate(w)
    if (allocated(f))          deallocate(f)
    if (allocated(fb))         deallocate(fb)
    if (allocated(p))          deallocate(p)
    if (allocated(pold))       deallocate(pold)
    if (allocated(pold2))      deallocate(pold2)
    if (allocated(pold3))      deallocate(pold3)
    if (allocated(pbold))      deallocate(pbold)
    if (allocated(pbold2))     deallocate(pbold2)
    if (allocated(ncocc))      deallocate(ncocc)
    if (allocated(ncvir))      deallocate(ncvir)
    if (allocated(nncf))       deallocate(nncf)
    if (allocated(nnce))       deallocate(nnce)
    if (allocated(icocc))      deallocate(icocc)
    if (allocated(icvir))      deallocate(icvir)
    if (allocated(kopt))       deallocate(kopt)
    if (allocated(ncf))        deallocate(ncf)
    if (allocated(nce))        deallocate(nce)
    if (allocated(pbold2))     deallocate(pbold2)
    if (allocated(ifact))      deallocate(ifact)
    if (allocated(cocc))       deallocate(cocc)
    if (allocated(cvir))       deallocate(cvir)
    if (allocated(q))          deallocate(q)
    if (allocated(dxyz))       deallocate(dxyz)
    if (allocated(errfn))      deallocate(errfn)
    if (allocated(pold))       deallocate(pold)
    if (allocated(p1))         deallocate(p1)
    if (allocated(p2))         deallocate(p2)
    if (allocated(p3))         deallocate(p3)
    if (allocated(ws))         deallocate(ws)
    if (allocated(fmo))        deallocate(fmo)
    if (allocated(ifmo))       deallocate(ifmo)
    if (allocated(partf))      deallocate(partf)
    if (allocated(partp))      deallocate(partp)
    if (allocated(parth))      deallocate(parth)
    if (allocated(part_dxyz))  deallocate(part_dxyz)
    if (allocated(nfmo))       deallocate(nfmo)
    if (allocated(eigs))       deallocate(eigs)
    if (allocated(idiag))      deallocate(idiag)
    if (allocated(iorbs))      deallocate(iorbs)
    if (allocated(nijbo))      deallocate(nijbo)
    return
  end subroutine delete_MOZYME_arrays
