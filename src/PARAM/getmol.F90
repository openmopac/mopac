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

subroutine getmol (iloop)
!
    use param_global_C, only : weight,nvars, locs, titls, comnts, keys, &
    ndeps, locpas, idepfs, locdes, depmuls, pas, pbs, nats, nfirss, nlasts, &
    numats, norbss, nlecss, nopens, nlm61, ncloss, fracts, &
    natmss, nas, nbs, ncs, lablss, heats, hions, dipls, Atom_pKas, pkas, &
    geos, names, molnam, refgeo, refhof, refips,  msdels, &
    refdip, wthof, wtips, wtgeo, wtdip, wtz, n2elecs, &
    refpka, wtpka, atom_pKa, atom_pKa_O, geotxt, l123s, nnalpha_open, &
    nnbeta_open


    use molkst_C, only : title, koment, natoms, norbs, mpack, numcal, &
    & ndep, nelecs, nvar, nclose, nopen, nalpha, nbeta, fract, lm61, &
    & keywrd, numat, n2elec,  atheat , uhf, id, msdel, refkey, l1u, &
    & l2u, l3u, l123, mol_weight, rhf, lxfac, nalpha_open, nbeta_open

    use common_arrays_C, only : nfirst, nlast, p, pa, pb, nat, labels, na, &
     nb, nc, geo, uspd, xparam, loc, nw, pdiag, atmass

    use parameters_C, only : ams

    use symmetry_C, only : locpar, idepfn, locdep, depmul
    implicit none
!----------------------------------------------------------------------
    integer, intent (in) :: iloop
    integer :: i, iatm, igeo, iloc, ips, isym, j
    save :: iatm, igeo, iloc, ips, isym
    if (iloop == 1) then
      igeo = 0
      isym = 0
      iloc = 0
      iatm = 0
      ips = 0
    end if
    molnam = names(iloop)(:40)
  !
  !     RESTORE ALL DATA RELATING TO THE MOLECULE.
  !
    numcal = iloop
  !#      MOLNUM=ILOOP
    natoms = natmss(iloop)
    numat = numats(iloop)
    ndep = ndeps(iloop)
    norbs = norbss(iloop)
    n2elec = n2elecs(iloop)
    mpack = (norbs*(norbs+1)) / 2
    nelecs = nlecss(iloop)
    nvar = nvars(iloop)
    keywrd = keys(iloop)
    refkey(1) = trim(keywrd)
    nclose = ncloss(iloop)
    nopen = nopens(iloop)
    nalpha_open = nnalpha_open(iloop)
    nbeta_open = nnbeta_open(iloop)
    msdel = msdels(iloop)
    l1u = l123s(1,iloop)
    l2u = l123s(2,iloop)
    l3u = l123s(3,iloop)
    l123 = (2*l1u + 1)*(2*l2u + 1)*(2*l3u + 1)
    rhf   = (Index (keywrd, " RHF") + Index(keywrd, " OPEN") + &
    index(keywrd, " MECI") /= 0 .or. Index (keywrd, " C.I.") /= 0)
    lxfac = (Index (keywrd, " XFAC") /= 0)
    if(index(keywrd," UHF") /=0)then
    nalpha = nclose
    nbeta  = nopen
    nclose = 0
    nopen  = 0
    uhf = .true.
    else
    nalpha = 0
    nbeta  = 0
    uhf = .false.
    end if
    if (.not. rhf .and. .not. uhf) then
      if (mod(nelecs,2) == 0) then
        rhf = .true.
      else
        nalpha = nclose
        nbeta  = nopen
        nclose = 0
        nopen  = 0
        uhf = .true.
      end if
    end if
    title = trim(titls(iloop))
    koment = comnts(iloop)
 !   nmos = nnmos(iloop)
    lm61 = nlm61(iloop)
    fract = fracts(iloop)
  !
  ! RESTORE GEOMETRIC DATA
  !
    do i = 1, natoms
      igeo = igeo + 1
      do j = 1, 3
        geo(j, i) = geos(j, igeo)
      end do
    end do
    do i = 1, ndep
      isym = isym + 1
      locpar(i) = locpas(isym)
      idepfn(i) = idepfs(isym)
      locdep(i) = locdes(isym)
      depmul(i) = depmuls(isym)
    end do
  !
  ! XPARAM is defined by the geo - it does NOT exist separately
  !
    do i = 1, nvar
      iloc = iloc + 1
      xparam(i) = geo(locs(2, iloc), locs(1, iloc))
      loc(1, i) = locs(1, iloc)
      loc(2, i) = locs(2, iloc)
    end do
    loc(1, nvar+1) = 0
  !
  !  RESTORE ATOMIC DATA.
  !
  mol_weight = 0.d0
    do i = 1, natoms
      iatm = iatm + 1
      na(i) = nas(iatm)
      nb(i) = nbs(iatm)
      nc(i) = ncs(iatm)
      nat(i) = nats(iatm)
      if ( i <= numat) then
        atmass(i) = ams(nat(i))
        mol_weight = mol_weight + atmass(i)
      end if
      nfirst(i) = nfirss(iatm)
      nlast(i) = nlasts(iatm)
      labels(i) = lablss(iatm)
    end do
    j = 1
    do i = 1, numat
      nw(i) = j
      j = j + ((nlast(i)-nfirst(i)+1)*(nlast(i)-nfirst(i)+2))/2
    end do
    id = 0
  !
  !  RESTORE ORBITAL DATA
  !
    call getusp (nat, nfirst, nlast, uspd, atheat)
    do i = 1, mpack
      pa(i) = pas(i+ips)
      pb(i) = pbs(i+ips)
      p(i) = pa(i) + pb(i)
    end do
    do i = 1, norbs
      pdiag(i) = p((i*(i+1))/2)
    end do
    ips = ips + mpack
  !
  !  RESTORE REFERENCE DATA
  !
    refhof = heats(iloop)
    refips = hions(iloop)
    refdip = dipls(iloop)
    refpka = pkas(iloop)
    atom_pKa = Atom_pKas(iloop,1)
    atom_pKa_O = Atom_pKas(iloop,2)
    do i = 1, 300
      refgeo(i) = geotxt(i, iloop)
    end do
    wthof = weight(1, iloop)
    wtdip = weight(2, iloop)
    wtips = weight(3, iloop)
    wtgeo = weight(4, iloop)
    wtz   = weight(5, iloop)
    wtpka = weight(6, iloop)
    if (atom_pKa /= 0) then  !  it's a pKa, so set pKa weight and set HoF weight to zero.
      wtpKa = wthof
      wthof = 0.d0
    else
      wtpKa = 0.d0
    end if
end subroutine getmol
subroutine getusp (nat, nfirst, nlast, uspd, atheat)
    use molkst_C, only : numat, norbs, keywrd, cutofp, id
    USE parameters_C, only : eisol, eheat, uss, upp, udd, eheat_sparkles, zs
    use funcon_C, only : fpc_9
    use common_arrays_C, only : geo, coord
    implicit none
    integer, dimension (numat), intent (in) :: nat, nfirst, nlast
    double precision, dimension (norbs), intent (out) :: uspd
    double precision, intent (out) :: atheat
    integer :: ia, ib, ii, j, k, k1, ni, i
    double precision :: eat
    double precision, external :: C_triple_bond_C
    save :: eat
!----------------------------------------------------------------
    atheat = 0.d0
    if (index(keywrd, " SPARKL") /= 0) then
      atheat = 0.d0
      do i = 1, numat
        if  (nat(i) > 56 .and. nat(i) < 72 .and. zs(nat(i)) < 0.1d0) then
          atheat = atheat + eheat_sparkles(nat(i))
        else
          atheat = atheat + eheat(nat(i))
        end if
      end do
    else
      atheat = sum(eheat(nat(:numat)))
    end if
    eat = sum(eisol(nat(:numat)))
    do ii = 1, numat
      ia = nfirst(ii)
      ib = nlast(ii)
      ni = nat(ii)
      if (ia <= ib) then
        uspd(ia) = uss(ni)
        if (ia /= ib) then
          k = ia + 1
          k1 = ia + 3
          do j = k, k1
            uspd(j) = upp(ni)
          end do
          if (k1 /= ib) then
            k = k1 + 1
            do j = k, ib
              uspd(j) = udd(ni)
            end do
          end if
        end if
      end if
    end do
    atheat = atheat - eat * fpc_9
    call gmetry (geo, coord)
    if (id == 0) then
      cutofp = 1.d10
    else
      cutofp = 30.d0
    end if
    call set_up_dentate
    atheat = atheat + C_triple_bond_C()
!
end subroutine getusp
