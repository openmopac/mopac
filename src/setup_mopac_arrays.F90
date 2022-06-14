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

  subroutine setup_mopac_arrays(n, mode)
!
!   Create and destroy allocatable arrays used by MOPAC
!
!  n == 0  destroy all arrays that have been created.
!  n >  0
!  mode == 1 Create the arrays essential for reading in data
!  mode /= 1 Create MOPAC arrays (holds electronic info.)
!
  use common_arrays_C, only : geo, coord, na, nb, nc, simbol, atmass, &
     & labels, loc, xparam, nat, nfirst, nlast, uspd, pdiag, xparef, &
     & h, w, p, pa, pb, f, c, eigs, dxyz, grad, eigb, txtatm, fb, cb, &
     & wk, errfn, aicorr, nbonds, ibonds, na_store, q, nw, lopt, &
     & hesinv, gnext1, gmin1, ifact, i1fact, ptot2, geoa, l_atom, coorda, &
       txtatm1, workmat1, workmat2, workmat3
!
  USE maps_C, only : react
!
  use iter_C, only : pold, pbold, pold2, pbold2, pold3, pbold3
!
  use symmetry_C, only : depmul, jelem, locpar, idepfn, locdep, jndex, &
  namo, ipo
!
  use molkst_C, only : mpack, n2elec, norbs, numat, nvar, uhf, &
  l123, natoms, mozyme, npulay
!
  USE chanel_C, only : iw
!
  use derivs_C, only : b, ab, fb_ci => fb, fmat, hmat, wmat, aidref, work2
!
  USE ef_C, ONLY: hess, ef_bmat => bmat, u, oldhss, oldu, pmat, uc, hessc, oldf, &
      & d, vmode
!
  use meci_C, only : microa, microb, &
      & cdiag, occa, nalmat, conf, &
      & vectci, ispin, eig, ispqr, deltap
!
  use drc_C, only: vref, vref0, allxyz, allvel, xyz3, vel3, allgeo, geo3, parref
!
  use cosmo_C, only : ipiden, idenat, gden, qdenet, &
  phinet, qscnet, qscat, bmat, abcmat, srad, iatsp, &
  nn, cosurf, xsp, nset, bh, qden, nsetf, isude, sude, &
  arat, amat, cmat
!
  use esp_C, only : cen, ex, iam, cc, &
      & cespm, espi, potpt, es, rnai, rnai2, pf0, pf1, pf2, &
      & pexs, pce, pexpn, ptd, pewcx, pewcy, pewcz, ird, indc,  &
      cequiv, exsr, itemp, dx_array, dy_array, dz_array, td, esp_array, qsc, &
      a2, ovl, cespm2, cespml, cesp, co, al, rad, b_esp, qesp, ind, temp, f0, &
      f1, ce, u_esp, exs, expn, ewcx, ewcy, ewcz, rnai1, fc
!
  use to_screen_C, only : cnorml
!
   use reimers_C, only : natm,nbf,ibf,iat,nbt,&
                 nprn,itrmet,&
                 natt,nstr,&
                 istr,nsym,ivv,npsn,mofrag,&
                 nfocc,nfvir,icifrag,imp2d,imp2s,&
                 isc,krefd,vv,iwk,&
                 istate,nseig,icif,jcif,ncif,&
                 nbtmo,iconf,nspn, &
                 vca,vcb,occfr,&
                 x,y,z,zcore,gamma, beta,betao,&
                 s,r,ppg,pg,dd,ff,aa,&
                 dm,xz,ci,evalmo,evalci,&
                 spintr,cimatr,dmci,&
                 stwt,e0,wk1,wk2,tr1,&
                 wk3,occ,qgs,dipgs,&
                 qbfcore,qbf,qcore,&
                 dip,dtmp, &
                 eec,ee2, &
                 aocc,bocc,aor1,bor1,&
                 aor,bor,aos,bos,aod,bod,&
                 ao1,bo1,cc0,wk0,q0
!
  implicit none
    integer :: n, mode, i, j
    double precision, external :: meci
    if (n > 0) then
      if (mode == 1) then
!
! Create essential arrays only at this point
!
        allocate(geo(3,n), coord(3,n), coorda(3,n), na(n), nb(n), nc(n))
        allocate(simbol(n*3), atmass(n), labels(n), loc(2, 3*n))
        allocate(xparam(3*n), nat(n), nfirst(n), nlast(n), uspd(9*n))
        allocate(pdiag(9*n), xparef(3*n), depmul(3*n), jelem(20, n))
        allocate(txtatm1(n), txtatm(n), locpar(3*n), idepfn(3*n), locdep(3*n))
        allocate(nbonds(n), ibonds(15,n), na_store(n), l_atom(n))
        na = 0
        na_store = 0
        nfirst = -9999
        nlast = -9999
        nbonds = 0
!
      else
!
!  Create main arrays - at this point, the array sizes are all known
!
        j = 0
        if (allocated(grad))     deallocate(grad)
        allocate(grad(nvar), stat=i)
          j = j + i
        if ( .not. mozyme) then
!
!  A conventional MOPAC job
!
          if (allocated(h))        deallocate(h)
          if (allocated(p))        deallocate(p)
          if (allocated(pa))       deallocate(pa)
          if (allocated(pb))       deallocate(pb)
          if (allocated(pold))     deallocate(pold)
          if (allocated(pold2))    deallocate(pold2)
          if (allocated(pold3))    deallocate(pold3)
          if (allocated(pbold))    deallocate(pbold)
          if (allocated(pbold2))   deallocate(pbold2)
          if (allocated(pbold3))   deallocate(pbold3)
          if (allocated(errfn))    deallocate(errfn)
          if (allocated(aicorr))   deallocate(aicorr)
          if (allocated(f))        deallocate(f)
          if (allocated(fb))       deallocate(fb)
          if (allocated(w))        deallocate(w)
          if (allocated(wk))       deallocate(wk)
          if (allocated(dxyz))     deallocate(dxyz)
          if (allocated(c))        deallocate(c)
          if (allocated(cb))       deallocate(cb)
          if (allocated(eigs))     deallocate(eigs)
          if (allocated(q))        deallocate(q)
          if (allocated(eigb))     deallocate(eigb)
          if (allocated(workmat1)) deallocate(workmat1)
          if (allocated(workmat2)) deallocate(workmat2)
          if (allocated(workmat3)) deallocate(workmat3)
          allocate(h(mpack), p(mpack), pa(mpack), pb(mpack), stat=i)
          j = j + i
          allocate(pold(npulay*mpack),  pold2(npulay*mpack), f(mpack), stat=i)
          j = j + i
          allocate(c(norbs, norbs), eigs(norbs + 1), q(numat), stat=i)
          j = j + i
          allocate(workmat1(norbs, norbs), workmat2(norbs, norbs), workmat3(norbs, norbs), stat=i)
          j = j + i
          allocate(eigb(norbs + 1), pold3(max(mpack, 400)), stat=i)
          j = j + i
          allocate(errfn(3*natoms*l123), aicorr(nvar), stat=i)
          j = j + i
          allocate(w(n2elec + 2025), stat=i)
          j = j + i
          allocate(dxyz(3*numat*l123), stat=i)
          j = j + i
          if (uhf) allocate(fb(mpack), cb(norbs, norbs), pbold(npulay*mpack), &
          pbold2(npulay*mpack), pbold3(max(mpack, npulay*npulay)), stat=i)
          j = j + i
          if (j /= 0) then
            call mopend("A problem occurred during memory assignment, most likely the system is too big to run. ")
            return
          end if
          p = 0.d0
          pold = 0.d0
          eigs = 0.d0
          eigb = 0.d0
          eigs = 0.d0
          if (uhf) pbold = 0.d0
          if (l123 > 1) then
            allocate(wk(n2elec + 2025), stat = i)
            if (i /= 0) then
              call mopend("Failed to allocate the two-electron integral array 'wk'")
              write(iw,"(10x,a,f9.2,a)") "  Size requested:",(n2elec*8.d0)/1.d6, "Mb"
              return
            end if
          end if
          dxyz = 0.d0
          errfn = 0.d0
          aicorr = 0.d0
          grad = 0.d0
        end if
      end if
    else
!
!  Delete all arrays
!
    mpack = 0
    numat = 0
    n2elec =0
    if ( meci() < - 1.d0) n2elec = 0
    if (allocated(f)) call fock2(f, f, f, w, w, w, numat, nfirst, nlast, mode)
    if (allocated(a2))         deallocate(a2, stat = i)
    if (allocated(ab))         deallocate(ab, stat = i)
    if (allocated(abcmat))     deallocate (abcmat, stat = i)
    if (allocated(aicorr))     deallocate (aicorr, stat = i)
    if (allocated(aidref))     deallocate(aidref, stat = i)
    if (allocated(al))         deallocate(al, stat = i)
    if (allocated(allgeo))     deallocate (allgeo, stat = i)
    if (allocated(allvel))     deallocate (allvel, stat = i)
    if (allocated(allxyz))     deallocate (allxyz, stat = i)
    if (allocated(amat))       deallocate (amat, stat = i)
    if (allocated(arat))       deallocate (arat, stat = i)
    if (allocated(atmass))     deallocate (atmass, stat = i)
    if (allocated(b))          deallocate(b, stat = i)
    if (allocated(b_esp))      deallocate(b_esp, stat = i)
    if (allocated(bh))         deallocate (bh, stat = i)
    if (allocated(bmat))       deallocate (bmat, stat = i)
    if (allocated(c))          deallocate (c, stat = i)
    if (allocated(cb))         deallocate (cb, stat = i)
    if (allocated(cc))         deallocate(cc, stat = i)
    if (allocated(ce))         deallocate(ce, stat = i)
    if (allocated(cen))        deallocate(cen, stat = i)
    if (allocated(cequiv))     deallocate(cequiv, stat = i)
    if (allocated(cesp))       deallocate(cesp, stat = i)
    if (allocated(cespm))      deallocate(cespm, stat = i)
    if (allocated(cespm2))     deallocate(cespm2, stat = i)
    if (allocated(cespml))     deallocate(cespml, stat = i)
    if (allocated(cmat))       deallocate (cmat, stat = i)
    if (allocated(cnorml))     deallocate(cnorml, stat = i)
    if (allocated(co))         deallocate(co, stat = i)
    if (allocated(coord))      deallocate (coord, stat = i)
    if (allocated(coorda))     deallocate (coorda, stat = i)
    if (allocated(cosurf))     deallocate (cosurf, stat = i)
    if (allocated(d))          deallocate (d, stat = i)
    if (allocated(deltap))     deallocate (deltap, stat = i)
    if (allocated(depmul))     deallocate (depmul, stat = i)
    if (allocated(dx_array))   deallocate(dx_array, stat = i)
    if (allocated(dxyz))       deallocate (dxyz, stat = i)
    if (allocated(dy_array))   deallocate(dy_array, stat = i)
    if (allocated(dz_array))   deallocate(dz_array, stat = i)
    if (allocated(ef_bmat))    deallocate (ef_bmat, stat = i)
    if (allocated(eigb))       deallocate (eigb, stat = i)
    if (allocated(eigs))       deallocate (eigs, stat = i)
    if (allocated(errfn))      deallocate (errfn, stat = i)
    if (allocated(es))         deallocate(es, stat = i)
    if (allocated(esp_array))  deallocate(esp_array, stat = i)
    if (allocated(espi))       deallocate(espi, stat = i)
    if (allocated(ewcx))       deallocate(ewcx, stat = i)
    if (allocated(ewcy))       deallocate(ewcy, stat = i)
    if (allocated(ewcz))       deallocate(ewcz, stat = i)
    if (allocated(ex))         deallocate(ex, stat = i)
    if (allocated(expn))       deallocate(expn, stat = i)
    if (allocated(exs))        deallocate(exs, stat = i)
    if (allocated(exsr))       deallocate(exsr, stat = i)
    if (allocated(f))          deallocate (f, stat = i)
    if (allocated(f0))         deallocate(f0, stat = i)
    if (allocated(f1))         deallocate(f1, stat = i)
    if (allocated(fb))         deallocate (fb, stat = i)
    if (allocated(fb_ci))      deallocate(fb_ci, stat = i)
    if (allocated(fc))         deallocate(fc, stat = i)
    if (allocated(fmat))       deallocate(fmat, stat = i)
    if (allocated(gden))       deallocate (gden, stat = i)
    if (allocated(geo))        deallocate (geo, stat = i)
    if (allocated(geo3))       deallocate (geo3, stat = i)
    if (allocated(geoa))       deallocate(geoa, stat = i)
    if (allocated(gmin1))      deallocate(gmin1, stat = i)
    if (allocated(gnext1))     deallocate(gnext1, stat = i)
    if (allocated(grad))       deallocate (grad, stat = i)
    if (allocated(h))          deallocate (h, stat = i)
    if (allocated(hesinv))     deallocate(hesinv, stat = i)
    if (allocated(hess))       deallocate (hess, stat = i)
    if (allocated(hessc))      deallocate (hessc, stat = i)
    if (allocated(hmat))       deallocate(hmat, stat = i)
    if (allocated(i1fact))     deallocate(i1fact, stat = i)
    if (allocated(iam))        deallocate(iam, stat = i)
    if (allocated(iatsp))      deallocate (iatsp, stat = i)
    if (allocated(ibonds))     deallocate (ibonds, stat = i)
    if (allocated(idenat))     deallocate (idenat, stat = i)
    if (allocated(idepfn))     deallocate (idepfn, stat = i)
    if (allocated(ifact))      deallocate(ifact, stat = i)
    if (allocated(ind))        deallocate(ind, stat = i)
    if (allocated(indc))       deallocate(indc, stat = i)
    if (allocated(ipiden))     deallocate (ipiden, stat = i)
    if (allocated(ipo))        deallocate(ipo, stat = i)
    if (allocated(ird))        deallocate(ird, stat = i)
    if (allocated(isude))      deallocate (isude, stat = i)
    if (allocated(itemp))      deallocate(itemp, stat = i)
    if (allocated(jelem))      deallocate (jelem, stat = i)
    if (allocated(jndex))      deallocate (jndex, stat = i)
    if (allocated(labels))     deallocate (labels, stat = i)
    if (allocated(loc))        deallocate (loc, stat = i)
    if (allocated(locdep))     deallocate (locdep, stat = i)
    if (allocated(locpar))     deallocate (locpar, stat = i)
    if (allocated(lopt))       deallocate (lopt, stat = i)
    if (allocated(na))         deallocate (na, stat = i)
    if (allocated(na_store))   deallocate (na_store, stat = i)
    if (allocated(l_atom))     deallocate (l_atom, stat = i)
    if (allocated(namo))       deallocate (namo, stat = i)
    if (allocated(nat))        deallocate (nat, stat = i)
    if (allocated(nb))         deallocate (nb, stat = i)
    if (allocated(nbonds))     deallocate (nbonds, stat = i)
    if (allocated(nc))         deallocate (nc, stat = i)
    if (allocated(nfirst))     deallocate (nfirst, stat = i)
    if (allocated(nlast))      deallocate (nlast, stat = i)
    if (allocated(nn))         deallocate (nn, stat = i)
    if (allocated(nset))       deallocate (nset, stat = i)
    if (allocated(nsetf))      deallocate (nsetf, stat = i)
    if (allocated(nw))         deallocate (nw, stat = i)
    if (allocated(oldf))       deallocate (oldf, stat = i)
    if (allocated(oldhss))     deallocate (oldhss, stat = i)
    if (allocated(oldu))       deallocate (oldu, stat = i)
    if (allocated(ovl))        deallocate(ovl, stat = i)
    if (allocated(p))          deallocate (p, stat = i)
    if (allocated(pa))         deallocate (pa, stat = i)
    if (allocated(parref))     deallocate (parref, stat = i)
    if (allocated(pb))         deallocate (pb, stat = i)
    if (allocated(pbold))      deallocate (pbold, stat = i)
    if (allocated(pbold2))     deallocate (pbold2, stat = i)
    if (allocated(pbold3))     deallocate (pbold3, stat = i)
    if (allocated(pbold3))     deallocate (pbold3, stat = i)
    if (allocated(pbold3))     deallocate (pbold3, stat = i)
    if (allocated(pce))        deallocate(pce, stat = i)
    if (allocated(pdiag))      deallocate (pdiag, stat = i)
    if (allocated(pewcx))      deallocate(pewcx, stat = i)
    if (allocated(pewcy))      deallocate(pewcy, stat = i)
    if (allocated(pewcz))      deallocate(pewcz, stat = i)
    if (allocated(pexpn))      deallocate(pexpn, stat = i)
    if (allocated(pexs))       deallocate(pexs, stat = i)
    if (allocated(pf0))        deallocate(pf0, stat = i)
    if (allocated(pf1))        deallocate(pf1, stat = i)
    if (allocated(pf2))        deallocate(pf2, stat = i)
    if (allocated(phinet))     deallocate (phinet, stat = i)
    if (allocated(pmat))       deallocate (pmat, stat = i)
    if (allocated(pold))       deallocate (pold, stat = i)
    if (allocated(pold2))      deallocate (pold2, stat = i)
    if (allocated(pold3))      deallocate (pold3,  stat = i)
    if (allocated(potpt))      deallocate(potpt, stat = i)
    if (allocated(ptd))        deallocate(ptd, stat = i)
    if (allocated(ptot2))      deallocate(ptot2, stat = i)
    if (allocated(q))          deallocate (q, stat = i)
    if (allocated(qden))       deallocate (qden, stat = i)
    if (allocated(qdenet))     deallocate (qdenet, stat = i)
    if (allocated(qesp))       deallocate(qesp, stat = i)
    if (allocated(qsc))        deallocate(qsc, stat = i)
    if (allocated(qscat))      deallocate (qscat, stat = i)
    if (allocated(qscnet))     deallocate (qscnet, stat = i)
    if (allocated(rad))        deallocate(rad, stat = i)
    if (allocated(react))      deallocate (react, stat = i)
    if (allocated(rnai))       deallocate(rnai, stat = i)
    if (allocated(rnai1))      deallocate(rnai1, stat = i)
    if (allocated(rnai2))      deallocate(rnai2, stat = i)
    if (allocated(simbol))     deallocate (simbol, stat = i)
    if (allocated(srad))       deallocate (srad, stat = i)
    if (allocated(sude))       deallocate (sude, stat = i)
    if (allocated(td))         deallocate(td, stat = i)
    if (allocated(temp))       deallocate(temp, stat = i)
    if (allocated(txtatm))     deallocate (txtatm, stat = i)
    if (allocated(txtatm1))    deallocate (txtatm1, stat = i)
    if (allocated(u))          deallocate (u, stat = i)
    if (allocated(u_esp))      deallocate(u_esp, stat = i)
    if (allocated(uc))         deallocate (uc, stat = i)
    if (allocated(uspd))       deallocate (uspd, stat = i)
    if (allocated(vel3))       deallocate (vel3, stat = i)
    if (allocated(vmode))      deallocate (vmode, stat = i)
    if (allocated(vref))       deallocate (vref, stat = i)
    if (allocated(vref0))      deallocate (vref0, stat = i)
    if (allocated(w))          deallocate (w, stat = i)
    if (allocated(wk))         deallocate (wk, stat = i) ! stat needed because the array can be faulty
    if (allocated(wmat))       deallocate(wmat, stat = i)
    if (allocated(work2))      deallocate(work2, stat = i)
    if (allocated(xparam))     deallocate (xparam, stat = i)
    if (allocated(xparef))     deallocate (xparef, stat = i)
    if (allocated(xsp))        deallocate (xsp, stat = i)
    if (allocated(xyz3))       deallocate (xyz3, stat = i)
    if (allocated(occa))       deallocate(occa, stat = i)
    if (allocated(microa))     deallocate(microa, stat = i)
    if (allocated(microb))     deallocate(microb, stat = i)
    if (allocated(cdiag))      deallocate(cdiag, stat = i)
    if (allocated(nalmat))     deallocate(nalmat, stat = i)
    if (allocated(eig))        deallocate(eig, stat = i)
    if (allocated(conf))       deallocate(conf, stat = i)
    if (allocated(vectci))     deallocate(vectci, stat = i)
    if (allocated(ispin))      deallocate(ispin, stat = i)
    if (allocated(ispqr))      deallocate(ispqr, stat = i)
    if (allocated(natm))       deallocate(natm, stat = i)
    if (allocated(nbf))        deallocate(nbf, stat = i)
    if (allocated(ibf))        deallocate(ibf, stat = i)
    if (allocated(iat))        deallocate(iat, stat = i)
    if (allocated(nbt))        deallocate(nbt, stat = i)
    if (allocated(nprn))       deallocate(nprn, stat = i)
    if (allocated(itrmet))     deallocate(itrmet, stat = i)
    if (allocated(natt))       deallocate(natt, stat = i)
    if (allocated(nstr))       deallocate(nstr, stat = i)
    if (allocated(istr))       deallocate(istr, stat = i)
    if (allocated(nsym))       deallocate(nsym, stat = i)
    if (allocated(ivv))        deallocate(ivv, stat = i)
    if (allocated(npsn))       deallocate(nspn, stat = i)
    if (allocated(mofrag))     deallocate(mofrag, stat = i)
    if (allocated(nfocc))      deallocate(nfocc, stat = i)
    if (allocated(nfvir))      deallocate(nfvir, stat = i)
    if (allocated(icifrag))    deallocate(icifrag, stat = i)
    if (allocated(imp2d))      deallocate(imp2d, stat = i)
    if (allocated(imp2s))      deallocate(imp2s, stat = i)
    if (allocated(isc))        deallocate(isc, stat = i)
    if (allocated(krefd))      deallocate(krefd, stat = i)
    if (allocated(vv))         deallocate(vv, stat = i)
    if (allocated(iwk))        deallocate(iwk, stat = i)
    if (allocated(istate))     deallocate(istate, stat = i)
    if (allocated(nseig))      deallocate(nseig, stat = i)
    if (allocated(icif))       deallocate(icif, stat = i)
    if (allocated(jcif))       deallocate(jcif, stat = i)
    if (allocated(ncif))       deallocate(ncif, stat = i)
    if (allocated(nbtmo))      deallocate(nbtmo, stat = i)
    if (allocated(iconf))      deallocate(iconf, stat = i)
    if (allocated(nspn))       deallocate(nspn, stat = i)
    if (allocated(vca))        deallocate(vca, stat = i)
    if (allocated(vcb))        deallocate(vcb, stat = i)
    if (allocated(occfr))      deallocate(occfr, stat = i)
    if (allocated(x))          deallocate(x, stat = i)
    if (allocated(y))          deallocate(y, stat = i)
    if (allocated(z))          deallocate(z, stat = i)
    if (allocated(zcore))      deallocate(zcore, stat = i)
    if (allocated(gamma))      deallocate(gamma, stat = i)
    if (allocated(beta))       deallocate(beta, stat = i)
    if (allocated(betao))      deallocate(betao, stat = i)
    if (allocated(s))          deallocate(s, stat = i)
    if (allocated(r))          deallocate(r, stat = i)
    if (allocated(ppg))        deallocate(ppg, stat = i)
    if (allocated(pg))         deallocate(pg, stat = i)
    if (allocated(dd))         deallocate(dd, stat = i)
    if (allocated(ff))         deallocate(ff, stat = i)
    if (allocated(aa))         deallocate(aa, stat = i)
    if (allocated(dm))         deallocate(dm, stat = i)
    if (allocated(xz))         deallocate(xz, stat = i)
    if (allocated(ci))         deallocate(ci, stat = i)
    if (allocated(evalmo))     deallocate(evalmo, stat = i)
    if (allocated(evalci))     deallocate(evalci, stat = i)
    if (allocated(spintr))     deallocate(spintr, stat = i)
    if (allocated(cimatr))     deallocate(cimatr, stat = i)
    if (allocated(dmci))       deallocate(dmci, stat = i)
    if (allocated(stwt))       deallocate(stwt, stat = i)
    if (allocated(e0))         deallocate(e0, stat = i)
    if (allocated(wk1))        deallocate(wk1, stat = i)
    if (allocated(wk2))        deallocate(wk2, stat = i)
    if (allocated(tr1))        deallocate(tr1, stat = i)
    if (allocated(wk3))        deallocate(wk3, stat = i)
    if (allocated(occ))        deallocate(occ, stat = i)
    if (allocated(qgs))        deallocate(qgs, stat = i)
    if (allocated(dipgs))      deallocate(dipgs, stat = i)
    if (allocated(qbfcore))    deallocate(qbfcore, stat = i)
    if (allocated(qbf))        deallocate(qbf, stat = i)
    if (allocated(qcore))      deallocate(qcore, stat = i)
    if (allocated(dip))        deallocate(dip, stat = i)
    if (allocated(dtmp))       deallocate(dtmp, stat = i)
    if (allocated(eec))        deallocate(eec, stat = i)
    if (allocated(ee2))        deallocate(ee2, stat = i)
    if (allocated(aocc))       deallocate(aocc, stat = i)
    if (allocated(bocc))       deallocate(bocc, stat = i)
    if (allocated(aor1))       deallocate(aor1, stat = i)
    if (allocated(bor1))       deallocate(bor1, stat = i)
    if (allocated(aor))        deallocate(aor, stat = i)
    if (allocated(bor))        deallocate(bor, stat = i)
    if (allocated(aos))        deallocate(aos, stat = i)
    if (allocated(bos))        deallocate(bos, stat = i)
    if (allocated(aod))        deallocate(aod, stat = i)
    if (allocated(bod))        deallocate(bod, stat = i)
    if (allocated(ao1))        deallocate(ao1, stat = i)
    if (allocated(bo1))        deallocate(bo1, stat = i)
    if (allocated(cc0))        deallocate(cc0, stat = i)
    if (allocated(wk0))        deallocate(wk0, stat = i)
    if (allocated(q0))         deallocate(q0, stat = i)
    if (allocated(workmat1))   deallocate(workmat1, stat = i)
    if (allocated(workmat2))   deallocate(workmat2, stat = i)
    if (allocated(workmat3))   deallocate(workmat3, stat = i)
    end if
  end subroutine setup_mopac_arrays
  subroutine memory_error(txt)
    USE chanel_C, only : iw
    character (len=*) :: txt
    write(iw,'(/10x,a,/)')"Unable to allocate memory in subroutine "//txt(:len_trim(txt))
    call mopend(txt(:len_trim(txt)))
    return
  end subroutine memory_error
