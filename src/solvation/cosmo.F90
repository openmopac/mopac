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

subroutine addfck (f, p)
    use molkst_C, only : lm61, numat, mpack
    use funcon_C, only : a0, ev
    use cosmo_C, only : ediel, fepsi, nps, &
     amat, bmat, iatsp, nsetf, phinet, qscnet, qdenet, ipiden, gden, qscat
    implicit none
    double precision, dimension(mpack), intent(in) :: p
    double precision, dimension(mpack), intent(inout) :: f
    integer :: i, iat, im, j
    double precision :: fcon, fim, phi, qsc3, s1, s3
    fcon = a0 * ev
    do i = 1, numat
      qscat(i) = 0.d0
    end do
   ! FIRST CALCULATE QDENEL FROM DENSITY MATRIX
    do i = 1, lm61
      qdenet(i, 2) = gden(i) * p(ipiden(i))
      qdenet(i, 3) = qdenet(i, 2) + qdenet(i, 1)
    end do
   !  NOW CALCULATE PHIEL FROM BMAT*QDENEL
    do i = 1, nps
      phi = 0.d0
      do j = 1, lm61
        phi = phi + bmat(j, i) * qdenet(j, 2)
      end do
      phinet(i, 2) = phi
      phinet(i, 3) = phinet(i, 1) + phi
    end do
   !  NOW CALCULATE QSCEL FROM A*QSCEL = -PHIEL
    call coscl2 (amat, nsetf, qscnet(1, 2), phinet(1, 2), nps)
    ediel = 0.d0
    s1 = 0.0d0
    s3 = 0.0d0
    do i = 1, nps
      iat = iatsp(i)
      qscnet(i, 2) = -fepsi * qscnet(i, 2)
      qsc3 = qscnet(i, 1) + qscnet(i, 2)
      qscnet(i, 3) = qsc3
      ediel = ediel + qsc3 * phinet(i, 3)
      s1 = s1 + qscnet(i, 1)
      s3 = s3 + qsc3
      qscat(iat) = qscat(iat) + qsc3
    end do
    ediel = ediel * fcon / 2
   ! NOW ADD BMAT*QSCEL TO FOCK MATRIX
    do i = 1, lm61
      im = ipiden(i)
      fim = 0.d0
      do j = 1, nps
        fim = fim + bmat(i, j) * qscnet(j, 2)
      end do
      f(im) = f(im) - fcon * fim
    end do
end subroutine addfck
subroutine addhcr
    use molkst_C, only : lm61
    use funcon_C, only : a0, ev
    use cosmo_C, only : nps, bmat, qscnet, ipiden
    use common_arrays_C, only : h
    implicit none
    integer :: i, im, j
    double precision :: fcon, him
    fcon = a0 * ev
    do i = 1, lm61
      im = ipiden(i)
      him = 0.d0
      do j = 1, nps
        him = him + bmat(i, j) * qscnet(j, 1)
      end do
      h(im) = h(im) - fcon * him
    end do
end subroutine addhcr
subroutine addnuc
    use molkst_C, only : numat, lm61, enuclr
    use funcon_C, only : a0, ev
    use common_arrays_c, only : nat
    use cosmo_C, only : nps, fepsi, &
    amat, bmat, nsetf, phinet, qscnet, qdenet, idenat
    use parameters_C, only : tore
    implicit none
    integer :: i, j
    double precision :: enclr, fcon, phi
   !
   ! FIRST CALCULATE QDENNUC
    fcon = a0 * ev
    do i = 1, lm61
      qdenet(i, 1) = 0
    end do
    do i = 1, numat
      qdenet(idenat(i), 1) = tore(nat(i))
    end do
   !  NOW CALCULATE PHINUC AS BMAT*QDENUC
    do i = 1, nps
      phi = 0.d0
      do j = 1, lm61
        phi = phi + bmat(j, i) * qdenet(j, 1)
      end do
      phinet(i, 1) = phi
    end do
   !  NOW CALCULATE QSCNUC FROM  AMAT*QSCNUC = PHINUC
    call coscl2 (amat, nsetf, qscnet(1, 1), phinet(1, 1), nps)
   ! SCALE QSCNUC AND CALCULATE INTERATION ENERGY
    enclr = 0.d0
    do i = 1, nps
      qscnet(i, 1) = -fepsi * qscnet(i, 1)
      enclr = enclr + qscnet(i, 1) * phinet(i, 1)
    end do
    enuclr = enuclr + fcon * enclr / 2
end subroutine addnuc
subroutine ansude (ra, rb, d, rs, ara, arb, aar, abr, arad, arbd, rinc)
    use funcon_C, only : pi
   !   THIS PROGROM CALCULATES THE AREA OF TWO INTERSECTING SPHERES
   !   WITH RADII RA AND RB AT A DISTANCE D AND A SOLVENT PROBE
   !   RADIUS RS. THE TWO AREAS ARE CALCULATED SEPARATELY (ARA, ARB).
   !   THE PART OF THE AREA ON THE CLOSURE PART IS AAR AND ABR.
   !   FOR BOTH AREAS ANALYTIC DERIVATIVES WITH RESPECT TO THE DISTANCE
   !   D ARE CALCULATED (ARAD, ARBD).    (WRITTEN BY ANDREAS KLAMT, 9/9/96)
    implicit none
    double precision, intent (in) :: d, ra, rb, rs
    double precision, intent (out) :: aar, abr, ara, arad, arb, arbd, rinc
    double precision :: ca, cad, cb, cbd, fza, fzad, fzb, fzbd, qa, qb, &
   & sa, sad, sb, sbd, ta, tad, tb, tbd, xa, xad, xb, xbd, ya, yad, yb, ybd, &
   & za, zad, zb, zbd
    qa = ra + rs
    qb = rb + rs
    ca = (qa**2 + d**2 - qb**2) / (2.d0*qa*d)
    cb = (qb**2 + d**2 - qa**2) / (2.d0*qb*d)
    sa = Sqrt (1.d0 - ca*ca)
    sb = Sqrt (1.d0 - cb*cb)
    ta = pi * sa
    tb = pi * sb
    fza = (1.d0 - Cos (ta)) / 2
    fzb = (1.d0 - Cos (tb)) / 2
    if (sa < 0 .or. sb < 0) then
      fza = 1.d0
    end if
    if (sa < 0 .or. sb < 0) then
      fzb = 1.d0
    end if
    xa = fzb ** 1 * rs * (ca + cb)
    xb = fza ** 1 * rs * (ca + cb)
    ya = ra * sa - fzb * rb * sb
    yb = rb * sb - fza * ra * sa
    za = Sqrt (xa*xa + ya*ya)
    zb = Sqrt (xb*xb + yb*yb)
    rinc = 0.5d0*(za + zb)/sqrt(rs*rs*(ca + cb)**2 + (ra*sa - rb*sb)**2)
    ara = pi * ra * (2.d0*(1.d0 + ca)*ra + sa*za)
    arb = pi * rb * (2.d0*(1.d0 + cb)*rb + sb*zb)
    aar = pi * ra * (sa*za)
    abr = pi * rb * (sb*zb)
   ! NOW DERIVATIVES
    cad = (qb**2 + d**2 - qa**2) / (2.d0*qa*d*d)
    cbd = (qa**2 + d**2 - qb**2) / (2.d0*qb*d*d)
    sad = -ca * cad / sa
    sbd = -cb * cbd / sb
    tad = pi * sad
    tbd = pi * sbd
    fzad = Sin (ta) * .5d0
    fzbd = Sin (tb) * .5d0
    if (sa < 0 .or. sb < 0) then
      fzad = 0.d0
    end if
    if (sa < 0 .or. sb < 0) then
      fzbd = 0.d0
    end if
    xad = rs * ((ca + cb)*fzbd*tbd + fzb*(cad + cbd))
    xbd = rs * ((ca + cb)*fzad*tad + fza*(cad + cbd))
    yad = ra * sad-fzbd * tbd * rb * sb - fzb * rb * sbd
    ybd = rb * sbd-fzad * tad * ra * sa - fza * ra * sad
    zad = (xa*xad + ya*yad) / za
    zbd = (xb*xbd + yb*ybd) / zb
    arad = pi * ra * (sad*za + sa*zad + 2.d0*ra*cad)
    arbd = pi * rb * (sbd*zb + sb*zbd + 2.d0*rb*cbd)
end subroutine ansude
subroutine ciint (c34, pq34)
    use molkst_C, only : lm61
    use cosmo_C, only : nps, cmat
    implicit none
    double precision, dimension (lm61), intent (in) :: c34
    double precision, dimension (lm61), intent (inout) :: pq34
    integer :: i, i0, j
    i0 = 0
    do i = 1, lm61
      pq34(i) = 0.d0
    end do
    if (nps < 0) return
    do i = 1, lm61
      do j = 1, i - 1
        i0 = i0 + 1
        pq34(j) = pq34(j) + cmat(i0) * c34(i)
        pq34(i) = pq34(i) + cmat(i0) * c34(j)
      end do
      i0 = i0 + 1
      pq34(i) = pq34(i) + cmat(i0) * c34(i)
    end do
end subroutine ciint
subroutine coscav
   !***********************************************************************
   !
   ! THIS ROUTINE CONSTRUCTS OR UPDATES THE SOLVENT - ACCESSIBLE
   ! SURFACE (SAS)
   !
   !***********************************************************************
    use molkst_C, only : numat
    use cosmo_C, only : lenabc, nps, area, isude, sude, &
       cosvol, rsolv, ioldcv, disex2, n0, dirvec, dirsm, srad, &
       cosurf, amat, iatsp, nar_csm, nsetf, phinet, arat
    use common_arrays_c, only : coord, nat
    use funcon_C, only : pi
    use chanel_C, only : iw
    implicit none
    integer :: i, i0, ik, ilipa, inset, ipm, ips, ix, j, jmax, jps, k, &
   & l, nara, narea, nfl1, nfl2, niter, nps0, maxrs
    double precision :: aa, d2, dist, dist1, dist2, dist3, fdiagr, r, ri, &
   & ri2, rj, rr, sininv, sp, spm, x1, x2, x3, x4, dists
    logical, dimension (10000) :: din ! Was 1082, but some calc'ns used more than 1082
    integer, dimension (1082) :: iseg
    double precision, dimension (3) :: xa, xb, xi, xj, xx
    double precision, dimension (:, :), allocatable :: dirtm
    double precision, dimension (:, :, :), allocatable :: finel
    double precision, dimension(:), allocatable :: rdat
    double precision, dimension(:, :), allocatable :: rsc
    double precision, dimension(:, :, :), allocatable :: tm
    integer, dimension(:), allocatable :: isort, ipsrs, nipsrs, nset, nipa, &
         & lipa
    integer, dimension(:, :), allocatable :: nn
    maxrs = 70 * numat
    allocate (rdat(numat), rsc(4, maxrs), tm(3, 3, numat), nn(3, numat), &
         & isort(maxrs), ipsrs(maxrs), nipsrs(lenabc), nset(1082*numat), &
         & nipa(numat), dirtm(3, 1082), finel(4, 500, 2))
    if (.not. allocated(srad)) allocate(srad(numat))
    if (.not. allocated(nar_csm))  allocate(nar_csm(lenabc + 1))
   ! MAKE COORDINATES A BIT ASYMMETRIC IN ORDER TO AVOID
   ! SYMMETRY PROBLEMS WITH CAVITY CONSTRUCTION
    do i = 1, numat
      do j = 1, 3
        coord(j, i) = coord(j, i) + Cos (i*j*.1d0) * 3.0d-9
      end do
    end do
   !
    ilipa = 0
    do i = 1, numat
      ri = srad(i)
      r = ri + rsolv
      rr = r + rsolv
      ri2 = ri * ri
      do ix = 1, 3
        xa(ix) = coord(ix, i)
      end do
      do j = 1, numat
        if (j /= i) then
          dist = 0.d0
          do ix = 1, 3
            dist = dist + (xa(ix) - coord(ix, j)) ** 2
          end do
          if (dist < (rr + srad(j))**2) then
            ilipa = ilipa + 1
          end if
        end if
      end do
    end do
    allocate (lipa(Max(1, ilipa)))
    fdiagr = 2.1d0 * Sqrt (pi)
    inset = 1
    ilipa = 0
    nps = 0
    area = 0.d0
    cosvol = 0.d0
   ! NOW A LARGE LOOP OVER ALL ATOMS STARTS, WHICH MAKES THE SEGMENTATION
   ! ON THE CONVEX (SPHERICAL) PART OF THE CORRESPONDING PART OF THE CAVITY
    do i = 1, numat
      nipa(i) = 0
      ri = srad(i)
      r = ri + rsolv
      rr = r + rsolv
      ri2 = ri * ri
      do ix = 1, 3
        xa(ix) = coord(ix, i)
      end do
      nps0 = nps + 1
      do j = 1, numat
        if (j /= i) then
          dist = 0.d0
          do ix = 1, 3
            dist = dist + (xa(ix) - coord(ix, j)) ** 2
          end do
          if (dist < (rr + srad(j))**2) then
            ilipa = ilipa + 1
            if (ilipa > maxrs) then
              write(iw, '(/10x, a)')"Solvent radius too large - reduce RSOLV"
              write(iw, '(/10x, a, f6.2, a, /)')"(Current value of RSOLV:", rsolv, " Angstroms)"
              call mopend ("Solvent radius too large - reduce RSOLV")
              return
            end if
            nipa(i) = nipa(i) + 1
            lipa(ilipa) = j
          end if
        end if
      end do
      ! SEARCH FOR 3 NEAREST NEIGHBOR ATOMS
      dist1 = 1.d20
      dist2 = 1.d20
      dist3 = 1.d20
      nn(1, i) = 0
      nn(2, i) = 0
      nn(3, i) = 0
      do j = 1, numat
        if (j /= i) then
          dist = 0.d0
          do ix = 1, 3
            dist = dist + (xa(ix) - coord(ix, j)) ** 2
          end do
          if (dist + 0.05d0 < dist3) then
            dist3 = dist
            nn(3, i) = j
          end if
          if (dist3 + 0.05d0 < dist2) then
            dist = dist2
            dist2 = dist3
            dist3 = dist
            nn(3, i) = nn(2, i)
            nn(2, i) = j
          end if
          if (dist2 + 0.05d0 < dist1) then
            dist = dist1
            dist1 = dist2
            dist2 = dist
            nn(2, i) = nn(1, i)
            nn(1, i) = j
          end if
        end if
      end do
      ! BUILD NEW TRANSFORMATION MATRIX
      if (nn(1, i) == 0) then
        tm(1, 1, i) = 1.d0
        tm(1, 2, i) = 0.d0
        tm(1, 3, i) = 0.d0
      else
        dist1 = 0.d0
        do ix = 1, 3
          dist1 = dist1 + (xa(ix) - coord(ix, nn(1, i))) ** 2
        end do
        dist = 1.d0 / Sqrt (dist1)
        tm(1, 1, i) = (coord(1, nn(1, i)) - xa(1)) * dist
        tm(1, 2, i) = (coord(2, nn(1, i)) - xa(2)) * dist
        tm(1, 3, i) = (coord(3, nn(1, i)) - xa(3)) * dist
      end if
      rr = -1.d0
      i_loop: do
        if (nn(2, i) == 0) then
          rr = rr + 1.d0
          dist = sqrt(9.d0*(tm(1, 2, I) + rr)**2 + 4.d0*(tm(1, 1, i) + rr**2)**2 + tm(1, 3, i)**2)
          xx(1) = 3.d0*(tm(1, 2, I) + rr)/dist
          xx(2) = tm(1, 3, I)/dist
          xx(3) = 2.d0*(tm(1, 1, I) + rr**2)/dist
        else
          dist2 = 0.d0
          do ix = 1, 3
              dist2 = dist2 + (xa(ix) - coord(ix, nn(2, i)))**2
          end do
          dist = 1.0d0/sqrt(dist2)
          xx(1) = (coord(1, nn(2, i)) - xa(1))*dist
          xx(2) = (coord(2, nn(2, i)) - xa(2))*dist
          xx(3) = (coord(3, nn(2, i)) - xa(3))*dist
        end if
        sp = xx(1)*tm(1, 1, i) + xx(2)*tm(1, 2, i) + xx(3)*tm(1, 3, i)
        if (sp*sp > 0.99d0) then
          nn(2, i) = nn(3, i)
          nn(3, i) = 0
          dist2 = dist3
          cycle i_loop
        else
          exit i_loop
        end if
      end do i_loop
      sininv = 1.d0 / Sqrt (1.d0 - sp*sp)
      tm(2, 1, i) = (xx(1) - sp*tm(1, 1, i)) * sininv
      tm(2, 2, i) = (xx(2) - sp*tm(1, 2, i)) * sininv
      tm(2, 3, i) = (xx(3) - sp*tm(1, 3, i)) * sininv
      tm(3, 1, i) = tm(1, 2, i) * tm(2, 3, i) - tm(2, 2, i) * tm(1, 3, i)
      tm(3, 2, i) = tm(1, 3, i) * tm(2, 1, i) - tm(2, 3, i) * tm(1, 1, i)
      tm(3, 3, i) = tm(1, 1, i) * tm(2, 2, i) - tm(2, 1, i) * tm(1, 2, i)
      ! TRANSFORM DIRVEC ACCORDING TO TM
      do j = 1, 1082
        xx(1) = dirvec(1, j)
        xx(2) = dirvec(2, j)
        xx(3) = dirvec(3, j)
        do ix = 1, 3
          dirtm(ix, j) = xx(1) * tm(1, ix, i) + xx(2) * tm(2, ix, i) + xx &
         & (3) * tm(3, ix, i)
        end do
      end do
      ! FIND THE POINTS OF THE BASIC GRID ON THE SAS
      narea = 0
      loop: do j = 1, 1082
        din(j) = .false.
        do ix = 1, 3
          xx(ix) = xa(ix) + dirtm(ix, j) * r
        end do
         !          - - - WE NEED ONLY TRY THOSE ATOMS INTERSECTING ATOM I
        do ik = ilipa - nipa(i) + 1, ilipa
          k = lipa(ik)
          dist = 0.d0
          do ix = 1, 3
            dist = dist + (xx(ix) - coord(ix, k)) ** 2
          end do
          dist = Sqrt (dist) - rsolv - srad(k)
          if (dist < 0) cycle loop
        end do
        narea = narea + 1
        cosvol = cosvol + ri2 * dirvec(4, j) * (dirtm(1, j)*xa(1) + dirtm(2, &
       & j)*xa(2) + dirtm(3, j)*xa(3) + ri)
        area = area + ri2 * dirvec(4, j)
        din(j) = .true.
      end do loop
      if (narea /= 0) then
         !
         !  IF HYDROGEN, THEN USE THE SMALLER SET OF POINTS (NORMALLY
         !  12 POINTS), OTHERWISE, USE THE LARGER SET (NORMALLY 42 POINTS)
         !
        i0 = 1
        if (nat(i) == 1) then
          i0 = 2
        end if
        jmax = n0(i0)
        i0 = (i0 - 1) * n0(1)
        do j = 1, jmax
          nps = nps + 1
          if (nps > lenabc) then
            call mopend ("NPS IS GREATER THAN LENABC - USE SMALLER NSPA")
            go to 100
          else
            iatsp(nps) = i
            do ix = 1, 3
              xx(ix) = dirsm(ix, i0 + j)
            end do
            do ix = 1, 3
              cosurf(ix, nps) = xx(1) * tm(1, ix, i) + xx(2) * tm(2, ix, &
             & i) + xx(3) * tm(3, ix, i)
            end do
          end if
        end do
        niter = 0
        do
          niter = niter + 1
          do ips = nps0, nps
            nar_csm(ips) = 0
            phinet(ips, 1) = 0.d0
            phinet(ips, 2) = 0.d0
            phinet(ips, 3) = 0.d0
          end do
          do j = 1, 1082
            if (din(j)) then
              spm = -1.d0
              ipm = 0
              x1 = dirtm(1, j)
              x2 = dirtm(2, j)
              x3 = dirtm(3, j)
              do ips = nps0, nps
                sp = x1 * cosurf(1, ips) + x2 * cosurf(2, ips) + x3 * cosurf &
               & (3, ips)
                if (sp >= spm) then
                  spm = sp * (1.d0 + 1.d-14)
                  ipm = ips
                end if
              end do
              iseg(j) = ipm
              nar_csm(ipm) = nar_csm(ipm) + 1
              do ix = 1, 3
                phinet(ipm, ix) = phinet(ipm, ix) + dirtm(ix, j) * dirvec &
               & (4, j)
              end do
            end if
          end do
          ips = nps0 - 1
          loop1: do
            ips = ips + 1
            do while (nar_csm(ips) == 0)
              niter = 1
              nps = nps - 1
              if (ips > nps) exit loop1
              do jps = ips, nps
                nar_csm(jps) = nar_csm(jps + 1)
                phinet(jps, 1) = phinet(jps + 1, 1)
                if (Abs(phinet(1, 1)) > 1.d-20) then
                  ips = ips
               end if
                phinet(jps, 2) = phinet(jps + 1, 2)
                phinet(jps, 3) = phinet(jps + 1, 3)
              end do
            end do
            dists = 0.d0
            do ix = 1, 3
              dists = dists + phinet(ips, ix) ** 2
            end do
            dists = Max (dists, 1.d-20)
            dist = 1.d0 / Sqrt (dists)
            do ix = 1, 3
              cosurf(ix, ips) = phinet(ips, ix) * dist
            end do
            if (ips >= nps) exit
          end do loop1
          if (niter >= 2) exit
        end do
         ! NOW ALL SEGMENTS ARE FINALLY DEFINED AND THE ASSOCIATED
         ! BASIC GRID POINTS ARE CLOSE - PACKED
        do ips = nps0, nps
          nsetf(ips) = inset
          inset = inset + nar_csm(ips)
          nar_csm(ips) = 0
          cosurf(4, ips) = 0.d0
          do ix = 1, 3
            cosurf(ix, ips) = cosurf(ix, ips) * ri + xa(ix)
          end do
        end do
        do j = 1, 1082
          if (din(j)) then
            ipm = iseg(j)
            nara = nar_csm(ipm)
            nset(nsetf(ipm) + nara) = j
            nar_csm(ipm) = nara + 1
            cosurf(4, ipm) = cosurf(4, ipm) + dirvec(4, j) * ri2
          end if
        end do
      end if
  !  Check lipa size
      if (ilipa >= 70*numat) then
        exit
      end if
    end do
     !
   ! HERE THE CONSTRUCTION FOR A SINGLE ATOM ENDS
    do i = 1, numat
      din(i) = .true.
    end do
    do j = 1, nps
      din(iatsp(j)) = .false.
    end do
   ! NOW THE SEGMENT FORMATION ON ALL ATOMS IS FINISHED
   ! NOW THE CLOSURE OF THE CONCAVE REGIONS OF THE SURFACE WILL BE DONE
    if (ioldcv == 0) then
      i = max(numat, 1082)
      call surclo (coord, nipa, lipa, din, i, rsc, isort, ipsrs, nipsrs, nat, &
     & srad, cosurf, iatsp, nar_csm, nsetf, isude, sude, maxrs)
    end if
    cosvol = cosvol / 3
    do i = 1, numat
      arat(i) = 0.d0
      rdat(i) = 0.d0
    end do
   ! FILLING AMAT
    do ips = 1, nps
      i = iatsp(ips)
      ri = srad(i)
      do ix = 1, 3
        xi(ix) = coord(ix, i)
        xa(ix) = cosurf(ix, ips)
      end do
      arat(i) = arat(i) + cosurf(4, ips)
      call mfinel (ips, 1, finel, nar_csm, nsetf, nset, rsc, nipsrs, dirvec, &
     & tm(1, 1, i), xi, ri, nfl1, ioldcv, maxrs, lenabc, numat)
      aa = 0.d0
      do k = 1, nfl1
        aa = aa + fdiagr * Sqrt (finel(4, k, 1)**3)
        x1 = finel(1, k, 1)
        x2 = finel(2, k, 1)
        x3 = finel(3, k, 1)
        x4 = finel(4, k, 1)
        do l = 1, k - 1
          aa = aa + 2 * x4 * finel(4, l, 1) / Sqrt ((x1 - finel(1, l, &
         & 1))**2 + (x2 - finel(2, l, 1))**2 + (x3 - finel(3, l, 1))**2)
        end do
      end do
      amat(((ips + 1)*ips)/2) = aa / cosurf(4, ips) ** 2
      rdat(i) = rdat(i) + aa
      do jps = 1, ips - 1
        j = iatsp(jps)
        d2 = 0.d0
        do ix = 1, 3
          xj(ix) = coord(ix, j)
          xb(ix) = cosurf(ix, jps)
          d2 = d2 + (xb(ix) - xa(ix)) ** 2
        end do
        if (d2 > disex2) then
          aa = 1.d0 / Sqrt (d2)
        else
          j = iatsp(jps)
          rj = srad(j)
          call mfinel (jps, 2, finel, nar_csm, nsetf, nset, rsc, nipsrs, &
         & dirvec, tm(1, 1, j), xj, rj, nfl2, ioldcv, maxrs, lenabc, numat)
          aa = 0.d0
          do k = 1, nfl1
            x1 = finel(1, k, 1)
            x2 = finel(2, k, 1)
            x3 = finel(3, k, 1)
            x4 = finel(4, k, 1)
            do l = 1, nfl2
              aa = aa + x4 * finel(4, l, 2) / Sqrt ((x1 - finel(1, l, &
             & 2))**2 + (x2 - finel(2, l, 2))**2 + (x3 - finel(3, l, 2))**2)
            end do
          end do
          aa = aa / cosurf(4, ips) / cosurf(4, jps)
        end if
        amat(((ips - 1)*ips)/2 + jps) = aa
        if (i == j) then
          rdat(i) = rdat(i) + 2 * aa * cosurf(4, ips) * cosurf(4, jps)
        end if
      end do
    end do
   !
   ! PERFORM CHOLESKY FACTORIZATION
   !
    call coscl1 (amat, nsetf, nps)
    deallocate (rdat, rsc, tm, nn, isort, ipsrs, nipsrs, nset, nipa, lipa, dirtm, finel)
    100 continue
end subroutine coscav
subroutine coscl1 (a, id, n)
   ! THIS ROUTINE PERFORMS A CHOLESKY FACTORIZATION
   ! INPUT:   A = PACKED LOWER TRIANGLE OF A
   !               SYMMETRIC POSITIVE DEFINITE N*N MATRIX
   ! OUTPUT:  A = LOWER TRIANGLE OF CHOLESKY MATRIX (INVERSE PIVOT ELEMEN
   !         ID = INTEGER VECTOR OF LENGTH N CONTAINING THE INDICES I(I - 1)
   !
    implicit none
    integer, intent (in) :: n
    double precision, dimension (*), intent (inout) :: a
    integer, dimension (n), intent (inout) :: id
    integer :: i, indi, info
    indi = 0
    do i = 1, n
      id(i) = indi
      indi = indi + i
    end do
    call dpptrf('U', n, a, info)
    do i = 1, n
      a(i+id(i)) = 1.0d0/a(i+id(i))
    end do
end subroutine coscl1
subroutine coscl2 (a, id, x, y, n)
   ! THIS ROUTINE SOLVES THE LINEAR SYSTEM CX = Y BASED ON
   ! CHOLESKY FACTORIZATION
   ! INPUT:   A = UPPER TRIANGLE OF CHOLESKY MATRIX
   ! INPUT:   Y = VECTOR OF LENGTH N
   !         ID = INTEGER VECTOR OF LENGTH N CONTAINING THE INDICES I(I - 1)
   ! OUTPUT:  X = VECTOR OF LENGTH N
    implicit none
    integer, intent (in) :: n
    double precision, dimension (*), intent (in) :: a
    integer, dimension (n), intent (in) :: id
    double precision, dimension (n), intent (out) :: x
    double precision, dimension (n), intent (in) :: y
    integer :: i, k, kk
    double precision :: summe
    x = 0.d0
    do k = 1, n
      summe = y(k)
      kk = id(k)
      do i = k - 1, 1, - 1
        summe = summe - a(i + kk) * x(i)
      end do
      x(k) = summe * a(k + kk)
    end do
    do k = n, 1, - 1
      summe = x(k)
      do i = k + 1, n
        summe = summe - a(k + id(i)) * x(i)
      end do
      x(k) = summe * a(k + id(k))
    end do
end subroutine coscl2
subroutine diegrd (dxyz)
    use molkst_C, only : numat
    use cosmo_C, only : nps, fepsi, nipc, &
     cosurf, iatsp, isude, sude, qscnet, qdenet, qscat, arat
    use parameters_C, only : dd, qq
    use common_arrays_c, only : coord, nfirst, nlast, nat
    use funcon_C, only : ev, a0, fpc_9
    implicit none
    integer :: i, ia, iak, ial, ib, idel, iden, ix, j, k, l, nati, iii
    double precision :: bsurf, ddi, deab, dist2, dx, fact, ff, ff0, &
   & qqi2, qsk, rm2, rm4, xxx
    double precision, dimension (3) :: xk, xl, xx
    double precision, dimension (0:3, 10) :: db
    double precision, dimension (3, numat), intent (inout) :: dxyz
    intrinsic Min, Sqrt
    do i = 1, 10
      do ix = 1, 3
        db(ix, i) = 0.d0
      end do
    end do
    db(0, 1) = 1.d0
    fact = -ev * a0 * fpc_9
    do k = 1, nps
      iak = iatsp(k)
      do ix = 1, 3
        xk(ix) = cosurf(ix, k)
      end do
      qsk = qscnet(k, 3)
      do l = 1, k - 1
        ial = iatsp(l)
        if (ial /= iak) then
          dist2 = 0.d0
          do ix = 1, 3
            xxx = cosurf(ix, l) - xk(ix)
            xl(ix) = xxx
            dist2 = dist2 + xxx * xxx
          end do
          ff = qsk * qscnet(l, 3) * fact * dist2 ** (- 1.5d0) / fepsi
          do ix = 1, 3
            dxyz(ix, iak) = dxyz(ix, iak) - xl(ix) * ff
            dxyz(ix, ial) = dxyz(ix, ial) + xl(ix) * ff
          end do
        end if
      end do
    end do
   ! NOW THE CONTRIBUTIONS FROM ANALYTIC SURFACE DERIVATIVES
   ! FOR THE SURFACE CLOSURE
    bsurf = 0.d0
    do i = 1, nipc
      ia = isude(1, i)
      ib = isude(2, i)
      deab = -0.25d0 * (qscat(ia)**2*sude(1, i)/arat(ia) + qscat(ib)**2*sude(2, &
     & i)/arat(ib) + bsurf*(sude(1, i) + sude(2, i)))
      xk(1) = coord(1, ib) - coord(1, ia)
      xk(2) = coord(2, ib) - coord(2, ia)
      xk(3) = coord(3, ib) - coord(3, ia)
      deab = deab / Sqrt (xk(1)**2 + xk(2)**2 + xk(3)**2)
      do ix = 1, 3
        dxyz(ix, ia) = dxyz(ix, ia) - xk(ix) * deab
        dxyz(ix, ib) = dxyz(ix, ib) + xk(ix) * deab
      end do
    end do
    do k = 1, nps
      iak = iatsp(k)
      do ix = 1, 3
        xk(ix) = cosurf(ix, k)
      end do
      qsk = qscnet(k, 3)
      iden = 0
      do i = 1, numat
        idel = nlast(i) + 1 - nfirst(i)
        if (i /= iak) then
          nati = nat(i)
          dist2 = 0.d0
          do ix = 1, 3
            xxx = xk(ix) - coord(ix, i)
            xx(ix) = xxx
            dist2 = dist2 + xxx * xxx
          end do
          ddi = dd(nati) * a0
          qqi2 = (a0*qq(nati)) ** 2
          ff0 = -qsk * fact * dist2 ** (- 1.5d0)
          if (idel  >  1) then
            rm2 = 1.d0 / dist2
            rm4 = rm2 ** 2
            db(0, 2) = ddi * 3 * xx(1) * rm2
            db(0, 4) = ddi * 3 * xx(2) * rm2
            db(0, 7) = ddi * 3 * xx(3) * rm2
            db(0, 3) = 1.d0 + qqi2 * (15*xx(1)**2*rm2 - 3.d0) * rm2
            db(0, 6) = 1.d0 + qqi2 * (15*xx(2)**2*rm2 - 3.d0) * rm2
            db(0, 10) = 1.d0 + qqi2 * (15*xx(3)**2*rm2 - 3.d0) * rm2
            db(0, 5) = qqi2 * 15 * xx(1) * xx(2) * rm4
            db(0, 8) = qqi2 * 15 * xx(1) * xx(3) * rm4
            db(0, 9) = qqi2 * 15 * xx(3) * xx(2) * rm4
            db(1, 2) = ddi
            db(2, 4) = db(1, 2)
            db(3, 7) = db(1, 2)
            db(1, 3) = 6 * qqi2 * xx(1) * rm2
            db(2, 6) = 6 * qqi2 * xx(2) * rm2
            db(3, 10) = 6 * qqi2 * xx(3) * rm2
            db(1, 5) = db(2, 6)
            db(2, 5) = db(1, 3)
            db(1, 8) = db(3, 10)
            db(3, 8) = db(1, 3)
            db(2, 9) = db(3, 10)
            db(3, 9) = db(2, 6)
          end if
          do j = 1, min(10, (idel*(idel + 1))/2)
            ff = -ff0 * qdenet(iden + j, 3)
            if(j == 1 .and. idel == 9) then
              do iii = 5, 9
                ff = ff - ff0 * qdenet(iden + (iii*(iii + 1))/2, 3)
              end do
            end if
            do ix = 1, 3
              dx = (xx(ix)*db(0, j) - db(ix, j)) * ff
              dxyz(ix, iak) = dxyz(ix, iak) + dx
              dxyz(ix, i) = dxyz(ix, i) - dx
            end do
          end do
        end if
        iden = iden + (idel*(idel + 1))/2
      end do
    end do
end subroutine diegrd
subroutine dmecip (coeffs, deltap, delta, eig, vectci, conf)
    use molkst_C, only : numat, norbs, lm61
    use meci_C, only : lab, nmos, nstate, nelec, microa, microb, occa, nalmat
    use cosmo_C, only : nps, fnsq, amat, bmat, nsetf, &
    phinet, qscnet, qdenet
    use common_arrays_C, only : nfirst, nlast
    use funcon_C, only : ev, a0
    implicit none
    double precision, dimension (lab), intent (inout) :: eig, vectci
    double precision, dimension (lab**2), intent (in) :: conf
    double precision, dimension (norbs, norbs), intent (in) :: coeffs
    double precision, dimension (norbs, nmos), intent (inout) :: delta
    double precision, dimension (nmos, nmos), intent (inout) :: deltap
    integer :: i, iat, id, iden, ii, ij, ist, ix, iy, j, jd, k
    double precision :: edie, fact, sum

    fact = -ev * a0 * fnsq / 2
    do ist = 1, lab
      ii = lab * (ist - 1)
      do j = 1, lab
        vectci(j) = conf(ii + j)
      end do
      !     INITIALIZE WITH THE OPPOSITE OF THE 'SCF' DENSITY.
      do i = 1, nmos
        deltap(i, i) = -occa(i) * 2.d0
        do j = 1, i - 1
          deltap(i, j) = 0.d0
        end do
      end do
      !
      !     ADD THE C.I. CORRECTION
      do id = 1, lab
        do jd = 1, id
            !     CHECK SPIN AGREEMENT
          if (nalmat(id) == nalmat(jd)) then
            ix = 0
            iy = 0
            do j = 1, nmos
              ix = ix + Abs (microa(j, id) - microa(j, jd))
              iy = iy + Abs (microb(j, id) - microb(j, jd))
            end do
               !     CHECK NUMBER OF DIFFERING M.O.
            if (ix + iy  <= 2) then
              if (ix == 2) then
                ! DETERMINANTS ID AND JD DIFFER BY M.O I IN ID AND M.O J IN JD:
                do i = 1, nmos
                  if (microa(i, id) /= microa(i, jd)) exit
                end do
                ij = microb(i, id)
                do j = i + 1, nmos
                  if (microa(j, id) /= microa(j, jd)) exit
                  ij = ij + microa(j, id) + microb(j, id)
                end do
                !        IJ GIVES THE SIGN OF THE PERMUTATION
                sum = 0.d0
                do k = 1, nstate
                  sum = sum + vectci(id + (k - 1)*lab) * vectci(jd + (k - 1)*lab)
                end do
                deltap(j, i) = deltap(j, i) + sum * Dble (1 - 2*Mod(ij, 2)) / &
               & nstate
              else if (iy == 2) then
                ! DETERMINANTS ID AND JD DIFFER BY M.O J IN ID AND M.O I IN JD:
                do i = 1, nmos
                  if (microb(i, id) /= microb(i, jd)) exit
                end do
                ij = 0
                do j = i + 1, nmos
                  if (microb(j, id) /= microb(j, jd)) exit
                  ij = ij + microa(j, id) + microb(j, id)
                end do
                ij = ij + microa(j, id)
                sum = 0.d0
                do k = 1, nstate
                  sum = sum + vectci(id + (k - 1)*lab) * vectci(jd + (k - 1)*lab)
                end do
                deltap(j, i) = deltap(j, i) + sum * Dble (1 - 2*Mod(ij, 2)) / &
                     & nstate
              else
                     !        DETERMINANTS ID AND JD ARE IDENTICAL:
                sum = 0.d0
                do k = 1, nstate
                  sum = sum + vectci(id + (k - 1)*lab) ** 2
                end do
                do i = 1, nmos
                  deltap(i, i) = deltap(i, i) &
                       & + (microa(i, id) + microb(i, id)) * sum / nstate
                end do
              end if
            end if
          end if
        end do
      end do
      !
      !     BACK TRANSFORM INTO A.O. BASIS.
      !    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !     P(C.I.) = P(SCF) + C * DELTAP * C'
      do i = 1, nmos
        do j = 1, i - 1
          deltap(j, i) = deltap(i, j)
        end do
      end do
      !     STEP 1: DELTAP = C * DELTAP
      call mxm (coeffs(1, nelec + 1), norbs, deltap, nmos, delta, nmos)
      !     STEP 2: P = P + DELTAP * C'
      iden = 0
      do iat = 1, numat
        do i = nfirst(iat), nlast(iat)
          do j = nfirst(iat), i
            sum = 0.d0
            do k = 1, nmos
              sum = sum + delta (i, k) * coeffs(j, nelec + k)
            end do
            iden = iden + 1
            qdenet(iden, 1) = sum
          end do
        end do
      end do
      ! FIRST CALCULATE THE CORRESPONDING POTENTIAL AS BMAT*QDENET
      do i = 1, nps
        sum = 0.d0
        do j = 1, lm61
          sum = sum + bmat(j, i) * qdenet(j, 1)
        end do
        phinet(i, 1) = sum
      end do
      ! NOW CALCULATE THE CORRESPONDING SCREENING CHARGES FROM A*QSCNET = -QDENE
      call coscl2 (amat, nsetf, qscnet, phinet, nps)
      edie = 0.d0
      ! NOW CALCULATE THE DIELECTRIC ENERGY AS QSCNET*PHINET
      do i = 1, nps
        edie = edie + qscnet(i, 1) * phinet(i, 1)
      end do
      edie = fact * edie
      eig(ist) = eig(ist) + edie
    end do
    do j = 1, lab
      vectci(j) = conf(j)
    end do
end subroutine dmecip
subroutine dvfill (nppa, dirvec)
   !***********************************************************************
   !
   !   DVFILL CALCULATES THE DIRECTION VECTORS.  THESE ARE PUT IN DIRVEC.
   !
   !  (THE VARIOUS SETS FORM ALMOST REGULAR POLYHEDRA WITH IH
   !   SYMMETRY.  A GOOD EXERCISE IS TO SKETCH THE POINTS USING,
   !   E.G. CHEM3D)
   !
   !***********************************************************************
    use chanel_C, only : iw
    use funcon_C, only : pi
    implicit none
    integer, intent (in) :: nppa
    double precision, dimension (4, nppa), intent (inout) :: dirvec
    integer :: i, ix, j, j1, j2, k, m, m2, na, nb, nc, nd
    double precision :: ar, beta, cphi, dist, dist2, h, r, sphi, sumar, t, &
   & xx, yy
    integer, dimension (2, 30) :: kset
    integer, dimension (3, 20) :: fset
    data kset / 1, 2, 1, 3, 1, 4, 1, 5, 1, 6, 12, 11, 12, 10, 12, 9, 12, 8, &
   & 12, 7, 2, 3, 3, 4, 4, 5, 5, 6, 6, 2, 7, 8, 8, 9, 9, 10, 10, 11, 11, 7, 2, &
   & 7, 7, 3, 3, 8, 8, 4, 4, 9, 9, 5, 5, 10, 10, 6, 6, 11, 11, 2 /
    data fset / 1, 2, 3, 1, 3, 4, 1, 4, 5, 1, 5, 6, 1, 6, 2, 12, 11, 10, 12, &
   & 10, 9, 12, 9, 8, 12, 8, 7, 12, 7, 11, 2, 3, 7, 3, 4, 8, 4, 5, 9, 5, 6, &
   & 10, 6, 2, 11, 7, 8, 3, 8, 9, 4, 9, 10, 5, 10, 11, 6, 11, 7, 2 /
    dirvec(1, 1) = -1.d0
    dirvec(2, 1) = 0.d0
    dirvec(3, 1) = 0.d0
    nd = 1
    r = Sqrt (0.8d0)
    h = Sqrt (0.2d0)
    do i = -1, 1, 2
      do j = 1, 5
        nd = nd + 1
        beta = 1.d0 + j * 0.4d0 * pi + (i + 1) * 0.1d0 * pi
        dirvec(2, nd) = r * Cos (beta)
        dirvec(3, nd) = r * Sin (beta)
        dirvec(1, nd) = i * h
      end do
    end do
    dirvec(1, 12) = 1.d0
    dirvec(2, 12) = 0.d0
    dirvec(3, 12) = 0.d0
    nd = 12
    cphi = Cos (1.d0)
    sphi = Sin (1.d0)
    do i = 1, 12
      xx = dirvec(1, i)
      yy = dirvec(2, i)
      dirvec(1, i) = cphi * xx + sphi * yy
      dirvec(2, i) = -sphi * xx + cphi * yy
    end do
   !   NPPA = 10*3**K*M**2 + 2
    m2 = (nppa - 2) / 10
    m = Nint (Sqrt(m2 + 0.d0))
    k = 0
    if (m2 /= m*m) then
      k = 1
      m2 = m2 / 3
      m = Nint (Sqrt(m2 + 0.d0))
    end if
    if (10*3**k*m**2 + 2 /= nppa) then
      write (iw, *) "VALUE OF NPPA NOT ALLOWED: IT MUST BE 10*3**K*M**2 + 2"
      return
    end if
   ! CREATE ON EACH EDGE M - 1 NEW POINTS
    do i = 1, 30
      na = kset(1, i)
      nb = kset(2, i)
      do j = 1, m - 1
        nd = nd + 1
        do ix = 1, 3
          dirvec(ix, nd) = dirvec(ix, na) * (m - j) + dirvec(ix, nb) * j
        end do
      end do
    end do
   ! CREATE POINTS WITHIN EACH TRIANGLE
    do i = 1, 20
      na = fset(1, i)
      nb = fset(2, i)
      nc = fset(3, i)
      do j1 = 1, m - 1
        do j2 = 1, m - j1 - 1
          nd = nd + 1
          do ix = 1, 3
            dirvec(ix, nd) = dirvec(ix, na) * (m - j1 - j2) + dirvec(ix, nb) * &
           & j1 + dirvec(ix, nc) * j2
          end do
        end do
      end do
      if (k /= 0) then
         ! CREATE TWO ADDITIONAL SUBGRIDS
        t = 1.d0 / 3
        do j1 = 0, m - 1
          do j2 = 0, m - j1 - 1
            nd = nd + 1
            do ix = 1, 3
              dirvec(ix, nd) = dirvec(ix, na) * (m - j1 - j2 - 2*t) + dirvec(ix, &
             & nb) * (j1 + t) + dirvec(ix, nc) * (j2 + t)
            end do
          end do
        end do
        t = 2.d0 / 3
        do j1 = 0, m - 2
          do j2 = 0, m - j1 - 2
            nd = nd + 1
            do ix = 1, 3
              dirvec(ix, nd) = dirvec(ix, na) * (m - j1 - j2 - 2*t) + dirvec(ix, &
             & nb) * (j1 + t) + dirvec(ix, nc) * (j2 + t)
            end do
          end do
        end do
      end if
    end do
   ! NORMALIZE ALL VECTORS
    sumar = 0.d0
    do i = 1, nppa
      dist = 0.d0
      do ix = 1, 3
        dist = dist + dirvec(ix, i) ** 2
      end do
      dist = 1.d0 / Sqrt (dist)
      dist2 = (m*dist) ** 2
      do ix = 1, 3
        dirvec(ix, i) = dirvec(ix, i) * dist
      end do
      if (i  <= 12) then
        ar = 5.d0
      else
        ar = 6.d0 * dist2
      end if
      dirvec(4, i) = ar
      sumar = sumar + ar
    end do
    sumar = 4 * pi / sumar
    do i = 1, nppa
      dirvec(4, i) = dirvec(4, i) * sumar
    end do
end subroutine dvfill
subroutine mfinel (ips, k, finel, nar_csm, nsetf, nset, rsc, nipsrs, dirvec, tm, &
& x, r, nfl, ioldcv, maxrs, lenabc, numat)
   !  THIS ROUTINE GENERATES THE LIST OF ALL BASIC GRID POINTS AND
   !  RING - SEGMENTS BELONGING TO SEGMENT IPS
    implicit none
    integer, intent (in) :: ioldcv, ips, k, lenabc, maxrs, numat
    integer, intent (out) :: nfl
    double precision, intent (in) :: r
    integer, dimension (1082*numat), intent (in) :: nset
    integer, dimension (lenabc), intent (in) :: nar_csm, nipsrs, nsetf
    double precision, dimension (3), intent (in) :: x
    double precision, dimension (3, 3), intent (in) :: tm
    double precision, dimension (4, 1082), intent (in) :: dirvec
    double precision, dimension (4, maxrs), intent (in) :: rsc
    double precision, dimension (4, 500, 2), intent (out) :: finel
    integer :: idir, irs, irs0, irs1, ix, l, nari
    double precision, dimension (3) :: y
   !
   ! FIRST THE TRANSFORMED BASIC GRID POINTS
    nfl = 0
    nari = nar_csm(ips)
    do l = nsetf(ips), nsetf(ips) + nari - 1
      idir = nset(l)
      nfl = nfl + 1
      do ix = 1, 3
        y(ix) = dirvec(ix, idir) * r
      end do
      finel(1, nfl, k) = y(1) * tm(1, 1) + y(2) * tm(2, 1) + y(3) * tm &
     & (3, 1) + x(1)
      finel(2, nfl, k) = y(1) * tm(1, 2) + y(2) * tm(2, 2) + y(3) * tm &
     & (3, 2) + x(2)
      finel(3, nfl, k) = y(1) * tm(1, 3) + y(2) * tm(2, 3) + y(3) * tm &
     & (3, 3) + x(3)
      finel(4, nfl, k) = dirvec(4, idir) * r * r
    end do
    if (ioldcv == 1) return
   ! NOW THE ASSOCIATED RINGSEGMENTS
    irs0 = 1
    if (ips > 1) then
      irs0 = nipsrs(ips - 1) + 1
    end if
    irs1 = nipsrs(ips)
    do irs = irs0, irs1
      nfl = nfl + 1
      do ix = 1, 4
        finel(ix, nfl, k) = rsc (ix, irs)
      end do
    end do
end subroutine mfinel
subroutine mkbmat
    use molkst_C, only : numat, method_indo
    use cosmo_C, only : nps, bmat, cosurf
    use parameters_C, only : dd, qq, dsd, dpd, ddd
    use funcon_C, only : a0
    use common_arrays_C, only : coord, nfirst, nlast, nat
    implicit none
    integer :: i, ia, idel, iden, ips, ix, nati, iii
    double precision :: ddi, dist, qqi2, rm1, rm3, rm5, dsdi, dpdi, dddi, atmp
    double precision, dimension (3) :: xa
   ! FILLING B - MATRIX
    iden = 0
    do i = 1, numat
      ia = nfirst(i)
      idel = nlast(i) - ia + 1
      nati = nat(i)
      ddi = dd(nati) * a0
      qqi2 = (a0*qq(nati)) ** 2
!
! replace dd and qq with terms that work for d orbitals
!
      dsdi = (dsd(nati)*a0)**2
      dpdi = dpd(nati)*a0
      dddi = (ddd(nati)*a0)**2
      do ips = 1, nps
        dist = 0.d0
        do ix = 1, 3
          xa(ix) = cosurf(ix, ips) - coord(ix, i)
          dist = dist + xa(ix) ** 2
        end do
        rm1 = 1.d0 / Sqrt (dist)
        bmat(iden + 1, ips) = rm1
        if (idel > 1) then
          rm3 = rm1 ** 3
          rm5 = rm1 ** 5
!
!   INDO off - diagonal terms are zero - preset here
!
          do iii = iden + 2, iden + 9
            bmat(iii, ips) = 0.d0
          end do
          bmat(iden + 3, ips) = rm1 + 3 * xa(1) ** 2 * qqi2 * rm5 - qqi2 * &
         & rm3
          bmat(iden + 6, ips) = rm1 + 3 * xa(2) ** 2 * qqi2 * rm5 - qqi2 * &
         & rm3
          bmat(iden + 10, ips) = rm1 + 3 * xa(3) ** 2 * qqi2 * rm5 - qqi2 &
         & * rm3
!
!  s - p terms (2, 4, 7) and p - p terms (5, 8, 9) for non - INDO methods
!
          if (.not. method_indo) then
            bmat(iden + 2, ips) = xa(1) * ddi * rm3
            bmat(iden + 4, ips) = xa(2) * ddi * rm3
            bmat(iden + 7, ips) = xa(3) * ddi * rm3
            bmat(iden + 5, ips) = 3 * xa(1) * xa(2) * qqi2 * rm5
            bmat(iden + 8, ips) = 3 * xa(1) * xa(3) * qqi2 * rm5
            bmat(iden + 9, ips) = 3 * xa(3) * xa(2) * qqi2 * rm5
          end if
! c now the d-orbitals in spherical symmetry
          if(idel > 4) then
            do iii = iden + 11, iden + 44
              bmat(iii, ips) = 0.d0
            end do
            do iii = 5, 9
              bmat(iden + (iii*(iii + 1))/2, ips) = rm1
            end do
! Replace spherically symmetric d orbitals with terms from
! J. Phys. Chem. B 2005, 109, 9799
            if (method_indo) then
              do iii = iden + 11, iden + 44
                bmat(iii, ips) = 0.d0
              end do
! Order of INDO d orbitals:  z^2, x^2 - y^2, xy, xz, yz
! d-d self interactions
              bmat(iden + 15, ips) = rm1 - dddi*rm3 + 3*dddi*xa(3)**2*rm5   ! z^2
              bmat(iden + 21, ips) = rm1 + dddi*rm3 - 3*dddi*xa(3)**2*rm5   ! x^2 - y^2
              bmat(iden + 28, ips) = rm1 + dddi*rm3 - 3*dddi*xa(3)**2*rm5   ! xy
              bmat(iden + 36, ips) = rm1 + dddi*rm3 - 3*dddi*xa(2)**2*rm5   ! xz
              bmat(iden + 45, ips) = rm1 + dddi*rm3 - 3*dddi*xa(1)**2*rm5   ! yz
            else
! Order of MOPAC d orbitals: x^2 - y^2, xz, z^2, yz, xy
              bmat(iden + 15, ips) = rm1 + dddi*rm3 - 3*dddi*xa(3)**2*rm5   ! x^2 - y^2
              bmat(iden + 21, ips) = rm1 + dddi*rm3 - 3*dddi*xa(2)**2*rm5   ! xz
              bmat(iden + 28, ips) = rm1 - dddi*rm3 + 3*dddi*xa(3)**2*rm5   ! z^2
              bmat(iden + 36, ips) = rm1 + dddi*rm3 - 3*dddi*xa(1)**2*rm5   ! yz
              bmat(iden + 45, ips) = rm1 + dddi*rm3 - 3*dddi*xa(3)**2*rm5   ! xy

!     s - d terms
              atmp = 2*xa(3)**2 - xa(1)**2 - xa(2)**2
              bmat(iden + 11, ips) = 3*dsdi*(xa(1)**2 - xa(2)**2)*rm5/2.d0    ! x^2 - y^2
              bmat(iden + 16, ips) = 3*dsdi*xa(1)*xa(3)*rm5                 ! xz
              bmat(iden + 22, ips) = sqrt(3.d0)*atmp*rm5/2.d0               ! z^2
              bmat(iden + 29, ips) = 3*dsdi*xa(2)*xa(3)*rm5                 ! yz
              bmat(iden + 37, ips) = 3*dsdi*xa(1)*xa(2)*rm5                 ! xy

!     p - d terms
!             dx^2 - y^2
              bmat(iden + 12, ips) = dpdi*xa(1)*rm3
              bmat(iden + 13, ips) = dpdi*xa(2)*rm3
              bmat(iden + 14, ips) = 0.d0
!             dxz
              bmat(iden + 17, ips) = dpdi*xa(3)*rm3
              bmat(iden + 18, ips) = 0.d0
              bmat(iden + 19, ips) = dpdi*xa(1)*rm3
!             dz^2
              bmat(iden + 23, ips) = -dpdi*xa(1)*rm3/sqrt(3.d0)
              bmat(iden + 24, ips) = -dpdi*xa(2)*rm3/sqrt(3.d0)
              bmat(iden + 25, ips) = 2*dpdi*xa(3)*rm3/sqrt(3.d0)
!             dyz
              bmat(iden + 30, ips) = 0.d0
              bmat(iden + 31, ips) = dpdi*xa(3)*rm3
              bmat(iden + 32, ips) = dpdi*xa(2)*rm3
!             dxy
              bmat(iden + 38, ips) = dpdi*xa(2)*rm3
              bmat(iden + 39, ips) = dpdi*xa(1)*rm3
              bmat(iden + 40, ips) = 0.d0

!     d-d non - self terms
              bmat(iden + 20, ips) = 3*dddi*xa(1)*xa(3)*rm5             ! dxz - dx^2 - y^2
              bmat(iden + 26, ips) = 3*dddi*rm5*(xa(1)**2 - xa(2)**2)   ! dz^2 - dx^2 - y^2
              bmat(iden + 27, ips) = sqrt(3.d0)*dddi*rm5*xa(1)*xa(3)    ! dz^2 - dxz
              bmat(iden + 33, ips) = 3*dddi*xa(2)*xa(3)*rm5             ! dyz - dx^2 - y^2
              bmat(iden + 34, ips) = 3*dddi*xa(1)*xa(2)*rm5             ! dyz - dxz
              bmat(iden + 35, ips) = sqrt(3.d0)*dddi*rm5*xa(2)*xa(3)    ! dyz - dz^2
              bmat(iden + 41, ips) = 0                                  ! dxy - dx^2 - y^2
              bmat(iden + 42, ips) = 3*dddi*xa(2)*xa(3)*rm5             ! dxy - dxz
              bmat(iden + 43, ips) = -2*sqrt(3.d0)*dddi*rm5*xa(1)*xa(2) ! dxy - dz^2
              bmat(iden + 44, ips) = 3*dddi*xa(1)*xa(3)*rm5             ! dxy - dyz
            end if
          end if
! c end d-orbitals
        end if
      end do
      iden = iden + (idel*(idel + 1))/2
    end do
end subroutine mkbmat
subroutine surclo (coord, nipa, lipa, din, dim_din, rsc, isort, ipsrs, nipsrs, nat, &
& srad, cosurf, iatsp, nar_csm, nsetf, isude, sude, maxrs)
   ! CREATES SEGMENTS WHICH CLOSE THE CONCAVE REGIONS OF THE CAVITY
    use cosmo_C, only : n0, lenabc, nipc, rsolv, area, &
       & cosvol, nps
    use molkst_C, only : numat
    use funcon_C, only : pi, twopi
    implicit none
    integer, intent (in) :: maxrs, dim_din
    logical, dimension (dim_din), intent (in) :: din
    integer, dimension (lenabc + 1), intent (inout) :: iatsp, nar_csm, nsetf
    integer, dimension (maxrs), intent (in) :: lipa
    integer, dimension (maxrs), intent (out) :: ipsrs, isort
    integer, dimension (numat), intent (in) :: nat, nipa
    integer, dimension (lenabc), intent (out) :: nipsrs
    integer, dimension (2, 30*numat), intent (out) :: isude
    double precision, dimension (numat), intent (in) :: srad
    double precision, dimension (2, 30*numat), intent (out) :: sude
    double precision, dimension (3, numat), intent (in) :: coord
    double precision, dimension (4, lenabc + 1), intent (inout) :: cosurf
    double precision, dimension (4, maxrs), intent (out) :: rsc
    integer :: i, i2, i3, ia, iat, iat0, ib, ic, ich, iib, iic, iii, &
         & ik, il, ilipa, ip1, ips, ips0, ips1, ipsmin, is, isum, &
         & isum2, ix, ja, jb, jc, k, l, nrs, nsa, nsab, ntrp, ntrp2
    double precision :: aa, aad, aar, ab, abd, abr, arseg, arsegn, ca, cj, &
         & cnrs, cosa, cosb, cphi, d2, d2max, d2min, da, dab, dabc, &
         & dabck, ddd, dist, dp, hh, htr, phi, phio, phiu, ra, rb, rc, &
         & sina, sj, sp
    double precision :: sphi, spn, spn2, spx, spy, sss, sumphi, &
         & xxx, dac, dbc, rinc, tar
    integer, dimension (50) :: iset
    double precision, dimension (2) :: fz, yx
    double precision, dimension (3) :: rvx, rvy, trp, tvx, tvy, xd, xja, &
         & xjb, xjc, xx
    double precision, dimension (50) :: phiset
    double precision, dimension (3, 3) :: ee, xta
    double precision, dimension (50) :: tarset
    htr = 0.d0
    tarset = 0.d0
    nipc = 0
    nrs = 0
    iset = 0
   ! GENERATION OF SEGMENTS ALONG THE INTERSECTION RINGS
    ilipa = 0
    do i = 1, numat
      nipsrs(i) = 0
    end do
    do ia = 1, numat - 1
      if (.not. din(ia)) then
        ra = srad(ia) + rsolv
        do ix = 1, 3
          xta(ix, 1) = coord(ix, ia)
        end do
        do iib = ilipa + 1, ilipa + nipa(ia)
          ib = lipa(iib)
          if (ib > ia) then
            if (.not. din(ib)) then
              rb = srad(ib) + rsolv
              dab = 0.d0
              nsab = 0
              do ix = 1, 3
                xta(ix, 2) = coord(ix, ib)
                xx(ix) = xta(ix, 2) - xta(ix, 1)
                dab = dab + xx(ix) ** 2
              end do
              dab = Sqrt (dab)
              call ansude (ra - rsolv, rb - rsolv, dab, rsolv, aa, ab, aar, abr, &
             & aad, abd, rinc)
              cosa = (ra**2 + dab**2 - rb**2) / (2*dab*ra)
              cosb = (rb**2 + dab**2 - ra**2) / (2*dab*rb)
              sina = Sqrt (1.d0 - cosa**2)
              da = ra * cosa
              hh = ra * sina
              ddd = rsolv * (cosa + cosb) / dab
              fz(1) = (1.d0 - Cos (hh*pi/ra)) / 2
              fz(2) = (1.d0 - Cos (hh*pi/rb)) / 2
              if (cosa*cosb < 0d0) then
                fz(1) = 1.d0
              end if
              if (cosa*cosb < 0d0) then
                fz(2) = 1.d0
              end if
              yx(1) = rsolv / ra
              yx(2) = rsolv / rb
              do ix = 1, 3
                xd(ix) = xta(ix, 1) + da * xx(ix) / dab
              end do
                  !  CREATE RING VECTORS
              rvx(1) = xx(2) * 3.d0 - xx(3) * 2.d0
              rvx(2) = xx(3) * 1.d0 - xx(1) * 3.d0
              rvx(3) = xx(1) * 2.d0 - xx(2) * 1.d0
              dist = Sqrt (rvx(1)**2 + rvx(2)**2 + rvx(3)**2)
              if(dist < 1d-9) then
                rvx(1) = xx(2)*4.d0 - xx(3)*2.d0
                rvx(2) = xx(3)*1.d0 - xx(1)*4.d0
                rvx(3) = xx(1)*2.d0 - xx(2)*1.d0
                dist = sqrt(rvx(1)**2 + rvx(2)**2 + rvx(3)**2)
              end if
              do ix = 1, 3
                rvx(ix) = hh * rvx(ix) / dist
              end do
              rvy(1) = (xx(2)*rvx(3) - xx(3)*rvx(2)) / dab
              rvy(2) = (xx(3)*rvx(1) - xx(1)*rvx(3)) / dab
              rvy(3) = (xx(1)*rvx(2) - xx(2)*rvx(1)) / dab
                  ! NOW ALL TRIPLE POINTS ON THE RING ARE SEARCHED
              ntrp = 0
              ntrp2 = 0
              do iic = ilipa + 1, ilipa + nipa(ia)
                ic = lipa(iic)
                if (ic /= ib) then
                  rc = srad(ic) + rsolv
                  dabc = 0.d0
                  sp = 0.d0
                  do ix = 1, 3
                    xxx = coord(ix, ic) - xd(ix)
                    xta(ix, 3) = coord(ix, ic)
                    sp = sp + xxx * xx(ix)
                    dabc = dabc + xxx ** 2
                  end do
                  dabc = Sqrt (dabc)
                  cosa = sp / dab / dabc
                           ! AVOID PROBLEMS WITH NUMERICAL ACCURACY
                  sina = Sqrt (Max(1.d-28, 1.d0 - cosa*cosa))
                  cj = (dabc*dabc + hh*hh - rc*rc) / (2*dabc*hh*sina)
                  if (cj < 1.d0) then
                    ntrp2 = ntrp2 + 2
                  end if
                  if (cj  <= 1.d0 .and. cj >= -1.d0) then
                    sj = Sqrt (1.d0 - cj*cj)
                              !
                    do ix = 1, 3
                      tvx(ix) = (xta(ix, 3) - xd(ix)) - cosa * dabc * xx(ix) / &
                     & dab
                    end do
                    dist = Sqrt (tvx(1)**2 + tvx(2)**2 + tvx(3)**2)
                    do ix = 1, 3
                      tvx(ix) = hh * tvx(ix) / dist
                    end do
                    tvy(1) = (xx(2)*tvx(3) - xx(3)*tvx(2)) / dab
                    tvy(2) = (xx(3)*tvx(1) - xx(1)*tvx(3)) / dab
                    tvy(3) = (xx(1)*tvx(2) - xx(2)*tvx(1)) / dab
                    loop: do l = -1, 1, 2
                      il = ntrp + 1
                      do ix = 1, 3
                        trp(ix) = xd(ix) + cj * tvx(ix) + sj * tvy(ix) * l
                      end do
                      do ik = ilipa + 1, ilipa + nipa(ia)
                        k = lipa(ik)
                        if (k /= ib .and. k /= ic) then
                          dabck = 0.d0
                          do ix = 1, 3
                            dabck = dabck + (trp(ix) - coord(ix, k)) ** 2
                          end do
                          dabck = Sqrt (dabck)
                          if (dabck < srad(k) + rsolv) cycle loop
                        end if
                      end do
                      ntrp = ntrp + 1
                      spx = 0.d0
                      spy = 0.d0
                      do ix = 1, 3
                        spx = spx + rvx(ix) * (trp(ix) - xd(ix))
                        spy = spy + rvy(ix) * (trp(ix) - xd(ix))
                      end do
                      phi = Acos (spx/(hh**2 + 1.d-10))
                      if (spy < 0.d0) then
                        phi = -phi
                      end if
                      phiset(il) = phi
                      sp = 0.d0
                      do ix = 1, 3
                        sp = sp + (- spy*rvx(ix) + spx*rvy(ix)) * &
                       & (trp(ix) - xta(ix, 3))
                      end do
                      iset(ntrp) = 1
                      if (sp < 0.d0) then
                        iset(ntrp) = -1
                      end if
! first calculate the edges of the triangle
                      sp = 0.d0
                      dac = 0.d0
                      dbc = 0.d0
                      do ix = 1, 3
                        ee(ix, 1) = trp(ix) + rsolv / srad(ia) &
                             & * (xta(ix, 1) - trp(ix))
                        ee(ix, 2) = trp(ix) + rsolv / srad(ib) &
                             & * (xta(ix, 2) - trp(ix))
                        ee(ix, 3) = trp(ix) + rsolv / srad(ic) &
                             & * (xta(ix, 3) - trp(ix))
                        sp = sp + (ee(ix, 1) - ee(ix, 3)) * (ee(ix, 2) - ee(ix, 3))
                        dac = dac + (ee(ix, 1) - ee(ix, 3)) &
                             & * (ee(ix, 1) - ee(ix, 3))
                        dbc = dbc + (ee(ix, 2) - ee(ix, 3)) &
                             & * (ee(ix, 2) - ee(ix, 3))
                      end do
                      tarset(il) = 0.8d0 * Sqrt (dac*dbc - sp*sp) / 12
                    end do loop
                  end if
                end if
              end do
                  ! SORT THE SET OF TRIPLE POINTS ON THE RING
              if (Mod(ntrp, 2) /= 0) then
                call mopend ("ODD NTRP")
              end if
              if (ntrp > 18) then
                call mopend ("NTRP TOO LARGE")
              end if
              if (ntrp + ntrp2 == 0) then
                phiset(1) = 0
                phiset(2) = twopi
                tarset(2) = 0.d0
                iset(1) = 1
                iset(2) = -1
                ntrp = 2
              end if
              do
                ic = 0
                do l = 2, ntrp
                  if (phiset(l) < phiset(l - 1)) then
                    phi = phiset(l)
                    iii = iset(l)
                    tar = tarset(l)
                    phiset(l) = phiset(l - 1)
                    iset(l) = iset(l - 1)
                    tarset(l) = tarset(l - 1)
                    phiset(l - 1) = phi
                    iset(l - 1) = iii
                    tarset(l - 1) = tar
                    ic = ic + 1
                  end if
                end do
                if (ic  <= 0) then
                  if (iset(1) == -1) then
                    phiset(1) = phiset(1) + twopi
                  else
                    exit
                  end if
                end if
              end do
                  ! NOW FOR EACH CONTINUOUS SECTION OF THE RING TRIANGLES ARE CREATED
              sumphi = 0.d0
              ips0 = nrs
              do l = 2, ntrp, 2
                k = l - 1
                phiu = phiset(k)
                phio = phiset(l)
                nsa = Int ((phio - phiu)/2/pi*20)
                nsa = Max (nsa + 1, 2)
                sumphi = sumphi + phio - phiu
                dp = (phio - phiu) / (nsa - 1)
                do ich = 1, 2
                  iat = ib
                  if (ich == 2) then
                    iat = ia
                  end if
                  if (iat == ia) then
                    htr = aar / twopi
                  end if
                  if (iat == ib) then
                    htr = abr / twopi
                  end if
                  do ja = ich, nsa, 2
                    jb = Max (ja - 1, 1)
                    jc = Min (ja + 1, nsa)
                    phi = phiu + (ja - 1) * dp
                    cphi = Cos (phi)
                    sphi = Sin (phi)
                    do ix = 1, 3
                      ca = xd(ix) + (cphi*rvx(ix) + sphi*rvy(ix)) * fz(ich)
                      ca = ca + (xta(ix, ich) - ca) * yx(ich)
                      xja(ix) = ca + (xta(ix, 3 - ich) - xta(ix, ich)) * ddd * &
                     & (1.d0 - fz(ich))
                    end do
                    phi = phiu + (jb - 1) * dp
                    cphi = Cos (phi)
                    sphi = Sin (phi)
                    do ix = 1, 3
                      ca = xd(ix) + cphi * rvx(ix) + sphi * rvy(ix)
                      xjb(ix) = ca + (xta(ix, 3 - ich) - ca) * yx(3 - ich)
                    end do
                    phi = phiu + (jc - 1) * dp
                    cphi = Cos (phi)
                    sphi = Sin (phi)
                    do ix = 1, 3
                      ca = xd(ix) + cphi * rvx(ix) + sphi * rvy(ix)
                      xjc(ix) = ca + (xta(ix, 3 - ich) - ca) * yx(3 - ich)
                    end do
                    nrs = nrs + 1
                    nsab = nsab + 1
                    sp = 0.d0
                    d2 = 0.d0
                    spn = 0.d0
                    spn2 = 0.d0
                    dist = 0.d0
                    ipsrs(nrs) = ib
                    if (ich == 2) then
                      ipsrs(nrs) = ia
                    end if
                    iat = ipsrs(nrs)
                    nipsrs(iat) = nipsrs(iat) + 1
                    do ix = 1, 3
                      rsc (ix, nrs) = (xja(ix)*.5d0 + xjb(ix) + xjc(ix)) / 2.5d0
                      i2 = Mod (ix, 3) + 1
                      i3 = Mod (i2, 3) + 1
                      cnrs = (xjc(i2) - xjb(i2)) * (xja(i3) - xjb(i3)) - &
                     & (xja(i2) - xjb(i2)) * (xjc(i3) - xjb(i3))
                      dist = dist + cnrs ** 2
                      spn = spn + cnrs * (rsc (ix, nrs) - coord(ix, iat))
                      spn2 = spn2 + cnrs * rsc (ix, nrs)
                    end do
                    dist = 1.d0 / Sqrt (dist)
                    if (spn < 0.d0) then
                      dist = -dist
                    end if
                    rsc(4, nrs) = (jc - jb)*dp*htr
                    if (ja  ==  1)   rsc (4, nrs) = rsc (4, nrs) + tarset(k) * rinc
                    if (ja  ==  nsa) rsc (4, nrs) = rsc (4, nrs) + tarset(l) * rinc
                    cosvol = cosvol + rsc (4, nrs) * spn2 * dist
                    do ik = ilipa + 1, ilipa + nipa(ia)
                      k = lipa(ik)
                      if (k  ==  ib) cycle
                      dabck = 0.d0
                      do ix = 1, 3
                        dabck = dabck + (rsc(ix, nrs) - coord(ix, k))**2
                      end do
                      dabck = sqrt(dabck)
                    end do
                    area = area + rsc (4, nrs)
                  end do
                end do
              end do
              if (sumphi > 1.d-10) then
                nipc = nipc + 1
                isude(1, nipc) = ia
                isude(2, nipc) = ib
                sumphi = sumphi / twopi
                sude(1, nipc) = aad * sumphi
                sude(2, nipc) = abd * sumphi
              end if
            end if
          end if
        end do
      end if
      ilipa = ilipa + nipa(ia)
    end do
    if (nrs > maxrs) then
      call mopend ("NRS .GT. MAXRS IN SURCLO")
    end if
   ! NOW SORT THE SEGMENTS WITH RESPECT TO ATOMS
    isum = 0
    do iat = 1, numat
      isum2 = isum + nipsrs(iat)
      nipsrs(iat) = isum
      isum = isum2
    end do
    do i = 1, nrs
      iat = ipsrs(i)
      isum2 = nipsrs(iat) + 1
      nipsrs(iat) = isum2
      isort(i) = isum2
    end do
    do i = 1, nrs
      do while (isort(i) /= i)
        is = isort(i)
        do ix = 1, 4
          sss = rsc (ix, i)
          rsc (ix, i) = rsc (ix, is)
          rsc (ix, is) = sss
        end do
        iii = ipsrs(i)
        ipsrs(i) = ipsrs(is)
        ipsrs(is) = iii
        isort(i) = isort(is)
        isort(is) = is
      end do
    end do
   ! NOW FIND FOR EACH RINGSEGMENT THE NEAREST SEGMENT OUT OF THE SAME ATOM
    do i = 1, nps
      nipsrs(i) = 0
    end do
    iat0 = 0
    ips0 = 0
    ips1 = 0
    d2max = 0.d0
    do i = 1, nrs
      iat = ipsrs(i)
      if (iat > iat0) then
        iat0 = iat
        d2max = 16 * srad(iat) ** 2 / n0(1)
        if (nat(iat) == 1) then
          d2max = d2max * n0(1) / n0(2)
        end if
        do ips = ips1 + 1, nps
          if (iatsp(ips) == iat) exit
        end do
        ips0 = ips
        do ips = ips0, nps - 1
          if (iatsp(ips + 1) > iat) exit
        end do
        ips1 = ips
      end if
      ipsmin = 0
      d2min = 1d6
      do ips = ips0, ips1
        d2 = 0.d0
        do ix = 1, 3
          d2 = d2 + (cosurf(ix, ips) - rsc (ix, i)) ** 2
        end do
        if (d2 < d2min) then
          ipsmin = ips
          d2min = d2
        end if
      end do
      if (d2min > d2max) then
         ! IF NO CLOSE ENOUGH SEGMENT IS PRESENT ADD A NEW ONE
        ips1 = ips1 + 1
        do ips = nps, ips1, - 1
          ip1 = ips + 1
          do ix = 1, 4
            cosurf(ix, ip1) = cosurf(ix, ips)
          end do
          iatsp(ip1) = iatsp(ips)
          nsetf(ip1) = nsetf(ips)
          nar_csm(ip1) = nar_csm(ips)
          nipsrs(ip1) = nipsrs(ips)
        end do
        ipsrs(i) = ips1
        nps = nps + 1
        do ix = 1, 4
          cosurf(ix, ips1) = rsc (ix, i)
        end do
        iatsp(ips1) = iatsp(ips1 - 1)
        nar_csm(ips1) = 0
        nipsrs(ips1) = 1
        nsetf(ips1) = nsetf(ips1 + 1)
      else
        ipsrs(i) = ipsmin
        nipsrs(ipsmin) = nipsrs(ipsmin) + 1
         ! UPDATE THE SEGMENT WITH RESPECT TO ADDED RING - SEGMENT
        arseg = cosurf(4, ipsmin)
        arsegn = arseg + rsc (4, i)
        do ix = 1, 3
          cosurf(ix, ipsmin) = (arseg*cosurf(ix, ipsmin) + rsc (4, i)*rsc (ix, &
         & i)) / arsegn
        end do
        cosurf(4, ipsmin) = arsegn
      end if
    end do
   ! NOW SORT THE RING - SEGMENTS WITH RESPECT TO PRIMARY SEGMENTS
    isum = 0
    do ips = 1, nps
      isum2 = isum + nipsrs(ips)
      nipsrs(ips) = isum
      isum = isum2
    end do
    do i = 1, nrs
      ips = ipsrs(i)
      isum2 = nipsrs(ips) + 1
      nipsrs(ips) = isum2
      isort(i) = isum2
    end do
    do i = 1, nrs
      do while (isort(i) /= i)
        is = isort(i)
        do ix = 1, 4
          sss = rsc (ix, i)
          rsc (ix, i) = rsc (ix, is)
          rsc (ix, is) = sss
        end do
        isort(i) = isort(is)
        isort(is) = is
      end do
    end do
end subroutine surclo
subroutine cosini(l_print)
    use cosmo_C, only : n0, ioldcv, fnsq, nps, rsolv, nspa, disex2, &
    dirsm, dirvec, srad, ipiden, gden, idenat, qdenet, amat, &
    cmat, lenabc, arat, sude, isude, bh, qden, nar_csm, nsetf, phinet, &
    qscnet, bmat, nset, xsp, abcmat, iatsp, nn, qscat, cosurf, nppa, &
    ffact, fepsi, nden
    use common_arrays_C, only : nat, nfirst, nlast
    use molkst_C, only : numat, keywrd, moperr, lm61, mozyme, line, jobnam
    use chanel_C, only : iw
    implicit none
    logical, intent (in) :: l_print
    integer :: i, i0, iat, idel, iden, incif, indise, inrsol, j, k, n1, &
   & n2, nfi, nfj, ns, nsp, nspd
    double precision :: disex, ri1, x
    double precision, dimension (107) :: rvdw, usevdw
    integer, external :: ijbo
    double precision, external :: reada
    data rvdw / &
   1.30d0,  1.64d0,  2.13d0,  2.19d0,  2.05d0,  2.00d0,  1.83d0,  1.72d0, &
   1.72d0,  1.80d0,  2.66d0,  2.02d0,  2.41d0,  2.46d0,  2.11d0,  2.16d0, &
   2.05d0,  2.20d0,  3.22d0,  2.54d0,  2.64d0,  2.64d0,  2.52d0,  2.40d0, &
   2.46d0,  2.41d0,  2.40d0,  1.91d0,  1.64d0,  1.63d0,  2.19d0,  2.46d0, &
   2.22d0,  2.22d0,  2.16d0,  2.36d0,  3.78d0,  3.44d0,  3.39d0,  3.33d0, &
   3.28d0,  2.57d0,  2.57d0,  2.57d0,  2.57d0,  1.91d0,  2.01d0,  1.85d0, &
   2.26d0,  2.54d0,  2.53d0,  2.41d0,  2.32d0,  2.53d0,  4.00d0,  3.47d0, &
   2.81d0,  2.81d0,  2.81d0,  2.81d0,  2.81d0,  2.81d0,  2.81d0,  2.81d0, &
   2.81d0,  2.81d0,  2.81d0,  2.81d0,  2.81d0,  2.81d0,  2.81d0,  2.57d0, &
   2.57d0,  2.57d0,  2.57d0,  2.57d0,  2.57d0,  2.05d0,  1.94d0,  1.81d0, &
   2.29d0,  2.36d0,  2.64d0,  2.64d0,  2.63d0,  2.69d0,  2.57d0,  2.57d0, &
   2.57d0,  2.57d0,  2.57d0,  2.18d0,  2.57d0,  2.57d0,  2.57d0,  2.57d0, &
   2.57d0,  2.57d0,  2.57d0,  2.57d0,  2.57d0,  2.57d0,  5*2.d0/
    lenabc = max(100, nspa*numat)
    ns = 0
    nsp = 0
    nspd = 0
    do i = 1, numat
      if (nlast(i) - nfirst(i) == 0) then
        ns = ns + 1
      else if (nlast(i) - nfirst(i) == 3) then
        nsp = nsp + 1
      else
        nspd = nspd + 1
      end if
    end do
    nden = 1*ns + 10*nsp + 45*nspd
    if (allocated(ipiden)) deallocate (ipiden)
    if (allocated(idenat)) deallocate (idenat)
    if (allocated(gden))   deallocate (gden)
    if (allocated(qdenet)) deallocate (qdenet)
    if (allocated(phinet)) deallocate (phinet)
    if (allocated(qscnet)) deallocate (qscnet)
    if (allocated(abcmat)) deallocate (abcmat)
    if (allocated(bmat))   deallocate (bmat)
    if (allocated(cmat))   deallocate (cmat)
    if (allocated(qscat))  deallocate (qscat)
    if (allocated(srad))   deallocate (srad)
    if (allocated(iatsp))  deallocate (iatsp)
    if (allocated(nar_csm))deallocate (nar_csm)
    if (allocated(nn))     deallocate (nn)
    if (allocated(cosurf)) deallocate (cosurf)
    if (allocated(xsp))    deallocate (xsp)
    if (allocated(nset))   deallocate (nset)
    if (allocated(bh))     deallocate (bh)
    if (allocated(qden))   deallocate (qden)
    if (allocated(nar_csm))    deallocate (nar_csm)
    if (allocated(nsetf))  deallocate (nsetf)
    if (allocated(isude))  deallocate (isude)
    if (allocated(sude))   deallocate (sude)
    if (allocated(arat))   deallocate (arat)
    if (allocated(amat))   deallocate (amat)
    allocate(ipiden(lm61), idenat(numat), gden(lm61), &
          qdenet(lm61, 3), phinet(lenabc + 1, 3), qscnet(lenabc + 1, 3), &
          qscat(numat), stat = i)
          qscat = 0.d0
    if (i /= 0) then
      call memory_error("COSINI (0) in Cosmo")
      return
    end if
    allocate(srad(numat), nn(3, numat), qden(lm61), &
          iatsp(lenabc + 1), isude(2, 30*numat), nar_csm(lenabc + 1), &
          arat(numat), sude(2, 30*numat), nsetf(lenabc + 1), cosurf(4, lenabc), stat = i)
    nsetf = 0 ! Added to fix gaps in nsetf's initialization, which might indicate a deeper flaw
    if (i /= 0) then
      call memory_error("COSINI (1) in Cosmo")
      return
    end if
    if (.not. mozyme) then
      allocate(abcmat(lenabc), xsp(3, lenabc), nset(nppa*numat), bh(lenabc), stat = j)
      allocate(bmat(lm61, lenabc), stat = i)
      j = j + i
      allocate(amat((lenabc*(lenabc + 1))/2), stat = i)
      j = j + i
      allocate(cmat((lm61*(lm61 + 1))/2), stat = i)
      j = j + i
      if (j /= 0) then
        if (l_print) then
          write(line, '(a, i5, a)')"Data set '"//trim(jobnam)//"' exists, "
          write(0, '(//10x, a)')trim(line)
          call mopend(trim(line))
          write(line, '(a)')"but is too large to run using COSMO."
          write(0, '(//10x, a)')trim(line)
          call mopend(trim(line))
          write(line, '(a, i5, a)')"(This job might run if MOZYME is used.)"
          write(0, '(/10x, a, //)')trim(line)
          write(iw, '(/10x, a, //)')trim(line)
        end if
        call memory_error("COSINI (2) in Cosmo")
        return
      end if
      cmat = 0.d0
    end if
    call extvdw (usevdw, rvdw)
    if (moperr) return
    rsolv = 1.30005d0
    ioldcv = 0
    inrsol = Index (keywrd, " RSOLV=")
    if (inrsol /= 0) then
      rsolv = reada (keywrd, inrsol)
    end if
    if (rsolv < 0.5d0) then
      write (iw, *) "RSOLV IS SET TO 0.5"
      return
    end if
    if (moperr) return
    ri1 = 2.d0
    incif = Index (keywrd, "N**2")
    if (incif /= 0) then
      ri1 = reada (keywrd, incif + 4)
      if (ri1 < 0.d0) then
        write (iw, *) "     N**2 CANNOT BE NEGATIVE"
        call mopend ("N**2 CANNOT BE NEGATIVE")
        return
      else if (0.d0 < ri1 .and. ri1 < 1.d0) then
        write (iw, *) "     N**2 IS SMALLER THAN 1.d0, OK?"
      end if
    end if
    fnsq = (ri1 - 1.d0) / (ri1 + .5d0)
!
!  Set ffact for later use
!
    if (fepsi > 0.d0) then
      ffact = fnsq/fepsi
    else
      ffact = 0.d0
    end if
    nps = 0
   ! NO KEYWORD DELSC ANYMORE. PEOPLE MAY USE EXPLICIT
   ! DEFINITION OF RADII TO CHANGE FROM DEFAULTS! RDS IS REPLACED BY RSOLV
    rsolv = Max (rsolv, 0.5d0)
    disex = 4.d0
    indise = Index (keywrd, " DISEX=")
    if (indise /= 0) then
      disex = reada (keywrd, indise)
    end if
   ! FILL THE COSMO - RADII (SRAD) AND INDEX - VECTORS IDENAT AND IPIDEN
    iden = 0
    do i = 1, numat
      nfi = nfirst(i)
      idenat(i) = iden + 1
      idel = nlast(i) + 1 - nfi
      if (mozyme) then
        i0 = ijbo (i, i)
        do j = 1, idel
          nfj = nfi - 1 + j
          do k = 1, j
            iden = iden + 1
            i0 = i0 + 1
            ipiden(iden) = i0
            gden(iden) = -2.0d0
          end do
          gden(iden) = -1.0d0
        end do
      else
        do j = 1, idel
          nfj = nfi - 1 + j
          i0 = (nfj*(nfj - 1)) / 2 + nfi - 1
          do k = 1, j
            iden = iden + 1
            ipiden(iden) = i0 + k
            gden(iden) = -2.0d0
          end do
          gden(iden) = -1.0d0
        end do
      end if
      iat = nat(i)
      srad(i) = usevdw(iat)
    end do

   !
   !  NUMBER OF SECTIONS PER ATOM
   !
    n0(1) = nspa
    if (nspa /= 42) then
      x = (nspa - 2) / 10.d0 + 1d-8
      n1 = 10 * Int (Sqrt(x)) ** 2 + 2
      n2 = 30 * Int (Sqrt(x/3)) ** 2 + 2
      n0(1) = Max (n1, n2)
      n0(1) = Max (12, n0(1))
      if (nspa /= n0(1)) then
        write (iw, "(A, I6)") "NSPA IS SET TO", n0 (1)
      end if
    end if
    nspa = n0(1)
    n1 = (n0(1) - 2) / 10
   !
   !   FOR HYDROGEN, USE ABOUT 1/3 THE NUMBER OF ATOMS.
   !
    n0(2) = 10 * n1 / 3 + 2
    if (Mod(n1, 3) /= 0) then
      n0(2) = 10 * Int (Sqrt(n1/3.d0)) ** 2 + 2
    end if
   !
   !   RESET N0(2) TO 12.  N0 = 42, 12
   !
    n0(2) = Max (n0(2), 12)
    if (n0(1) + n0(2) > 1082) then
      call mopend ("CHOSE NSPA < 0.75*1082")
    end if
    call dvfill (n0(1), dirsm)
    call dvfill (n0(2), dirsm(1, n0(1) + 1))
    disex2 = 4 * (1.7d0*disex) ** 2 / nspa
    call dvfill (1082, dirvec)
end subroutine cosini
subroutine extvdw (vdw, refvdw)
    use molkst_C, only : numat, keywrd, line
    use chanel_C, only : iw
    use common_arrays_C, only : nat
    use elemts_C, only : cap_elemnt, elemnt
    implicit none
    double precision, dimension (107), intent (out) :: vdw
    double precision, dimension (107), intent (in) :: refvdw
    integer :: i, j, k
    double precision, external :: reada
   !
   ! ... Executable Statements ...
   !
   !  Modify Van der Waal's radii of various atoms.  The format of
   ! the keyword is: VDW(chemical symbol = n.nn:chemical symbol = n.nn ...)
   ! e.g. VDW(:H = 1.0:Cl = 1.7)
   !
   i = index(keywrd, " VDW(")
    if (i == 0) then
      line = " "
    else
      if (keywrd(i + 5:i + 5) /= ";") keywrd(i + 5:) = ";"//keywrd(i + 5:)
      j = index(keywrd(i + 5:), ")") + i + 5
      do k = i + 5, j
         if (keywrd(k:k) == ":") keywrd(k:k) = ";"
         if (keywrd(k:k) == ", ") keywrd(k:k) = ";"
      end do
      line = keywrd(i + 5:j - 1)
    end if
    vdw(:107) = refvdw(:107)
    if (line /= " ") then
      do i = 1, 107
        j = 2
        if (cap_elemnt(i) (2:2) == " ") j = 1
        j = index(line, ";"//cap_elemnt(i)(:j)//"=")
        if (j > 0) vdw(i) = reada(line, j)
      end do
    end if
    !
    !  Verify that all radii that will be used are, in fact, set correctly
    !
    do i = 1, numat
      if (nat(i) > 102) cycle
      if (vdw(nat(i)) > 900.d0) then
        write (line, *) "MISSING VAN DER WAALS RADIUS " // elemnt (nat(i))
        call mopend (trim(line))
        j = 2
        if (elemnt(nat(i)) (2:2) == " ") j = 1
        write (iw, "(2x, 3a)") "To correct this, add keyword 'VDW(", &
       & elemnt (nat(i)) (1:j), " = n.nn)'"
        return
      end if
    end do
end subroutine extvdw
