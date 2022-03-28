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

! ____________________________________________________________________________
! _________________________ OVERLAP ROUTINES _________________________________
! ____________________________________________________________________________

      subroutine ovlap (s, x, y, z)

      use reimers_C, only : n, na, natm, iat, nbt, nprn, nbeta, &
          zeta, zetad, zetawt, fssig, fpsig, fppi, fdsig, fdpi, &
          fddel, pp, bincoe, fact, d, e, fspdf, kflag, jflag
      USE molkst_C, ONLY: mpack, keywrd

      implicit none
      double precision ::  x(na), y(na), z(na), s(mpack), &
                        one, two, zero, &
                        g, abc, acu, amu, bcu, bmu, bin, &
                        edel, ephi, epi, esig, ggg, &
                        pig, delg, phig, &
                        reada, rij, sigg, ss, sp, sd, sf, &
                        sss, ssp, ssd, ssf, &
                        u1, u2, u3
      integer ::        langl(16), i, j, k, &
                        i1, ind, iss, issmax, &
                        ja, jb, jj, jl, jlp1, jm, jnat, jnp, jss, jssmax, &
                        ka, kb, kl, klp1, knat, knp, la, lh, minkm, nm
      character*9       atomk, atomkx, atomj, atomjx
      logical           agovlp
      double precision  agfact

      data one/1.d0/, two/2.d0/, zero/0.0d0/
      data langl/0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3/

      edel = 0.d0
      epi = 0.d0
      pig = 0.d0
      phig = 0.d0
      delg = 0.d0
      ss = 0.d0
      sp = 0.d0
      sd = 0.d0
      sf = 0.d0
      sigg = 0.d0
      jm = 0

! RMG - add scaling factor
      agovlp = .False.
      agfact = 0.d0
      if (index(keywrd, ' AGOVLP = ') /= 0) then
          agovlp = .True.
          agfact = reada(keywrd, index(keywrd, ' AGOVLP = '))
      end if

      if (nbeta == 4) then
!	**** Tomono betas, set overlap to nearest neighbour indicator ****
        nm = 0
        do i = 1, n
          do j = 1, i
            nm = nm + 1
            s(nm) = zero
            if (r(i, j) < 1.7d0) s(nm) = one
          end do
        end do
        return
      end if

!     *********************************************************
!     ******* general overlap code ****************************
!     *********************************************************

!     **** n - factorial ****
      fact(1) = 1
      do i = 2, 30
        fact(i) = fact(i - 1) * (i - 1)
      end do

!
!     ********* form binomial coeffecient *******
!
!      (n )
!      (m )  index = n*(n + 1)/2 + m + 1
!
      bincoe(1) = 1.d0
      ind = 2
      do 601 i = 1, 29
      i1 = i + 1
      do 602 j = 1, i1
      jj = i1 - j + 1
      if (j < 26)go to 66
      bin = 1.d0
      if (j >= jj)then
      do 62 k = j, i
      bin = bin*dfloat(k)
62    continue
      bincoe(ind) = bin/fact(jj)
      else
      do 64 k = jj, i
      bin = bin*dfloat(k)
64    continue
      bincoe(ind) = bin/fact(j)
      end if
      go to 60
66    bincoe(ind) = fact(i1)/fact(j)/fact(jj)
60    ind = ind + 1
602   continue
601   continue

!  introduce atomk and atomj to test if treating the same shell
!
!    kflag(1) - - > k atom no
!    kflag(2) - - > k pqn
!    kflag(3) - - > k lqn
!    jflag(1) - - > j atom no
!    jflag(2) - - > j pqn
!    jflag(3) - - > j lqn
!
!  initialise atomkx and atomjx to zero
      atomkx = '  0  0  0'
      atomjx = ' 0 0 0'
!
!     fspdf(i, j) is weighting factor.
!     i = angular quantum number + 1,
!     j = 1, 2, 3, 4 for sigma, pi, delta, phi.
      do 1 i = 1, 4
      do 11 j = 1, 4
      fspdf(i, j) = one
   11 continue
    1 continue
      fspdf(1, 1) = fssig
      fspdf(2, 1) = fpsig
      fspdf(2, 2) = fppi
      fspdf(3, 1) = fdsig
      fspdf(3, 2) = fdpi
      fspdf(3, 3) = fddel
!
!  set pp (in common pptt) to - 10.d0 so as to initialise ovlap
!
      pp = -10.d0
!
      nm = 0
      do 100 i = 1, n
        ka = iat(i)
        kb = nbt(i)
        knat = natm(ka)
        knp = nprn(i)
        kl = langl(kb + 1)
        klp1 = kl + 1
        kflag(1) = ka
        kflag(2) = knp
        kflag(3) = kl
!
!      do an internal write to convert the saved flag values to
!      character strings.
!
        write(atomk, '(3i3)') kflag
        do 90 j = 1, i
          nm = nm + 1
          if (j == i) then
            s(nm) = 1.d0
            goto 90
          end if
          ja = iat(j)
          jb = nbt(j)
          jnat = natm(ja)
          jnp = nprn(j)
          jl = langl(jb + 1)
          jlp1 = jl + 1
          jflag(1) = ja
          jflag(2) = jnp
          jflag(3) = jl
          write (atomj, '(3i3)') jflag
          rij = r(ja, ka)
!
!    one center overlaps
!
         if (ja == ka) then
            s(nm) = zero
            go to 90
         end if

         if (kb <= jb) then
            la = jb
            lh = kb
         else
           la = kb
           lh = jb
           rij = -rij
         end if

!	  **** geometrical matrices ****

          if (abs(rij) < 1.0D-4) then
            u1 = one
            u2 = one
            u3 = zero
          else
            u1 = (z(ja) - z(ka))/rij
            g = one - u1**2
            if (g  <=  1.D-7) then
              u2 = one
              if (u1 < zero) u2 = -one
              u3 = zero
            else
              g = one/dsqrt(g)
              u2 = (x(ja) - x(ka))*g/rij
              u3 = (y(ja) - y(ka))*g/rij
            end if
          end if

      call geome (u1, u2, u3, lh, d)
      call geome (u1, u2, u3, la, e)

!.......................................................
!
!     get phase for e transform
!     phase (- 1)**(l - m) on la >= the other l value
!
      if (la == jb) jm = jlp1
      if (la == kb) jm = klp1
      go to (27, 28, 29, 30), jm
!
!  f - orbs
!
   30 ephi = one
      edel = -one
!
!  p - orbs
!
   28 epi = one
      esig = -one
      go to 31
!
!  d-orbs
!
   29 edel = one
      epi = -one
!
!  s - orbs
!
   27 esig = one
   31 continue
!.......................................................
      rij = abs (rij)
      if (atomk == atomkx.and.atomj == atomjx) go to 46
      atomkx = atomk
      atomjx = atomj
      ss = zero
      sp = zero
      sd = zero
      sf = zero
!
!     do double zeta for each atom.  when exponent is zero, quit.
!
      issmax = 1
      if (nbt(i) > 3) issmax = 2
      jssmax = 1
      if (nbt(j) > 3) jssmax = 2

      do 38 iss = 1, issmax
         if (nbt(i) > 3)  then
           acu = zetawt(iss, knat)
           amu = zetad(iss, knat)
         else
           acu = 1.d0
           amu = zeta(knat)
         end if
         if (amu < 1.0d-4) go to 38

         do 32 jss = 1, jssmax
            if (nbt(j) > 3)  then
              bcu = zetawt(jss, jnat)
              bmu = zetad(jss, jnat)
            else
              bcu = 1.d0
              bmu = zeta(jnat)
            end if
            if (bmu < 1.0d-4) go to 32

            abc = acu*bcu
            call ovlaap(knp, kl, amu, jnp, jl, bmu, rij, sss, ssp, ssd, ssf)
            ss = ss + sss*abc
            sp = sp + ssp*abc
            sd = sd + ssd*abc
            sf = sf + ssf*abc
   32       continue
   38    continue

!     **** this is some type of damping for low r, uses ZINDOs alt single STO set
      ggg = 0.d0
      if (kl /= jl) then
!	**** set interaction factors to one for cross sp, sd, pd interactions ****
        sigg = one
        pig = one
        delg = one
        phig = one
      else
!	**** ss, pp, or dd interaction, damped ****
        sigg = (fspdf(klp1, 1) + fspdf(jlp1, 1))/two
        pig = (fspdf(klp1, 2) + fspdf(jlp1, 2))/two
        delg = (fspdf(klp1, 3) + fspdf(jlp1, 3))/two
        phig = (fspdf(klp1, 4) + fspdf(jlp1, 4))/two
        sigg = sigg + (one - sigg)*ggg
        pig = pig + (one - pig )*ggg
        delg = delg + (one - delg)*ggg
        phig = phig + (one - phig)*ggg
      end if

!.......................................................
!
!  form sigma overlap
!
 46   s(nm) = e(1)*d(1)*ss*esig*sigg
!
!  check how many terms in overlap
!
      minkm = min0(kl, jl) + 1
      go to (50, 47, 48, 49), minkm
!
!  phi overlap
!
   49 s(nm) = s(nm) + (e(6)*d(6) + e(7)*d(7))*sf*ephi*phig
!
!  delta overlap
!
   48 s(nm) = s(nm) + (e(4)*d(4) + e(5)*d(5))*sd*edel*delg
!
!  pi overlap
!
   47 s(nm) = s(nm) + (e(2)*d(2) + e(3)*d(3))*sp*epi*pig
!
!.......................................................
   50 continue

! RMG - Ag overlap scaling
      if (agovlp) then
        if ((knat == 47.or.jnat == 47).and.knat /= jnat) then
          s(nm) = s(nm)*agfact
        end if
      end if

   90 continue
  100 continue

      return

      contains
            double precision function r(i, j)
                  implicit none
                  integer :: i, j
                  r = sqrt((x(j) - x(i))**2 + (y(j) - y(i))**2 + (z(j) - z(i))**2 )
            end function r

      end subroutine ovlap

!     *************************************************************************

      subroutine ovlaap (n1, l1, amu, n2, l2, bmu, r, s, sp, sd, sf)
!
!     evaluates overlap between two slater type orbitals.
!     n1, n2, are princ. q.no. - s, l1, l2, are secondary q.n. - s, , amu and
!     bmu are the exponential constants. r is the separation in angstrom
!     s is the sigma, sp the pi, sd the delta and sf the f compoents of
!     the overlap, to be put together using lh and sub. geom.
!     fact(i + 1) = factorial i.
!
      use reimers_C, only : au2ang, pp, fact, ss

      implicit none
      double precision ::  zero, one, two, &
                        f, p, r, s, t, &
                        amu, bmu, rr, sp, sd, sf, ssss, tt
      integer ::           i, ii, l1, l2, lmin, mm, n1, n2, nn
      save :: tt

      data              zero/0.d0/, one/1.d0/, two/2.d0/

      if (abs(pp + 10.d0) < 1.d-10)then
      tt = zero
      end if
      rr = r/au2ang
      p = (amu + bmu)*rr/two
      t = (amu - bmu)/(amu + bmu)
      if (r <= 1.d-3) go to 300
      if(dabs (t) <= 1.d-4) t = zero
      if(abs(p - pp) < 1.d-10 .and. abs(t - tt) < 1.d-10) go to 5
      pp = p
      tt = t
      call aux (p, t)
    5 continue
      ss(2) = zero
      ss(3) = zero
      ss(4) = zero
      lmin = l1
      if (l2 < l1) lmin = l2
      lmin = lmin + 1
      do 10 ii = 1, lmin
      i = ii - 1
      call molpab (n1, n2, l1, l2, i, i, amu, bmu, rr, ssss)
      ss(ii) = ssss
   10 continue
      s = ss(1)
      sp = ss(2)
      sd = ss(3)
      sf = ss(4)
      go to 500
!     one center overlap.
  300 s = zero
      sp = zero
      sd = zero
      sf = zero
      if (l1 /= l2) go to 500
      nn = 2*n1 + 1
      mm = 2*n2 + 1
      f = (one - t)**mm*(one + t)**nn
      s = fact(n1 + n2 + 1)*dsqrt(f/(fact(nn)*fact(mm)))
      if (l1 > 0) sp = (- one)**(l1 + 1)*s
      if (l1 > 1) sd = (- one)**(l1 + 2)*s
      if (l1 > 2) sf = (- one)**(l1 + 3)*s
      s = (- one)**l1*s
  500 continue
      return
      end subroutine ovlaap

!     *************************************************************************

      subroutine aux (pp, tt)

!     this subroutine calculates a and b fns. for molecular integrals.
!     b - fns. of argument pp*tt and with indices of 0 to ix, a fns. of
!     argument pp and indices from 0 to ix.
! **  references:
! **  (a) mullikan/rieke/orloff/orloff (1949) jcp 17, 1248;
! **      pp (= p or = rho) = 1/2*[mu(a) + mu(b)]*r/a(h)
! **      tt (= t or = tau) = [mu(a) - mu(b)]/[mu(a) + mu(b)]
! **  (b) kotani/amemiya/simose (1940) proc. phys. mat. soc. japan 22, 1

      use reimers_C, only : a, b
      implicit none
!      implicit          double precision :: (a - h, o - z)
      double precision ::  pp, tt, upplim, zero, one, two, &
                        c, d, h, r, t, ra, rho2, tr
      integer ::           i, j, k, il, in, ir, is, ix, ixm, ixs, jy

      data upplim       /1.0D36/
      data zero, one, two /0.0d0, 1.0d0, 2.d0/

        ix = 32
        rho2 = pp*tt
        ir = int(abs(rho2 + rho2))
        if (ir > 170) go to 100
        is = min0(ir + 1, 15)
!
!     first the more difficult b - fn.
!
        if (dabs(rho2) <  1.d-3) go to 35
           if (pp > 40.0d0) go to 100
              d = dexp(rho2)
              h = dexp(- rho2)
              r = d-h
!
!     if the value of rho is too small the sinh must be obtained by
!     summing the infinite series rather than by addition of two
!     exponentials
!
        if (dabs(r) - 1.d-1 >= 0.d0) go to 28
           ra = zero
           t = rho2
              do 27 i = 2, 25
              if (dabs(t) < 1.d-18) go to 127
              t = t*rho2*rho2/dfloat((i + i - 1)*(i + i - 2))
              ra = ra + t
27            continue
127     continue
        r = (ra + rho2)*two
28      b(1) = r/rho2
!
!     as many successive b functions are generated from b(0) by the
!     recursion formula as accuracy will permit.
!
        ixs = ix
            do 51 i = 2, ixs, is
            if (ir == 0) go to 40

            il = is - 1
            if ((i + il) > ix) il = ix - i + 1
              do 31 j = 1, il
              k = i + j - 1
                 if (mod(k, 2) <= 0) go to 30
                 b(k) = (r + dfloat(k - 1)*b(k - 1))/rho2
                 go to 31
30               b(k) = -(d + h - dfloat(k - 1)*b(k - 1))/rho2
31               continue
40               in = i + is - 1
                 if (in - ix <= 0) go to 39
                 if (in - ix > 0) go to 38
!
!     after the recurrence formula has been applied an appropriate no.
!     of times the next b function is obtained by summing the infinite
!     series.
!
39      if (mod(in, 2) > 0) go to 44
        tr = rho2
        b(in) = -(tr + tr)/dfloat(in + 1)
          do 43 j = 1, 500
          tr = tr*rho2*rho2/dfloat((j + j)*(j + j + 1))
!      note accuracy criterion
!      if (dabs(tr/b(in)) - 0.000011 ) 51, 51, 43
          if (dabs(tr) - 1.d-7*dabs(b(in)) <= 0.d0) go to 51
          b(in) = b(in) - (tr + tr)/dfloat(in + 1 + j + j)
43        continue
44        tr = one
          b(in) = (tr + tr)/dfloat(in)
          do 46 j = 1, 500
          tr = tr*rho2*rho2/dfloat((j + j)*(j + j - 1))
!      note accuracy criterion
!      if (dabs(tr/b(in)) - 0.00001  ) 51, 51, 46
          if (dabs(tr) - 1.d-7*dabs(b(in)) <= 0.d0) go to 51
            b(in) = b(in) + (tr + tr)/dfloat(in + j + j)
46        continue
51        continue
        go to 38
!
!     if the argument of the b - fn is zero a separate formula is used.
!
35      continue
        jy = ix/2
           do 36 i = 1, jy
              b(i + i) = zero
           b(i + i - 1) = two/dfloat(i + i - 1)
36      continue
38      continue
!
!     now the a - functions
!
        c = dexp(- pp)
        a(1) = c/pp
        ixm = ix - 3
            do 15 i = 2, ixm
              if (a(i - 1) >= upplim) then
                  a(i) = a(i - 1)
                     go to 15
              end if
              a(i) = (dfloat(i - 1)*a(i - 1) + c)/pp
15          continue
            go to 500
100     continue
          do 120 i = 1, ix
             a(i) = 0.0d0
             b(i) = 0.0d0
120       continue
500       return
          end subroutine aux

!     *************************************************************************

       subroutine molpab (n1, n2, l1, l2, m1, m2, sk1, sk2, r, vest)
!
!     this subroutine calculates two centered overlap integrals
!
      use reimers_C, only : fact

!      implicit          double precision :: (a - h, o - z)
      implicit none
      double precision ::   sk1, sk2, r, vest, &
                         f1, f2, f3, f4, f5, f6, &
                         f11, f12, f13, f14, f15, f16, f17, f18, &
                         q1, q2, rhoa, rhob, strad, strad1, term, value
      integer ::            j, k, iff, jend, ju, kend, ku, &
                         l1, l2, m1, m2, n1, n2

      vest = 0.0d0
      if (m1 /= m2)  go to 500
      strad = 0.0d0
      rhoa = r*sk1
      rhob = r*sk2
      f1 = fact(2*n1 + 1)
      f2 = fact(2*n2 + 1)
      iff = l1 - m1 + 1
      f3 = fact(iff)
      iff = l2 - m2 + 1
      f4 = fact(iff)
      iff = l1 + m1 + 1
      f5 = fact(iff)
      iff = l2 + m2 + 1
      f6 = fact(iff)
      q1 = (sk1**n1/dsqrt(f1))*r**n1
      q2 = (sk2**n2/dsqrt(f2))*r**n2
      term = dsqrt(dfloat((2*l1 + 1)*(2*l2 + 1))*&
     &  (f3/f5)*(f4/f6)*rhoa*rhob)/2.d0**(l1 + l2 + 1)
      term = q1*term*q2
      jend = 1 + ((l1 - m1)/2)
      kend = 1 + ((l2 - m2)/2)
      do 501 j = 1, jend
      ju = j - 1
      iff = 2*l1 - 2*ju + 1
      f11 = fact(iff)
      iff = l1 - m1 - 2*ju + 1
      f13 = fact(iff)
      f15 = fact(ju + 1)
      iff = l1 - ju + 1
      f17 = fact(iff)
      do 502 k = 1, kend
      ku = k - 1
      iff = 2*l2 - 2*ku + 1
      f12 = fact(iff)
      iff = l2 - m2 - 2*ku + 1
      f14 = fact(iff)
      f16 = fact(ku + 1)
      iff = l2 - ku + 1
      f18 = fact(iff)
      call cfunct (n1 - l1 + 2*ju, n2 - l2 + 2*ku, l1 - m1 - 2*ju, l2 - m2 - 2*ku, m1, value)
      strad1 = value * (f11/f13/f15/f17) * (f12/f14/f16/f18)
      if (mod(ju + ku, 2) == 1) strad1 = -strad1
      strad = strad + strad1
  502 continue
  501 continue
      vest = term*strad
  500 return
      end subroutine molpab

!     *************************************************************************

      subroutine cfunct (ia, ib, ic, id, ie, snag)

      use reimers_C, only : bincoe, a, b
!      implicit          double precision :: (a - h, o - z)
      implicit none
      double precision ::   count, b1, b2, b3, b4, b5, b6, snag, term
      integer ::            i, i1, i2, i3, i4, i5, i6, &
                         ia, ib, ic, id, ie, ip, ir, &
                         iab, ibb, icb, idb, ieb, &
                         inda, indb, indc, indd, inde, nin

      nin(i) = i * (i - 1) / 2

      count = 0.0d0
      iab = ia + 1
      ibb = ib + 1
      icb = ic + 1
      idb = id + 1
      ieb = ie + 1
      inde = nin(ieb)
      indd = nin(idb)
      indc = nin(icb)
      indb = nin(ibb)
      inda = nin(iab)
      do 901 i6 = 1, ieb
      b6 = bincoe(inde + i6)
      do 902 i5 = 1, ieb
      b5 = bincoe(inde + i5)
      do 903 i4 = 1, idb
      b4 = bincoe(indd + i4)
      do 904 i3 = 1, icb
      b3 = bincoe(indc + i3)
      do 905 i2 = 1, ibb
      b2 = bincoe(indb + i2)
      do 906 i1 = 1, iab
      b1 = bincoe(inda + i1)
      term = b1*b2*b3*b4*b5*b6
      if (mod(i2 + i5 + i6 + i4 + ie + id, 2) == 1) term = -term
      ir = i1 + i2 - i3 - i4 + ie + ie - i6 - i6 + ic + id + 3
      ip = ia - i1 + ib - i2 + ie + ie - i5 - i5 + ic - i3 + id-i4 + 7
      count = count + a(ip)*b(ir)*term
  906 continue
  905 continue
  904 continue
  903 continue
  902 continue
  901 continue
      snag = count
      end subroutine cfunct

!     *************************************************************************

      subroutine geome (ct, cp, sp, m, h)
!
!     this subroutine gives the transformation matrix necessary for the
!     expression of p, d, and f tensors in the simplest canonical terms.
!
!     m = 0, 1, 2, 3, 4, 5, 6, 7, 8 for s, x, y, z, 3z2 - r2, x2 - y2, xy, xz, yz respectively
!     ____________ deleted to save space ______________
!     m = 9, 10, 11, 12, 13, 14, 15 for z3, xz2, yz2, z(x2 - y2), xyz, x(x2 - 3y2), and
!        y(3x2 - y2)
!     m = 51, 52, 53 for x2, y2, z2
!
!     h(i) returns as sigma, pi - x, pi - y, delta - (x**2 - y**2), delta - xy, phi -
!          x(x2 - 3y2), phi - y(3x2 - y2), and s component if any.
!
!     the phases of the real fns. are those of frank harris in the upps
!     ala monograph.
!     spherical coordinates. ct = cos(theta), cp = cos(phi), sp = sin(phi).
!     h returns as trans elements of 3 - dim. space for p orb. and of 5 - di
!     m space for the d orb., etc.

!      implicit          double precision :: (a - h, o - z)
      implicit none
      double precision ::   h(8), ct, cp, sp, dsqrt3, one, two, three, zero, &
                         r, st, z
      integer ::            m

      data dsqrt3/1.732050807d0/one/1.d0/zero/0.d0/two/2.d0/three/3.d0/

      h(1) = one
      h(2) = zero
      h(3) = zero
      h(4) = zero
      h(5) = zero
      h(6) = zero
      h(7) = zero
      h(8) = zero
      if (m == 0) go to 50
      continue
      r = two*sp*cp
      st = dabs (one - ct**2)
      if (1.D-7 - st <= 0.d0) go to 22
      st = zero
      go to 24
   22 continue
      st = dsqrt(st)
   24 continue
      z = cp**2 - sp**2
      go to (1, 2, 3, 4, 5, 6, 7, 8), m
!     p type fn.
    1 h(1) = cp*st
      h(2) = cp*ct
      h(3) = -sp
      go to 50
    2 h(1) = sp*st
      h(2) = sp*ct
      h(3) = cp
      go to 50
    3 h(1) = ct
      h(2) = -st
      h(3) = zero
      go to 50
!     d type fn .
    4 h(1) = (three*ct*ct - one)/two
      h(2) = -dsqrt3*ct*st
      h(3) = zero
      h(4) = dsqrt3*(one - ct**2)/two
      h(5) = zero
      go to 50
    5 h(1) = dsqrt3*z*(one - ct*ct)/two
      h(2) = ct*z*st
      h(3) = -st*r
      h(4) = z*(one + ct*ct)/two
      h(5) = -r*ct
      go to 50
    6 h(1) = dsqrt3*cp*sp*(one - ct*ct)
      h(2) = ct*st*r
      h(3) = z*st
      h(4) = cp*sp*(one + ct*ct)
      h(5) = ct*z
      go to 50
    7 h(1) = dsqrt3*cp*ct*st
      h(2) = cp*(two*ct*ct - one)
      h(3) = -sp*ct
      h(4) = -h(1)/dsqrt3
      h(5) = sp*st
      go to 50
    8 h(1) = dsqrt3*sp*ct*st
      h(2) = sp*(two*ct*ct - one)
      h(3) = cp*ct
      h(4) = -h(1)/dsqrt3
      h(5) = -cp*st
   50 continue
      return
      end subroutine geome
