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

      subroutine analyt(psum, palpha, pbeta, coord, nat, jja, jjd, iia, iid, &
        eng) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE funcon_C, only : a0, fpc_9
      USE molkst_C, only : keywrd, numcal, mpack
      USE analyt_C, only : ds, dg, dr, g, nztype
      USE parameters_C, only : alp, natorb, tore, guess1, guess2, guess3, &
      betas, betap 
      use chanel_C, only : irot
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: jja 
      integer , intent(in) :: jjd 
      integer , intent(in) :: iia 
      integer , intent(in) :: iid 
      integer , intent(in) :: nat(2) 
      double precision , intent(in) :: psum(mpack) 
      double precision , intent(in) :: palpha(mpack) 
      double precision , intent(in) :: pbeta(mpack) 
      double precision  :: coord(3,2) 
      double precision , intent(inout) :: eng(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icalcn, jd, ja, id, ia, j, i, ni, istart, nj, jstart&
        , ix, isp, iol, k, ka, kg, l, la, lg, is, ig, m, mn&
        , n, kk, ll, kl, mk, nk, ml, nl, i22
      double precision, dimension(3) :: eaa, eab, enuc 
      double precision, dimension(4) :: bi, bj 
      double precision ::  r2, rij, r0, rr, del1, termaa, termab&
        , del2, del3, c1, termnc, &
        part1, part2, part3, f3, dd, anam1, termk, bb, aa 
      logical :: am1

      save icalcn, am1
!-----------------------------------------------
!***********************************************************************
!                                                                      *
!         CALCULATION OF ANALYTICAL DERIVATIVES                        *
!                                                                      *
!***********************************************************************
      data icalcn/ 0/   
      if (icalcn /= numcal) then 
        icalcn = numcal 
        am1 = index(keywrd,'AM1') + index(keywrd,'PM3') /= 0 
      end if 
      jd = jjd - jja + 1 
      ja = 1 
      id = iid - iia + 1 + jd 
      ia = jd + 1 
      eaa = 0.0D0 
      eab = 0.0D0 
      enuc = 0.0D0 
      eng = 0.0D0 
      i = 2 
      ni = nat(i) 
      istart = nztype(ni)*4 - 3 
      j = 1 
      nj = nat(j) 
      jstart = nztype(nj)*4 - 3 
      r2 = (coord(1,i)-coord(1,j))**2 + (coord(2,i)-coord(2,j))**2 + (coord(3,i&
        )-coord(3,j))**2 
      rij = sqrt(r2) 
      r0 = rij/a0 
      rr = r2/(a0*a0) 
      do ix = 1, 3 
        del1 = coord(ix,i) - coord(ix,j) 
        termaa = 0.0D0 
        termab = 0.0D0 
        isp = 0 
        iol = 0 
!   THE FIRST DERIVATIVES OF OVERLAP INTEGRALS
        do k = ia, id 
          ka = k - ia 
          kg = istart + ka 
          do l = ja, jd 
            la = l - ja 
            lg = jstart + la 
            iol = iol + 1 
            ds(iol) = 0.0D0 
            if (ka==0 .and. la==0) then 
!   (S/S) TERM
              if (abs(del1) <= 1.0D-6) cycle  
              is = 1 
            else if (ka==0 .and. la>0) then 
!   (S/P) TERM
              is = 3 
              if (ix == la) go to 20 
              if (abs(del1) <= 1.0D-6) cycle  
              is = 2 
              del2 = coord(la,i) - coord(la,j) 
            else if (ka>0 .and. la==0) then 
!   (P/S) TERM
              is = 5 
              if (ix == ka) go to 20 
              if (abs(del1) <= 1.0D-6) cycle  
              is = 4 
              del2 = coord(ka,i) - coord(ka,j) 
            else 
!   (P/P) TERM
              if (ka == la) then 
!    P/P
                is = 9 
                if (ix == ka) go to 20 
                if (abs(del1) <= 1.0D-6) cycle  
!    P'/P'
                is = 8 
                del2 = coord(ka,i) - coord(ka,j) 
              else if (ix/=ka .and. ix/=la) then 
!    P'/P"
                if (abs(del1) <= 1.0D-6) cycle  
                is = 7 
                del2 = coord(ka,i) - coord(ka,j) 
                del3 = coord(la,i) - coord(la,j) 
              else 
!    P/P' OR P'/P
                del2 = coord(ka+la-ix,i) - coord(ka+la-ix,j) 
                is = 6 
              end if 
            end if 
!
!        CALCULATE OVERLAP DERIVATIVES, STORE RESULTS IN DS
!
   20       continue 
            call ders (kg, lg, rr, del1, del2, del3, is, iol) 
          end do 
        end do 
        if (ix==1) read (irot) (g(i22),i22=1,22) 
        call delri (dg, ni, nj, r0, del1) 
        call delmol (coord, i, j, ni, nj, ia, id, ja, jd, ix, rij, del1, isp)   
!
!   THE FIRST DERIVATIVE OF NUCLEAR REPULSION TERM
!
!      CORE-CORE TERMS, MNDO AND AM1
!
!
!  SPECIAL TREATMENT FOR N-H AND O-H TERMS
!
          if (rij<1.D0 .and. natorb(ni)*natorb(nj)==0) then 
            termnc = 0.D0 
            go to 50 
          end if 
          c1 = tore(ni)*tore(nj) 
          if (ni==1 .and. (nj==7 .or. nj==8)) then 
            f3 = 1.0D0 + exp((-alp(1)*rij)) + rij*exp((-alp(nj)*rij)) 
            dd = (dg(1)*f3-g(1)*(del1/rij)*(alp(1)*exp((-alp(1)*rij))+(alp(&
              nj)*rij-1.0D0)*exp((-alp(nj)*rij))))*c1 
          else if ((ni==7 .or. ni==8) .and. nj==1) then 
            f3 = 1.0D0 + exp((-alp(1)*rij)) + rij*exp((-alp(ni)*rij)) 
            dd = (dg(1)*f3-g(1)*(del1/rij)*(alp(1)*exp((-alp(1)*rij))+(alp(&
              ni)*rij-1.0D0)*exp((-alp(ni)*rij))))*c1 
          else 
!
!  SPECIAL CASE OF TWO SPARKLES
!
            part1 = dg(1)*c1 
            part2 = -(g(1)*(del1/rij)*(alp(ni)*exp((-alp(ni)*rij))+alp(nj)*&
              exp((-alp(nj)*rij))))*abs(c1) 
            part3 = dg(1)*(exp((-alp(ni)*rij))+exp((-alp(nj)*rij)))*abs(c1) 
            dd = part1 + part2 + part3 
!
!   THE GENERAL CASE
!
          end if 
          termnc = dd 
!
!   ****   START OF THE AM1 SPECIFIC DERIVATIVE CODE   ***
!
!      ANALYT=-A*(1/(R*R)+2.D0*B*(R-C)/R)*EXP(-B*(R-C)**2)
!
!    ANALYTICAL DERIVATIVES
!
        if (am1) then 
          anam1 = 0.D0 
          do ig = 1, 4 
            if (abs(guess1(ni,ig)) > 0.D0) anam1 = anam1 + guess1(ni,ig)*(1.D0/(rij*&
              rij) + 2.D0*guess2(ni,ig)*(rij-guess3(ni,ig))/rij)*exp(max(-30.D0,(-guess2&
              (ni,ig)*(rij-guess3(ni,ig))**2))) 
            if (abs(guess1(nj,ig)) <= 0.D0) cycle  
            anam1 = anam1 + guess1(nj,ig)*(1.D0/(rij*rij) + 2.D0*guess2(nj,ig)*(rij-&
              guess3(nj,ig))/rij)*exp(max(-30.D0,(-guess2(nj,ig)*(rij-guess3(nj,ig))**2)&
              )) 
          end do 
          anam1 = anam1*tore(ni)*tore(nj) 
          termnc = termnc - anam1*del1/rij 
        end if 
!
!   ****   END OF THE AM1 SPECIFIC DERIVATIVE CODE   ***
!
   50   continue 
!
!   COMBINE TOGETHER THE OVERLAP DERIVATIVE PARTS
!
          bi(1) = betas(ni) 
          bi(2) = betap(ni) 
          bi(3) = bi(2) 
          bi(4) = bi(2) 
          bj(1) = betas(nj) 
          bj(2) = betap(nj) 
          bj(3) = bj(2) 
          bj(4) = bj(2) 
        iol = 0 
        do k = ia, id 
          if (jd - ja + 1 > 0) then 
            termk = bi(k-ia+1) 
            termab = termab + sum((termk + bj(:jd-ja+1))*psum(ja+k*(k-1)/2:jd+k&
              *(k-1)/2)*ds(iol+1:jd-ja+1+iol)) 
            iol = jd - ja + 1 + iol 
          end if 
        end do 
!
!        FIRST, CORE-ELECTRON ATTRACTION DERIVATIVES (MNDO AND AM1)
!
!          ATOM CORE I AFFECTING A.O.S ON J
          isp = 0 
          do m = ja, jd 
            bb = 1.D0 
            do n = m, jd 
              mn = m + (n*(n - 1))/2 
              isp = isp + 1 
              termab = termab - bb*tore(ni)*psum(mn)*dr(isp) 
              bb = 2.D0 
            end do 
          end do 
!          ATOM CORE J AFFECTING A.O.S ON I
          k = max(jd - ja + 1,1) 
          k = (k*(k + 1))/2 
          isp = (-k) + 1 
          do m = ia, id 
            bb = 1.D0 
            do n = m, id 
              mn = m + (n*(n - 1))/2 
              isp = isp + k 
              termab = termab - bb*tore(nj)*psum(mn)*dr(isp) 
              bb = 2.D0 
            end do 
          end do 
          isp = 0 
!
!   NOW FOR COULOMB AND EXCHANGE TERMS (MNDO AND AM1)
!
          do k = ia, id 
            aa = 1.D0 
            kk = (k*(k - 1))/2 
            do l = k, id 
              ll = (l*(l - 1))/2 
              do m = ja, jd 
                bb = 1.D0 
                do n = m, jd 
                  isp = isp + 1 
                  kl = k + ll 
                  mn = m + (n*(n - 1))/2 
!
!    COULOMB TERM
!
                  termaa = termaa + aa*bb*psum(kl)*psum(mn)*dr(isp) 
                  mk = m + kk 
                  nk = n + kk 
                  ml = m + ll 
                  nl = n + ll 
!
!    EXCHANGE TERM
!
                  termaa = termaa - 0.5D0*aa*bb*(palpha(mk)*palpha(nl)+palpha(&
                    nk)*palpha(ml)+pbeta(mk)*pbeta(nl)+pbeta(nk)*pbeta(ml))*dr(&
                    isp) 
                  bb = 2.D0 
                end do 
              end do 
              aa = 2.D0 
            end do 
          end do 
        eaa(ix) = eaa(ix) + termaa 
        eab(ix) = eab(ix) + termab 
        enuc(ix) = enuc(ix) + termnc 
      end do 
      eng = eaa + eab + enuc 
      eng = -eng*fpc_9 
      return  
      end subroutine analyt 


      subroutine delmol(coord, i, j, ni, nj, ia, id, ja, jd, ix, rij, tomb, isp&
        ) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE analyt_C, only : dr, dg, tdx, tdy, tdz, g, tx, ty, tz
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: i 
      integer  :: j 
      integer , intent(in) :: ni 
      integer , intent(in) :: nj 
      integer , intent(in) :: ia 
      integer , intent(in) :: id 
      integer , intent(in) :: ja 
      integer , intent(in) :: jd 
      integer  :: ix 
      integer , intent(inout) :: isp 
      double precision  :: rij 
      double precision  :: tomb 
      double precision  :: coord(3,*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ib, jb, k, kk, l, ll, m, mm, n, nn 
      double precision :: temp1, temp2 
!-----------------------------------------------
      if (ni>1 .or. nj>1) call rotat (coord, i, j, ix, rij, tomb, 2) 
      ib = max0(ia,id) 
      jb = max0(ja,jd) 
   do k = ia, ib
      kk = k - ia
      do l = k, ib
        ll = l - ia
        do m = ja, jb
          mm = m - ja
          do n = m, jb
            nn = n - ja
            isp = isp + 1
            if (nn == 0) then
              if (ll == 0) then
                     !   (SS/SS)
                dr(isp) = dg(1)
              else if (kk == 0) then
                     !   (SP/SS)
                dr(isp) = dg(2) * tx(ll) + g(2) * tdx(ll)
              else
                     !   (PP/SS)
                dr(isp) = dg(3) * tx(kk) * tx(ll) + g(3) * &
               & (tdx(kk)*tx(ll)+tx(kk)*tdx(ll)) + dg(4) * &
               & (ty(kk)*ty(ll)+tz(kk)*tz(ll)) + g(4) * &
               & (tdy(kk)*ty(ll)+ty(kk)*tdy(ll)+tdz(kk)*tz(ll)+tz(kk)*tdz(ll))
              end if
            else if (mm == 0) then
              if (ll == 0) then
                     !   (SS/SP)
                dr(isp) = dg(5) * tx(nn) + g(5) * tdx(nn)
              else if (kk == 0) then
                     !   (SP/SP)
                dr(isp) = dg(6) * tx(ll) * tx(nn) + g(6) * &
               & (tdx(ll)*tx(nn)+tx(ll)*tdx(nn)) + dg(7) * &
               & (ty(ll)*ty(nn)+tz(ll)*tz(nn)) + g(7) * &
               & (tdy(ll)*ty(nn)+ty(ll)*tdy(nn)+tdz(ll)*tz(nn)+tz(ll)*tdz(nn))
              else
                     !   (PP/SP)
                dr(isp) = dg(8) * tx(kk) * tx(ll) * tx(nn) + g(8) * &
               & (tdx(kk)*tx(ll)*tx(nn)+tx(kk)*tdx(ll)*tx(nn)+&
               & tx(kk)*tx(ll)*tdx(nn)) + dg(9) * (ty(kk)*ty(ll)+&
               & tz(kk)*tz(ll)) * tx(nn) + g(9) * ((tdy(kk)*ty(ll)+&
               & ty(kk)*tdy(ll)+tdz(kk)*tz(ll)+tz(kk)*tdz(ll))*tx(nn)+&
               & (ty(kk)*ty(ll)+tz(kk)*tz(ll))*tdx(nn)) + dg(10) * &
               & (tx(kk)*(ty(ll)*ty(nn)+tz(ll)*tz(nn))+tx(ll)*(ty(kk)*ty(nn)+&
               & tz(kk)*tz(nn))) + g(10) * (tdx(kk)*(ty(ll)*ty(nn)+&
               & tz(ll)*tz(nn))+tdx(ll)*(ty(kk)*ty(nn)+tz(kk)*tz(nn))+&
               & tx(kk)*(tdy(ll)*ty(nn)+ty(ll)*tdy(nn)+tdz(ll)*tz(nn)+&
               & tz(ll)*tdz(nn))+tx(ll)*(tdy(kk)*ty(nn)+ty(kk)*tdy(nn)+&
               & tdz(kk)*tz(nn)+tz(kk)*tdz(nn)))
              end if
            else if (ll == 0) then
                  !   (SS/PP)
              dr(isp) = dg(11) * tx(mm) * tx(nn) + g(11) * &
             & (tdx(mm)*tx(nn)+tx(mm)*tdx(nn)) + dg(12) * &
             & (ty(mm)*ty(nn)+tz(mm)*tz(nn)) + g(12) * &
             & (tdy(mm)*ty(nn)+ty(mm)*tdy(nn)+tdz(mm)*tz(nn)+tz(mm)*tdz(nn))
            else if (kk == 0) then
                  !   (SP/PP)
              dr(isp) = dg(13) * tx(ll) * tx(mm) * tx(nn) + g(13) * &
             & (tdx(ll)*tx(mm)*tx(nn)+tx(ll)*tdx(mm)*tx(nn)+&
             & tx(ll)*tx(mm)*tdx(nn)) + dg(14) * tx(ll) * (ty(mm)*ty(nn)+&
             & tz(mm)*tz(nn)) + g(14) * (tdx(ll)*(ty(mm)*ty(nn)+&
             & tz(mm)*tz(nn))+tx(ll)*(tdy(mm)*ty(nn)+ty(mm)*tdy(nn)+&
             & tdz(mm)*tz(nn)+tz(mm)*tdz(nn))) + dg(15) * &
             & (ty(ll)*(ty(mm)*tx(nn)+ty(nn)*tx(mm))+tz(ll)*(tz(mm)*tx(nn)+&
             & tz(nn)*tx(mm))) + g(15) * (tdy(ll)*(ty(mm)*tx(nn)+&
             & ty(nn)*tx(mm))+tdz(ll)*(tz(mm)*tx(nn)+tz(nn)*tx(mm))+&
             & ty(ll)*(tdy(mm)*tx(nn)+ty(mm)*tdx(nn)+tdy(nn)*tx(mm)+&
             & ty(nn)*tdx(mm))+tz(ll)*(tdz(mm)*tx(nn)+tz(mm)*tdx(nn)+&
             & tdz(nn)*tx(mm)+tz(nn)*tdx(mm)))
            else
                  !   (PP/PP)
              dr(isp) = dg(16) * tx(kk) * tx(ll) * tx(mm) * tx(nn) + g &
             & (16) * (tdx(kk)*tx(ll)*tx(mm)*tx(nn)+&
             & tx(kk)*tdx(ll)*tx(mm)*tx(nn)+tx(kk)*tx(ll)*tdx(mm)*tx(nn)+&
             & tx(kk)*tx(ll)*tx(mm)*tdx(nn)) + dg(17) * (ty(kk)*ty(ll)+&
             & tz(kk)*tz(ll)) * tx(mm) * tx(nn) + g(17) * ((tdy(kk)*ty(ll)+&
             & ty(kk)*tdy(ll)+tdz(kk)*tz(ll)+tz(kk)*tdz(ll))*tx(mm)*tx(nn)+&
             & (ty(kk)*ty(ll)+tz(kk)*tz(ll))*(tdx(mm)*tx(nn)+tx(mm)*tdx(nn))) +&
             &  dg(18) * tx(kk) * tx(ll) * (ty(mm)*ty(nn)+tz(mm)*tz(nn)) + &
             & g(18) * ((tdx(kk)*tx(ll)+tx(kk)*tdx(ll))*(ty(mm)*ty(nn)+&
             & tz(mm)*tz(nn))+tx(kk)*tx(ll)*(tdy(mm)*ty(nn)+ty(mm)*tdy(nn)+&
             & tdz(mm)*tz(nn)+tz(mm)*tdz(nn)))
              dr(isp) = dr(isp) + dg(19) * (ty(kk)*ty(ll)*ty(mm)*ty(nn)+&
             & tz(kk)*tz(ll)*tz(mm)*tz(nn)) + g(19) * &
             & (tdy(kk)*ty(ll)*ty(mm)*ty(nn)+ty(kk)*tdy(ll)*ty(mm)*ty(nn)+&
             & ty(kk)*ty(ll)*tdy(mm)*ty(nn)+ty(kk)*ty(ll)*ty(mm)*tdy(nn)+&
             & tdz(kk)*tz(ll)*tz(mm)*tz(nn)+tz(kk)*tdz(ll)*tz(mm)*tz(nn)+&
             & tz(kk)*tz(ll)*tdz(mm)*tz(nn)+tz(kk)*tz(ll)*tz(mm)*tdz(nn)) + dg &
             & (20) * (tx(kk)*(tx(mm)*(ty(ll)*ty(nn)+tz(ll)*tz(nn))+&
             & tx(nn)*(ty(ll)*ty(mm)+tz(ll)*tz(mm)))+&
             & tx(ll)*(tx(mm)*(ty(kk)*ty(nn)+tz(kk)*tz(nn))+&
             & tx(nn)*(ty(kk)*ty(mm)+tz(kk)*tz(mm))))
                  !      TO AVOID COMPILER DIFFICULTIES THIS IS DIVIDED
              temp1 = tdx(kk) * (tx(mm)*(ty(ll)*ty(nn)+tz(ll)*tz(nn))+&
             & tx(nn)*(ty(ll)*ty(mm)+tz(ll)*tz(mm))) + tdx(ll) * &
             & (tx(mm)*(ty(kk)*ty(nn)+tz(kk)*tz(nn))+tx(nn)*(ty(kk)*ty(mm)+&
             & tz(kk)*tz(mm))) + tx(kk) * (tdx(mm)*(ty(ll)*ty(nn)+&
             & tz(ll)*tz(nn))+tdx(nn)*(ty(ll)*ty(mm)+tz(ll)*tz(mm))) + tx(ll) &
             & * (tdx(mm)*(ty(kk)*ty(nn)+tz(kk)*tz(nn))+tdx(nn)*(ty(kk)*ty(mm)+&
             & tz(kk)*tz(mm)))
              temp2 = tx(kk) * (tx(mm)*(tdy(ll)*ty(nn)+ty(ll)*tdy(nn)+&
             & tdz(ll)*tz(nn)+tz(ll)*tdz(nn))+tx(nn)*(tdy(ll)*ty(mm)+&
             & ty(ll)*tdy(mm)+tdz(ll)*tz(mm)+tz(ll)*tdz(mm))) + tx(ll) * &
             & (tx(mm)*(tdy(kk)*ty(nn)+ty(kk)*tdy(nn)+tdz(kk)*tz(nn)+&
             & tz(kk)*tdz(nn))+tx(nn)*(tdy(kk)*ty(mm)+ty(kk)*tdy(mm)+&
             & tdz(kk)*tz(mm)+tz(kk)*tdz(mm)))
              dr(isp) = dr(isp) + g(20) * (temp1+temp2)
              dr(isp) = dr(isp) + dg(21) * (ty(kk)*ty(ll)*tz(mm)*tz(nn)+&
             & tz(kk)*tz(ll)*ty(mm)*ty(nn)) + g(21) * &
             & (tdy(kk)*ty(ll)*tz(mm)*tz(nn)+ty(kk)*tdy(ll)*tz(mm)*tz(nn)+&
             & ty(kk)*ty(ll)*tdz(mm)*tz(nn)+ty(kk)*ty(ll)*tz(mm)*tdz(nn)+&
             & tdz(kk)*tz(ll)*ty(mm)*ty(nn)+tz(kk)*tdz(ll)*ty(mm)*ty(nn)+&
             & tz(kk)*tz(ll)*tdy(mm)*ty(nn)+tz(kk)*tz(ll)*ty(mm)*tdy(nn))
              dr(isp) = dr(isp) + dg(22) * (ty(kk)*tz(ll)+tz(kk)*ty(ll)) * &
             & (ty(mm)*tz(nn)+tz(mm)*ty(nn)) + g(22) * ((tdy(kk)*tz(ll)+&
             & ty(kk)*tdz(ll)+tdz(kk)*ty(ll)+tz(kk)*tdy(ll))*(ty(mm)*tz(nn)+&
             & tz(mm)*ty(nn))+ (ty(kk)*tz(ll)+tz(kk)*ty(ll))*(tdy(mm)*tz(nn)+&
             & ty(mm)*tdz(nn)+tdz(mm)*ty(nn)+tz(mm)*tdy(nn)))
            end if
          end do
        end do
      end do
    end do
      return  
      end subroutine delmol 


      subroutine delri(dg, ni, nj, rr, del1) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      use funcon_C, only : a0, ev
      use molkst_C, only : numcal
      use parameters_C, only : natorb, am, ad, aq, dd, qq
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ni 
      integer , intent(in) :: nj 
      double precision , intent(in) :: rr 
      double precision , intent(in) :: del1 
      double precision , intent(out) :: dg(22) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icalcn 
      double precision :: term, da, db, qa, qb, aee, ee, ade, aqe, dze, &
        qzze, qxxe, aed, aeq, edz, eqzz, eqxx, add, adq, aqd, aqq, dxdx, dzdz, &
        dzqxx, qxxdz, dzqzz, qzzdz, qxxqxx, qxxqyy, qxxqzz, qzzqxx, qzzqzz, &
        dxqxz, qxzdx, qxzqxz 

      save icalcn 
!-----------------------------------------------
!***********************************************************************
!                                                                      *
!    ON INPUT NI = ATOMIC NUMBER OF FIRST ATOM                         *
!             NJ = ATOMIC NUMBER OF SECOND ATOM                        *
!             RR = INTERATOMIC DISTANCE IN BOHRS                       *
!                                                                      *
!***********************************************************************
      data icalcn/ 0/  
      dze = 0.d0
      qxxe = 0.d0
      qzze = 0.d0
      if (icalcn /= numcal) then 
        icalcn = numcal  
      end if 
      term = (ev*del1)/(rr*a0*a0) 
      da = dd(ni) 
      db = dd(nj) 
      qa = qq(ni) 
      qb = qq(nj) 
!   HYDROGEN-HYDROGEN
      aee = 0.25D0*(1.0D0/am(ni)+1.0D0/am(nj))**2 
      ee = -rr/sqrt(rr**2 + aee)**3 
      dg(1) = term*ee 
      if (natorb(ni)<=2 .and. natorb(nj)<=2) return  
      if (natorb(ni) > 2) then 
!   HEAVY ATOM-HYDROGEN
        ade = 0.25D0*(1.0D0/ad(ni)+1.0D0/am(nj))**2 
        aqe = 0.25D0*(1.0D0/aq(ni)+1.0D0/am(nj))**2 
        dze = (rr + da)/sqrt((rr + da)**2 + ade)**3 - (rr - da)/sqrt((rr - da)&
          **2 + ade)**3 
        qzze = (-(rr + 2.0D0*qa)/sqrt((rr + 2.0D0*qa)**2 + aqe)**3) - (rr - &
          2.0D0*qa)/sqrt((rr - 2.0D0*qa)**2 + aqe)**3 + (2.0D0*rr)/sqrt(rr**2&
           + aqe)**3 
        qxxe = (-(2.0D0*rr)/sqrt(rr**2 + 4.0D0*qa**2 + aqe)**3) + (2.0D0*rr)/&
          sqrt(rr**2 + aqe)**3 
        dg(2) = -(term*dze)/2.0D0 
        dg(3) = term*(ee + qzze/4.0D0) 
        dg(4) = term*(ee + qxxe/4.0D0) 
        if (natorb(nj) <= 2) return  
      end if 
!   HYDROGEN-HEAVY ATOM
      aed = 0.25D0*(1.0D0/am(ni)+1.0D0/ad(nj))**2 
      aeq = 0.25D0*(1.0D0/am(ni)+1.0D0/aq(nj))**2 
      edz = (rr - db)/sqrt((rr - db)**2 + aed)**3 - (rr + db)/sqrt((rr + db)**2&
         + aed)**3 
      eqzz = (-(rr - 2.0D0*qb)/sqrt((rr - 2.0D0*qb)**2 + aeq)**3) - (rr + 2.0D0&
        *qb)/sqrt((rr + 2.0D0*qb)**2 + aeq)**3 + (2.0D0*rr)/sqrt(rr**2 + aeq)**&
        3 
      eqxx = (-(2.0D0*rr)/sqrt(rr**2 + 4.0D0*qb**2 + aeq)**3) + (2.0D0*rr)/&
        sqrt(rr**2 + aeq)**3 
      dg(5) = -(term*edz)/2.0D0 
      dg(11) = term*(ee + eqzz/4.0D0) 
      dg(12) = term*(ee + eqxx/4.0D0) 
      if (natorb(ni) <= 2) return  
!   HEAVY ATOM-HEAVY ATOM
      add = 0.25D0*(1.D0/ad(ni)+1.D0/ad(nj))**2 
      adq = 0.25D0*(1.D0/ad(ni)+1.D0/aq(nj))**2 
      aqd = 0.25D0*(1.D0/aq(ni)+1.D0/ad(nj))**2 
      aqq = 0.25D0*(1.D0/aq(ni)+1.D0/aq(nj))**2 
      dxdx = (-(2.D0*rr)/sqrt(rr**2 + (da - db)**2 + add)**3) + (2.D0*rr)/sqrt(&
        rr**2 + (da + db)**2 + add)**3 
      dzdz = (-(rr + da - db)/sqrt((rr + da - db)**2 + add)**3) - (rr - da + db&
        )/sqrt((rr - da + db)**2 + add)**3 + (rr - da - db)/sqrt((rr - da - db)&
        **2 + add)**3 + (rr + da + db)/sqrt((rr + da + db)**2 + add)**3 
      dzqxx = 2.D0*(rr + da)/sqrt((rr + da)**2 + 4.D0*qb**2 + adq)**3 - 2.D0*(&
        rr - da)/sqrt((rr - da)**2 + 4.D0*qb**2 + adq)**3 - 2.D0*(rr + da)/&
        sqrt((rr + da)**2 + adq)**3 + 2.D0*(rr - da)/sqrt((rr - da)**2 + adq)**&
        3 
      qxxdz = 2.D0*(rr - db)/sqrt((rr - db)**2 + 4.D0*qa**2 + aqd)**3 - 2.D0*(&
        rr + db)/sqrt((rr + db)**2 + 4.D0*qa**2 + aqd)**3 - 2.D0*(rr - db)/&
        sqrt((rr - db)**2 + aqd)**3 + 2.D0*(rr + db)/sqrt((rr + db)**2 + aqd)**&
        3 
      dzqzz = (rr + da - 2.D0*qb)/sqrt((rr + da - 2.D0*qb)**2 + adq)**3 - (rr&
         - da - 2.D0*qb)/sqrt((rr - da - 2.D0*qb)**2 + adq)**3 + (rr + da + &
        2.D0*qb)/sqrt((rr + da + 2.D0*qb)**2 + adq)**3 - (rr - da + 2.D0*qb)/&
        sqrt((rr - da + 2.D0*qb)**2 + adq)**3 + 2.D0*(rr - da)/sqrt((rr - da)**&
        2 + adq)**3 - 2.D0*(rr + da)/sqrt((rr + da)**2 + adq)**3 
      qzzdz = (rr + 2.D0*qa - db)/sqrt((rr + 2.D0*qa - db)**2 + aqd)**3 - (rr&
         + 2.D0*qa + db)/sqrt((rr + 2.D0*qa + db)**2 + aqd)**3 + (rr - 2.D0*qa&
         - db)/sqrt((rr - 2.D0*qa - db)**2 + aqd)**3 - (rr - 2.D0*qa + db)/&
        sqrt((rr - 2.D0*qa + db)**2 + aqd)**3 - 2.D0*(rr - db)/sqrt((rr - db)**&
        2 + aqd)**3 + 2.D0*(rr + db)/sqrt((rr + db)**2 + aqd)**3 
      qxxqxx = (-(2.D0*rr)/sqrt(rr**2 + 4.D0*(qa - qb)**2 + aqq)**3) - (2.D0*rr&
        )/sqrt(rr**2 + 4.D0*(qa + qb)**2 + aqq)**3 + (4.D0*rr)/sqrt(rr**2 + &
        4.D0*qa**2 + aqq)**3 + (4.D0*rr)/sqrt(rr**2 + 4.D0*qb**2 + aqq)**3 - (&
        4.D0*rr)/sqrt(rr**2 + aqq)**3 
      qxxqyy = (-(4.D0*rr)/sqrt(rr**2 + 4.D0*qa**2 + 4.D0*qb**2 + aqq)**3) + (&
        4.D0*rr)/sqrt(rr**2 + 4.D0*qa**2 + aqq)**3 + (4.D0*rr)/sqrt(rr**2 + &
        4.D0*qb**2 + aqq)**3 - (4.D0*rr)/sqrt(rr**2 + aqq)**3 
      qxxqzz = (-2.D0*(rr - 2.D0*qb)/sqrt((rr - 2.D0*qb)**2 + 4.D0*qa**2 + aqq)&
        **3) - 2.D0*(rr + 2.D0*qb)/sqrt((rr + 2.D0*qb)**2 + 4.D0*qa**2 + aqq)**&
        3 + 2.D0*(rr - 2.D0*qb)/sqrt((rr - 2.D0*qb)**2 + aqq)**3 + 2.D0*(rr + &
        2.D0*qb)/sqrt((rr + 2.D0*qb)**2 + aqq)**3 + (4.D0*rr)/sqrt(rr**2 + 4.D0&
        *qa**2 + aqq)**3 - (4.D0*rr)/sqrt(rr**2 + aqq)**3 
      qzzqxx = (-2.D0*(rr + 2.D0*qa)/sqrt((rr + 2.D0*qa)**2 + 4.D0*qb**2 + aqq)&
        **3) - 2.D0*(rr - 2.D0*qa)/sqrt((rr - 2.D0*qa)**2 + 4.D0*qb**2 + aqq)**&
        3 + 2.D0*(rr + 2.D0*qa)/sqrt((rr + 2.D0*qa)**2 + aqq)**3 + 2.D0*(rr - &
        2.D0*qa)/sqrt((rr - 2.D0*qa)**2 + aqq)**3 + (4.D0*rr)/sqrt(rr**2 + 4.D0&
        *qb**2 + aqq)**3 - (4.D0*rr)/sqrt(rr**2 + aqq)**3 
      qzzqzz = (-(rr + 2.D0*qa - 2.D0*qb)/sqrt((rr + 2.D0*qa - 2.D0*qb)**2 + &
        aqq)**3) - (rr + 2.D0*qa + 2.D0*qb)/sqrt((rr + 2.D0*qa + 2.D0*qb)**2 + &
        aqq)**3 - (rr - 2.D0*qa - 2.D0*qb)/sqrt((rr - 2.D0*qa - 2.D0*qb)**2 + &
        aqq)**3 - (rr - 2.D0*qa + 2.D0*qb)/sqrt((rr - 2.D0*qa + 2.D0*qb)**2 + &
        aqq)**3 + 2.D0*(rr - 2.D0*qa)/sqrt((rr - 2.D0*qa)**2 + aqq)**3 + 2.D0*(&
        rr + 2.D0*qa)/sqrt((rr + 2.D0*qa)**2 + aqq)**3 + 2.D0*(rr - 2.D0*qb)/&
        sqrt((rr - 2.D0*qb)**2 + aqq)**3 + 2.D0*(rr + 2.D0*qb)/sqrt((rr + 2.D0*&
        qb)**2 + aqq)**3 - (4.D0*rr)/sqrt(rr**2 + aqq)**3 
      dxqxz = 2.D0*(rr - qb)/sqrt((rr - qb)**2 + (da - qb)**2 + adq)**3 - 2.D0*&
        (rr + qb)/sqrt((rr + qb)**2 + (da - qb)**2 + adq)**3 - 2.D0*(rr - qb)/&
        sqrt((rr - qb)**2 + (da + qb)**2 + adq)**3 + 2.D0*(rr + qb)/sqrt((rr + &
        qb)**2 + (da + qb)**2 + adq)**3 
      qxzdx = 2.D0*(rr + qa)/sqrt((rr + qa)**2 + (qa - db)**2 + aqd)**3 - 2.D0*&
        (rr - qa)/sqrt((rr - qa)**2 + (qa - db)**2 + aqd)**3 - 2.D0*(rr + qa)/&
        sqrt((rr + qa)**2 + (qa + db)**2 + aqd)**3 + 2.D0*(rr - qa)/sqrt((rr - &
        qa)**2 + (qa + db)**2 + aqd)**3 
      qxzqxz = (-2.D0*(rr + qa - qb)/sqrt((rr + qa - qb)**2 + (qa - qb)**2 + &
        aqq)**3) + 2.D0*(rr + qa + qb)/sqrt((rr + qa + qb)**2 + (qa - qb)**2 + &
        aqq)**3 + 2.D0*(rr - qa - qb)/sqrt((rr - qa - qb)**2 + (qa - qb)**2 + &
        aqq)**3 - 2.D0*(rr - qa + qb)/sqrt((rr - qa + qb)**2 + (qa - qb)**2 + &
        aqq)**3 + 2.D0*(rr + qa - qb)/sqrt((rr + qa - qb)**2 + (qa + qb)**2 + &
        aqq)**3 - 2.D0*(rr + qa + qb)/sqrt((rr + qa + qb)**2 + (qa + qb)**2 + &
        aqq)**3 - 2.D0*(rr - qa - qb)/sqrt((rr - qa - qb)**2 + (qa + qb)**2 + &
        aqq)**3 + 2.D0*(rr - qa + qb)/sqrt((rr - qa + qb)**2 + (qa + qb)**2 + &
        aqq)**3 
      dg(6) = (term*dzdz)/4.0D0 
      dg(7) = (term*dxdx)/4.0D0 
      dg(8) = -term*(edz/2.0D0 + qzzdz/8.0D0) 
      dg(9) = -term*(edz/2.0D0 + qxxdz/8.0D0) 
      dg(10) = -(term*qxzdx)/8.0D0 
      dg(13) = -term*(dze/2.0D0 + dzqzz/8.0D0) 
      dg(14) = -term*(dze/2.0D0 + dzqxx/8.0D0) 
      dg(15) = -(term*dxqxz)/8.0D0 
      dg(16) = term*(ee + eqzz/4.0D0 + qzze/4.0D0 + qzzqzz/16.0D0) 
      dg(17) = term*(ee + eqzz/4.0D0 + qxxe/4.0D0 + qxxqzz/16.0D0) 
      dg(18) = term*(ee + eqxx/4.0D0 + qzze/4.0D0 + qzzqxx/16.0D0) 
      dg(19) = term*(ee + eqxx/4.0D0 + qxxe/4.0D0 + qxxqxx/16.0D0) 
      dg(20) = (term*qxzqxz)/16.0D0 
      dg(21) = term*(ee + eqxx/4.0D0 + qxxe/4.0D0 + qxxqyy/16.0D0) 
      dg(22) = term*(qxxqxx - qxxqyy)/32.0D0 
      return  
      end subroutine delri 


      subroutine rotat(coord, i, j, ix, rij, del1, idx) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE analyt_C, only : tx, ty, tz,  tdx, tdy, tdz 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: i 
      integer , intent(in) :: j 
      integer , intent(in) :: ix 
      integer , intent(in) :: idx 
      double precision , intent(in) :: rij 
      double precision , intent(in) :: del1 
      double precision , intent(in) :: coord(3,25) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      double precision :: xd, yd, zd, rxy, ryz, rzx, term, tmp 
!-----------------------------------------------
      xd = coord(1,i) - coord(1,j) 
      yd = coord(2,i) - coord(2,j) 
      zd = coord(3,i) - coord(3,j) 
      rxy = sqrt(xd*xd + yd*yd) 
      ryz = sqrt(yd*yd + zd*zd) 
      rzx = sqrt(zd*zd + xd*xd) 
      tx = 0.0D0 
      ty = 0.0D0 
      tz = 0.0D0 
      tdx = 0.0D0 
      tdy = 0.0D0 
      tdz = 0.0D0 
      if (rxy < 1.0D-4) then 
!   MOLECULAR Z AXIS IS PARALLEL TO DIATOMIC Z AXIS
        tx(3) = 1.0D0 
        if (zd < 0.0D0) tx(3) = -1.0D0 
        ty(2) = 1.0D0 
        tz(1) = tx(3) 
        if (idx == 1) return  
        if (ix == 1) tdx(1) = 1.0D0/rij 
        if (ix == 2) tdx(2) = 1.0D0/rij 
        if (ix == 1) tdz(3) = -1.0D0/rij 
        if (ix == 2) tdy(3) = -tx(3)/rij 
      else if (ryz < 1.0D-4) then 
!   MOLECULAR X AXIS IS PARALLEL TO DIATOMIC Z AXIS
        tx(1) = 1.0D0 
        if (xd < 0.0D0) tx(1) = -1.0D0 
        ty(2) = tx(1) 
        tz(3) = 1.0D0 
        if (idx == 1) return  
        if (ix == 2) tdx(2) = 1.0D0/rij 
        if (ix == 3) tdx(3) = 1.0D0/rij 
        if (ix == 2) tdy(1) = -1.0D0/rij 
        if (ix == 3) tdz(1) = -tx(1)/rij 
      else if (rzx < 1.0D-4) then 
!   MOLECULAR Y AXIS IS PARALLEL TO DIATOMIC Z AXIS
        tx(2) = 1.0D0 
        if (yd < 0.0D0) tx(2) = -1.0D0 
        ty(1) = -tx(2) 
        tz(3) = 1.0D0 
        if (idx == 1) return  
        if (ix == 1) tdx(1) = 1.0D0/rij 
        if (ix == 3) tdx(3) = 1.0D0/rij 
        if (ix == 1) tdy(2) = 1.0D0/rij 
        if (ix == 3) tdz(2) = -tx(2)/rij 
      else 
        tx(1) = xd/rij 
        tx(2) = yd/rij 
        tx(3) = zd/rij 
        tz(3) = rxy/rij 
        ty(1) = -tx(2)*sign(1.0D0,tx(1))/tz(3) 
        ty(2) = abs(tx(1)/tz(3)) 
        ty(3) = 0.0D0 
        tz(1) = -tx(1)*tx(3)/tz(3) 
        tz(2) = -tx(2)*tx(3)/tz(3) 
        if (idx == 1) return  
        term = del1/(rij*rij) 
        select case (ix)  
        case (1)  
          tdx(1) = 1.0D0/rij - tx(1)*term 
          tdx(2) = -tx(2)*term 
          tdx(3) = -tx(3)*term 
          tdz(3) = tx(1)/rxy - tz(3)*term 
        case (2)  
          tdx(1) = -tx(1)*term 
          tdx(2) = 1.0D0/rij - tx(2)*term 
          tdx(3) = -tx(3)*term 
          tdz(3) = tx(2)/rxy - tz(3)*term 
        case (3)  
          tdx(1) = -tx(1)*term 
          tdx(2) = -tx(2)*term 
          tdx(3) = 1.0D0/rij - tx(3)*term 
          tdz(3) = -tz(3)*term 
        end select 
        tmp = tdz(3)/tz(3)**2 
        tdy(1) = (-tdx(2)/tz(3)) + tx(2)*tmp 
        if (tx(1) < 0.0D0) tdy(1) = -tdy(1) 
        tdy(2) = tdx(1)/tz(3) - tx(1)*tmp 
        if (tx(1) < 0.0D0) tdy(2) = -tdy(2) 
        tdy(3) = 0.0D0 
        tdz(1) = tx(1)*tx(3)*tmp - (tx(3)*tdx(1)+tx(1)*tdx(3))/tz(3) 
        tdz(2) = tx(2)*tx(3)*tmp - (tx(3)*tdx(2)+tx(2)*tdx(3))/tz(3) 
      end if 
      return  
      end subroutine rotat 


      subroutine ders(m, n, rr, del1, del2, del3, is, iol) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE overlaps_C, only :  ccc,  zzz
      use funcon_C, only : a0
      USE analyt_C, only : ds 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: m 
      integer , intent(in) :: n 
      integer , intent(in) :: is 
      integer , intent(in) :: iol 
      double precision , intent(in) :: rr 
      double precision , intent(in) :: del1 
      double precision , intent(in) :: del2 
      double precision , intent(in) :: del3 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j 
      double precision, dimension(6,6) :: ss 
      double precision :: apb, amb, adb, adr, abn 
!-----------------------------------------------
!***********************************************************************
!                                                                      *
!    ON INPUT M    = INDEX OF FIRST ATOMIC ORBITAL                     *
!             N    = INDEX OF SECOND ATOMIC ORBITAL                    *
!             RR   = SQUARE IF INTERATOMIC DIATANCE (IN BOHR)          *
!             DEL1 = CATERSIAN DISTANCE IN DERIVATIVE DIRECTION        *
!             DEL2 = CARTESIAN DISTANCE IN M A.O.'S DIRECTION          *
!             DEL3 = CARTESIAN DISTANCE IN N A.O.'S DIRECTION          *
!             IS   = INDICATES TYPE OF A.O.-A.O. INTERACTION           *
!                  = 1 S/S, 2 S/P', 3 S/P, 4 P'/S, 5 P/S, 6 P/P',      *
!                    7 P'/P", 8 P'P', 9 P/P                            *
!             IOL  = INDEX FOR STORING DERIVATIVES IN DS               *
!                                                                      *
!***********************************************************************
      do i = 1, 6 
        do j = 1, 6 
          ss(i,j) = 0.0D0 
          apb = zzz(m,i)*zzz(n,j) 
          amb = zzz(m,i) + zzz(n,j) 
          adb = apb/amb 
          adr = min(adb*rr,35.D0) 
          select case (is)  
          case default 
            abn = -2.0D0*adb*del1/a0**2 
          case (2)  
            abn = -4.0D0*adb**2*del1*del2/(sqrt(zzz(n,j))*a0**3) 
          case (3)  
            abn = (2.0D0*adb/(sqrt(zzz(n,j))*a0))*(1.0D0 - 2.0D0*adb*del1**2/a0&
              **2) 
          case (4)  
            abn = 4.0D0*adb**2*del1*del2/(sqrt(zzz(m,i))*a0**3) 
          case (5)  
            abn = -(2.0D0*adb/(sqrt(zzz(m,i))*a0))*(1.0D0 - 2.0D0*adb*del1**2/a0&
              **2) 
          case (6)  
            abn = -(4.0D0*adb**2*del2/(sqrt(apb)*a0**2))*(1.0D0 - 2.0D0*adb*&
              del1**2/a0**2) 
          case (7)  
            abn = 8.0D0*adb**3*del1*del2*del3/(sqrt(apb)*a0**4) 
          case (8)  
            abn = -(8.0D0*adb**2*del1/(sqrt(apb)*a0**2))*(0.5D0 - adb*del2**2/&
              a0**2) 
          case (9)  
            abn = -(8.0D0*adb**2*del1/(sqrt(apb)*a0**2))*(1.5D0 - adb*del1**2/&
              a0**2) 
          end select  
          ss(i,j) = sqrt((2.0D0*sqrt(apb)/amb)**3)*exp((-adr))*abn 
        end do 
      end do 
      do i = 1, 6 
        ds(iol) = ds(iol) + sum(ss(i,:)*ccc(m,i)*ccc(n,:)) 
      end do 
      return  
      end subroutine ders 
