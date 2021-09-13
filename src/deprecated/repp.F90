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

      subroutine repp(ni, nj, rij, ri, core, cutoff, a0, ev, ev1, ev2, ev3, ev4) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      use parameters_C, only : natorb, am, ad, aq, dd, qq, tore

!***********************************************************************
!
!..VECTOR VERSION WRITTEN BY ERNEST R. DAVIDSON, INDIANA UNIVERSITY
!
!
!  REPP CALCULATES THE TWO-ELECTRON REPULSION INTEGRALS AND THE
!       NUCLEAR ATTRACTION INTEGRALS.
!
!     ON INPUT RIJ     = INTERATOMIC DISTANCE
!              NI      = ATOM NUMBER OF FIRST ATOM
!              NJ      = ATOM NUMBER OF SECOND ATOM
!    (REF)     ADD     = ARRAY OF GAMMA, OR TWO-ELECTRON ONE-CENTER,
!                        INTEGRALS.
!    (REF)     TORE    = ARRAY OF NUCLEAR CHARGES OF THE ELEMENTS
!    (REF)     DD      = ARRAY OF DIPOLE CHARGE SEPARATIONS
!    (REF)     QQ      = ARRAY OF QUADRUPOLE CHARGE SEPARATIONS
!              CUTOFF  = MAXIMUM DISTANCE FOR 'NORMAL' INTERACTION.
!                        BEYOND CUTOFF ALL INTERACTIONS BECOME
!                        POINT-CHARGE, AND ARE BROUGHT IN TO
!                        CUTOFF.
!
!
!    ON OUTPUT RI      = ARRAY OF TWO-ELECTRON REPULSION INTEGRALS
!              CORE    = 4 X 2 ARRAY OF ELECTRON-CORE ATTRACTION
!                        INTEGRALS
!
!
! *** THIS ROUTINE COMPUTES THE TWO-CENTRE REPULSION INTEGRALS AND THE
! *** NUCLEAR ATTRACTION INTEGRALS.
! *** THE TWO-CENTRE REPULSION INTEGRALS (OVER LOCAL COORDINATES) ARE
! *** STORED AS FOLLOWS (WHERE P-SIGMA = O,  AND P-PI = P AND P* )
!     (SS/SS)=1,   (SO/SS)=2,   (OO/SS)=3,   (PP/SS)=4,   (SS/OS)=5,
!     (SO/SO)=6,   (SP/SP)=7,   (OO/SO)=8,   (PP/SO)=9,   (PO/SP)=10,
!     (SS/OO)=11,  (SS/PP)=12,  (SO/OO)=13,  (SO/PP)=14,  (SP/OP)=15,
!     (OO/OO)=16,  (PP/OO)=17,  (OO/PP)=18,  (PP/PP)=19,  (PO/PO)=20,
!     (PP/P*P*)=21,   (P*P/P*P)=22.
! *** THE STORAGE OF THE NUCLEAR ATTRACTION INTEGRALS  CORE(KL/IJ) IS
!     (SS/)=1,   (SO/)=2,   (OO/)=3,   (PP/)=4
!     WHERE IJ=1 IF THE ORBITALS CENTRED ON ATOM I,  =2 IF ON ATOM J.
! *** NI AND NJ ARE THE ATOMIC NUMBERS OF THE TWO ELEMENTS.
!
!***********************************************************************
!***********************************************************************
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ni 
      integer , intent(in) :: nj 
      double precision , intent(in) :: rij 
      double precision , intent(in) :: cutoff 
      double precision , intent(in) :: a0 
      double precision , intent(in) :: ev 
      double precision , intent(in) :: ev1 
      double precision , intent(in) :: ev2 
      double precision , intent(in) :: ev3 
      double precision , intent(in) :: ev4 
      double precision , intent(out) :: ri(22) 
      double precision , intent(out) :: core(4,2) 

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i 
      double precision, dimension(72) :: arg, sqr 
      double precision :: td, pp, r, aee, da, qa, ade, aqe, rsq, xxx, ee, db, qb, &
        aed, aeq, axx, adq, aqd, aqq, yyy, zzz, www, dze, qzze, qxxe, edz, eqzz&
        , eqxx, dxdx, dzdz, dzqxx, qxxdz, dzqzz, qzzdz, qxxqxx, qxxqyy, qxxqzz&
        , qzzqxx, qzzqzz, dxqxz, qxzdx, qxzqxz 
      logical :: si, sj 

      save td, pp 
!-----------------------------------------------
      data td/ 2.D00/  
      data pp/ 0.5D00/  
      ri = 0.D0 
      core(:,1) = 0.D0 
      core(:,2) = 0.D0 
      r = min(cutoff,rij)/a0 
        si = natorb(ni) >= 3 
        sj = natorb(nj) >= 3 
        aee = pp/am(ni) + pp/am(nj) 
        aee = aee*aee 
        if (.not.si .and. .not.sj) then 
!
!     HYDROGEN - HYDROGEN  (SS/SS)
!         
          ri(1) = ev/sqrt(r*r + aee) 
          core(1,1) = tore(nj)*ri(1) 
          core(1,2) = tore(ni)*ri(1) 
!
        else if (si .and. .not.sj) then 
!
!     HEAVY ATOM - HYDROGEN
!
          da = dd(ni) 
          qa = qq(ni)*td 
          ade = pp/ad(ni) + pp/am(nj) 
          ade = ade*ade 
          aqe = pp/aq(ni) + pp/am(nj) 
          aqe = aqe*aqe 
          rsq = r*r 
          arg(1) = rsq + aee 
          xxx = r + da 
          arg(2) = xxx*xxx + ade 
          xxx = r - da 
          arg(3) = xxx*xxx + ade 
          xxx = r + qa 
          arg(4) = xxx*xxx + aqe 
          xxx = r - qa 
          arg(5) = xxx*xxx + aqe 
          arg(6) = rsq + aqe 
          arg(7) = arg(6) + qa*qa 
          do i = 1, 7 
            sqr(i) = sqrt(arg(i)) 
          end do 
          ee = ev/sqr(1) 
          ri(1) = ee 
          ri(2) = ev1/sqr(2) - ev1/sqr(3) 
          ri(3) = ee + ev2/sqr(4) + ev2/sqr(5) - ev1/sqr(6) 
          ri(4) = ee + ev1/sqr(7) - ev1/sqr(6) 
          if (rij > cutoff) then 
!
!   Convert RI to point-charge system
!
            ri = 0.D0 
            ri(1) = ee 
            ri(3) = ee 
            ri(4) = ee 
          end if 
          core(1,1) = tore(nj)*ri(1) 
          core(1,2) = tore(ni)*ri(1) 
          core(2,1) = tore(nj)*ri(2) 
          core(3,1) = tore(nj)*ri(3) 
          core(4,1) = tore(nj)*ri(4) 
!
        else if (.not.si .and. sj) then 
!
!     HYDROGEN - HEAVY ATOM
!
          db = dd(nj) 
          qb = qq(nj)*td 
          aed = pp/am(ni) + pp/ad(nj) 
          aed = aed*aed 
          aeq = pp/am(ni) + pp/aq(nj) 
          aeq = aeq*aeq 
          rsq = r*r 
          arg(1) = rsq + aee 
          xxx = r - db 
          arg(2) = xxx*xxx + aed 
          xxx = r + db 
          arg(3) = xxx*xxx + aed 
          xxx = r - qb 
          arg(4) = xxx*xxx + aeq 
          xxx = r + qb 
          arg(5) = xxx*xxx + aeq 
          arg(6) = rsq + aeq 
          arg(7) = arg(6) + qb*qb 
          do i = 1, 7 
            sqr(i) = sqrt(arg(i)) 
          end do 
          ee = ev/sqr(1) 
          ri(1) = ee 
          ri(5) = ev1/sqr(2) - ev1/sqr(3) 
          ri(11) = ee + ev2/sqr(4) + ev2/sqr(5) - ev1/sqr(6) 
          ri(12) = ee + ev1/sqr(7) - ev1/sqr(6) 
          if (rij > cutoff) then 
!
!   Convert RI to point-charge system
!
            ri = 0.D0 
            ri(1) = ee 
            ri(11) = ee 
            ri(12) = ee 
          end if 
          core(1,1) = tore(nj)*ri(1) 
          core(1,2) = tore(ni)*ri(1) 
          core(2,2) = tore(ni)*ri(5) 
          core(3,2) = tore(ni)*ri(11) 
          core(4,2) = tore(ni)*ri(12) 
!
        else 
!
!     HEAVY ATOM - HEAVY ATOM
!
!     DEFINE CHARGE SEPARATIONS.
          da = dd(ni) 
          db = dd(nj) 
          qa = qq(ni)*td 
          qb = qq(nj)*td 
!
          ade = pp/ad(ni) + pp/am(nj) 
          ade = ade*ade 
          aqe = pp/aq(ni) + pp/am(nj) 
          aqe = aqe*aqe 
          aed = pp/am(ni) + pp/ad(nj) 
          aed = aed*aed 
          aeq = pp/am(ni) + pp/aq(nj) 
          aeq = aeq*aeq 
          axx = pp/ad(ni) + pp/ad(nj) 
          axx = axx*axx 
          adq = pp/ad(ni) + pp/aq(nj) 
          adq = adq*adq 
          aqd = pp/aq(ni) + pp/ad(nj) 
          aqd = aqd*aqd 
          aqq = pp/aq(ni) + pp/aq(nj) 
          aqq = aqq*aqq 
          rsq = r*r 
          arg(1) = rsq + aee 
          xxx = r + da 
          arg(2) = xxx*xxx + ade 
          xxx = r - da 
          arg(3) = xxx*xxx + ade 
          xxx = r - qa 
          arg(4) = xxx*xxx + aqe 
          xxx = r + qa 
          arg(5) = xxx*xxx + aqe 
          arg(6) = rsq + aqe 
          arg(7) = arg(6) + qa*qa 
          xxx = r - db 
          arg(8) = xxx*xxx + aed 
          xxx = r + db 
          arg(9) = xxx*xxx + aed 
          xxx = r - qb 
          arg(10) = xxx*xxx + aeq 
          xxx = r + qb 
          arg(11) = xxx*xxx + aeq 
          arg(12) = rsq + aeq 
          arg(13) = arg(12) + qb*qb 
          xxx = da - db 
          arg(14) = rsq + axx + xxx*xxx 
          xxx = da + db 
          arg(15) = rsq + axx + xxx*xxx 
          xxx = r + da - db 
          arg(16) = xxx*xxx + axx 
          xxx = r - da + db 
          arg(17) = xxx*xxx + axx 
          xxx = r - da - db 
          arg(18) = xxx*xxx + axx 
          xxx = r + da + db 
          arg(19) = xxx*xxx + axx 
          xxx = r + da 
          arg(20) = xxx*xxx + adq 
          arg(21) = arg(20) + qb*qb 
          xxx = r - da 
          arg(22) = xxx*xxx + adq 
          arg(23) = arg(22) + qb*qb 
          xxx = r - db 
          arg(24) = xxx*xxx + aqd 
          arg(25) = arg(24) + qa*qa 
          xxx = r + db 
          arg(26) = xxx*xxx + aqd 
          arg(27) = arg(26) + qa*qa 
          xxx = r + da - qb 
          arg(28) = xxx*xxx + adq 
          xxx = r - da - qb 
          arg(29) = xxx*xxx + adq 
          xxx = r + da + qb 
          arg(30) = xxx*xxx + adq 
          xxx = r - da + qb 
          arg(31) = xxx*xxx + adq 
          xxx = r + qa - db 
          arg(32) = xxx*xxx + aqd 
          xxx = r + qa + db 
          arg(33) = xxx*xxx + aqd 
          xxx = r - qa - db 
          arg(34) = xxx*xxx + aqd 
          xxx = r - qa + db 
          arg(35) = xxx*xxx + aqd 
          arg(36) = rsq + aqq 
          xxx = qa - qb 
          arg(37) = arg(36) + xxx*xxx 
          xxx = qa + qb 
          arg(38) = arg(36) + xxx*xxx 
          arg(39) = arg(36) + qa*qa 
          arg(40) = arg(36) + qb*qb 
          arg(41) = arg(39) + qb*qb 
          xxx = r - qb 
          arg(42) = xxx*xxx + aqq 
          arg(43) = arg(42) + qa*qa 
          xxx = r + qb 
          arg(44) = xxx*xxx + aqq 
          arg(45) = arg(44) + qa*qa 
          xxx = r + qa 
          arg(46) = xxx*xxx + aqq 
          arg(47) = arg(46) + qb*qb 
          xxx = r - qa 
          arg(48) = xxx*xxx + aqq 
          arg(49) = arg(48) + qb*qb 
          xxx = r + qa - qb 
          arg(50) = xxx*xxx + aqq 
          xxx = r + qa + qb 
          arg(51) = xxx*xxx + aqq 
          xxx = r - qa - qb 
          arg(52) = xxx*xxx + aqq 
          xxx = r - qa + qb 
          arg(53) = xxx*xxx + aqq 
          qa = qq(ni) 
          qb = qq(nj) 
          xxx = da - qb 
          xxx = xxx*xxx 
          yyy = r - qb 
          yyy = yyy*yyy 
          zzz = da + qb 
          zzz = zzz*zzz 
          www = r + qb 
          www = www*www 
          arg(54) = xxx + yyy + adq 
          arg(55) = xxx + www + adq 
          arg(56) = zzz + yyy + adq 
          arg(57) = zzz + www + adq 
          xxx = qa - db 
          xxx = xxx*xxx 
          yyy = qa + db 
          yyy = yyy*yyy 
          zzz = r + qa 
          zzz = zzz*zzz 
          www = r - qa 
          www = www*www 
          arg(58) = zzz + xxx + aqd 
          arg(59) = www + xxx + aqd 
          arg(60) = zzz + yyy + aqd 
          arg(61) = www + yyy + aqd 
          xxx = qa - qb 
          xxx = xxx*xxx 
          arg(62) = arg(36) + td*xxx 
          yyy = qa + qb 
          yyy = yyy*yyy 
          arg(63) = arg(36) + td*yyy 
          arg(64) = arg(36) + td*(qa*qa + qb*qb) 
          zzz = r + qa - qb 
          zzz = zzz*zzz 
          arg(65) = zzz + xxx + aqq 
          arg(66) = zzz + yyy + aqq 
          zzz = r + qa + qb 
          zzz = zzz*zzz 
          arg(67) = zzz + xxx + aqq 
          arg(68) = zzz + yyy + aqq 
          zzz = r - qa - qb 
          zzz = zzz*zzz 
          arg(69) = zzz + xxx + aqq 
          arg(70) = zzz + yyy + aqq 
          zzz = r - qa + qb 
          zzz = zzz*zzz 
          arg(71) = zzz + xxx + aqq 
          arg(72) = zzz + yyy + aqq 
          do i = 1, 72 
            sqr(i) = sqrt(arg(i)) 
          end do 
          ee = ev/sqr(1) 
          dze = (-ev1/sqr(2)) + ev1/sqr(3) 
          qzze = ev2/sqr(4) + ev2/sqr(5) - ev1/sqr(6) 
          qxxe = ev1/sqr(7) - ev1/sqr(6) 
          edz = (-ev1/sqr(8)) + ev1/sqr(9) 
          eqzz = ev2/sqr(10) + ev2/sqr(11) - ev1/sqr(12) 
          eqxx = ev1/sqr(13) - ev1/sqr(12) 
          dxdx = ev1/sqr(14) - ev1/sqr(15) 
          dzdz = ev2/sqr(16) + ev2/sqr(17) - ev2/sqr(18) - ev2/sqr(19) 
          dzqxx = ev2/sqr(20) - ev2/sqr(21) - ev2/sqr(22) + ev2/sqr(23) 
          qxxdz = ev2/sqr(24) - ev2/sqr(25) - ev2/sqr(26) + ev2/sqr(27) 
          dzqzz = (-ev3/sqr(28)) + ev3/sqr(29) - ev3/sqr(30) + ev3/sqr(31) - &
            ev2/sqr(22) + ev2/sqr(20) 
          qzzdz = (-ev3/sqr(32)) + ev3/sqr(33) - ev3/sqr(34) + ev3/sqr(35) + &
            ev2/sqr(24) - ev2/sqr(26) 
          qxxqxx = ev3/sqr(37) + ev3/sqr(38) - ev2/sqr(39) - ev2/sqr(40) + ev2/&
            sqr(36) 
          qxxqyy = ev2/sqr(41) - ev2/sqr(39) - ev2/sqr(40) + ev2/sqr(36) 
          qxxqzz = ev3/sqr(43) + ev3/sqr(45) - ev3/sqr(42) - ev3/sqr(44) - ev2/&
            sqr(39) + ev2/sqr(36) 
          qzzqxx = ev3/sqr(47) + ev3/sqr(49) - ev3/sqr(46) - ev3/sqr(48) - ev2/&
            sqr(40) + ev2/sqr(36) 
          qzzqzz = ev4/sqr(50) + ev4/sqr(51) + ev4/sqr(52) + ev4/sqr(53) - ev3/&
            sqr(48) - ev3/sqr(46) - ev3/sqr(42) - ev3/sqr(44) + ev2/sqr(36) 
          dxqxz = (-ev2/sqr(54)) + ev2/sqr(55) + ev2/sqr(56) - ev2/sqr(57) 
          qxzdx = (-ev2/sqr(58)) + ev2/sqr(59) + ev2/sqr(60) - ev2/sqr(61) 
          qxzqxz = ev3/sqr(65) - ev3/sqr(67) - ev3/sqr(69) + ev3/sqr(71) - ev3/&
            sqr(66) + ev3/sqr(68) + ev3/sqr(70) - ev3/sqr(72) 
          ri(1) = ee 
          ri(2) = -dze 
          ri(3) = ee + qzze 
          ri(4) = ee + qxxe 
          ri(5) = -edz 
          ri(6) = dzdz 
          ri(7) = dxdx 
          ri(8) = (-edz) - qzzdz 
          ri(9) = (-edz) - qxxdz 
          ri(10) = -qxzdx 
          ri(11) = ee + eqzz 
          ri(12) = ee + eqxx 
          ri(13) = (-dze) - dzqzz 
          ri(14) = (-dze) - dzqxx 
          ri(15) = -dxqxz 
          ri(16) = ee + eqzz + qzze + qzzqzz 
          ri(17) = ee + eqzz + qxxe + qxxqzz 
          ri(18) = ee + eqxx + qzze + qzzqxx 
          ri(19) = ee + eqxx + qxxe + qxxqxx 
          ri(20) = qxzqxz 
          ri(21) = ee + eqxx + qxxe + qxxqyy 
          ri(22) = pp*(qxxqxx - qxxqyy) 
          if (rij > cutoff) then 
!
!   Convert RI to point-charge system
!
            ri = 0.D0 
            ri(1) = ee 
            ri(3) = ee 
            ri(4) = ee 
            ri(11) = ee 
            ri(12) = ee 
            ri(16) = ee 
            ri(17) = ee 
            ri(18) = ee 
            ri(19) = ee 
            ri(21) = ee 
          end if 
!
!     CALCULATE CORE-ELECTRON ATTRACTIONS.
!
          core(1,1) = tore(nj)*ri(1) 
          core(2,1) = tore(nj)*ri(2) 
          core(3,1) = tore(nj)*ri(3) 
          core(4,1) = tore(nj)*ri(4) 
          core(1,2) = tore(ni)*ri(1) 
          core(2,2) = tore(ni)*ri(5) 
          core(3,2) = tore(ni)*ri(11) 
          core(4,2) = tore(ni)*ri(12) 
!
        end if 
!
        return  
      end subroutine repp 
