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

      subroutine enpart()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use common_arrays_C, only : h, p, pa, pb, nfirst, nlast, nat, &
      & coord, w, set_a, set_b
      use molkst_C, only : numat, uhf, keywrd, E_disp, E_hb, numcal, itemp_1
      use elemts_C, only : elemnt
      use chanel_C, only : iw
      USE funcon_C, only : fpc
      USE parameters_C, only : uss, upp, udd
!----------------------------------------------------------*
!
!     SUBROUTINE ENPART,  MODIFIED BY TSUNEO HIRANO 1986/6/3/
!
!----------------------------------------------------------*
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k,ii, ia, ib, ni, j, kl, n, iminus, jj, ja, &
      jb, ka, ik, iss, nj, kinc, jap1, kc, kk, jss, l, iap1, ij, jij, jkl, &
      im,km,kik,lm, kil, kjk, i3, jm, kjl, kb, k3, jk, l3, il, jl, j3, numat1, j1, j4, &
      parts(10), nparts, iu, ju
      double precision, dimension(45) :: e1b, e2a
      double precision, dimension(:), allocatable :: w2
      double precision, allocatable, dimension(:,:) :: ea, e, ex
      double precision :: t, sum1, aa, sum, eau, eae, tone, oneii, onejj, g, &
        bb, pij, eabr, eabx, eabee, eaben, eabnn, eabrx, &
        eabe, ttwo, et
      logical :: si, sj, large
      double precision, external :: reada
!-----------------------------------------------
!--- DEFINED HERE, AND TO BE USED FOR ENPART-PRINT ONLY ---*
!--- END OF DIMENSION DEFINITION ----------------- BY TH --*
!***********************************************************************
!
! *** ENERGY PARTITIONING WITHIN THE UMNDO SCHEME
!     ROUTINE WRITTEN BY S.OLIVELLA, BARCELONA NOV. 1979.
!     EXTENDED TO AM1 AND PM3 BY JJPS.
!
!   ON INPUT UHF     = .TRUE. IF A U.H.F. CALCULATION.
!            H       = ONE-ELECTRON MATRIX.
!            pa      = ALPHA ELECTRON DENSITY.
!            pb      = BETA ELECTRON DENSITY.
!            P       = TOTAL ELECTRON DENSITY.
!            Q       = ATOM ELECTRON DENSITIES.
!
!    NOTHING IS CHANGED ON EXIT.
!
!***********************************************************************
! *** RECALCULATE THE DENSITY MATRICES IN THE UHF SCHEME
!
      allocate (ea((numat*(numat + 1))/2,2), e((numat*(numat + 1))/2,4), &
      ex((numat*(numat + 1))/2,3), w2(8100), stat = i)
      if (i /= 0) then
        write(iw,*)" A problem occurred during memory assignment in ENPART. "
        call mopend("A problem occurred during memory assignment in ENPART.")
        return
      end if
      ea = 0.d0
      if (.not.uhf) pb = pa
      large = (index(keywrd, " LARGE") /= 0)
!
! *** ONE-CENTER ENERGIES
!
      k = 0
      do ii = 1, numat
        ia = nfirst(ii)
        ib = nlast(ii)
        ni = nat(ii)
        ea(ii,1) = 0.0D0
        go to (60,40,40,40,20,20,20,20,20) ib - ia + 1
   20   continue
        t = udd(ni)
        do j = ia + 4, ib
          ea(ii,1) = ea(ii,1) + p((j*(j+1))/2)*t
        end do
   40   continue
        t = upp(ni)
        do j = ia + 1, ia + 3
          ea(ii,1) = ea(ii,1) + p((j*(j+1))/2)*t
        end do
   60   continue
        ea(ii,1) = ea(ii,1) + p((ia*(ia+1))/2)*uss(ni)
      end do
      if (large) write (iw, '(3/,10X,''TOTAL ENERGY PARTITIONING'')')
      if (large) write (iw, '(/10X,''ALL ENERGIES ARE IN ELECTRON VOLTS'')')
      kl = 0
      k = kl + 1
      kl = kl + 10
      kl = min0(kl,numat)
      do while(numat > kl)
        k = kl + 1
        kl = kl + 10
        kl = min0(kl,numat)
      end do
      eau = 0.0D0
      eae = 0.0D0
      do i = 1, numat
        eau = eau + ea(i,1)
        eae = eae + ea(i,2)
      end do
      tone = eau + eae
! *** TWO-CENTER ENERGIES
!     RESONANCE (E(N,1)) TERMS
      n = 1
      do ii = 2, numat
        e(n,1) = 0.0D0
        ia = nfirst(ii)
        ib = nlast(ii)
        iminus = ii - 1
        oneii = 1.D0
        if (nat(ii) == 102) oneii = 0.D0
        do jj = 1, iminus
          n = n + 1
          ja = nfirst(jj)
          jb = nlast(jj)
          onejj = 1.D0
          if (nat(jj) == 102) onejj = 0.D0
          e(n,1) = 0.0D0
          do i = ia, ib
            ka = (i*(i - 1))/2
            do k = ja, jb
              ik = ka + k
              e(n,1) = e(n,1) + 2.0D0*p(ik)*h(ik)*oneii*onejj
            end do
          end do
        end do
        n = n + 1
      end do
!
!
!     CORE-CORE REPULSION (E(N,2)) AND CORE-ELEC. ATTRACTION (E(N,3)).
      n = 0
      kk = 0
      outer_ii_loop: do ii = 1, numat
        ia = nfirst(ii)
        ib = nlast(ii)
        si = (ib < ia)
        ni = nat(ii)
        iss = (ia*(ia+1)) / 2
        iminus = ii - 1
        outer_jj_loop: do jj = 1, iminus
          n = n + 1
          e(n, 2) = 0.0d0
          e(n, 3) = 0.0d0
          ja = nfirst(jj)
          jb = nlast(jj)
          sj = (jb < ja)
          nj = nat(jj)
          jss = (ja*(ja+1)) / 2
!
!  Use rotate here because g, e1b, and e2a need to be re-calculated.
!
          call rotate (ni, nj, coord(1,ii), coord(1,jj), w2, i, e1b, e2a, g)
          if(w2(1) < -1.d20) return ! Dummy use of w2
          if (ib >= ia .and. jb >= ja) then
            kk = kk + 1
!
!  Calculate the nuclear-nuclear term (e(n,2).
!
            e(n, 2) = g
!
!  Calculate the electron-nuclear stabilization. (e(n,3))
!
            if (ib >= ia ) e(n, 3) =  e(n, 3) + (pa(iss) + pb(iss)) * e1b(1)
            if (jb >= ja ) e(n, 3) =  e(n, 3) + (pa(jss) + pb(jss)) * e2a(1)
          else  ! add in sparkle contribution
            e(n, 2) = g
            if (ib >= ia ) e(n, 3) =  (pa(iss) + pb(iss)) * e1b(1)
            if (jb >= ja ) e(n, 3) =  (pa(jss) + pb(jss)) * e2a(1)
          end if
!
!  Calculate the electron-nuclear stabilization. (e(n,3))
!
!
!  The s-s term has already been calculated.  Now do the
!  <ZZ|other terms>  Go along the side of the square matrix of
!  W integrals (in steps of 1).
!  Note:  The first term is the <ZZ|sz> integral.
!
          kinc = ((jb-ja+1)*(jb-ja+2)) / 2 - 1
          jap1 = ja + 1
          i = 1
          do k = jap1, jb
            kc = (k*(k-1)) / 2
            do l = ja, k
              kl = kc + l
              bb = 2.0d0
              if (k == l) then
                bb = 1.0d0
              end if
              if ( si ) then
                i = i + 1  ! Sparkle contribution to A.O. of atom j
                e(n, 3) = e(n, 3) + (pa(kl) + pb(kl)) * bb * e2a(i)
              else
                i = i + 1
                kk = kk + 1
                e(n, 3) = e(n, 3) + (pa(kl) + pb(kl)) * bb * e2a(i)
              end if
            end do
          end do
         !
         !  Calculate the electron-nuclear stabilization.
         !  The s-s term has already been calculated.  Now do the
         !  <ZZ|other terms>  Go down the side of the square matrix of
         !  W integrals (in steps of kinc+1).
         !  Note: The first term is the <sz|ZZ> integral.
         !
          iap1 = ia + 1
          l = 1
          do i = iap1, ib
            ka = (i*(i-1)) / 2
            do j = ia, i
              ij = ka + j
              aa = 2.0d0
              if (i == j) then
                aa = 1.0d0
              end if
              if ( sj ) then
                l = l + 1
                e(n, 3) = e(n, 3) +  (pa(ij) + pb(ij)) * aa * e1b(l)
              else
                l = l + 1
                kk = kk + 1
                e(n, 3) = e(n, 3) + (pa(ij) + pb(ij)) * aa * e1b(l)
               !
               !  KINC is the number of terms on a side, less 1
               !
                kk = kk + kinc
              end if
            end do
          end do
        end do outer_jj_loop
        sum1 = 0.d0
       !
       !   One-center coulomb and exchange terms for atom II.
       !
       !  F(i,j)=F(i,j)+sum(k,l)((PA(k,l)+PB(k,l))*<i,j|k,l>
       !                        -(PA(k,l)        )*<i,k|j,l>), k,l on atom II.
       !
        do i = ia, ib
          aa = 2.d0
          do j = ia, i
            if (i == j) then
              aa = 1.d0
            end if
               !
               !   'J' address in P
               !
            jij = (i*(i-1)) / 2 + j
            sum = 0.d0
            do k = ia, ib
              bb = 2.d0
              do l = ia, k
                if (k == l) then
                  bb = 1.d0
                end if
                     !
                     !   'J' address in P
                     !
                jkl = (k*(k-1)) / 2 + l
                     !
                     !   'K' addresses in P
                     !
                im = Max (i, k)
                km = Min (i, k)
                kik = (im*(im-1)) / 2 + km
                im = Max (i, l)
                lm = Min (i, l)
                kil = (im*(im-1)) / 2 + lm
                jm = Max (j, k)
                km = Min (j, k)
                kjk = (jm*(jm-1)) / 2 + km
                jm = Max (j, l)
                lm = Min (j, l)
                kjl = (jm*(jm-1)) / 2 + lm
                     !
                     !   The term itself
                     !
                kk = kk + 1
                sum = sum + aa * bb * w(kk) * ((pa(jij)+&
               & pb(jij))*(pa(jkl)+pb(jkl))-(pa(kik)*pa(kjl)+&
               & pa(kil)*pa(kjk)+pb(kik)*pb(kjl)+&
               & pb(kil)*pb(kjk))*0.5d0)
              end do
            end do
            sum1 = sum1 + sum
          end do
        end do
        ea(ii, 2) = sum1 * 0.5d0
        n = n + 1
        e(n, 2) = 0.0d0
        e(n, 3) = 0.0d0
      end do outer_ii_loop
      eau = 0.0d0
      eae = 0.0d0
      do i = 1, numat
        eau = eau + ea(i, 1)
        eae = eae + ea(i, 2)
      end do
      tone = eau + eae
!
!  Two center Coulomb and exchange terms
!
! (All indices start over)
!
      !     COULOMB (E(N,4)) AND EXCHANGE (EX(N)) TERMS
      n = 1
      ia = nfirst(1)
      ib = nlast(1)
      kk = (((ib-ia+1)*(ib-ia+2))/2) ** 2
      do ii = 2, numat
        e(n, 4) = 0.0d0
        ex(n, 1) = 0.0d0
        ia = nfirst(ii)
        ib = nlast(ii)
        iminus = ii - 1
        do jj = 1, iminus
          ja = nfirst(jj)
          jb = nlast(jj)
          n = n + 1
          e(n, 4) = 0.0d0
          ex(n, 1) = 0.0d0
          i3 = 0
          do i = ia, ib
            i3 = i3 + 1
            ka = (i*(i-1)) / 2
            j3 = 0
            do j = ia, i
              j3 = j3 + 1
              kb = (j*(j-1)) / 2
              ij = ka + j
              aa = 2.0d0
              if (i == j) then
                aa = 1.0d0
              end if
              pij = (pa(ij)+pb(ij))
              k3 = 0
              do k = ja, jb
                k3 = k3 + 1
                kc = (k*(k-1)) / 2
                ik = ka + k
                jk = kb + k
                l3 = 0
                do l = ja, k
                  l3 = l3 + 1
                  il = ka + l
                  jl = kb + l
                  kl = kc + l
                  bb = 2.0d0
                  if (k == l) then
                    bb = 1.0d0
                  end if
                  kk = kk + 1
                  g = w(kk)
                  e(n, 4) = e(n, 4) + aa * bb * g * pij * &
                 & (pa(kl)+pb(kl))
                  ex(n, 1) = ex(n, 1) - 0.5d0 * aa * bb * g * &
                 & (pa(ik)*pa(jl)+pa(il)*pa(jk)+pb(ik)*pb(jl)+&
                 & pb(il)*pb(jk))
                end do
              end do
            end do
          end do
        end do
        kk = kk + (((ib-ia+1)*(ib-ia+2))/2) ** 2
        n = n + 1
      end do
      numat1 = (numat*(numat + 1))/2
      e(numat1,:) = 0.0D0
      ex(numat1,:) = 0.0D0
!@ --------------------------*
!-----PRINT OUT ONE AND TWO CENTER ENERGIES
!
!     E(I,1):     RESONANCE ENERGY
!     E(I,2):     NUCLEAR-NUCLEAR REPULSION ENERGY
!     E(I,3):     ELECTRON-NUCLEAR ATTRACTION ENERGY
!     E(I,4):     ELECTRON-ELECTRON REPULSION ENERGY
!     EX(I,1):    EXCHANGE  ENERGY
!     EX(I,2):    EXCHANGE + RESONANCE ENERGY
      ex(:numat1,2) = e(:numat1,1) + ex(:numat1,1)
!
!   ADD IN MONOCENTRIC EXCHANGE AND COULOMBIC TERM
!
      do i = 1, numat
        ex((i*(i+1))/2,2) = ea(i,2)
      end do
!
!
      do i = 1, numat
        e((i*(i+1))/2,3) = ea(i,1)
      end do
!
!
      ex(:numat1,3) = e(:numat1,4) + e(:numat1,3) + e(:numat1,2)
!     PRINT OUT OF TOTAL COULOMB TERM
!     PRINT OUT OF TWO-CENTER SUM(OFF-DIAGONAL) +
!                  ONE-CENTER SUM(DIAGONAL).
      if (large) then
        write (iw, '(/,8(10X,A,/))') '  ONE-CENTER TERMS ', &
          'E-E:  ELECTRON-ELECTRON REPULSION', &
          'E-N:  ELECTRON-NUCLEAR ATTRACTION'
        write (iw, '(/,''   ATOM      E-E       E-N    (E-E + E-N)'')')
        do i = 1, numat
          j = (i*(i + 1))/2
          write (iw, '(2X,A2,I4,1X,2F10.4,F10.4)') &
          elemnt(nat(i)), i, ex(j,2), e(j,3), ex(j,2) + e(j,3)
        end do
        write (iw, '(/,8(10X,A,/))') '    TWO-CENTER TERMS', ' ', &
          'J:   RESONANCE ENERGY          E-E: ELECTRON-ELECTRON REPULSION', &
          'K:   EXCHANGE ENERGY           E-N: ELECTRON-NUCLEAR ATTRACTION', &
          '                               N-N: NUCLEAR-NUCLEAR REPULSION', &
          'C:   COULOMBIC INTERACTION = E-E + E-N + N-N', &
          'EE:  TOTAL OF ELECTRONIC AND NUCLEAR ENERGIES'
        write (iw, &
        '(/,''     ATOM          J        K       E-E       E-N      N-N      C        EE'')')
        write (iw, '(''     PAIR'')')
        ij = 0
        do i = 1, numat
          if (i<6 .or. i==numat) then
            do j = 1, i
              ij = ij + 1
              if (i /= j) then
                write (iw, &
        '(1X,A2,I5,1X,A2,I5,1X,2F9.4,F9.4,F10.4,F9.4,F8.4,F9.4)') elemnt(nat(i)), i, &
         elemnt(nat(j)), j, e(ij,1), ex(ij,1), e(ij,4), e(ij,3), &
         e(ij,2), ex(ij,3), ex(ij,2) + ex(ij,3)
              else
                write (iw, *)
              end if
            end do
          else
            do j = 1, i
              ij = ij + 1
              if (i /= j) then
                write (iw, &
        '(1X,A2,I5,1X,A2,I5,1X,2F9.4,F9.4,F10.4,F9.4,F8.4,F9.4)') elemnt(nat(i)), i, &
        elemnt(nat(j)), j, e(ij,1), ex(ij,1), e(ij,4), e(ij,3), e(ij,2), ex(ij,3), ex(ij,2) + ex(ij,3)
              else
                write (iw, &
        '(/,''   ATOM          J        K       E-E       E-N      N-N      C        EE'')')
                write (iw, '(''   PAIR'')')
              end if
            end do
          end if
        end do
      else
        write(iw,'(/12x,a)')"For more detail of energy partitioning, add keyword 'LARGE'"
      end if
!
!  If the system is to be split into parts, do so now.
!
      if (index(keywrd,"ENPART(") /= 0) then
        i = index(keywrd,"ENPART(")
        j = index(keywrd(i:),")") + i - 1
        keywrd(j:j) = ","
        k = 0
        do nparts = 1,10
           parts(nparts) = nint(reada(keywrd,i)) + k
           k = parts(nparts)
           i = index(keywrd(i + 1:),",") + i
           if (i >= j) exit
        end do
        if (k > numat) then
          write(iw,"(a)")"  Number of atoms in parts exceeds total number of atoms"
          write(iw,"(a,i5)")"  Number of atoms in parts:", k
          write(iw,"(a,i5)")"  Number of atoms:         ", numat
        end if
        keywrd(j:j) = ")"
        if (k <= numat) then
          write(iw,"(/,11x,a,i1,a)")" In the energy partitioning, the system is to be split into ",nparts, " parts"
          write(iw,"(15x,a)")" Part     No. of atoms     Atoms"
          k = 0
          do i = 1, nparts
            write(iw,"(16x,i3,i12,i10,' to',i4)")i, parts(i) - k, k + 1, parts(i)
            k = parts(i)
          end do
          if (allocated(set_a))  deallocate(set_a)
          if (allocated(set_b))  deallocate(set_b)
          allocate (set_a(numat), set_b(numat))   ! two fragments
          write(iw,'(40x,a,10x,a)')"     Contribution from:    H-bonds    +  Dispersion   =   Total"
          il = 1
          do nparts  = 1, nparts
            iu = parts(nparts)
!
!  Evaluate hydrogen bonding and dispersion contributions
!  First set: atoms il - iu, second part: jl - ju
!
              itemp_1 =  -99
              set_a = 0
              set_b = 0
              j = 0
              do i = il, iu
                j = j + 1
                set_a(j) = i
                set_b(j) = i
              end do
              numcal = numcal + 1
              call post_scf_corrections(sum, .false.)
!
!   Calculate energy term due to the part
!
            sum = 0.d0
            do i = il, iu
              k = (i*(i-1))/2 + il - 1
              do j = il, i - 1
                k = k + 1
                sum = sum + ex(k,2) + ex(k,3)
              end do
              k = k + 1
              sum = sum + ex(k,2) + e(k,3)
            end do
            write(iw,"(a,i1,a,f14.4,' eV =',f13.3,' Kcal/mol', f10.3, f16.3,f15.3,a)") &
            " Part ", nparts," self - energy:",sum, sum*fpc(9), E_hb, E_disp, E_disp + E_hb, " Kcal/mol"
!
!  Calculate energy term due to interaction of parts
!
            jl = 1
            do l = 1, nparts - 1
              ju = parts(l)
              sum = 0.d0
              do i = il, iu
                k = (i*(i - 1))/2 + jl - 1
                do j = jl, ju
                  k = k + 1
                  sum = sum + ex(k,2) + ex(k,3)
                end do
              end do
!
!  Evaluate hydrogen bonding and dispersion contributions
!  First set: atoms il - iu, second part: jl - ju
!
              itemp_1 =  -99
              set_a = 0
              set_b = 0
              j = 0
              do i = il, iu
                j = j + 1
                set_a(j) = i
              end do
              j = 0
              do i = jl, ju
                j = j + 1
                set_b(j) = i
              end do
              numcal = numcal + 1
              call post_scf_corrections(aa, .false.)
              write(iw,"(a,i1,a,i1,a,f11.4,' eV =',f13.3,' Kcal/mol', f10.3, f16.3,f15.3,a )") &
              " Parts ", nparts, " - ",l," interaction:",sum, sum*fpc(9), E_hb, E_disp, E_disp + E_hb, " Kcal/mol"
              jl = ju + 1
            end do
            il = iu + 1
          end do
        end if
      end if
      itemp_1 =  0
      numcal = numcal + 1
      call post_scf_corrections(sum, .false.)
      write(iw,'(/,a, f10.3, a,f9.3,a,f10.3,a,/)') &
        "                    Total contribution from hydrogen bonds:    ", &
        E_hb," disp.:",E_disp," Tot:",E_disp + E_hb," kcal/mol"

!
!     ++++   TOTALS   ++++
!
      eabr = 0.0D0
      eabx = 0.0D0
      eabee = 0.0D0
      eaben = 0.0D0
      eabnn = 0.0D0
      e((/(j1,j1=1,numat)/)*(/(j4,j4=2,numat+1)/)/2,3) = 0.D0
      do i = 1, numat1
        eabr = eabr + e(i,1)
        eabx = eabx + ex(i,1)
        eabee = eabee + e(i,4)
        eaben = eaben + e(i,3)
        eabnn = eabnn + e(i,2)
      end do
      eabrx = eabr + eabx
      eabe = eabee + eaben + eabnn
      ttwo = eabrx + eabe
      et = tone + ttwo
!@ ***************************************************************
      write (iw, 450)
  450 format(/,/,/,'***  SUMMARY OF ENERGY PARTITION  ***')
      write (iw, 460)
  460 format(' ','---------------------------------------')
      write (iw, '(''     ONE-CENTER TERMS'')')
      write (iw, 470) eau
  470 format(/,' ELECTRON-NUCLEAR  (ONE-ELECTRON) ',f17.4,' EV')
      write (iw, 480) eae
  480 format(' ELECTRON-ELECTRON (TWO-ELECTRON) ',f17.4,' EV')
      write (iw, 490) tone
  490 format(/,' TOTAL OF ONE-CENTER TERMS ',9x,f15.4,' EV')
      write (iw, 460)
      write (iw, '(''     TWO-CENTER TERMS'')')
      write (iw, 500) eabr
  500 format(/,' RESONANCE ENERGY',19x,f15.4,' EV')
      write (iw, 510) eabx
  510 format(' EXCHANGE ENERGY ',19x,f15.4,' EV')
      write (iw, 520) eabrx
  520 format(/,' EXCHANGE + RESONANCE ENERGY:       ',f15.4,' EV')
      write (iw, 530) eabee
  530 format(/,' ELECTRON-ELECTRON REPULSION',f23.4,' EV')
      write (iw, 540) eaben
  540 format(' ELECTRON-NUCLEAR ATTRACTION',f23.4,' EV')
      write (iw, 550) eabnn
  550 format(' NUCLEAR-NUCLEAR REPULSION  ',f23.4,' EV')
      write (iw, 560) eabe
  560 format(/,' TOTAL ELECTROSTATIC INTERACTION ',f18.4,' EV',/)
      write (iw, 570) ttwo
  570 format(' GRAND TOTAL OF TWO-CENTER TERMS',f19.4,' EV')
      write (iw, 460)
      write (iw, 580) et
  580 format(' ETOT (EONE + ETWO) ',14x,f17.4,' EV'/,/)
      deallocate (ea, e, ex, w2)
      return
      end subroutine enpart
