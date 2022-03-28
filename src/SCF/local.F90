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

      subroutine local(c, nocc, eig, iprint, txt)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : norbs, numat, keywrd, nbeta
      use common_arrays_C, only : nat, nfirst, nlast, p, pa, pb
      use chanel_C, only : iw
      use symmetry_C, only : namo, jndex
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: nocc, iprint
      double precision  :: c(norbs,norbs)
      double precision  :: eig(norbs)
      character*2 :: txt
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(20) :: iel
      integer :: niter, i, j, iter, k, k1, kl, ku, il, iu, i1, ii
      double precision, dimension(norbs) :: eig1, psi1, psi2, cii, refeig
      double precision :: eps, sum, xijjj, xjiii, xiiii, xjjjj, xijij, xiijj, dij, &
        dii, djj, aij, bij, ca, sa, sum1, x, co
      character :: elemnt(99)*2, num_1*1, num_2*1
      double precision, allocatable  :: cold(:,:)

      save elemnt
!-----------------------------------------------
!**********************************************************************
!
!   LOCALISATION SUBROUTINE
! ON INPUT
!        C = EIGENVECTORS IN AN norbs*norbs MATRIX
!        NOCC = NUMBER OF FILLED LEVELS
!        NORBS = NUMBER OF ORBITALS
!        NUMAT = NUMBER OF ATOMS
!        NLAST   = INTEGER ARRAY OF ATOM ORBITAL COUNTERS
!        NFIRST   = INTEGER ARRAY OF ATOM ORBITAL COUNTERS
!
!       SUBROUTINE MAXIMIZES (PSI)**4
!       REFERENCE_
!       A NEW RAPID METHOD FOR ORBITAL LOCALISATION, P.G. PERKINS AND
!       J.J.P. STEWART, J.C.S. FARADAY (II) 78, 285-296, (1982).
!
!       MODIFIED AND CORRECTED TO AVOID SIGMA-PI ORBITAL MIXING BY
!       JUAN CARLOS PANIAGUA, UNIVERSITY OF BARCELONA, MAY 1983.
!
!**********************************************************************
      data elemnt/ 'H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F', 'NE', 'NA', &
        'MG', 'AL', 'SI', 'P', 'S', 'CL', 'AR', 'K', 'CA', 'SC', 'TI', 'V', &
        'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA', 'GE', 'AS', 'SE', 'BR'&
        , 'KR', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG'&
        , 'CD', 'IN', 'SN', 'SB', 'TE', 'I', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR'&
        , 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', &
        'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB', &
        'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH', 'PA', 'U', 'NP', 'PU', &
        'AM', 'CM', 'BK', 'CF', 'XX'/
!
! Set all symmetry names to "a", the lowest symmetry, because LMOs do not normally have any symmetry.
!
      namo(:norbs) = "a   "
      do i = 1, norbs
        jndex(i) = i
      end do
      niter = 100
      eps = 1.0D-10
      refeig(:norbs) = eig(:norbs)
      allocate (cold(norbs,norbs), stat = i)
       if (i /= 0) then
        call memory_error ("Unable to allocate memory in LOCAL")
        return
      end if
      cold(:norbs,:norbs) = c(:norbs,:norbs)
      iter = 0
   20 continue
      sum = 0.D0
      iter = iter + 1
      do i = 1, nocc
        do j = 1, nocc
          if (j == i) cycle
          xijjj = 0.0D0
          xjiii = 0.0D0
          xiiii = 0.0D0
          xjjjj = 0.0D0
          xijij = 0.0D0
          xiijj = 0.0D0
          k = 1
          if (norbs > 0) then
            psi1(:norbs) = c(:norbs,i)
            psi2(:norbs) = c(:norbs,j)
            k = norbs + 1
          end if
! NOW FOLLOWS THE RATE-DETERMINING STEP FOR THE CALCULATION
          do k1 = 1, numat
            kl = nfirst(k1)
            ku = nlast(k1)
            dij = 0.D0
            dii = 0.D0
            djj = 0.D0
            k = kl
            if (ku - kl + 1 > 0) then
              do k = 1, ku - kl + 1
                dij = dij + psi1(k-1+kl)*psi2(k-1+kl)
                dii = dii + psi1(k-1+kl)*psi1(k-1+kl)
                djj = djj + psi2(k-1+kl)*psi2(k-1+kl)
              end do
              k = ku + 1
            end if
            xijjj = xijjj + dij*djj
            xjiii = xjiii + dij*dii
            xiiii = xiiii + dii*dii
            xjjjj = xjjjj + djj*djj
            xijij = xijij + dij*dij
            xiijj = xiijj + dii*djj
          end do
          aij = xijij - (xiiii + xjjjj - 2.0D0*xiijj)/4.0D0
          bij = xjiii - xijjj
          ca = sqrt(aij*aij + bij*bij)
          sa = aij + ca
          if (sa < 1.0D-14) cycle
          sum = sum + sa
          ca = -aij/ca
          ca = (1.0D0 + sqrt((1.0D0 + ca)/2.0D0))/2.0D0
          if ((2.0D0*ca - 1.0D0)*bij < 0.0D0) ca = 1.0D0 - ca
          sa = sqrt(1.0D0 - ca)
          ca = sqrt(ca)
          k = 1
          if (norbs > 0) then
            c(:norbs,i) = ca*psi1(:norbs) + sa*psi2(:norbs)
            c(:norbs,j) = (-sa*psi1(:norbs)) + ca*psi2(:norbs)
            k = norbs + 1
          end if
        end do
      end do
      sum1 = 0.D0
      psi1(:nocc) = 0.d0
      do i = 1, nocc
        do j = 1, numat
          il = nfirst(j)
          iu = nlast(j)
          x = 0.D0
          do k = il, iu
            x = x + c(k,i)**2
          end do
          sum1 = sum1 + x*x
          psi1(i) = psi1(i) + x*x
        end do
      end do
      if (sum > eps .and. iter < niter) go to 20
!
!   Check for LMOs that involve the same atom(s).  Resolve any
!   ill-definition.
!
      do i = 1, nocc
        psi1(i) = 1.d0/psi1(i)
      end do
      call resolv (c, cold, norbs, eig, nocc, psi1)
!
!   Work out LMO energy levels
!
      do i = 1, nocc
        sum = 0.D0
        do j = 1, nocc
          co = 0.D0
          do k = 1, norbs
            co = co + cold(k,j)*c(k,i)
          end do
          sum = sum + co*co*eig(j)
        end do
        eig1(i) = sum
      end do
!
!  Sort into increasing energy order
!
     do i = 1, nocc
        x = 100.D0
        i1 = 0
        do j = i, nocc
          if (x < eig1(j)) cycle
          x = eig1(j)
          i1 = j
        end do
        eig(i) = eig1(i1)
        x = eig1(i1)
        eig1(i1) = eig1(i)
        eig1(i) = x
        do j = 1, norbs
          x = c(j,i1)
          c(j,i1) = c(j,i)
          c(j,i) = x
        end do
      end do
      if (iprint == 1) then
        num_1 = char(ichar("1") +int(log10(numat + 0.05)))
        write (iw, "(/,10x,'NUMBER OF ITERATIONS =',i4,/,10x,'LOCALIZATION VALUE =',f14.9,/)") iter, sum1
        write (iw, "("//num_1//"x,'NUMBER OF CENTERS  LMO ENERGY     COMPOSITION OF ORBITALS ')")
        write (iw, "("//num_1//"x,34x,'(AS PERCENT OF THE LMO)',/)")
        num_1 = char(ichar("3") +int(log10(numat + 0.05)))
        num_2 = char(ichar("1") +int(log10(norbs + 0.05)))
        do i = 1, nocc
          x = 0.D0
          do k1 = 1, numat
            kl = nfirst(k1)
            ku = nlast(k1)
            dii = 0.D0
            do k = kl, ku
              dii = dii + c(k,i)**2
            end do
            x = x + dii*dii
            psi1(k1) = dii*100.D0
          end do
          x = 1.D0/x
          do ii = 1, numat
            sum = 0.D0
            k = 0
            do j = 1, numat
              if (psi1(j) < sum) cycle
              sum = psi1(j)
              k = j
            end do
            psi1(k) = 0.D0
            cii(ii) = sum
            iel(ii) = k
            if (sum >= 1.D0) cycle
            exit
          end do
          ii = ii - 1
          if (ii == 1) then
            if (cii(1) > 99.949d0) then
              write (iw, '(i'//num_1//',f10.4,f17.5, 3x,a2,i'//num_2//',f6.1)') i, x,eig(i), elemnt(nat(iel(1))),iel(1),cii(1)
            else
              write (iw, '(i'//num_1//',f10.4,f17.5, 3x,a2,i'//num_2//',f6.2)') i, x,eig(i), elemnt(nat(iel(1))),iel(1),cii(1)
            end if
          else
            if (ii < 6) then
              write (iw, '(i'//num_1//',f10.4,f17.5, 5(3x,a2,i'//num_2//',f6.2))') &
              i, x,eig(i), (elemnt(nat(iel(k))),iel(k),cii(k),k=1,ii)
            else if (ii < 11) then
              write (iw, '(i'//num_1//',f10.4,f17.5, 5(3x,a2,i'//num_2//',f6.2),/31x,5(3x,a2,i'//num_2//',f6.2))') &
              i, x,eig(i), (elemnt(nat(iel(k))),iel(k),cii(k),k=1,ii)
            else
              write (iw, '(i'//num_1//',f10.4,f17.5, 5(3x,a2,i'//num_2//',f6.2),2(/31x,5(3x,a2,i'//num_2//',f6.2)))') &
              i, x,eig(i), (elemnt(nat(iel(k))),iel(k),cii(k),k=1,ii)
            end if
          end if
        end do
        call phase_lock(c, norbs)
        if (numat > 25 .and. index(keywrd, "LARGE") == 0) then
          write(iw,'(/10x,a)')"Localized orbitals for systems of over 25 atoms are not printed by default"
          write(iw,'(10x,a)')"(To print localized orbitals, add keyword ""LARGE"")"
        else
          write (iw, "(/,/,20x,' LOCALIZED ORBITALS ',/,/)")
          call matout (c, eig, nocc, norbs, norbs)
        end if
        call to_screen("To_file: LMOs")
      end if
      if (txt == "c ") then
        if (nbeta == 0) then
          write (iw, '(/10X,''BONDING CONTRIBUTION OF EACH M.O.'',/)')
          call molval (c, p, 2.D0)
          call to_screen("To_file: Localized RHF M.O.s")
        else
          write (iw, '(/10X,''BONDING CONTRIBUTION OF EACH ALPHA M.O.'',/)')
          call molval (c, pa, 2.D0)
          call to_screen("To_file: Localized alpha M.O.s")
        end if
      else
        write (iw, '(/10X,''BONDING CONTRIBUTION OF EACH BETA  M.O.'',/)')
        call molval (c, pb, 2.D0)
        call to_screen("To_file: Localized beta M.O.s")
      end if
      if (index(keywrd,' GRAPH') == 0) then
        eig(:nocc) = refeig(:nocc)
        c(:norbs,:nocc) = cold(:norbs,:nocc)
      end if
      return
      end subroutine local
