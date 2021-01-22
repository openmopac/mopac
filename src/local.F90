      subroutine local(c, nocc, eig, iprint, txt) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use molkst_C, only : norbs, numat, keywrd, nbeta
      use common_arrays_C, only : nat, nfirst, nlast, p, pa, pb
      use chanel_C, only : iw
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  10:47:23  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use resolv_I 
      use matout_I 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: nocc, iprint
      real(double)  :: c(norbs,norbs) 
      real(double)  :: eig(norbs) 
      real(double)  :: cold(norbs,norbs) 
      character*2 :: txt
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(20) :: iel 
      integer :: niter, i, j, iter, k, k1, kl, ku, il, iu, i1, ii 
      real(double), dimension(norbs) :: eig1, psi1, psi2, cii, refeig 
      real(double) :: eps, sum, xijjj, xjiii, xiiii, xjjjj, xijij, xiijj, dij, &
        dii, djj, aij, bij, ca, sa, sum1, x, co 
      character, dimension(99) :: elemnt*2 

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
!       J.J.P. STEWART, J.C.S. FARADAY (II) 77, 000, (1981).
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
      niter = 100 
      eps = 1.0D-7 
      refeig(:norbs) = eig(:norbs) 
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
          endif 
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
            endif 
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
          endif 
        end do 
      end do 
      sum1 = 0.D0 
      do i = 1, nocc 
        do j = 1, numat 
          il = nfirst(j) 
          iu = nlast(j) 
          x = 0.D0 
          do k = il, iu 
            x = x + c(k,i)**2 
          end do 
          sum1 = sum1 + x*x 
        end do 
      end do 
      if (sum > eps .and. iter < niter) go to 20 
!
!   Check for LMOs that involve the same atom(s).  Resolve any
!   ill-definition.
!
      call resolv (c, cold, norbs, eig, nocc) 
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
        write (iw, 110) iter, sum1 
  110 format(/,10x,'NUMBER OF ITERATIONS =',i4,/,10x,'LOCALISATION VALUE =',f&
        14.9,/) 
        write (iw, 120) 
  120 format(3x,'NUMBER OF CENTERS',14x,'(COMPOSITION OF ORBITALS)'/,/)   
       
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
          write (iw, 240) x, (elemnt(nat(iel(k))),iel(k),cii(k),k=1,ii) 
    240   format(f10.4,4(5(3x,a2,i3,f6.2),/,10x)) 
        end do 
  260 format(/,/,20x,' LOCALIZED ORBITALS ',/,/) 
        call phase_lock(c, norbs)
        write (iw, 260) 
        call matout (c, eig, nocc, norbs, norbs) 
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
      endif  
      if (index(keywrd,' GRAPH') == 0) then
        eig(:nocc) = refeig(:nocc) 
        c(:norbs,:nocc) = cold(:norbs,:nocc) 
      end if
      return  
      end subroutine local 
