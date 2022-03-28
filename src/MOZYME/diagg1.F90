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

subroutine diagg1 (fao, nocc, nvir, eigv, ws, latoms, ifmo, fmo, fmo_dim, nij, idiagg,  avir, aocc, aov)
   !**********************************************************************
   !
   !  FAO:      FOCK MATRIX OVER ATOMIC ORBITALS
   !  FMO:      FOCK MATRIX OVER MOLECULAR ORBITALS
   !  VECTOR:   EIGENVECTORS OF M.O.S
   !  FILLED:   OCCUPIED M.O.S
   !  eigs:     EIGENVALUES OF OCCUPIED M.O.S
   !  EMPTY:    VIRTUAL M.O.S
   !  EIGV:     EIGENVALUES OF VIRTUAL M.O.S
   !  NOCC:     NUMBER OF OCCUPIED M.O.S
   !  EIG:      ALL THE EIGENVALUES
   !  N:        NUMBER OF M.O.S = NUMBER OF A.O.S*
   !
   !**********************************************************************
   !
   !  ALL OCCUPIED AND UNOCCUPIED LMO ENERGY LEVELS ARE CALCULATED, AS WELL
   !  AS ALL MATRIX ELEMENTS WHICH CAN BE NON ZERO CONNECTING THE OCCUPIED
   !  AND VIRTUAL SETS OF LMO'S.
   !
   !**********************************************************************
   !
    use molkst_C, only: numat, norbs, mpack, numcal, keywrd
    use MOZYME_C, only : nvirtual, icocc_dim, &
       & nfmo, lijbo, nijbo, &
       tiny, sumt, ijc, ovmax, ncf, nce, nncf, nnce, ncocc, ncvir, &
     & iorbs, icocc, icvir, cocc, cvir
    use common_arrays_C, only : eigs, nfirst, nlast, p
    implicit none
    integer, intent (in) :: idiagg, nocc, nvir, fmo_dim
    integer, intent (inout) :: nij
    logical, dimension (numat), intent (out) :: latoms
    integer, dimension (2, fmo_dim), intent (inout) :: ifmo
    double precision, dimension (fmo_dim), intent (out) :: fmo
    double precision, dimension (numat), intent (out) :: aov
    double precision, dimension (norbs), intent (out) :: avir, ws
    double precision, dimension (mpack), intent (in) :: fao
    double precision, dimension (nvirtual), intent (inout) :: eigv
    double precision, dimension (icocc_dim), intent (out) :: aocc
    integer, save :: nf, icalcn = 0, mydisp = 0, ij0
    integer :: i, i1, i2, i4, ii, j, j1, j2, j4, &
         & jj, jl, jx, k, k1, kk, kl, l, loopi, loopj, kj, i1j1, &
         & i1j2, i2j1, i2j2, k1j1
    logical :: lij
    logical, save :: times
    double precision :: cutlim, flim
    double precision, save :: fref, oldlim, safety
    double precision :: cutoff, sum, sum1
    integer, external :: ijbo
    if (numcal /= icalcn) then
      icalcn = numcal
      fref = 10.0d0
      safety = 1.0d0
      oldlim = 0.0d0
      nf = 0
      times = (Index (keywrd, " TIMES") /= 0)
      if (Index (keywrd, " OLDENS") /= 0) then
        fref = 0.d0
      end if
    end if
    !
    !
    aocc(:) = 0.d0
    !
    !   IF THE CONTRIBUTION OF AN ATOM IN AN OCCUPIED LMO IS VERY SMALL
    !   THEN DO NOT USE THAT ATOM IN CALCULATING THE OCCUPIED-VIRTUAL
    !   INTERACTION.  PUT THE CONTRIBUTIONS INTO AN ARRAY 'AOCC'.
    !
    do j = 1, nocc
        loopj = ncocc(j)
        kl = 0
        do kk = nncf(j) + 1, nncf(j) + ncf(j)
          k1 = icocc(kk)
          sum = 0.d0
          do k = nfirst(k1), nlast(k1)
            kl = kl + 1
            sum = sum + cocc(kl+loopj) ** 2
          end do
          !
          !   AOCC(KK) HOLDS THE SQUARE OF THE CONTRIBUTION OF THE KK'TH ATOM
          !   IN THE OCCUPIED SET.  NOTE:  THIS IS NOT ATOM KK.
          !
          aocc(kk) = sum
          !
        end do
      !
    end do
    !
    !
    !    CUTLIM    PRECISION OF PL
    !
    !    1.D-6      0.004
    !    1.D-7      0.00004
    !
    cutlim = 1.d-8
    cutoff = Max (cutlim, tiny*10.d0*cutlim)
    flim = Min (3.d0, fref*0.5d0)
    fref = 0.d0
    if (idiagg <= 5) then
      cutoff = cutlim
    end if
    sumt = 0.d0
    ijc = 0
    tiny = 0.d0
    !
    !
    do i = 1, nvir
      !
      loopi = ncvir(i)
      latoms(:) = .false.
      l = 0
      do j = nnce(i) + 1, nnce(i) + nce(i)
        j1 = icvir(j)
        sum = 0.d0
        do k = l + 1, l + iorbs(j1)
          sum = sum + cvir(k+loopi) ** 2
        end do
        l = l + iorbs(j1)
        !
        !   AVIR(J1) HOLDS THE SQUARE OF THE CONTRIBUTION OF THE ATOM J1
        !   IN THE VIRTUAL LMO 'I'.
        !
        avir(j1) = sum
        latoms(icvir(j)) = .true.
      end do
      !
      if (lijbo) then
        do jj = nnce(i) + 1, nnce(i) + nce(i)
          j1 = icvir(jj)
          do jx = 1, iorbs(j1)
            ws(nfirst(j1)+jx-1) = 0.0d00
          end do
          !
          kl = loopi
          do kk = nnce(i) + 1, nnce(i) + nce(i)
            k1 = icvir(kk)
            kj = nijbo(k1, j1)
            if (kj >= 0) then
              if (avir(k1)*p(kj+1) > cutoff) then
                !
                !  EXTRACT THE ATOM-ATOM INTERSECTION OF FAO
                !
                if (iorbs(k1) .eq. 1 .and. iorbs(j1) .eq. 1) then
                  ws(nfirst(j1)) = ws(nfirst(j1)) + fao(kj+1) * cvir(kl+1)
                else
                  if (k1 > j1) then
                    ii = kj
                    do i4 = 1, iorbs(k1)
                      do jx = 1, iorbs(j1)
                        ii = ii + 1
                        ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) &
                             & * cvir(kl+i4)
                      end do
                    end do
                  else if (k1 < j1) then
                    ii = kj
                    do jx = 1, iorbs(j1)
                      do i4 = 1, iorbs(k1)
                        ii = ii + 1
                        ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) &
                             & * cvir(kl+i4)
                      end do
                    end do
                  else
                    do jx = 1, iorbs(j1)
                      do i4 = 1, iorbs(j1)
                        if (i4 > jx) then
                          ii = kj + (i4*(i4-1)) / 2 + jx
                        else
                          ii = kj + (jx*(jx-1)) / 2 + i4
                        end if
                        ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) &
                             & * cvir(kl+i4)
                      end do
                    end do
                  end if
                end if
              end if
            end if
            kl = kl + iorbs(k1)
          end do
        end do
      else
        do jj = nnce(i) + 1, nnce(i) + nce(i)
          j1 = icvir(jj)
          do jx = 1, iorbs(j1)
            ws(nfirst(j1)+jx-1) = 0.0d00
          end do
          !
          kl = loopi
          do kk = nnce(i) + 1, nnce(i) + nce(i)
            k1 = icvir(kk)
            kj = ijbo (k1, j1)
            if (kj >= 0) then
              if (avir(k1)*p(kj+1) > cutoff) then
                !
                !  EXTRACT THE ATOM-ATOM INTERSECTION OF FAO
                !
                if (iorbs(k1) == 1 .and. iorbs(j1) == 1) then
                  ws(nfirst(j1)) = ws(nfirst(j1)) + fao(kj+1) * cvir(kl+1)
                else
                  if (k1 > j1) then
                    ii = kj
                    do i4 = 1, iorbs(k1)
                      do jx = 1, iorbs(j1)
                        ii = ii + 1
                        ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) &
                             & * cvir(kl+i4)
                      end do
                    end do
                  else if (k1 < j1) then
                    ii = kj
                    do jx = 1, iorbs(j1)
                      do i4 = 1, iorbs(k1)
                        ii = ii + 1
                        ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) &
                             & * cvir(kl+i4)
                      end do
                    end do
                  else
                    do jx = 1, iorbs(j1)
                      do i4 = 1, iorbs(j1)
                        if (i4 > jx) then
                          ii = kj + (i4*(i4-1)) / 2 + jx
                        else
                          ii = kj + (jx*(jx-1)) / 2 + i4
                        end if
                        ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) &
                             & * cvir(kl+i4)
                      end do
                    end do
                  end if
                end if
              end if
            end if
            kl = kl + iorbs(k1)
          end do
        end do
      end if
      !
      do j = 1, numat
        if (latoms(j)) then
          sum = 0.d0
          do k = nfirst(j), nlast(j)
            sum = sum + ws(k) ** 2
          end do
          !
          !   AOV(J) HOLDS THE SQUARE OF THE ENERGY CONTRIBUTION OF THE
          !   ATOM J IN THE VIRTUAL LMO 'I'.
          !
          aov(j) = sum
        else
          aov(j) = 0.d0
        end if
      end do
      !
      !   EVALUATE THE VIRTUAL ENERGY LEVELS
      !
      sum = 0.d0
      kl = loopi
      do kk = nnce(i) + 1, nnce(i) + nce(i)
        k1 = icvir(kk)
        if (aov(k1)*avir(k1) > cutoff) then
          do k = nfirst(k1), nlast(k1)
            kl = kl + 1
            sum = sum + ws(k) * cvir(kl)
          end do
        else
          kl = kl + iorbs(k1)
        end if
      end do
      !
      eigv(i) = sum
      !
      !  EVALUATE THE OCCUPIED-VIRTUAL LMO INTERACTION ENERGIES
      !
      if (idiagg <= 5 .or. Mod (idiagg, 2) == 0) then
        nf = 0
        i1 = icvir(nnce(i)+1)
        if (nce(i) > 1) then
          i2 = icvir(nnce(i)+2)
        else
          i2 = i1
        end if
        if (lijbo) then
          do j = 1, nocc
            if (ijc /= nij) then
              !
              !  FAST TEST TO SEE IF THE INTEGRAL IS WORTH EVALUATING
              !
              j1 = icocc(nncf(j)+1)
              if (ncf(j) > 1) then
                j2 = icocc(nncf(j)+2)
              else
                j2 = j1
              end if
              !
              i1j1 = nijbo(i1, j1)
              i1j2 = nijbo(i1, j2)
              i2j1 = nijbo(i2, j1)
              i2j2 = nijbo(i2, j2)
              !
              if (i1j1 >= 0 .or. i1j2 >= 0 .or. i2j1 >= 0 .or. i2j2 >= 0) then
                sum = 0.d0
                if (i1j1 >= 0) then
                  sum = Abs (fao(i1j1+1))
                end if
                if (i1j2 >= 0) then
                  sum = sum + Abs (fao(i1j2+1))
                end if
                if (i2j1 >= 0) then
                  sum = sum + Abs (fao(i2j1+1))
                end if
                if (i2j2 >= 0) then
                  sum = sum + Abs (fao(i2j2+1))
                end if
                if (sum >= flim) then
                  loopj = ncocc(j)
                  lij = .false.
                  sum = 0.d0
                  kl = 0
                  do kk = nncf(j) + 1, nncf(j) + ncf(j)
                    k1 = icocc(kk)
                    if (aocc(kk)*aov(k1) < cutoff .or. .not. latoms(k1)) then
                      kl = kl + iorbs(k1)
                    else
                      lij = .true.
                      do k = nfirst(k1), nlast(k1)
                        kl = kl + 1
                        sum = sum + ws(k) * cocc(kl+loopj)
                      end do
                    end if
                  end do
                  sumt = sumt + Abs (sum)
                  tiny = Max (tiny, Abs (sum))
                  if (lij) then
                    if (Abs (sum) > oldlim) then
                      nf = nf + 1
                      ijc = ijc + 1
                      ifmo(1, ijc) = i
                      ifmo(2, ijc) = j
                      fmo(ijc) = sum
                    end if
                  end if
                end if
              end if
            end if
          end do
        else
          do j = 1, nocc
            if (ijc /= nij) then
              !
              !  FAST TEST TO SEE IF THE INTEGRAL IS WORTH EVALUATING
              !
              j1 = icocc(nncf(j)+1)
              if (ncf(j) > 1) then
                j2 = icocc(nncf(j)+2)
              else
                j2 = j1
              end if
              !
              i1j1 = ijbo (i1, j1)
              i1j2 = ijbo (i1, j2)
              i2j1 = ijbo (i2, j1)
              i2j2 = ijbo (i2, j2)
              !
              if (i1j1 >= 0 .or. i1j2 >= 0 .or. i2j1 >= 0 .or. i2j2 >= 0) &
                   & then
                sum = 0.d0
                if (i1j1 >= 0) then
                  sum = Abs (fao(i1j1+1))
                end if
                if (i1j2 >= 0) then
                  sum = sum + Abs (fao(i1j2+1))
                end if
                if (i2j1 >= 0) then
                  sum = sum + Abs (fao(i2j1+1))
                end if
                if (i2j2 >= 0) then
                  sum = sum + Abs (fao(i2j2+1))
                end if
                if (sum >= flim) then
                  loopj = ncocc(j)
                  lij = .false.
                  sum = 0.d0
                  kl = 0
                  do kk = nncf(j) + 1, nncf(j) + ncf(j)
                    k1 = icocc(kk)
                    if (aocc(kk)*aov(k1) < cutoff .or. .not. latoms(k1)) then
                      kl = kl + iorbs(k1)
                    else
                      lij = .true.
                      do k = nfirst(k1), nlast(k1)
                        kl = kl + 1
                        sum = sum + ws(k) * cocc(kl+loopj)
                      end do
                    end if
                  end do
                  sumt = sumt + Abs (sum)
                  tiny = Max (tiny, Abs (sum))
                  if (lij) then
                    if (Abs (sum) > oldlim) then
                      nf = nf + 1
                      ijc = ijc + 1
                      ifmo(1, ijc) = i
                      ifmo(2, ijc) = j
                      fmo(ijc) = sum
                    end if
                  end if
                end if
              end if
            end if
          end do
        end if
        nfmo(i) = nf
      else
        ij0 = ijc + mydisp
        do jj = 1, nfmo(i)
          if (ijc+mydisp == nij) exit
          j = ifmo(2, ij0+jj)
          loopj = ncocc(j)
          sum = 0.d0
          kl = 0
          do kk = nncf(j) + 1, nncf(j) + ncf(j)
            k1 = icocc(kk)
            if (aov(k1)*aocc(kk) < cutoff .or. .not. latoms(k1)) then
              kl = kl + iorbs(k1)
            else
              do k = nfirst(k1), nlast(k1)
                kl = kl + 1
                sum = sum + ws(k) * cocc(kl+loopj)
              end do
            end if
          end do
          sumt = sumt + Abs (sum)
          tiny = Max (tiny, Abs (sum))
          ijc = ijc + 1
          fmo(ijc) = sum
        end do
      end if
    end do
    !
    if (ijc == nij) then
      !
      !   THERE WAS NOT ENOUGH STORAGE TO HOLD ALL THE INTEGRALS.
      !   THEREFORE, ON THE NEXT ITERATION, CALCULATE FEWER INTEGRALS.
      !
      safety = safety * 2.d0
    else
      !
      !  THERE IS ENOUGH STORAGE FOR ALL THE INTEGRALS.  IF NECESSARY,
      !  CALCULATE MORE INTEGRALS.
      !
      safety = Max (safety*0.5d0, 1.d0)
    end if

    nij = ijc
    if (idiagg > 2 .and. Mod (idiagg, 4) /= 0) then
      fref = tiny ** 4
      if (times) then
        call timer (" AFTER DIAGG1 IN ITER")
      end if
        !
    else
      !
      !  EVALUATE THE OCCUPIED ENERGY LEVELS
      !
        !
      do i = 1, nocc
        !
          loopi = ncocc(i)
          l = 0
          do j = nncf(i) + 1, nncf(i) + ncf(i)
            j1 = icocc(j)
            sum = 0.d0
            do k = l + 1, l + iorbs(j1)
              sum = sum + cocc(k+loopi) ** 2
            end do
            l = l + iorbs(j1)
            !
            !   AOCC(KK) HOLDS THE SQUARE OF THE CONTRIBUTION OF THE KK'TH ATOM
            !   IN THE OCCUPIED SET.  NOTE:  THIS IS NOT ATOM KK.
            !
            !   AVIR(J1) HOLDS THE SQUARE OF THE CONTRIBUTION OF ATOM J1
            !   IN THE OCCUPIED SET.  NOTE:  THIS IS NOT THE J1'th ATOM.
            !
            avir(j1) = sum
            !
            !   AOCC(KK) HOLDS THE SQUARE OF THE CONTRIBUTION OF THE KK'TH ATOM
            !   IN THE OCCUPIED SET.  NOTE:  THIS IS NOT ATOM KK.
            !
            !   AVIR(J1) HOLDS THE SQUARE OF THE CONTRIBUTION OF ATOM J1
            !   IN THE OCCUPIED SET.  NOTE:  THIS IS NOT THE J1'th ATOM.
            !
          end do
          sum = 0.d0
          jl = loopi
          !
          if (lijbo) then
            do j = nncf(i) + 1, nncf(i) + ncf(i)
              j1 = icocc(j)
              kl = loopi
              do k = nncf(i) + 1, nncf(i) + ncf(i)
                k1 = icocc(k)
                k1j1 = nijbo(k1, j1)
                if (k1j1 >= 0) then
                  if (avir(k1)*p(k1j1+1)*aocc(k) >= cutoff) then
                    !
                    !  EXTRACT THE ATOM-ATOM INTERSECTION OF FAO
                    !
                    if (k1 > j1) then
                      !
                      !   LOWER TRIANGLE
                      !
                      do jx = 1, iorbs(j1)
                        sum1 = 0.d0
                        do i4 = 1, iorbs(k1)
                          ii = k1j1 + (i4-1) * iorbs(j1) + jx
                          sum1 = sum1 + fao(ii) * cocc(kl+i4)
                        end do
                        sum = sum + cocc(jl+jx) * sum1
                      end do
                    else if (k1 < j1) then
                      !
                      !   UPPER TRIANGLE
                      !
                      do jx = 1, iorbs(j1)
                        sum1 = 0.d0
                        do i4 = 1, iorbs(k1)
                          ii = k1j1 + (jx-1) * iorbs(k1) + i4
                          sum1 = sum1 + fao(ii) * cocc(kl+i4)
                        end do
                        sum = sum + cocc(jl+jx) * sum1
                      end do
                    else
                      !
                      !   DIAGONAL TERM
                      !
                      do jx = 1, iorbs(j1)
                        sum1 = 0.d0
                        do j4 = 1, jx
                          ii = k1j1 + (jx*(jx-1)) / 2 + j4
                          sum1 = sum1 + fao(ii) * cocc(kl+j4)
                        end do
                        ii = k1j1 + (jx*(jx+1)) / 2
                        do i4 = jx + 1, iorbs(k1)
                          ii = k1j1 + (i4*(i4-1)) / 2 + jx
                          sum1 = sum1 + fao(ii) * cocc(kl+i4)
                        end do
                        sum = sum + cocc(jl+jx) * sum1
                      end do
                    end if
                  end if
                end if
                kl = kl + iorbs(k1)
              end do
              !
              jl = jl + iorbs(j1)
            end do
          else
            do j = nncf(i) + 1, nncf(i) + ncf(i)
              j1 = icocc(j)
              kl = loopi
              do k = nncf(i) + 1, nncf(i) + ncf(i)
                k1 = icocc(k)
                k1j1 = ijbo (k1, j1)
                if (k1j1 >= 0) then
                  if (avir(k1)*p(k1j1+1)*aocc(k) >= cutoff) then
                    !
                    !  EXTRACT THE ATOM-ATOM INTERSECTION OF FAO
                    !
                    if (k1 > j1) then
                      !
                      !   LOWER TRIANGLE
                      !
                      do jx = 1, iorbs(j1)
                        sum1 = 0.d0
                        do i4 = 1, iorbs(k1)
                          ii = k1j1 + (i4-1) * iorbs(j1) + jx
                          sum1 = sum1 + fao(ii) * cocc(kl+i4)
                        end do
                        sum = sum + cocc(jl+jx) * sum1
                      end do
                    else if (k1 < j1) then
                      !
                      !   UPPER TRIANGLE
                      !
                      do jx = 1, iorbs(j1)
                        sum1 = 0.d0
                        do i4 = 1, iorbs(k1)
                          ii = k1j1 + (jx-1) * iorbs(k1) + i4
                          sum1 = sum1 + fao(ii) * cocc(kl+i4)
                        end do
                        sum = sum + cocc(jl+jx) * sum1
                      end do
                    else
                      !
                      !   DIAGONAL TERM
                      !
                      do jx = 1, iorbs(j1)
                        sum1 = 0.d0
                        do j4 = 1, jx
                          ii = k1j1 + (jx*(jx-1)) / 2 + j4
                          sum1 = sum1 + fao(ii) * cocc(kl+j4)
                        end do
                        ii = k1j1 + (jx*(jx+1)) / 2
                        do i4 = jx + 1, iorbs(k1)
                          ii = k1j1 + (i4*(i4-1)) / 2 + jx
                          sum1 = sum1 + fao(ii) * cocc(kl+i4)
                        end do
                        sum = sum + cocc(jl+jx) * sum1
                      end do
                    end if
                  end if
                end if
                kl = kl + iorbs(k1)
              end do
              !
              jl = jl + iorbs(j1)
            end do
          end if
          !
          eigs(i) = sum
          !
      end do
            oldlim = tiny * safety * 1.d-3
      fref = tiny ** 4
      if (times) then
        call timer (" AFTER DIAGG1 IN ITER")
      end if
    end if

    ovmax = tiny
end subroutine diagg1
