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

  subroutine diagg2 (nocc, nvir,  eigv, iused, latoms, &
     & nij, idiagg, storei, storej)
   !***********************************************************************
   !
   !   DIAGG2  PERFORMS A SIMPLE JACOBIAN ANNIHILATION OF THE ENERGY TERMS
   !   CONNECTING THE OCCUPIED LMOS AND THE VIRTUAL LMOS.  THE ENERGY TERMS
   !   ARE IN THE ARRAY FMO, AND THE INDICES OF THE LMOS ARE IN IFMO.
   !
   !***********************************************************************
    use molkst_C, only: numat, norbs, numcal, keywrd, id
    use MOZYME_C, only : nvirtual, icocc_dim, shift, &
       & icvir_dim, cocc_dim, cvir_dim, ipad2, thresh, tiny, sumb, &
       ncf, nce, nncf, nnce, ncocc, ncvir, iorbs, cocc, cvir, icocc, icvir, &
       ifmo, fmo
!
    use common_arrays_C, only : eigs, nat
    use parameters_C, only: main_group
    use chanel_C, only: iw
    implicit none
    integer, intent (in) :: idiagg, nij, nocc, nvir
    logical, dimension (numat), intent (out) :: latoms
    integer, dimension (numat), intent (out) :: iused
    double precision, dimension (norbs), intent (out) :: storei, storej
    double precision, dimension (nvirtual), intent (in) :: eigv
    logical :: bug = .false.
    logical :: retry
    logical, save :: debug, times
    integer :: i, ii, jur, l
    integer, save :: icalcn = 0
    integer :: ij, ilr, incv, iur, j, jlr, jncf, k, le, lf, lij, loopi, loopj, &
   & mie, mle, mlee, mlf, mlff, ncei, ncfj, nrej
    double precision :: a, alpha, b, beta, biglim, c, d, e, sum
    double precision, save :: const, eps, eta, bigeps
    double precision, external :: reada
    integer, dimension (2) :: nrejct
    data nrejct / 2 * 0 /
    if (numcal /= icalcn) then
      icalcn = numcal
      times = (Index (keywrd, " TIMES") /= 0)
      debug = (Index (keywrd, " DIAGG2") /= 0)
      !
      !   IF THE SYSTEM IS A SOLID, THEN DAMP ROTATION OF VECTORS,
      !   IN AN ATTEMPT TO PREVENT AUTOREGENERATIVE CHARGE OSCILLATION.
      !
      i = Index (keywrd, " DAMP")
      if (i /= 0) then
        const = reada (keywrd, i+5)
      else if (id == 3) then
        const = 0.5d0
      else
        a = 1.d0               !  If a transition metal, set DAMP to 0.5d0
        do i = 1, numat        !
          if (.not. main_group(nat(i)))  a = 0.5d0
        end do
        const = a
      end if
      !
      !   EPS IS THE SMALLEST NUMBER WHICH, WHEN ADDED TO 1.D0, IS NOT
      !   EQUAL TO 1.D0
      call epseta (eps, eta)
      !
      !   INCREASE EPS TO ALLOW FOR A LOT OF ROUND-OFF
      !
      bigeps = 50.d0 * Sqrt (eps)
    end if
   !
   !  RETRY IS .TRUE. IF THE NUMBER OF REJECTED ANNIHILATIONS IS IN THE
   !         RANGE 1 TO 20  AND THE SAME ON TWO ITERATIONS.  THIS WILL
   !         OCCUR NEAR THE END OF A SCF CALCULATION, WHEN ONLY A FEW
   !         LMOS ARE BADLY BEHAVED.
   !
    retry = (nrejct(1) == nrejct(2) .and. nrejct(1) /= 0 .and. nrejct(1) < 20)
    if (Mod(idiagg, 5) == 0 .or. idiagg <= 5) then
      tiny = -1.d0
      biglim = -1.d0
    else
      tiny = 0.01d0 * tiny
      biglim = bigeps
    end if
   !***********************************************************************
   !
   !   DO A CRUDE 2 BY 2 ROTATION TO "ELIMINATE" SIGNIFICANT ELEMENTS
   !
   !***********************************************************************
    iused(:) = -1
    latoms(:) = .false.
    if (debug) then
      write (iw,*)
      write (iw,*) "            SIZE OF OCCUPIED ARRAYS IN DIAGG2"
      write (iw,*)
      write (iw,*) "    LMO    NNCF     NCF   SPACE  ", " NCOCC    SIZE   SPACE"
      do i = 1, nocc - 1
        l = ncocc(i)
        do j = nncf(i) + 1, nncf(i) + ncf(i)
          l = l + iorbs(icocc(j))
        end do
        write (iw, "(7I8)") i, nncf (i), ncf (i), nncf (i+1) - &
       & nncf(i) - ncf(i), ncocc(i), l - ncocc(i), ncocc(i+1) - l
      end do
      i = nocc
      l = ncocc(i)
      do j = nncf(i) + 1, nncf(i) + ncf(i)
        l = l + iorbs(icocc(j))
      end do
      write (iw, "(7I8)") i, nncf (i), ncf (i), icocc_dim - nncf (i) - ncf &
     & (i), ncocc(i), l - ncocc(i), cocc_dim - l
      write (iw,*)
      write (iw,*) "            SIZE OF VIRTUAL ARRAYS IN DIAGG2"
      write (iw,*)
      write (iw,*) "    LMO    NNCE     NCE   SPACE  ", " NCVIR    SIZE   SPACE"
      do i = 1, nvir - 1
        l = ncvir(i)
        do j = nnce(i) + 1, nnce(i) + nce(i)
          l = l + iorbs(icvir(j))
        end do
        write (iw, "(7I8)") i, nnce (i), nce (i), nnce (i+1) - nnce (i) - nce &
       & (i), ncvir(i), l - ncvir(i), ncvir(i+1) - l
      end do
      i = nvir
      l = ncvir(i)
      do j = nnce(i) + 1, nnce(i) + nce(i)
        l = l + iorbs(icvir(j))
      end do
      write (iw, "(7I8)") i, nnce (i), nce (i), icvir_dim - nnce (i) - nce (i), &
     & ncvir(i), l - ncvir(i), cvir_dim - l
      if (bug) then
        write (iw,*)
        write (iw, "(A,I3,A)") " THIS FAULT CAN PROBABLY BE " // &
                             & "CORRECTED BY USE OF KEYWORD 'NLMO=", &
                             & ipad2 + 50, "'"
        write (iw,*)
        call mopend ("VALUE OF NLMO IS TOO SMALL")
      end if
    end if
    sumb = 0.d0
    nrej = 0
    lij = 0
    outer_loop: do ij = 1, nij
      i = ifmo(1, ij)
      j = ifmo(2, ij)
      if (Abs (fmo(ij)) >= tiny) then
        c = fmo(ij) * const
        d = eigs(j) - eigv(i) - shift
        if (Abs (c/d) >= biglim) then
          ncfj = ncf(j)
          ncei = nce(i)
      !
      !  STORE LMOS FOR POSSIBLE REJECTION, IF LMOS EXPAND TOO MUCH.
      !
          jlr = ncocc(j) + 1
          if (j /= nocc) then
            jur = ncocc(j+1)
            jncf = nncf(j+1)
          else
            jur = cocc_dim
            jncf = icocc_dim
          end if
          jur = Min (jlr+norbs-1, jur)
          ilr = ncvir(i) + 1
          if (i /= nvir) then
            iur = ncvir(i+1)
            incv = nnce(i+1)
          else
            iur = cvir_dim
            incv = icvir_dim
          end if
          iur = Min (ilr+norbs-1, iur)
          l = 0
          do k = jlr, jur
            l = l + 1
            storej(l) = cocc(k)
          end do
          l = 0
          do k = ilr, iur
            l = l + 1
            storei(l) = cvir(k)
          end do
          !
          !   STORAGE DONE.
          !
          lij = lij + 1
          e = Sign (Sqrt(4.d0*c*c+d*d), d)
          alpha = Sqrt (0.5d0*(1.d0+d/e))
          do
            beta = -Sign (Sqrt(1.d0-alpha*alpha), c)
            sumb = sumb + Abs (beta)
            !
            ! IDENTIFY THE ATOMS IN THE OCCUPIED LMO.  ATOMS NOT USED ARE
            ! FLAGGED BY '-1' IN IUSED.
            !
            mlf = 0
            !
            do lf = nncf(j) + 1, nncf(j) + ncf(j)
              ii = icocc(lf)
              iused(ii) = mlf
              mlf = mlf + iorbs(ii)
            end do
            loopi = ncvir(i)
            loopj = ncocc(j)
            mle = 0
         !
         !      ROTATION OF PSEUDO-EIGENVECTORS
         !
            do le = nnce(i) + 1, nnce(i) + nce(i)
              mie = icvir(le)

              latoms(mie) = .true.
              mlff = iused(mie) + loopj
              if (iused(mie) >= 0) then
                !
                !  TWO BY TWO ROTATION OF ATOMS WHICH ARE COMMON
                !  TO OCCUPIED LMO J AND VIRTUAL LMO I
                !
                do mlee = mle + 1 + loopi, mle + iorbs(mie) + loopi
                  mlff = mlff + 1
                  a = cocc(mlff)
                  b = cvir(mlee)
                  cocc(mlff) = alpha * a + beta * b
                  cvir(mlee) = alpha * b - beta * a
                end do
              else
                !
                !   FILLED  LMO ATOM 'MIE' DOES NOT EXIST.
                !   CHECK IF IT SHOULD EXIST
                !
                sum = 0.d0
                do mlee = mle + 1 + loopi, mle + iorbs(mie) + loopi
                  sum = sum + (beta*cvir(mlee)) ** 2
                end do
                if (sum > thresh) then
                  !
                  if (nncf(j)+ncf(j) >= jncf) go to 1000
                  if (mlf+iorbs(mie)+loopj > jur) go to 1000
                  !
                  !  YES, OCCUPIED LMO ATOM 'MIE' SHOULD EXIST
                  !
                  ncf(j) = ncf(j) + 1
                  icocc(nncf(j)+ncf(j)) = mie
                  !
                  iused(mie) = mlf
                  mlf = mlf + iorbs(mie)
                  !
                  !   PUT INTENSITY INTO OCCUPIED LMO ATOM 'MIE'
                  !
                  mlff = iused(mie) + loopj
                  do mlee = mle + 1 + loopi, mle + iorbs(mie) + loopi
                    mlff = mlff + 1
                    cocc(mlff) = beta * cvir(mlee)
                    cvir(mlee) = alpha * cvir(mlee)
                  end do
                end if
              end if
              mle = mle + iorbs(mie)
            end do
            !
            !  NOW CHECK ALL ATOMS WHICH WERE IN THE OCCUPIED LMO
            !  WHICH ARE NOT IN THE VIRTUAL LMO, TO SEE IF THEY
            !  SHOULD BE IN THE VIRTUAL LMO.
            !
            do lf = nncf(j) + 1, nncf(j) + ncf(j)
              ii = icocc(lf)

              if ( .not. latoms(ii)) then
                sum = 0.d0
                do mlff = iused(ii) + loopj + 1, iused(ii) + loopj + &
                     & iorbs(ii)
                  sum = sum + (beta*cocc(mlff)) ** 2
                end do
                if (sum > thresh) then
                  if (nnce(i)+nce(i) >= incv) go to 1000
                  if (mle+iorbs(ii)+loopi > iur) go to 1000
                  !
                  !  YES, VIRTUAL  LMO ATOM 'II' SHOULD EXIST
                  !
                  nce(i) = nce(i) + 1
                  icvir(nnce(i)+nce(i)) = ii
                  latoms(ii) = .true.
                  !
                  !   PUT INTENSITY INTO VIRTUAL  LMO ATOM 'II'
                  !
                  mlff = iused(ii) + loopj
                  do mlee = mle + 1 + loopi, mle + iorbs(ii) + loopi
                    mlff = mlff + 1
                     !
                    cvir(mlee) = -beta * cocc(mlff)
                    cocc(mlff) = alpha * cocc(mlff)
                  end do
                  mle = mle + iorbs(ii)
                end if
              end if
            end do
          exit
      1000  continue
            nrej = nrej + 1
              !
              !   THE ARRAY BOUNDS WERE GOING TO BE EXCEEDED.
              !   TO PREVENT THIS, RESET THE LMOS.
              !
            l = 0
            do k = jlr, jur
              l = l + 1
              cocc(k) = storej(l)
            end do
            l = 0
            do k = ilr, iur
              l = l + 1
              cvir(k) = storei(l)
            end do
            ncf(j) = ncfj
            nce(i) = ncei
            do k = 1, numat
              iused(k) = -1
              latoms(k) = .false.
            end do
            if (retry) then
                !
                !   HALF THE ROTATION ANGLE.  WILL THIS PREVENT THE
                !   ARRAY BOUND FROM BEING EXCEEDED?
                !
              alpha = 0.5d0 * (alpha+1.d0)
            else
              cycle outer_loop
            end if
          end do
        !
        !  RESET COUNTERS WHICH HAVE BEEN SET.
        !
          do le = nnce(i) + 1, nnce(i) + nce(i)
            mie = icvir(le)
            latoms(mie) = .false.
          end do
        !
          do lf = nncf(j) + 1, nncf(j) + ncf(j)
            iused(icocc(lf)) = -1
          end do
        end if
      end if
    end do outer_loop
    nrejct(2) = nrejct(1)
    nrejct(1) = nrej
    if (times) then
      call timer (" AFTER DIAGG2 IN ITER")
    end if
  end subroutine diagg2
