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

      subroutine freqcy(fmatrx, freq, travel, force_const, eorc, deldip, ff, oldf, ts)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE molkst_C, only : numat, keywrd, nvar, mers, id
      USE common_arrays_C, only : atmass, coord, tvec
      USE funcon_C, only : fpc_10, fpc_8, fpc_6, pi
      USE chanel_C, only : iw, brillouin_fn, ibrz
      use to_screen_C, only : redmas, cnorml
      implicit none
      double precision  :: ff(*)
      logical , intent(in) :: eorc, ts
      double precision  :: fmatrx((3*numat*(3*numat+1))/2)
      double precision  :: freq(3*numat)
      double precision , intent(out) :: travel(3*numat), force_const(3*numat)
      double precision  :: deldip(3,3*numat)
      double precision , intent(out) :: oldf((3*numat*(3*numat+1))/2)
!
      integer :: loop, i, j, ij, iu, il, im1, ju, jl, ii, jj, l, linear, jii, k
      double precision, dimension(numat*3) :: wtmass
      double precision :: fact, c2pi, sumerr, sum, err, weight, summ, sum1
      logical :: bcc
!********************************************************************
!
!  FREQCY CALCULATES THE FORCE CONSTANTS AND VIBRATIONAL FREQUENCIES
!       FOR A MOLECULE.  IT USES THE ISOTOPIC MASSES TO WEIGHT THE
!       FORCE MATRIX
!
! ON INPUT   FMATRX   =  FORCE MATRIX, OF SIZE NUMAT*3*(NUMAT*3+1)/2.
!
!********************************************************************
      fact = fpc_10
!
!    CONVERSION FACTOR FOR SPEED OF LIGHT AND 2 PI.
!
      c2pi = 1.D0/(fpc_8*pi*2.D0)
! NOW TO CALCULATE THE VIBRATIONAL FREQUENCIES
!
!   FIND CONVERSION CONSTANTS FOR MASS WEIGHTED SYSTEM
      if (eorc .and. .not. ts .and. index(keywrd,' NOSYM') == 0) then
!
!     The diagonal terms of the Force Matrix should
!     be equal to minus the sum of the off-diagonal terms.
!
        do loop = 1, 20
          sumerr = 0.D0
          do i = 1, nvar
            sum = 0.D0
            err = 0.D0
            do j = 1, i - 1
              sum = sum + fmatrx((i*(i-1))/2+j)
            end do
            do j = i + 1, nvar
              sum = sum + fmatrx((j*(j-1))/2+i)
            end do
            err = err + fmatrx((i*(i+1))/2) + sum
            sumerr = sumerr + abs(err)
            fmatrx((i*(i+1))/2) = (-sum) - err*0.5D0
          end do
          call symt (fmatrx, deldip, ff)
          if (sumerr >= 1.D-6) cycle
          exit
        end do
      end if
      if (index(keywrd,' FREQCY') /= 0) then
        write (iw, '(A)') ' SYMMETRIZED HESSIAN MATRIX'
!
!   THE FORCE MATRIX IS PRINTED AS AN ATOM-ATOM MATRIX RATHER THAN
!   AS A 3N*3N MATRIX, AS THE 3N MATRIX IS VERY CONFUSING!
!
        ij = 0
        iu = 0
        do i = 1, numat
          il = iu + 1
          iu = il + 2
          im1 = i - 1
          ju = 0
          do j = 1, im1
            jl = ju + 1
            ju = jl + 2
            sum = 0.D0
            do ii = il, iu
              do jj = jl, ju
                sum = sum + fmatrx((ii*(ii-1))/2+jj)**2
              end do
            end do
            ij = ij + 1
            cnorml(ij) = sqrt(sum)
          end do
          ij = ij + 1
          cnorml(ij) = sqrt(fmatrx(((il+0)*(il+1))/2)**2 + &
                            fmatrx(((il+1)*(il+2))/2)**2 + &
                            fmatrx(((il+2)*(il+3))/2)**2 + &
                            2.D0*(fmatrx(((il+1)*(il+2))/2-1)**2 + &
                                  fmatrx(((il+2)*(il+3))/2-2)**2 + &
                                  fmatrx(((il+2)*(il+3))/2-1)**2))
        end do
        i = -numat
        call vecprt (cnorml, i)
      end if
      l = 0
      do i = 1, numat
        weight = 1.D0/sqrt(atmass(i))
        wtmass(l+1) = weight
        wtmass(l+2) = weight
        wtmass(l+3) = weight
        l = l + 3
        wtmass(l) = weight
      end do
!    CONVERT TO MASS WEIGHTED FMATRX
      linear = 0
      do i = 1, nvar
        if (i > 0) then
          oldf(linear+1:i+linear) = fmatrx(linear+1:i+linear)*1.D5
          fmatrx(linear+1:i+linear) = fmatrx(linear+1:i+linear)* &
            wtmass(i)*wtmass(:i)
          linear = i + linear
        end if
      end do
!
!    1.D5 IS TO CONVERT FROM MILLIDYNES/ANGSTROM TO DYNES/CM.

      if (Index (keywrd, " MERS") /= 0) then
!
!  Write out information for BZ.  Do NOT stop the run - the user might want
!  vibrational frequencies and thermodynamics of the cluster.
!
        bcc = (Index (keywrd, " BCC") /= 0)
        open (unit=ibrz, file=brillouin_fn, status="UNKNOWN", form="FORMATTED")
        rewind(ibrz)
        write (ibrz,*) numat*3, mers, bcc
        write (ibrz,*) (fmatrx(i), i=1, linear)
        write (ibrz,*) tvec, id, numat, ((coord(j, i)-coord(j,1), j=1, 3), i=1,numat)
        write (ibrz,*) ((i-1)*3+1,i*3, i=1, numat)
      end if
!
!    DIAGONALIZE
!
      if (.not. ts) call frame (fmatrx, numat, 1)
      call rsp (fmatrx, nvar, freq, cnorml)
      call phase_lock(cnorml, nvar)
      if (eorc .and. nvar == 3*numat) call symtrz (cnorml, freq, 2, .TRUE.)
      do i = 1, nvar
        j = int((freq(i)+50.D0)*0.01D0)
        freq(i) = freq(i) - dble(j*100)
      end do
      freq(:nvar) = freq(:nvar)*1.D5
!
!    CALCULATE REDUCED MASSES, STORE IN REDMAS
!
      do i = 1, nvar
        ii = (i - 1)*nvar
        summ = 0.D0
        do j = 1, nvar/3
          summ = summ + (cnorml(ii+j*3-2)**2 + &
                         cnorml(ii+j*3-1)**2 + &
                         cnorml(ii+j*3  )**2)**2*atmass(j)
        end do
        sum = 0.D0
        do j = 1, nvar
          jii = j + ii
          jj = (j*(j - 1))/2
          do k = 1, j
            sum = sum + cnorml(jii)*oldf(jj+k)*cnorml(k+ii)
          end do
          do k = j + 1, nvar
            sum = sum + cnorml(jii)*oldf((k*(k-1))/2+j)*cnorml(k+ii)
          end do
        end do
        sum = sum*0.5d0
        sum1 = sum*2.D0
        if (abs(freq(i)) > abs(sum)*1.D-20) then
          sum = 1.D0*sum/freq(i)
        else
          sum = 0.D0
        end if
      !
      !  Given   E = (h/2pi)freqency, then
      !          Energy = SQRT(FREQ(I)) = sqrt(force constant/reduced mass) and
      !          Reduced mass = (force constant)/(FREQ(I))
      !
      !   At this point, FREQ(I) is the square root of the energy.
      !
        redmas(i,1) = summ
        if (Abs (freq(i)) > Abs (sum)*1.d-20) then
          redmas(i, 2) = Abs (sum1/freq(i))
          sum = sum / freq(i)
        else
          redmas(i, 2) = 0.d0
          sum = 0.d0
        end if
        force_const(i) = freq(i)*redmas(i,1)*1.d-5
        freq(i) = sign(sqrt(fact*abs(freq(i)))*c2pi,freq(i))
!
! Convert frequency into SI units (ergs)
!
        sum = freq(i)*fpc_6*fpc_8
!
! Travel, in Angstroms
!
        if (force_const(i) == 0.D0) then
          travel(i) = 0.D0
        else
          travel(i) = sqrt(2.d0*sum/(force_const(i)*1.d5))*1.d8
        end if
        if (travel(i) > 1.D0) travel(i) = 0.D0
        if (abs(freq(i)) < abs(sum1)*1.D+20) then
          sum1 = sqrt(abs(freq(i)/(sum1*1.D-5)))
        else
          sum1 = 0.D0
        end if
      end do
      if (eorc) then
!
!    CONVERT NORMAL VECTORS TO CARTESIAN COORDINATES
!    (DELETED) AND NORMALIZE SO THAT THE TOTAL MOVEMENT IS 1.0 ANGSTROM.
!
        ij = 0
        do i = 1, nvar
          sum = 0.D0
          j = 0
          do jj = 1, nvar/3
            sum1 = 0.D0
            cnorml(ij+1) = cnorml(ij+1)*wtmass(j+1)
            sum1 = sum1 + cnorml(ij+1)**2
!
            cnorml(ij+2) = cnorml(ij+2)*wtmass(j+2)
            sum1 = sum1 + cnorml(ij+2)**2
!
            cnorml(ij+3) = cnorml(ij+3)*wtmass(j+3)
            sum1 = sum1 + cnorml(ij+3)**2
!
            j = j + 3
            ij = ij + 3
            sum = sum + sqrt(sum1)
          end do
          sum = 1.D0/sum
          ij = ij - nvar
          cnorml(ij+1:nvar+ij) = cnorml(ij+1:nvar+ij)*sum
          ij = nvar + ij
        end do
!
!          RETURN HESSIAN IN MILLIDYNES/ANGSTROM IN FMATRX
!
        fmatrx(:linear) = oldf(:linear)*1.D-5
      else
!
!  RETURN HESSIAN AS MASS-WEIGHTED FMATRIX
        linear = 0
!
        do i = 1, nvar
          if (i > 0) then
            fmatrx(linear+1:i+linear) = oldf(linear+1:i+linear)*1.D-5* &
              wtmass(i)*wtmass(:i)
            linear = i + linear
          end if
        end do
      end if
      return
      end subroutine freqcy
