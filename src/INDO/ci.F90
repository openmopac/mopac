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

      subroutine inttr (i, j, a, b, c, gamma, nbtmo, j12, k12, ifstor)

!     **** generates J and K integrals in MO basis for CI routines   ****
!     **** by index transformation on the AO - basis integrals	     ****
!     **** is set up for notation (ij||ab) = <ia||jb> of Georges notes****

!     **** J = [ij|ab] = <phi_i(1) phi_j(1) | 1/r | phi_a(2) phi_b(2)> ****
!     **** K = [ia|jb] = <phi_i(1) phi_a(1) | 1/r | phi_j(2) phi_b(2)> ****

!     **** integrals in the AO basis set are:				   ****
!     **** (mu, mu|nu, nu) = gamma(mu, nu) = CNDO Coulomb integrals		   ****
!     **** (mu, nu|mu, nu) = gamma(nu, mu) = INDO one - centre Exchange integrals ****
!     **** (mu, nu|nu, mu) = gamma(nu, mu) = INDO one - centre Exchange integrals ****
!     **** (ig1, ig2|ig3, ig4) INDO trans metal irregular integrals	   ****

!     **** nb, in an alternate notation used in Georges notes, we have *
!     **** J = <ib|ja> = <phi_i(1) phi_b(2) | 1/r | phi_j(1) phi_a(2)> ****
!     **** K = <ib|aj> = <phi_i(1) phi_b(2) | 1/r | phi_a(1) phi_j(2)> ****
!     **** ie, <ib||ja> = (ij||ab), and parameter definition reverses ****

!     **** J = <ij|ab> = <phi_i(1) phi_j(2) | 1/r | phi_a(1) phi_b(2)> ****
!     **** K = <ij|ba> = <phi_i(1) phi_j(2) | 1/r | phi_b(1) phi_a(2)> ****
!     **** ie, <ia||jb> = (ia||jb), and parameter definition reverses ****
      use reimers_C, only : n, matind, nham, ibf, iat, nslwr, nsupr, &
          istr, nsym, ixprd, ig1, ig2, ig3, ig4, g, ntrmet, itrmet, &
          nirreg, natt, ncore, edef, ee2, wk1, wk2, isc
!      USE chanel_C, only : iw

      implicit none
      double precision ::  c(n, n), gamma(n, n), aa, gk
      integer ::           nbtmo(n), a, b, i, j, k, &
                        ifstor, imu, imu0, ir, irreg, &
                        k1, k2, k3, k4, ka, kt, mu, mu1, &
                        nm1, nm2, nmc, nme, nu, nu1, nu1l
      logical*1 ::         diagon
      double precision ::  j12, k12, jk12
      double precision ::            edefs
      nu1l = 0
      nmc = 0
      nme = 0

      if (allocated(wk1)) deallocate(wk1)
      if (allocated(wk2)) deallocate(wk2)
      allocate(wk1(n))
      allocate(wk2(n))
      if (ifstor == 1) then
!	**** see if element already evaluated and stored ****
        nm1 = matind(max(i, j) - ncore) + min(i, j) - ncore
        nm2 = matind(max(a, b) - ncore) + min(a, b) - ncore
        nmc = matind(max(nm1, nm2)) + min(nm1, nm2)

        nm1 = matind(max(i, a) - ncore) + min(i, a) - ncore
        nm2 = matind(max(j, b) - ncore) + min(j, b) - ncore
        nme = matind(max(nm1, nm2)) + min(nm1, nm2)

        edefs = edef / 2.
        if (ee2(nme) < edefs .and. ee2(nmc) < edefs) then
!	  **** retriev from store ****
          j12 = ee2(nmc)
          k12 = ee2(nme)
          return
        end if
      end if

!     ********** evaluate integrals **************

      diagon = i == j .and. a == b
      j12 = 0.d0
      k12 = 0.d0

!     **** elements of symmetry i x j ****

      if (iand (nbtmo(i), nbtmo(j)) /= 0) then
        ir = ixprd(nsym(i), nsym(j))
        imu0 = 0
        do mu1 = nslwr(ir), nsupr(ir)
          mu = istr(1, mu1)
          wk1(mu) = c(i, mu)*c(j, mu)
          wk2(mu) = c(a, mu)*c(b, mu)
          j12 = j12 + wk1(mu)*wk2(mu) * gamma(mu1, mu1)
          do nu1 = nslwr(ir), mu1 - 1
!	     **** CNDO coulomb integrals (mu, mu|nu, nu) ****
             nu = istr(1, nu1)
             j12 = j12 + (wk1(mu)*wk2(nu) + wk1(nu)*wk2(mu)) *&
     &                                gamma(mu1, nu1)
          end do
          if (nham == 2) then
!	    **** INDO exchange terms (mu, nu|mu, nu) ****
            imu = iat(mu)
            if (imu /= imu0) then
              nu1l = mu1
              imu0 = imu
            else
              do nu1 = nu1l, mu1 - 1
                nu = istr(1, nu1)
                k12 = k12 + (wk1(mu)*wk2(nu) + wk1(nu)*wk2(mu))&
     &                        * gamma(nu1, mu1)
              end do
            end if
          end if
        end do
      end if

!     **** integrals of symmetry i x a ****

      if (iand (nbtmo(i), nbtmo(a)) /= 0) then
        ir = ixprd(nsym(i), nsym(a))
        jk12 = 0.d0
        imu0 = 0
        do mu1 = nslwr(ir), nsupr(ir)
          mu = istr(1, mu1)
          wk1(mu) = c(i, mu)*c(a, mu)
          wk2(mu) = c(j, mu)*c(b, mu)
          k12 = k12 + wk1(mu)*wk2(mu) * gamma(mu1, mu1)
          do nu1 = nslwr(ir), mu1 - 1
!	    **** CNDO coulomb integrals (mu, mu|nu, nu) ****
            nu = istr(1, nu1)
            k12 = k12 + (wk2(mu)*wk1(nu) + wk2(nu)*wk1(mu)) * &
     &                                gamma(mu1, nu1)
          end do
          if (nham == 2) then
!	    **** INDO exchange terms (mu, nu|mu, nu) ****
            imu = iat(mu)
            if (imu /= imu0) then
              nu1l = mu1
              imu0 = imu
            else
              do nu1 = nu1l, mu1 - 1
                nu = istr(1, nu1)
                jk12 = jk12 + (wk1(mu)*wk2(nu) + wk1(nu)*wk2(mu))&
     &                                * gamma(nu1, mu1)
              end do
            end if
          end if
        end do
        if (diagon) then
!	  **** this trick adds also the i x b terms for determinate energy ****
          j12 = j12 + jk12 * 2.d0
          k12 = k12 + jk12
        else
!	  **** this is the i x a term only ****
          j12 = j12 + jk12
        end if
      end if

!     **** integrals (mu, nu|nu, mu) from INDO, symmetry = i x b ****

      if (.not.diagon .and. nham == 2 .and.&
     &                iand (nbtmo(i), nbtmo(b)) /= 0) then
        jk12 = 0.d0
        imu0 = 0
        ir = ixprd(nsym(i), nsym(b))
        do mu1 = nslwr(ir), nsupr(ir)
          mu = istr(1, mu1)
          imu = iat(mu)
          wk1(mu) = c(i, mu)*c(b, mu)
          wk2(mu) = c(j, mu)*c(a, mu)
          if (imu /= imu0) then
            nu1l = mu1
            imu0 = imu
          else
            do nu1 = nu1l, mu1 - 1
              nu = istr(1, nu1)
              jk12 = jk12 + (wk2(mu)*wk1(nu) + wk2(nu)*wk1(mu))&
     &                                * gamma(nu1, mu1)
            end do
          end if
        end do
        j12 = j12 + jk12
        k12 = k12 + jk12
      end if

!      **** Zerners irregular integrals for transition metals ****

      do kt = 1, ntrmet
        ka = itrmet(kt)
        k = natt(ka)
!	**** select unique transition metal atom ****
        if (isc(ka) == 1) then
          mu = ibf(ka)
          do irreg = 1, nirreg
            k1 = ig1(irreg, k)
            k2 = ig2(irreg, k)
            k3 = ig3(irreg, k)
            k4 = ig4(irreg, k)
            gk = g(irreg, k)

!	    **** terms (k1, k2|k3, k4) ****
            aa = c(i, mu + k1) * c(b, mu + k4) * gk
            k12 = k12 + c(a, mu + k2) * c(j, mu + k3) * aa
            j12 = j12 + c(a, mu + k3) * c(j, mu + k2) * aa

!	    **** terms (k1, k2|k4, k3) ****
            aa = c(i, mu + k1) * c(b, mu + k3) * gk
            k12 = k12 + c(a, mu + k2) * c(j, mu + k4) * aa
            j12 = j12 + c(a, mu + k4) * c(j, mu + k2) * aa

!	    **** terms (k3, k4|k1, k2) ****
            aa = c(i, mu + k3) * c(b, mu + k2) * gk
            k12 = k12 + c(a, mu + k4) * c(j, mu + k1) * aa
            j12 = j12 + c(a, mu + k1) * c(j, mu + k4) * aa

!	    **** terms (k4, k3|k1, k2) ****
            aa = c(i, mu + k4) * c(b, mu + k2) * gk
            k12 = k12 + c(a, mu + k3) * c(j, mu + k1) * aa
            j12 = j12 + c(a, mu + k1) * c(j, mu + k3) * aa

            if (k1 /= k2) then
!	      **** terms (k2, k1|k3, k4) ****
              aa = c(i, mu + k2) * c(b, mu + k4) * gk
              k12 = k12 + c(a, mu + k1) * c(j, mu + k3) * aa
              j12 = j12 + c(a, mu + k3) * c(j, mu + k1) * aa

!	      **** terms (k2, k1|k4, k3) ****
              aa = c(i, mu + k2) * c(b, mu + k3) * gk
              k12 = k12 + c(a, mu + k1) * c(j, mu + k4) * aa
              j12 = j12 + c(a, mu + k4) * c(j, mu + k1) * aa

!	      **** terms (k3, k4|k2, k1) ****
              aa = c(i, mu + k3) * c(b, mu + k1) * gk
              k12 = k12 + c(a, mu + k4) * c(j, mu + k2) * aa
              j12 = j12 + c(a, mu + k2) * c(j, mu + k4) * aa

!	      **** terms (k4, k3|k2, k1) ****
              aa = c(i, mu + k4) * c(b, mu + k1) * gk
              k12 = k12 + c(a, mu + k3) * c(j, mu + k2) * aa
              j12 = j12 + c(a, mu + k2) * c(j, mu + k3) * aa
            end if

          end do
        end if
      end do

!     **** store matrix elements ****

      if (ifstor == 1) then
        ee2(nmc) = sngl (j12)
        ee2(nme) = sngl (k12)
        j12 = ee2(nmc)
        k12 = ee2(nme)
      end if

      return
      end

!     *************************************************************************

      function eec1 (i1, a1, c, gamma, nbtmo)

!     **** evaluates sum over core j of (i, a|j, j) with i, a in x space ****
      use reimers_C, only : n, matind, ncore, edef, eec
!      USE chanel_C, only : iw

      implicit none
      double precision ::  c(n, n), gamma(n, n), eec1
      integer ::           nbtmo(n)
      integer ::           a, a1, i, i1, j, nm
      double precision ::  j12, k12
      double precision ::            edefs

      i = i1 - ncore
      a = a1 - ncore
      nm = matind(max(i, a)) + min(i, a)
      edefs = edef / 2.
      if (eec(nm) > edefs) then
!	**** value not previously determined, calculate it ****
        eec(nm) = 0.0
        do j = 1, ncore
          call inttr (i1, a1, j, j, c, gamma, nbtmo, j12, k12, 0)
          eec(nm) = eec(nm) + 2.d0*j12 - k12
        end do
      end if

      eec1 = eec(nm)
      return

      end function eec1

!     *************************************************************************

      function align (i, j, spn, aocc, bocc)

!     **** calculates the number of column interchanges necessary to give ****
!     **** optimum alignment of two determinates.			  ****
!     **** returns 1.d0 if even nber, - 1.d0 if odd number		  ****
!     **** I, J are orbs that are different in spin type SPN		  ****
!
!     **** ELECTRON ORDERING IS: alpha(1) beta(1) alpha(2) beta(2) ...	  ****
!     **** ie, any electron for orbital i appears before any for orb i + 1  ****
!     **** which is different to Georges book which has all alpha els	  ****
!     **** followed by all beta els.  The consequence is that signs of    ****
!     **** linear combinations change, eg, closed shell i - >a singlet state ***
!     **** is |i, a_bar> - |i_bar, a> wheras in Goerges book its written as ****
!     ****    |i, a_bar> + |a, i_bar>					  ****
!     **** This choice has consequences only here and in S2 matrix gen	  ****
      use reimers_C, only : nov

      implicit none
      double precision ::  align
      integer ::           i, j, k, i1, j1, ninter
      logical*1 ::         aocc(nov), bocc(nov), spn

      i1 = min(i, j)
      j1 = max(i, j)
      ninter = 0
      do k = i1 + 1, j1 - 1
        if (aocc(k)) ninter = ninter + 1
        if (bocc(k)) ninter = ninter + 1
      end do
      if (spn) then
!	**** electron is alpha spin, move past a beta spin in orb i1 ****
        if (bocc(i1)) ninter = ninter + 1
      else
!	**** electron is beta spin, move past an alpha spin in orb j1 ****
        if (aocc(j1)) ninter = ninter + 1
      end if

      align = 1.d0
      if (mod(ninter, 2) == 1) align = -1.d0
      return

      end function align

!     **********************************************************************

      subroutine cimat (c, gamma, beta, nbtmo, istsym, e, aocc, bocc, &
     &                spintr, nspn, ci, ndump)
      use reimers_C, only : n, nb2, matind, nconf, multci, ncore, nol, nvl, &
           nvh, nov, fastci, edef, eec, ee2, iwk, wk0, nex, cc0, mspn
      use cosmo_C, only : iseps, diagsl

      implicit none
      double precision ::  e(nconf), beta(nb2), ci(nconf*(nconf + 1)/2), &
                        c(n, n), gamma(n, n), spintr(mspn, nconf), &
                        ciener, e1, e2, xxx, edefs
      integer ::      istsym(nconf), nspn(nconf), io(nconf), b, &
                        iv(nconf), i, j, k, iout, ioutt, iw = 0, &
                        k1, k2, ndump, nm, nxs, nxss
      double precision ::  j12, k12
      logical*1 ::         aocc(nov, mspn, nex), bocc(nov, mspn, nex), same
      integer ::           a, nbtmo(n)

!     ****  RHF: generates singlet or triplet CI matrix h ****
!     **** ROHF: generates CI matrix in (S) with spin Sz = MULTCI   ****
!     ****	 and hence S**2 components > = Sz.  The S**2 matrix ****
!     ****	 is written into matrix (T).			   ****

!     **** zero entire matrix ****
      do i = 1, matind(nconf) + nconf
        ci(i) = 0.d0
      end do

      if (fastci) then
!	************ RHF calc, singles excitations only ****************

!	**** CI matrix element of excitn i - >a with excitn j - >b ****
        do k2 = 2, nconf
!	  **** first, determine occ and vir levels in excitation ****
          io(k2) = 1
          do while (aocc(io(k2), 1, k2))
            io(k2) = io(k2) + 1
          end do
          iv(k2) = nvl - ncore
          do while (.not.aocc(iv(k2), 1, k2))
            iv(k2) = iv(k2) + 1
          end do

          nm = matind(k2) + k2
          ci(nm) = ci(nm) + e(k2)
          do k1 = 2, k2 - 1
            if (istsym(k2) == istsym(k1)) then
!	      **** matrix elements, J = [ij|ab], K = [ia|jb] ****
              call inttr (ncore + io(k2), ncore + io(k1), ncore + iv(k2), &
     &                    ncore + iv(k1), c, gamma, nbtmo, j12, k12, 0)
              nm = matind(k2) + k1
              ci(nm) = ci(nm) - j12
              if (multci == 1) ci(nm) = ci(nm) + 2.d0 * k12
            end if
          end do

        end do
      else

!	**************** General open shell matrix elements *************

        do k2 = 1, nconf
          nm = matind(k2) + k2
          ci(nm) = ci(nm) + e(k2)
          do k1 = 1, k2 - 1
            e1 = 0.d0
            if (istsym(k2) == istsym(k1)) then


!	      **** each config is a spin-adapted linear combs of dets	****
!	      **** if both sets of dets are the same, then E = 0 as these ****
!	      **** are both eigenvectors of the same orb class in SPCLAS ***
              same = nspn(k1)  ==  nspn(k2)
              i = 0
              do while (same .and. i < nov)
                i = i + 1
                same = (aocc(i, 1, k1) .eqv. aocc(i, 1, k2)) .and.&
     &                (bocc(i, 1, k1) .eqv. bocc(i, 1, k2))
              end do
              if (.not. same) then
                do i = 1, nspn(k1)
                  e2 = 0.d0
                  do j = 1, nspn(k2)
                    xxx = ciener (aocc(1, i, k1), &
     &                 bocc(1, i, k1), aocc(1, j, k2), bocc(1, j, k2), &
     &                 c, gamma, nbtmo, beta)
                    e2 = e2 + spintr(j, k2) * xxx
                  end do
                  e1 = e1 + spintr(i, k1) * e2
                end do
              end if
              nm = matind(k2) + k1
              ci(nm) = ci(nm) + e1

            end if
          end do
        end do

      end if

! RMG - add COSMO correction to CI diagonal matrix elements
      if (iseps) then
        do i = 1, n
          do j = 1, n
            cc0(i, j) = c(i, j)
          end do
        end do
        if (allocated(diagsl)) deallocate(diagsl)
        allocate(diagsl(nconf))
        call diagci(ci, aocc, bocc, diagsl)
      end if

!     ********************** dump output ******************

      i = 4
      if (ndump > 0) i = 3
      nxs = nov*(nov + 1)/2
      nxss = nxs*(nxs + 1)/2
      edefs = edef / 2.d0

      if (ndump >= 2) then
!       **** dump matrix elements differing by 1 electron ****
        write (iw, "(/' Contributions to integrals differing by 1 electron:', &
     &         ' from BETA, from sum over CI core, and total'/)")
        iout = 0
        ioutt = 0
        nm = 0
        do i = nol, nvh
          do a = nol, i
            nm = nm + 1
            if (eec(nm) < edefs) then
              ioutt = ioutt + 1
              iout = iout + 1
              iwk(1, iout) = i
              iwk(2, iout) = a
              wk0(1, iout) = beta(matind(i) + a)
              wk0(2, iout) = eec(nm)
            end if
            if (iout == 5 .or. iout > 0 .and. nm == nxs) then
              write (iw, '(5(2i4, 2f6.0, f6.1))') (iwk(1, k), iwk(2, k), wk0(1, k), wk0(2, k), &
     &            wk0(1, k) + wk0(2, k), k = 1, iout)
              iout = 0
            end if
          end do
        end do
        write (iw, "(1x, a3, ' used = ', i6, ' maximum = ', i6, ' fraction = ', f6.2)") &
          'eec', ioutt, nxs, ioutt/float(nxs)

!       **** dump matrix elements differing by 2 electrons ****
        write (iw, "(/' Non - Zero Contributions to integrals differing by 2', &
     &         ' electrons'/)")
        iout = 0
        ioutt = 0
        nm = 0
        do i = nol, nvh
         do j = nol, i
          do a = nol, j
           do b = nol, a
            nm = nm + 1
            if (ee2(nm) < edefs) then
              ioutt = ioutt + 1
              if (ee2(nm) /= 0.E0) then
                iout = iout + 1
                iwk(1, iout) = i
                iwk(2, iout) = j
                iwk(3, iout) = a
                iwk(4, iout) = b
                wk0(1, iout) = ee2(nm)
              end if
            end if
            if (iout == 6 .or. iout > 0 .and. nm == nxss) then
              write (iw, '(6(4i4, f6.2))') (iwk(1, k), iwk(2, k), iwk(3, k), iwk(4, k), &
     &            wk0(1, k), k = 1, iout)
              iout = 0
            end if
           end do
          end do
         end do
        end do
        write (iw, "(1x, a3, ' used = ', i6, ' maximum = ', i6, ' fraction = ', f6.2)") &
          'ee2', ioutt, nxss, ioutt/float(nxss)
      end if

      return
      end subroutine cimat

!     ***********************************************************************

      function ciener (aoi, boi, aoj, boj, c, gamma, nbtmo, beta)

!     **** determines the interaction energy between two necessarily	****
!     **** different dets. first is aoi, boi; second is aoj, boj		****
      use cosmo_C, only : bmat, nden, gden, ipiden, nps, qscnet, useps
      use funcon_C, only : a0, ev
      use reimers_C, only : n, nb2, matind, ncore, nol, nvh, nov, &
          ndiff, io, iv, cc0
      use molkst_c, only : mpack

      implicit none
      integer ::  nbtmo(n)
      double precision ::  beta(nb2), c(n, n), gamma(n, n), &
                        j12, k12, ciener, e, eec1, align
      integer ::           i, j, ivo, j1, ndiffo, ndiffv, a
      logical*1 ::         aoi(nov), boi(nov), aoj(nov), boj(nov), &
     &                  ao, bo, spno(2), spnv(2), ao1(nov), bo1(nov)
      double precision  ::  ediel, fcon, phi(nps), qden(nden), p(mpack)

      e = 0.d0

!      ************ determine nber of different electrons ************
!      **** ndiff = total, o = occ in k1 not k2; v = occ in k2 not k1****
!      **** spn = true if alpha ****
      ndiff = 0
      ndiffo = 0
      ndiffv = 0
      do j1 = nol, nvh
        j = j1 - ncore

!        **** check alpha occupancy ****
        if (aoi(j) .neqv. aoj(j)) then
          if (ndiff == 4) goto 120
          ndiff = ndiff + 1
          if (aoi(j)) then
            ndiffo = ndiffo + 1
            if (ndiffo > 2) goto 120
            io(ndiffo) = j1
            spno(ndiffo) = .true.
          else
            ndiffv = ndiffv + 1
            if (ndiffv > 2) goto 120
            iv(ndiffv) = j1
            spnv(ndiffv) = .true.
          end if
        end if

!       **** check beta occupancy ****
        if (boi(j) .neqv. boj(j)) then
          if (ndiff == 4) goto 120
          ndiff = ndiff + 1
          if (boi(j)) then
            ndiffo = ndiffo + 1
            if (ndiffo > 2) goto 120
            io(ndiffo) = j1
            spno(ndiffo) = .false.
          else
            ndiffv = ndiffv + 1
            if (ndiffv > 2) goto 120
            iv(ndiffv) = j1
            spnv(ndiffv) = .false.
          end if
        end if

      end do

      if (ndiff == 2) then
!	******* differ in one electron, excitn io(1) - >iv(1) ********

!	**** core term sum_j in core (io, iv|j, j) ****
        e = e + eec1 (io(1), iv(1), c, gamma, nbtmo)

!	**** (io, iv|j, j) terms with j in xcitn space ****
        do j1 = nol, nvh
          j = j1 - ncore
          ao = aoi(j)
          bo = boi(j)
          if (j1 == io(1) .or. j1 == iv(1)) then
!	    **** remove the electon thats changing from sum ****
            if (spno(1)) then
              ao = .false.
            else
              bo = .false.
            end if
          end if

          if (ao .or. bo) then
!	    **** calc <i, j||a, j> via J = [i, a|j, j] and K = [i, j|a, j] ****
            call inttr (io(1), iv(1), j1, j1, &
     &          c, gamma, nbtmo, j12, k12, 0)
            if (ao .and. bo) then
!		      **** orb k has both electrons occupied ****
              e = e + 2.d0*j12 - k12
            else if (ao .eqv. spno(1)) then
!              **** one electron, same spin as moving one ****
              e = e + j12 - k12
            else
!              **** one electron, opp spin to moving one ****
              e = e + j12
            end if
          end if
         end do

!        **** add one - electron term ****
         ivo = matind (max(io(1), iv(1))) + min(io(1), iv(1))
         e = e + beta(ivo)

!        **** determine sign due to reorg of det to line it all up ****
! Add solvent correction for electron - induced surface charges to beta(ivo)
         if (useps) then
           ediel = 0.d0
           fcon = a0 * ev
           a = 0
! Construct electron density for the product of orbitals
           do i = 1, n
             do j = 1, i
               a = a + 1
               p(a) = cc0(iv(1), i)*cc0(io(1), j)
             end do
           end do
! Set up charge density for solvation
           do i = 1, nden
             qden(i) = gden(i) * p(ipiden(i))
           end do
! Compute electrostatic potential at each surface point
           do i = 1, nps
             phi(i) = 0.d0
             do j = 1, nden
               phi(i) = phi(i) + bmat(j, i) * qden(j)
             end do
           end do
           do i = 1, nps
             ediel = ediel + qscnet(i, 2) * phi(i)
           end do
           e = e + ediel * fcon
          end if

!        **** determine sign due to reorg of det to line it all up ****
         e = e * align (iv(1) - ncore, io(1) - ncore, spno(1), aoi, boi)


       else
!        ********* dets differ by two electrons ***********

         if (spno(1) .neqv. spnv(1) ) then
!          **** ensure io(1) and iv(1) are same spin ****
           j = iv(1)
           iv(1) = iv(2)
           iv(2) = j
         end if
!	 *** calc <ij||ab> = <ib||aj> via J = [i, a|b, j] and K = [ib|aj] ****
        call inttr (io(1), iv(1), iv(2), io(2), c, gamma, nbtmo, &
      &              j12, k12, 0)
        e = j12
        if (spno(1) .eqv. spno(2)) e = e - k12

!       **** determine sign due to reorg of det, electron 1 ****
        e = e * align (iv(1) - ncore, io(1) - ncore, spno(1), aoi, boi)

!       **** determine change to det k1 after this rearrangement ****
        do i = 1, nov
          ao1(i) = aoi(i)
          bo1(i) = boi(i)
        end do
        if (spno(1)) then
          ao1(io(1) - ncore) = .false.
          ao1(iv(1) - ncore) = .true.
        else
          bo1(io(1) - ncore) = .false.
          bo1(iv(1) - ncore) = .true.
        end if

!       **** determine sign due to reorg of det, electron 2 ****
        e = e * align (iv(2) - ncore, io(2) - ncore, spno(2), ao1, bo1)

      end if

120   continue
      ciener = e
      return
      end function ciener

!     *************************************************************************

      subroutine cidiag (ciin, ci, aii, aocc, bocc)
      use reimers_C, only : matind, nconf, nci, nr, nov, nex, mspn
      use cosmo_C, only : iseps, diagsl

      implicit none
      double precision ::  ciin(nconf*(nconf + 1)/2), ci(nconf, nconf), &
                        aii(nconf), wk0(nconf), emin, xx
      integer ::           i, j, k, i0, ier, nm, nconf1
      logical*1 ::         aocc(nov, mspn, nex), bocc(nov, mspn, nex)


!     **** expand symmetric matrix into full storage mode ****
      nm = matind(nconf) + nconf
      do i = nconf, 1, - 1
        do j = i, 1, - 1
          ci(j, i) = ciin(nm)
          nm = nm - 1
        end do
      end do
      do i = nconf, 1, - 1
        do j = i, 1, - 1
          ci(i, j) = ci(j, i)
        end do
      end do

!     **** block - wise diagonalization ****
      i0 = 1
      do i = 1, nr
        if (nci(i) > 0) then

          call tred2e (nconf, nci(i), ci(i0, i0), aii(i0), wk0, ci(i0, i0))
          call tql2e  (nconf, nci(i),         aii(i0), wk0, ci(i0, i0), ier)
          i0 = i0 + nci(i)
        end if
      end do

! RMG - correct eigenvalues if using solvent
      if (iseps) then
        call corrci(ciin, ci, aocc, bocc, diagsl, aii)
      end if

!     **** order in terms of increasing eigenvalue ****
      do nconf1 = 1, nconf
        k = 0
        emin = 1.d30
        do i = nconf1, nconf
          if (emin > aii(i)) then
            k = i
            emin = aii(i)
          end if
        end do
        xx = aii(nconf1)
        aii(nconf1) = aii(k)
        aii(k) = xx
        do i = 1, nconf
          xx = ci(i, nconf1)
          ci(i, nconf1) = ci(i, k)
          ci(i, k) = xx
        end do
      end do
!     **** transpose eigenvectors ****
      do i = 2, nconf
        do j = 1, min(i - 1, nconf)
          xx = ci(i, j)
          ci(i, j) = ci(j, i)
          ci(j, i) = xx
        end do
      end do

      return
      end subroutine cidiag

!     *************************************************************************

      subroutine ciout (ci, evalci, istsym, dmci, aocc, bocc, tot)
      use reimers_C, only : n, matind, au2ev, au2ang, au2cm, debye, &
          nconf, ndtype, nmrep, nptg, nr, ixprd, nov, &
          n2phot, nese, nciout, nciouv, filenm, lenf, irrtyp, occ, &
          istate, nex, mspn
      USE chanel_C, only : iw
      use molkst_C, only : keywrd

      implicit none
      double precision :: mom, evalci(nconf), ci(nconf, nconf), &
                        dmci(nconf*(nconf + 1)/2, 3), &
                        beta, deldip, dipch, dlo, dmtot, eng, feps = 0.d0, &
                        fosc, foscf, foscm, fref = 0.d0, freq, &
                        rmom, solva, solve, tot, &
                        wl, cisum
      integer ::           istsym(nconf), &
                        i, j, k, i0, il, kk, &
                        nm, nm0, nmd, nst, nsx, nsy
      logical*1 ::       solv = .false.
      logical*1 ::         aocc(nov, mspn, nex), bocc(nov, mspn, nex)
      double precision ::  wrtconf
      double precision, external :: reada
      integer, allocatable :: nelc(:, :)
      allocate (nelc(n, nconf))
 1005 format (a5, 1x, i4, 1x, f7.4, '     CI coeff  CI percent')
 1006 format (4x, a7, 2x, i5, 1x, f12.8, f12.8)
 1007 format (a24, 7x, f12.8)
 1150 format (/' Depression of ground-state after CI=', &
     &         f12.0, ' cm**-1', ' Energy=', f14.7, ' eV' /&
     &         ' Dipole moment=', 3f10.6, ' tot=', f10.6, ' Debyes'/)
 1170 format (i4, i2, 1x, a3, f8.1, f10.6, f8.3)
      nr = 1
      nciouv = nciout
      nese = 1
      ndtype = 1
      if(allocated(occ)) deallocate(occ)
      allocate(occ(nov))

!     **** determine CI state symmetry ****
      if(allocated(istate)) deallocate(istate)
      allocate(istate(nconf))
      do il = 1, nconf
        j = 1
        do while (abs(ci(il, j)) < 1.d-2)
          j = j + 1
        end do
        nsy = istsym(j)
        istate(il) = nsy
      end do

!     **** determine MO occupancy CI states ****
      do j = 1, nconf
!	**** determine nber of electrons in MO i in all configs j ****
        do i = 1, nov
          nelc(i, j) = 0
          if (aocc(i, 1, j)) nelc(i, j) = 1
          if (bocc(i, 1, j)) nelc(i, j) = nelc(i, j) + 1
        end do
      end do
      do j = 1, nciout
!	**** weight these by CI coeffs ****
        do i = 1, nov
          occ(i) = 0.d0
          do k = 1, nconf
            occ(i) = occ(i) + nelc(i, k) * ci(j, k)**2
          end do
        end do
      end do

!     **** ground state properties ****
      dmtot = sqrt (dmci(1, 1)**2 + dmci(1, 2)**2 + dmci(1, 3)**2) * debye
      dlo = evalci(1)/au2ev*au2cm
      write (iw, 1150) dlo, evalci(1) + tot, (debye*dmci(1, kk), kk = 1, 3), dmtot

!     **** read solvent shift parameters ****

      solva = 0.d0
      solve = 0.d0
      foscf = 2.d0/3.d0/au2ang**2/au2cm
      foscm = 2.d0/3.d0*au2ang**2*au2cm

!     *************** output excitations from CI states ****************
! RMG - add option to write transition dipoles between excited states
      if (index(keywrd, ' TDIP') /= 0) then
        filenm(lenf:lenf + 4) = '.tdip'
        open (14, file = filenm, status = 'unknown')
        write(14, *) nciout
      end if
! End RMG
      do i0 = 1, nese

        nsx = istate(i0)
        write (iw, "(/'  CI trans.  energy frequency wavelength oscillator', &
     & '--------- polarization---------  dipole  ', &
     & '------ components-----')")
        write (iw, "(' st.  symm.    eV      cm - 1       nm      strength ', &
     & '      x          y         z       moment     x       y  ', &
     & '     z ' /)")
        call polzro
        nm0 = matind(i0) + i0

        do il = 1, nciout
         if (il /= i0) then
          eng = evalci(il) - evalci(i0)
          freq = eng / au2ev * au2cm
! RMG - add option to write transition dipoles between excited states
          if (index(keywrd, ' TDIP') /= 0) then
            write(14, *) eng
          end if
! End RMG

          wl = 0.d0
          if (abs(freq) > 1.d0) wl = 1.D7/abs(freq)
          nsy = istate(il)
          nst = ixprd(nsx, nsy)

!         **** dipole moments and solvent shifts ****
          nmd = matind(il) + il
          dmtot = sqrt (dmci(nmd, 1)**2 + dmci(nmd, 2)**2 + dmci(nmd, 3)**2)&
     &                * debye
          if (i0 == 1 .and. solv) then
            beta = 0.d0
            solva = 0.d0
            solve = 0.d0
            do kk = 1, 3
              deldip = dmci(nmd, kk) - dmci(nm0, kk)
              beta = beta + deldip**2
              solva = solva + dmci(nm0, kk)*deldip
              solve = solve + dmci(nmd, kk)*deldip
            end do
            beta =  fref * beta * 0.5
            solva = -feps * solva - beta
            solve = -feps * solve + beta

!	    **** Jeffs special code to evaluate polarizability after shift ****
            evalci(il) = evalci(il) + solva * au2ev / au2cm
            eng = evalci(il) - evalci(i0)
            freq = eng / au2ev * au2cm
            wl = 0.d0
            if (abs(freq) > 1.d0) wl = 1.D7/abs(freq)
          end if

!	  **** transition moments and oscillator strengths ****
          nm = matind(max(i0, il)) + min(i0, il)
          mom = dmci(nm, 1)**2 + dmci(nm, 2)**2 + dmci(nm, 3)**2
!	  write (iy, '(3i4, 4f10.5)') i0, il, nm, sqrt(mom)*debye,
!     $				   (dmci(nm, k)*debye, k = 1, 3)
          if (ndtype /= 2) then
!	    **** f = 2/3*m**2*e with m, e in au; with m in Ang and e in cm** - 1
            fosc = foscf * freq * mom
            if (ndtype == 3) fosc = fosc * 1000.d0
          else
!	    **** f = 2/3*m**2/e with m, e in au; with m in Ang** - 1 and e in cm** - 1
            fosc = foscm / freq * mom
          end if
          call polizn (dmci(nm, 1), dmci(nm, 2), dmci(nm, 3), freq)

          dipch = sqrt ((dmci(nmd, 1) - dmci(1, 1))**2 + &
     &      (dmci(nmd, 2) - dmci(1, 2))**2 + (dmci(nmd, 3) - dmci(1, 3))**2)
          if (irrtyp > 0 .and. freq < 4.D5) write (40, 1170)&
     &      il - 1, nst, nmrep(nst, nptg), freq, fosc, dipch

          if (il > i0 .and. fosc > 0.00000001d0) then

!	    **** main output to log file ****
            rmom = 1.d0 / sqrt(mom)
            if (ndtype <= 1) then
              write (iw, "(i4, 1x, a3, f11.7, f10.0, f10.2, f10.6, 3f11.6, f10.6, 3f8.3)") &
                il, nmrep(nst, nptg), eng, freq, wl, fosc, &
     &                (dmci(nm, kk)*rmom, kk = 1, 3), &
     &                dmtot, (debye*dmci(nmd, kk), kk = 1, 3)
            else
              write (iw, "(i4, 1x, a3, f11.7, f10.0, f10.2, f10.6, 3f11.6, f10.6, 3f8.3)") &
                il, nmrep(nst, nptg), eng, freq, wl, fosc, &
     &                (dmci(nm, kk)*rmom, kk = 1, 3)
            end if

          else if (il > i0) then
!	    **** main output to log file, no osc strength ****
            if (ndtype <= 1) then
              write (iw, '(i4, 1x, a3, f11.7, f10.0, f10.2, f10.6, 33x,  f10.6, 3f8.3)') &
                il, nmrep(nst, nptg), eng, freq, wl, fosc, &
     &                dmtot, (debye*dmci(nmd, kk), kk = 1, 3)
            else
              write (iw, '(i4, 1x, a3, f11.7, f10.0, f10.2, f10.6, 33x,  f10.6, 3f8.3)') &
                il, nmrep(nst, nptg), eng, freq, wl, fosc
            end if
          end if

         end if
         end do

        if (ndtype <= 1) call polout
        if (n2phot > 0) call twopho (i0, evalci, dmci, istate)
      end do

!     ************************ write out ci matrix **********************
!     RMG add-print major state contributions for all states

      wrtconf = 0.04d0
      if (index(keywrd, ' WRTCONF=') /= 0) then
        wrtconf = reada(keywrd, index(keywrd, ' WRTCONF='))
      end if

      write (iw, *)
      write (iw, *)
      write (iw, *) 'Major CI contributions to CI states'
      write (iw, *)

      do il = 1, nciout
        eng = evalci(il) - evalci(1)
        write (iw, 1005) 'State ', il, eng

!       Find configurations with large contributions
!
        cisum = 0
        do i = 1, nconf
          if (abs(ci(il, i)) > wrtconf) then
            write(iw, 1006) 'Config ', i, ci(il, i), ci(il, i)**2
            cisum = cisum + ci(il, i)**2
          end if
        end do
        write (iw, 1007) 'Total coeff printed    ', cisum
        write (iw, *) ' '
      end do
!     End RMG

! RMG - add option to write transition dipoles between excited states
      if (index(keywrd, ' TDIP') /= 0) then
        do i = 1, nciout
          do j = 1, i
            nm = matind(max(i, j)) + min(i, j)
            write(14, '(i5, i5, 3f12.6)') i, j, dmci(nm, 1)*debye, dmci(nm, 2)*debye, dmci(nm, 3)*debye
          end do
         end do
        close(14)
      end if

! End RMG

      return

      end subroutine ciout

!     *************************************************************************

      subroutine foscil (dmci, dm, istsym, e, aocc, bocc, spintr, nspn, xz, zcore)

!     **** uses dipole moment in MO basis to construct dipole integrals ****
!     **** between CI configurations					****
      use reimers_C, only : na, nb2, matind, au2ev, au2ang, au2cm, debye, &
          icifrag, nconf, multci, ncore, nov, fastci, n2phot, nciouv, &
          nmrep, nptg, aor1, bor1, dipsym, ndtype, nex, mspn
      USE chanel_C, only : iw

      implicit none
      double precision ::  dmci(nconf*(nconf + 1)/2, 3), dm(nb2, 3), &
                        e(nconf), spintr(mspn, nconf), xz(na, 3), &
                        zcore(na), dcore(3), sqrt2, foscf, foscm, &
                        xx, align, emin, freq, diplen, tmsq, f
      integer ::           istsym(nconf), nspn(nconf), &
                        i, j, k, i1, ilow, iov, ivo, j1, k1, k2, kk, &
                        maxx, ndiff, nel, nio, niv, nm, nspni
      logical*1 ::         aocc(nov, mspn, nex), bocc(nov, mspn, nex), &
     &                  same, adiff, bdiff, abdiff
      character*4       lo(9), lv(9), ll
      character*3       rep

      write (iw, "(' The lowest ', i4, ' spin-adapted configurations of ', &
     &        'multiplicity=', i3 /&
     &        5x, 'sym   eV   cm**-1 -dets- dipole oscilator X FRAG ....', &
     &        'Excitations named from first reference determinate' /&
     &        23x, 'tot  #  Debye  strength')") nconf, multci
      nciouv = nconf

      if (ndtype == 3) write (iw, "(45x, ' x 1000')")
      write (iw, *)
      call polzro
      sqrt2 = sqrt (2.d0)
      foscf = 2.d0/3.d0/au2ang**2/au2cm
      foscm = 2.d0/3.d0*au2ang**2*au2cm

      if (ndtype == 4) then
!	**** dont do any dipole or transition moment calcs ****
        do i = 1, 3
          do j = 1, nconf*(nconf + 1)/2
            dmci(j, i) = 0.d0
          end do
        end do
        goto 500
      end if

!     **** determine dipole moment of (nucleas + inner shells) + valence core ****

      do kk = 1, 3
        dcore(kk) = 0.d0
         if (abs(dipsym - 1.d0) < 1.d-10) then
          do i = 1, na
            dcore(kk) = dcore(kk) + xz(i, kk) * zcore(i)
          end do
          do i = 1, ncore
            dcore(kk) = dcore(kk) + dm(matind(i) + i, kk)*2.d0
          end do
        end if
      end do

!     ***** calc dipole and transition moments ****

      nm = 0
      do i = 1, nconf
!	**** transition moments ****
        do j = 1, i - 1
          nm = nm + 1
          dmci(nm, 1) = 0.d0
          dmci(nm, 2) = 0.d0
          dmci(nm, 3) = 0.d0

!	  **** loop over each det in each config ****
          do i1 = 1, nspn(i)
            do j1 = 1, nspn(j)

!	      **** count the nber of electrons differing, proceed if 1 only ****
              ndiff = 0
              k = 0
              do while (ndiff <= 3 .and. k < nov)
                k = k + 1
                adiff = aocc(k, i1, i) .neqv. aocc(k, j1, j)
                bdiff = bocc(k, i1, i) .neqv. bocc(k, j1, j)
                if (adiff .and. bdiff) then
                  ndiff = 100
                else if (adiff .or. bdiff) then
                  ndiff = ndiff + 1
                  k2 = k
                  if (ndiff == 1) k1 = k
                  abdiff = adiff
                end if
              end do

              if (ndiff == 2) then
!		**** determine sign due to reorg of det to line it all up ****
                xx = align (k1, k2, abdiff, aocc(1, i1, i), bocc(1, j1, j)) *&
     &                spintr(i1, i) * spintr(j1, j)
!		**** check if dip ints real or imag; get symm of matrix ****
                ivo = matind(ncore + k2) + ncore + k1
                do kk = 1, 3
                  dmci(nm, kk) = dmci(nm, kk) + xx * dm(ivo, kk)
                end do
              end if

            end do
          end do

          if (fastci .and. j == 1) then
!	    **** GS interaction; apply directly + / - linear combs of det pair ****
            do kk = 1, 3
              dmci(nm, kk) = sqrt2 * dmci(nm, kk)
              if (multci == 3) dmci(nm, kk) = 0.d0
            end do
          end if
        end do

!	**** dipole moment of CI state i ****
        nm = nm + 1
        dmci(nm, 1) = dcore(1)
        dmci(nm, 2) = dcore(2)
        dmci(nm, 3) = dcore(3)
        if (abs(dipsym - 1.d0) < 1.d-10) then
          do k = 1, nov
            nel = 0
            if (aocc(k, 1, i)) nel = 1
            if (bocc(k, 1, i)) nel = nel + 1
            if (nel > 0) then
              iov = matind(ncore + k) + ncore + k
              do kk = 1, 3
                dmci(nm, kk) = dmci(nm, kk) + nel * dm(iov, kk)
              end do
            end if
          end do
        end if
      end do
500   continue

!     *************** generate names and output single - dets **************

!     **** first, calc max nber of excitations for format spec ****

      maxx = 0
      do i = 1, nconf
        nio = 0
        niv = 0
         do j = 1, nov
           if (aor1(j) .neqv. aocc(j, 1, i)) then
            if (aor1(j)) then
              nio = nio + 1
            else
              niv = niv + 1
            end if
          end if
          if (bor1(j) .neqv. bocc(j, 1, i)) then
            if (bor1(j)) then
              nio = nio + 1
            else
              niv = niv + 1
            end if
          end if
        end do
        if (nio /= niv) then
          write (iw, *) 'Cant interpret excitation', i
          write (iw, *) 'ref = ', (aor1(j), j = 1, nov), '  ', (bor1(j), j = 1, nov)
          write (iw, *) 'det = ', (aocc(j, 1, i), j = 1, nov), '  ', &
     &                      (bocc(j, 1, i), j = 1, nov)
          stop 'DET FUNNY in FOSCIL'
        end if
        maxx = max (maxx, nio)
      end do
      if (maxx > 9) stop 'EXCITE - OUTPUT IS FOR MAX OF 9 EXCITNS !'

!     *** determine state of lowest energy; calc transition moments from it ****

      emin = 1.D30
      do i = 1, nconf
        if (emin > e(i)) then
          emin = e(i)
          ilow = i
        end if
      end do

!     ***** actually do output *****

      do i = 1, nconf
!	**** first reset flags for Occ and Vir orbs ****
        nio = 0
        niv = 0
        do j = 1, maxx
          lo(j) = '    '
          lv(j) = '    '
        end do

!	**** loop over orbitals, check diffs from first ref det ****
        do j = 1, nov

          if (aor1(j) .neqv. aocc(j, 1, i)) then
!	    **** alpha differs from first ref det ****
            if (aor1(j)) then
              nio = nio + 1
              write (ll, '(i4)') j + ncore
              lo(nio) = ll
            else
              niv = niv + 1
              write (ll, '(i4)') j + ncore
              lv(niv) = ll
            end if
          end if

          if (bor1(j) .neqv. bocc(j, 1, i)) then
!	    **** beta differs from first ref det ****
            if (bor1(j)) then
              nio = nio + 1
              write (ll, '(i4)') j + ncore
              lo(nio) = ll
            else
              niv = niv + 1
              write (ll, '(i4)') j + ncore
              lv(niv) = ll
            end if
          end if
        end do

!	**** check if det list is same as any previous config ****
        nspni = 1
        do k = 1, i - 1
          same = nspn(k)  ==  nspn(i)
          j = 0
          do while (j < nov .and. same)
            j = j + 1
            same = (aocc(j, 1, i) .eqv. aocc(j, 1, k)) .and.&
     &          (bocc(j, 1, i) .eqv. bocc(j, 1, k))
          end do
          if (same) nspni = nspni + 1
        end do

!	**** symmetry, wavenumber, dip length, tr mom to ref state ****
        rep = nmrep(istsym(i), nptg)
        freq = (e(i) - e(ilow))/au2ev*au2cm
        nm = matind(i) + i
        diplen = sqrt (dmci(nm, 1)**2 + dmci(nm, 2)**2 + dmci(nm, 3)**2)

!	**** calculate ground-state polarizability ****
        nm = matind(max(i, ilow)) + min(i, ilow)
        if (ndtype <= 1)&
     &        call polizn (dmci(nm, 1), dmci(nm, 2), dmci(nm, 3), freq)
!	**** calculate oscillator strength ****
        tmsq = dmci(nm, 1)**2 + dmci(nm, 2)**2 + dmci(nm, 3)**2
        if (ndtype /= 2) then
!	  **** f = 2/3*m**2*e with m, e in au; with m in Ang and e in cm** - 1
          f = foscf * tmsq * freq
          if (ndtype == 3) f = f * 1000.d0
        else
!	  **** velocity formalism ****
!	  **** f = 2/3*m**2/e with m, e in au; with m in Ang** - 1 and e in cm** - 1
          f = foscm * tmsq / freq
        end if

        write (iw, "(i4, 1x, a3, f7.3, f8.0, 2i3, f8.4, f9.6, i2, i3, &
     &        ' (', 20a4)") i, rep, e(i), freq, nspn(i), nspni, diplen*debye, &
     &   f, nio, icifrag(i), (lo(j), j = 1, maxx), ')->(', (lv(j), j = 1, maxx), ')'

      end do

!     **** output of estimated GS polarizability ****

!     **** two - photon intensities ****
      if (n2phot > 0) call twopho (ilow, e, dmci, istsym)

      return
      end subroutine foscil

!     *************************************************************************

      subroutine twopho (i0, ener, dm, istsym)

!     **** this calculates 2 - photon relative einstein - B coefficients	****
!     **** for absorption from state i0 which must be non - degenerate	****
!     **** D.P. Craig and T. Thirunamachandran Molecular Quantum	****
!     **** Electrodynamics (Academic, London, 1984) P109		****
      use reimers_C, only : matind, nconf, nptg, nmrep, ixprd, nciout
      USE chanel_C, only : iw

      implicit none
      double precision ::   ener(nconf), dm(nconf*(nconf + 1)/2, 3), &
                         alp11, alp21, alp22, alp31, alp32, alp33, &
                         d, trace, od, planar, circ
      integer ::            istsym(nconf), j, k, i0, isym, nj, nk
      write (iw, "('1Single - Beam 2 - photon relative Einstein B coefficients', &
     &        ' exciting from state', i4 /&
     &        ' transition  energy  lin    cir!       xx        xy', &
     &        '        yy        xz        yz        zz' / )") i0

!     **** loop over all outputted states j ****

      do j = 1, nciout
        if (j /= i0) then
          alp11 = 0.d0
          alp21 = 0.d0
          alp22 = 0.d0
          alp31 = 0.d0
          alp32 = 0.d0
          alp33 = 0.d0

!	  **** loop over all intermediate states k ****
          do k = 1, nconf
            if (k /= i0 .and. k /= j) then
              nj = matind (max(i0, k)) + min(i0, j)
              nk = matind (max(k, j)) + min(k, j)
              d = 1.d0 / (ener(k) - ener(j))
              alp11 = alp11 + dm(nj, 1)*dm(nk, 1) * d
              alp22 = alp22 + dm(nj, 2)*dm(nk, 2) * d
              alp33 = alp33 + dm(nj, 3)*dm(nk, 3) * d
              alp21 = alp21 + (dm(nj, 2)*dm(nk, 1) + dm(nj, 1)*dm(nk, 2)) * d
              alp31 = alp31 + (dm(nj, 3)*dm(nk, 1) + dm(nj, 1)*dm(nk, 3)) * d
              alp32 = alp32 + (dm(nj, 3)*dm(nk, 2) + dm(nj, 2)*dm(nk, 3)) * d
            end if
          end do

!	  **** double count diag terms, calc single - beam intensities ****
          alp11 = alp11 * 2.d0
          alp22 = alp22 * 2.d0
          alp33 = alp33 * 2.d0
          trace = alp11 + alp22 + alp33
          od = alp11**2 + alp22**2 + alp33**2 + 2.d0 * &
     &       (alp21**2 + alp31**2 + alp32**2)
          planar = trace**2 + 2.d0 * od
          circ = -trace**2 + 3.d0 * od
          isym = ixprd (istsym(i0), istsym(j))
          write (iw, '(i4, 1x, a3, f8.3, 2f8.3, 6f10.7)') j, nmrep(isym, nptg), &
            ener(j) - ener(i0), planar, circ, alp11, alp21, alp22, alp31, alp32, alp33
        end if
      end do

      return
      end subroutine twopho

!     *************************************************************************
