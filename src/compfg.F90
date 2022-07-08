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

      subroutine compfg(xparam, int, escf, fulscf, grad, lgrad)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!
      USE funcon_C, only : fpc_9
!
      USE chanel_C, only : iw
!
      use elemts_C, only : elemnt
!
      USE molmec_C, only : nnhco, nhco, htype
!
      use MOZYME_C, only : partf
!
      use common_arrays_C, only : nat, loc, geo, na, nb, nc, geoa, &
      & coord, xparef, aicorr, tvec, labels, pa, p, pdiag, f, nlast, h, nfirst, c, eigs
!
      USE molkst_C, ONLY: numat, norbs, nelecs, nclose, nopen, fract, natoms, numcal, &
      & ndep, nvar, elect, enuclr, keywrd, moperr, emin, mozyme, method_PM7, lxfac, &
      atheat, id, pressure, method_pm6, density, stress, N_3_present, Si_O_H_present, &
      use_ref_geo, hpress, nsp2_corr, Si_O_H_corr, sum_dihed, method_PM6_D3H4X, method_PM6_D3H4, &
      method_PM6_D3, method_pm7_minus, method_pm6_dh_plus, method_pm7_hh, method_pm8, &
      method_indo, mpack, e_disp
!
      use cosmo_C, only : iseps, useps, noeps, solv_energy
!
      use linear_cosmo, only : ini_linear_cosmo, coscavz
      use reimers_C, only: x, y, z, xz, zcore, beta, gamma, s, betao, ibf, &
      & natm, r, nbf, nbt, nprn, iat, natt, zcorea, betaa, matind, n, &
      & nprin, vnn, dm, ef, dd, ff, cc0, aa, dtmp, nb2, ppg, pg, nsym
      use mopac_interface_flags, only : reset_compfg_L
!
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision , intent(out) :: escf
      logical , intent(in) :: int
      logical, intent(in)  :: fulscf
      logical , intent(in) :: lgrad
      double precision , intent(in) :: xparam(nvar)
      double precision  :: grad(nvar)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icalcn, i, j, k, l
      double precision, dimension(3) :: degree
      double precision :: angle, atheat_store, sum, store_e_disp
      double precision, external ::  nsp2_correction, Si_O_H_correction
      double precision, external :: helecz
      logical :: debug, print, large, usedci, force, times, aider, &
        dh, l_locate_ts
      double precision, external :: xfac_value, reada
      double precision, external :: ddot, dot, volume
      character :: tmpkey*241

      save debug, print, large, usedci, force, times, aider, degree, &
        icalcn, dh, l_locate_ts
!***********************************************************************
!
!   COMPFG CALCULATES (A) THE HEAT OF FORMATION OF THE SYSTEM, AND
!                     (B) THE GRADIENTS, IF LGRAD IS .TRUE.
!
!   ON INPUT  XPARAM = ARRAY OF PARAMETERS TO BE USED IN INTERNAL COORDS
!             LGRAD  = .TRUE. IF GRADIENTS ARE NEEDED, .FALSE. OTHERWISE
!             INT    = .TRUE. IF HEAT OF FORMATION IS TO BE CALCULATED
!             FULSCF = .TRUE. IF FULL SCF TO BE DONE, .FALSE. OTHERWISE.
!
!   ON OUTPUT ESCF  = HEAT OF FORMATION.
!             GRAD   = ARRAY OF GRADIENTS, IF LGRAD = .TRUE.
!
!***********************************************************************
      data icalcn/ 0/
      if (lxfac) then
        if (index(keywrd," POP") /= 0) then
          pa = 0.d0
          i = index(keywrd," POP") + 4
          pa(1)  = reada(keywrd,i+1)    !  "s"-population
          pa(3)  = reada(keywrd,i+3)/3  !  "p"-population
          pa(6)  = pa(3)
          pa(10) = pa(3)
          pa(15) = reada(keywrd,i+5)/5  !  "d"-population
          pa(21) = pa(15)
          pa(28) = pa(15)
          pa(36) = pa(15)
          pa(45) = pa(15)
          p = pa
          pa = pa*0.5d0
          pdiag(1)   = p(1)
          pdiag(2:4) = p(3)
          pdiag(5:9) = p(15)
        else
          escf = xfac_value()
          return
        end if
      end if
      if (icalcn /= numcal .or. reset_compfg_L) then
         reset_compfg_L = .false.
        icalcn = numcal
        hpress = 0.d0
        nsp2_corr = 0.d0
        Si_O_H_corr = 0.d0
        sum_dihed = 0.d0
        if (iseps) then
          iseps = .true.
          noeps = .true.
          call cosini(.true.)
          if (moperr) return
          mozyme = (index(keywrd," MOZ") + index(keywrd," LOCATE-TS") + index(keywrd," RAPID") /= 0)
          if (mozyme .and. numat == 1) then
            call mopend("MOZYME cannot be used for systems composed of only one atom!")
            return
          end if
          if (mozyme) call ini_linear_cosmo
        else
          iseps = .false.
          noeps = .false.
        end if
        aider = index(keywrd,'AIDER') /= 0
        times = index(keywrd,'TIMES') /= 0
        usedci = nclose/=nopen .and. Abs(fract - 2.d0) > 1.d-20 .and. &
          fract > 1.d-20 .or. index(keywrd,'C.I.')/=0
        force = index(keywrd,'FORCE') /= 0
        large = index(keywrd,'LARGE') /= 0
        print = index(keywrd,'COMPFG') /= 0
        l_locate_ts = index(keywrd,'LOCATE-TS') /= 0
        debug = index(keywrd,'DEBUG')/=0 .and. print
        dh = (method_pm6_d3h4x .or. method_pm6_d3h4    .or. &
              method_pm6_d3    .or. method_pm6_dh_plus .or. &
              method_pm7_hh    .or. method_pm7_minus   .or. &
              method_PM7 .or. method_pm8)
        emin = 0.D0
        xparef(:nvar) = xparam(:nvar)
      end if
!
! SET UP COORDINATES FOR CURRENT CALCULATION
!
!       PLACE THE NEW VALUES OF THE VARIABLES IN THE ARRAY GEO.
!       MAKE CHANGES IN THE GEOMETRY.
      do i = 1, nvar
        k = loc(1,i)
        l = loc(2,i)
        geo(l,k) = xparam(i)
      end do
!      IMPOSE THE SYMMETRY CONDITIONS + COMPUTE THE DEPENDENT-PARAMETERS
      if (ndep /= 0) call symtry
!      NOW COMPUTE THE ATOMIC COORDINATES.
      if (debug) then
        if (large) then
          k = natoms
        else
          k = min(5,natoms)
        end if
        write(iw,"(a)")" COORDINATES IN ARRAY 'GEO'"
        degree(1) = 1.d0
        do i = 1, k
          if (na(i) > 0) then
            degree(2:3) = 57.29577951308232D0
          else
            degree(2:3) = 1.d0
          end if
          write (iw, "(i4,3x,a2,3x,3F14.5,3i5)") i,elemnt(labels(i)), &
          (geo(j, i)*degree(j), j=1, 3), na(i), nb(i), nc(i)
        end do
      end if
      call gmetry (geo, coord)
      if (moperr) return
      if (debug) then
        if (large) then
          k = numat
        else
          k = min(5,numat)
        end if
        write (iw, '('' CARTESIAN COORDINATES'',/10000(/,i4,3x,a2,3x,3F16.9))') &
          (i,elemnt(nat(i)),(coord(j,i),j=1,3),i=1,k)
      end if
      if (iseps) then
      ! The following routine constructs the dielectric screening surface
      if (mozyme) then
          call coscavz(coord, nat)
        else
          call coscav
          call mkbmat
        end if
        if (moperr) return
        if (noeps) useps = .false.
      end if
      if (index(keywrd,' HCORE') /= 0) call prtpar
      if (times) call timer ('BEFORE HCORE')
      if (mozyme) then
        if (iseps) useps = .true.
        if (l_locate_ts .or. int) call hcore_for_MOZYME ()
        if (moperr) return
      else if (method_indo) then
! Set up Reimers data
        if (allocated(x))     deallocate(x)
        if (allocated(y))     deallocate(y)
        if (allocated(z))     deallocate(z)
        if (allocated(xz))    deallocate(xz)
        if (allocated(zcore)) deallocate(zcore)
        if (allocated(beta))  deallocate(beta)
        if (allocated(gamma)) deallocate(gamma)
        if (allocated(s))     deallocate(s)
        if (allocated(betao)) deallocate(betao)
        if (allocated(ibf))   deallocate(ibf)
        if (allocated(natm))  deallocate(natm)
        if (allocated(r))     deallocate(r)
        if (allocated(nbf))   deallocate(nbf)
        if (allocated(nbt))   deallocate(nbt)
        if (allocated(nprn))  deallocate(nprn)
        if (allocated(iat))   deallocate(iat)
        if (allocated(natt))  deallocate(natt)
        if (allocated(dd))    deallocate(dd)
        if (allocated(cc0))   deallocate(cc0)
        if (allocated(aa))    deallocate(aa)
        if (allocated(ff))    deallocate(ff)
        if (allocated(dtmp))  deallocate(dtmp)
        if (allocated(dd))    deallocate(dd)
        if (allocated(ppg))   deallocate(ppg)
        if (allocated(pg))    deallocate(pg)
        if (allocated(nsym))  deallocate(nsym)
        allocate(x(numat))
        allocate(y(numat))
        allocate(z(numat))
        allocate(xz(numat,3))
        allocate(zcore(numat))
        allocate(beta(mpack))
        allocate(gamma(norbs,norbs))
        gamma = 0.d0
        allocate(s(mpack))
        allocate(betao(norbs))
        allocate(ibf(numat))
        allocate(natm(numat))
        allocate(r(numat,numat))
        allocate(nbf(numat))
        allocate(nbt(norbs))
        allocate(nprn(norbs))
        allocate(iat(norbs))
        allocate(natt(norbs))


        matind(1) = 0
        do i=2,50000
          matind(i) = matind(i-1) + i-1
        end do
        n = norbs
     !   na = numat

        do i=1,numat
          x(i) = coord(1,i)
          y(i) = coord(2,i)
          z(i) = coord(3,i)
          xz(i,1) = coord(1,i)
          xz(i,2) = coord(2,i)
          xz(i,3) = coord(3,i)
          zcore(i) = zcorea(labels(i))
          ibf(i) = nfirst(i)
          natm(i) = labels(i)
          natt(i) = labels(i)
          nbf(i) = nlast(i) - nfirst(i) + 1
          do j=nfirst(i),nlast(i)
            nbt(j) = j - nfirst(i)
            iat(j) = i
            if (j-nfirst(i).eq.0) then
              betao(j) = betaa(1,labels(i))
              nprn(j) = nprin(labels(i))
            else if (j-nfirst(i).le.3) then
              betao(j) = betaa(2,labels(i))
              nprn(j) = nprin(labels(i))
            else
              betao(j) = betaa(3,labels(i))
              nprn(j) = nprin(labels(i))-1
            end if
          end do
        end do
! Call Reimers INDO routine
        call replsn ()

        enuclr = vnn
        call ovlap  (s,x,y,z)
        call beta1  (s,betao,beta)
! Add electric field correction - borrow code from hcore
        tmpkey = trim(keywrd)
        i = index(tmpkey,' FIELD(') + index(tmpkey,' FIELD=(')
        if (i /= 0) then
!   ERASE ALL TEXT FROM TMPKEY EXCEPT FIELD DATA
          tmpkey(:i) = ' '
          tmpkey(index(tmpkey,')'):) = ' '
!   READ IN THE EFFECTIVE FIELD IN X,Y,Z COORDINATES
          ef(1) = reada(tmpkey,i)
          i = index(tmpkey,',')
          if (i /= 0) then
            tmpkey(i:i) = ' '
            ef(2) = reada(tmpkey,i)
            i = index(tmpkey,',')
            if (i /= 0) then
              tmpkey(i:i) = ' '
              ef(3) = reada(tmpkey,i)
            end if
          end if
          write (iw,'(/10X,''THE ELECTRIC FIELD IS'',3F10.5,&
             &''VOLTS/ANGSTROM'',/)') ef
! Reimers code uses E field in V/A
          do i = 1,3
            ef(i) = -ef(i)
          end do
! Back to Reimers-specific stuff
          nb2 = mpack
          if(.not. allocated(dm))     allocate(dm(nb2,3))
          call dipol (x,y,z,dm)
          call efmods (beta,zcore,dm)
        end if
! Put Reimers data back into h
        do i= 1,mpack
          h(i) = beta(i)
        end do
! Solvent correction
        if (useps) then
          call addnuc ()
          call addhcr ()
!     Put h back into beta
          do i= 1,mpack
            beta(i) = h(i)
          end do
        end if
      else
        if (int) call hcore ()
        if (moperr) return
      end if
       atheat_store = atheat
      if (times) call timer ('AFTER  HCORE')
!
! COMPUTE THE HEAT OF FORMATION.
!
      if (norbs > 0 .and. nelecs > 0) then
!
!  Put any ad-hoc corrections to the HoF here, so they will show up if PL
!  is used
!
        hpress = 0.d0
        if (Abs (pressure) > 1.d-4) then
          if (id == 1) then
            hpress = -pressure * Sqrt (dot(tvec(1, 1), tvec(1, 1), 3))
          else if (id == 3) then
            hpress = -pressure * volume (tvec, 3)
          end if
          atheat = atheat + hpress
        end if
        if (useps .and. .not. mozyme) atheat = atheat + solv_energy * fpc_9
!
!  Add in any molecular-mechanics type corrections here
!
        if (method_pm6 .and. N_3_present) then
          nsp2_corr = nsp2_correction()
          atheat = atheat + nsp2_corr
        end if
        if (method_pm7 .and. Si_O_H_present) then
          Si_o_H_corr = Si_O_H_correction()
          atheat = atheat + Si_O_H_corr
        end if
        call setup_nhco(i)
        sum_dihed = 0.d0
        do i = 1, nnhco
          call dihed (coord, nhco(1,i), nhco(2,i), nhco(3,i), nhco(4,i), angle)
          sum_dihed = sum_dihed + htype*sin(angle)**2
        end do
        atheat = atheat + sum_dihed
        stress = 0.d0
        if(use_ref_geo) then
          do i = 1, numat
            do j = 1,3
              stress = stress + (geo(j,i) - geoa(j,i))**2
            end do
          end do
        end if
        if (dh) then
          call post_scf_corrections(sum, .false.)
          if (moperr) return
          atheat =  sum + atheat
        end if
        if (times) call timer ('BEFORE ITER')
        if (int) then
          if (mozyme) then
            call iter_for_MOZYME (elect)
          else
            call iter (elect, fulscf, .TRUE.)
          end if
          if (moperr) return
          if (noeps) then
            noeps = .false.
            useps = .true.
            if ( .not. mozyme) then
              if (.not. method_indo) call hcore ()
              if (method_indo) then
                call addnuc ()
                call addhcr ()
                do i= 1,mpack
                  beta(i) = h(i)
                end do
              end if
              call iter (elect, fulscf, .TRUE.)
            end if
          end if
        else
          if (mozyme) then
             call buildf (f, partf, 0)
             elect = helecz()
           end if
        end if
        stress = stress*density
        atheat = atheat + stress
        if (moperr) return
        if (times) call timer ('AFTER  ITER')
      else
        elect = 0.D0
      end if
      escf = (elect + enuclr)*fpc_9 + atheat
      if (useps .and. mozyme) then
            escf = escf + solv_energy * fpc_9
      end if
      if (.not. dh) then
        call post_scf_corrections(sum, .false.)
        if (moperr) return
        escf =  sum + escf
      end if
      atheat = atheat_store
      if (escf < emin .or. emin == 0.D0) emin = escf
      if (method_indo) then
        call output (c,eigs)
      end if
!
! FIND DERIVATIVES IF DESIRED
!
      if (lgrad) then
        store_e_disp = e_disp
        if (times) call timer ('Before DERIV')
        if (nelecs > 0) call deriv (geo, grad)
        if (moperr) return
        if (times) call timer ('AFTER  DERIV')
        e_disp = store_e_disp
      end if
      if (aider) then
!
!  ADD IN AB INITIO CORRECTION
!
        escf = escf + ddot(nvar,xparam(:nvar)-xparef(:nvar),1,aicorr(:nvar),1)
      end if
      if (int .and. print) write (iw, '(/1X,'' HEAT OF FORMATION'',G30.17)') &
        escf
      if (print .and. lgrad) then
        write (iw, '('' PARAMETERS     '',8F8.3,(/10F8.3))') (xparam(i),i=1,nvar)
        write (iw, '('' GRADIENT       '',8F8.3,(/10F8.3))') (grad(i),i=1,nvar)
      end if
!
! REFORM DENSITY MATRIX, IF A C.I. DONE AND EITHER THE LAST SCF OR A
! FORCE CALCULATION
!
      if (usedci .and. force .and. .not. method_indo) call mecip ()
      return
      end subroutine compfg
!
!
!
      double precision function xfac_value()
      use parameters_C, only :  tore, guess1, guess2, guess3, alpb, xfac, po, pocord, &
        uss, upp, udd
      use common_arrays_C, only : nat, coord, pa, p, h, f, w
      use funcon_C, only : a0, ev, fpc_9
      USE molkst_C, ONLY: keywrd
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ni, nj, i
      double precision :: enuc, abond, fff, scale, ax, r, &
        gab, enuclr, point, const
      double precision, external :: reada
        if (index(keywrd," POP") /= 0) then
          p = 0.d0
          i = index(keywrd," POP") + 4
          p(1)  = reada(keywrd,i+1)  !  "s"-population
          p(3)  = reada(keywrd,i+3)/3  !  "p"-population
          p(6)  = p(3)
          p(10) = p(3)
          p(15) = reada(keywrd,i+5)/5  !  "d"-population
          p(21) = p(15)
          p(28) = p(15)
          p(36) = p(15)
          p(45) = p(15)
          pa = p*0.5d0
          h = 0.d0
          h(1) = uss(nat(1))
          h(3) = upp(nat(1))
          h(6) = h(3)
          h(10) = h(3)
          h(15) = udd(nat(1))
          h(21) = h(15)
          h(28) = h(15)
          h(36) = h(15)
          h(45) = h(15)
          f = h
          call fock1(f, p, pa, 45, w, i, 1, 9, 45)
          xfac_value = 0.0d0 ! dummy return value for unaccessed branch
          return
        end if
        r = coord(1,2)
        ni = nat(1)
        nj = nat(2)
        if (pocord(ni) > 1.D-5) po(9,ni) = pocord(ni)
        if (pocord(nj) > 1.D-5) po(9,nj) = pocord(nj)
        gab = eV/sqrt((r/a0)**2 + (po(9,ni) + po(9,nj))**2)
        call to_point(r, point, const)
        gab = gab*const + (1.d0 - const)*point
        enuc = tore(ni)*tore(nj)*gab
        abond = alpb(ni,nj)
        if (abond  > 1.d-3) then
          fff = xfac(ni,nj)
          scale = 2.d0 * fff * Exp (-abond*(r + 0.0003*r**6))
          enuclr = enuc * scale
          scale = 0.d0
          ax = guess2(ni,1)*(r - guess3(ni,1))**2
          if (ax < 25.D0) scale = scale + tore(ni)*tore(nj)/r*guess1(ni,1)*exp((-ax))
          ax = guess2(nj,1)*(r - guess3(nj,1))**2
          if (ax < 25.D0) scale = scale + tore(ni)*tore(nj)/r*guess1(nj,1)*exp((-ax))
          enuclr = enuclr + scale
          ax = r/(ni**0.3333d0 + nj**0.3333d0)
          if (ax < 3.d0) then
            scale = 1.d-8/ax**12
            enuclr = enuclr + min(scale, 1.d5)
          end if
        else
          scale = 10.d0*exp((-2.18d0*r)) ! This is a generic core-core term.
          enuclr = abs(scale*enuc)
      end if
        xfac_value = enuclr*fpc_9
      end function xfac_value
