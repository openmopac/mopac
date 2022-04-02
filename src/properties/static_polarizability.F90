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

      subroutine static_polarizability
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : keywrd, numat, natoms, nvar, escf, ndep
      use chanel_C, only : iw
      use common_arrays_C, only : geo, coord, nat, labels, na
      use elemts_C, only : elemnt
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k, l
      double precision, dimension(1) :: grad
      double precision, dimension(3,3) :: rotvec
      double precision :: sum, sumx, sumy, sumz, summax
      logical :: let
!-----------------------------------------------
!**********************************************************************
!
!   POLAR SETS UP THE CALCULATION OF THE MOLECULAR ELECTRIC RESPONSE
!   PROPERTIES BY FFHPOL.
!
!**********************************************************************
      let = index(keywrd,' LET') /= 0
      write (iw, 10)
   10 format(' ',20('*'),' FINITE-FIELD POLARIZABILITIES ',20('*'),/,/,&
        '    THE F-F METHOD IS PERFORMED USING BOTH AN ENERGY',/,&
        '    AND DIPOLE MOMENT EXPANSION.  THESE RESULTS ARE',/,&
        '    LISTED BELOW AS "H.o.F." AND "Dipole", RESPECTIVELY.')
      call gmetry (geo, coord)
!
!  Orient the molecule with the moments of inertia.
!  This is done to ensure a unique, reproduceable set of directions.
!  If LET is specified, the input orientation will be used.
!
      if (.not.let) then
        call axis (sumx, sumy, sumz, rotvec)
        write (iw, 20)
   20   format(/,' ROTATION MATRIX FOR ORIENTATION OF MOLECULE:'/)
        do i = 1, 3
          write (iw, 30) (rotvec(i,j),j=1,3)
   30     format(5x,3f12.6)
        end do
!
!  ROTATE ATOMS
!
        do i = 1, numat
          do j = 1, 3
            sum = 0.0D00
            do k = 1, 3
              sum = sum + coord(k,i)*rotvec(k,j)
            end do
            geo(j,i) = sum
          end do
        end do
        coord(:,:numat) = geo(:,:numat)
        write (iw, '(2/10X,''CARTESIAN COORDINATES '',/)')
        write (iw, &
      '(4X,''NO.'',7X,''ATOM'',9X,''X'',                       9X,''Y'',9X,''Z'&
      &',/)')
        l = 0
        do i = 1, numat
          if (nat(i)==99 .or. nat(i)==107) cycle
          l = l + 1
          write (iw, '(I6,8X,A2,4X,3F10.4)') l, elemnt(nat(i)), (coord(j,l),j=1,&
            3)
        end do
      end if
!
!
!  SET UP THE VARIABLES IN XPARAM AND LOC, THESE ARE IN CARTESIAN
!  COORDINATES.
!
      numat = 0
      sumx = 0.D0
      sumy = 0.D0
      sumz = 0.D0
      do i = 1, natoms
        if (labels(i)==99 .or. labels(i)==107) cycle
        numat = numat + 1
        labels(numat) = labels(i)
        sumx = sumx + coord(1,numat)
        sumy = sumy + coord(2,numat)
        sumz = sumz + coord(3,numat)
        geo(:,numat) = coord(:,numat)
      end do
      sumx = sumx/numat
      sumy = sumy/numat
      sumz = sumz/numat
      summax = 0.D0
      do i = 1, numat
        geo(1,i) = geo(1,i) - sumx
        summax = dmax1(abs(geo(1,i)),summax)
        geo(2,i) = geo(2,i) - sumy
        summax = dmax1(abs(geo(2,i)),summax)
        geo(3,i) = geo(3,i) - sumz
        summax = dmax1(abs(geo(3,i)),summax)
      end do
!
      ndep = 0
      natoms = numat
      nvar = 0
      na = 0
      call compfg (geo, .TRUE., escf, .TRUE., grad, .FALSE.)
      write (iw,"(/,/,' ENERGY OF ""REORIENTED"" SYSTEM WITHOUT FIELD:',f15.5,a)") &
      escf, " Kcal/mol"
!
      call ffhpol ()
!
      return
      end subroutine static_polarizability
!
      subroutine ffhpol()
      use molkst_C, only : keywrd, numat, efield,escf
      use common_arrays_C, only : geo, nat
      use parameters_C, only : tore
      use funcon_C, only : fpc_9, eV,a0, fpc_1, fpc_8
      use chanel_C, only : iw
      implicit none
!
      integer :: nbdip, nbcnt, ngcnt, id, ivl, kd, kvl, &
        idm1, jd, i3, counter
      double precision, dimension(3) :: eigs
      double precision, dimension(9) :: vectrs
      double precision, dimension(3) :: dipe4
      double precision, dimension(6) :: apole4
      double precision, dimension(3) :: dipdp
      double precision, dimension(6) :: apoldp
      double precision, dimension(3) :: dip1p, dip1m, dip2p, dip2m
      double precision :: autokc, autodb, autovm, efval, sfe&
        , hnuc, heat1p, grad(1), heat1m, heat2p, heat2m, eterm, dmu, &
        ae, aki, hnucj, hpp, hpm, hmm, hmp, &
        h2pp, h2pm, h2mm, h2mp, aterm, aij, &
        dipe4t, dipe4d, dipdpt, dipdpd, avgpe4, avga3, avgesu, avgpdp, &
        avga3d, avgesd, heat0
      logical :: debug
      character, dimension(3) :: axis
      character :: ch_xyz(3)*1
      double precision, external :: pol_vol
      data ch_xyz/"X","Y","Z"/
!***********************************************************************
!  SUBROUTINE FOR THE FINITE FIELD CALCULATION OF ELECTRIC RESPONSE
!  PROPERTIES (DIPOLE MOMENT, POLARIZABILITY, AND 1ST AND 2ND
!  HYPERPOLARIZABILITY.
!
!  HENRY A. KURTZ, DEPARTMENT OF CHEMISTRY
!                  MEMPHIS STATE UNIVERSITY
!                  MEMPHIS, TN   38152
!
!***********************************************************************
!
!
!     DIPE4 AND DIPDP HOLD THE CALCULATED DIPOLE MOMENTS
!
!     APOLE4 AND APOLDP HOLD THE POLARIZABILITY TENSOR AS
!                                A PACKED ARRAY XX,XY,YY,XZ,YZ,ZZ
!
! Energy: a.u. to kcal/mole
      autokc = fpc_9*eV
! Dipole: a.u. to debye
      autodb = a0*fpc_8*fpc_1*1.d-10
! Electric Field: a.u. to volt/meter
      autovm = eV/a0 ! = 51.42
      nbdip = 1
      nbcnt = 4
      ngcnt = 4
      heat0 = escf
!
      data axis/ 'X', 'Y', 'Z'/
      debug = index(keywrd,'DEBUG') /= 0
!
!  FIELD STRENGTH IN A.U.
!
      efval = 0.002D0
      write (iw, "(/,/,' APPLIED ELECTRIC FIELD MAGNITUDE: ',f15.5,a)") efval*ev/a0, &
      " VOLTS PER ANGSTROM"
      sfe = 1.D00/efval
!.......................................................................
!  CALCULATE THE POLARIZABILITY AND HYPERPOLARIZABILITIES ALONG
!  THE THREE PRINCIPLE AXES.  (THESE AXES DEPEND ON YOUR ARBITRARY
!  ORIENTATION AND MAY NOT BE THE TRUE PRINCIPLE AXES.)
!.......................................................................
      counter = 0
      do id = 1, 3
        if (debug) then
          write (iw, 40) axis(id)
   40     format(/,/,' ****** ',a1,' DIRECTION *****',/)
        end if
!
! ZERO THE FIELD
!
        efield = 0.0D00
        hnuc = sum(efval*geo(id,:numat)*tore(nat(:numat))*autovm)
        hnuc = hnuc*fpc_9
! +E(ID)
        efield(id) = efval
        call compfg (geo, .TRUE., heat1p, .TRUE., grad, .FALSE.)
        counter = counter + 1
        write(iw,'(a,i3,a)')"    Step", counter, " of 36 done (+"//ch_xyz(id)//")"
        call dipind (dip1p)
! -E(ID)
        efield(id) = -efval
        call compfg (geo, .TRUE., heat1m, .TRUE., grad, .FALSE.)
        counter = counter + 1
        write(iw,'(a,i3,a)')"    Step", counter, " of 36 done (-"//ch_xyz(id)//")"
        call dipind (dip1m)
! +2E(ID)
        efield(id) = 2.0D00*efval
        call compfg (geo, .TRUE., heat2p, .TRUE., grad, .FALSE.)
        counter = counter + 1
        write(iw,'(a,i3,a)')"    Step", counter, " of 36 done (+2"//ch_xyz(id)//")"
        call dipind (dip2p)
! -2E(ID)
        efield(id) = -2.0D00*efval
        call compfg (geo, .TRUE., heat2m, .TRUE., grad, .FALSE.)
        counter = counter + 1
        write(iw,'(a,i3,a)')"    Step", counter, " of 36 done (-2"//ch_xyz(id)//")"
        call dipind (dip2m)
        if (debug) then
          write (iw, 70)
   70     format(' ENERGIES AT: ',5x,'F',21x,'2F',/)
          write (iw, 80) heat1p, heat2p, heat1m, heat2m
   80     format('   + ',2(f20.10,3x),/,'   - ',2(f20.10,3x))
        end if
!
! DIPOLE
!
        eterm = (1.0D00/12.D00)*(heat2p - heat2m) - &
          (2.0D00/3.0D00)*(heat1p - heat1m)
        dipe4(id) = eterm*sfe/autokc
!
! ALPHA
!
        ivl = (id*(id + 1))/2
        eterm = 2.5D00*heat0 - (4.D00/3.D00)*(heat1p + heat1m) + &
          (1.D00/12.0D00)*(heat2p + heat2m)
        apole4(ivl) = eterm*sfe*sfe/autokc
!
! DIPOLE CALCULATIONS
!
        dmu = (2.0D00/3.0D00)*(dip1p(id)+dip1m(id)) - &
          (1.D00/6.0D00)*(dip2p(id) + dip2m(id))
        dipdp(id) = dmu/autodb
        ae = (2.0D00/3.0D00)*(dip1p(id)-dip1m(id)) - &
          (1.0D00/12.D00)*(dip2p(id) - dip2m(id))
        apoldp(ivl) = ae*sfe/autodb
        do kd = 1, 3
          if (kd < id) then
            kvl = (id*(id - 1))/2 + kd
            aki = (2.0D00/3.0D00)*(dip1p(kd)-dip1m(kd)) - &
              (1.0D00/12.0D00)*(dip2p(kd)-dip2m(kd))
            apoldp(kvl) = aki*sfe/autodb
          end if
          if (kd == id) cycle
          nbdip = nbdip + 1
        end do
!.......................................................................
!
!  NOW CALCULATE THE OFF AXIS RESULTS.
!
!.......................................................................
        idm1 = id - 1
        do jd = 1, idm1
          hnucj = sum(efval*geo(jd,:numat)*tore(nat(:numat))*eV/a0)
          hnucj = hnucj*fpc_9
          efield = 0.0D00
!
! DIAGONAL FIELDS WITH COMPONENTS EQUAL TO EFVAL
!
          efield(id) = efval
          efield(jd) = efval
          call compfg (geo, .TRUE., hpp, .TRUE., grad, .FALSE.)
          counter = counter + 1
          write(iw,'(a,i3,a)')"    Step", counter, " of 36 done (+"//ch_xyz(id)//"+"//ch_xyz(jd)//")"
          call dipind (dip1p)
          efield(jd) = -efval
          call compfg (geo, .TRUE., hpm, .TRUE., grad, .FALSE.)
          counter = counter + 1
          write(iw,'(a,i3,a)')"    Step", counter, " of 36 done (+"//ch_xyz(id)//"-"//ch_xyz(jd)//")"
          call dipind (dip1p)
          efield(id) = -efval
          call compfg (geo, .TRUE., hmm, .TRUE., grad, .FALSE.)
          counter = counter + 1
          write(iw,'(a,i3,a)')"    Step", counter, " of 36 done (-"//ch_xyz(id)//"-"//ch_xyz(jd)//")"
          call dipind (dip1p)
          efield(jd) = efval
          call compfg (geo, .TRUE., hmp, .TRUE., grad, .FALSE.)
          counter = counter + 1
          write(iw,'(a,i3,a)')"    Step", counter, " of 36 done (-"//ch_xyz(id)//"+"//ch_xyz(jd)//")"
          call dipind (dip1p)
          hpp = hpp + hnuc + hnucj
          hpm = hpm + hnuc - hnucj
          hmm = hmm - hnuc - hnucj
          hmp = hmp - hnuc + hnucj
          if (debug) then
            write (iw, 120)
  120       format(/,'  ',12x,'+,+',15x,'+,-',15x,'-,+',15x,'-,-')
            write (iw, 130) hpp, hpm, hmp, hmm
  130       format('  E ',4f18.6)
          end if
!
!  DIAGONAL FIELDS WITH COMPONENTS EQUAL TO 2*EFVAL
!
          efield(id) = efval*2.D00
          efield(jd) = efval*2.D00
          call compfg (geo, .TRUE., h2pp, .TRUE., grad, .FALSE.)
          counter = counter + 1
          write(iw,'(a,i3,a)')"    Step", counter, " of 36 done (+2"//ch_xyz(id)//"+2"//ch_xyz(jd)//")"
          efield(jd) = -efval*2.D00
          call compfg (geo, .TRUE., h2pm, .TRUE., grad, .FALSE.)
          counter = counter + 1
          write(iw,'(a,i3,a)')"    Step", counter, " of 36 done (+2"//ch_xyz(id)//"-2"//ch_xyz(jd)//")"
          efield(id) = -efval*2.D00
          call compfg (geo, .TRUE., h2mm, .TRUE., grad, .FALSE.)
          counter = counter + 1
          write(iw,'(a,i3,a)')"    Step", counter, " of 36 done (-2"//ch_xyz(id)//"-2"//ch_xyz(jd)//")"
          efield(jd) = efval*2.D00
          call compfg (geo, .TRUE., h2mp, .TRUE., grad, .FALSE.)
          counter = counter + 1
          write(iw,'(a,i3,a)')"    Step", counter, " of 36 done (-2"//ch_xyz(id)//"+2"//ch_xyz(jd)//")"
          h2pp = h2pp + 2.0D00*(hnuc + hnucj)
          h2pm = h2pm + 2.0D00*(hnuc - hnucj)
          h2mm = h2mm - 2.0D00*(hnuc + hnucj)
          h2mp = h2mp - 2.0D00*(hnuc - hnucj)
          if (debug) then
            write (iw, 140) h2pp, h2pm, h2mp, h2mm
  140       format(' 2E ',4f18.6)
          end if
!
          aterm = (1.0D00/48.0D00)*(h2pp - h2pm - h2mp + h2mm) - &
            (1.0D00/3.0D00)*(hpp - hpm - hmp + hmm)
          aij = aterm*sfe*sfe/autokc
          ivl = (id*(id - 1))/2 + jd
          apole4(ivl) = aij
          nbcnt = nbcnt + 1
          nbcnt = nbcnt + 1
          ngcnt = ngcnt + 1
        end do
!
      end do
!-----------------------------------------------------------------------
!  SUMMARIZE THE RESULTS
!-----------------------------------------------------------------------
      write (iw, 170)
  170 format(/,/,' ',30('*'),' DIPOLE ',30('*'),/,/)
      dipe4t = sqrt(dipe4(1)*dipe4(1)+dipe4(2)*dipe4(2)+dipe4(3)*dipe4(3))
      dipe4d = dipe4t*autodb
      dipdpt = sqrt(dipdp(1)*dipdp(1)+dipdp(2)*dipdp(2)+dipdp(3)*dipdp(3))
      dipdpd = dipdpt*autodb
      write (iw, 180)
  180 format(21x,'H.o.F.         Dipole',/)
      write (iw, 190) 'X', dipe4(1), dipdp(1)
      write (iw, 190) 'Y', dipe4(2), dipdp(2)
      write (iw, 190) 'X', dipe4(3), dipdp(3)
  190 format(5x,a1,7x,2f15.6)
      write (iw, 200) dipe4t, dipdpt, dipe4d, dipdpd
  200 format(/,/,' MAGNITUDE:  ',2f15.6,'  (A.U.)',/,' ',12x,2f15.6,'  (DEBYE)')
!
! FIND EIGENVALUES AND EIGENVECTORS OF POLARIZATION MATRIX.
!
      write (iw, 210)
  210 format(/,/,' ',30('*'),' POLARIZABILITY ',20('*'),/,/)
      write (iw, 220)
  220 format(/,' Polarizability Tensor from Heat of Formation:')
      apole4(1) = pol_vol(apole4(1))
      apole4(3) = pol_vol(apole4(3))
      apole4(6) = pol_vol(apole4(6))
      i3 = 3
      call vecprt (apole4, (-i3))
      call rsp (apole4, i3, eigs, vectrs)
      call matout (vectrs, eigs, i3, i3, i3)
      avgpe4 = (eigs(1)+eigs(2)+eigs(3))/3.D0
      avga3 = avgpe4*0.14818D00
      avgesu = avgpe4*0.296352D-24
      write (iw, 230)
  230 format(/,' Polarizability Tensor from Dipole Moment:')
      apoldp(1) = pol_vol(apoldp(1))
      apoldp(3) = pol_vol(apoldp(3))
      apoldp(6) = pol_vol(apoldp(6))
      call vecprt (apoldp, (-i3))
      call rsp (apoldp, i3, eigs, vectrs)
      call matout (vectrs, eigs, i3, i3, i3)
      avgpdp = (eigs(1)+eigs(2)+eigs(3))/3.D0
      avga3d = avgpdp*0.14818D00
      avgesd = avgpdp*0.296352D-24
      write (iw, 240) avgpe4, avgpdp, avga3, avga3d, avgesu, avgesd
  240 format(/,/,' Average Polarizability from:     H.o.F         Dipole',/,' ',24x,2f15.6,&
        '  A.U.',/,' ',24x,2f15.6,'  ANG.**3',/,' ',24x,2(1p,d15.6),'  ESU')

!
      return
      end subroutine ffhpol
      subroutine dipind(dipvec)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use common_arrays_C, only : p, nfirst, nlast, nat, geo
      use molkst_C, only : numat, mozyme
      use funcon_C, only : a0, fpc_8, fpc_1
      use parameters_C, only : ams, tore, dd
!...............................................................
!  MODIFICATION OF DIPOLE SUBROUTINE FOR USE IN THE CALUCLATION
!  OF THE INDUCED DIPOLES FOR POLARIZABILITIES.
!...............................................................
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision , intent(out) :: dipvec(3)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, l, j, ni, ia, k
      double precision, dimension(numat) :: q, q2
      double precision, dimension(3) :: center
      double precision, dimension(3,numat) :: coord
      double precision, dimension(4,3) :: dip
      double precision :: sum, hyfsp
      double precision, save :: wtmol = 0.0d0
      logical, save :: first, chargd
      integer, external :: ijbo
!-----------------------------------------------
!
!***********************************************************************
!     DIPOLE CALCULATES DIPOLE MOMENTS
!
!  ON INPUT P     = DENSITY MATRIX
!           Q     = TOTAL ATOMIC CHARGES, (NUCLEAR + ELECTRONIC)
!           NUMAT = NUMBER OF ATOMS IN MOLECULE
!           NAT   = ATOMIC NUMBERS OF ATOMS
!           NFIRST= START OF ATOM ORBITAL COUNTERS
!           COORD = COORDINATES OF ATOMS
!
!  OUTPUT  DIPOLE = DIPOLE MOMENT
!***********************************************************************
!
!     IN THE ZDO APPROXIMATION, ONLY TWO TERMS ARE RETAINED IN THE
!     CALCULATION OF DIPOLE MOMENTS.
!     1. THE POINT CHARGE TERM (INDEPENDENT OF PARAMETERIZATION).
!     2. THE ONE-CENTER HYBRIDIZATION TERM, WHICH ARISES FROM MATRIX
!     ELEMENTS OF THE FORM <NS/R/NP>. THIS TERM IS A FUNCTION OF
!     THE SLATER EXPONENTS (ZS,ZP) AND IS THUS DEPENDENT ON PARAMETER-
!     IZATION. THE HYBRIDIZATION FACTORS (HYF(I)) USED IN THIS SUB-
!     ROUTINE ARE CALCULATED FROM THE FOLLOWING FORMULAE.
!     FOR SECOND ROW ELEMENTS <2S/R/2P>
!     HYF(I)= 469.56193322*(SQRT(((ZS(I)**5)*(ZP(I)**5)))/
!           ((ZS(I) + ZP(I))**6))
!     FOR THIRD ROW ELEMENTS <3S/R/3P>
!     HYF(I)=2629.107682607*(SQRT(((ZS(I)**7)*(ZP(I)**7)))/
!           ((ZS(I) + ZP(I))**8))
!     FOR FOURTH ROW ELEMENTS AND UP :
!     HYF(I)=2*(2.10716)*DD(I)
!     WHERE DD(I) IS THE CHARGE SEPARATION IN ATOMIC UNITS
!
!
!     REFERENCES:
!     J.A.POPLE & D.L.BEVERIDGE: APPROXIMATE M.O. THEORY
!     S.P.MCGLYNN, ET AL: APPLIED QUANTUM CHEMISTRY
!
      data first/ .TRUE./, chargd /.FALSE./

!
!  SETUP FOR DIPOLE CALCULATION
!
      call chrge (p, q2)
      q(:numat) = tore(nat(:numat)) - q2(:numat)
      call gmetry (geo, coord)
!
      if (first) then
        wtmol = 0.D0
        sum = 0.D0
        do i = 1, numat
          wtmol = wtmol + ams(nat(i))
          sum = sum + q(i)
        end do
        chargd = abs(sum) > 0.5D0
        first = .FALSE.
      end if
      if (chargd) then
!
!   NEED TO RESET ION'S POSITION SO THAT THE CENTER OF MASS IS AT THE
!   ORIGIN.
!
        center = 0.D0
        do i = 1, 3
          do j = 1, numat
            center(i) = center(i) + ams(nat(j))*coord(i,j)
          end do
        end do
        center = center/wtmol
        do i = 1, 3
          coord(i,:numat) = coord(i,:numat) - center(i)
        end do
      end if
      dip = 0.0D00
      do i = 1, numat
        ni = nat(i)
        ia = nfirst(i)
        l = Min(nlast(i) - ia, 3)
        hyfsp = 2.0D0*dd(ni)*a0*fpc_8*fpc_1*1.D-10
        if (mozyme) then
          do j = 1, l
            k = ijbo (i, i) + 1 + (j*(j+1)) / 2
            dip(j, 2) = dip(j, 2) - hyfsp * p(k)
          end do
        else
          do j = 1, l
            k = ((ia + j)*(ia + j - 1))/2 + ia
            dip(j,2) = dip(j,2) - hyfsp*p(k)
          end do
        end if
        dip(:3,1) = dip(:3,1) + 4.803D00*q(i)*coord(:,i)
      end do
      dip(:3,3) = dip(:3,2) + dip(:3,1)
      do j = 1, 3
        dip(4,j) = sqrt(dip(1,j)**2+dip(2,j)**2+dip(3,j)**2)
      end do
      dipvec(1) = -dip(1,3)
      dipvec(2) = -dip(2,3)
      dipvec(3) = -dip(3,3)
      return
!
      end subroutine dipind
