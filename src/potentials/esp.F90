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

      subroutine esp()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : keywrd, moperr
      use chanel_C, only : iw
      use esp_C, only : scale, den, nesp
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: n, i
      double precision :: scincr, time1
      double precision, external :: reada, seconds
!-----------------------------------------------
      time1 = seconds (1)
      if (index(keywrd, " CUBE") /= 0 .or. index(keywrd, " ESPGRID") /= 0) then
        call new_esp()
        return
      end if
!
!     SET STANDARD PARAMETERS FOR THE SURFACE GENERATION
!
      if (index(keywrd,'SCALE=') /= 0) then
        scale = reada(keywrd,index(keywrd,'SCALE='))
      else
        scale = 1.4D0
      end if
!
      if (index(keywrd,'DEN=') /= 0) then
        den = reada(keywrd,index(keywrd,'DEN='))
      else
        den = 1.0D0
      end if
!
      if (index(keywrd,'SCINCR=') /= 0) then
        scincr = reada(keywrd,index(keywrd,'SCINCR='))
      else
        scincr = 0.20D0
      end if
!
      if (index(keywrd,'NSURF=') /= 0) then
        n = nint(reada(keywrd,index(keywrd,'NSURF=')))
      else
        n = 4
      end if
!
      time1 = seconds(1)
!
!     NOW CALCULATE THE SURFACE POINTS
!
      call setup_esp(1)
      nesp = 0
      if (index(keywrd,'WILLIAMS') /= 0) then
        call pdgrid ()
        if (moperr) return
      else
        do i = 1, n
          call surfac ()
          if (moperr) return
          scale = scale + scincr
        end do
      end if
!
!     NEXT CALCULATE THE ESP AT THE POINTS CALCULATED BY SURFAC
!
      call potcal ()
!
!     END OF CALCULATION
!
      time1 = seconds(1) - time1
      write (iw, 20) 'TIME TO CALCULATE ESP:', time1, ' SECONDS'
   20 format(/,9x,a,f8.2,a)
      return
      end subroutine esp
       subroutine pdgrid()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : numat, natoms
      use common_arrays_C, only : geo, coord, labels, nat
      use esp_C, only : nesp, potpt, co
      use chanel_C, only : iw
!
!     ROUTINE TO CALCULATE WILLIAMS SURFACE
!
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icntr, i, j, ia, npnt, l, jz
      double precision, dimension(53) :: vderw
      double precision, dimension(100) :: dist
      double precision, dimension(3) :: xmin, xmax
      double precision :: shell, grid, closer, vdmax, xstart, ystart, zstart, zgrid&
        , ygrid, xgrid
!-----------------------------------------------
!
      data vderw/ 53*0.0D0/
      vderw(1) = 2.4D0
      vderw(5) = 3.0D0
      vderw(6) = 2.9D0
      vderw(7) = 2.7D0
      vderw(8) = 2.6D0
      vderw(9) = 2.55D0
      vderw(15) = 3.1D0
      vderw(16) = 3.05D0
      vderw(17) = 3.0D0
      vderw(35) = 3.15D0
      vderw(53) = 3.35D0
      shell = 1.2D0
      nesp = 0
      grid = 0.8D0
      closer = 0.D0
!     CHECK IF VDERW IS DEFINED FOR ALL ATOMS
!
!     CONVERT INTERNAL TO CARTESIAN COORDINATES
!
      call gmetry (geo, coord)
!
!     STRIP COORDINATES AND ATOM LABEL FOR DUMMIES (I.E. 99)
!
      icntr = 0
      do i = 1, natoms
        co(:,i) = coord(:,i)
        if (labels(i) == 99) cycle
        icntr = icntr + 1
        nat(icntr) = labels(i)
      end do
      numat = icntr
!
      do i = 1, numat
        j = nat(i)
        if (vderw(j) == 0.0D0) go to 40
      end do
      go to 50
   40 continue
      write (iw, *) 'VAN DER WAALS'' RADIUS NOT DEFINED FOR ATOM', i
      write (iw, *) 'IN WILLIAMS SURFACE ROUTINE PDGRID!'
      call mopend (&
      'VAN DER WAALS'' RADIUS NOT DEFINED FOR ATOM IN WILLIAMS SURFACE ROUTINE &
      &PDGRID!')
      return
!     NOW CREATE LIMITS FOR A BOX
   50 continue
      xmin = 100000.0D0
      xmax = -100000.0D0
      do ia = 1, numat
        where (co(:,ia) - xmin < 0.D0)
          xmin = co(:,ia)
        end where
        where (co(:,ia) - xmax > 0.D0)
          xmax = co(:,ia)
        end where
      end do
!     ADD (OR SUBTRACT) THE MAXIMUM VDERW PLUS SHELL
      vdmax = 0.0D0
      vdmax = dmax1(maxval(vderw),vdmax)
      xmin = xmin - vdmax - shell
      xmax = xmax + vdmax + shell
! STEP GRID BACK FROM ZERO TO FIND STARTING POINTS
      xstart = 0.0D0
      xstart = xstart - grid
      do while(xstart > xmin(1))
        xstart = xstart - grid
      end do
      ystart = 0.0D0
      ystart = ystart - grid
      do while(ystart > xmin(2))
        ystart = ystart - grid
      end do
      zstart = 0.0D0
      zstart = zstart - grid
      do while(zstart > xmin(3))
        zstart = zstart - grid
      end do
      npnt = 0
      zgrid = zstart
  160 continue
      ygrid = ystart
  170 continue
      xgrid = xstart
  180 continue
      do l = 1, numat
        jz = nat(l)
        dist(l) = sqrt((co(1,l)-xgrid)**2+(co(2,l)-ygrid)**2+(co(3,l)-zgrid)**2)
!     REJECT GRID POINT IF ANY ATOM IS TOO CLOSE
        if (dist(l) < vderw(jz) - closer) go to 220
      end do
! BUT AT LEAST ONE ATOM MUST BE CLOSE ENOUGH
      do l = 1, numat
        jz = nat(l)
        if (dist(l) > vderw(jz) + shell) cycle
        go to 210
      end do
      go to 220
  210 continue
      npnt = npnt + 1
      nesp = nesp + 1
      potpt(1,nesp) = xgrid
      potpt(2,nesp) = ygrid
      potpt(3,nesp) = zgrid
  220 continue
      xgrid = xgrid + grid
      if (xgrid <= xmax(1)) go to 180
      ygrid = ygrid + grid
      if (ygrid <= xmax(2)) go to 170
      zgrid = zgrid + grid
      if (zgrid <= xmax(3)) go to 160
      return
      end subroutine pdgrid
      subroutine surfac()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use common_arrays_C, only :geo, coord, labels, nat
      use molkst_C, only : numat, natoms
      use esp_C, only : nesp, rad, dens, scale, potpt, co, den
      use chanel_C, only : iw
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(numat) :: ias
      integer :: icntr, i, ipoint, iatom, nnbr, jatom, ncon
      double precision, dimension(100) :: vander
      double precision, dimension(3,1000) :: con
      double precision, dimension(3,200) :: cnbr
      double precision, dimension(200) :: rnbr
      double precision, dimension(3) :: ci, temp0
      double precision, dimension(3,2) :: cw
      double precision :: pi, rw, ri, d2
      logical :: si
      double precision, external :: distance
      logical, external :: collid
!-----------------------------------------------
!***********************************************************************
!
!      THIS SUBROUTINE CALCULATES THE MOLECULAR SURFACE OF A MOLECULE
!      GIVEN THE COORDINATES OF ITS ATOMS.  VAN DER WAALS' RADII FOR
!      THE ATOMS AND THE PROBE RADIUS MUST ALSO BE SPECIFIED.
!
!      ON INPUT    SCALE = INITIAL VAN DER WAALS' SCALE FACTOR
!                  DENS  = DENSITY OF POINTS PER UNIT AREA
!
!      THIS SUBROUTINE WAS LIFTED FROM MICHAEL CONNOLLY'S SURFACE
!      PROGRAM FOR UCSF GRAPHICS SYSTEM BY U.CHANDRA SINGH AND
!      P.A.KOLLMAN AND MODIFIED FOR USE IN QUEST. K.M.MERZ
!      ADAPTED AND CLEANED UP THIS PROGRAM FOR USE IN AMPAC/MOPAC
!      IN FEB. 1989 AT UCSF.
!
!***********************************************************************
!
!
!     NEIGHBOR ARRAYS
!
!     THIS SAME DIMENSION FOR THE MAXIMUM NUMBER OF NEIGHBORS
!     IS USED TO DIMENSION ARRAYS IN THE LOGICAL FUNCTION COLLID
!
!#      DIMENSION INBR(200)
!
!     ARRAYS FOR ALL ATOMS
!
!     IATOM, JATOM AND KATOM COORDINATES
!
!
!     GEOMETRIC CONSTRUCTION VECTORS
!
!
!     LOGICAL VARIABLES
!
!
!     LOGICAL FUNCTIONS
!
!
!     DATA FOR VANDER VALL RADII
!
      data vander/ 1.20D0, 1.20D0, 1.37D0, 1.45D0, 1.45D0, 1.50D0, 1.50D0, &
        1.40D0, 1.35D0, 1.30D0, 1.57D0, 1.36D0, 1.24D0, 1.17D0, 1.80D0, 1.75D0&
        , 1.70D0, 17*0.0D0, 2.3D0, 65*0.0D0/
!
      pi = 4.D0*atan(1.D0)
!     INSERT VAN DER WAAL RADII FOR ZINC
      vander(30) = 1.00D0
!
!     CONVERT INTERNAL TO CARTESIAN COORDINATES
!
      call gmetry (geo, coord)
!
!     STRIP COORDINATES AND ATOM LABEL FOR DUMMIES (I.E. 99)
!
      icntr = 0
      co(:,:numat) = coord(:,:numat)
      do i = 1, natoms
        if (labels(i) == 99) cycle
        icntr = icntr + 1
        nat(icntr) = labels(i)
      end do
!
!     ONLY VAN DER WAALS' TYPE SURFACE IS GENERATED
!
!#      NESP=0
      rw = 0.0D0
      numat = icntr
      dens = den
      do i = 1, numat
        ipoint = nat(i)
        rad(i) = vander(ipoint)*scale
        if (rad(i) < 0.01D0) write (iw, &
      '(T2,''VAN DER WAALS'''' RADIUS FOR ATOM '',I3,'' IS ZERO, SUPPLY A VALUE IN SUBROUTINE SURFAC)''                )')
        ias(i) = 2
      end do
!
!     BIG LOOP FOR EACH ATOM
!
      do iatom = 1, numat
        if (ias(iatom) == 0) cycle
!
!     TRANSFER VALUES FROM LARGE ARRAYS TO IATOM VARIABLES
!
        ri = rad(iatom)
        si = ias(iatom) == 2
        ci = co(:,iatom)
!
!     GATHER THE NEIGHBORING ATOMS OF IATOM
!
        nnbr = 0
        do jatom = 1, numat
          if (iatom==jatom .or. ias(jatom)==0) cycle
          d2 = distance(iatom,jatom)
          if (d2 >= (2*rw + ri + rad(jatom))**2) cycle
!
!     WE HAVE A NEW NEIGHBOR
!     TRANSFER ATOM COORDINATES, RADIUS AND SURFACE REQUEST NUMBER
!
          nnbr = nnbr + 1
          if (nnbr > 200) then
            write (iw, '(''ERROR'',2X,''TOO MANY NEIGHBORS:'',I5)') nnbr
            call mopend ('ERROR.  TOO MANY NEIGHBORS')
            return
          end if
!#            INBR(NNBR) = JATOM
          cnbr(:,nnbr) = co(:,jatom)
          rnbr(nnbr) = rad(jatom)
        end do
!
!     CONTACT SURFACE
!
        if (.not.si) cycle
        ncon = int((4*pi*ri**2)*den)
        ncon = min0(1000,ncon)
!
!     THIS CALL MAY DECREASE NCON SOMEWHAT
!
        if (ncon == 0) then
          write (iw, '(T2,''VECTOR LENGTH OF ZERO IN SURFAC'')')
          call mopend ('VECTOR LENGTH OF ZERO IN SURFAC')
          return
        end if
        call genun (con, ncon)
!#         AREA = (4 * PI * RI ** 2) / NCON
!
!     CONTACT PROBE PLACEMENT LOOP
!
        do i = 1, ncon
          cw(:,1) = ci + (ri + rw)*con(:,i)
!
!     CHECK FOR COLLISION WITH NEIGHBORING ATOMS
!
          if (collid(cw(1,1),rw,cnbr,rnbr,nnbr,1)) cycle
          temp0 = ci + ri*con(:,i)
!
!     STORE POINT IN POTPT AND INCREMENT NESP
!
          nesp = nesp + 1
          potpt(1,nesp) = temp0(1)
          potpt(2,nesp) = temp0(2)
          potpt(3,nesp) = temp0(3)
        end do
      end do
      return
      end subroutine surfac

     subroutine potcal()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : keywrd, numat, ux, uy, uz, method_mndo
      USE common_arrays_C, only : q, nat
      use funcon_C, only : a0
      use elemts_C, only : elemnt
      use esp_C, only : cf, cequiv, co, qsc, nesp, &
      & rms, rrms, dx, dy, dz, iz, idip
      use chanel_C, only : iw
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, ieq
      double precision :: dipx, dipy, dipz, slope, dip
      double precision, external :: reada
!-----------------------------------------------
!***********************************************************************
!
!     THIS SUBROUTINE CALCULATES THE TOTAL ELECTROSTATIC POTENTIAL
!     THE NUCLEAR CONTRIBUTION IS EVALUATED BY NUCPOT
!     THE ELECTRONIC CONTRIBUTION IS EVALUATED BY ELESP
!     ESPFIT FITS THE QUANTUM POTENTIAL TO A CLASSICAL POINT CHARGE
!     MODEL.
!     THIS SUBROUTINE WAS WRITTEN BY B.H.BESLER AND K.M.MERZ IN FEB.
!     1989 AT UCSF
!
!***********************************************************************
      data dipx, dipy, dipz/ 3*0.D0/
!
!
      call elesn ()
      rms = 0.D0
      rrms = 0.D0
      slope = 0.d0
!
!     NOW FIT THE ELECTROSTATIC POTENTIAL
!
      write (iw, '(2/12X,''ELECTROSTATIC POTENTIAL CHARGES'',/)')
      iz = 0
      if (index(keywrd,'CHARGE=') /= 0) iz = nint(reada(keywrd,index(keywrd,&
        'CHARGE=')))
!
!     DIPOLAR CONSTRAINTS IF DESIRED
!
      if (index(keywrd,'DIPOLE') /= 0) then
        idip = 1
        if (iz /= 0) then
          idip = 0
          write (iw, '(/12X,''  DIPOLE CONSTRAINTS NOT USED'')')
          write (iw, '(12X,''        CHARGED MOLECULE'',/)')
        end if
      else
        idip = 0
      end if
      if (idip == 1) write (iw, '(/12X,''DIPOLE CONSTRAINTS WILL BE USED'',/)')
!
!     GET X,Y,Z DIPOLE COMPONENTS IF DESIRED
!
      if (index(keywrd,'DIPX=') /= 0) then
        dx = reada(keywrd,index(keywrd,'DIPX='))
      else
        dx = ux
      end if
      if (index(keywrd,'DIPY=') /= 0) then
        dy = reada(keywrd,index(keywrd,'DIPY='))
      else
        dy = uy
      end if
      if (index(keywrd,'DIPZ=') /= 0) then
        dz = reada(keywrd,index(keywrd,'DIPZ='))
      else
        dz = uz
      end if
      call espfit ()
!
!     WRITE OUT OUR RESULTS TO CHANNEL 6
!     THE CHARGES ARE SCALED TO REPRODUCE 6-31G* CHARGES FOR MNDO ONLY
!     AM1 CHARGES ARE NOT SCALED DUE TO THE LOW COORELATION
!     COEFFICIENT. SEE BESLER,MERZ,KOLLMAN IN J. COMPUT. CHEM.
!     (IN PRESS)
!
      if (.not. method_mndo) then
        write (iw, '(15X,''ATOM NO.    TYPE    CHARGE'')')
        do i = 1, numat
          write (iw, '(17X,I4,9X,A2,1X,F10.4)') i, elemnt(nat(i)), q(i)
        end do
        call to_screen("To_file: Esp")
      else
!
!     MNDO CALCULATION-SCALE THE CHARGES. TEST FOR SLOPE KEYWORD
!
        if (index(keywrd,'SLOPE=') /= 0) then
          slope = reada(keywrd,index(keywrd,'SLOPE='))
        else
          slope = 1.422D0
        end if
        qsc(:numat) = slope*q(:numat)
        write (iw, '(7X,''ATOM NO.    TYPE    CHARGE   SCALED CHARGE'')')
        do i = 1, numat
          write (iw, '(9X,I4,9X,A2,1X,F10.4,2X,F10.4)') i, elemnt(nat(i)), q(i), qsc(i)
        end do
!     Temporarily scale charges in q for AUX printing
        q(:numat) = qsc(:numat)
        call to_screen("To_file: Esp")
        q(:numat) = qsc(:numat)/slope
      end if
      write (iw, '(/12X,A,4X,I6)') 'THE NUMBER OF POINTS IS:', nesp
      write (iw, '(12X,A,4X,F9.4)') 'THE RMS DEVIATION IS:', rms
      write (iw, '(12X,A,3X,F9.4)') 'THE RRMS DEVIATION IS:', rrms
!
!     CALCULATE DIPOLE MOMENT IF NEUTRAL MOLECULE
!
      if (iz /= 0) then
        go to 60
      else
        write (iw, 40)
   40   format(/,/,5x,'DIPOLE MOMENT EVALUATED FROM THE POINT CHARGES',/)
        dipx = 0.d0
        dipy = 0.d0
        dipz = 0.d0
        dipx = dipx + sum(co(1,:numat)*q(:numat)/a0)
        dipy = dipy + sum(co(2,:numat)*q(:numat)/a0)
        dipz = dipz + sum(co(3,:numat)*q(:numat)/a0)
        dip = sqrt(dipx**2 + dipy**2 + dipz**2)
        write (iw, '(12X,'' X        Y        Z       TOTAL'')')
        write (iw, '(8X,4F9.4)') dipx*cf, dipy*cf, dipz*cf, dip*cf
      end if
   60 continue
      if (index(keywrd,'SYMAVG') /= 0) then
        allocate (cequiv(numat,numat))
        do i = 1, numat
          do j = 1, numat
            cequiv(i,j) = .FALSE.
            if (abs(abs(q(i))-abs(q(j))) >= 1.D-5) cycle
            cequiv(i,j) = .TRUE.
          end do
        end do
        do i = 1, numat
          ieq = 0
          qsc(i) = 0.D0
          do j = 1, numat
            if (.not.cequiv(i,j)) cycle
            qsc(i) = qsc(i) + abs(q(j))
            ieq = ieq + 1
          end do
          q(i) = q(i)/abs(q(i))*qsc(i)/ieq
        end do
        write (iw, *) ' '
        write (iw, *) '   ELECTROSTATIC POTENTIAL CHARGES AVERAGED FOR'
        write (iw, *) '   SYMMETRY EQUIVALENT ATOMS'
        write (iw, *) ' '
        if (index(keywrd,'AM1')/=0 .or. index(keywrd,'PM3')/=0) then
          write (iw, '(7X,''ATOM NO.    TYPE    CHARGE'')')
          do i = 1, numat
            write (iw, '(9X,I4,9X,A2,1X,F10.4)') i, elemnt(nat(i)), q(i)
          end do
        else
          write (iw, &
      '(7X,''ATOM NO.    TYPE     CHARGE   SCALED CHARGE'')')
          do i = 1, numat
            write (iw, '(9X,I4,9X,A2,1X,F10.4,2X,F10.4)') i, elemnt(nat(i)), q(i), q(i)*slope
          end do
        end if
      end if
      return
      end subroutine potcal


      subroutine elesn()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : numat, keywrd, nopen, fract, nclose
      use overlaps_C, only : allc, allz
      use parameters_C, only : tore, zs, zp
      use chanel_C, only : iesp, esp_fn
      use common_arrays_C, only : c, nat
      use esp_C, only : b_esp, ovl, fac, &
      & cespm2, esp_array, cesp, cespm, es, cespml, potpt, &
      & cc, dex, fv, ex, cen, iam, ind, temp, itemp, nesp, co, fc, indc, icd, &
      & npr, ncc, is_esp, ip, ipe, ipx, isc
      use funcon_C, only : a0
!
      implicit none
!
      integer :: j1, i, m, k, j, nqn, ipc,  in, l
      double precision :: norm
      double precision :: pi,  t, fval, sum, ra, rij
      logical :: potwrt, sto3g
      double precision, external :: dex2
      pi = 4.D0*atan(1.D0)
!
!     PUT STO-6G BASIS SET ON ATOM CENTERS
!
      do i = -1, 10
        dex(i) = dex2(i)
      end do
      if (fac(4) > 1.D0) then
        fac(:7) = 1.D0/fac(:7)
      end if
      do m = 0, 8
        k = 1
        fv(m,1) = 1.D0/(2.D0*m + 1.D0)
        t = 0.D0
        do j1 = 1, int(41.D0/0.05D0)
          t = t + 0.05D0
          k = k + 1
          call fsub (m, t, fval)
          fv(m,k) = fval
        end do
        t = t + 0.05D0
      end do
!
!     LOAD BASIS FUNCTIONS INTO ARRAYS
!
      sto3g = index(keywrd,'STO3G') /= 0
      if (sto3g) then
        icd = 3
        call setup3
      else
        icd = 6
        call setupg
      end if
      ncc = 0
      npr = 0
      do i = 1, numat
        if (nat(i) <= 2) then
          cc(npr+1:icd+npr) = allc(:icd,1,1)
          ex(npr+1:icd+npr) = allz(:icd,1,1)*zs(1)**2
          cen(npr+1:icd+npr,1) = co(1,i)/a0
          cen(npr+1:icd+npr,2) = co(2,i)/a0
          cen(npr+1:icd+npr,3) = co(3,i)/a0
          iam(npr+1:icd+npr,1) = 0
          iam(npr+1:icd+npr,2) = 0
          fc(npr+1:icd+npr) = i
          ncc = ncc + 1
          npr = npr + icd
        else
!        DETERMINE PRINCIPAL QUANTUM NUMBER(NQN)
!        OF ORBITALS TO BE USED
!
          nqn = 2
          if (nat(i)>10 .and. nat(i)<=18) nqn = 3
          if (nat(i)>18 .and. nat(i)<=36) nqn = 4
          if (nat(i)>36 .and. nat(i)<=54) nqn = 5
!
          cc(npr+1:icd+npr) = allc(:icd,nqn,1)
          ex(npr+1:icd+npr) = allz(:icd,nqn,1)*zs(nat(i))**2
          cen(npr+1:icd+npr,1) = co(1,i)/a0
          cen(npr+1:icd+npr,2) = co(2,i)/a0
          cen(npr+1:icd+npr,3) = co(3,i)/a0
          iam(npr+1:icd+npr,1) = 0
          iam(npr+1:icd+npr,2) = 0
          ncc = ncc + 1
          npr = npr + icd
          do k = 1, 3
            cc(npr+1:icd+npr) = allc(:icd,nqn,2)
            ex(npr+1:icd+npr) = allz(:icd,nqn,2)*zp(nat(i))**2
            cen(npr+1:icd+npr,1) = co(1,i)/a0
            cen(npr+1:icd+npr,2) = co(2,i)/a0
            cen(npr+1:icd+npr,3) = co(3,i)/a0
            iam(npr+1:icd+npr,1) = 1
            iam(npr+1:icd+npr,2) = k
            ncc = ncc + 1
            npr = npr + icd
          end do
        end if
      end do
!
!     CALCULATE NORMALIZATION CONSTANTS AND INCLUDE
!     THEM IN THE CONTRACTION COEFFICIENTS
!
      do i = 1, npr
        norm = (2.D0*ex(i)/pi)**0.75D0*(4.D0*ex(i))**(iam(i,1)/2.D0)/sqrt(dex(2&
          *iam(i,1)-1))
        cc(i) = cc(i)*norm
      end do
!
!     PERFORM SORT OF PRIMITIVES BY ANGULAR MOMENTUM
!
      is_esp = 0
      isc = 0
      j = 0
      do i = 1, npr
        if (iam(i,1) /= 0) cycle
        is_esp = is_esp + 1
        ind(is_esp) = i
      end do
      ip = is_esp
      do i = 1, npr
        if (iam(i,1)/=1 .or. iam(i,2)/=1) cycle
        ip = ip + 1
        ind(ip) = i
      end do
      do i = 1, npr
        if (iam(i,1)/=1 .or. iam(i,2)/=2) cycle
        ip = ip + 1
        ind(ip) = i
      end do
      do i = 1, npr
        if (iam(i,1)/=1 .or. iam(i,2)/=3) cycle
        ip = ip + 1
        ind(ip) = i
      end do
      do i = 1, ncc
        in = i*icd - icd + 1
        if (iam(in,1) /= 0) cycle
        isc = isc + 1
        indc(isc) = i
      end do
      ipc = isc
      do i = 1, ncc
        in = i*icd - icd + 1
        if (iam(in,1)/=1 .or. iam(in,2)/=1) cycle
        ipc = ipc + 1
        indc(ipc) = i
      end do
      do i = 1, ncc
        in = i*icd - icd + 1
        if (iam(in,1)/=1 .or. iam(in,2)/=2) cycle
        ipc = ipc + 1
        indc(ipc) = i
      end do
      do i = 1, ncc
        in = i*icd - icd + 1
        if (iam(in,1)/=1 .or. iam(in,2)/=3) cycle
        ipc = ipc + 1
        indc(ipc) = i
      end do
      temp(:npr) = cc(ind(:npr))
      cc(:npr) = temp(:npr)
      temp(:npr) = ex(ind(:npr))
      ex(:npr) = temp(:npr)
      temp(:npr) = cen(ind(:npr),1)
      cen(:npr,1) = temp(:npr)
      temp(:npr) = cen(ind(:npr),2)
      cen(:npr,2) = temp(:npr)
      temp(:npr) = cen(ind(:npr),3)
      cen(:npr,3) = temp(:npr)
      itemp(:npr) = iam(ind(:npr),1)
      iam(:npr,1) = itemp(:npr)
      itemp(:npr) = iam(ind(:npr),2)
      iam(:npr,2) = itemp(:npr)
!     CALCULATE OVERLAP MATRIX OF STO-6G FUNCTIONS
!
      do j = 1, ncc
        call ovlp (j)
      end do
!
      do j = 1, ncc
        do k = 1, ncc
          cespm2(indc(j),indc(k)) = ovl(j,k)
        end do
      end do
      ovl(:ncc,:ncc) = cespm2(:ncc,:ncc)
      l = 0
      do i = 1, ncc
        if (i > 0) then
          cesp(l+1:i+l) = ovl(i,:i)
          l = i + l
        end if
      end do
!
!     DEORTHOGONALIZE THE COEFFICIENTS AND REFORM THE DENSITY MATRIX
!
      call rsp (cesp, ncc, temp, cespml)
      do i = 1, ncc ! TODO: to be parallelized
        do j = 1, i
          sum = 0.D0
          do k = 1, ncc
            sum = sum + cespml(i+(k-1)*ncc)/dsqrt(temp(k))*cespml(j+(k-1)*ncc)
            cesp(i+(j-1)*ncc) = sum
            cesp(j+(i-1)*ncc) = sum
          end do
        end do
      end do
      call mult (c, cesp, cespml, ncc)
      call densit (cespml, ncc, ncc, nclose, 2.d0, nopen, fract, cesp, 2)


!
!     NOW CALCULATE THE ELECTRONIC CONTRIBUTION TO THE ELECTROSTATIC POT
!
      l = 0
      do i = 1, ncc ! TODO: to be parallelized
        do j = 1, i
          l = l + 1
          cespm(i,j) = cesp(l)
          cespm(j,i) = cesp(l)
        end do
      end do
      ipx = (npr - is_esp)/3
      ipe = is_esp + ipx
      es(:nesp) = 0.D0
      call naicas () !TODO: to be parallelized
      call naicap () !TODO: to be parallelized
!     CALCULATE TOTAL ESP AND FORM ARRAYS FOR ESPFIT
      b_esp = 0.d0
      do i = 1, nesp
        esp_array(i) = 0.D0
        do j = 1, numat
          ra = sqrt((co(1,j)-potpt(1,i))**2+(co(2,j)-potpt(2,i))**2+(co(3,j) - potpt(3,i))**2)
          esp_array(i) = esp_array(i) + tore(nat(j))/(ra/a0)
        end do
        esp_array(i) = esp_array(i) - es(i)
        do j = 1, numat
          rij = sqrt((co(1,j)-potpt(1,i))**2+(co(2,j)-potpt(2,i))**2+(co(3,j) - potpt(3,i))**2)/a0
          b_esp(j) = b_esp(j) + esp_array(i)*1.D0/rij
        end do
      end do
!
!     IF REQUESTED WRITE OUT ELECTRIC POTENTIAL DATA TO
!     UNIT 21
!
      potwrt = index(keywrd,'POTWRT') /= 0
      if (potwrt) then
        open(iesp, file=esp_fn, status='UNKNOWN', position='asis')
        write (iesp, '(I5)') nesp
        do i = 1, nesp
          write (iesp, 420) esp_array(i), potpt(1,i)/a0, potpt(2,i)/a0, potpt(3,i)/a0
        end do
        close(unit=iesp, status='KEEP')
      end if
  420 format(1x,4e16.7)
      return
      end subroutine elesn
      subroutine espfit()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE esp_C, only : nesp, cf, co, potpt, idip, iz, dx, dy, dz, &
      & b_esp, esp_array, rms, rrms
      use common_arrays_C, only : q
      use molkst_C, only : numat
      use funcon_C, only : a0, fpc_1, fpc_8
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, j, i, info
      double precision ::rik, rij, espc, atmp
      integer, allocatable :: ipiv(:)
      double precision, allocatable :: a(:,:)
    !
    !.. Local Scalars ..
    integer ::  istart, iend
    cf = fpc_1 * a0 * fpc_8 * 1.d-10
     allocate (a(numat+4,numat+4))
    !
    !     THE FOLLOWING SETS UP THE LINEAR EQUATION A*Q=B
    !     SET UP THE A(J,K) ARRAY
    !


    istart = 1
    iend = nesp
    do k = 1, numat
      do j = 1, numat
        atmp = 0
        do i = istart, iend
          rik = (co(1, k) - potpt(1, i)) ** 2 &
               & + (co(2, k) - potpt(2, i)) ** 2 &
               & + (co(3, k) - potpt(3, i)) ** 2
          rij = (co(1, j) - potpt(1, i)) ** 2 &
               & + (co(2, j) - potpt(2, i)) ** 2 &
               & + (co(3, j) - potpt(3, i)) ** 2
          atmp = atmp + 1 / Sqrt (rij * rik)
        end do
        a(j, k) = atmp * a0**2
      end do
    end do
      a(numat+1, 1:numat) = 1.d0 !row
      a(numat+1:, numat+1:) = 0.d0 !corner
      a(1:numat, numat+1) = 1.d0 !column
      if (idip == 1) then
        a(numat+2:numat+4, 1:numat) = co(1:3, 1:numat) / a0
        a(1:numat, numat+2:numat+4) = Transpose (co(1:3, 1:numat) / a0)
        a(numat+2:numat+4, numat+2:numat+4) = 0.d0
      else
        a(:, numat+2:) = 0.d0
        a(numat+2:, 1:numat+1) = 0.d0
      end if
    !
    !     INSERT CHARGE AND DIPOLAR (IF DESIRED) CONSTRAINTS
    !
    b_esp(numat+1) = iz
    b_esp(numat+2) = dx / cf
    b_esp(numat+3) = dy / cf
    b_esp(numat+4) = dz / cf
    !
    if (allocated(q)) deallocate(q)
    allocate(q(numat + 4))
    q(1:numat+4) = b_esp(1:numat+4)
    !
    allocate (ipiv(numat+4))
    !
    if (idip == 1) then
      call dgesv (numat+4, 1, a, numat+4, ipiv, q, numat+4, info)
    else
      call dgesv (numat+1, 1, a, numat+4, ipiv, q, numat+4, info)
    end if
    !
    deallocate (ipiv)
    !
    !     CALCULATE ROOT MEAN SQUARE FITS AND RELATIVE ROOT MEAN SQUARE FITS
    !
    do i = 1, nesp
      espc = 0.d0
      do j = 1, numat
        rij = Sqrt ((co(1, j) - potpt(1, i)) ** 2 &
             & + (co(2, j) - potpt(2, i)) ** 2 &
             & + (co(3, j) - potpt(3, i)) ** 2) / a0
        espc = espc + q(j) / rij
      end do
      rms = rms + (espc-esp_array(i)) ** 2
      rrms = rrms + esp_array(i) ** 2
    end do
    rms = Sqrt (rms/nesp)
    rrms = rms / Sqrt (rrms/nesp)
    rms = rms * 627.51d0
  end subroutine espfit


      subroutine naicas()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use chanel_C, only : iw, iesr, restart_fn
      use funcon_C, only : a0
      use molkst_C, only : keywrd
      use esp_C, only : cen, exsr, ex, tf, fv, iam, cc, ip, npr, fac, &
      & ncc, ipe, cespm, espi, potpt, es, dx_array, dy_array, dz_array, f1, td, &
      & u_esp, expn, rnai, ewcx, ewcy, ewcz, f0, rnai1, ipx, is_esp, isc, &
      & icd, exs, ce, indc, nesp
!***********************************************************************
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ipx2, i, jstart, iesps, np, ic, ipr, istart, j, iesp, iref, k, ips, jps
      double precision :: pi, potp1, potp2, potp3, ref, res, term, f, ts, term1, fi, ts1
!-----------------------------------------------
      a0 = a0
!
!     CALCULATE DISTANCE ARRAYS
!
      write (iw, *)
      pi = 4.D0*atan(1.D0)
      ipx2 = 2*ipx
!     IF THIS IS A RESTART RUN, READ IN RESTART INFO
      if (index(keywrd,'ESPRST') /= 0) then
        open(unit=iesr, file=restart_fn, status='OLD', form=&
          'UNFORMATTED', position='asis')
        read (iesr) jstart, iesps
        if (jstart == isc*2) then
          close(iesr)
          return
        end if
        do i = 1, nesp
          read (iesr) es(i)
        end do
        close(iesr)
!
        jstart = jstart + 1
      else
        jstart = 1
      end if
      np = is_esp + 1
      do ic = jstart, isc
        ipr = ic*icd - icd + 1
        istart = ipr
        dx_array(istart:ipe) = cen(ipr,1) - cen(istart:ipe,1)
        dy_array(istart:ipe) = cen(ipr,2) - cen(istart:ipe,2)
        dz_array(istart:ipe) = cen(ipr,3) - cen(istart:ipe,3)
        td(istart:ipe) = dx_array(istart:ipe)**2 + dy_array(istart:ipe)**2 + dz_array(istart:ipe)**2
!
!     CALCULATE EXPONENT SUM
!
        do i = istart, ipe
          do j = 1, icd
            exsr(i,j) = ex(ipr+j-1) + ex(i)
            exs(i,j) = 1.D0/exsr(i,j)
            ce(i,j) = ex(ipr+j-1)*ex(i)*exs(i,j)
            expn(i,j) = exp((-ce(i,j)*td(i)))
          end do
        end do
!
!     CALCULATE EXPONENT WEIGHTED CENTERS
!
        do i = istart, ipe
          ewcx(i,:icd) = (ex(i)*cen(i,1)+ex(ipr:icd-1+ipr)*cen(ipr:icd-1+ipr,1)&
            )*exs(i,:icd)
          ewcy(i,:icd) = (ex(i)*cen(i,2)+ex(ipr:icd-1+ipr)*cen(ipr:icd-1+ipr,2)&
            )*exs(i,:icd)
          ewcz(i,:icd) = (ex(i)*cen(i,3)+ex(ipr:icd-1+ipr)*cen(ipr:icd-1+ipr,3)&
            )*exs(i,:icd)
        end do
!
!     BEGIN LOOP OVER ESP POINTS
!
        do iesp = 1, nesp
          potp1 = potpt(1,iesp)/a0
          potp2 = potpt(2,iesp)/a0
          potp3 = potpt(3,iesp)/a0
!
!     BEGIN LOOP OVER COMPONENTS OF CONTRACTED FUNCTION IC
!
          do j = 1, icd
!
!     CALCULATE DISTANCE BETWEEN EXPONENT WEIGHTED AND PROBE POINT
!
            do i = istart, ipe
              u_esp(i,j) = ((ewcx(i,j)-potp1)**2+(ewcy(i,j)-potp2)**2+(ewcz(i,j)-&
                potp3)**2)*exsr(i,j)
              rnai(i,j) = sqrt(pi/u_esp(i,j))
            end do
!
!     CALCULATE ESP INTEGRALS
!
            do i = istart, ipe
              if (u_esp(i,j) <= tf(0)) then
                iref = nint(u_esp(i,j)*20.D0)
                ref = 0.05D0*iref
                res = u_esp(i,j) - ref
                term = 1.D0
                f0(i,j) = 0.D0
                do k = 0, 6
                  f = fv(k,iref+1)
                  ts = f*term*fac(k)
                  term = -term*res
                  f0(i,j) = f0(i,j) + ts
                end do
              else
                f0(i,j) = rnai(i,j)*0.5D0
              end if
            end do
            do i = np, ipe
              if (u_esp(i,j) <= tf(1)) then
                iref = nint(u_esp(i,j)*20.D0)
                ref = 0.05D0*iref
                res = u_esp(i,j) - ref
                term1 = 1.D0
                f1(i,j) = 0.D0
                do k = 0, 6
                  fi = fv(k+1,iref+1)
                  ts1 = fi*term1*fac(k)
                  term1 = -term1*res
                  f1(i,j) = f1(i,j) + ts1
                end do
              else
                f1(i,j) = rnai(i,j)*0.25D0/u_esp(i,j)
              end if
            end do
            u_esp(istart:is_esp,j) = 2.D0*pi*exs(istart:is_esp,j)*expn(istart:is_esp,j)* &
            f0(istart:is_esp,j)
            np = is_esp + 1
            rnai(np:ipe,j) = 2.D0*pi*exs(np:ipe,j)*expn(np:ipe,j)*f0(np:ipe,j)
            rnai1(np:ipe,j) = 2.D0*pi*exs(np:ipe,j)*expn(np:ipe,j)*f1(np:ipe,j)
!
!     CALCULATE (S||P) ESP INTEGRALS
!
            if (iam(ipr,1)/=0 .or. is_esp==ip) cycle
            u_esp(np:ipe,j) = (ewcx(np:ipe,j)-cen(np:ipe,1))*rnai(np:ipe,j) - &
            (ewcx(np:ipe,j)-potp1)*rnai1(np:ipe,j)
            u_esp(ipe+1:1+ipx+ipe,j) = (ewcy(ipe+1-ipx:1+ipe,j)-cen(ipe+1-ipx:1+ipe,2))* &
            rnai(ipe+1-ipx:1+ipe,j) - (ewcy(ipe+1-ipx:1+ipe,j)-potp2)*&
              rnai1(ipe+1-ipx:1+ipe,j)
            u_esp(ipe+1+ipx:npr,j) = (ewcz(ipe+1+ipx-ipx2:npr-ipx2,j) - &
            cen(ipe+1+ipx-ipx2:npr-ipx2,3))*rnai(ipe+1+ipx-ipx2:npr-ipx2,j) - &
            (ewcz(ipe+1+ipx-ipx2:npr-ipx2,j)-potp3)*rnai1(ipe+1+ipx-ipx2:npr-ipx2,j)
          end do
          ips = ic*icd - icd + 1
          do i = ic, ncc
            jps = i*icd - icd + 1
            espi(i,ic) = 0.D0
            do j = jps, jps + icd - 1
              espi(i,ic) = espi(i,ic) + sum(cc(j)*cc(ips:icd-1+ips)*u_esp(j,:icd))
            end do
            es(iesp) = es(iesp) + 2.D0*cespm(indc(i),indc(ic))*espi(i,ic)
          end do
          es(iesp) = es(iesp) - cespm(indc(ic),indc(ic))*espi(ic,ic)
        end do
!     WRITE OUT RESTART INFORMATION
        open(unit=iesr, file=restart_fn, status='UNKNOWN', form='UNFORMATTED', position='asis')
        iesps = 0
        write (iesr) ic, iesps
        do i = 1, nesp
          write (iesr) es(i)
        end do
        close(iesr)
!
        write (iw, '(A,F6.2,A)') 'NAICAS DUMPED: ', 100.D0/isc*ic, &
          ' PERCENT COMPLETE'
      end do
      return
      end subroutine naicas


      subroutine ovlp(ic)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use esp_C, only : cen, ex, iam, cc, is_esp, ip, npr, ncc, icd, ovl, &
      &  dx_array, dy_array, dz_array, td, ce, exs, expn, rnai, ewcx, ewcy, ewcz
!***********************************************************************
      implicit none
      integer :: ic
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ipr, istart, i, j, np, ips, jps
      double precision :: pi
!-----------------------------------------------
!
!     CALCULATE DISTANCE ARRAYS
!
      pi = 4.D0*atan(1.D0)
      ipr = ic*icd - icd + 1
      istart = ipr
      dx_array(istart:npr) = cen(ipr,1) - cen(istart:npr,1)
      dy_array(istart:npr) = cen(ipr,2) - cen(istart:npr,2)
      dz_array(istart:npr) = cen(ipr,3) - cen(istart:npr,3)
      td(istart:npr) = dx_array(istart:npr)**2 + dy_array(istart:npr)**2 + dz_array(istart:npr)**2
!
!     CALCULATE EXPONENT SUM
!
      do i = istart, npr
        exs(i,:icd) = 1.D0/(ex(ipr:icd-1+ipr)+ex(i))
        ce(i,:icd) = ex(ipr:icd-1+ipr)*ex(i)*exs(i,:icd)
!
!     CALCULATE EXPONENT WEIGHTED CENTERS
!
        ewcx(i,:icd) = (ex(i)*cen(i,1)+ex(ipr:icd-1+ipr)*cen(ipr:icd-1+ipr,1))*&
          exs(i,:icd)
        ewcy(i,:icd) = (ex(i)*cen(i,2)+ex(ipr:icd-1+ipr)*cen(ipr:icd-1+ipr,2))*&
          exs(i,:icd)
        ewcz(i,:icd) = (ex(i)*cen(i,3)+ex(ipr:icd-1+ipr)*cen(ipr:icd-1+ipr,3))*&
          exs(i,:icd)
      end do
      do i = 1, npr
        do j = 1, icd
          expn(i,j) = exp((-ce(i,j)*td(i)))
          rnai(i,j) = (pi*exs(i,j))**1.5D0*expn(i,j)
          expn(i,j) = rnai(i,j)
        end do
      end do
!
!     CALCULATE (S||P) ESP INTEGRALS
!
      if (iam(ipr,1)==0 .and. is_esp/=ip) then
        np = is_esp + 1
        do i = np, npr
          do j = 1, icd
            select case (iam(i,2))
            case default
              rnai(i,j) = (ewcx(i,j)-cen(i,1))*expn(i,j)
            case (2)
              rnai(i,j) = (ewcy(i,j)-cen(i,2))*expn(i,j)
            case (3)
              rnai(i,j) = (ewcz(i,j)-cen(i,3))*expn(i,j)
            end select
          end do
        end do
      end if
!
!     CALCULATE (P||S) ESP INTEGRALS
!
      if (iam(ipr,1)==1 .and. is_esp/=ip) then
        np = is_esp + 1
        do i = istart, npr
          do j = 1, icd
            select case (iam(ipr+j-1,2))
            case default
              rnai(i,j) = (ewcx(i,j)-cen(ipr+j-1,1))*expn(i,j)
            case (2)
              rnai(i,j) = (ewcy(i,j)-cen(ipr+j-1,2))*expn(i,j)
            case (3)
              rnai(i,j) = (ewcz(i,j)-cen(ipr+j-1,3))*expn(i,j)
            end select
          end do
        end do
        do i = istart, npr
          do j = 1, icd
            select case (iam(i,2))
            case default
              rnai(i,j) = (ewcx(i,j)-cen(i,1))*rnai(i,j)
              if (iam(ipr+j-1,2) == iam(i,2)) rnai(i,j) = rnai(i,j) + exs(i,j)*&
                0.5D0*expn(i,j)
            case (2)
              rnai(i,j) = (ewcy(i,j)-cen(i,2))*rnai(i,j)
              if (iam(ipr+j-1,2) == iam(i,2)) rnai(i,j) = rnai(i,j) + exs(i,j)*&
                0.5D0*expn(i,j)
            case (3)
              rnai(i,j) = (ewcz(i,j)-cen(i,3))*rnai(i,j)
              if (iam(ipr+j-1,2) == iam(i,2)) rnai(i,j) = rnai(i,j) + exs(i,j)*&
                0.5D0*expn(i,j)
            end select
          end do
        end do
      end if
      ips = ic*icd - icd + 1
      do i = ic, ncc
        jps = i*icd - icd + 1
        ovl(ic,i) = 0.D0
        do j = jps, jps + icd - 1
          ovl(ic,i) = ovl(ic,i) + sum(cc(j)*cc(ips:icd-1+ips)*rnai(j,:icd))
        end do
        ovl(i,ic) = ovl(ic,i)
      end do
      return
      end subroutine ovlp
      subroutine setup3
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE overlaps_C, only : allc, allz
!***********************************************************************
      implicit none
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!     SET-UP THE STEWART'S STO-3G EXPANSIONS
!     FROM J. CHEM. PHYS. 52 431.
!                                            1S
      allz(1,1,1) = 2.227660584D00
      allz(2,1,1) = 4.057711562D-01
      allz(3,1,1) = 1.098175104D-01
!
      allc(1,1,1) = 1.543289673D-01
      allc(2,1,1) = 5.353281423D-01
      allc(3,1,1) = 4.446345422D-01
!                                      2S
      allz(1,2,1) = 2.581578398D00
      allz(2,2,1) = 1.567622104D-01
      allz(3,2,1) = 6.018332272D-02
!
      allc(1,2,1) = -5.994474934D-02
      allc(2,2,1) = 5.960385398D-01
      allc(3,2,1) = 4.581786291D-01
!                                     2P
      allz(1,2,2) = 9.192379002D-01
      allz(2,2,2) = 2.359194503D-01
      allz(3,2,2) = 8.009805746D-02
!
      allc(1,2,2) = 1.623948553D-01
      allc(2,2,2) = 5.661708862D-01
      allc(3,2,2) = 4.223071752D-01
!                                      3S
      allz(1,3,1) = 5.641487709D-01
      allz(2,3,1) = 6.924421391D-02
      allz(3,3,1) = 3.269529097D-02
!
      allc(1,3,1) = -1.782577972D-01
      allc(2,3,1) = 8.612761663D-01
      allc(3,3,1) = 2.261841969D-01
!                                     3P
      allz(1,3,2) = 2.692880368D00
      allz(2,3,2) = 1.489359592D-01
      allz(3,3,2) = 5.739585040D-02
!
      allc(1,3,2) = -1.061945788D-02
      allc(2,3,2) = 5.218564264D-01
      allc(3,3,2) = 5.450015143D-01
!                                      4S
      allz(1,4,1) = 2.267938753D-01
      allz(2,4,1) = 4.448178019D-02
      allz(3,4,1) = 2.195294664D-02
!
      allc(1,4,1) = -3.349048323D-01
      allc(2,4,1) = 1.056744667D00
      allc(3,4,1) = 1.256661680D-01
!                                     4P
      allz(1,4,2) = 4.859692220D-01
      allz(2,4,2) = 7.430216918D-02
      allz(3,4,2) = 3.653340923D-02
!
      allc(1,4,2) = -6.147823411D-02
      allc(2,4,2) = 6.604172234D-01
      allc(3,4,2) = 3.932639495D-01
!                                      5S
      allz(1,5,1) = 1.080198458D-01
      allz(2,5,1) = 4.408119382D-02
      allz(3,5,1) = 2.610811810D-02
!
      allc(1,5,1) = -6.617401158D-01
      allc(2,5,1) = 7.467595004D-01
      allc(3,5,1) = 7.146490945D-01
!                                     5P
      allz(1,5,2) = 2.127482317D-01
      allz(2,5,2) = 4.729648620D-02
      allz(3,5,2) = 2.604865324D-02
!
      allc(1,5,2) = -1.389529695D-01
      allc(2,5,2) = 8.076691064D-01
      allc(3,5,2) = 2.726029342D-01
!
      return
      end subroutine setup3

!****************************************************************
      double precision function dex2 (m)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: m
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------
      if (m < 2) then
        dex2 = 1
      else
        dex2 = 1
        do i = 1, m, 2
          dex2 = dex2*i
        end do
      end if
      return
      end function dex2

!***********************************************************************
      subroutine fsub(n, x, fval)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      double precision , intent(in) :: x
      double precision , intent(out) :: fval
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ku, k, i
      double precision, dimension(21) :: ff
      double precision, dimension(200) :: term
      double precision :: a0, a1s2, pie4, a1, xsw, e, fac0, fact, term0, sum, suma, &
        sum1
!-----------------------------------------------
!***********************************************************************
!
!     CALCULATE THE FM(T). KINDLY SUPPLIED BY RUS PITZER AND CLEANED UP
!     BY K.M.MERZ IN FEB. 1989 AT UCSF.
!
!     ON INPUT:  N = INDEX
!                X = EXPONENT
!     ON OUTPUT: FVAL = VALUE OF THE FUNCTION
!
!     FOR MORE DETAILS SEE: OBARA AND SAIKA J. CHEM. PHYS. 1986,84,3963
!***********************************************************************
      data a0, a1s2, pie4, a1/ 0.0D0, 0.5D0, 0.7853981633974483096156608D0, &
        1.0D0/
      data xsw/ 24.0D0/
      e = a1s2*exp((-x))
      fac0 = n
      fac0 = fac0 + a1s2
      if (x > xsw) go to 50
!
!     USE POWER SERIES
!
   10 continue
      fact = fac0
      term0 = e/fact
      sum = term0
      ku = int(x - fac0)
      if (ku >= 1) then
!
!     SUM INCREASING TERMS FORWARDS
!
        do k = 1, ku
          fact = fact + a1
          term0 = term0*x/fact
          sum = sum + term0
        end do
      end if
      i = 1
      fact = fact + a1
      term(1) = term0*x/fact
      suma = sum + term(1)
      if (Abs(sum - suma) < 1.d-20) go to 90
      i = i + 1
      fact = fact + a1
      term(i) = term(i-1)*x/fact
      sum1 = suma
      suma = suma + term(i)
      do while(sum1 - suma /= 0.D0)
        i = i + 1
        fact = fact + a1
        term(i) = term(i-1)*x/fact
        sum1 = suma
        suma = suma + term(i)
      end do
      go to 90
!
!     USE ASYMPTOTIC SERIES
!
   50 continue
      sum = sqrt(pie4/x)
      if (n /= 0) then
        fact = -a1s2
        do k = 1, n
          fact = fact + a1
          sum = sum*fact/x
        end do
      end if
      i = 1
      term(1) = -e/x
      suma = sum + term(1)
      if (Abs(sum - suma) < 1.d-20) go to 90
      fact = fac0
      ku = int(x + fac0 - a1)
      do i = 2, ku
        fact = fact - a1
        term(i) = term(i-1)*fact/x
        sum1 = suma
        suma = suma + term(i)
        if (Abs(sum - suma) < 1.d-20) go to 90
      end do
!
!     XSW SET TOO LOW. USE POWER SERIES.
!
      go to 10
!
!     SUM DECREASING TERMS BACKWARDS
!
   90 continue
      sum1 = a0
      do k = 1, i
        sum1 = sum1 + term(i+1-k)
      end do
      ff(n+1) = sum + sum1
!
!     USE RECURRENCE RELATION
!
      if (n /= 0) then
        do k = 1, n
          fac0 = fac0 - a1
          ff(n+1-k) = (e + x*ff(n+2-k))/fac0
        end do
      end if
      fval = ff(n+1)
      return
      end subroutine fsub

!****************************************************************
      subroutine naicap()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use funcon_C, only : a0
      use chanel_C, only : iw, iesr, esr_fn
      use esp_C, only : cen, ex, tf, fv, iam, cc, isc, is_esp, ip, npr, ncc, ipe, &
      & ipx, icd, cespm, espi, potpt, es, rnai, rnai2, pf0, pf1, pf2, fac, &
      & pexs, pce, pexpn, ptd, pewcx, pewcy, pewcz, ird, indc, nesp
      use molkst_C, only : keywrd
!***********************************************************************
      implicit none

!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer,  allocatable, dimension (:) :: id
      integer :: idn, idc, ipx2, np, i, l, j, jstart, iesps, iesp, il&
        , iref, k, ic, ipr, istart, in, ir, ir2, ips, jps, mptd
      double precision :: pi, potp1, potp2, potp3, ref, res, term, f, ts, &
        term1, fi, ts1, term2, fii, ts2
!-----------------------------------------------
      idn = 10
      mptd = ((ipe - is_esp)*(ipe - is_esp + 1))/2
      allocate(id(ipe - is_esp))
      allocate (pf0(mptd), pf1(mptd), pf2(mptd), pexs(mptd), pce(mptd), &
         & pexpn(mptd), ptd(mptd), pewcx(mptd), pewcy(mptd), pewcz(mptd))
!
      idc = 0
      write (iw, *)
      ipx2 = 2*ipx
      pi = 4.D0*atan(1.D0)
      np = is_esp + 1
!     SETUP INDEX ARRAY
      do i = np, ipe
        ird(i) = i - is_esp
        ird(i+ipx) = i - is_esp
        ird(i+ipx2) = i - is_esp
      end do
!
!     CALCULATE QUANTITIES INVARIANT WITH ESP POINT FOR
!     (P|P) ESP INTEGRALS
!
!      IL=L
      l = 0
      do i = np, ipe
        do j = i, ipe
          l = l + 1
          ptd(l) = (cen(i,1)-cen(j,1))**2 + (cen(i,2)-cen(j,2))**2 + (cen(i,3)-&
            cen(j,3))**2
          pexs(l) = 1.D0/(ex(i)+ex(j))
          pce(l) = ex(i)*ex(j)*pexs(l)
          pexpn(l) = exp((-pce(l)*ptd(l)))
          pewcx(l) = (ex(i)*cen(i,1)+ex(j)*cen(j,1))*pexs(l)
          pewcy(l) = (ex(i)*cen(i,2)+ex(j)*cen(j,2))*pexs(l)
          pewcz(l) = (ex(i)*cen(i,3)+ex(j)*cen(j,3))*pexs(l)
        end do
!
!     SET UP OTHER INDEX ARRAY FOR PACKED SYMMETRIC ARRAY
!     STORAGE
!
        id(i-is_esp) = l - ipx
      end do
!
!     READ IN RESTART INFORMATION IF THIS IS A RESTART
!
      if (index(keywrd,'ESPRST') /= 0) then
        open(unit=iesr, file=esr_fn, status='UNKNOWN', form=&
          'UNFORMATTED', position='asis')
        read (iesr) jstart, iesps
        if (jstart /= isc*2) then
          iesps = 0
          close(iesr)
          go to 50
        end if
        do i = 1, nesp
          read (iesr) es(i)
        end do
        close(iesr)
        idc = int((iesps/float(nesp))*10)
      else
        iesps = 0
      end if
   50 continue
      do iesp = iesps + 1, nesp
        potp1 = potpt(1,iesp)/a0
        potp2 = potpt(2,iesp)/a0
        potp3 = potpt(3,iesp)/a0
!     CALCULATE QUANTITY U
!
        l = 0
        do i = np, ipe
          do j = i, ipe
            l = l + 1
            ptd(l) = ((pewcx(l)-potp1)**2+(pewcy(l)-potp2)**2+(pewcz(l)-potp3)**2)/pexs(l)
            pce(l) = sqrt(pi/ptd(l))
          end do
        end do
!
!     CALCULATE F0, F1, AND F2(U) USING TAYLOR SERIES
!     OR ASYMPTOTIC EXPANSION
!
        il = l
        l = 0
        do i = 1, il
          if (ptd(i) <= tf(0)) then
            iref = nint(ptd(i)*20.D0)
            ref = 0.05D0*iref
            res = ptd(i) - ref
            term = 1.D0
            pf0(i) = 0.D0
            do k = 0, 6
              f = fv(k,iref+1)
              ts = f*term*fac(k)
              term = -term*res
              pf0(i) = pf0(i) + ts
            end do
          else
            pf0(i) = pce(i)*0.5D0
          end if
          if (ptd(i) <= tf(1)) then
            iref = nint(ptd(i)*20.D0)
            ref = 0.05D0*iref
            res = ptd(i) - ref
            term1 = 1.D0
            pf1(i) = 0.D0
            do k = 0, 6
              fi = fv(k+1,iref+1)
              ts1 = fi*term1*fac(k)
              term1 = -term1*res
              pf1(i) = pf1(i) + ts1
            end do
          else
            pf1(i) = pce(i)*0.25D0/ptd(i)
          end if
          if (ptd(i) <= tf(2)) then
            iref = nint(ptd(i)*20.D0)
            ref = 0.05D0*iref
            res = ptd(i) - ref
            term2 = 1.D0
            pf2(i) = 0.D0
            do k = 0, 6
              fii = fv(k+2,iref+1)
              ts2 = fii*term2*fac(k)
              term2 = -term2*res
              pf2(i) = pf2(i) + ts2
            end do
          else
            pf2(i) = pce(i)*0.375D0/(ptd(i)*ptd(i))
          end if
        end do
!
!     CALCULATE (S||S) TYPE INTEGRALS
!
        pf0(:il) = 2.D0*pi*pexs(:il)*pexpn(:il)*pf0(:il)
        ptd(:il) = pf0(:il)
        pf1(:il) = 2.D0*pi*pexs(:il)*pexpn(:il)*pf1(:il)
        pf2(:il) = 2.D0*pi*pexs(:il)*pexpn(:il)*pf2(:il)
!
        do ic = isc + 1, ncc
          ipr = ic*icd - icd + 1
          istart = ipr
          do j = 1, icd
!
!     CALCULATE (P||S) ESP INTEGRALS
!
            if (iam(ipr,1)==1 .and. is_esp/=ip) then
              do i = istart, npr
                in = ipr + j - 1
                ir = ird(i) + id(ird(in))
                ir2 = id(ird(i)) + ird(in)
                ir = min0(ir2,ir)
                select case (iam(in,2))
                case default
                  rnai2(i,j) = (pewcx(ir)-cen(in,1))*pf1(ir) - pf2(ir)*(pewcx(ir)-potp1)
                  rnai(i,j) = (pewcx(ir)-cen(in,1))*pf0(ir) - pf1(ir)*(pewcx(ir)-potp1)
                case (2)
                  rnai2(i,j) = (pewcy(ir)-cen(in,2))*pf1(ir) - pf2(ir)*(pewcy(ir)-potp2)
                  rnai(i,j) = (pewcy(ir)-cen(in,2))*pf0(ir) - pf1(ir)*(pewcy(ir)-potp2)
                case (3)
                  rnai2(i,j) = (pewcz(ir)-cen(in,3))*pf1(ir) - pf2(ir)*(pewcz(ir)-potp3)
                  rnai(i,j) = (pewcz(ir)-cen(in,3))*pf0(ir) - pf1(ir)*(pewcz(ir)-potp3)
                end select
              end do
            end if
!
!     CALCULATE (P||P) ESP INTEGRALS
!
            if (iam(ipr,1)/=1 .or. is_esp==ip) cycle
            do i = istart, npr
              in = ipr + j - 1
              ir = ird(i) + id(ird(in))
              ir2 = id(ird(i)) + ird(in)
              ir = min0(ir2,ir)
              select case (iam(i,2))
              case default
                rnai(i,j) = (pewcx(ir)-cen(i,1))*rnai(i,j) - (pewcx(ir)-potp1)*rnai2(i,j)
                if (iam(in,2) == iam(i,2)) rnai(i,j) = rnai(i,j) + pexs(ir)*0.5D0*(ptd(ir)-pf1(ir))
              case (2)
                rnai(i,j) = (pewcy(ir)-cen(i,2))*rnai(i,j) - (pewcy(ir)-potp2)*rnai2(i,j)
                if (iam(in,2) == iam(i,2)) rnai(i,j) = rnai(i,j) + pexs(ir)*0.5D0*(ptd(ir)-pf1(ir))
              case (3)
                rnai(i,j) = (pewcz(ir)-cen(i,3))*rnai(i,j) - (pewcz(ir)-potp3)*rnai2(i,j)
                if (iam(in,2) == iam(i,2)) rnai(i,j) = rnai(i,j) + pexs(ir)*0.5D0*(ptd(ir)-pf1(ir))
              end select
            end do
          end do
!
!     FORM INTEGRALS OVER CONTRACTED FUNCTIONS
!
          ips = ic*icd - icd + 1
          do i = ic, ncc
            jps = i*icd - icd + 1
            espi(i,ic) = 0.D0
            do j = jps, jps + icd - 1
              espi(i,ic) = espi(i,ic) + sum(cc(j)*cc(ips:icd-1+ips)*rnai(j,:icd))
            end do
            es(iesp) = es(iesp) + 2.D0*cespm(indc(i),indc(ic))*espi(i,ic)
          end do
          es(iesp) = es(iesp) - cespm(indc(ic),indc(ic))*espi(ic,ic)
        end do
!
!     WRITE OUT RESTART INFORMATION EVERY NESP/10 POINTS
!
        if (mod(iesp,nesp/idn) /= 0) cycle
        open(unit=iesr, file=esr_fn, status='UNKNOWN', form=&
          'UNFORMATTED', position='asis')
        jstart = isc*2
        write (iesr) jstart, iesp
        do i = 1, nesp
          write (iesr) es(i)
        end do
        close(iesr)
        idc = idc + 1
        write (iw, '(A,F6.2,A)') 'NAICAP DUMPED: ', 100.D0/idn*idc, &
          ' PERCENT COMPLETE'
      end do
      deallocate(id, pf0, pf1, pf2, pexs, pce, pexpn, ptd, pewcx, pewcy, pewcz)
      return
      end subroutine naicap
  subroutine setup_esp(mode)
    use esp_C, only : ind, itemp, ird, iam, cc, ex, temp, rad, es, &
      & b_esp, esp_array, cesp, cespml, al, qesp, td, dx_array, dy_array, dz_array, &
      & cespm, cen, fc, co, indc, &
      & potpt, ovl, cespm2, a2, exs, ce, expn, ewcx, ewcy, ewcz, f0, f1, u_esp, &
      & rnai, rnai1, rnai2, espi,  exsr, qsc
    use molkst_C, only : norbs, numat
    implicit none
    integer :: i, mode
    integer, parameter :: mesp=5000
    if (mode == 1) then
      i = 6*norbs
      if (allocated(a2))        deallocate(a2)
      if (allocated(ovl))       deallocate(ovl)
      if (allocated(cespm2))    deallocate(cespm2)
      if (allocated(cespml))    deallocate(cespml)
      if (allocated(cesp))      deallocate(cesp)
      if (allocated(potpt))     deallocate(potpt)
      if (allocated(es))        deallocate(es)
      if (allocated(co))        deallocate(co)
      if (allocated(al))        deallocate(al)
      if (allocated(rad))       deallocate(rad)
      if (allocated(b_esp))     deallocate(b_esp)
      if (allocated(qesp))      deallocate(qesp)
      if (allocated(cespm))     deallocate(cespm)
      if (allocated(cc))        deallocate(cc)
      if (allocated(cen))       deallocate(cen)
      if (allocated(iam))       deallocate(iam)
      if (allocated(ind))       deallocate(ind)
      if (allocated(ex))        deallocate(ex)
      if (allocated(temp))      deallocate(temp)
      if (allocated(f0))        deallocate(f0)
      if (allocated(f1))        deallocate(f1)
      if (allocated(ce))        deallocate(ce)
      if (allocated(u_esp))     deallocate(u_esp)
      if (allocated(exs))       deallocate(exs)
      if (allocated(expn))      deallocate(expn)
      if (allocated(cesp))      deallocate(cesp)
      if (allocated(rnai))      deallocate(rnai)
      if (allocated(ewcx))      deallocate(ewcx)
      if (allocated(ewcy))      deallocate(ewcy)
      if (allocated(ewcz))      deallocate(ewcz)
      if (allocated(rnai1))     deallocate(rnai1)
      if (allocated(rnai2))     deallocate(rnai2)
      if (allocated(fc))        deallocate(fc)
      if (allocated(exsr))      deallocate(exsr)
      if (allocated(indc))      deallocate(indc)
      if (allocated(itemp))     deallocate(itemp)
      if (allocated(dx_array))  deallocate(dx_array)
      if (allocated(dy_array))  deallocate(dy_array)
      if (allocated(dz_array))  deallocate(dz_array)
      if (allocated(td))        deallocate(td)
      if (allocated(esp_array)) deallocate(esp_array)
      if (allocated(qsc))       deallocate(qsc)
      if (allocated(espi))      deallocate(espi)
      if (allocated(ird))       deallocate(ird)
      allocate(a2((numat+4), (numat+4)), ovl(norbs, norbs), cespm2(norbs, norbs), cespml(norbs**2), &
      & cesp(norbs**2), potpt(3, mesp), es(mesp), co(3,numat), al((numat+4)**2),  &
      & rad(numat), b_esp(numat+4), qesp(numat+4), cespm(norbs, norbs))
      allocate( cc(i), cen(i,3), iam(i,2), ind(i), ex(i), temp(i), f0(i,6), &
      & f1(i,6), ce(i,6), u_esp(i,6), exs(i,6), expn(i,6), rnai(i,6), ewcx(i,6), ewcy(i,6), &
      & ewcz(i,6), rnai1(i,6), rnai2(i,6), fc(36*norbs), exsr(i,6), indc(norbs))
      allocate( itemp(i), dx_array(i), dy_array(i), dz_array(i), td(i), esp_array(mesp), &
      & qsc(numat), espi(norbs, norbs), &
      & ird(i))
       es = 1.d9
       a2 = 1.d9
       ovl = 1.d9
       cespm2 = 1.d9
       cesp = 1.d9
       potpt = 0.d0
       co = 1.d9
       al = 1.d9
       rad = 1.d9
       b_esp = 1.d9
       qesp = 1.d9
       cespm = 1.d9
       cc = 1.d9
       rnai1 = 0.d0


    end if
    end subroutine setup_esp
     subroutine getattrib(xmin, xmax)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
          use molkst_C, only : numat, natoms
          use common_arrays_C, only : geo, coord, labels, nat
          use esp_C, only : nesp, co
          use chanel_C, only : iw
    !
    !     ROUTINE TO CALCULATE WILLIAMS SURFACE
    !
    !***********************************************************************
    !-----------------------------------------------
    !   I n t e r f a c e   B l o c k s
    !-----------------------------------------------
          implicit none
        !.. Formal Arguments ..
          double precision, dimension (3), intent (out) :: xmin
          double precision, dimension (3), intent (out) :: xmax
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
          integer :: icntr, i, j, ia
          double precision, dimension(53) :: vderw
          double precision :: shell, vdmax
    !-----------------------------------------------
    !
          data vderw/ 53*0.0D0/
          vderw(1) = 2.4D0
          vderw(5) = 3.0D0
          vderw(6) = 2.9D0
          vderw(7) = 2.7D0
          vderw(8) = 2.6D0
          vderw(9) = 2.55D0
          vderw(15) = 3.1D0
          vderw(16) = 3.05D0
          vderw(17) = 3.0D0
          vderw(35) = 3.15D0
          vderw(53) = 3.35D0
          shell = 1.2D0
          nesp = 0
    !     CHECK IF VDERW IS DEFINED FOR ALL ATOMS
    !
    !     CONVERT INTERNAL TO CARTESIAN COORDINATES
    !
          call gmetry (geo, coord)
    !
    !     STRIP COORDINATES AND ATOM LABEL FOR DUMMIES (I.E. 99)
    !
          icntr = 0
          do i = 1, natoms
            co(:,i) = coord(:,i)
            if (labels(i) == 99) cycle
            icntr = icntr + 1
            nat(icntr) = labels(i)
          end do
          numat = icntr
    !
          do i = 1, numat
            j = nat(i)
            if (vderw(j) == 0.0D0) go to 40
          end do
          go to 50
       40 continue
          write (iw, *) 'VAN DER WAALS'' RADIUS NOT DEFINED FOR ATOM', i
          write (iw, *) 'IN WILLIAMS SURFACE ROUTINE PDGRID!'
          call mopend (&
          'VAN DER WAALS'' RADIUS NOT DEFINED FOR ATOM IN WILLIAMS SURFACE ROUTINE PDGRID!')
          return
    !     NOW CREATE LIMITS FOR A BOX
       50 continue
          xmin = 100000.0D0
          xmax = -100000.0D0
          do ia = 1, numat
            where (co(:,ia) - xmin < 0.D0)
              xmin = co(:,ia)
            end where
            where (co(:,ia) - xmax > 0.D0)
              xmax = co(:,ia)
            end where
          end do
    !     ADD (OR SUBTRACT) THE MAXIMUM VDERW PLUS SHELL
          vdmax = 0.0D0
          vdmax = dmax1(maxval(vderw),vdmax)
          xmin = xmin - vdmax - shell
          xmax = xmax + vdmax + shell
    end subroutine getattrib

    subroutine espplane(iplane, xmin, step, ngridpts2, ngridpts3)
      use esp_C, only : nesp, potpt
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      integer :: iplane, ngridpts2, ngridpts3
      double precision, dimension (3), intent (in) :: xmin, step
      !
      !.. Local Scalars ..
      integer :: ix, iy
      double precision :: ygrid, zgrid

      if (nesp == 0) then
          do iy = 1, ngridpts2
              ygrid = xmin(2) + step(2) * (iy - 1)
              do ix = 1, ngridpts3
                  nesp = nesp + 1
                  potpt(1, nesp) = xmin(1) + step(1) * (ix - 1)
                  potpt(2, nesp) = ygrid
              end do
          end do
      end if

      zgrid = xmin(3) + step(3) * (iplane - 1)
      do ix = 1, nesp
          potpt(3, ix) = zgrid
      end do
    end subroutine espplane
