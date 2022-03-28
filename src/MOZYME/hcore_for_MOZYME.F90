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

subroutine hcore_for_MOZYME ()
  use molkst_C, only: numat, norbs, id, numcal, n2elec, moperr, &
     & enuclr, l1u, l2u, l3u, keywrd, efield, mpack, cutofp
  use cosmo_C, only : useps,  phinet, qscnet, qdenet
  use linear_cosmo, only : addnucz
  use chanel_C, only: iw
  use funcon_C, only: a0, ev
  use overlaps_C, only : cutof1, cutof2
  use parameters_C, only: tore, dd, natorb
  use common_arrays_C, only: h, coord, nat, w, wj => w, wk, uspd, tvec
  use MOZYME_C, only : semidr, direct, cutofs, parth, &
    iorbs, jopt, mode, numred, refnuc
  implicit none
!
  character (len=248) :: tmpkey
  logical :: calci, calcij, calcj, fldon
  logical, save :: debug
  integer :: i1, i2, ii, im1, io1, ione, ired, j, j1, jj, jo1, jred, k, &
 & krmax, kro, ks, mm, nj, kr, i, ni, itemp, imol=0
  double precision :: const, enuc, fldcon, fnuc, half, hterme, xf, yf, zf, xj(3)
  double precision, parameter :: eps = 1.d-10
  double precision, dimension (45) :: e1b, e2a
  double precision, dimension (2025) :: wjd, wkd
  double precision, dimension (9, 9) :: di, dibits
  double precision, external :: reada
  integer, external :: ijbo
  fldcon = 0.d0
  fnuc = 0.d0
  debug = (Index (keywrd, " HCORE") /= 0)
  call add_more_interactions()
  if (moperr) return
  if (imol /= numcal) then
    if (index(keywrd, " SILENT") == 0 .and. (debug .or. id /= 0)) then
      write (iw, "(A,31X,F12.6,A)") " Overlap Cutoff Distance:", Sqrt(cutofs), " CUTOFS"
      write (iw, "(A,F12.6,A)") &
       & " Cutoff for quadrupolar and higher 2-electron integrals:", Sqrt(cutof2), " CUTOF2"
      write (iw, "(A,26X,F12.6,A)") " Cutoff for dipolar integrals:", Sqrt(cutof1), " CUTOF1"
      if (id /= 0) &
          write (iw, "(A,39X,F12.6,A)") " Madelung cutoff:", cutofp, " CUTOFP"
    end if
    imol = numcal
    xf = 0.d0
    yf = 0.d0
    zf = 0.d0
    tmpkey = trim(keywrd)
    i = Index (tmpkey, " FIELD(") + Index (tmpkey, " FIELD=(")
    if (i /= 0) then
       !
       !   ERASE ALL TEXT FROM TMPKEY EXCEPT FIELD DATA
       !
      tmpkey (:i) = " "
      itemp = Index (tmpkey, ")")
      tmpkey (itemp:) = " "
       !
       !   READ IN THE EFFECTIVE FIELD IN X,Y,Z COORDINATES
       !
      xf = reada (tmpkey, i)
      i = Index (tmpkey, ",")
      if (i /= 0) then
        tmpkey (i:i) = " "
        yf = reada (tmpkey, i)
        i = Index (tmpkey, ",")
        if (i /= 0) then
          tmpkey (i:i) = " "
          zf = reada (tmpkey, i)
        end if
      end if
      write (iw, "(/10X,'THE ELECTRIC FIELD IS',3F10.5, ' VOLTS/ANGSTROM',&
       &/)") xf, yf, zf
    end if
    !
    !        CONST = Ao/(8h)  (h=Hartree = eV/atomic unit)
    !
    const = a0 / ev
    !
    efield(1) = xf * const
    efield(2) = yf * const
    efield(3) = zf * const
  end if
  ione = 1
  if (id /= 0) ione = 0
  if (mode == -1) then
    h(:mpack) = -h(:mpack)
    enuclr = -enuclr
  else if (mode == 0) then
    enuclr = 0.d0
  else
    h = parth
    enuclr = refnuc
  end if
  krmax = n2elec + 101
  fldon = .false.
  if (Abs (efield(1)) > eps .or. Abs (efield(2)) > eps .or. Abs (efield(3)) > eps) then
    fldcon = eV/a0 ! = 51.42
    fldon = .true.
  end if
 !
  kr = 1
 !
  mm = 0
  ired = 1
  if (mode == 0) then
    !
    !   ZERO OUT THE H MATRIX
    !
    h(1:mpack) = 0.d0
  end if
  do i = 1, numat

    calci = (jopt(ired) == i)
    if (calci .and. ired < numred) then
      ired = ired + 1
    end if
    ni = nat(i)
    if (mode == 0) then
        !
        ! FILL THE DIAGONALS, AND OFF-DIAGONALS ON THE SAME ATOM
        !
      i2 = ijbo (i, i)
      do i1 = 1, iorbs(i)
        do j1 = 1, i1
          i2 = i2 + 1
          h(i2) = 0.d0
          if (fldon) then
            io1 = i1 - 1
            jo1 = j1 - 1
            if ((jo1 == 0) .and. (io1 == 1)) then
              hterme = -a0 * dd(ni) * efield(1) * fldcon
              h(i2) = hterme
            end if
            if ((jo1 == 0) .and. (io1 == 2)) then
              hterme = -a0 * dd(ni) * efield(2) * fldcon
              h(i2) = hterme
            end if
            if ((jo1 == 0) .and. (io1 == 3)) then
              hterme = -a0 * dd(ni) * efield(3) * fldcon
              h(i2) = hterme
            end if
          end if
        end do
        mm = mm + 1
        h(i2) = uspd(mm)
        if (fldon) then
          fnuc = -(efield(1)*coord(1, i) + efield(2)*coord(2, i) + &
               & efield(3)*coord(3, i)) * fldcon
          h(i2) = h(i2) + fnuc
        end if
      end do
    end if
    if (fldon) then
      enuclr = enuclr - fnuc * tore(nat(i))
    end if
    !
    !   FILL THE ATOM-OTHER ATOM ONE-ELECTRON MATRIX<PSI(LAMBDA)|PSI(SIGMA)>
    !
    jred = 1
    im1 = i - ione
    do j = 1, im1
      half = 1.d0
      if (i == 46 .and. j < 888) then
    enuclr = enuclr
   end if
      if (i == j) half = 0.5d0
      calcj = (jopt(jred) == j)
      if (calcj .and. jred < numred) jred = jred + 1
      calcij = (calci .or. calcj .or. mode == 0)
      nj = nat(j)
      if (id == 0) then
        !
        !   Molecular system
        !
        if (ijbo(i, j) >= 0) then
          if (calcij) then
            call h1elec (ni, nj, coord(1, i), coord(1, j), di)
            ii = ijbo (i, j)
            if (i == j) then
              do i1 = 1, iorbs(i)
                do j1 = 1, i1
                  ii = ii + 1
                  h(ii) = h(ii) + di(i1, j1)
                end do
              end do
            else
              do i1 = 1, iorbs(i)
                do j1 = 1, iorbs(j)
                  ii = ii + 1
                  h(ii) = h(ii) + di(i1, j1)
                end do
              end do
            end if
            !
            ! CALCULATE THE TWO-ELECTRON INTEGRALS, W; THE ELECTRON
            ! NUCLEAR TERMS E1B AND E2A; AND THE NUCLEAR-NUCLEAR TERM ENUC.
            !
            call rotate (ni, nj, coord(1, i), coord(1, j), w(kr), i1, e1b, e2a, enuc)
            enuclr = enuclr + enuc
           !
          else if ( .not. direct) then
            kr = kr + (natorb(ni)*(natorb(ni)+1))/2 * (natorb(nj)*(natorb(nj)+1))/2
          end if
        else if (ijbo(i, j) ==-2) then
          if (calcij) then
            call outer2 (ni, nj, coord(1, i), coord(1, j), w(kr), kr, e1b, e2a, enuc, id, semidr)
            enuclr = enuclr + enuc
          else if ( .not. semidr) then
            if (natorb(ni)*natorb(nj) > 0) then
              if (natorb(ni) > 1) then
                if (natorb(nj) > 1) then
                  kr = kr + 7
                else
                  kr = kr + 4
                end if
              else if (natorb(nj) > 1) then
                kr = kr + 4
              else
                kr = kr + 1
              end if
            end if
          end if
        else if (calcij) then
          call outer1 (ni, nj, coord(1, i), coord(1, j), w(kr), kr, e1b, e2a, enuc, 0, semidr)
          enuclr = enuclr + enuc
        else if ( .not. semidr) then
          if (natorb(ni)*natorb(nj) > 0) kr = kr + 1
        end if
        if (calcij) then
          !
          !   ADD ON THE ELECTRON-NUCLEAR ATTRACTION TERM FOR ATOM I.
          !
          ii = ijbo (i, i)
          j1 = (iorbs(i)*(iorbs(i)+1)) / 2
          do i1 = 1, j1
            ii = ii + 1
            h(ii) = h(ii) + e1b(i1) * half
          end do
          !
          !   ADD ON THE ELECTRON-NUCLEAR ATTRACTION TERM FOR ATOM J.
          !
          ii = ijbo (j, j)
          j1 = (iorbs(j)*(iorbs(j)+1)) / 2
          do i1 = 1, j1
            ii = ii + 1
            h(ii) = h(ii) + e2a(i1) * half
          end do
        end if
        if (kr > krmax-100 .and. i /= numat) then
          write (iw,*) kr, krmax
          write (iw, "(' Running out of storage for W in HCORE ')")
          write (iw,*) " NUMBER OF ATOMS CALCULATED FOR 'W':", i
          write (iw,*) " NUMBER OF ATOMS IN SYSTEM:", numat
          call mopend ("Running out of storage for W in HCORE")
          return
        end if
      else
        !
        !   Solid-state system
        !
        if (ijbo(i, j) >= 0) then
          if (calcij) then
            di = 0.D0
            do ii = -l1u, l1u
              do jj = -l2u, l2u
                do k = -l3u, l3u
                  xj = coord(:,j) + tvec(:,1)*ii + tvec(:,2)*jj + tvec(:,3)*k
                  call h1elec (ni, nj, coord(1,i), xj, dibits)
                  di = di + dibits
                end do
              end do
            end do
            ii = ijbo (i, j)
            if (i == j) then
              do i1 = 1, iorbs(i)
                do j1 = 1, i1
                  ii = ii + 1
                  h(ii) = h(ii) + di(i1, j1)
                end do
              end do
            else
              do i1 = 1, iorbs(i)
                do j1 = 1, iorbs(j)
                  ii = ii + 1
                  h(ii) = h(ii) + di(i1, j1)
                end do
              end do
            end if
            !
            ! CALCULATE THE TWO-ELECTRON INTEGRALS, W;
            ! THE ELECTRON NUCLEAR TERMS E1B AND E2A;
            ! AND THE NUCLEAR-NUCLEAR TERM ENUC.
            !
            kro = kr
            call solrot (ni, nj, coord(1,i), coord(1,j), wjd, wkd, kr, e1b, e2a, enuc)
            w(kro:kr - 1) = wjd(:kr-kro)
            wk(kro:kr - 1) = wkd(:kr-kro)
            enuclr = enuclr + enuc
          else if (natorb(ni) == 1) then
            if (natorb(nj) == 1) then
              kr = kr + 1
            else
              kr = kr + 10
            end if
          else if (natorb(nj) == 1) then
            kr = kr + 10
          else
            kr = kr + 100
          end if
        else if (ijbo(i, j) ==-2) then
          if (calcij) then
            ks = kr
            call outer2 (ni, nj, coord(1, i), coord(1, j), wj(kr), kr, &
                 & e1b, e2a, enuc, id, semidr)
            do k = ks, kr - 1
              wk(k) = 0.d0
            end do
            enuclr = enuclr + enuc
          else if (natorb(ni)*natorb(nj) > 0) then
            if (natorb(ni) > 1) then
              if (natorb(nj) > 1) then
                kr = kr + 7
              else
                kr = kr + 4
              end if
            else if (natorb(nj) > 1) then
              kr = kr + 4
            else
              kr = kr + 1
            end if
          end if
        else if (calcij) then
          wk(kr) = 0.d0
          call outer1 (ni, nj, coord(1, i), coord(1, j), wj(kr), kr, &
               & e1b, e2a, enuc, id, semidr)
          enuclr = enuclr + enuc
        else if (natorb(ni)*natorb(nj) > 0) then
          kr = kr + 1
        end if
        if (calcij) then
           !
           !   ADD ON THE ELECTRON-NUCLEAR ATTRACTION TERM FOR ATOM I.
           !
          ii = ijbo (i, i)
          j1 = (iorbs(i)*(iorbs(i)+1)) / 2
          do i1 = 1, j1
            ii = ii + 1
            h(ii) = h(ii) + e1b(i1) * half
          end do
           !
           !   ADD ON THE ELECTRON-NUCLEAR ATTRACTION TERM FOR ATOM J.
           !
          ii = ijbo (j, j)
          j1 = (iorbs(j)*(iorbs(j)+1)) / 2
          do i1 = 1, j1
            ii = ii + 1
            h(ii) = h(ii) + e2a(i1) * half
          end do
        end if
        if (kr > krmax-100 .and. i /= numat) then
          write (iw,*) kr, krmax
          write (iw, "(' Running out of storage for W in HCORE ')")
          write (iw,*) " NUMBER OF ATOMS CALCULATED FOR 'W':", i
          write (iw,*) " NUMBER OF ATOMS IN SYSTEM:", numat
          call mopend ("Running out of storage for W in HCORE")
        end if
      end if
    end do
    ii = iorbs(i)
    ii = (ii*(ii+1)) / 2
    if (id /= 0) then
      do i1 = kr, kr + ii * ii - 1
        wk(i1) = 0.d0
      end do
    end if
    if (ii /= 0) then
      call wstore (w(kr), kr, ni, ii)
    end if
  end do
 !
 !
  if (mode == -1) then
    parth(:mpack) = -h(:mpack)
    refnuc = -enuclr
  end if
  if (useps) then
    ! In the following routine the dielectric correction to the core-core-
    ! interaction is added to ENUCLR (just set arrays to zero)

    call addnucz (phinet, qscnet, qdenet)

  end if
  kr = kr - 1
 !
 !
  if (debug) then
    write (iw, "(//10X,'ONE-ELECTRON MATRIX FROM HCORE')")
    if (mode ==-1) then
      write (iw, "(10X,A)") " AFTER REMOVAL OF TERMS FOR MOVING ATOMS"
    end if
    if (mode == 1) then
      write (iw, "(10X,A)") " AFTER ADDITION OF TERMS FOR MOVING ATOMS"
    end if
    call vecprt_for_MOZYME (h, norbs)
    if (kr > 2000) then
      write (iw,*) " THE TWO-ELECTRON MATRIX IS TOO LARGE TO PRINT"
    end if
    j = Min (kr, 2000)
    if (id == 0) then
      write (iw, "(//10X,'TWO-ELECTRON MATRIX IN HCORE'/)")
10000   format (10 f8.4)
      write (iw, 10000) (w(i), i=1, j)
    else
      write (iw, "(//10X,'TWO-ELECTRON J MATRIX IN HCORE'/)")
      write (iw, 10000) (wj(i), i=1, j)
      write (iw, "(//10X,'TWO-ELECTRON K MATRIX IN HCORE'/)")
      write (iw, 10000) (wk(i), i=1, j)
    end if
  end if
end subroutine hcore_for_MOZYME
