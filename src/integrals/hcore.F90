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

      subroutine hcore()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE molkst_C, only : numcal, numat, norbs, id, l1u, l2u, l3u, &
      & keywrd, enuclr, n2elec, efield, mpack, line
      use common_arrays_C, only : nfirst, nlast, nat, uspd, &
      & coord, h, w, wk, tvec
      use cosmo_C, only : useps
      USE funcon_C, only : a0, ev, fpc_9
      USE parameters_C, only : tore, dd
      use MOZYME_C, only : cutofs
      use overlaps_C, only : cutof1, cutof2
      USE chanel_C, only : iw

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icalcn, ione
      integer :: i, j, k, kr, ia, ib, ni, i1, i2, j1, io1, jo1, im1, &
        ja, jb, nj, ii, jj, kro, q_at, link_at, mol_at
      double precision, dimension(45) :: e1b(45), e2a(45), di(9,9), wjd(2025), wkd(2025), xj(3), dibits(9,9)
      double precision, allocatable:: vqc(:)
      double precision ::  xf, yf, zf, const, fldcon, hterme, fnuc, half, enuc
      logical :: fldon, first, debug, lmolaris_qmmm, exists
      character :: tmpkey*241
      double precision, external :: reada
      save first, debug, icalcn, ione, vqc, lmolaris_qmmm
!***********************************************************************
!
!   HCORE GENERATES THE ONE-ELECTRON MATRIX AND TWO ELECTRON INTEGRALS
!         FOR A GIVEN MOLECULE WHOSE GEOMETRY IS DEFINED IN CARTESIAN
!         COORDINATES.
!
!  ON INPUT  COORD   = COORDINATES OF THE MOLECULE.
!
!  ON OUTPUT  H      = ONE-ELECTRON MATRIX.
!             W      = TWO-ELECTRON INTEGRALS.
!             ENUCLR = NUCLEAR ENERGY
!***********************************************************************
      data icalcn/ 0/
      fnuc = 0.d0
      fldcon = 0.d0
      first = icalcn /= numcal
      icalcn = numcal
      if (first) then
        ione = 1
        cutofs = 1.d4
        cutof2 = 1.D10
        cutof1 = 225.d0 ! (Cutoff distance)**2 for overlap integrals
        if (id /= 0) then
          if (allocated(wk)) deallocate (wk)
          allocate(wk(n2elec))
        end if
        if (id /= 0) ione = 0
        debug = (index(keywrd,'HCORE') /= 0)
        lmolaris_qmmm = (index(keywrd,'QMMM') /= 0)
        if (lmolaris_qmmm) then
!
!  Read in the energy in kcal/mol that an electron on each atom would have, arising from
!  the partial charges on all atoms.
!
          line="mol.in"
          call add_path(line)
          inquire (file=trim(line), exist = exists)
          if (.not. exists) then
             call mopend("A file named '"//trim(line)//"' was expected, but was not found.")
             return
          end if
          open(85, file=line, status='unknown', form='formatted')
          read(85, *, iostat = i)
          if (i /= 0) then
             call mopend("A file named '"//trim(line)//"' exists, but is empty.")
             return
          end if
          read(85, *, iostat = i)q_at, link_at
          if (i /= 0) then
             call mopend("A file named '"//trim(line)//"' exists, but is corrupt.")
             return
          end if
          mol_at = q_at + link_at
          if(mol_at /= numat) then
            write(line,'(a,i5,a,i5,a)') &
              " The number of atoms in 'mol.in':",mol_at," and in the MOPAC data set:",numat," are different."
            call mopend(trim(line))
            write(iw,'(10x,a)')" Correct fault and resubmit."
            return
          end if
          if (allocated(vqc)) deallocate(vqc)
          allocate(vqc(numat))
          do i = 1, numat
            read(85,*)line,line,line,line,vqc(i)
            if (debug) write(iw,'(a, i4, a, f9.3)')"ATOM No.",i," VQC(I)",vqc(i)
          end do
          close(85)
        end if
        xf = 0.D0
        yf = 0.D0
        zf = 0.D0
        tmpkey = trim(keywrd)
        i = index(tmpkey,' FIELD(') + index(tmpkey,' FIELD=(')
        if (i /= 0) then
!
!   ERASE ALL TEXT FROM TMPKEY EXCEPT FIELD DATA
!
          tmpkey(:i) = ' '
          tmpkey(index(tmpkey,')'):) = ' '
!
!   READ IN THE EFFECTIVE FIELD IN X,Y,Z COORDINATES
!
          xf = reada(tmpkey,i)
          i = index(tmpkey,',')
          if (i /= 0) then
            tmpkey(i:i) = ' '
            yf = reada(tmpkey,i)
            i = index(tmpkey,',')
            if (i /= 0) then
              tmpkey(i:i) = ' '
              zf = reada(tmpkey,i)
            end if
          end if
          write (iw, '(/10X,''THE ELECTRIC FIELD IS'',3F10.5,'' VOLTS/ANGSTROM'',/)') xf, yf, zf
        end if
        const = a0/ev
!
        efield(1) = xf*const
        efield(2) = yf*const
        efield(3) = zf*const
      end if
      fldon = .FALSE.
      if (efield(1) /= 0.0D00 .or. efield(2) /= 0.0D00 .or. efield(3) /= 0.0D00) then
!
!   FLDCON = h/Ao
!
        fldcon = ev/a0
        fldon = .TRUE.
      end if
      enuclr = 0.d0
      h(:mpack) = 0.d0
      kr = 1
      do i = 1, numat
        ia = nfirst(i)
        ib = nlast(i)
        ni = nat(i)
!
! FIRST WE FILL THE DIAGONALS, AND OFF-DIAGONALS ON THE SAME ATOM
!
        if (.not.fldon) then
          do i1 = ia, ib
            i2 = i1*(i1 - 1)/2 + ia - 1
            if (i1 - ia + 1 > 0) then
              h(i2+1:i1-ia+1+i2) = 0.D0
              i2 = i1 - ia + 1 + i2
            end if
            h(i2) = uspd(i1)
            if(lmolaris_qmmm) then
              if (debug) write(iw,'(''OLD 1e MATRIX ELEMENT '',i5,f12.5,'' i'',i5,&
            & '' vqc(i)'',f12.5)')i2,h(i2),i,-vqc(i)/fpc_9
              h(i2) = h(i2) - vqc(i)/fpc_9
              if (debug) write(iw,'(''UPD 1e MATRIX ELEMENT '',i5,f12.5)')i2,h(i2)
            end if
            cycle
          end do
        else
          do i1 = ia, ib
            i2 = i1*(i1 - 1)/2 + ia - 1
            do j1 = ia, i1
              i2 = i2 + 1
              h(i2) = 0.D0
              io1 = i1 - ia
              jo1 = j1 - ia
              if (jo1==0 .and. io1==1) then
                hterme = -a0*dd(ni)*efield(1)*fldcon
                h(i2) = hterme
              end if
              if (jo1==0 .and. io1==2) then
                hterme = -a0*dd(ni)*efield(2)*fldcon
                h(i2) = hterme
              end if
              if (jo1/=0 .or. io1/=3) cycle
              hterme = -a0*dd(ni)*efield(3)*fldcon
              h(i2) = hterme
            end do
            h(i2) = uspd(i1)
            fnuc = -(efield(1)*coord(1,i) + efield(2)*coord(2,i) + efield(3)*coord(3,i))*fldcon
            h(i2) = h(i2) + fnuc
          end do
        end if
        if (fldon) enuclr = enuclr - fnuc*tore(nat(i))
        if (lmolaris_qmmm) enuclr = enuclr + vqc(i)/fpc_9*tore(nat(i))
!
!   FILL THE ATOM-OTHER ATOM ONE-ELECTRON MATRIX<PSI(LAMBDA)|PSI(SIGMA)>
!
        im1 = i - ione
        do j = 1, im1
          half = 1.D0
          if (i == j) half = 0.5D0
          ja = nfirst(j)
          jb = nlast(j)
          nj = nat(j)
          if (id == 0) then
            call h1elec (ni, nj, coord(1,i), coord(1,j), di)
          else
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
          end if
          i2 = 0
          do i1 = ia, ib
            ii = i1*(i1 - 1)/2 + ja - 1
            i2 = i2 + 1
            jj = min(i1,jb)
            h(ii+1:jj-ja+1+ii) = h(ii+1:jj-ja+1+ii) + di(i2,:jj-ja+1)
          end do
!
!   CALCULATE THE TWO-ELECTRON INTEGRALS, W; THE ELECTRON NUCLEAR TERMS
!   E1B AND E2A; AND THE NUCLEAR-NUCLEAR TERM ENUC.
!
          if (id == 0) then
            call rotate (ni, nj, coord(1,i), coord(1,j), w(kr), kr, e1b, e2a, enuc)
          else
            kro = kr
            call solrot (ni, nj, coord(1,i), coord(1,j), wjd, wkd, kr, e1b, e2a, enuc)
            jj = 0
            w(kro:kr - 1) = wjd(:kr-kro)
            wk(kro:kr - 1) = wkd(:kr-kro)
          end if
          enuclr = enuclr + enuc
!
!   ADD ON THE ELECTRON-NUCLEAR ATTRACTION TERM FOR ATOM I.
!
          i2 = 0
          do i1 = ia, ib
            ii = i1*(i1 - 1)/2 + ia - 1
            if (i1 - ia + 1 > 0) then
              h(ii+1:i1-ia+1+ii) = h(ii+1:i1-ia+1+ii) + e1b(i2+1:i1-ia+1+i2)*half
              i2 = i1 - ia + 1 + i2
            end if
          end do
!
!   ADD ON THE ELECTRON-NUCLEAR ATTRACTION TERM FOR ATOM J.
!
          i2 = 0
          do i1 = ja, jb
            ii = i1*(i1 - 1)/2 + ja - 1
            if (i1 - ja + 1 > 0) then
              h(ii+1:i1-ja+1+ii) = h(ii+1:i1-ja+1+ii) + e2a(i2+1:i1-ja+1+i2)*half
              i2 = i1 - ja + 1 + i2
            end if
          end do
        end do
        ii = ib - ia + 1
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
      if (useps) then
! In the following routine the dielectric correction to the core-core-
! interaction is added to ENUCLR
        call addnuc ()
! The following routine adds the dielectric correction for the electron-
! interaction to the diagonal elements of H
        call addhcr ()
      end if
! end of COSMO change
      if (debug) then
        write (iw, '(2/10X,''ONE-ELECTRON MATRIX FROM HCORE'')')
        call vecprt (h, norbs)
        j = min(400,kr - 1)
        j = kr - 1
        if (id == 0) then
          write (iw, '(2/10X,''TWO-ELECTRON MATRIX IN HCORE''/)')
          write (iw, 200) (w(i),i=1,j)
        else
          write (iw, '(2/10X,''TWO-ELECTRON J MATRIX IN HCORE''/)')
          write (iw, 200) (w(i),i=1,j)
          write (iw, '(2/10X,''TWO-ELECTRON K MATRIX IN HCORE''/)')
          write (iw, 200) (wk(i),i=1,j)
        end if
  200   format(10f8.4)
      end if
      return
      end subroutine hcore
