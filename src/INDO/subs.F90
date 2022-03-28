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

! ____________________________________________________________________________
! _________________________GENERALLY USED ROUTINES ___________________________
! ____________________________________________________________________________

      subroutine upcas (s)

      implicit none
      integer ::          i
      character*(*)    s

!     **** converts string to upper case ****

      do i = 1, len(s)
        if (s(i:i) <= 'z' .and. s(i:i) >= 'a')&
     &          s(i:i) = char(ichar(s(i:i)) - 32)
      end do

      return
      end subroutine upcas

!     ***********************************************************************

      subroutine st (a, b, ihalf)
      use reimers_C, only : n, stwt, nr, nslwr, nsupr, nstr, istr

      implicit none
      double precision ::   a(n, n), b(n, n), bb, bbb
      integer ::            i, j, ihalf, i1, ib, ir, j1, jb

!     **** this does the symmetry - orbital transform on a matrix in ao basis ***
!     ****			b = T(transp) . a . T			    ***
!     **** ihalf = 1 lower, = 2 upper, triangle of matrices only		    ***

      do i = 1, n
        do j = i, n
          if (ihalf == 2) then
            b(i, j) = 0.d0
          else
            b(j, i) = 0.d0
          end if
        end do
      end do

      do ir = 1, nr
        do i = nslwr(ir), nsupr(ir)
          do j = i, nsupr(ir)
            bbb = 0.d0
            do i1 = 1, nstr(i)
              bb = 0.d0
              ib = abs(istr(i1, i))
              do j1 = 1, nstr(j)
                jb = abs(istr(j1, j))
                if (jb > ib .and. ihalf == 2 .or.&
     &              jb <= ib .and. ihalf /= 2) then
                  bb = bb + isign(1, istr(j1, j)) * a(ib, jb) * stwt(j)
                else
                  bb = bb + isign(1, istr(j1, j)) * a(jb, ib) * stwt(j)
                end if
              end do
              bbb = bbb + bb * isign(1, istr(i1, i)) * stwt(i)
            end do
            if (ihalf == 2) then
              b(i, j) = bbb
            else
              b(j, i) = bbb
            end if
          end do
        end do
      end do

      return
      end subroutine st

!     *************************************************************************
      subroutine stgamm (gamma, g1)
      use reimers_C, only : n, nham

      implicit none
      double precision ::  g1(n, n), gamma(n, n)
      integer ::           i, j

!     **** symmetry - trandforms the gamma (i, i|j, j) and (i, j|i, j) matrices ****

!     **** first copy gamma into g1 ****
      do i = 1, n
        do j = 1, n
          g1(i, j) = gamma(i, j)
        end do
      end do

!     **** transform integrals (i, j|i, j) in upper half ****
      if (nham == 2) call st (g1, gamma, 2)

!     **** transform integrals (i, i|j, j) in lower half ****
      call st (g1, gamma, 1)

      return
      end subroutine stgamm

!     *************************************************************************

      subroutine ao2mo1 (a, f, c, d, i1, i2, j1, j2, weight)

!     **** this adds C(transp) . F . C onto A over MO range with weighting ****
!     **** as specified by Edwards and Zerner, eqn 24			   ****
!     **** in MO basis, the proj operators simply limit the range of i, j   ****
      use reimers_C, only : n, nb2, matind

      implicit none
      double precision ::  a(nb2), f(nb2), c(n, n), d(nb2, 2), wk(n, n), weight, &
                        xx
      integer ::           i, j, k, i1, i2, j1, j2, nm

! RMG edit
      do i = 1, n
        wk(i, 1) = d(i, 1)
        do j = 2, n
          wk(i, j) = 0.d0
        end do
      end do
! End RMG

      do i = 1, n
        do j = j1, j2
          xx = 0.d0
          do k = 1, n
            nm = matind (max(i, k)) + min(i, k)
            xx = xx + f(nm) * c(j, k)
          end do
          wk(i, j) = xx
        end do
      end do

      do j = j1, j2
        do i = i1, min(j, i2)
          xx = 0.d0
          do k = 1, n
            xx = xx + c(i, k) * wk(k, j)
          end do
          nm = matind(max(i, j)) + min(i, j)
          a(nm) = a(nm) + xx * weight
        end do
      end do

! RMG edit
      do i = 1, n
         d(i, 1) = wk(i, 1)
      end do
! End RMG

      return
      end subroutine ao2mo1

!     *************************************************************************

      subroutine ao2mo (a, c, wk, mm, nn, n1, n2)

!     **** this transforms a matrix from AO to MO basis ****
!     **** or from configuration state basis to CI states ****
!     **** calculating output for states n1 to n2 only ****
      use reimers_C, only : matind

      implicit none
      integer ::           i, j, k, mm, nn, n1, n2, nm
      double precision ::  a(*), c(mm, nn), wk(mm, nn), xx

      do i = 1, nn
        do j = n1, n2
          xx = 0.d0
          do k = 1, nn
            nm = matind (max(i, k)) + min(i, k)
            xx = xx + a(nm) * c(j, k)
          end do
          wk(i, j) = xx
        end do
      end do

      do j = n1, n2
        nm = matind(j)
        do i = n1, j
          nm = nm + 1
          xx = 0.d0
          do k = 1, nn
            xx = xx + c(i, k) * wk(k, j)
          end do
          a(nm) = xx
        end do
      end do

      return
      end subroutine ao2mo

!     *************************************************************************

      subroutine mo2ao (a, c, d, nn)
      use reimers_C, only : nb2, matind

!     **** this transforms a matrix A from MO to AO basis set ****
!     **** or from CI state basis to spin - adapted configuration basis
      implicit none
      integer ::           i, j, k, nm, nn
      double precision ::  a(*), c(nn, nn), wk(nn, nn), d(nb2, 2), xx

! RMG edit
      do i = 1, nn
        wk(i, 1) = d(i, 1)
        do j = 2, nn
          wk(i, j) = 0.d0
        end do
      end do
! End RMG
      do i = 1, nn
        do j = 1, nn
          xx = 0.d0
          do k = 1, nn
            nm = matind(max(i, k)) + min(i, k)
            xx = xx + a(nm) * c(k, j)

          end do
          wk(i, j) = xx
        end do
      end do

      nm = 0
      do j = 1, nn
        do i = 1, j
          nm = nm + 1
          xx = 0.d0
          do k = 1, nn
            xx = xx + c(k, i) * wk(k, j)
          end do
          a(nm) = xx
        end do
      end do
! RMG edit
      do i = 1, nn
         d(i, 1) = wk(i, 1)
      end do
! End RMG
      return
      end subroutine mo2ao

!     *************************************************************************

      subroutine dump (a, nn, ibtype, name)
      use reimers_C, only : matind, nham, natm, iat, nbt, nci, iatsym, kind, &
          nptg, nr, nmrep
      USE chanel_C, only : iw

      implicit none
      double precision ::  a(*), amax
      integer ::           i, j, k, i0, i1, i2, ibtype, ip, ir, k1, k2, nm, nn, nnb2
      character*(*) ::     name
      character*21 ::      f1150
      data f1150        /'(1x,2i3,a2,a3,15f8.?)'/

!     **** dumps a symmetric matrix					    ****
!     **** ibtype = 1 atom by atom; = 2 orbital by orbital; = 3 state by state;   ****
!     **** = 4 state by state max 15, 	 = 0 is basic			    ****

      write (iw, *)
      write (iw, *) name

!     **** find largest elements ****
      nnb2 = nn*(nn + 1)/2
      amax = 1.D-10
      do i = 1, nnb2
        amax = max (amax, abs(a(i)))
      end do
      ip = int(log10(amax) + 0.99d0)
      ip = max (0, min (4, 6 - ip))

      k2 = 0
      do while (k2 < nn)
        k1 = k2 + 1
        k2 = min (k1 + 14, nn)
        write (iw, *)
        write (iw, '(10x, 15i8)') (i, i = k1, k2)
        if (ibtype < 3) then
          do j = k1, nn
            nm = matind(j)
            if (ibtype == 0) then
!	      **** generic matrix ****
              write (iw, '(1x, i6, 1x, a2, 2x, 15f8.4)') j, '  ', (a(nm + k), k = k1, min(j, k2))
            else if (ibtype == 1) then
!	      **** atom by atom matrix ****
              write (iw, '(1x, i6, 1x, a2, 2x, 15f8.4)') j, iatsym(natm(j)), (a(nm + k), k = k1, min(j, k2))
            else
!	      **** orbital by orbital matrix ****
              f1150(20:20) = char (ip + ichar('0'))
              i = iat(j)
              write (iw, f1150) j, i, iatsym(natm(i)), kind(nbt(j)), &
     &                                (a(nm + k), k = k1, min(j, k2))
              if (nham /= 1 .and. nham /= 3 .and. (j == nn .or.&
     &           i /= iat(j + 1)) ) write (iw, *)
            end if
          end do

        else
!	  **** configuration by configuration matrix, symmetry - block wize ****

          i0 = 1
          do ir = 1, nr
            if (nci(ir) > 0) then
              if (nptg > 1) write (iw, "(' Block of symmetry ', a4)") nmrep(ir, nptg)
              i1 = i0 - 1
              i2 = i1 + nci(ir)
              if (ibtype == 4) i2 = min (i2, i1 + 14)
              do i = i1 + 1, i1 + nci(ir)
                nm = matind(i) + i1
                write (iw, '(1x, i6, 5x, 15f8.4)') i, (a(nm + k), k = 1, min(15, i - i1))
              end do
              i0 = i0 + nci(ir)
            end if
          end do
          k2 = nn

        end if
      end do

      return
      end subroutine dump

!     *************************************************************************

      subroutine polzro

      use reimers_C, only : pol
      implicit none
      integer ::           j

      do 100 j = 1, 6
        pol(j) = 0.d0
100   continue

      return
      end subroutine polzro

!     *************************************************************************

      subroutine polizn (tmx, tmy, tmz, dele)
      use reimers_C, only : pol

      implicit none
      double precision ::  tmx, tmy, tmz, dele

      if (abs(dele) < 1.d1) then
        return
      end if
      pol(1) = pol(1) + tmx**2 / dele
      pol(2) = pol(2) + tmx*tmy/ dele
      pol(3) = pol(3) + tmy**2 / dele
      pol(4) = pol(4) + tmx*tmz/ dele
      pol(5) = pol(5) + tmy*tmz/ dele
      pol(6) = pol(6) + tmz**2 / dele

      return
      end subroutine polizn

!     *************************************************************************

      subroutine polout
      use reimers_C, only : au2ang, au2cm, pol
      USE chanel_C, only : iw

      implicit none
      double precision ::  conv, toang
      integer ::           j

      conv = 2.d0 * au2cm / au2ang**2
      toang = au2ang**3
      write (iw, *)
      write (iw, 800) 'au  ', (pol(j)*conv, j = 1, 6)
      write (iw, 800) 'A**3', (pol(j)*conv*toang, j = 1, 6)

      return
800   format (' Polarizability (', a4, ') xx=', f8.2, &
     & ' xy=', f8.2, ' yy=', f8.2, ' xz=', f8.2, ' yz=', f8.2, ' zz=', f8.2)
      end subroutine polout
