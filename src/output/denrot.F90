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

      subroutine denrot()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : norbs, numat, mpack, id, l1u, l2u, l3u, &
      mozyme
!
      use chanel_C, only : iw
      use elemts_C, only : elemnt
      use common_arrays_C, only : nat, geo, nfirst, nlast, coord, p, tvec
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(norbs) :: natom
      integer , dimension(5,35) :: irot
      integer , dimension(9) :: isp
      integer :: iprt, i, if, il, ipq, ii, i1, j1, j, mm, nn, kk, jf, jl, jpq, &
        jj, k, ll, l, l1, l2, ij, limit, na, m, ma, n
      double precision, dimension(9,9) :: arot
      double precision, dimension(3,5,5) :: c
      double precision, dimension(9,9) :: pab, vect
      double precision, dimension(:), allocatable :: b
      double precision :: delx, dely, delz, rmin, r, delxm, delym, delzm, sum
      character , dimension(21) :: line*6
      character , dimension(9) :: atorbs*7
      character, dimension(norbs) :: itext*7, jtext*2

      save irot, isp, atorbs
!-----------------------------------------------
!***********************************************************************
!
! DENROT PRINTS THE DENSITY MATRIX AS (S-SIGMA, P-SIGMA, P-PI) RATHER
!        THAN (S, PX, PY, PZ).
!
!***********************************************************************
      data atorbs/ 'S-SIGMA', 'P-SIGMA', '  P-PI ', '  P-PI ', 'D-SIGMA', &
        '  D-PI ', '  D-PI ', ' D-DELL', ' D-DELL'/
!**********************************************************************
! IROT IS A MAPPING LIST. FOR EACH ELEMENT OF AROT 5 NUMBERS ARE
! NEEDED. THESE ARE, IN ORDER, FIRST AND SECOND SUBSCRIPTS OF AROT,
! AND FIRST,SECOND, AND THIRD SUBSCRIPTS OF C, THUS THE FIRST
! LINE OF IROT DEFINES AROT(1,1)=C(1,3,3)
!
!**********************************************************************
      data irot/ 1, 1, 1, 3, 3, 2, 2, 2, 4, 3, 3, 2, 2, 2, 3, 4, 2, 2, 3, 3, 2&
        , 3, 2, 4, 2, 3, 3, 2, 2, 2, 4, 3, 2, 3, 2, 2, 4, 2, 4, 4, 3, 4, 2, 2, &
        4, 4, 4, 2, 3, 4, 5, 5, 3, 1, 5, 6, 5, 3, 4, 3, 7, 5, 3, 3, 3, 8, 5, 3&
        , 2, 3, 9, 5, 3, 5, 3, 5, 6, 3, 1, 2, 6, 6, 3, 4, 2, 7, 7, 3, 3, 2, 8, &
        6, 3, 2, 2, 9, 6, 3, 5, 2, 5, 7, 3, 1, 4, 6, 7, 3, 4, 4, 7, 7, 3, 3, 4&
        , 8, 7, 3, 2, 4, 9, 7, 3, 5, 4, 5, 8, 3, 1, 1, 6, 8, 3, 4, 1, 7, 8, 3, &
        3, 1, 8, 8, 3, 2, 1, 9, 8, 3, 5, 1, 5, 9, 3, 1, 5, 6, 9, 3, 4, 5, 7, 9&
        , 3, 3, 5, 8, 9, 3, 2, 5, 9, 9, 3, 5, 5/
      data isp/ 1, 2, 3, 3, 4, 5, 5, 6, 6/
      if (mozyme) then
        call denrot_for_MOZYME()
        return
      end if
      allocate(b(mpack))
      call gmetry (geo, coord)
      iprt = 0
      j1 = 0
      do i = 1, numat
        if = nfirst(i)
        il = nlast(i)
        ipq = il - if - 1
        ii = ipq + 2
        if (ii == 0) cycle
        do i1 = 1, ii
          j1 = iprt + isp(i1)
          itext(j1) = atorbs(i1)
          jtext(j1) = elemnt(nat(i))
          natom(j1) = i
        end do
        iprt = j1
        if (ipq /= 2) ipq = min(max(ipq,1),3)
        do j = 1, i
          delx = coord(1,j) - coord(1,i)
          dely = coord(2,j) - coord(2,i)
          delz = coord(3,j) - coord(3,i)
          if (id /= 0) then
            rmin = 100.D0
            delxm = 0.d0
            delym = 0.d0
            delzm = 0.d0
            do mm = -l1u, l1u
              do nn = -l2u, l2u
                do kk = -l3u, l3u
                  r = (tvec(1,1)*mm+tvec(1,2)*nn+tvec(1,3)*kk+delx)**2 + (tvec(&
                    2,1)*mm+tvec(2,2)*nn+tvec(2,3)*kk+dely)**2 + (tvec(3,1)*mm+&
                    tvec(3,2)*nn+tvec(3,3)*kk+delz)**2
                  if (r >= rmin) cycle
                  rmin = r
                  delxm = tvec(1,1)*mm + tvec(1,2)*nn + tvec(1,3)*kk + delx
                  delym = tvec(2,1)*mm + tvec(2,2)*nn + tvec(2,3)*kk + dely
                  delzm = tvec(3,1)*mm + tvec(3,2)*nn + tvec(3,3)*kk + delz
                end do
              end do
            end do
            delx = delxm
            dely = delym
            delz = delzm
          end if
          jf = nfirst(j)
          jl = nlast(j)
          jpq = jl - jf - 1
          jj = jpq + 2
          if (jj == 0) cycle
          if (jpq /= 2) jpq = min(max(jpq,1),3)
          pab = 0.D0
          kk = 0
          do k = if, il
            kk = kk + 1
            ll = 0
            pab(kk,:jl-jf+1) = p(jf + (k*(k-1))/2:jl + (k*(k-1))/2)
          end do
          call coe (delx, dely, delz, ipq, jpq, c, r)
          arot = 0.D0
          do i1 = 1, 35
            arot(irot(1,i1),irot(2,i1)) = c(irot(3,i1),irot(4,i1),irot(5,i1))
          end do
          l1 = isp(ii)
          l2 = isp(jj)
          vect = -1.D0
          j1 = 10
          j1 = 1
          if (l2 > 0) then
            vect(:l1,:l2) = 0.D0
            j1 = l2 + 1
          end if
          if (i /= j) then
            ij = max(ii,jj)
            do i1 = 1, ii
              do j1 = 1, jj
                sum = 0.D0
                do l1 = 1, ij
                  do l2 = 1, ij
                    sum = sum + arot(l1,i1)*pab(l1,l2)*arot(l2,j1)
                  end do
                end do
                vect(isp(i1),isp(j1)) = vect(isp(i1),isp(j1)) + sum**2
              end do
            end do
          end if
          k = 0
          do i1 = if, il
            k = k + 1
            l = 0
            do j1 = jf, jl
              l = l + 1
              if (j1 > i1) cycle
              b(j1+(i1*(i1-1))/2) = vect(k,l)
            end do
          end do
        end do
      end do
!
! NOW TO REMOVE ALL THE DEAD SPACE IN P, CHARACTERIZED BY -1.0
!
      l = 0
      do i = 1, mpack
        if (b(i) <= (-0.1D0)) cycle
        l = l + 1
        b(l) = b(i)
      end do
!
!   PUT ATOMIC ORBITAL VALENCIES ONTO THE DIAGONAL
!
      do i = 1, iprt
        sum = 0.D0
        ii = (i*(i - 1))/2
        do j = 1, i
          sum = sum + b(j+ii)
        end do
        do j = i + 1, iprt
          sum = sum + b((j*(j-1))/2+i)
        end do
        b((i*(i+1))/2) = sum
      end do
      line = '------'
      limit = (iprt*(iprt + 1))/2
      kk = 8
      na = 1
  190 continue
      ll = 0
      m = min0(iprt + 1 - na,6)
      ma = 2*m + 1
      m = na + m - 1
      write (iw, '(/16X,10(1X,A7,3X))') (itext(i),i=na,m)
      write (iw, '(15X,10(2X,A2,I3,4X))') (jtext(i),natom(i),i=na,m)
      write (iw, '(20A6)') (line(k),k=1,ma)
      do i = na, iprt
        ll = ll + 1
        k = (i*(i - 1))/2
        l = min0(k + m,k + i)
        k = k + na
        if (kk + ll > 50) then
          write (iw, '(/17X,10(1X,A7,3X))') (itext(n),n=na,m)
          write (iw, '( 17X,10(2X,A2,I3,4X))') (jtext(n),natom(n),n=na,m)
          write (iw, '(20A6)') (line(n),n=1,ma)
          kk = 4
          ll = 0
        end if
        write (iw, '(1X,A7,1X,A2,I3,10F11.6)') itext(i), jtext(i), natom(i), (b&
          (n),n=k,l)
      end do
      if (l >= limit) go to 220
      kk = kk + ll + 4
      na = m + 1
      if (kk + iprt + 1 - na <= 50) go to 190
      kk = 4
      go to 190
  220 continue
      deallocate (b)
      return
      end subroutine denrot
