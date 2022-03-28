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

      subroutine ijkl(cp, cf, nelec, nmos, dijkl, cij, ckl, wcij, xy)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE molkst_C, only : norbs, numat, lm61
      USE cosmo_C, only : iseps
      use common_arrays_C, only : nfirst, nlast, w
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nelec
      integer , intent(in) :: nmos
      double precision , intent(in) :: cp(norbs,nmos)
      double precision , intent(in) :: cf(norbs,norbs)
      double precision , intent(inout) :: dijkl(norbs,nmos,(nmos*(nmos + 1))/2)
      double precision  :: cij(lm61)
      double precision  :: ckl(lm61)
      double precision  :: wcij(lm61)
      double precision , intent(out) :: xy(nmos,nmos,nmos,nmos)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ij, i, j, ipq, ii, ip, k, l, kk
      double precision :: sum
!-----------------------------------------------
!***********************************************************************
!
!   IJKL FILLS THE XY ARRAY.  XY HOLDS THE TWO-ELECTRON INTEGRALS OVER
!        MOLECULAR ORBITALS IN THE ACTIVE SPACE.
!        XY(I,J,K,L) = <IJ|1/R(1,2)|KL>
!
!           ON INPUT
!
! CP     = M.O.'S OVER C.I. ACTIVE SPACE (NORMALLY 1 TO 5 M.O.S)
! CF     = ALL M.O.'S, INCLUDING THOSE IN CP
! NORBS  = NUMBER OF ATOMIC ORBITALS
! NELEC  = NUMBER OF OCCUPIED M.O.S NOT INVOLVED IN THE C.I.
! NMOS   = NUMBER OF M.O.S INVOLVED IN THE C.I. (NORMALLY 1 TO 5 M.O.S)
!          ALSO CALLED THE ACTIVE SPACE OF THE C.I.
!
!  NOTE: THIS SUBROUTINE IS UNUSUAL IN THAT ONE FUNCTION IS TO
!        FILL THE ARRAY XY WHICH IS NOT PASSED AS AN ARGUMENT,
!        INSTEAD IT IS PASSED VIA COMMON BLOCK XYIJKL.
!
!***********************************************************************
!
!  CALCULATE TWO-ELECTRON INTEGRALS FOR THE SET DIJKL(K,L,IJ)
!  THE INDEX K RUNS OVER ALL M.O.'S, L OVER ACTIVE-SPACE M.O.'S,
!  AND IJ OVER LOWER-HALF TRIANGLE OF ACTIVE-SPACE M.O.'S, J FASTER THAN
!  I.
!  ALL ACTIVE-SPACE INTERACTIONS ARE COPIED INTO THE ARRAY XY
!
      ij = 0
      do i = 1, nmos
        do j = 1, i
          ij = ij + 1
          ipq = 0
          do ii = 1, numat
            do ip = nfirst(ii), nlast(ii)
              if (ip - nfirst(ii) + 1 > 0) then
                cij(ipq+1:ip-nfirst(ii)+1+ipq) = cp(ip,i)*cp(nfirst(ii):ip,j)&
                   + cp(ip,j)*cp(nfirst(ii):ip,i)
                ipq = ip - nfirst(ii) + 1 + ipq
              end if
            end do
          end do
!
!  CIJ HOLDS THE DENSITY DISTRIBUTION PSI(I)*PSI(J) OVER ATOMIC BASES
!  I AND J ARE M.O. INDICES WITHIN THE ACTIVE SPACE.  CIJ(M,N) IS FOR
!  THE ATOMIC BASES M AND N FOR M.O.'S I AND J.
!
          call partxy (cij, wcij, w)
!
! WCIJ HOLDS THE KET PART OF THE INTEGRAL <K,L|1/R(1,2)|I,J>
! THAT IS, |1/R(1,2)|I,J>.  WCIJ(M,N) IS FOR THE ATOMIC BASES M AND N
! FOR M.O.'S K AND L.
!
          do k = 1, norbs
            do l = 1, nmos
!
!  ABOUT TO CALCULATE <I,J|1/R(1,2)|K,L>
!
              ipq = 0
              do ii = 1, numat
                do ip = nfirst(ii), nlast(ii)
                  if (ip - nfirst(ii) + 1 > 0) then
                    ckl(ipq+1:ip-nfirst(ii)+1+ipq) = cf(ip,k)*cp(nfirst(ii):ip,&
                      l) + cp(ip,l)*cf(nfirst(ii):ip,k)
                    ipq = ip - nfirst(ii) + 1 + ipq
                  end if
                end do
              end do
!
! CKL HOLDS THE DENSITY DISTRIBUTION PSI(K)*PSI(L) OVER ATOMIC BASES.
! K IS THE INDEX OF A M.O.; L IS AN INDEX OF A M.O. IN THE ACTIVE SPACE.
!
              sum = 0.D0
              do ii = 1, ipq
                sum = sum + ckl(ii)*wcij(ii)
              end do
!
!  SUM IS THE INTEGRAL <I,J|1/R(1,2)|K,L>
!
              dijkl(k,l,ij) = sum
            end do
          end do
        end do
      end do
!
!  NOW SPREAD THE INTEGRALS OVER THE XY ARRAY.  XY IS ENTIRELY
!  IN ACTIVE SPACE
!
      do k = 1, nmos
        kk = nelec + k
!
!  K IS A M.O. INDEX IN ACTIVE SPACE
! KK IS A M.O. INDEX
!
        do l = 1, nmos
          ij = 0
          do i = 1, nmos
            do j = 1, i
              ij = ij + 1
              sum = dijkl(kk,l,ij)
              xy(i,j,k,l) = sum
              xy(i,j,l,k) = sum
              xy(j,i,k,l) = sum
              xy(j,i,l,k) = sum
              xy(k,l,i,j) = sum
              xy(k,l,j,i) = sum
              xy(l,k,i,j) = sum
            end do
          end do
        end do
      end do
      if (.not.iseps) return
      do i = 1, nmos
        ipq = 0
        do ii = 1, numat
          do ip = nfirst(ii), nlast(ii)
            if (ip - nfirst(ii) + 1 > 0) then
              cij(ipq+1:ip-nfirst(ii)+1+ipq) = cp(ip,i)*cp(nfirst(ii):ip,i)
              ipq = ip - nfirst(ii) + 1 + ipq
            end if
          end do
        end do
        call ciint (cij, wcij)
        do l = 1, i
          ipq = 0
          do ii = 1, numat
            do ip = nfirst(ii), nlast(ii)
              if (ip - nfirst(ii) + 1 > 0) then
                ckl(ipq+1:ip-nfirst(ii)+1+ipq) = cp(ip,l)*cp(nfirst(ii):ip,l)
                ipq = ip - nfirst(ii) + 1 + ipq
              end if
            end do
          end do
        end do
      end do
      return
      end subroutine ijkl
