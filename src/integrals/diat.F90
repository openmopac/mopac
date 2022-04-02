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

      subroutine diat(ni, nj, xj, di)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use molkst_C, only : numcal
      use funcon_C, only : a0
      use overlaps_C, only : cutof1
      use parameters_C, only : natorb, zs, zp, zd, npq
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)  :: ni
      integer, intent(in)  :: nj
      double precision, intent(in)  :: xj(3)
      double precision, intent(out)  :: di(9,9)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: a, pq2, b, pq1
      integer , dimension(3,5) :: ival
      integer :: icalcn
      integer :: i, j, ia, ib, newk, nk1, iss, jss, k, kss, kmin, kmax, lmin, &
        lmax, l, ii, jj, pi, pj
      double precision, dimension(3,3,3) :: s
      double precision, dimension(3) :: ul1, ul2
      double precision, dimension(3,5,5) :: c
      double precision, dimension(27) :: slin
      double precision, dimension(3,5) :: c1, c2, c3, c4, c5
      double precision, dimension(3,3) :: s1, s2, s3
      double precision :: x2, y2, z2, r, aa, bb
      double precision, external :: ss
      logical :: use_diat2(107)

      save ival, icalcn, use_diat2
!-----------------------------------------------
!***********************************************************************
!
!   DIAT CALCULATES THE DI-ATOMIC OVERLAP INTEGRALS BETWEEN ATOMS
!        OF ATOMIC NUMBER NI AND NJ, WHERE NJ IS AT POSITION
!        XJ RELATIVE TO NI.
!
!   ON INPUT NI  = ATOMIC NUMBER OF THE FIRST ATOM.
!            NJ  = ATOMIC NUMBER OF THE SECOND ATOM.
!            XJ  = CARTESIAN COORDINATES OF THE SECOND ATOM,
!                  RELATIVE TO THE FIRST ATOM.
!
!  ON OUTPUT DI  = DIATOMIC OVERLAP, IN A 9 * 9 MATRIX. LAYOUT OF
!                  ATOMIC ORBITALS IN DI IS
!                  1   2   3   4   5            6     7       8     9
!                  S   PX  PY  PZ  D(X**2-Y**2) D(XZ) D(Z**2) D(YZ)D(XY)
!
!   LIMITATIONS:  IN THIS FORMULATION, NI AND NJ MUST BE LESS THAN 107
!         EXPONENTS ARE ASSUMED TO BE PRESENT IN COMMON BLOCK EXPONT.
!
!***********************************************************************
      equivalence (slin(1), s(1,1,1))
      equivalence (c1(1,1), c(1,1,1)), (c2(1,1), c(1,1,2)), (c3(1,1), c(1,1,3)), &
        (c4(1,1), c(1,1,4)), (c5(1,1), c(1,1,5)), (s1(1,1), s(1,1,1)), (s2(1,1), &
        s(1,1,2)), (s3(1,1), s(1,1,3))
      data ival/ 1, 0, 9, 1, 3, 8, 1, 4, 7, 1, 2, 6, 0, 0, 5/
      data icalcn/ 0/
      if (icalcn /= numcal) then
        icalcn = numcal
        use_diat2 = .false.
        do i = 1, 17
          use_diat2(i) = (natorb(i) < 5)
        end do
        use_diat2(2)  = .false.
        use_diat2(10) = .false.
      end if
      x2 = xj(1)
      y2 = xj(2)
      z2 = xj(3)
      pq1 = npq(ni,1)
      pq2 = npq(nj,1)
      di = 0.0D0
      r = x2**2 + y2**2 + z2**2
      if (pq1==0 .or. pq2==0 .or. r>=cutof1) return
      if (natorb(ni)==0 .or. natorb(nj)==0) return
      call coe (x2, y2, z2, natorb(ni), natorb(nj), c, r)
      if (r < 0.001D0) return
      ia = min(pq1 + 1,3)
      ib = min(pq2 + 1,3)
      a = ia - 1
      b = ib - 1
      if (use_diat2(ni) .and. use_diat2(nj)) then
        call diat2 (ni, zs(ni), zp(ni), r, nj, zs(nj), zp(nj), s, a0)
      else
        ul1(1) = zs(ni)
        ul2(1) = zs(nj)
        ul1(2) = zp(ni)
        ul2(2) = zp(nj)
        ul1(3) = max(zd(ni),0.3D0)
        ul2(3) = max(zd(nj),0.3D0)
        slin = 0.0D0
        newk = min(a,b)
        nk1 = newk + 1
        do i = 1, ia
          iss = i
          ib = b + 1
          pq1 = npq(ni,i)
          do j = 1, ib
            jss = j
            pq2 = npq(nj,j)
            do k = 1, nk1
              if (k>i .or. k>j) cycle
              kss = k
              pi = max(pq1,iss)
              pj = max(pq2,jss)
              s(i,j,k) = ss(pi,pj,iss,jss,kss,ul1(i),ul2(j),r,a0)
            end do
          end do
        end do
      end if
      do i = 1, ia ! i goes over  L = s, p, d for atom a
        kmin = 4 - i
        kmax = 2 + i
        do j = 1, ib ! j goes over  L = s, p, d for atom b
          if (j == 2) then
            aa = -1.D0
            bb = 1.D0
          else
            aa = 1.D0
            if (j == 3) then
              bb = -1.D0
            else
              bb = 1.D0
            end if
          end if
          lmin = 4 - j
          lmax = 2 + j
          do k = kmin, kmax   ! k goes over M_L = -del, -pi, sigma, pi, del for atom a
            do l = lmin, lmax ! l goes over M_L = -del, -pi, sigma, pi, del for atom b
              ii = ival(i,k)  !                     1    2     3     4   5
              jj = ival(j,l)
              di(ii,jj) = s1(i,j)*(c3(i,k)*c3(j,l))*aa + &
                          s2(i,j)*(c4(i,k)*c4(j,l)+c2(i,k)*c2(j,l))*bb + &
                          s3(i,j)*(c5(i,k)*c5(j,l)+c1(i,k)*c1(j,l))
            end do
          end do
        end do
      end do
      return
      end subroutine diat
      subroutine diat2(na, esa, epa, r12, nb, esb, epb, s, a0)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE overlaps_C, only : sa, sb, a, b, isp, ips
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: na
      integer  :: nb
      double precision  :: esa
      double precision  :: epa
      double precision , intent(in) :: r12
      double precision  :: esb
      double precision  :: epb
      double precision , intent(in) :: a0
      double precision , intent(inout) :: s(3,3,3)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(17) :: inmb
      integer , dimension(78) :: iii
      integer :: jmax, jmin, nbond, ii
      double precision :: rab, rab4, w, rt3, d, e, rab6

      save inmb, iii
!-----------------------------------------------
!***********************************************************************
!
! OVERLP CALCULATES OVERLAPS BETWEEN ATOMIC ORBITALS FOR PAIRS OF ATOMS
!        IT CAN HANDLE THE ORBITALS 1S, 2S, 3S, 2P, AND 3P.
!
!***********************************************************************
      data inmb/ 1, 0, 2, 2, 3, 4, 5, 6, 7, 0, 8, 8, 8, 9, 10, 11, 12/
!     NUMBERING CORRESPONDS TO BOND TYPE MATRIX GIVEN ABOVE
!      THE CODE IS
!
!     III=1      FIRST - FIRST  ROW ELEMENTS
!        =2      FIRST - SECOND
!        =3      FIRST - THIRD
!        =4      SECOND - SECOND
!        =5      SECOND - THIRD
!        =6      THIRD - THIRD
      data iii/  1, 2, 4, 2, 4, 4, 2, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, &
        2, 4, 4, 4, 4, 4, 4, 3, 5, 5, 5, 5, 5, 5, 6, 3, 5, 5, 5, 5, 5, 5, 6, 6, &
        3, 5, 5, 5, 5, 5, 5, 6, 6, 6, 3, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 3, 5, &
        5, 5, 5, 5, 5, 6, 6, 6, 6, 6/
!
!      ASSIGNS BOND NUMBER
!
      jmax = max0(inmb(na),inmb(nb))
      jmin = min0(inmb(na),inmb(nb))
      nbond = (jmax*(jmax - 1))/2 + jmin
      ii = iii(nbond)
      s = 0.D0
      rab = r12/a0
      select case (ii)
!
!     ------------------------------------------------------------------
! *** THE ORDERING OF THE ELEMENTS WITHIN S IS
! *** S(1,1,1)=(S(B)/S(A))
! *** S(1,2,1)=(P-SIGMA(B)/S(A))
! *** S(2,1,1)=(S(B)/P-SIGMA(A))
! *** S(2,2,1)=(P-SIGMA(B)/P-SIGMA(A))
! *** S(2,2,2)=(P-PI(B)/P-PI(A))
!     ------------------------------------------------------------------
! *** FIRST ROW - FIRST ROW OVERLAPS
!
      case default
        call set (esa, esb, na, nb, rab, ii)
        s(1,1,1) = .25D0*sqrt((sa*sb*rab*rab)**3)*(a(3)*b(1)-b(3)*a(1))
        return
!
! *** FIRST ROW - SECOND ROW OVERLAPS
!
      case (2)
        call set (esa, esb, na, nb, rab, ii)
        rab4 = rab**4*0.125D00
        w = sqrt(sa**3*sb**5)*rab4
        s(1,1,1) = sqrt(1.D00/3.D00)
        s(1,1,1) = w*s(1,1,1)*(a(4)*b(1)-b(4)*a(1)+a(3)*b(2)-b(3)*a(2))
        if (na > 1) call set (epa, esb, na, nb, rab, ii)
        if (nb > 1) call set (esa, epb, na, nb, rab, ii)
        w = sqrt(sa**3*sb**5)*rab4
        s(isp,ips,1) = w*(a(3)*b(1)-b(3)*a(1)+a(4)*b(2)-b(4)*a(2))
        return
!
! *** FIRST ROW - THIRD ROW OVERLAPS
!
      case (3)
        call set (esa, esb, na, nb, rab, ii)
        rab4 = rab**5*0.0625D00
        w = sqrt(sa**3*sb**7/22.5D00)*rab4
        s(1,1,1) = w*(a(5)*b(1)-b(5)*a(1)+(a(4)*b(2)-b(4)*a(2))*2.D0)
        if (na > 1) call set (epa, esb, na, nb, rab, ii)
        if (nb > 1) call set (esa, epb, na, nb, rab, ii)
        w = sqrt(sa**3*sb**7/7.5D00)*rab4
        s(isp,ips,1) = w*(a(4)*(b(1)+b(3))-b(4)*(a(1)+a(3))+b(2)*(a(3)+a(5))-a(2)*(b(3)+b(5)))
        return
!
! *** SECOND ROW - SECOND ROW OVERLAPS
!
      case (4)
        call set (esa, esb, na, nb, rab, ii)
        rab4 = rab**5*0.0625D00
        w = sqrt((sa*sb)**5)*rab4
        s(1,1,1) = w*(a(5)*b(1)+b(5)*a(1)-2.0D00*a(3)*b(3))/3.0D00
        call set (esa, epb, na, nb, rab, ii)
        if (na > nb) call set (epa, esb, na, nb, rab, ii)
        w = sqrt((sa*sb)**5)*rab4
        rt3 = 1.D00/sqrt(3.D00)
        d = a(4)*(b(1)-b(3)) - a(2)*(b(3)-b(5))
        e = b(4)*(a(1)-a(3)) - b(2)*(a(3)-a(5))
        s(isp,ips,1) = w*rt3*(d + e)
        call set (epa, esb, na, nb, rab, ii)
        if (na > nb) call set (esa, epb, na, nb, rab, ii)
        w = sqrt((sa*sb)**5)*rab4
        d = a(4)*(b(1)-b(3)) - a(2)*(b(3)-b(5))
        e = b(4)*(a(1)-a(3)) - b(2)*(a(3)-a(5))
        s(ips,isp,1) = w*rt3*(d - e)
        call set (epa, epb, na, nb, rab, ii)
        w = sqrt((sa*sb)**5)*rab4
        s(2,2,1) = -w*(b(3)*(a(5)+a(1))-a(3)*(b(5)+b(1)))
        s(2,2,2) = 0.5D0*w*(a(5)*(b(1)-b(3))-b(5)*(a(1)-a(3))-a(3)*b(1)+b(3)*a(1))
        return
!
! *** SECOND ROW - THIRD ROW OVERLAPS
!
      case (5)
        call set (esa, esb, na, nb, rab, ii)
        rab6 = rab**6*0.03125D0/sqrt(7.5D0)
        w = sqrt(sa**5*sb**7)*rab6
        rt3 = 1.D00/sqrt(3.D00)
        s(1,1,1) = w*(a(6)*b(1)+a(5)*b(2)-2.D0*(a(4)*b(3)+a(3)*b(4))+a(2)*b(5)+a(1)*b(6))/3.D00
!
        call set (esa, epb, na, nb, rab, ii)
        if (na > nb) call set (epa, esb, na, nb, rab, ii)
        w = sqrt(sa**5*sb**7)*rab6
        s(isp,ips,1) = w*rt3*(a(6)*b(2)+a(5)*b(1)-2.D0*(a(4)*b(4)+a(3)*b(3))+a(2)*b(6)+a(1)*b(5))
!
        call set (epa, esb, na, nb, rab, ii)
        if (na > nb) call set (esa, epb, na, nb, rab, ii)
        w = sqrt(sa**5*sb**7)*rab6
        s(ips,isp,1) = -w*rt3*(a(5)*(2.D0*b(3)-b(1))-b(5)*(2.D0*a(3)-a(1))-a(2) &
          *(b(6)-2.D0*b(4))+b(2)*(a(6)-2.D0*a(4)))
!
        call set (epa, epb, na, nb, rab, ii)
        w = sqrt(sa**5*sb**7)*rab6
        s(2,2,1) = -w*(b(4)*(a(1)+a(5))-a(4)*(b(1)+b(5))+b(3)*(a(2)+a(6))-a(3)* &
          (b(2)+b(6)))
        s(2,2,2) = 0.5D0*w*(a(6)*(b(1)-b(3))-b(6)*(a(1)-a(3))+a(5)*(b(2)-b(4)) - &
          b(5)*(a(2)-a(4))-a(4)*b(1)+b(4)*a(1)-a(3)*b(2)+b(3)*a(2))
        return
!
! *** THIRD ROW - THIRD ROW OVERLAPS
!
      case (6)
        call set (esa, esb, na, nb, rab, ii)
        rab4 = rab**7/480.D00
        w = sqrt((sa*sb)**7)*rab4
        rt3 = 1.D00/sqrt(3.D00)
        s(1,1,1) = w*(a(7)*b(1)-3.D00*(a(5)*b(3)-a(3)*b(5))-a(1)*b(7))/3.D00
        call set (esa, epb, na, nb, rab, ii)
        if (na > nb) call set (epa, esb, na, nb, rab, ii)
        w = sqrt((sa*sb)**7)*rab4
        d = a(6)*(b(1)-b(3)) - 2.D00*a(4)*(b(3)-b(5)) + a(2)*(b(5)-b(7))
        e = b(6)*(a(1)-a(3)) - 2.D00*b(4)*(a(3)-a(5)) + b(2)*(a(5)-a(7))
        s(isp,ips,1) = w*rt3*(d - e)
        call set (epa, esb, na, nb, rab, ii)
        if (na > nb) call set (esa, epb, na, nb, rab, ii)
        w = sqrt((sa*sb)**7)*rab4
        d = a(6)*(b(1)-b(3)) - 2.D00*a(4)*(b(3)-b(5)) + a(2)*(b(5)-b(7))
        e = b(6)*(a(1)-a(3)) - 2.D00*b(4)*(a(3)-a(5)) + b(2)*(a(5)-a(7))
        s(ips,isp,1) = -w*rt3*((-d) - e)
        call set (epa, epb, na, nb, rab, ii)
        w = sqrt((sa*sb)**7)*rab4
        d = a(3)*(b(7)+b(3)+b(3)) - a(5)*(b(1)+b(5)+b(5)) - b(5)*a(1) + a(7)*b(3)
        s(2,2,1) = -w*d
        d = a(7)*(b(1)-b(3)) + b(7)*(a(1)-a(3))
        e = a(5)*(b(5)-b(3)-b(1)) + b(5)*(a(5)-a(3)-a(1)) + 2.D0*a(3)*b(3)
        s(2,2,2) = 0.5D0*w*(d + e)
        return
      end select
!
      end subroutine diat2
      double precision function ss (na, nb, la1, lb1, m1, ua, ub, r1, a0)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use overlaps_C, only : fact
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: na
      integer , intent(in) :: nb
      integer , intent(in) :: la1
      integer , intent(in) :: lb1
      integer , intent(in) :: m1
      double precision , intent(in) :: ua
      double precision , intent(in) :: ub
      double precision , intent(in) :: r1
      double precision , intent(in) :: a0
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: m, lb, la, i, i1, j, n, lam1, lbm1, ia, ic, ib, id, iab, k1, &
        k2, k3, k4, k5, iaf, k6, ibf
      double precision, dimension(0:2,0:2,0:2) :: aff
      double precision, dimension(0:19) :: af, bf
      double precision, dimension(0:12,0:12) :: bi
      double precision :: r, p, b, quo, sum, sum1
      logical :: first

      save first, aff, bi
!-----------------------------------------------
      data first/ .TRUE./
      data aff/ 27*0.D0/
      m = m1 - 1
      lb = lb1 - 1
      la = la1 - 1
      r = r1/a0
      if (first) then
        first = .FALSE.
!
!           INITIALISE SOME CONSTANTS
!
!                  BINOMIALS
!
        do i = 0, 12
          bi(i,0) = 1.D0
          bi(i,i) = 1.D0
        end do
        do i = 0, 11
          i1 = i - 1
          bi(i+1,1:i1+1) = bi(i,1:i1+1) + bi(i,:i1)
        end do
        aff(0,0,0) = 1.D0
        aff(1,0,0) = 1.D0
        aff(1,1,0) = sqrt(0.5D0)
        aff(2,0,0) = 1.5D0
        aff(2,1,0) = sqrt(1.5D0)
        aff(2,2,0) = sqrt(0.375D0)
!
!   AFF(2,0,2) CORRESPONDS TO C(2,0,1) IN THE MANUAL
!
        aff(2,0,2) = -0.5D0
      end if
      p = (ua + ub)*r*0.5D0
      b = (ua - ub)*r*0.5D0
      quo = 1/p
      af(0) = quo*exp((-p))
      do n = 1, 19
        af(n) = n*quo*af(n-1) + af(0)
      end do
      call bfn (b, bf)
      sum = 0.D0
      lam1 = la - m
      lbm1 = lb - m
!
!          START OF OVERLAP CALCULATION PROPER
!
      do i = 0, lam1, 2
        ia = na + i - la
        ic = la - i - m
        do j = 0, lbm1, 2
          ib = nb + j - lb
          id = lb - j - m
          sum1 = 0.D0
          iab = ia + ib
!
!   In the Manual ka = K6
!                 kb = K5
!                 Pa = K1
!                 Pb = K2
!                 qa = K3
!                 qb = K4
!
          do k1 = 0, ia
            do k2 = 0, ib
              do k3 = 0, ic
                do k4 = 0, id
                  do k5 = 0, m
                    iaf = iab - k1 - k2 + k3 + k4 + 2*k5
                    do k6 = 0, m
                      ibf = k1 + k2 + k3 + k4 + 2*k6
                      sum1 = sum1 + bi(id,k4)*bi(ic,k3)*bi(ib,k2)*bi(ia,k1)*bi(m,k5)*bi(m,k6)* &
                      (1 - 2*mod(m + k2 + k4 + k5 + k6,2))*af(iaf)*bf(ibf)
                    end do
                  end do
                end do
              end do
            end do
          end do
          sum = sum + sum1*aff(la,m,i)*aff(lb,m,j)
        end do
      end do
      ss = sum*r**(na + nb + 1)*ua**na*ub**nb/2.D0* &
      sqrt(ua*ub/(fact(na+na)*fact(nb+ nb))*((la+la+1)*(lb+lb+1)))
      return
      end function ss
