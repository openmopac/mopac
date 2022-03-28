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

  module YGXX_SimpleGaussianInts
  implicit none
  public :: GGJPP_IntQ, GPJPS_Int, GPJPS_IntQ, SSJPS_IntQ, SSJPP_Int, SPJPS_Int, GSJPS_IntQ

  private

  integer, parameter :: XSN = 3
  double precision, parameter :: XSZ(3) = (/ 0.4068833884920483d0, 0.09179674341953627d0, 3.515225231758639d0 /)
  double precision, parameter :: XSC(3) = (/ 0.3620527755096057d0, 0.6262050528632612d0, 0.01174217162750757d0 /)
  double precision, parameter :: sqrt_pi = 1.772453850905516027298167483341d0

  contains

  subroutine GSJPS_IntQ(Rab, na, Ca, Za, dZadQa, Zb, dZbdQb, X, Xa, Xb)
    implicit none
    double precision, intent(in) :: Rab(3)
    integer, intent(in) :: na
    double precision, intent(in) :: Za(na), Ca(na), dZadQa(na), Zb, dZbdQb
    double precision, intent(OUT) :: X(3, 1), Xa(3, 1), Xb(3, 1)
!
!  Local variables
!
    double precision :: Gb(XSN)
    double precision :: dGb(XSN)
    Gb = XSZ * Zb**2
    dGb = XSZ * 2*Zb * dZbdQb
    call GGJPS_IntQ(Rab, na, Ca, Za, dZadQa, XSN, XSC, Gb, dGb, X, Xa, Xb)
  end subroutine GSJPS_IntQ

  subroutine SSJPS_IntQ(Rab, Za, dZadQa, Zb, dZbdQb, X, Xa, Xb)
    implicit none
    double precision, intent(in) :: Rab(3), Za, Zb, dZadQa, dZbdQb
    double precision, intent(OUT) :: X(3, 1), Xa(3, 1), Xb(3, 1)
!
!  Local variables
!
    double precision :: Ga(XSN), Gb(XSN)
    double precision :: dGa(XSN), dGb(XSN)
    Ga = XSZ * Za**2
    Gb = XSZ * Zb**2
    dGa = XSZ * 2*Za * dZadQa
    dGb = XSZ * 2*Zb * dZbdQb
    call GGJPS_IntQ(Rab, XSN, XSC, Ga, dGa, XSN, XSC, Gb, dGb, X, Xa, Xb)
  end subroutine SSJPS_IntQ

  subroutine GGJPS_IntQ(Rab, na, Ca, Za, dZadQa, nb, Cb, Zb, dZbdQb, X, Xa, Xb)
    implicit none
    double precision, intent(in) :: Rab(3)
    integer, intent(in) :: na, nb
    double precision, intent(in) :: Ca(na), Za(na), dZadQa(na), Cb(nb), Zb(nb), dZbdQb(nb)
    double precision, intent(OUT) :: X( 3, 1 ), Xa( 3, 1 ), Xb( 3, 1 )
!
!  Local variables
!
    integer :: i, j
    double precision :: r2, r, b, br, cc, j00, d00dr, j10, e, U(3), dbdQa, dbdQb, d00db, &
      dedb, d00drdb, d10db, d10dQa, d10dQb
    X = 0.d0
    Xa = 0.d0
    Xb = 0.d0
    r2 = Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3)
    if ( r2 < 1.d-25 ) then
       X = 0.d0
       Xa = 0.d0
       Xb = 0.d0
    else
       r = sqrt(r2)
       U = Rab/r
       j10 = 0.d0
       d10dQa = 0.d0
       d10dQb = 0.d0
       do i=1, na
          do j=1, nb
             b  = sqrt( Za(i) * Zb(j) / ( Za(i) + Zb(j) ) )
             dbdQa = dZadQa(i) * ( Zb(j)/( Za(i)+Zb(j) ) )**2  / (2.d0*b)
             dbdQb = dZbdQb(j) * ( Za(i)/( Za(i)+Zb(j) ) )**2  / (2.d0*b)
             br = b * r
             cc = Ca(i) * Cb(j)
             e = (2.d0 / sqrt_PI) * b * EXP(-br*br)
             j00    = erf(br)/r
             d00dr  = e/r - j00/r
             d00db = (2.d0 / sqrt_PI) * EXP(-br*br)
             dedb  = (2.d0 / sqrt_PI) * (1.d0-2.d0*br*br) * EXP(-br*br)
             d00drdb = dedb/r - d00db/r
             j10    = j10 + cc * d00dr
             d10db  = cc * d00drdb
             d10dQa = d10dQa + d10db * dbdQa
             d10dQb = d10dQb + d10db * dbdQb
          end do
       end do
       X(1, 1) = X(1, 1) + j10 * U(3)
       X(2, 1) = X(2, 1) + j10 * U(1)
       X(3, 1) = X(3, 1) + j10 * U(2)
       Xa(1, 1) = Xa(1, 1) + d10dQa * U(3)
       Xa(2, 1) = Xa(2, 1) + d10dQa * U(1)
       Xa(3, 1) = Xa(3, 1) + d10dQa * U(2)
       Xb(1, 1) = Xb(1, 1) + d10dQb * U(3)
       Xb(2, 1) = Xb(2, 1) + d10dQb * U(1)
       Xb(3, 1) = Xb(3, 1) + d10dQb * U(2)
    end if
  end subroutine GGJPS_IntQ

  subroutine SSJPP_Int(Rab, Za, Zb, X)
    implicit none
    double precision, intent(in) :: Rab(3)
    double precision, intent(in) :: Za, Zb
    double precision, intent(OUT) :: X(3, 3)
!
!  Local variables
!
    double precision :: Ga(XSN), Gb(XSN)
    Ga = XSZ * Za**2
    Gb = XSZ * Zb**2
    call GGJPP_Int(Rab, XSN, XSC, Ga, XSN, XSC, Gb, X)
  end subroutine SSJPP_Int




  subroutine GGJPP_Int(Rab, na, Ca, Za, nb, Cb, Zb, X)
    implicit none
    double precision, intent(in) :: Rab(3)
    integer, intent(in) :: na, nb
    double precision, intent(in) :: Ca(na), Za(na)
    double precision, intent(in) :: Cb(nb), Zb(nb)
    double precision, intent(OUT) :: X( 3, 3 )
!
!  Local variables
!

    integer :: i, j
    double precision :: r2, r
    double precision :: b, br, cc
    double precision :: j00, d00dr, d00dr2
    double precision :: j10
    double precision :: j11
    double precision :: e, dedr
    double precision :: U(3)
    double precision :: dUda(3, 3)
    integer :: ix, iy, jx, jy
    X = 0.d0
    r2 = Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3)
    if ( r2 < 1.d-25 ) then
       e = 0.d0
       do i=1, na
          do j=1, nb
             b  = sqrt( Za(i) * Zb(j) / ( Za(i) + Zb(j) ) )
             cc = Ca(i) * Cb(j)
             e = e + cc * 4.d0 * b**3 / ( 3.d0 * sqrt_PI )
          end do
       end do
       do i=1, 3
          X(i, i) = e
       end do
    else
       r = sqrt(r2)
       U = Rab/r
       do i=1, 3
          do j=1, 3
             dUda(i, j) = -U(i)*U(j)/r
          end do
          dUda(i, i) = dUda(i, i) + 1.d0/r
       end do
       j10 = 0.d0
       j11 = 0.d0
       do i=1, na
          do j=1, nb
             b  = sqrt( Za(i) * Zb(j) / ( Za(i) + Zb(j) ) )
             br = b * r
             cc = Ca(i) * Cb(j)
             e     =  (2.d0 / sqrt_PI) * b * EXP(-br*br)
             dedr  = -2*b*br * e
             j00    = erf(br)/r
             d00dr  = e/r - j00/r
             d00dr2 = dedr/r - e/r2 - d00dr/r + j00/r2
             j10    = j10 + cc * d00dr
             j11    = j11 + cc * d00dr2
          end do
       end do
       do iy=1, 3
          jy = MOD(iy+1, 3)+1
          do ix=1, 3
             jx = MOD(ix+1, 3)+1
             X(ix, iy) =   j11 * U(jx) * (-U(jy)) &
                  &     - j10 * dUda(jx, jy)
          end do
       end do
    end if
  end subroutine GGJPP_Int

  subroutine GGJPP_IntQ(Rab, na, Ca, Za, dZadQa, nb, Cb, Zb, dZbdQb, X, Xa, Xb)
    implicit none
    double precision, intent(in) :: Rab(3)
    integer, intent(in) :: na, nb
    double precision, intent(in) :: Ca(na), Za(na), dZadQa(na)
    double precision, intent(in) :: Cb(nb), Zb(nb), dZbdQb(nb)
    double precision, intent(OUT) :: X( 3, 3 )
    double precision, intent(OUT) :: Xa( 3, 3 )
    double precision, intent(OUT) :: Xb( 3, 3 )
!
!  Local variables
!
    integer :: i, j
    double precision :: r2, r
    double precision :: b, br, cc
    double precision :: j00, d00dr, d00dr2
    double precision :: j10
    double precision :: j11
    double precision :: e, dedr
    double precision :: U(3)
    double precision :: dUda(3, 3)
    integer :: ix, iy, jx, jy
    double precision :: dbdQa, dbdQb
    double precision :: dedb, dedrdb, d00db, d00drdb, d00dr2db
    double precision :: dedQa, dedQb
    double precision :: d10dQa, d10dQb
    double precision :: d11dQa, d11dQb
    X  = 0.d0
    Xa = 0.d0
    Xb = 0.d0
    r2 = Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3)
    if ( r2 < 1.d-25 ) then
       e = 0.d0
       dedQa = 0.d0
       dedQb = 0.d0
       do i=1, na
          do j=1, nb
             b  = sqrt( Za(i) * Zb(j) / ( Za(i) + Zb(j) ) )
             dbdQa = dZadQa(i) * ( Zb(j)/( Za(i)+Zb(j) ) )**2  / (2.d0*b)
             dbdQb = dZbdQb(j) * ( Za(i)/( Za(i)+Zb(j) ) )**2  / (2.d0*b)
             cc = Ca(i) * Cb(j)
             e = e + cc * 4.d0 * b**3 / ( 3.d0 * sqrt_PI )
             dedb =  cc * 12.d0 * b**2 / ( 3.d0 * sqrt_PI )
             dedQa = dedQa + dedb * dbdQa
             dedQb = dedQb + dedb * dbdQb
          end do
       end do
       do i=1, 3
          X(i, i) = e
          Xa(i, i) = dedQa
          Xb(i, i) = dedQb
       end do
    else
       r = sqrt(r2)
       U = Rab/r
       do i=1, 3
          do j=1, 3
             dUda(i, j) = -U(i)*U(j)/r
          end do
          dUda(i, i) = dUda(i, i) + 1.d0/r
       end do
       j10 = 0.d0
       j11 = 0.d0
       d10dQa = 0.d0
       d10dQb = 0.d0
       d11dQa = 0.d0
       d11dQb = 0.d0
       do i=1, na
          do j=1, nb
             b  = sqrt( Za(i) * Zb(j) / ( Za(i) + Zb(j) ) )
             dbdQa = dZadQa(i) * ( Zb(j)/( Za(i)+Zb(j) ) )**2  / (2.d0*b)
             dbdQb = dZbdQb(j) * ( Za(i)/( Za(i)+Zb(j) ) )**2  / (2.d0*b)
             br = b * r
             cc = Ca(i) * Cb(j)
             e     =  (2.d0 / sqrt_PI) * b * EXP(-br*br)
             dedr  = -2*b*br * e
             dedb = (2.d0 / sqrt_PI) * (1.d0-2.d0*br*br) * EXP(-br*br)
             dedrdb = -4*br*e - 2*b*br*dedb
             j00    = erf(br)/r
             d00dr  = e/r - j00/r
             d00dr2 = dedr/r - e/r2 - d00dr/r + j00/r2
             d00db    = (2.d0 / sqrt_PI) * EXP(-br*br)
             d00drdb  = dedb/r - d00db/r
             d00dr2db = dedrdb/r - dedb/r2 - d00drdb/r + d00db/r2
             j10    = j10 + cc * d00dr
             j11    = j11 + cc * d00dr2
             d10dQa = d10dQa + cc * d00drdb * dbdQa
             d10dQb = d10dQb + cc * d00drdb * dbdQb
             d11dQa = d11dQa + cc * d00dr2db * dbdQa
             d11dQb = d11dQb + cc * d00dr2db * dbdQb
          end do
       end do
       do iy=1, 3
          jy = MOD(iy+1, 3)+1
          do ix=1, 3
             jx = MOD(ix+1, 3)+1
             X(ix, iy) =   j11 * U(jx) * (-U(jy))  - j10 * dUda(jx, jy)
             Xa(ix, iy) =   d11dQa * U(jx) * (-U(jy))  - d10dQa * dUda(jx, jy)
             Xb(ix, iy) =   d11dQb * U(jx) * (-U(jy))  - d10dQb * dUda(jx, jy)
          end do
       end do
    end if
  end subroutine GGJPP_IntQ

  subroutine SPJPS_Int(Rab, Za, X)
    implicit none
    double precision, intent(in) :: Rab(3), Za
    double precision, intent(OUT) :: X(3, 1)
!
!  Local variables
!
    double precision :: Ga(XSN)
    Ga = XSZ * Za**2
    call GPJPS_Int(Rab, XSN, XSC, Ga, X)
  end subroutine SPJPS_Int




  subroutine GPJPS_Int(Rab, na, Ca, Za, X)
    implicit none
    double precision, intent(in) :: Rab(3)
    integer, intent(in) :: na
    double precision, intent(in) :: Ca(na), Za(na)
    double precision, intent(OUT) :: X( 3, 1 )
!
!  Local variables
!
    integer :: i
    double precision :: r2, r, b, br, cc, j00, d00dr, j10, e, U(3)
    X = 0.d0
    r2 = Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3)
    if ( r2 < 1.d-25 ) then
       X = 0.d0
    else
       r = sqrt(r2)
       U = Rab/r
       j10 = 0.d0
       do i=1, na
          b  = sqrt( Za(i)  )
          br = b * r
          cc = Ca(i)
          e = (2.d0 / sqrt_PI) * b * EXP(-br*br)
          j00    = erf(br)/r
          d00dr  = e/r - j00/r
          j10    = j10 + cc * d00dr
       end do
       X(1, 1) = X(1, 1) + j10 * U(3)
       X(2, 1) = X(2, 1) + j10 * U(1)
       X(3, 1) = X(3, 1) + j10 * U(2)
    end if
  end subroutine GPJPS_Int



  subroutine GPJPS_IntQ(Rab, na, Ca, Za, dZadQa, X, Xa)
    implicit none
    double precision, intent(in) :: Rab(3)
    integer, intent(in) :: na
    double precision, intent(in) :: Ca(na), Za(na), dZadQa(na)
    double precision, intent(OUT) :: X( 3, 1 ), Xa( 3, 1 )
!
!  Local variables
!
    integer :: i
    double precision :: r2, r, b, br, cc, j00, d00dr, j10, e, U(3), d10dQa, dbdQa, &
      dedb, d00db, d00drdb
    X = 0.d0
    Xa = 0.d0
    r2 = Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3)
    if ( r2 < 1.d-25 ) then
       X = 0.d0
       Xa = 0.d0
    else
       r = sqrt(r2)
       U = Rab/r
       j10 = 0.d0
       d10dQa = 0.d0
       do i=1, na
          b  = sqrt( Za(i) )
          dbdQa = dZadQa(i) / ( 2.d0*b )
          br = b * r
          cc = Ca(i)
          e = (2.d0 / sqrt_PI) * b * EXP(-br*br)
          dedb = (2.d0 / sqrt_PI) * (1.d0-2.d0*br*br) * EXP(-br*br)
          j00    = erf(br)/r
          d00dr  = e/r - j00/r
          d00db = (2.d0 / sqrt_PI) * EXP(-br*br)
          d00drdb = dedb/r - d00db/r
          j10    = j10 + cc * d00dr
          d10dQa = d10dQa + cc * d00drdb * dbdQa
       end do
       X(1, 1) = X(1, 1) + j10 * U(3)
       X(2, 1) = X(2, 1) + j10 * U(1)
       X(3, 1) = X(3, 1) + j10 * U(2)
       Xa(1, 1) = Xa(1, 1) + d10dQa * U(3)
       Xa(2, 1) = Xa(2, 1) + d10dQa * U(1)
       Xa(3, 1) = Xa(3, 1) + d10dQa * U(2)
    end if
  end subroutine GPJPS_IntQ
  end module YGXX_SimpleGaussianInts

  subroutine Invert_Symmetric_Matrix( N, A, AI )
    implicit none
    integer, intent(in) :: N
    double precision, intent(in) :: A(N, N)
    double precision, intent(OUT) :: AI(N, N)
!
!  Local variables
!
    integer :: INFO, i, j
    AI = A
    call dpotrf( "U", N, AI(1, 1), N, INFO )
    call dpotri( "U", N, AI(1, 1), N, INFO )
    do j=2, N
      do i=1, j-1
          AI(j, i) = AI(i, j)
      end do
    end do
  end subroutine Invert_Symmetric_Matrix

 subroutine CPE_QdepDipole_Contribution(N_qm, Id_qm, Crd_qm, Q_qm, Zeta_qm, E, Pot, C )
    implicit none
    integer, intent(in) :: N_qm          ! Number of qm atoms
    integer, intent(in) :: Id_qm( N_qm ) ! Parameter indices
    double precision, intent(in) :: Crd_qm( 3, N_qm ) ! qm coordinates (xyz)
    double precision, intent(in) :: Q_qm( N_qm )      ! qm electric partial charges
    double precision, intent(in) :: Zeta_qm( N_qm )   ! qm aux slater exponent
    double precision, intent(OUT) :: E                ! response energy contribution
    double precision, intent(OUT) :: Pot( N_qm )      ! effective monopolar potential correction to be mapped into the fock matrix
    double precision, intent(OUT) :: C( 3 * N_qm )    ! Response dipole coefs (zxy)
!
!  Local variables
!
    integer :: N
    double precision, ALLOCATABLE :: Eta(:, :), dEtadQa(:, :), dEtadQb(:, :), EtaInv(:, :), &
      M(:), Mmat(:, :), dMdQa(:), Q0(:)
    integer :: i, i1, i2, j, j1, j2
    Pot = 0.d0
    N = 3 * N_qm
    allocate( Eta(N, N), dEtadQa(N, N), dEtadQb(N, N), EtaInv(N, N), M(N), Mmat( N, N_qm ) )
    allocate( dMdQa(N), Q0( N_qm) )
    Q0( 1:N_qm ) = Q_qm( 1:N_qm )
    call Cpt_SecondOrderMatrix( N_qm, Id_qm, Crd_qm, Q_qm, Eta, dEtadQa, dEtadQb )
    call Cpt_FirstOrderMatrix( N_qm, Id_qm, Crd_qm, Q_qm, Zeta_qm, Mmat, dMdQa )
    call Invert_Symmetric_Matrix( N, Eta, EtaInv )
    M = matmul( Mmat, Q0 )
    C = -matmul(EtaInv, M )
    E =  dot_product( C, M ) + 0.5d0 * dot_product( C, matmul( Eta, C ) )
    do i=1, N_qm
       i1 = (i-1)*3 + 1
       i2 = i1 + 2
       Pot(i) = Pot(i) + dot_product( C, Mmat(:, i) )
       Pot(i) = Pot(i) + dot_product( C(i1:i2), dMdQa(i1:i2) )
    end do
    do i=1, N_qm
       i1 = (i-1)*3 + 1
       i2 = i1 + 2
       do j=1, N_qm
          j1 = (j-1)*3 + 1
          j2 = j1 + 2
          Pot(i) = Pot(i) + 0.5d0 * dot_product(C(i1:i2), matmul(dEtadQa(i1:i2, j1:j2), C(j1:j2)))
          Pot(j) = Pot(j) + 0.5d0 * dot_product(C(i1:i2), matmul(dEtadQb(i1:i2, j1:j2), C(j1:j2)))
       end do
    end do
  end subroutine CPE_QdepDipole_Contribution

  subroutine Cpt_FirstOrderMatrix( N_qm, Id_qm, Crd_qm, Q_qm, Zeta_qm, Mmat, dMdQa )
    use YGXX_SimpleGaussianInts, only : GSJPS_IntQ
    USE parameters_C, only : CPE_Z0, CPE_B, CPE_Xlo, CPE_Xhi
    implicit none
    integer, intent(in) :: N_qm          ! Number of qm atoms
    integer, intent(in) :: Id_qm( N_qm ) ! Parameter indices
    double precision, intent(in) :: Crd_qm( 3, N_qm ) ! qm coordinates (xyz)
    double precision, intent(in) :: Q_qm( N_qm )      ! qm electric partial charges
    double precision, intent(in) :: Zeta_qm( N_qm )   ! qm aux slater exponent
    double precision, intent(OUT) :: Mmat( 3*N_qm, N_qm), dMdQa( 3*N_qm )
!
!  Local variables
!
    double precision :: Rab(3), X(3, 1), Xa(3, 1), Xb(3, 1), Za, dZadQa, Zb, dZbdQb
    integer :: i, i1, i2, ia, j, ib
    double precision :: rlo, rhi, r, s
    Mmat  = 0.d0
    dMdQa = 0.d0
    do i=1, N_qm
       i1 = (i-1)*3 + 1
       i2 = i1 + 2
       ia = Id_qm(i)
       Za = CPE_Z0(ia) * EXP( CPE_B(ia) * Q_qm(i) )
       dZadQa = CPE_B(ia) * Za
       do j=1, N_qm
          ib = Id_qm(j)
          Rab = Crd_qm(1:3, i) - Crd_qm(1:3, j)
          r   = sqrt( Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3) )
          Zb = Zeta_qm(j)
          dZbdQb = 0.d0
          call GSJPS_IntQ(Rab, 1, (/1.d0/), (/Za/), (/dZadQa/), Zb, dZbdQb, X, Xa, Xb)
          rlo = CPE_Xlo(ia) + CPE_Xlo(ib)
          rhi = CPE_Xhi(ia) + CPE_Xhi(ib)
          call SwitchOn(r, rlo, rhi, s)
          X  = s * X
          Xa = s * Xa
          Mmat( i1:i2, j ) = Mmat( i1:i2, j ) + X(1:3, 1)
          dMdQa( i1:i2 ) = dMdQa( i1:i2 ) + Xa(1:3, 1) * Q_qm(j)
       end do
    end do
    if (xb(1, 1) > 1.d9) return
  end subroutine Cpt_FirstOrderMatrix

  subroutine Cpt_SecondOrderMatrix( N_qm, Id_qm, Crd_qm, Q_qm, Eta, dEtadQa, dEtadQb )
    use YGXX_SimpleGaussianInts, only : GGJPP_IntQ
    USE parameters_C, only : CPE_Z0, CPE_B
    implicit none
    integer, intent(in) :: N_qm          ! Number of qm atoms
    integer, intent(in) :: Id_qm( N_qm ) ! Parameter indices
    double precision, intent(in) :: Crd_qm( 3, N_qm ) ! qm coordinates (xyz)
    double precision, intent(in) :: Q_qm( N_qm )      ! qm electric partial charges

    double precision, intent(OUT) :: Eta( 3*N_qm, 3*N_qm ), dEtadQa( 3*N_qm, 3*N_qm ), &
      dEtadQb( 3*N_qm, 3*N_qm )
!
!  Local variables
!
    double precision :: X(3, 3), Xa(3, 3), Xb(3, 3), Za, dZadQa, Zb, dZbdQb, Rab(3)
    integer :: i, i1, i2, j, j1, j2, ia, ib
    Eta     = 0.d0
    dEtadQa = 0.d0
    dEtadQb = 0.d0
    do i=1, N_qm
       i1 = (i-1)*3 + 1
       i2 = i1 + 2
       ia = Id_qm(i)
       Za = CPE_Z0(ia) * EXP( CPE_B(ia) * Q_qm(i) )
       dZadQa = CPE_B(ia) * Za
       do j=1, N_qm
          j1 = (j-1)*3 + 1
          j2 = j1 + 2
          ib = Id_qm(j)
          Zb = CPE_Z0(ib) * EXP( CPE_B(ib) * Q_qm(j) )
          dZbdQb = CPE_B(ib) * Zb
          Rab = Crd_qm(1:3, i) - Crd_qm(1:3, j)
          call GGJPP_IntQ(Rab, 1, (/1.d0/), (/Za/), (/dZadQa/), 1, (/1.d0/), (/Zb/), (/dZbdQb/), X, Xa, Xb)
          Eta( i1:i2, j1:j2 ) = X
          dEtadQa( i1:i2, j1:j2 ) = Xa
          dEtadQb( i1:i2, j1:j2 ) = Xb
       end do
    end do
  end subroutine Cpt_SecondOrderMatrix

  subroutine SwitchOn(x, xlo, xhi, f)
    implicit none
    double precision, intent(in) :: x, xlo, xhi
    double precision, intent(OUT) :: f
!
!  Local variables
!
    double precision :: s, sf
!
! switch off
!
    if ( x >= xhi ) then
       f = 0.d0
    else if ( x <= xlo ) then
       f = 1.d0
    else
       s  = (xhi-x)/(xhi-xlo)
       sf =    10.d0 * s**3 &
            & -15.d0 * s**4 &
            & + 6.d0 * s**5
       f = sf
    end if
!
! switch on
!
    f = 1.d0 - f
  end subroutine SwitchOn

  Subroutine CPE_energy(Ene, Pot, C)
!
! Evaluates the total energy due to CPE (chemical-potential equalization)
!
    use molkst_C, only : numat
    use Common_arrays_C, only : nat, coord, chrg
    use funcon_C, only : a0, eV, fpc_9
    USE parameters_C, only : CPE_Zeta
    implicit none
    double precision :: Pot(numat)! Effective monopolar potential correction to be mapped into the Fock matrix
    double precision :: C(3*numat)! Response dipole coefs (zxy)
    double precision :: Ene       ! Response energy contribution
!
!  Local variables
!
    double precision, allocatable :: Crd_qm(:,:), Zeta_qm(:)
!=============================================================================
!  Initialization
!=============================================================================
    allocate(Zeta_qm(numat), Crd_qm(3, numat))
    Crd_qm(:,:numat) = coord(:,:numat)/a0 !  Convert geometry from Angstroms to Bohrs
    Zeta_qm(:numat) = CPE_Zeta(nat(:numat))
!=============================================================================
! Initialization complete.  Now do the calculation.
!=============================================================================
    call CPE_QdepDipole_Contribution( numat, nat, Crd_qm, chrg, Zeta_qm, Ene, Pot, C)
!
!  Convert energy into kcal.mol^(-1), Pot and C into eV. (This works in MOPAC, but does not work "stand-alone")
!
    Ene = Ene * eV*fpc_9
    Pot(:numat) = Pot(:numat)*eV
    C(:3*numat) = C(:3*numat)
    return
  end subroutine CPE_energy
