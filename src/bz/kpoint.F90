! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

subroutine kpoint (sec_det, nvecs, xk, yk, zk, v, eigs)
!***********************************************************************
!                                                                      *
!  Calculates the eigenvalues and complex eigenvectors for the         *
!  k-point (XK,YK,ZK). On exit, the eigenvalues are in EIGS, and the   *
!  complex eigenvectors are in V.                                      *
!                                                                      *
!***********************************************************************
  use common_common, only : phonon, ncell, keywrd, iw_new
  implicit none
  integer, intent(in) :: nvecs
  double precision, intent(in) :: xk, yk, zk
  complex, dimension(nvecs, nvecs), intent(inout) :: v
  double precision, dimension(2*nvecs), intent(inout) :: eigs
  double precision, dimension(nvecs, nvecs*ncell), intent(in) :: sec_det
! 
  integer :: i, j, l, m, nprt, k
  complex :: sum
  double precision :: fact = 6.0221367d23, &
  & c2pi = 1/(2.99792458d10*2*3.14159265358979d0) 
  complex, dimension(ncell) :: phase
  complex, dimension(nvecs,nvecs) :: complex_fmat
  double precision, allocatable, dimension (:,:) :: real_fmat, real_vecs
  double precision, allocatable, dimension (:) :: real_1d_fmat
! 
!  Get the complex phase for each unit cell at the point (xk, yk, zk)
!
!
  call get_phase(phase, xk, yk, zk)
!
!  Build the complex secular determinant for the point (xk, yk, zk)
!
  do i = 1, nvecs
    do j = i, nvecs
      sum = 0.0
      if (i == 25 .and. j == 25) then
        sum = 0.d0
      end if
      do m = 1, ncell
        sum = sum + Real(sec_det(i, j + m*nvecs - nvecs))*phase(m)
   !     if (i == 1 .and. j == 10) then
   !       write(iw_new,'(i3, f12.3, f8.3,f13.3, f8.3, f12.3, f12.3,f8.3, i7,2i3)')m, sum, &
   !       real(sec_det(i, j + m*nvecs - nvecs))*phase(m), &
   !         sec_det(i, j + m*nvecs - nvecs), phase(m), nijk(:,m)
   !     end if
      end do
      complex_fmat(i,j) = sum
      complex_fmat(j,i) = Conjg(sum)
    end do
  end do
  if (Index (keywrd, "SECDET") /= 0) then
    nprt = Min (25, nvecs)
    write(iw_new, '(a)') " Complex Secular Determinant, Real Part"
    do i = 1, nvecs
      write(iw_new, "(25f8.4)") (dble(complex_fmat(j,i)), j = 1,nprt)
    end do
    write(iw_new, '(a)') " Complex Secular Determinant, Imaginary Part"
    do i = 1, nvecs
      write(iw_new, "(25f8.4)") (Aimag (complex_fmat(j,i)), j = 1,nprt)
    end do
  end if
  allocate (real_fmat(2*nvecs,2*nvecs), real_vecs(2*nvecs,2*nvecs))
  allocate (real_1d_fmat(4*nvecs*nvecs)) 
!
!  Store complex_fmat as a real matrix, for use with the real diagonalizer.
!
    do i = 1,nvecs
      do j = 1,i
        real_fmat(j,i) = dble(complex_fmat(j,i))
        real_fmat(j+nvecs,i+nvecs) = real_fmat(j,i)
        real_fmat(j,i+nvecs) = Aimag(complex_fmat(j,i))
        real_fmat(i,j+nvecs) = -real_fmat(j,i+nvecs)
      end do
    end do
    if (nvecs < 30) then
!
! Solve for real eigenvalues and complex eigenvectors using a subroutine that
! diagonalizes a complex Hermitian matrix.
!
! The lower half triangle of the input complex square matrix is set to zero, to make debugging easier.
!
      do i = 1, nvecs
        do j = i + 1, nvecs
          complex_fmat(j,i) = 0.0
        end do
      end do
      call HEigensystem(nvecs, complex_fmat, eigs, v, 1)
    else
!
!  Convert to lower-half triangular form
!
       l = 0
      do i = 1,2*nvecs
        do j = 1,i
        l = l + 1
        real_1d_fmat(l) = real_fmat(j,i)
        end do
      end do
!
! The complex diagonalizer has not generated good eigenvectors.  
! Get eigenvectors from real diagonalizer.
!
      call rsp(real_1d_fmat, 2*nvecs, 2*nvecs, eigs, real_vecs)
      
      do k = 1,nvecs
        eigs(k) = eigs(2*k)
        do j = 1,nvecs
          v(j,k) = Cmplx(real_vecs(j,2*k),real_vecs(j+nvecs,2*k))
        end do
      end do
    end if
!
! Lock phase so that the largest coefficient is positive real
!
    call phase_lock(v,nvecs)
    if ( phonon ) eigs(1:nvecs) = Sign(sqrt(fact*Abs(eigs(1:nvecs)*1.d5))*c2pi,eigs(1:nvecs))
  if (Index (keywrd, " VECTORS") /= 0) then
    m = Min (25, nvecs)
    write(iw_new, '(a,i5,24i10)') " Eigenvalues",(i, i = 2,m)    
    write(iw_new, "(25F10.4)") (eigs(i), i = 1, m)
    write(iw_new, '(a)') "Complex Eigenvectors, Before Operation, Real part"
    do i = 1, nvecs
      write(iw_new, "(25F10.4)") (real(v(i, j)), j = 1, m)
    end do
    write(iw_new, '(a)') "Complex Eigenvectors, ", "Before Operation, Imaginary part"
    do i = 1, nvecs
      write(iw_new, "(25F10.4)") (Aimag (v(i, j)), j = 1, m)
    end do
  end if
end subroutine kpoint
!
!
!
subroutine get_phase(phase, xk, yk, zk)
!
!   Put the complex phase factor for each unit cell into the array "phase"
!
  use common_common, only : ncell, nijk
  implicit none
  double precision :: xk, yk, zk
  complex :: phase(ncell)
  integer :: loop, i, j ,k
  double precision :: point, twopi = 6.28318530717959d0
   do loop = 1, ncell
    i = nijk(1, loop)
    j = nijk(2, loop)
    k = nijk(3, loop)
    point = (xk*i + yk*j + zk*k)*twopi
    phase(loop) = Cmplx (Cos(point), Sin(point))
  end do
end subroutine get_phase

subroutine phase_lock(vecs,n)
!
!  Rotate the complex eigenvectors vecs(:,i) so that the largest coefficient is
!  real and positive
!
  implicit none
!
  integer :: n
  complex, dimension (n,n) :: vecs
!  Local
  integer :: i, j
  real :: asum
  complex :: ssum
!
!  Work out the (complex) rotation phase ssum
!
  do i = 1, n
    asum = 0.0
    do j = 1, n
      if (Abs(vecs(j,i)) > asum) then
        asum = Abs(vecs(j,i))
        ssum = vecs(j,i)
      end if
    end do
    ssum = ssum/asum
!
!  Rotate the eigenvector 
!
    do j = 1, n
      vecs(j,i) = cmplx(real(ssum)*real(vecs(j,i)) + aimag(ssum)*aimag(vecs(j,i)), &
      real(ssum)*aimag(vecs(j,i)) - aimag(ssum)*real(vecs(j,i)))
    end do
  end do
  return
  end subroutine phase_lock

