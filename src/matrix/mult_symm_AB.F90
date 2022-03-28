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

! ---------------------------------------------------------------------------------------
      subroutine mult_symm_AB(a, b, alpha, ndim, mdim, c, beta, iopc)

!        Use mod_vars_cuda, only: ngpus
        Use iso_c_binding
#ifdef GPU
        Use call_gemm_cublas
        Use mamult_cuda_i
        use common_arrays_C, only : ifact
#endif
        implicit none
        Integer :: iopc,ndim,mdim
        Integer :: i
#ifdef GPU
        integer :: igrid, iblock
        real :: tt
#endif
#if CC12
        real(c_float), dimension(mdim) :: a, b, c
#else
        real(c_double), dimension(mdim) :: a, b, c
#endif
        double precision, allocatable, dimension(:,:) :: xa, xb, xc

        double precision :: alpha,beta

        continue

! here, performs matrix multiplications

        Select case (iopc)

          case (1) ! mamult
            call mamult (a, b, c, ndim, beta)
#ifdef GPU
          case (2) ! mamult_gpu
            igrid = 512 ; iblock = 512
            tt = 0.0
            call mamult_gpu(a, b, c, ndim, mdim, ifact, beta, igrid, iblock, tt, 0)
#endif
          case (3) ! dgemm
            allocate (xa(ndim,ndim), xb(ndim,ndim), xc(ndim,ndim),stat=i)

!            forall (i=1:ndim,j=1:ndim)
!              xa(i,j) = 0.d0
!              xb(i,j) = 0.d0
!            endforall

            call dtpttr( 'u', ndim, a, xa, ndim, i )
            if (i /= 0) stop 'error in dtpttr'

            call dtpttr( 'u', ndim, b, xb, ndim, i )
            if (i /= 0) stop 'error in dtpttr'

            do i = 1,ndim-1
               call dcopy(ndim-i,xa(i,i+1),ndim,xa(i+1,i),1)
               call dcopy(ndim-i,xb(i,i+1),ndim,xb(i+1,i),1)
            end do

            if (.not.(beta == 0.d0)) then
!              forall (i=1:ndim,j=1:ndim)
!                xc(i,j) = 0.d0
!              endforall

              call dtpttr( 'u', ndim, c, xc, ndim, i )
              if (i /= 0) stop 'error in dtpttr'
              do i = 1,ndim-1
                 call dcopy(ndim-i,xc(i,i+1),ndim,xc(i+1,i),1)
              end do
            end if

            call dgemm ("N", "N", ndim, ndim, ndim, alpha, xa, ndim, xb, ndim, beta, xc, &
                       & ndim)

            call dtrttp('u', ndim, xc, ndim, c, i )

            deallocate (xa,xb,xc,stat=i)

#ifdef GPU
          case (4) ! dgemm_gpu

            allocate (xa(ndim,ndim), xb(ndim,ndim), xc(ndim,ndim),stat=i)

!            forall (i=1:ndim,j=1:ndim)
!              xa(i,j) = 0.d0
!              xb(i,j) = 0.d0
!            endforall

            call dtpttr( 'u', ndim, a, xa, ndim, i )
            if (i /= 0) stop 'error in dtpttr'

            call dtpttr( 'u', ndim, b, xb, ndim, i )
            if (i /= 0) stop 'error in dtpttr'

            do i = 1,ndim-1
               call dcopy(ndim-i,xa(i,i+1),ndim,xa(i+1,i),1)
               call dcopy(ndim-i,xb(i,i+1),ndim,xb(i+1,i),1)
            end do

            if (beta /= 0.d0) then
!              forall (i=1:ndim,j=1:ndim)
!                xc(i,j) = 0.d0
!              endforall

              call dtpttr( 'u', ndim, c, xc, ndim, i )
              if (i /= 0) stop 'error in dtpttr'
              do i = 1,ndim-1
                 call dcopy(ndim-i,xc(i,i+1),ndim,xc(i+1,i),1)
              end do
            end if

!            if (ngpus > 1 .and. ndim > 100) then
!
!               call gemm_cublas_mgpu ("N", "N", ndim, ndim, ndim, alpha, xa, ndim, xb, ndim, beta, xc, &
!                           & ndim)
!
!            else
               call gemm_cublas ("N", "N", ndim, ndim, ndim, alpha, xa, ndim, xb, ndim, beta, xc, &
                           & ndim)

!            end if

            call dtrttp('u', ndim, xc, ndim, c, i )

            deallocate (xa,xb,xc,stat=i)
#endif
        end select

        continue

      return
      end subroutine mult_symm_AB

! ---------------------------------------------------------------------------------------
