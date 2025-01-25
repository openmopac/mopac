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

subroutine density_for_GPU (c, fract, ndubl, nsingl, occ, mpack, norbs, mode, pp, iopc)
#ifdef GPU
      Use mod_vars_cuda, only: real_cuda, prec, nthreads_gpu, nblocks_gpu
      Use iso_c_binding
      Use density_cuda_i
      Use call_gemm_cublas
      Use call_syrk_cublas
#endif
      implicit none
      Integer :: ndubl, nsingl, mode, mpack, norbs, nl1, nl2, nu1, nu2, i, j, l, &
	           & nl21, nl11, iopc
      double precision,allocatable :: xmat(:,:)
      double precision :: c(norbs,norbs), pp(mpack)
      double precision :: cst, sign, fract, frac, occ, sum1, sum2
#ifdef GPU
      double precision, allocatable :: pdens(:)
#endif
      if (ndubl /= 0 .and. nsingl > (norbs/2) .and. mode /= 2) then
        !
        !    TAKE POSITRON EQUIVALENT
        !
        sign = -1.d0
        frac = occ - fract
        cst = occ
        nl2 = nsingl + 1
        nu2 = norbs
        nl1 = ndubl + 1
        nu1 = nsingl
      else
        !
        !    TAKE ELECTRON EQUIVALENT
        !
        sign = 1.d0
        frac = fract
        cst = 0.d0
        nl2 = 1
        nu2 = ndubl
        nl1 = ndubl + 1
        nu1 = nsingl
      end if
      Select case (iopc)
        case(2)   ! Option to use dgemm from CUBLAS
#ifdef GPU

          nl21 = Min (norbs, nl2)
          nl11 = Min (norbs, nl1)
          allocate(xmat(norbs,norbs),stat = i)
          call gemm_cublas ('N', 'T', norbs, norbs, nu2-nl2+1, 2.0_prec*sign, c(1:norbs,nl21:norbs),&
                        &   norbs, c(1:norbs,nl21:norbs), norbs, 0.0_prec, xmat, norbs)
          call gemm_cublas ('N', 'T', norbs, norbs, nu1-nl1+1, frac*sign, c(1:norbs,nl11:norbs), &
                        &   norbs, c(1:norbs,nl11:norbs), norbs, 1.0_prec, xmat, norbs)
          forall (i=1:norbs)
             xmat(i,i) = xmat(i,i) + cst
          endforall
          call dtrttp('u', norbs, xmat, norbs, pp, i )
          deallocate (xmat,stat=i)
#endif
        case(3)   ! Option to use dgemm from BLAS

          nl21 = Min (norbs, nl2)
          nl11 = Min (norbs, nl1)

          allocate(xmat(norbs,norbs),stat = i)
          if (norbs < 0) forall (j=1:norbs,i=1:norbs) xmat(i,j) = 0.d0  ! Dummy statement to 'fool' FORTRAN checks

          call dgemm ('N', 'T', norbs, norbs, nu2-nl2+1, 2.0d0*sign, c(1:norbs,nl21:norbs),&
                        &   norbs, c(1:norbs,nl21:norbs), norbs, 0.0d0, xmat, norbs)
          call dgemm ('N', 'T', norbs, norbs, nu1-nl1+1, frac*sign, c(1:norbs,nl11:norbs), &
                        &   norbs, c(1:norbs,nl11:norbs), norbs, 1.0d0, xmat, norbs)

          forall (i=1:norbs)
             xmat(i,i) = xmat(i,i) + cst
          endforall

          call dtrttp('u', norbs, xmat, norbs, pp, i )

          deallocate (xmat,stat=i)
        case(4)   ! Option to use dsyrk from CUBLAS
#ifdef GPU
          allocate(xmat(norbs,norbs),stat = i)
          forall (j = 1:norbs, i=1:norbs) xmat(i, j) = 0.d0
          call syrk_cublas ('U','N',norbs,ndubl, &
                 & occ,c(1:norbs,1:ndubl),norbs, &
                 & 0.d0,xmat,norbs)
          call dtrttp('u', norbs, xmat, norbs, pp, i )
          deallocate(xmat,stat=i)
#endif
        case(5)   ! Option to use dsyrk from BLAS
          if (fract < 1.d-2) then
            allocate(xmat(norbs,norbs),stat = i)
            forall (j = 1:norbs, i=1:norbs) xmat(i, j) = 0.d0
            call dsyrk ('u', 'n', norbs, ndubl, occ, c(1:norbs,1:ndubl), norbs, 0.d0, xmat, norbs) ! For RHF	      	      	
            call dtrttp('u', norbs, xmat, norbs, pp, i )
            deallocate (xmat,stat=i)
          else
!
! The following block should be re-cast in a modern style "someday"
! It's used only infrequently, so updating it is not urgent.
!
            l = 0
            do i = 1, norbs
              do j = 1, i
                l = l + 1
                sum1 = 0.D0
                sum2 = sum(c(i,nl2:nu2)*c(j,nl2:nu2))
                sum2 = sum2*occ
                sum1 = sum(c(i,nl1:nu1)*c(j,nl1:nu1))
                pp(l) = (sum2 + sum1*frac)*sign
              end do
              pp(l) = cst + pp(l)
            end do
          end if
    End select
    continue
    return
End subroutine density_for_GPU
