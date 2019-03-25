
module density_cuda_i
    interface density_cuda
        subroutine density_gpu(c, norbs, p, nn, nl1, nl2, nu1, nu2, cst, frac, sign, xtime, &
            & nThreads,nblocks) bind(c, &
            & name='density_cuda_driver')
            use iso_c_binding
            
            implicit none
            
            integer(c_int), value :: norbs, nn, nl1, nl2, nu1, nu2, nThreads,nblocks
            
            real(c_float) :: xtime
            real(c_double), dimension(nn) :: p
            real(c_double), dimension(norbs, norbs) :: c
            real(c_double), value :: cst, frac, sign
            
        end subroutine
    end interface
end module  

module call_rot_cuda
    interface diag2_cuda
      subroutine rot_cuda(fmo,eig,vector,ci0,ca0,nocc,lumo,n,bigeps,tiny) &
                         & bind(c, name='diag2GPU_Driver')
    	use iso_c_binding
    	implicit none
    	integer(c_int), value  :: nocc, n,lumo
	    real(c_double)  :: fmo(n*n), eig(n), vector(n,n),ci0(n,nocc),ca0(n,n-nocc)
	    real(c_double), value :: bigeps, tiny
	  end subroutine    
    end interface 

    interface diag2_cuda_2gpu
      subroutine rot_cuda_2gpu(fmo,eig,vector,ci0,ca0,nocc,lumo,n,bigeps,tiny) &
                         & bind(c, name='diag2GPU_Driver_2gpu')
        use iso_c_binding
        implicit none
        integer(c_int), value  :: nocc, n,lumo
        real(c_double)  :: fmo(n*n), eig(n), vector(n,n),ci0(n,nocc),ca0(n,n-nocc)
        real(c_double), value :: bigeps, tiny
      end subroutine
    end interface
end module

