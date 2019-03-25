module mamult_cuda_i 
    interface mamult_cuda
        subroutine mamult_gpu(a, b, c, n, nn, ipoint, cst, gridx, blockx, tempo, job) bind(c, name='mamult_driver')
            use iso_c_binding
            implicit none
            integer(c_int), value :: n, nn, gridx, blockx, job
            real(c_float) :: tempo
#if CC12
            real(c_float), dimension(nn) :: a, b, c
            real(c_float),value :: cst           
#else
            real(c_double), dimension(nn) :: a, b, c
            real(c_double),value :: cst   
#endif
            integer(c_int), dimension(n) :: ipoint
        end subroutine
    end interface

end module
