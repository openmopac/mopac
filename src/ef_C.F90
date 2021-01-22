      module ef_C 
      USE vast_kind_param, ONLY:  double 
      integer :: ef_mode, nstep, negreq, iprnt, iloop
      real(double) :: ddx, rmin, rmax, omin = 0.d0, xlamd, xlamd0, skal, x0, x1, x2
      real(double), dimension(:), allocatable :: oldf, d, vmode, pmat, uc, hessc
      real(double), dimension(:,:), allocatable :: hess, bmat, u, oldhss, oldu, alparm
      end module ef_C 
      module drc_C 
        USE vast_kind_param, ONLY:  double 
        real(double), dimension(:), allocatable :: vref, vref0, now  
        real(double), dimension(:,:), allocatable :: allxyz, allvel, xyz3, vel3, allgeo, geo3 
        real(double), dimension(:), allocatable :: parref
        double precision :: time
      end module drc_C 
      module derivs_C
        USE vast_kind_param, ONLY:  double 
        real(double), dimension(:), allocatable :: wmat, hmat, fmat
        real(double), dimension(:,:), allocatable  :: b, ab, fb
          real(double), dimension(:), allocatable :: aidref, work2
      end module derivs_C

