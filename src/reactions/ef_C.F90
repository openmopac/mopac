      module ef_C 
      integer :: ef_mode, nstep, negreq, iprnt, iloop
      double precision :: ddx, rmin, rmax, omin = 0.d0, xlamd, xlamd0, skal, x0, x1, x2
      double precision, dimension(:), allocatable :: oldf, d, vmode, pmat, uc, hessc
      double precision, dimension(:,:), allocatable :: hess, bmat, u, oldhss, oldu, alparm
      end module ef_C 
      module drc_C 
        double precision, dimension(:), allocatable :: vref, vref0, now  
        double precision, dimension(:,:), allocatable :: allxyz, allvel, xyz3, vel3, allgeo, geo3 
        double precision, dimension(:), allocatable :: parref
        double precision :: time
      end module drc_C 
      module derivs_C
        double precision, dimension(:), allocatable :: wmat, hmat, fmat
        double precision, dimension(:,:), allocatable  :: b, ab, fb
          double precision, dimension(:), allocatable :: aidref, work2
      end module derivs_C

