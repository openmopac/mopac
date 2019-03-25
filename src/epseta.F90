subroutine epseta (eps, eta)
   !
   !.. Implicit Declarations ..
    implicit none
   !
   !.. Formal Arguments ..
    double precision, intent (out) :: eps, eta
   !
   !.. Intrinsic Functions ..
    intrinsic Max, Tiny, Epsilon

   ! ... Executable Statements ...
   !
   !     RETURN ETA, THE SMALLEST REPRESENTABLE NUMBER,
   !     AND EPS, THE SMALLEST NUMBER FOR WHICH 1+EPS.NE.1.
   !
    eps = Max (Tiny (1.d0), 1.d-39)
    eta = Max (Epsilon (1.d0), 1.d-17)
end subroutine epseta
