      MODULE spline_I   
      INTERFACE
      subroutine spline(x, f, df, xhigh, xlow,  xmin, npnts)
      USE vast_kind_param, ONLY:  double 
      integer :: npnts
      real(double), dimension(12) :: x, f, df 
      real(double) :: xlow, xhigh, xmin, fmin, dfmin
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
