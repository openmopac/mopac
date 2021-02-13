      MODULE spline_I   
      INTERFACE
      subroutine spline(x, f, df, xhigh, xlow,  xmin, npnts)
      integer :: npnts
      double precision, dimension(12) :: x, f, df 
      double precision :: xlow, xhigh, xmin, fmin, dfmin
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
