      subroutine bfn(x, bf) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      use overlaps_C, only : fact
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision , intent(in) :: x 
      double precision , intent(out) :: bf(13) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, io, last, i, m 
      double precision :: absx, expx, expmx, y, xf 
!-----------------------------------------------
!**********************************************************************
!
!     BINTGS FORMS THE "B" INTEGRALS FOR THE OVERLAP CALCULATION.
!
!**********************************************************************
      k = 12 
      io = 0 
      absx = abs(x) 
      if (absx <= 3.D00) then 
        if (absx > 2.D00) then 
          last = 15 
          go to 60 
        end if 
        if (absx > 1.D00) then 
          last = 12 
          go to 60 
        end if 
        if (absx > 0.5D00) then 
          last = 7 
          go to 60 
        end if 
        if (absx <= 1.D-6) go to 90 
        last = 6 
        go to 60 
      end if 
      expx = exp(x) 
      expmx = 1.D00/expx 
      bf(1) = (expx - expmx)/x 
      do i = 1, k 
        bf(i+1) = (i*bf(i)+(-1.D00)**i*expx-expmx)/x 
      end do 
      go to 110 
   60 continue 
      do i = io, k 
        y = 0.0D00 
        do m = io, last 
          xf = 1.0D00 
          if (m /= 0) xf = fact(m) 
          y = y + (-x)**m*(2*mod(m + i + 1,2))/(xf*(m + i + 1)) 
        end do 
        bf(i+1) = y 
      end do 
      go to 110 
   90 continue 
      do i = io, k 
        bf(i+1) = (2*mod(i + 1,2))/(i + 1.D0) 
      end do 
  110 continue 
      return  
!
      end subroutine bfn 
