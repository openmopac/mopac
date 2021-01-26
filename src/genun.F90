      subroutine genun(u, n) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!****************************************************************
!
!     GENERATE UNIT VECTORS OVER SPHERE. USED BY SURFAC ONLY.
!
!****************************************************************
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(inout) :: n 
      double precision , intent(out) :: u(3,n) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nequat, nvert, nu, i, nhor, j 
      double precision :: pi, fi, z, xy, fj, x, y 
!-----------------------------------------------
      pi = 4.D0*atan(1.D0) 
      nequat = int(sqrt(n*pi)) 
      nvert = nequat/2 
      nu = 0 
      l20: do i = 1, nvert + 1 
        fi = (pi*(i - 1))/nvert 
        z = cos(fi) 
        xy = sin(fi) 
        nhor = int(nequat*xy) 
        nhor = max0(1,nhor) 
        do j = 1, nhor 
          fj = (2.D0*pi*(j - 1))/nhor 
          x = dcos(fj)*xy 
          y = dsin(fj)*xy 
          if (nu >= n) exit  l20 
          nu = nu + 1 
          u(1,nu) = x 
          u(2,nu) = y 
          u(3,nu) = z 
        end do 
      end do l20 
      n = nu 
      return  
      end subroutine genun 


      logical function collid (cw, rw, cnbr, rnbr, nnbr, ishape) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!****************************************************************
!
!     COLLISION CHECK OF PROBE WITH NEIGHBORING ATOMS
!     USED BY SURFAC ONLY.
!
!****************************************************************
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nnbr 
      integer , intent(in) :: ishape 
      double precision , intent(in) :: rw 
      double precision , intent(in) :: cw(3) 
      double precision , intent(in) :: cnbr(3,200) 
      double precision , intent(in) :: rnbr(200) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i 
      double precision :: sumrad, vect1, vect2, vect3, sr2, dd2 
!-----------------------------------------------
      if (nnbr > 0) then 
!
!     CHECK WHETHER PROBE IS TOO CLOSE TO ANY NEIGHBOR
!
        if (ishape == 3) then 
        else 
          do i = 1, nnbr 
            sumrad = rw + rnbr(i) 
            vect1 = dabs(cw(1)-cnbr(1,i)) 
            if (vect1 >= sumrad) cycle  
            vect2 = dabs(cw(2)-cnbr(2,i)) 
            if (vect2 >= sumrad) cycle  
            vect3 = dabs(cw(3)-cnbr(3,i)) 
            if (vect3 >= sumrad) cycle  
            sr2 = sumrad**2 
            dd2 = vect1**2 + vect2**2 + vect3**2 
            if (dd2 < sr2) go to 30 
          end do 
        endif 
      endif 
      collid = .FALSE. 
      go to 40 
   30 continue 
      collid = .TRUE. 
   40 continue 
      return  
      end function collid 
