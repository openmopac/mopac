      subroutine prttim(tleft, tprt, txt) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision , intent(in) :: tleft 
      double precision , intent(out) :: tprt 
      character , intent(out) :: txt 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      tprt = tleft 
      txt = 'S' 
      if (tprt >= 604800D0) then 
        tprt = tprt/604800.D0 
        txt = 'W' 
      else if (tprt >= 86400.D0) then 
        tprt = tprt/86400.D0 
        txt = 'D' 
      else if (tprt >= 3600.D0) then 
        tprt = tprt/3600.D0 
        txt = 'H' 
      else if (tprt >= 60.D0) then 
        tprt = tprt/60.D0 
        txt = 'M' 
      end if 
      return  
      end subroutine prttim 
