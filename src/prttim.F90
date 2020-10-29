      subroutine prttim(tleft, tprt, txt) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:35:39  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double) , intent(in) :: tleft 
      real(double) , intent(out) :: tprt 
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
      endif 
      return  
      end subroutine prttim 
