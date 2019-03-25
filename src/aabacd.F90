      real(kind(0.0d0)) function aabacd (iocca1, ioccb1, iocca2, ioccb2, nmos, xy) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  10:46:59  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nmos 
      integer , intent(in) :: iocca1(nmos) 
      integer , intent(in) :: ioccb1(nmos) 
      integer , intent(in) :: iocca2(nmos) 
      integer , intent(in) :: ioccb2(nmos) 
      real(double) , intent(in) :: xy(nmos,nmos,nmos,nmos) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ij, i, j, k, l 
      real(double) :: sum 
!-----------------------------------------------
!**********************************************************************
!
! AABACD EVALUATES THE C.I. MATRIX ELEMENT FOR TWO MICROSTATES DIFFERING
!       BY TWO ALPHA MOS. ONE MICROSTATE HAS ALPHA ELECTRONS IN
!       M.O.S PSI(I) AND PSI(J) FOR WHICH THE OTHER MICROSTATE HAS
!       ELECTRONS IN PSI(K) AND PSI(L)
!
!**********************************************************************
      ij = 0 
      do i = 1, nmos 
        if (iocca1(i) >= iocca2(i)) cycle  
        exit  
      end do 
      do j = i + 1, nmos 
        if (iocca1(j) < iocca2(j)) exit  
        ij = ij + iocca2(j) + ioccb2(j) 
      end do 
      do k = 1, nmos 
        if (iocca1(k) <= iocca2(k)) cycle  
        exit  
      end do 
      do l = k + 1, nmos 
        if (iocca1(l) > iocca2(l)) exit  
        ij = ij + iocca1(l) + ioccb1(l) 
      end do 
      ij = ij + ioccb2(i) + ioccb1(k) 
      sum = xy(i,k,j,l) - xy(i,l,k,j) 
      if (mod(ij,2) == 1) sum = -sum 
      aabacd = sum 
      return  
      end function aabacd 
