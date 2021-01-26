      subroutine dfield() 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE parameters_C, only : tore  
      USE molkst_C, only :numat, efield
      USE funcon_C, only : ev, a0, fpc_9
      use common_arrays_C, only : nat, p, dxyz
!***********************************************************************
!DECK MOPAC
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use chrge_I 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!----------------------------------------------- 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      double precision, dimension(numat) :: q2 
      double precision :: fldcon 
!-----------------------------------------------
      call chrge (p, q2) 
      q2(:numat) = tore(nat(:numat)) - q2(:numat) 
!
!   FLDCON=(h/Ao)*(eV to Kcal/mol)
!
      fldcon = ev/a0*fpc_9 
      do i = 1, numat 
        dxyz(3*(i-1) + 1) = dxyz(3*(i-1) + 1) + efield(1)*q2(i)*fldcon 
        dxyz(3*(i-1) + 2) = dxyz(3*(i-1) + 2) + efield(2)*q2(i)*fldcon 
        dxyz(3*(i-1) + 3) = dxyz(3*(i-1) + 3) + efield(3)*q2(i)*fldcon 
      end do 
      return  
      end subroutine dfield 
