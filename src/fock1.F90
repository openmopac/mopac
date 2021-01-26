      subroutine fock1(f, ptot, pa, pb) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      use common_arrays_C, only : nat, nlast, nfirst
      use molkst_C, only : numat
      USE parameters_C, only : gpp, gp2, hsp, gss, gsp
!***********************************************************************
!DECK MOPAC
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision , intent(inout) :: f(*) 
      double precision  :: ptot(*) 
      double precision  :: pa(*) 
      double precision , intent(in) :: pb(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ii, ia, ib, ni, nbases, ka, l, j, m 
      double precision :: ptpop, papop 
!-----------------------------------------------
! *********************************************************************
!
! *** COMPUTE THE REMAINING CONTRIBUTIONS TO THE ONE-CENTER ELEMENTS.
!
! *********************************************************************
      do ii = 1, numat 
        ia = nfirst(ii) 
        ib = nlast(ii) 
        ni = nat(ii) 
        ptpop = 0.D0 
        papop = 0.D0 
        nbases = ib - ia + 1 
        go to (1002,20,10,10,10,10,10,10,10,10) nbases + 1 
   10   continue 
        ptpop = ptot((ib*(ib+1))/2) + ptot(((ib-1)*ib)/2) + ptot(((ib-2)*(ib-1)&
          )/2) 
        papop = pa((ib*(ib+1))/2) + pa(((ib-1)*ib)/2) + pa(((ib-2)*(ib-1))/2) 
   20   continue 
!
!     F(S,S)
!
        ka = (ia*(ia + 1))/2 
        f(ka) = f(ka) + pb(ka)*gss(ni) + ptpop*gsp(ni) - papop*hsp(ni) 
        if (nbases /= 1) then 
          l = ka 
          do j = ia + 1, ib 
            m = l + ia 
            l = l + j 
!
!     F(P,P)
!
            f(l) = f(l) + ptot(ka)*gsp(ni) - pa(ka)*hsp(ni) + pb(l)*gpp(ni) + (&
              ptpop - ptot(l))*gp2(ni) - 0.5D0*(papop - pa(l))*(gpp(ni)-gp2(ni)&
              ) 
!
!     F(S,P)
!
            f(m) = f(m) + 2.D0*ptot(m)*hsp(ni) - pa(m)*(hsp(ni)+gsp(ni)) 
          end do 
!
!     F(P,P*)
!
          do j = ia + 1, ib - 1 
            do l = j + 1, ib 
              m = (l*(l - 1))/2 + j 
              f(m) = f(m) + ptot(m)*(gpp(ni)-gp2(ni)) - 0.5D0*pa(m)*(gpp(ni)+&
                gp2(ni)) 
            end do 
          end do 
        endif 
 1002   continue 
      end do 
      return  
      end subroutine fock1 
