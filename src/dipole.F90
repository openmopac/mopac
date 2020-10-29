      real(kind(0.0d0)) function dipole (p, coord, dipvec, mode) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double  
      USE molkst_C, only : numat, keywrd, numcal, ux, uy, uz, mpack, &
      mol_weight
      USE parameters_C, only : ddp, dd
      use to_screen_C, only : dip
      USE funcon_C, only : fpc_8, fpc_1, a0
      USE chanel_C, only : iw
      use common_arrays_C, only : nfirst, nlast, nat, atmass, q
!...Translated by Pacific-Sierra Research 77to90  4.4G  09:21:03  03/16/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: mode 
      real(double) , intent(in) :: p(mpack) 
      real(double) , intent(inout) :: coord(3,*) 
      real(double) , intent(out) :: dipvec(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icalcn, i, j, ni, ia, ib, l, k, ll 
      real(double), dimension(3) :: center 
      real(double) :: sum, hyfsp, xt, hyfpd, dx, dy, dz, factdm 
      logical :: first, chargd, force 

      save first, chargd, force, icalcn 
!-----------------------------------------------
!***********************************************************************
!     DIPOLE CALCULATES DIPOLE MOMENTS
!
!  ON INPUT P     = DENSITY MATRIX
!           Q     = TOTAL ATOMIC CHARGES, (NUCLEAR + ELECTRONIC)
!           NUMAT = NUMBER OF ATOMS IN MOLECULE
!           NAT   = ATOMIC NUMBERS OF ATOMS
!           NFIRST= START OF ATOM ORBITAL COUNTERS
!           COORD = COORDINATES OF ATOMS
!
!  OUTPUT  DIPOLE = DIPOLE MOMENT
!***********************************************************************
!
!     IN THE ZDO APPROXIMATION, ONLY TWO TERMS ARE RETAINED IN THE
!     CALCULATION OF DIPOLE MOMENTS.
!     1. THE POINT CHARGE TERM (INDEPENDENT OF PARAMETERIZATION).
!     2. THE ONE-CENTER HYBRIDIZATION TERM, WHICH ARISES FROM MATRIX
!     ELEMENTS OF THE FORM <NS/R/NP>. THIS TERM IS A FUNCTION OF
!     THE SLATER EXPONENTS (ZS,ZP) AND IS THUS DEPENDENT ON PARAMETER-
!     IZATION.
!
!
!     REFERENCES:
!     J.A.POPLE & D.L.BEVERIDGE: APPROXIMATE M.O. THEORY
!     S.P.MCGLYNN, ET AL: APPLIED QUANTUM CHEMISTRY
!  
      data icalcn/ 0/  
      first = icalcn /= numcal 
      icalcn = numcal 
!
!   CONST = c.C.a0.2
!
      if (first) then 
        sum = 0.D0 
        do i = 1, numat  
          sum = sum + q(i) 
        end do 
        chargd = abs(sum) > 0.5D0 
        force = index(keywrd,'FORCE') + index(keywrd," THERMO") + index(keywrd,'IRC') /= 0 
      endif 
      if (.not.force .and. chargd) then 
!
!   NEED TO RESET ION'S POSITION SO THAT THE CENTER OF MASS IS AT THE
!   ORIGIN.
!
        center = 0.D0 
        do i = 1, 3 
          do j = 1, numat 
            center(i) = center(i) + atmass(j)*coord(i,j) 
          end do 
        end do 
        center = center/mol_weight 
        do i = 1, 3 
          coord(i,:numat) = coord(i,:numat) - center(i) 
        end do 
      endif 
      dip = 0.0D00 
      do i = 1, numat 
        ni = nat(i) 
        ia = nfirst(i) 
        ib = nlast(i) 
        l = ib - ia 
        if (l > 0) then  
            hyfsp = 2.0D0*dd(ni)*a0*fpc_8*fpc_1*1.D-10 
          do j = 1, 3 
            k = ((ia + j)*(ia + j - 1))/2 + ia 
            dip(j,2) = dip(j,2) - hyfsp*p(k) 
          end do 
!
!   PD ONE-CENTER TERM
          if (ib - ia == 8) then 
            xt = 1.D0/sqrt(3.D00) 
            hyfpd = 2.0D0*ddp(5,ni)*a0*fpc_8*fpc_1*1.D-10 
            ll = (ia + 5)*(ia + 4)/2 + ia + 3 
            dx = p(ll) 
            ll = (ia + 4)*(ia + 3)/2 + ia + 1 
            dx = dx + p(ll) 
            ll = (ia + 8)*(ia + 7)/2 + ia + 2 
            dx = dx + p(ll) 
            ll = (ia + 6)*(ia + 5)/2 + ia + 1 
            dx = dx - xt*p(ll) 
!
            ll = (ia + 7)*(ia + 6)/2 + ia + 3 
            dy = p(ll) 
            ll = (ia + 4)*(ia + 3)/2 + ia + 2 
            dy = dy - p(ll) 
            ll = (ia + 8)*(ia + 7)/2 + ia + 1 
            dy = dy + p(ll) 
            ll = (ia + 6)*(ia + 5)/2 + ia + 2 
            dy = dy - xt*p(ll) 
!
            ll = (ia + 5)*(ia + 4)/2 + ia + 1 
            dz = p(ll) 
            ll = (ia + 7)*(ia + 6)/2 + ia + 2 
            dz = dz + p(ll) 
            ll = (ia + 6)*(ia + 5)/2 + ia + 3 
            dz = dz + 2.D0*xt*p(ll) 
            dip(1,2) = dip(1,2) - dx*hyfpd 
            dip(2,2) = dip(2,2) - dy*hyfpd 
            dip(3,2) = dip(3,2) - dz*hyfpd 
          endif 
        endif 
!
!  FPC(8)=SPEED OF LIGHT,  FPC(1)=CHARGE ON ELECTRON.
!
        factdm = fpc_8*fpc_1*1.D-10 
        dip(:3,1) = dip(:3,1) + q(i)*coord(:,i)*factdm 
      end do 
      dip(:3,3) = dip(:3,2) + dip(:3,1) 
      do j = 1, 3 
        dip(4,j) = sqrt(dip(1,j)**2+dip(2,j)**2+dip(3,j)**2) 
      end do 
      if (force) then 
        dipvec(1) = dip(1,3) 
        dipvec(2) = dip(2,3) 
        dipvec(3) = dip(3,3) 
      endif 
      if (mode == 1) write (iw, 140) ((dip(i,j),i=1,4),j=1,3)
!     STORE DIPOLE MOMENT COMPONENTS IN UX,UY,UZ FOR USE IN
!     ASSIGNING CHARGES DETERMINED FROM THE ESP. BHB
      ux = dip(1,3) 
      uy = dip(2,3) 
      uz = dip(3,3) 
      dipole = dip(4,3) 
      if (.not.force .and. chargd) then 
!
!   NEED TO RESET ION'S POSITION AGAIN, SO THAT IT IS BACK WHERE IT STARTED
!
        do i = 1, 3 
          coord(i,:numat) = coord(i,:numat) + center(i) 
        end do 
      endif 
      return  
!
  140 format(' DIPOLE',11x,'X         Y         Z       TOTAL',/,' POINT-CHG.',&
        4f10.3,/,' HYBRID',4x,4f10.3,/,' SUM',7x,4f10.3)  
!
      end function dipole 
