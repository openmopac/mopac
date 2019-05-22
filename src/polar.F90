      subroutine polar() 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use common_arrays_C, only : nat, na, labels, &
      & geo, coord, w, c, tvec
      USE parameters_C, only : polvol
      USE funcon_C, only : ev
      USE elemts_C, only : elemnt 
      USE molkst_C, ONLY: numat, norbs, natoms, ndep, nvar, keywrd, last, &
      limscf, id, moperr
      USE polar_C, only : omega 
      USE chanel_C, only : iw
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use gmetry_I 
      use axis_I 
      use compfg_I 
      use wrdkey_I 
      use reada_I 
      use openda_I 
      use alphaf_I 
      use betaf_I 
      use nonor_I 
      use nonope_I 
      use nonbet_I 
      use beopor_I 
      use ngamtg_I 
      use ngefis_I 
      use ngidri_I 
      use ngoke_I 
      implicit none
      real(double), allocatable, dimension (:,:) :: x1, x2, x3, x4, x5, x6, x7, x8, &
        x19, x10, x11, x12, x13
      integer :: i, j, k, l, iwflb, ibet, igam, maxitu, maxita, &
        nfreq, iwfla, ref_na(natoms)
      real(double), dimension(nvar) :: grad 
      real(double), dimension(3,3) :: rotvec, tempv 
      real(double), dimension(10) :: dataev 
      real(double), dimension(3,natoms) :: refgeo 
      real(double) :: a, b, sum, sumx, sumy, sumz, summax, atpol, heat0, &
        atol, btol, hartr, omegau, cc
      logical :: let 
      character :: polkey*241 
!-----------------------------------------------
!
!   THIS COMPUTER PROGRAM HAS BEEN PLACED IN THE PUBLIC DOMAIN
!   BY THE AUTHORS PROF HENRY KURTZ AND PRAKASHAN KORAMBATH
!   DEPARTMENT OF CHEMISTRY, MEMPHIS STATE UNIVERSITY, MEMPHIS TN
!   1992
!.. 6/13/91
!**********************************************************************
!
!   POLAR SETS UP THE CALCULATION OF THE MOLECULAR ELECTRIC RESPONSE
!   PROPERTIES BY FFHPOL.
!
!**********************************************************************
!..
    allocate(x1(norbs,norbs), x2(norbs,norbs), x3(norbs,norbs), x4(norbs,norbs), &
      x5(norbs,norbs), x6(norbs,norbs), x7(norbs,norbs), x8(norbs,norbs), &
      x19(norbs,norbs), x10(norbs,norbs), x11(norbs,norbs), x12(norbs,norbs), x13(norbs,norbs), &
      stat=i)
    if (i /= 0) then
      call memory_error ("polar")
      return
    end if
    x1 = 0.d0
    x2 = 0.d0
    x3 = 0.d0
    x4 = 0.d0
    x5 = 0.d0
    x6 = 0.d0
    x7 = 0.d0
    x8 = 0.d0
    x19 = 0.d0
    x10 = 0.d0
    x11 = 0.d0
      limscf = .FALSE. 
      let = index(keywrd,'LET') /= 0 
      write (iw, 10) 
   10 format(' ',20('*'),' TDHF POLARIZABILITIES ',20('*'),/,/) 
   refgeo(:,:natoms) = geo(:,:natoms) 
   ref_na(:natoms) = na(:natoms)
      call gmetry (geo, coord) 
!
!  ORIENT THE MOLECULE WITH THE MOMENTS OF INERTIA.
!  THIS IS DONE TO ENSURE A UNIQUE, REPRODUCEABLE SET OF DIRECTIONS.
!  IF LET IS SPECIFIED, THE INPUT ORIENTATION WILL BE USED.
!
      if (.not.let) then 
        call axis (a, b, cc, rotvec) 
        if (Int(a + b + cc) == -1000) return ! dummy use of arguments of axis
        write (iw, 20) 
   20   format(/,' ROTATION MATRIX FOR ORIENTATION OF MOLECULE:'/) 
        do i = 1, 3 
          write (iw, 30) (rotvec(i,j),j=1,3) 
   30     format(5x,3f12.6) 
        end do 
!
!  ROTATE ATOMS
!
        do i = 1, numat 
          do j = 1, 3 
            sum = 0.0D00 
            do k = 1, 3 
              sum = sum + coord(k,i)*rotvec(k,j) 
            end do 
            geo(j,i) = sum 
          end do 
        end do 
        coord(:,:numat) = geo(:,:numat) 
        write (iw, '(2/10X,''CARTESIAN COORDINATES '',/)') 
        write (iw, &
      '(4X,''NO.'',7X,''ATOM'',9X,''X'',                      9X,''Y'',9X,''Z''&
      &,/)') 
        l = 0 
        do i = 1, numat 
          if (nat(i)==99 .or. nat(i)==107) cycle  
          l = l + 1 
          write (iw, '(I6,8X,A2,4X,3F10.4)') l, elemnt(nat(i)), (coord(j,l),j=1&
            ,3) 
        end do 
!
!  IF POLYMER, ROTATE TVEC
!  (BEWARE:  THE POLYMER SECTIONS MAY NOT WORK YET)
!
        if (id > 0) then 
          do i = 1, id 
            do j = 1, 3 
              sum = 0.0D00 
              do k = 1, 3 
                sum = sum + tvec(k,i)*rotvec(k,j) 
              end do 
              tempv(j,i) = sum 
            end do 
          end do 
          tvec(:,:id) = tempv(:,:id) 
          write (iw, 160) ((tvec(j,i),j=1,3),i=1,id) 
  160     format(/,' NEW TRANSLATION VECTOR:'/,' ',3(3f15.5)) 
        endif 
      endif 
!
      last = 1 
      na = 0
!
!  SET UP THE VARIABLES IN XPARAM AND LOC, THESE ARE IN CARTESIAN
!  COORDINATES.
!
      ndep = 0 
      numat = 0 
      sumx = 0.D0 
      sumy = 0.D0 
      sumz = 0.D0 
      do i = 1, natoms 
        if (labels(i)==99 .or. labels(i)==107) cycle  
        numat = numat + 1 
        labels(numat) = labels(i) 
        sumx = sumx + coord(1,numat) 
        sumy = sumy + coord(2,numat) 
        sumz = sumz + coord(3,numat) 
        geo(:,numat) = coord(:,numat) 
      end do 
      sumx = sumx/numat 
      sumy = sumy/numat 
      sumz = sumz/numat 
      summax = 0.D0 
      atpol = 0.D0 
      do i = 1, numat 
        if (labels(i) /= 107) atpol = atpol + polvol(labels(i)) 
        geo(1,i) = geo(1,i) - sumx 
        summax = dmax1(abs(geo(1,i)),summax) 
        geo(2,i) = geo(2,i) - sumy 
        summax = dmax1(abs(geo(2,i)),summax) 
        geo(3,i) = geo(3,i) - sumz 
        summax = dmax1(abs(geo(3,i)),summax) 
      end do 
!
      nvar = 0 
      natoms = numat 
      call compfg (geo, .TRUE., heat0, .TRUE., grad, .FALSE.) 
      write (iw, 200) heat0 
  200 format(/,/,' ENERGY OF "REORIENTED" SYSTEM WITHOUT FIELD:',f20.10) 
!...............................................................
!
!  VARIABLES USED FOR TIME-DEPENDENT CALCULATIONS
!
!    OMEGA .........  FREQUENCY OF LIGHT (ACTUALLY INPUT AS ENERGY
!                     IN EV'S.
!    IWFLA .........  TYPE OF ALPHA CALCULATION FOR STORING MATRICES
!                     0 = STATIC
!                     1 = OMEGA
!                     2 = 2*OMEGA
!                     3 = 3*OMEGA
!    IWFLB .........  TYPE OF BETA CALCULATION FOR STORING MATRICES
!                     0 = (0,0)
!                     1 = (W,W) (SHG)
!                     2 = (0,W) (EOPE)
!                     3 = (W,-W) (OR)
!
!  INPUT NUMBER OF FREQENCIES TO RUN
!
!     IBET = 0  NO BETA CALC
!            1  ITERATIVE BETA
!           -1  NONITER BETA (SHG)
!           -2  NONITER EOPE
!           -3  NONITER OR
!
!     IGAM = 0  NO GAMMA CALC
!            1 THIRD HARMONIC GENERATION INPUT N,0,1,1
!            2 DC-EFISHG INPUT N,0,1,2
!            3 IDRI N,0,1,3
!            4 OKE N,0,1,4
!            5 DC EFIOR (NOT AVAILABLE)
!
      iwflb = nint(wrdkey(keywrd,'POLAR(',6,'IWFLB',5,0.D0)) 
      ibet = nint(wrdkey(keywrd,'POLAR(',6,'BETA',4,1.D0)) 
      igam = nint(wrdkey(keywrd,'POLAR(',6,'GAMMA',5,1.D0)) 
      atol = wrdkey(keywrd,'POLAR(',6,'TOL',3,0.001D0) 
      maxitu = nint(wrdkey(keywrd,'POLAR(',6,'MAXITU',6,500.D0)) 
      maxita = nint(wrdkey(keywrd,'POLAR(',6,'MAXITA',6,150.D0)) 
      btol = wrdkey(keywrd,'POLAR(',6,'BTOL',4,0.001D0) 
      i = index(keywrd,'POLAR(') 
      if (i /= 0) then 
        polkey = keywrd(i:) 
        i = index(polkey,') ') 
        polkey(i:) = ' ' 
      else 
        polkey = ' ' 
      endif 
      i = index(polkey,'E=(') 
      if (i /= 0) then 
        j = index(polkey(i:),')') 
        polkey(:i) = ' ' 
        polkey(i+j:) = ' ' 
        do i = 1, 10 
          dataev(i) = reada(polkey,1) 
          j = index(polkey,',') 
          if (j /= 0) then 
            polkey(:j) = ' ' 
          else 
            nfreq = i 
            go to 230 
          endif 
        end do 
      else 
        do i = 1, 3 
          dataev(i) = (i - 1)*0.25D0 
        end do 
        nfreq = 3 
      endif 
  230 continue 
      if (igam /= 0) ibet = 1 
      write (iw, 240) nfreq, iwflb, ibet, igam 
  240 format(/,/,'  NFREQ=',i3,'  IWFLB=',i3,'  IBET=',i3,'  IGAM=',i3) 
!
! ATOL IS THE MAXIMUM TOLERANCE IN MAKEUF AND BTOL IS THAT IN BMAKUF
! MAXITU IS THE MAXIMUM ITERATION IN BETAF AND MAXITA IS THE MAXIMUM
! ITERATION IN ALPHAF
!
      write (iw, 250) atol, btol, maxitu, maxita 
  250 format('  ATOL=',d12.5,'  BTOL=',d12.5,'    MAXITU=',i5,'    MAXITA=',i5) 
!
! SET UP DIRECT ACCESS FILE FOR T-D MATRICES
      call openda (0) 
!
! CALCULATE ALPHA AT STATIC VALUES
!
      if (iwflb==2 .or. igam==2 .or. igam==4 .or. ibet<=(-2)) then 
        iwfla = 0 
        omega = 0.0D00 
        call alphaf (iwfla, atol, maxita, x1, x2, x3, x4, x5, x6, x7, x8, c, w) 
        if (moperr) return
!
      endif 
      if (igam == 4) then 
        iwflb = 0 
        call betaf (iwflb, maxitu, btol, x1, x2, x3, x4, x5, x6, x7, x8, x19, &
          x10, x11, x12, x13, c, w) 
      endif 
!
! CALCULATE FREQUENCY DEPENDENT VALUES
!
      hartr = ev 
      do i = 1, nfreq 
!
!  READ IN FREQ:  ACTUALLY READ IN AS ENERGY IN EV.
!
        omega = dataev(i) 
        omegau = omega/hartr 
        if (omega < 1.0D-8) then 
!#           WRITE(IW,401) OMEGA
          write (iw, 260) 
  260     format(/,/,' ',65('*'),/,' CALCULATION OF STATIC FIELD QUANTITIES',/,&
            ' ',65('*')) 
        else 
          write (iw, 270) omega, omegau, 1239.8424D0/omega, 8065.541D0*omega 
  270     format(/,/,' ',70('*'),/,' CALCULATION FOR A FREQUENCY OF ',f10.5,&
            ' EV  =',f14.5,' A.U. '/,18x,'WAVELENGTH OF ',f10.2,' NM  =',f14.5,&
            ' CM(-1)',/,' ',70('*')) 
        endif 
!
!  CALCULATE ALPHA(W)
!
        iwfla = 1 
        call alphaf (iwfla, atol, maxita, x1, x2, x3, x4, x5, x6, x7, x8, c, w) 
        if (moperr) return
!
!  PERFORM NONITERATIVE BETA CALCULATIONS
!
!   OPTICAL RECTIFICATION
        if (ibet == (-3)) call nonor (x1, x2, x3, x4, x5, x6, x7, x8, x19, x10&
          , x11, x12) 
!   ELECTROPTIC POCKELS EFFECT
        if (ibet == (-2)) call nonope (x1, x2, x3, x4, x5, x6, x7, x8, x19, x10&
          , x11, x12) 
!   SECOND HARMONIC GENERATION
        if (ibet == (-1)) then 
          iwfla = 2 
          omega = omega*2.0D00 
          call alphaf (iwfla, atol, maxita, x1, x2, x3, x4, x5, x6, x7, x8, c, w) 
          if (moperr) return
!
          omega = omega/2.0D00 
          call nonbet (x1, x2, x3, x4, x5, x6, x7, x8, x19, x10, x11, x12) 
        endif 
!
!  PERFORM ITERATIVE BETA (SHG AND STATIC)CALCULATIONS
!
!         IF ((IBET.GT.0) .AND.(IGAM .EQ. 0)) THEN
        if (ibet==1 .and. iwflb<=1 .and. igam==0) then 
          call betaf (iwflb, maxitu, btol, x1, x2, x3, x4, x5, x6, x7, x8, x19&
            , x10, x11, x12, x13, c, w) 
!
! PERFORM ITERATIVE BETA (EOPE AND OR) CALCULATIONS
!
        else if (ibet==1 .and. iwflb>1 .and. igam==0) then 
          call beopor (iwflb, maxitu, btol, x1, x2, x3, x4, x5, x6, x7, x8, x19&
            , x10, x11, x12, x13, c, w) 
        endif 
!.......................................................................
! CALCULATE GAMMA VALUES
!.......................................................................
        if (ibet>0 .and. igam<=3 .and. igam/=0) then 
          iwflb = 1 
          call betaf (iwflb, maxitu, btol, x1, x2, x3, x4, x5, x6, x7, x8, x19&
            , x10, x11, x12, x13, c, w) 
        endif 
! THIRD HARMONIC GENERATION
        if (igam == 1) then 
          iwfla = 3 
          omega = omega*3.0D00 
          call alphaf (iwfla, atol, maxita, x1, x2, x3, x4, x5, x6, x7, x8, c, w)
          if (moperr) return
!
          omega = omega/3.0D00 
          call ngamtg (x1, x2, x3, x4, x5, x6, x7, x8, x19) 
        endif 
! DC-EFISHG
        if (igam == 2) then 
          iwfla = 2 
          omega = 2.0D00*omega 
          call alphaf (iwfla, atol, maxita, x1, x2, x3, x4, x5, x6, x7, x8, c, w) 
          if (moperr) return
!
          omega = omega/2.0D00 
          iwflb = 2 
          call beopor (iwflb, maxitu, btol, x1, x2, x3, x4, x5, x6, x7, x8, x19&
            , x10, x11, x12, x13, c, w) 
          call ngefis (x1, x2, x3, x4, x5, x6, x7, x8, x19) 
        endif 
! IDRI
        if (igam == 3) then 
          iwflb = 3 
          call beopor (iwflb, maxitu, btol, x1, x2, x3, x4, x5, x6, x7, x8, x19&
            , x10, x11, x12, x13, c, w) 
          call ngidri (x1, x2, x3, x4, x5, x6, x7, x8, x19) 
        endif 
! OKE
        if (igam /= 4) cycle  
        iwflb = 2 
        call beopor (iwflb, maxitu, btol, x1, x2, x3, x4, x5, x6, x7, x8, x19, &
          x10, x11, x12, x13, c, w) 
        call ngoke (igam, x1, x2, x3, x4, x5, x6, x7, x8, x19) 
      end do 
!
      geo(:,:natoms) = refgeo(:,:natoms) 
      na(:natoms) = ref_na(:natoms) 
      return  
      end subroutine polar 



      subroutine tf(ua, ga, ub, gb, t, norbs) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use zerom_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: norbs 
      real(double) , intent(in) :: ua(norbs,norbs) 
      real(double) , intent(in) :: ga(norbs,norbs) 
      real(double) , intent(in) :: ub(norbs,norbs) 
      real(double) , intent(in) :: gb(norbs,norbs) 
      real(double)  :: t(norbs,norbs) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j
      real(double) :: sum1, sum2 
!-----------------------------------------------
!
!  THIS SUBROUTINE CREATES THE NEW T MATRIX
!
!
!  ZERO MATRIX INITIALLY
!
      call zerom (t, norbs) 
!
! CALCULATE T (IJ)(W,W)= SUM(GA(IK)(W)*UB(KJ)(W)+
! GB(IK)(W)*UA(KJ)(W)-UA(IK)(W)GB(KJ)(W)-UB(IK)(W)GA(KJ)(W)
!
      do i = 1, norbs 
        do j = 1, norbs 
          sum1 = 0.0D0 
          sum2 = 0.0D0 
! CALCULATE FOR (W,W), (0,W) VALUES
!
          sum1 = sum(ga(i,:norbs)*ub(:norbs,j)+gb(i,:norbs)*ua(:norbs,j)-ua(i,:&
            norbs)*gb(:norbs,j)-ub(i,:norbs)*ga(:norbs,j)) 
          sum2 = sum(ga(j,:norbs)*ub(:norbs,i)+gb(j,:norbs)*ua(:norbs,i)-ua(j,:&
            norbs)*gb(:norbs,i)-ub(j,:norbs)*ga(:norbs,i)) 
          t(i,j) = sum1 
          t(j,i) = sum2 
        end do 
      end do 
!
      return  
      end subroutine tf 


      subroutine transf(f, g, c, norb) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use molkst_C, only : norbs
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: norb 
      real(double) , intent(in) :: f(norb,norb) 
      real(double) , intent(out) :: g(norb,norb) 
      real(double) , intent(in) :: c(norb,norb) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j
      real(double), dimension(norbs) :: work 
      real(double) :: term, term2 
!-----------------------------------------------
!
!  THIS SUBROUTINE FORMS THE G MATRIX BY TRANSFORMING F WITH C
!
!
      do i = 1, norb 
        do j = 1, norb 
          term = 0.0D00 
          term = sum(f(j,:)*c(:,i)) 
          work(j) = term 
        end do 
        do j = 1, norb 
          term2 = 0.0D00 
          term2 = sum(work(:norb)*c(:,j)) 
          g(j,i) = term2 
        end do 
      end do 
      return  
      end subroutine transf 


      real(kind(0.0d0)) function trsub (ul, x, ur, l1, lm, ndim) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
! THIS PROGRAM CALCULATES TRACES OF MATRICES
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: l1 
      integer , intent(in) :: lm 
      integer , intent(in) :: ndim 
      real(double) , intent(in) :: ul(ndim,ndim) 
      real(double) , intent(in) :: x(ndim,ndim) 
      real(double) , intent(in) :: ur(ndim,ndim) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k, l 
      real(double) :: sum, suml 
!-----------------------------------------------
!
      sum = 0.0D00 
      do i = 1, l1 
        do k = 1, lm 
          suml = 0.0D00 
          do l = 1, lm 
            suml = suml + x(k,l)*ur(l,i) 
          end do 
          sum = sum + suml*ul(i,k) 
        end do 
      end do 
      trsub = 2.0D00*sum 
      return  
      end function trsub 


      real(kind(0.0d0)) function trudgu (ul, x, ur, l1, lm, ndim) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: l1 
      integer , intent(in) :: lm 
      integer , intent(in) :: ndim 
      real(double) , intent(in) :: ul(ndim,ndim) 
      real(double) , intent(in) :: x(ndim,ndim) 
      real(double) , intent(in) :: ur(ndim,ndim) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k, l 
      real(double) :: sum, suml 
!-----------------------------------------------
!
      sum = 0.0D00 
      do i = 1, l1 
        do k = 1, lm 
          suml = 0.0D00 
          do l = 1, lm 
            suml = suml + x(k,l)*ur(l,i) 
          end do 
          sum = sum + suml*ul(k,i) 
        end do 
      end do 
      trudgu = 2.0D00*sum 
      return  
      end function trudgu 


      real(kind(0.0d0)) function trugdu (ul, x, ur, l1, lm, ndim) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: l1 
      integer , intent(in) :: lm 
      integer , intent(in) :: ndim 
      real(double) , intent(in) :: ul(ndim,ndim) 
      real(double) , intent(in) :: x(ndim,ndim) 
      real(double) , intent(in) :: ur(ndim,ndim) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k, l 
      real(double) :: sum, suml 
!-----------------------------------------------
!
      sum = 0.0D00 
      do i = 1, l1 
        do k = 1, lm 
          suml = 0.0D00 
          do l = 1, lm 
            suml = suml + x(l,k)*ur(l,i) 
          end do 
          sum = sum + suml*ul(i,k) 
        end do 
      end do 
      trugdu = 2.0D00*sum 
      return  
      end function trugdu 


      real(kind(0.0d0)) function trugud (ul, x, ur, l1, lm, ndim) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: l1 
      integer , intent(in) :: lm 
      integer , intent(in) :: ndim 
      real(double) , intent(in) :: ul(ndim,ndim) 
      real(double) , intent(in) :: x(ndim,ndim) 
      real(double) , intent(in) :: ur(ndim,ndim) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k, l 
      real(double) :: sum, suml 
!-----------------------------------------------
!
      sum = 0.0D00 
      do i = 1, l1 
        do k = 1, lm 
          suml = 0.0D00 
          do l = 1, lm 
            suml = suml + x(k,l)*ur(i,l) 
          end do 
          sum = sum + suml*ul(i,k) 
        end do 
      end do 
      trugud = 2.0D00*sum 
      return  
      end function trugud 


      real(kind(0.0d0)) function wrdkey (keywrd, key, nk, refkey, nr, def) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use reada_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nk 
      integer , intent(in) :: nr 
      real(double) , intent(in) :: def 
      character  :: keywrd*241 
      character , intent(in) :: key*(*) 
      character , intent(in) :: refkey*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k 
!-----------------------------------------------
!***********************************************************************
!
!   WRDKEY is a way to assign values using keywords.  If the keyword is
!          present, it is used, otherwise the default is used.
!
!   Format:
!          KEY    =  Start of keyword
!          NK     =  Number of characters in (start of keyword)
!          REFKEY =  Variable in keyword
!          NR     =  Number of characters in (variable in keyword)
!
!   For example, if call is WRDKEY(KEYWRD,'VDW(',4,':ZN',3,1.5D0)
!   then if VWD(:ZN=1.4) exists, WRDKEY will be 1.4, otherwise 1.5
!   Note: the "=" and ")" are not necessary, but make the line easier
!   to read.
!***********************************************************************
      i = index(keywrd,key(:nk)) 
      if (i /= 0) then 
        j = index(keywrd(i:),') ') 
        k = index(keywrd(i:i+j),refkey(:nr)) 
        if (k /= 0) then 
          wrdkey = reada(keywrd,i + k) 
        else 
          wrdkey = def 
        endif 
      else 
        wrdkey = def 
      endif 
      return  
      end function wrdkey 


      subroutine zerom(x, m) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: m 
      real(double) , intent(out) :: x(m,m) 
!-----------------------------------------------
!
!  ZEROM ZEROS THE MATRIX X
!
      x = 0.0D00 
      return  
      end subroutine zerom 


      subroutine alphaf(iwfla, atol, maxita, u, f, g, uold, h1, d, da, w2, c, w) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use molkst_C, only : norbs, keywrd, nclose, numat
      use chanel_C, only : iw
      use funcon_C, only: a0, ev
      use common_arrays_C, only : eigs, nat, nfirst, nlast, coord
      USE polar_C, only : omega, alpavg
!
!  SUBROUTINE FOR THE CALCULATION OF THE FREQUENCY DEPENDENT FIRST-ORDER
!  RESPONCE MATRICIES UA AND DENSITIES DA.
!  USED TO COMPUTE THE FREQUENCY DEPENDENT POLARIZABILITY AND FOR
!  SOLVING THE SECOND-ORDER PROBLEM.
!
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use second_I 
      use hmuf_I 
      use copym_I 
      use zerom_I 
      use transf_I 
      use makeuf_I 
      use densf_I 
      use aval_I 
      use ffreq2_I 
      use ffreq1_I 
      use hplusf_I 
      use dawrit_I 
      use pol_vol_I
      use to_screen_I
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: iwfla 
      integer , intent(in) :: maxita 
      real(double)  :: atol 
      real(double)  :: u(norbs,norbs) 
      real(double)  :: f(norbs,norbs) 
      real(double)  :: g(norbs,norbs) 
      real(double)  :: uold(norbs,norbs) 
      real(double)  :: h1(norbs,norbs) 
      real(double)  :: d(norbs,norbs) 
      real(double)  :: da(norbs,norbs) 
      real(double)  :: w2(*) 
      real(double)  :: c(*) 
      real(double)  :: w(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nsqr, iposu, iposg, id, icount, ic, j, iexp
      real(double), dimension(3,3) :: allalp 
      real(double) :: cmptim, alpold, diff, alphaw, dela, &
      const
      logical :: last, debug 
      character, dimension(3) :: alab 

      save alab 
!-----------------------------------------------
!
      data alab/ 'X', 'Y', 'Z'/  
!
      debug = index(keywrd,' DEBUG') /= 0 
      nsqr = norbs*norbs 
      alpavg = 0.0D00 
! COMPUTE OFFSETS FOR U AND G MATRICES
      iposu = 1 + 6*iwfla 
      iposg = 4 + 6*iwfla 
      write (iw, 10) omega 
   10 format(/,' +++++ ALPHA AT ',1f13.5,' EV.') 
!
!  CHOOSE A  COMPONENT
!  X: ID=1   Y: ID=2   Z: ID=3
!
      do id = 1, 3 
        cmptim = second(1) 
        last = .FALSE. 
!
!  CALCULATE THE DIPOLE MATRIX.
!
        call hmuf (h1, id, coord, nfirst, nlast, nat, norbs, numat) 
        call copym (h1, f, norbs) 
        
!
!
!  INITIALIZE UOLD TO ZERO
!
        call zerom (uold, norbs) 
!.................................................................
!  LOOP STARTS HERE
!.................................................................
        icount = 0 
        alpold = 0.0D00 
        icount = icount + 1 
        if (icount > maxita) last = .TRUE. 
!
!  CREATE G MATRIX.
!
        call transf (f, g, c, norbs) 
!
!  FORM U MATRIX
!
        call makeuf (u, uold, g, eigs, last, norbs, nclose, diff, atol) 
!
!  FORM NEW DENSITY MATRIX
!
        call densf (u, c, d, da, norbs, nclose, w2) 
!
! COMPUTE TEST ALPHA TO BE USED FOR A CONVERGENCE TEST
!
        alphaw = aval(h1,d,norbs) 
        dela = dabs(alpold - alphaw) 
        alpold = alphaw 
!
!  CREATE NEW FOCK MATRIX
!
        f = 0.d0
        call ffreq2 (f, d, w) 
        call ffreq1 (f, d, da, da, norbs)  
        f = h1 + f/ev
        do while(.not.last) 
          icount = icount + 1 
          if (icount > maxita) last = .TRUE. 
!
!  CREATE G MATRIX.
!
          call transf (f, g, c, norbs) 
!
!  FORM U MATRIX
!
          call makeuf (u, uold, g, eigs, last, norbs, nclose, diff, atol) 
!
!  FORM NEW DENSITY MATRIX
!
          call densf (u, c, d, da, norbs, nclose, w2) 
!
! COMPUTE TEST ALPHA TO BE USED FOR A CONVERGENCE TEST
!
          alphaw = aval(h1,d,norbs) 
          dela = dabs(alpold - alphaw) 
          alpold = alphaw 
!
!  CREATE NEW FOCK MATRIX
!
          call zerom (f, norbs) 
          call ffreq2 (f, d, w) 
          call ffreq1 (f, d, da, da, norbs) 
          if (Abs(f(1,1)) > 1.d6) then
            call mopend("Calculation of polarizability failed (This part of MOPAC is fragile)")
            return
          end if
          call hplusf (f, h1, norbs) 
!..............................................................
        end do 
        cmptim = second(1) - cmptim 
        if (debug) write (iw, 30) icount, cmptim, diff, dela 
   30   format(/,' CONVERGED IN',i4,' ITERATIONS IN',f10.2,' SECONDS',/,&
          '           DENSITY CONVERG. TO ',1p,d12.5,/,&
          '             ALPHA CONVERG. TO ',1p,d12.5,/) 
!
! COMPUTE ALPHA
!
        alphaw = aval(h1,d,norbs) 
        allalp(id,id) = alphaw 
        if (debug) write (iw, 40) alab(id), alab(id), alphaw 
   40   format('      ALPHA(',a1,',',a1,') = ',1p,d14.7) 
        alpavg = alpavg + alphaw 
!
!  WRITE OUT U AND G FOR FUTURE USE
!
        call dawrit (u, nsqr, iposu + id) 
        call dawrit (g, nsqr, iposg + id) 
!
!  COMPUTE OTHER COMPONENTS
!
        do ic = 1, 3 
          if (ic == id) cycle  
          call hmuf (h1, ic, coord, nfirst, nlast, nat, norbs, numat) 
          alphaw = aval(h1,d,norbs) 
          allalp(ic,id) = alphaw 
          if (debug) write (iw, 50) alab(ic), alab(id), alphaw 
   50     format('      ALPHA(',a1,',',a1,') = ',1p,d14.7) 
        end do 
      end do 
      
      allalp(1,1) = pol_vol(allalp(1,1))
      allalp(2,2) = pol_vol(allalp(2,2))
      allalp(3,3) = pol_vol(allalp(3,3))
      alpavg = (allalp(1,1) + allalp(2,2) + allalp(3,3))/3.d0
      if (alpavg > 1.d6) then
        iexp = Nint (Log10(alpavg)) - 5
        const = 10.d0 ** (-iexp)
        write (iw, "(/27X, 'COMPONENTS OF ALPHA*10**(',I3,')',/)") - iexp
      else
        const = 1.d0
        write (iw, "(/27X, 'COMPONENTS OF ALPHA',/)")
      end if
    write (iw, "(21X,'*X',11X,'*Y',11X,'*Z')")
    write (iw, "(10X,' *=X:',3F13.4)") (allalp(1, j)*const, j=1, 3)
    write (iw, "(10X,' *=Y:',3F13.4)") (allalp(2, j)*const, j=1, 3)
    write (iw, "(10X,' *=Z:',3F13.4)") (allalp(3, j)*const, j=1, 3)
     if (const < 0.9d0) then
      write (iw, 10035) - iexp, alpavg * const,alpavg * const*a0**3
10035 format (/, "  ISOTROPIC AVERAGE ALPHA*10**(", I3, ") = ", f13.5, " A.U. =",&
  & f13.5, " ANG.**3")
    else
      write (iw, 10040) alpavg, alpavg*a0**3
10040 format (/, "  ISOTROPIC AVERAGE ALPHA = ", 1 f13.5, " A.U. =", &
           & f13.5, " ANG.**3")
    end if
    call to_screen("To_file: POLAR")
!
      return  
      end subroutine alphaf 


      real(kind(0.0d0)) function aval (h, d, norbs) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: norbs 
      real(double) , intent(in) :: h(norbs,norbs) 
      real(double) , intent(in) :: d(norbs,norbs) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j 
      real(double) :: sum 
!-----------------------------------------------
!.................................................................
!  COMPUTE POLARIZABILITY AS TRACE OF H*D
!.................................................................
      sum = 0.0D00 
      do i = 1, norbs 
        do j = 1, norbs 
          sum = sum + h(i,j)*d(j,i) 
        end do 
      end do 
      aval = -sum 
      return  
      end function aval 


      subroutine bdenin(bdcon, ua, ub, c, norbs, nclose) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: norbs 
      integer , intent(in) :: nclose 
      real(double) , intent(out) :: bdcon(norbs,norbs) 
      real(double) , intent(in) :: ua(norbs,norbs) 
      real(double) , intent(in) :: ub(norbs,norbs) 
      real(double) , intent(in) :: c(norbs,norbs) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, l, k 
      real(double), dimension(norbs) :: w1, w2, w3, w4 
      real(double) :: sum1, sum2, sum3, sum4, s3, s4 
!-----------------------------------------------
!
!  THIS SUBROUTINE IS USED TO COMPUTE THE FIRST-ORDER DENSITY
!
!
! CALCULATE
!
      do i = 1, norbs 
        do j = 1, norbs 
!
          do l = nclose + 1, norbs 
            sum1 = 0.D0 
            sum2 = 0.D0 
            sum1 = sum(ub(l,:nclose)*c(j,:nclose)) 
            sum2 = sum(ua(l,:nclose)*c(j,:nclose)) 
            w1(l) = sum1 
            w2(l) = sum2 
          end do 
          do k = 1, nclose 
            sum3 = 0.D0 
            sum3 = sum((ua(k,nclose+1:norbs)*w1(nclose+1:norbs)+ub(k,nclose+1:&
              norbs)*w2(nclose+1:norbs))) 
            w3(k) = sum3 
          end do 
!
          do l = 1, nclose 
            sum1 = 0.D0 
            sum2 = 0.D0 
            sum1 = sum(ub(l,nclose+1:norbs)*c(j,nclose+1:norbs)) 
            sum2 = sum(ua(l,nclose+1:norbs)*c(j,nclose+1:norbs)) 
            w1(l) = sum1 
            w2(l) = sum2 
          end do 
          do k = nclose + 1, norbs 
            sum4 = 0.D0 
            sum4 = sum((ua(k,:nclose)*w1(:nclose)+ub(k,:nclose)*w2(:nclose))) 
            w4(k) = sum4 
          end do 
          s3 = 0.D0 
          s3 = sum(w3(:nclose)*c(i,:nclose)) 
          s4 = 0.D0 
          s4 = sum(w4(nclose+1:norbs)*c(i,nclose+1:norbs)) 
          bdcon(i,j) = s3 - s4 
        end do 
      end do 
!
      return  
      end subroutine bdenin 


      subroutine bdenup(bdcon, uab, c, d, da, norbs, nclose) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use zerom_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: norbs 
      integer , intent(in) :: nclose 
      real(double) , intent(in) :: bdcon(norbs,norbs) 
      real(double) , intent(in) :: uab(norbs,norbs) 
      real(double) , intent(in) :: c(norbs,norbs) 
      real(double)  :: d(norbs,norbs) 
      real(double) , intent(inout) :: da(norbs,norbs) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, k, l, i 
      real(double), dimension(norbs) :: w1 
      real(double) :: sum, s1, s2 
!-----------------------------------------------
!
!  THIS SUBROUTINE IS USED TO COMPUTE THE FIRST-ORDER DENSITY
!
!
! FORM DENSITY MATRIX
!
!
      call zerom (d, norbs) 
!
! CALCULATE
!
      do j = 1, norbs 
        do k = 1, norbs 
          sum = 0.D0 
          do l = 1, nclose 
            sum = sum + uab(k,l)*c(j,l) 
          end do 
          da(k,j) = sum 
        end do 
      end do 
      do i = 1, norbs 
        do k = 1, norbs 
          sum = 0.D0 
          do l = 1, nclose 
            sum = sum + c(i,l)*uab(l,k) 
          end do 
          w1(k) = sum 
        end do 
        do j = 1, norbs 
          s1 = 0.0D00 
          s2 = 0.0D00 
          do k = 1, norbs 
            s1 = s1 + c(i,k)*da(k,j) 
            s2 = s2 + w1(k)*c(j,k) 
          end do 
          d(i,j) = 2.0D00*(s1 - s2 + bdcon(i,j)) 
        end do 
      end do 
      da(:norbs,:norbs) = d(:norbs,:norbs)*0.5D0 
!
      return  
      end subroutine bdenup 


      subroutine beopor(iwflb, maxitu, btol, ua, ub, f, ga, gb, t, h1, d, da, &
        uab, uold1, g, x, c, w) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use molkst_C, only : norbs, numat, nclose, line
      use common_arrays_C, only : eigs, nat, nfirst, nlast, coord
      use chanel_C, only : iw
      USE polar_C, only : omega 
!
! THIS SUBROUTINE CALCULATES ITERATIVE BETA VALUES FOR
! THE ELECTROOPTIC POCKELS EFFECT AND OPTICAL RECTIFICATION
!
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use second_I 
      use hmuf_I 
      use zerom_I 
      use daread_I 
      use fhpatn_I 
      use tf_I 
      use bdenin_I 
      use bdenup_I 
      use aval_I 
      use ffreq2_I 
      use ffreq1_I 
      use hplusf_I 
      use transf_I 
      use bmakuf_I 
      use epsab_I 
      use dawrit_I 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: iwflb 
      integer , intent(in) :: maxitu 
      real(double)  :: btol 
      real(double)  :: ua(norbs,norbs) 
      real(double)  :: ub(norbs,norbs) 
      real(double)  :: f(norbs,norbs) 
      real(double)  :: ga(norbs,norbs) 
      real(double)  :: gb(norbs,norbs) 
      real(double)  :: t(norbs,norbs) 
      real(double)  :: h1(norbs,norbs) 
      real(double)  :: d(norbs,norbs) 
      real(double)  :: da(norbs,norbs) 
      real(double)  :: uab(norbs,norbs) 
      real(double)  :: uold1(norbs,norbs) 
      real(double)  :: g(norbs,norbs) 
      real(double)  :: x(norbs,norbs) 
      real(double)  :: c(*) 
      real(double)  :: w(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(9) :: ida, idb 
      integer :: maxsq, iposu, iposg, ipose, iposum, id, ia, ib, jpu, jpg, &
        icount, ic 
      real(double) :: maxu, one, bavx, bavy, bavz, cmptim, betaw, diff, bvec 
      logical :: last 
      character, dimension(3) :: alab 

      save alab, ida, idb 
!-----------------------------------------------
      data alab/ 'X', 'Y', 'Z'/  
      data ida/ 1, 1, 1, 2, 2, 2, 3, 3, 3/  
      data idb/ 1, 2, 3, 1, 2, 3, 1, 2, 3/  
      one = 1.0D00 
      maxsq = norbs*norbs 
      if (iwflb == 2) then 
        iposu = 73 
      else 
        iposu = 109 
      endif 
      iposg = iposu + 9 
      ipose = iposg + 9 
      iposum = ipose + 9 
      if (iwflb == 0) then 
        write (iw, 10) omega 
   10   format(/,' +++++ BETA (STATIC) AT ',1f15.5,' EV.'/) 
      else if (iwflb == 2) then 
        write (iw, 20) omega 
   20   format(/,' +++++ BETA',' (ELECTROOPTIC POCKELS EFFECT) AT ',1f15.5,&
          ' EV.'/) 
      else 
        write (iw, 30) omega 
   30   format(/,' +++++ BETA',' (OPTICAL RECTIFICATION) AT ',1f15.5,' EV.'/) 
      endif 
!
!  LOOP OVER COMPONENTS
!
      bavx = 0.0D+00 
      bavy = 0.0D+00 
      bavz = 0.0D+00 
      do id = 1, 9 
        cmptim = second(1) 
        ia = ida(id) 
        ib = idb(id) 
        last = .FALSE. 
!
!  CALCULATE THE DIPOLE MATRIX.
!
        call hmuf (h1, ia, coord, nfirst, nlast, nat, norbs, numat) 
!
!  INITIALIZE ZERO ARRAYS
!
        call zerom (uold1, norbs) 
        call zerom (uab, norbs) 
        call zerom (f, norbs) 
!
!  INPUT U AND GA FROM ALPHA CALCULATIONS
!
        if (iwflb==2 .or. iwflb==0) then 
!  UA CONTAINS UA(0)
          jpu = 1 + ia 
          call daread (ua, maxsq, jpu) 
!  GA CONTAINS GA(0)
          jpg = 4 + ia 
          call daread (ga, maxsq, jpg) 
        else 
!  UA CONTAINS UA(W)
          jpu = 7 + ia 
          call daread (ua, maxsq, jpu) 
!  GA CONTAINS GA(W)
          jpg = 10 + ia 
          call daread (ga, maxsq, jpg) 
        endif 
!
! READ VALUES FOR (W,-W) CALCULATION  :  OR
!
        if (iwflb == 3) then 
!  UB CONTAINS UB(-W) = -UB+(W)
          jpu = 7 + ib 
          call daread (x, maxsq, jpu) 
          call fhpatn (ub, x, norbs, 2, (-one)) 
!  GB CONTAINS GB(-W) = GB+(W)
          jpg = 10 + ib 
          call daread (x, maxsq, jpg) 
          call fhpatn (gb, x, norbs, 2, one) 
!
! READ VALUES FOR (0,W) CALCULATION  :  OKE
!
        else if (iwflb == 0) then 
!  UB CONTAINS UB(0)
          jpu = 1 + ib 
          call daread (ub, maxsq, jpu) 
!  GB CONTAINS GB(0)
          jpg = 4 + ib 
          call daread (gb, maxsq, jpg) 
        else 
!  UB CONTAINS UB(W)
          jpu = 7 + ib 
          call daread (ub, maxsq, jpu) 
!  GB CONTAINS GB(W)
          jpg = 10 + ib 
          call daread (gb, maxsq, jpg) 
        endif 
!
!  CONSTRUCT T-MATRIX ONE TIME
!
        call tf (ua, ga, ub, gb, t, norbs) 
!
!  CALCULATE INITIAL DENSITY AND BETA VALUE
!
        call bdenin (x, ua, ub, c, norbs, nclose) 
        call bdenup (x, uab, c, d, da, norbs, nclose) 
        betaw = aval(h1,d,norbs) 
!
! INITIALIZE FOCK MATRIX
!
        call ffreq2 (f, d, w) 
        call ffreq1 (f, d, da, da, norbs) 
        call zerom (da, norbs) 
        call hplusf (f, da, norbs) 
!.................................................................
!  LOOP STARTS HERE
!.................................................................
        icount = 0 
        icount = icount + 1 
        if (icount >= maxitu) last = .TRUE. 
!
!  CREATE G MATRIX.
!
        call transf (f, g, c, norbs) 
!
!  FORM U MATRIX
!
        call bmakuf (ua, ub, uab, t, uold1, g, eigs, last, norbs, nclose, diff&
          , iwflb, maxu, btol) 
!
!  FORM NEW DENSITY MATRIX
!
        call bdenin (x, ua, ub, c, norbs, nclose) 
        call bdenup (x, uab, c, d, da, norbs, nclose) 
!...
! COMPUTE TEST BETA
!
        betaw = aval(h1,d,norbs) 
!         IF (LAST.OR.(ICOUNT.GT.(MAXITU-5))) THEN
!             WRITE(IW,1500) ICOUNT,DELA,MAXU,DIFF
! 1500        FORMAT(' ',I4,'  DELTA BETA = ', D12.5,
!     X       ' MAXU = ', D12.5, '  UDIFF = ', D12.5)
!         ENDIF
!
!  CREATE NEW FOCK MATRIX
!
        call zerom (f, norbs) 
        call ffreq2 (f, d, w) 
        call ffreq1 (f, d, da, da, norbs) 
        call zerom (da, norbs) 
        call hplusf (f, da, norbs) 
        do while(.not.last) 
          icount = icount + 1 
          if (icount >= maxitu) last = .TRUE. 
!
!  CREATE G MATRIX.
!
          call transf (f, g, c, norbs) 
!
!  FORM U MATRIX
!
          call bmakuf (ua, ub, uab, t, uold1, g, eigs, last, norbs, nclose, &
            diff, iwflb, maxu, btol) 
!
!  FORM NEW DENSITY MATRIX
!
          call bdenin (x, ua, ub, c, norbs, nclose) 
          call bdenup (x, uab, c, d, da, norbs, nclose) 
!...
! COMPUTE TEST BETA
!
          betaw = aval(h1,d,norbs) 
!         IF (LAST.OR.(ICOUNT.GT.(MAXITU-5))) THEN
!             WRITE(IW,1500) ICOUNT,DELA,MAXU,DIFF
! 1500        FORMAT(' ',I4,'  DELTA BETA = ', D12.5,
!     X       ' MAXU = ', D12.5, '  UDIFF = ', D12.5)
!         ENDIF
!
!  CREATE NEW FOCK MATRIX
!
          call zerom (f, norbs) 
          call ffreq2 (f, d, w) 
          call ffreq1 (f, d, da, da, norbs) 
          call zerom (da, norbs) 
          call hplusf (f, da, norbs) 
!..............................................................
        end do 
        cmptim = second(1) - cmptim 
        write (iw, 50) icount, cmptim 
   50   format(/,' CONVERGED IN',i4,' ITERATIONS IN',f10.2,' SECONDS') 
        write (iw, 60) maxu, diff 
   60   format(' MAXIMUM UAB ELEMENT =',1f15.5,',  MAXIMUM DIFFERENCE =',1f15.5&
          ,/) 
!
!  COMPUTE OTHER COMPONENTS
!
        do ic = 1, 3 
          call hmuf (h1, ic, coord, nfirst, nlast, nat, norbs, numat) 
          betaw = aval(h1,d,norbs) 
          write (iw, 70) alab(ic), alab(ia), alab(ib), betaw 
   70     format('      BETA(',a1,',',a1,',',a1,') = ',1f15.5) 
! CALCULATES THE AVERAGE VALUE OF BETA
!
          if (id==1 .and. ic==1) then 
            bavx = bavx + 3.0D0*betaw 
          else if ((id==5 .or. id==9) .and. ic==1) then 
            bavx = bavx + betaw 
          else if ((id==2 .or. id==4) .and. ic==2) then 
            bavx = bavx + betaw 
          else if ((id==3 .or. id==7) .and. ic==3) then 
            bavx = bavx + betaw 
          endif 
! CALCULATES AVERAGE BETA IN Y-DIRECTION
!
          if (id==5 .and. ic==2) then 
            bavy = bavy + 3.0D0*betaw 
          else if ((id==2 .or. id==4) .and. ic==1) then 
            bavy = bavy + betaw 
          else if ((id==1 .or. id==9) .and. ic==2) then 
            bavy = bavy + betaw 
          else if ((id==6 .or. id==8) .and. ic==3) then 
            bavy = bavy + betaw 
          endif 
! CALCULATES AVERAGE BETA IN THE Z-DIRECTION
!
          if (id==9 .and. ic==3) then 
            bavz = bavz + 3.0D0*betaw 
          else if ((id==3 .or. id==7) .and. ic==1) then 
            bavz = bavz + betaw 
          else if ((id==6 .or. id==8) .and. ic==2) then 
            bavz = bavz + betaw 
          else if ((id==1 .or. id==5) .and. ic==3) then 
            bavz = bavz + betaw 
          endif 
        end do 
!
! CALL SUBROUTINE TO CALCULATE EPSILON AND UMINUS OMEGA,OMEGA
!  EPSILON IN H1 AND UMINUS IN DA
        call epsab (h1, eigs, g, ga, gb, ua, ub, uab, da, norbs, nclose, iwflb) 
        call dawrit (uab, maxsq, iposu + id) 
        call dawrit (g, maxsq, iposg + id) 
        call dawrit (h1, maxsq, ipose + id) 
        call dawrit (da, maxsq, iposum + id) 
      end do 
      bavx = bavx/5.0D+00 
      bavy = bavy/5.0D+00 
      bavz = bavz/5.0D+00 
      bvec = (bavx*bavx + bavy*bavy + bavz*bavz)**0.5D+00 
!
      if (bvec < 1.d10) then
        line = "(' AVERAGE BETA',a,' VALUE AT',f10.5,' EV = ',f11.4)"
      else if (bvec < 1.d15) then
        line = "(' AVERAGE BETA',a,' VALUE AT',f10.5,' EV = ',f16.4)"
      else
        line = "(' AVERAGE BETA',a,' VALUE AT',f10.5,' EV = ',f21.4)"
      end if
      write (iw,"(/)")
      write (iw, trim(line)) "X", omega, bavx 
      write (iw, trim(line)) "Y", omega, bavy 
      write (iw, trim(line)) "Z", omega, bavz 
      write (iw,"(2/)")
      write (iw, trim(line)) " ", omega, bvec 
      write (iw,"(/)")
      return  
      end subroutine beopor 


      subroutine betaf(iwflb, maxitu, btol, ua, ub, f, ga, gb, t, h1, d, da, &
        uab, uold1, g, x, c, w) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use common_arrays_C, only : eigs, nat, nfirst, nlast, coord
      use chanel_C, only : iw
      use molkst_C, only : norbs, nclose, numat, keywrd, line
      USE polar_C, only : omega 
!
! THIS SUBROUTINE CALCULATES ITERATIVE BETA VALUES FOR SECOND HARMONIC
! GENERATION.
!
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use second_I 
      use hmuf_I 
      use zerom_I 
      use daread_I 
      use fhpatn_I 
      use tf_I 
      use bdenin_I 
      use bdenup_I 
      use aval_I 
      use ffreq2_I 
      use ffreq1_I 
      use hplusf_I 
      use transf_I 
      use bmakuf_I 
      use epsab_I 
      use dawrit_I 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: iwflb 
      integer , intent(in) :: maxitu 
      real(double)  :: btol 
      real(double)  :: ua(norbs,norbs) 
      real(double)  :: ub(norbs,norbs) 
      real(double)  :: f(norbs,norbs) 
      real(double)  :: ga(norbs,norbs) 
      real(double)  :: gb(norbs,norbs) 
      real(double)  :: t(norbs,norbs) 
      real(double)  :: h1(norbs,norbs) 
      real(double)  :: d(norbs,norbs) 
      real(double)  :: da(norbs,norbs) 
      real(double)  :: uab(norbs,norbs) 
      real(double)  :: uold1(norbs,norbs) 
      real(double)  :: g(norbs,norbs) 
      real(double)  :: x(norbs,norbs) 
      real(double)  :: c(*) 
      real(double)  :: w(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(6) :: ida, idb 
      integer :: maxsq, iposu, iposg, ipose, iposum, id, ia, ib, jpu, jpg, &
        icount, ic, i, j 
      real(double) :: maxu 
      real(double), dimension(3,3,3) :: allbet 
      real(double) :: one, bavx, bavy, bavz, cmptim, betaw, diff, bvec, &
        au_to_esu = 8.639418d-33
      logical :: last, debug 
      character, dimension(3) :: alab 

      save alab, ida, idb 
!-----------------------------------------------
!
!
      data alab/ 'X', 'Y', 'Z'/  
      data ida/ 1, 1, 1, 2, 2, 3/  
      data idb/ 1, 2, 3, 2, 3, 3/  
!
      debug = index(keywrd,' BETAF ') /= 0 
      one = 1.0D00 
      maxsq = norbs*norbs 
      iposu = 25 + 24*iwflb 
      iposg = iposu + 6 
      ipose = iposg + 6 
      iposum = ipose + 6 
!
      if (iwflb == 0) then 
        write (iw, 10) omega 
   10   format(/,' +++++ BETA (STATIC) AT ',1f15.5,' EV.'/) 
      else 
        write (iw, 20) omega 
   20   format(/,' +++++ BETA',' (SECOND HARMONIC GENERATION) AT ',1f13.5,&
          ' EV.'/) 
      endif 
!
!  CHOOSE A  COMPONENT
!  X: ID=1   Y: ID=2   Z: ID=3
!
      bavx = 0.0D+00 
      bavy = 0.0D+00 
      bavz = 0.0D+00 
      do id = 1, 6 
        cmptim = second(1) 
        ia = ida(id) 
        ib = idb(id) 
        last = .FALSE. 
!
!  CALCULATE THE DIPOLE MATRIX.
!
        call hmuf (h1, ia, coord, nfirst, nlast, nat, norbs, numat) 
!
!  INITIALIZE ZERO ARRAYS
!
        call zerom (uold1, norbs) 
        call zerom (uab, norbs) 
        call zerom (f, norbs) 
!
!  INPUT U AND GA FROM ALPHA CALCULATIONS
!
        if (iwflb==2 .or. iwflb==0) then 
          jpu = 1 + ia 
          call daread (ua, maxsq, jpu) 
          jpg = 4 + ia 
          call daread (ga, maxsq, jpg) 
        else 
          jpu = 7 + ia 
          call daread (ua, maxsq, jpu) 
          jpg = 10 + ia 
          call daread (ga, maxsq, jpg) 
        endif 
! READ VALUES FOR (W,-W)
        if (iwflb == 3) then 
          jpu = 7 + ib 
          call daread (x, maxsq, jpu) 
          call fhpatn (ub, x, norbs, 2, (-one)) 
          jpg = 10 + ib 
          call daread (x, maxsq, jpg) 
          call fhpatn (gb, x, norbs, 2, one) 
! READ VALUES FOR OKE
!
        else if (iwflb == 0) then 
          jpu = 1 + ib 
          call daread (ub, maxsq, jpu) 
          jpg = 4 + ib 
          call daread (gb, maxsq, jpg) 
        else 
          jpu = 7 + ib 
          call daread (ub, maxsq, jpu) 
          jpg = 10 + ib 
          call daread (gb, maxsq, jpg) 
        endif 
!
!  CONSTRUCT T-MATRIX ONE TIME
!
        call tf (ua, ga, ub, gb, t, norbs) 
!
!  CALCULATE INITIAL DENSITY AND BETA VALUE
!
        call bdenin (x, ua, ub, c, norbs, nclose) 
        call bdenup (x, uab, c, d, da, norbs, nclose) 
        betaw = aval(h1,d,norbs) 
!
! INITIALIZE FOCK MATRIX
!
        call ffreq2 (f, d, w) 
        call ffreq1 (f, d, da, da, norbs) 
        call zerom (da, norbs) 
        call hplusf (f, da, norbs) 
!.................................................................
!  LOOP STARTS HERE
!.................................................................
        icount = 0 
        icount = icount + 1 
        if (icount >= maxitu) last = .TRUE. 
!
!  CREATE G MATRIX.
!
        call transf (f, g, c, norbs) 
!
!  FORM U MATRIX
!
        call bmakuf (ua, ub, uab, t, uold1, g, eigs, last, norbs, nclose, diff&
          , iwflb, maxu, btol) 
!
!  FORM NEW DENSITY MATRIX
!
        call bdenup (x, uab, c, d, da, norbs, nclose) 
!...
! COMPUTE TEST BETA
!
        betaw = aval(h1,d,norbs) 
!         IF (LAST.OR.(ICOUNT.GT.(MAXITU-5))) THEN
!             WRITE(IW,1500) ICOUNT,DELA,MAXU,DIFF
! 1500        FORMAT(' ',I4,'  DELTA BETA = ', D12.5,
!     X       ' MAXU = ', D12.5, '  UDIFF = ', D12.5)
!         ENDIF
!
!  CREATE NEW FOCK MATRIX
!
        call zerom (f, norbs) 
        call ffreq2 (f, d, w) 
        call ffreq1 (f, d, da, da, norbs) 
        call zerom (da, norbs) 
        call hplusf (f, da, norbs) 
        do while(.not.last) 
          icount = icount + 1 
          if (icount >= maxitu) last = .TRUE. 
!
!  CREATE G MATRIX.
!
          call transf (f, g, c, norbs) 
!
!  FORM U MATRIX
!
          call bmakuf (ua, ub, uab, t, uold1, g, eigs, last, norbs, nclose, &
            diff, iwflb, maxu, btol) 
!
!  FORM NEW DENSITY MATRIX
!
          call bdenup (x, uab, c, d, da, norbs, nclose) 
!...
! COMPUTE TEST BETA
!
          betaw = aval(h1,d,norbs) 
!         IF (LAST.OR.(ICOUNT.GT.(MAXITU-5))) THEN
!             WRITE(IW,1500) ICOUNT,DELA,MAXU,DIFF
! 1500        FORMAT(' ',I4,'  DELTA BETA = ', D12.5,
!     X       ' MAXU = ', D12.5, '  UDIFF = ', D12.5)
!         ENDIF
!
!  CREATE NEW FOCK MATRIX
!
          call zerom (f, norbs) 
          call ffreq2 (f, d, w) 
          call ffreq1 (f, d, da, da, norbs) 
          call zerom (da, norbs) 
          call hplusf (f, da, norbs) 
!..............................................................
        end do 
        cmptim = second(1) - cmptim 
        if (debug) write (iw, 40) icount, cmptim 
   40   format(/,' CONVERGED IN',i4,' ITERATIONS IN',f10.2,' SECONDS') 
        if (debug) write (iw, 50) maxu, diff 
   50   format(' MAXIMUM UAB ELEMENT =',1f15.5,',  MAXIMUM DIFFERENCE =',1f15.5&
          ,/) 
!
! COMPUTE BETA
!
!        CALL HMUF(H1,ID,COORD,NFIRST,NLAST,NAT,NORBS,NUMAT)
!        BETAW = AVAL(H1,D,NORBS)
!        WRITE(IW,2000) ALAB(ID),ALAB(ID),ALAB(ID),BETAW
!2000    FORMAT('BETA(',A1,',',A1,','A1,') = ',D12.5)
!
!  COMPUTE OTHER COMPONENTS
!
        do ic = 1, 3 
          call hmuf (h1, ic, coord, nfirst, nlast, nat, norbs, numat) 
          betaw = aval(h1,d,norbs) 
          allbet(ic,ia,ib) = betaw 
          if (debug) write (iw, 60) alab(ic), alab(ia), alab(ib), betaw 
   60     format('      BETA(',a1,',',a1,',',a1,') = ',1f15.5) 
!
! CALCULATE AVERAGE BETA IN THE X-DIRECTION
!
          if (id==1 .and. ic==1) then 
            bavx = bavx + 3.0D0*betaw 
          else if (id==2 .and. ic==2) then 
            bavx = bavx + 2.0D0*betaw 
          else if (id==3 .and. ic==3) then 
            bavx = bavx + 2.0D0*betaw 
          else if ((id==4 .or. id==6) .and. ic==1) then 
            bavx = bavx + betaw 
          endif 
! CALCULATES AVERAGE BETA IN THE Y-DIRECTION
          if (id==4 .and. ic==2) then 
            bavy = bavy + 3.0D0*betaw 
          else if (id==2 .and. ic==1) then 
            bavy = bavy + 2.0D0*betaw 
          else if (id==5 .and. ic==3) then 
            bavy = bavy + 2.0D0*betaw 
          else if ((id==1 .or. id==6) .and. ic==2) then 
            bavy = bavy + betaw 
          endif 
! CALCULATES AVERAGE BETA IN THE Z-DIRECTION
          if (id==6 .and. ic==3) then 
            bavz = bavz + 3.0D0*betaw 
          else if (id==3 .and. ic==1) then 
            bavz = bavz + 2.0D0*betaw 
          else if (id==5 .and. ic==2) then 
            bavz = bavz + 2.0D0*betaw 
          else if ((id==4 .or. id==1) .and. ic==3) then 
            bavz = bavz + betaw 
          endif 
        end do 
!
!
! CALL SUBROUTINE TO CALCULATE EPSILON AND UMINUS OMEGA,OMEGA
!  EPSILON IN H1 AND UMINUS IN DA
        call epsab (h1, eigs, g, ga, gb, ua, ub, uab, da, norbs, nclose, iwflb) 
        call dawrit (uab, maxsq, iposu + id) 
        call dawrit (g, maxsq, iposg + id) 
        call dawrit (h1, maxsq, ipose + id) 
        call dawrit (da, maxsq, iposum + id) 
      end do 
      write (iw, '(30X,''COMPONENTS OF BETA'')')        
      bavx = bavx/5.0D+00 
      bavy = bavy/5.0D+00 
      bavz = bavz/5.0D+00 
! CALCULATES AVERAGE BETA
      bvec = (bavx*bavx + bavy*bavy + bavz*bavz)**0.5D+00 
      if (bvec < 0.999d5) then
        line = "(3X,a,6F12.5)"
      else 
        line = "(3X,a,6g12.4)"
      end if
      write (iw, &
      '(12X,''*XX'', 9X,''*XY'', 9X,''*YY'', 9X,               ''*XZ'', 9X,''*Y&
      &Z'', 9X,''*ZZ'')') 
      write (iw, trim(line)) "*=X", ((allbet(1,i,j),i=1,j),j=1,3) 
      write (iw, trim(line)) "*=Y", ((allbet(2,i,j),i=1,j),j=1,3) 
      write (iw, trim(line)) "*=Z", ((allbet(3,i,j),i=1,j),j=1,3) 
!
!
      if (bvec < 1.d10) then
        line = "(' AVERAGE BETA',a,'(SHG) VALUE AT',f10.5,' EV = ',f15.4,a,g14.6,a)"
      else 
        line = "(' AVERAGE BETA',a,'(SHG) VALUE AT',f10.5,' EV = ',g13.6,a,g14.6,a)"
      end if
      write (iw,"(/)")
      write (iw, trim(line)) "X", omega, bavx, " a.u. =", bavx*au_to_esu, " ESU"
      write (iw, trim(line)) "Y", omega, bavy, " a.u. =", bavy*au_to_esu, " ESU"
      write (iw, trim(line)) "Z", omega, bavz, " a.u. =", bavz*au_to_esu, " ESU"
      write (iw,"(2/)")
      write (iw, trim(line)) " ", omega, bvec, " a.u. =", bvec*au_to_esu, " ESU"
      write (iw, "(a)") "   (1 a.u. = 8.639418X10-33 esu)"
      write (iw,"(/)")
      return  
      end subroutine betaf 


      subroutine betal1(u0a, g0a, u1b, g1b, u1c, g1c, nclose, norbs, term) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!
! THIS SUBROUTINE CALCULATES THE TRACE OF UGU MATRICES
!
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use trugud_I 
      use trudgu_I 
      use trugdu_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: nclose 
      integer  :: norbs 
      real(double) , intent(out) :: term 
      real(double)  :: u0a(norbs,norbs) 
      real(double)  :: g0a(norbs,norbs) 
      real(double)  :: u1b(norbs,norbs) 
      real(double)  :: g1b(norbs,norbs) 
      real(double)  :: u1c(norbs,norbs) 
      real(double)  :: g1c(norbs,norbs) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double) :: t1a, t2a, t3a, t4a, t5a, t6a, t1b, t2b, t3b, t4b, t5b, &
        t6b 
!-----------------------------------------------
      t1a = trugud(u0a,g1b,u1c,nclose,norbs,norbs) 
      t2a = trudgu(u1c,g1b,u0a,nclose,norbs,norbs) 
      t3a = trugdu(u1b,g1c,u0a,nclose,norbs,norbs) 
      t4a = trugdu(u0a,g1c,u1b,nclose,norbs,norbs) 
      t5a = trudgu(u1c,g0a,u1b,nclose,norbs,norbs) 
      t6a = trugud(u1b,g0a,u1c,nclose,norbs,norbs) 
      t1b = trugud(u0a,g1b,u1c,norbs,nclose,norbs) 
      t2b = trudgu(u1c,g1b,u0a,norbs,nclose,norbs) 
      t3b = trugdu(u1b,g1c,u0a,norbs,nclose,norbs) 
      t4b = trugdu(u0a,g1c,u1b,norbs,nclose,norbs) 
      t5b = trudgu(u1c,g0a,u1b,norbs,nclose,norbs) 
      t6b = trugud(u1b,g0a,u1c,norbs,nclose,norbs) 
      term = t1b - t1a + t2b - t2a + t3a - t3b + t4a - t4b + t5b - t5a + t6b - &
        t6a 
      return  
      end subroutine betal1 


      subroutine betall(u2a, g2a, u1b, g1b, u1c, g1c, nclose, norbs, term) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!
! THIS SUBROUTINE CALCULATES TRACE OF UGU MATRICES
! WHEN A,B,C DIRECTIONS ARE DIFFERENT
!
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use trudgu_I 
      use trugud_I 
      use trugdu_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: nclose 
      integer  :: norbs 
      real(double) , intent(out) :: term 
      real(double)  :: u2a(norbs,norbs) 
      real(double)  :: g2a(norbs,norbs) 
      real(double)  :: u1b(norbs,norbs) 
      real(double)  :: g1b(norbs,norbs) 
      real(double)  :: u1c(norbs,norbs) 
      real(double)  :: g1c(norbs,norbs) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double) :: t1a, t2a, t3a, t4a, t5a, t6a, t1b, t2b, t3b, t4b, t5b, &
        t6b 
!-----------------------------------------------
      t1a = trudgu(u2a,g1b,u1c,nclose,norbs,norbs) 
      t2a = trugud(u1c,g1b,u2a,nclose,norbs,norbs) 
      t3a = trugud(u1b,g1c,u2a,nclose,norbs,norbs) 
      t4a = trudgu(u2a,g1c,u1b,nclose,norbs,norbs) 
      t5a = trugdu(u1c,g2a,u1b,nclose,norbs,norbs) 
      t6a = trugdu(u1b,g2a,u1c,nclose,norbs,norbs) 
      t1b = trudgu(u2a,g1b,u1c,norbs,nclose,norbs) 
      t2b = trugud(u1c,g1b,u2a,norbs,nclose,norbs) 
      t3b = trugud(u1b,g1c,u2a,norbs,nclose,norbs) 
      t4b = trudgu(u2a,g1c,u1b,norbs,nclose,norbs) 
      t5b = trugdu(u1c,g2a,u1b,norbs,nclose,norbs) 
      t6b = trugdu(u1b,g2a,u1c,norbs,nclose,norbs) 
      term = t1b - t1a + t2b - t2a + t3b - t3a + t4b - t4a + t5a - t5b + t6a - &
        t6b 
      return  
      end subroutine betall 


      subroutine betcom(u1, g1, u2, g2, nclose, norbs, term) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!
! THIS SUBROUTINE CALCULATES TRACE OF UGU MATRICES
!
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use trudgu_I 
      use trugud_I 
      use trugdu_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: nclose 
      integer  :: norbs 
      real(double) , intent(out) :: term 
      real(double)  :: u1(norbs,norbs) 
      real(double)  :: g1(norbs,norbs) 
      real(double)  :: u2(norbs,norbs) 
      real(double)  :: g2(norbs,norbs) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double) :: t1a, t2a, t3a, t1b, t2b, t3b 
!-----------------------------------------------
      t1a = trudgu(u2,g1,u1,nclose,norbs,norbs) 
      t2a = trugud(u1,g1,u2,nclose,norbs,norbs) 
      t3a = trugdu(u1,g2,u1,nclose,norbs,norbs) 
      t1b = trudgu(u2,g1,u1,norbs,nclose,norbs) 
      t2b = trugud(u1,g1,u2,norbs,nclose,norbs) 
      t3b = trugdu(u1,g2,u1,norbs,nclose,norbs) 
      term = 2.0D0*(t1b - t1a + t2b - t2a + t3a - t3b) 
      return  
      end subroutine betcom 


      subroutine bmakuf(ua, ub, uab, t, uold1, gab, eigs, last, norbs, nclose, &
        diff, iwflb, maxu, btol) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use funcon_C, only : ev
      USE polar_C, only : omega 
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: norbs 
      integer , intent(in) :: nclose 
      integer , intent(in) :: iwflb 
      real(double) , intent(out) :: diff 
      real(double) , intent(out) :: maxu 
      real(double) , intent(in) :: btol 
      logical , intent(out) :: last 
      real(double) , intent(in) :: ua(norbs,norbs) 
      real(double) , intent(in) :: ub(norbs,norbs) 
      real(double) , intent(inout) :: uab(norbs,norbs) 
      real(double) , intent(in) :: t(norbs,norbs) 
      real(double) , intent(inout) :: uold1(norbs,norbs) 
      real(double) , intent(in) :: gab(norbs,norbs) 
      real(double) , intent(in) :: eigs(norbs) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, kll, kul, k, l 
      real(double) :: hartr, sum, udif 
!-----------------------------------------------
!
!  THIS SUBROUTINE CREATES THE NEW TRANSFORMATION MATRIX U
!  AND THEN CHECKS FOR CONVERGENCE
!
      hartr = ev 
!
!  ZERO MATRIX INITIALLY
!      CALL ZEROM(UAB,NORBS)
!
!  CREATE DIAGONAL BLOCKS (OCC,OCC) AND (UNOCC,UNOCC)
!
      do i = 1, norbs 
        do j = 1, i 
          sum = 0.0D00 
          if (i <= nclose) then 
            kll = nclose + 1 
            kul = norbs 
          else if (i>nclose .and. j>nclose) then 
            kll = 1 
            kul = nclose 
          endif 
          do k = kll, kul 
            sum = sum + ua(i,k)*ub(k,j) + ub(i,k)*ua(k,j) 
          end do 
          uab(i,j) = sum*0.5D00 
          uab(j,i) = sum*0.5D00 
        end do 
      end do 
!
!  CREATE OFF-DIAGONAL BLOCKS
!
      do k = nclose + 1, norbs 
        do l = 1, nclose 
          select case (iwflb)  
! CALCULATE FOR (W,W) VALUES
!
          case default 
            uab(k,l) = hartr*((gab(k,l)+t(k,l))/((eigs(l)-eigs(k))-2.0D00*omega&
              )) 
            uab(l,k) = hartr*((gab(l,k)+t(l,k))/((eigs(k)-eigs(l))-2.0D00*omega&
              )) 
! CALCULATE FOR (0,W) VALUES
!
          case (2)  
            uab(k,l) = hartr*((gab(k,l)+t(k,l))/((eigs(l)-eigs(k))-omega)) 
            uab(l,k) = hartr*((gab(l,k)+t(l,k))/((eigs(k)-eigs(l))-omega)) 
! CALCULATE FOR (W,-W) VALUES
!
          case (3)  
            uab(k,l) = hartr*((gab(k,l)+t(k,l))/(eigs(l)-eigs(k))) 
            uab(l,k) = hartr*((gab(l,k)+t(l,k))/(eigs(k)-eigs(l))) 
          end select  
        end do 
      end do 
!
!  CHECK FOR CONVERGENCE
!
      diff = 0.0D00 
      maxu = -1000 
      do i = 1, norbs 
        do j = 1, norbs 
          udif = uab(i,j) - uold1(i,j) 
          diff = dmax1(abs(udif),diff) 
          maxu = dmax1(uab(i,j),maxu) 
        end do 
      end do 
      if (diff < btol) last = .TRUE. 
!
      uold1 = uab 
!
      return  
      end subroutine bmakuf 


      subroutine copym(h, f, m) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: m 
      real(double) , intent(in) :: h(m,m) 
      real(double) , intent(out) :: f(m,m) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!  COPYM COPIES MATRIX H INTO F
!

! GBR_new_addition
! ORG      f = h 
      call mkl_domatcopy('c', 'n', m, m, 1.d0, h, m, f, m)
      !call dlacpy('u', m, m, h, m, f, m )
!        
      return  
      end subroutine copym 


      subroutine darea1(v, len, idaf, ns) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: len 
      integer , intent(in) :: idaf 
      integer , intent(in) :: ns 
      real(double)  :: v(len) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!       READ A PHYSICAL RECORD FROM THE DAF
!
      read (unit=idaf, rec=ns) v 
      return  
      end subroutine darea1 


      subroutine daread(v, len, nrec) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use chanel_C, only : iw, idaf, ioda, irecln
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use darea1_I 
      use mopend_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: len 
      integer , intent(in) :: nrec 
      real(double)  :: v(len) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: n, is, ns, lent, if, nsp, lenw 
!-----------------------------------------------
!
      n = ioda(nrec) 
      if (n /= (-1)) then 
        is = (-irecln) + 1 
        ns = n 
        lent = len 
        is = is + irecln 
        if = is + lent - 1 
        if (if - is + 1 > irecln) if = is + irecln - 1 
        nsp = ns 
        lenw = if - is + 1 
        call darea1 (v(is), lenw, idaf, nsp) 
        lent = lent - irecln 
        ns = ns + 1 
        n = ns 
        do while(lent >= 1) 
          is = is + irecln 
          if = is + lent - 1 
          if (if - is + 1 > irecln) if = is + irecln - 1 
          nsp = ns 
          lenw = if - is + 1 
          call darea1 (v(is), lenw, idaf, nsp) 
          lent = lent - irecln 
          ns = ns + 1 
          n = ns 
        end do 
        return  
      endif 
!
      write (iw, 30) nrec, len 
      call mopend (&
        '*** ERROR ***, ATTEMPT TO READ A DAF RECORD THAT WAS NEVER WRITTEN') 
!
   30 format(1x,'*** ERROR ***, ATTEMPT TO READ A DAF RECORD',&
        ' THAT WAS NEVER WRITTEN. NREC,LEN=',i5,i10) 
      return  
      end subroutine daread 


      subroutine dawrit(v, len, nrec) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use chanel_C, only : ioda, ifilen, irecst, idaf, irecln, iw
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dawrt1_I 
      use mopend_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: len 
      integer , intent(in) :: nrec 
      real(double)  :: v(len) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: n, ist, ns, lent, if, nsp, lenw 
      logical :: newrec 
!-----------------------------------------------
      n = ioda(nrec) 
      if (n<=0 .or. len==ifilen(nrec)) then 
        newrec = .FALSE. 
        if (n <= 0) then 
          ioda(nrec) = irecst 
          ifilen(nrec) = len 
          newrec = .TRUE. 
          irecst = irecst + (len - 1)/irecln + 1 
          n = ioda(nrec) 
        endif 
        ist = (-irecln) + 1 
        ns = n 
        lent = len 
        ist = ist + irecln 
        if = ist + lent - 1 
        if (if - ist + 1 > irecln) if = ist + irecln - 1 
        nsp = ns 
        lenw = if - ist + 1 
        call dawrt1 (v(ist), lenw, idaf, nsp) 
        lent = lent - irecln 
        ns = ns + 1 
        n = ns 
        do while(lent >= 1) 
          ist = ist + irecln 
          if = ist + lent - 1 
          if (if - ist + 1 > irecln) if = ist + irecln - 1 
          nsp = ns 
          lenw = if - ist + 1 
          call dawrt1 (v(ist), lenw, idaf, nsp) 
          lent = lent - irecln 
          ns = ns + 1 
          n = ns 
        end do 
        if (newrec) write (unit=idaf, rec=1) irecst, ioda, ifilen 
        return  
      endif 
!
      write (iw, 40) nrec, len, ifilen(nrec) 
      call mopend (&
      'DAWRIT HAS REQUESTED A RECORD WITH LENGTH DIFFERENT THAN BEFORE - ABORT &
      &FORCED.') 
!
   40 format(1x,'DAWRIT HAS REQUESTED A RECORD WITH LENGTH',1x,&
        'DIFFERENT THAN BEFORE - ABORT FORCED.'/,1x,'DAF RECORD ',i5,&
        ' NEW LENGTH =',i5,' OLD LENGTH =',i5) 
      return  
      end subroutine dawrit 


!*MODULE IOLIB   *DECK DAWRT1
      subroutine dawrt1(v, len, idaf, ns) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: len 
      integer , intent(in) :: idaf 
      integer , intent(in) :: ns 
      real(double) , intent(in) :: v(len) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!     ----- WRITE A PHYSICAL RECORD ON THE DAF -----
!
      write (unit=idaf, rec=ns) v 
      return  
      end subroutine dawrt1 


      subroutine densf(u, c, d, da, norbs, nclose, w2) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: norbs 
      integer , intent(in) :: nclose 
      real(double) , intent(in) :: u(norbs,norbs) 
      real(double) , intent(in) :: c(norbs,norbs) 
      real(double) , intent(out) :: d(norbs,norbs) 
      real(double) , intent(out) :: da(norbs,norbs) 
      real(double) , intent(inout) :: w2(norbs,norbs) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, k, l, i 
      real(double), dimension(norbs) :: w1 
      real(double) :: sum 
!-----------------------------------------------
!
!  THIS SUBROUTINE IS USED TO COMPUTE THE FIRST-ORDER DENSITY
!  FROM CA = C*U
!
!
!  FORM DENSITY MATRIX      CA*N*C+ + C*N*CA+
!
      do j = 1, norbs 
        do k = 1, norbs 
          sum = 0.D0 
          do l = 1, nclose 
            sum = sum + u(k,l)*c(j,l) 
          end do 
          w2(k,j) = sum 
        end do 
      end do 
      do i = 1, norbs 
        do k = 1, norbs 
          sum = 0.D0 
          do l = 1, nclose 
            sum = sum + c(i,l)*u(l,k) 
          end do 
          w1(k) = sum 
        end do 
        do j = 1, norbs 
          sum = 0.0D00 
          do k = 1, norbs 
            sum = sum + c(i,k)*w2(k,j) - w1(k)*c(j,k) 
          end do 
          d(i,j) = 2.0D00*sum 
          da(i,j) = sum 
        end do 
      end do 
!
      return  
      end subroutine densf 


      subroutine epsab(eigsab, eigs, gab, ga, gb, ua, ub, uab, udms, norbs, &
        nclose, iwflb) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use funcon_C, only : ev
      USE polar_C, only : omega 
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use zerom_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: norbs 
      integer , intent(in) :: nclose 
      integer , intent(in) :: iwflb 
      real(double)  :: eigsab(norbs,norbs) 
      real(double) , intent(in) :: eigs(norbs) 
      real(double) , intent(in) :: gab(norbs,norbs) 
      real(double) , intent(in) :: ga(norbs,norbs) 
      real(double) , intent(in) :: gb(norbs,norbs) 
      real(double) , intent(in) :: ua(norbs,norbs) 
      real(double) , intent(in) :: ub(norbs,norbs) 
      real(double) , intent(in) :: uab(norbs,norbs) 
      real(double)  :: udms(norbs,norbs) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j
      real(double) :: hartr, omval, s1, s2 
!-----------------------------------------------
!
!  THIS SUBROUTINE CREATES THE NEW EPSILON MATRIX AND UDMS MATRIX
!
      hartr = ev 
!
!  ZERO EPSILON OMEGA OMEGA MATRIX INITIALLY
!
      call zerom (eigsab, norbs) 
!
!  ZERO UAB MINUS OMEGA,OMEGA MATRIX INITIALLY
!
      call zerom (udms, norbs) 
!
      if (iwflb==0 .or. iwflb==1) then 
        omval = 2.0D00*omega 
      else if (iwflb == 3) then 
        omval = 0.0D00 
      else if (iwflb == 2) then 
        omval = omega 
      endif 
      do i = 1, nclose 
        do j = 1, nclose 
          s1 = 0.0D00 
!
! CALCULATION FOR EPSAB
!
          s1 = sum(ga(i,nclose+1:norbs)*ub(nclose+1:norbs,j)+gb(i,nclose+1:&
            norbs)*ua(nclose+1:norbs,j)) 
!
          eigsab(i,j) = gab(i,j) + s1 + uab(i,j)*(eigs(i)-eigs(j)+omval)/hartr 
        end do 
      end do 
!
! CALCULATION FOR UMS
!
      do i = 1, norbs 
        do j = 1, norbs 
          s2 = 0.0D00 
          s2 = sum(ua(i,:norbs)*ub(:norbs,j)+ub(i,:norbs)*ua(:norbs,j)) 
!
          udms(i,j) = s2 - uab(i,j) 
        end do 
      end do 
!
      return  
      end subroutine epsab 


      subroutine ffreq1(f, ptot, pa, pb, ndim) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      USE parameters_C, only : gss, gsp, gpp, gp2, hsp
      use common_arrays_C, only : nat, nfirst, nlast
      use molkst_C, only : numat
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ndim 
      real(double) , intent(inout) :: f(ndim,ndim) 
      real(double) , intent(in) :: ptot(ndim,ndim) 
      real(double) , intent(in) :: pa(ndim,ndim) 
      real(double) , intent(in) :: pb(ndim,ndim) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ii, ia, ib, ni, iplus, j, iminus, icc
      real(double) :: ptpop, papop 
      logical :: first 

      save first 
!-----------------------------------------------
! *********************************************************************
!
! *** COMPUTE THE REMAINING CONTRIBUTIONS TO THE ONE-CENTRE ELEMENTS.
!
! *********************************************************************
      data first/ .TRUE./  
      if (first) first = .FALSE. 
      do ii = 1, numat 
        ia = nfirst(ii) 
        ib = min(nlast(ii), ia + 3)
        ni = nat(ii) 
        ptpop = 0.D0 
        papop = 0.D0 
        go to (1002,20,10,10,10,10,10,10,10,10) ib - ia + 2 
   10   continue 
        ptpop = ptot(ib,ib) + ptot(ib-1,ib-1) + ptot(ib-2,ib-2) 
        papop = pa(ib,ib) + pa(ib-1,ib-1) + pa(ib-2,ib-2) 
   20   continue 
!
!     F(S,S)
!
        f(ia,ia) = f(ia,ia) + pb(ia,ia)*gss(ni) + ptpop*gsp(ni) - papop*hsp(ni) 
        if (ni >= 3) then 
          iplus = ia + 1 
          do j = iplus, ib 
!
!     F(P,P)
!
            f(j,j) = f(j,j) + ptot(ia,ia)*gsp(ni) - pa(ia,ia)*hsp(ni) + pb(j,j)&
              *gpp(ni) + (ptpop - ptot(j,j))*gp2(ni) - 0.5D0*(papop - pa(j,j))*&
              (gpp(ni)-gp2(ni)) 
!
!     F(S,P)
!
            f(ia,j) = f(ia,j) + 2.D0*ptot(ia,j)*hsp(ni) - pa(ia,j)*(hsp(ni)+gsp&
              (ni)) 
            f(j,ia) = f(j,ia) + 2.D0*ptot(j,ia)*hsp(ni) - pa(j,ia)*(hsp(ni)+gsp&
              (ni)) 
          end do 
!
!     F(P,P*)
!
          iminus = ib - 1 
          do j = iplus, iminus 
            icc = j + 1 
            f(j,icc:ib) = f(j,icc:ib) + ptot(j,icc:ib)*(gpp(ni)-gp2(ni)) - &
              0.5D0*pa(j,icc:ib)*(gpp(ni)+gp2(ni)) 
            f(icc:ib,j) = f(icc:ib,j) + ptot(icc:ib,j)*(gpp(ni)-gp2(ni)) - &
              0.5D0*pa(icc:ib,j)*(gpp(ni)+gp2(ni)) 
          end do 
        endif 
 1002   continue 
      end do 
      return  
      end subroutine ffreq1 


      subroutine ffreq2(f, ptot, w) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      USE molkst_C, only : numat, norbs
      use common_arrays_C, only : nfirst, nlast
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double) , intent(inout) :: f(norbs,norbs) 
      real(double) , intent(in) :: ptot(norbs,norbs) 
      real(double) , intent(in) :: w(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ii, iim1, ia, ib, jj, ja, jb, i, j, k, l, kk 
      real(double) :: fij, fkl, a, aint 
!-----------------------------------------------
!*******************************************************************
!
!  TDHF FORMS TWO ELECTRON TWO CENTER REPULSION PART OF THE FOCK
!  MATRIX
! ON INPUT PTOT = TOTAL DENSITY MATRIX
!          P    = ALPHA OR BETA DENSITY MATRIX
!          W    = TWO ELECTRON INTEGRAL MATRIX
!
!  ON OUTPUT F = PARTIAL FOCK MATRIX
!
!********************************************************************
  
        kk = 0 
!
        do ii = 1, numat 
          iim1 = ii - 1 
          ia = nfirst(ii) 
          ib = nlast(ii) 
          do jj = 1, iim1 
            ja = nfirst(jj) 
            jb = nlast(jj) 
            do i = ia, ib 
              do j = ia, i 
                fij = 1.0D00 
                if (i == j) fij = 0.5D00 
                do k = ja, jb 
                  do l = ja, k 
                    fkl = 1.0D00 
                    if (k == l) fkl = 0.5D00 
                    kk = kk + 1 
                    a = w(kk) 
                    aint = a*fkl*fij 
                    f(i,j) = f(i,j) + aint*(ptot(k,l)+ptot(l,k)) 
                    f(j,i) = f(j,i) + aint*(ptot(k,l)+ptot(l,k)) 
                    f(k,l) = f(k,l) + aint*(ptot(i,j)+ptot(j,i)) 
                    f(l,k) = f(l,k) + aint*(ptot(i,j)+ptot(j,i)) 
                    aint = aint*0.5D00 
                    f(i,l) = f(i,l) - aint*ptot(j,k) 
                    f(l,i) = f(l,i) - aint*ptot(k,j) 
                    f(k,j) = f(k,j) - aint*ptot(l,i) 
                    f(j,k) = f(j,k) - aint*ptot(i,l) 
                    f(i,k) = f(i,k) - aint*ptot(j,l) 
                    f(k,i) = f(k,i) - aint*ptot(l,j) 
                    f(j,l) = f(j,l) - aint*ptot(i,k) 
                    f(l,j) = f(l,j) - aint*ptot(k,i) 
                  end do 
                end do 
              end do 
            end do 
          end do 
          kk = kk + (((ib - ia +1 )*(ib - ia + 2))/2) **2
        end do 
      return  
      end subroutine ffreq2 


      subroutine fhpatn(a, b, norbs, itw, sign) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: norbs 
      integer , intent(in) :: itw 
      real(double) , intent(in) :: sign 
      real(double) , intent(out) :: a(norbs,norbs) 
      real(double) , intent(in) :: b(norbs,norbs) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
! THIS SUBROUTINE CONVERTS THE MATRICES INTO ITS ADJOINTS
!
      go to (10,40,40,10) itw 
   10 continue 
      a = b 
      go to 70 
   40 continue 
      a = transpose(sign*b) 
   70 continue 
      return  
      end subroutine fhpatn 


      subroutine hmuf(h1, id, coord, nfirst, nlast, nat, norbs, numat) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      USE parameters_C, only : dd 
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use zerom_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: id 
      integer  :: norbs 
      integer , intent(in) :: numat 
      integer , intent(in) :: nfirst(numat) 
      integer , intent(in) :: nlast(numat) 
      integer , intent(in) :: nat(numat) 
      real(double)  :: h1(norbs,norbs) 
      real(double) , intent(in) :: coord(3,*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ia, ib, ni, i1, j1, io1, jo1 
!-----------------------------------------------
!
!  FORM THE DIPOLE MOMENT MATRIX FOR COMPONENT ID
!
!
!  ZERO H1 MATRIX
!
      call zerom (h1, norbs) 
!
!  FORM DIPOLE MATRIX
!
      do i = 1, numat 
        ia = nfirst(i) 
        ib = min(nlast(i), ia + 3)
        ni = nat(i) 
        do i1 = ia, ib 
          do j1 = ia, i1 
            h1(i1,j1) = 0.0D00 
            io1 = i1 - ia 
            jo1 = j1 - ia 
            if (id==1 .and. jo1==0 .and. io1==1) then 
              h1(i1,j1) = dd(ni) 
              h1(j1,i1) = dd(ni) 
            endif 
            if (id==2 .and. jo1==0 .and. io1==2) then 
              h1(i1,j1) = dd(ni) 
              h1(j1,i1) = dd(ni) 
            endif 
            if (.not.(id==3 .and. jo1==0 .and. io1==3)) cycle  
            h1(i1,j1) = dd(ni) 
            h1(j1,i1) = dd(ni) 
          end do 
          h1(i1,i1) = 0.0D00 
!.. ADDED FOR TRANSLATION OF CENTER FROM ORIGIN
          h1(i1,i1) = h1(i1,i1) + 1.8897262D0*coord(id,i) 
        end do 
      end do 
!
      return  
      end subroutine hmuf 


      subroutine hplusf(f, h, norbs) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use funcon_C, only : ev
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: norbs 
      real(double) , intent(inout) :: f(norbs,norbs) 
      real(double) , intent(in) :: h(norbs,norbs) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double) :: hartr
!-----------------------------------------------
!
! HPLUSF ADDS THE 1 AND 2-ELECTRON PARTS OF THE FOCK MATRIX
!
      hartr = ev 
      f = h + f/hartr 
      return  
      end subroutine hplusf 


      subroutine makeuf(u, uold, g, eigs, last, norbs, nclose, diff, atol) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use funcon_C, only : ev
       USE polar_C, only : omega 
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use zerom_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: norbs 
      integer , intent(in) :: nclose 
      real(double) , intent(out) :: diff 
      real(double) , intent(in) :: atol 
      logical , intent(out) :: last 
      real(double)  :: u(norbs,norbs) 
      real(double) , intent(inout) :: uold(norbs,norbs) 
      real(double) , intent(in) :: g(norbs,norbs) 
      real(double) , intent(in) :: eigs(norbs) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, i, j 
      real(double) :: hartr, udif 
!-----------------------------------------------
!
!  THIS SUBROUTINE CREATES THE NEW TRANSFORMATION MATRIX U
!  AND THEN CHECKS FOR CONVERGENCE
!
      hartr = ev 
!
!  ZERO MATRIX INITIALLY
!
      call zerom (u, norbs) 
!
!  CREATE OFF-DIAGONAL BLOCKS
!
     do k = nclose + 1, norbs
      do i = 1, nclose
        u(i, k) = hartr * g(i, k) / (eigs(k)-eigs(i)-omega)
        u(k, i) = hartr * g(k, i) / (eigs(i)-eigs(k)-omega)
      end do
    end do
!
!  CHECK FOR CONVERGENCE
!
      diff = 0.0D00 
      do i = 1, norbs 
        do j = 1, norbs 
          udif = abs(u(i,j)-uold(i,j)) 
          diff = dmax1(udif,diff) 
        end do 
      end do 
      if (diff < atol) last = .TRUE. 
!
      uold(:norbs,:norbs) = u(:norbs,:norbs) 
!
      return  
      end subroutine makeuf 


      subroutine ngamtg(x, gd3, ud3, g1, u1, gs, usmd, eps, us) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use molkst_C, only : norbs, nclose
      use chanel_C, only : iw
       USE polar_C, only : omega 
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use daread_I 
      use fhpatn_I 
      use trsub_I 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double)  :: x(norbs,norbs) 
      real(double)  :: gd3(norbs,norbs) 
      real(double)  :: ud3(norbs,norbs) 
      real(double)  :: g1(norbs,norbs) 
      real(double)  :: u1(norbs,norbs) 
      real(double)  :: gs(norbs,norbs) 
      real(double)  :: usmd(norbs,norbs) 
      real(double)  :: eps(norbs,norbs) 
      real(double)  :: us(norbs,norbs)  
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(9) :: ida, idb, idc, idd 
      integer , dimension(3,3) :: ipair 
      integer :: msq, jgarc, juarc, jurec, jgrec, jg2rec, ju2rec, ju2mrc, &
        jeprec, ie, ia, ib, ic, id, icd, ibd, ibc, imove, j2, j34 
      real(double), dimension(9) :: gamma 
      real(double) :: one, gav, yy, gave 
      character, dimension(3) :: alab 

      save alab, ida, idb, idc, idd, ipair 
!-----------------------------------------------
!.....................................................................
!  CALCULATE GAMMA(THG) IN A NONITERATIVE FASHION
!..................................................................... 
      data alab/ 'X', 'Y', 'Z'/  
      data ida/ 1, 2, 3, 1, 1, 2, 2, 3, 3/  
      data idb/ 1, 2, 3, 1, 1, 2, 2, 3, 3/  
      data idc/ 1, 2, 3, 2, 3, 1, 3, 1, 2/  
      data idd/ 1, 2, 3, 2, 3, 1, 3, 1, 2/  
      data ipair/ 1, 2, 3, 2, 4, 5, 3, 5, 6/  
      one = 1.D0 
      msq = norbs*norbs 
      write (iw, 20) omega 
   20 format(/,/,' GAMMA (THIRD HARMONIC GENERATION) AT ',f10.5,' EV.'/,/) 
!
! IGAM=1 (THIRD HARMONIC GENERATION)
!
      jgarc = 22 
      juarc = 19 
      jurec = 7 
      jgrec = 10 
      jg2rec = 55 
      ju2rec = 49 
      ju2mrc = 67 
      jeprec = 61 
!
! LOOP BEGINS FOR THE CALCULATION OF GAMMA(ABCD)
!
      gav = 0.0D+00 
      do ie = 1, 9 
!
        ia = ida(ie) 
        ib = idb(ie) 
        ic = idc(ie) 
        id = idd(ie) 
        icd = ipair(ic,id) 
        ibd = ipair(ib,id) 
        ibc = ipair(ib,ic) 
!
!  READ IN THE FIRST ORDER U3 OMEGA AND G3 OMEGA IN THE DIRECTION A
!
! MAKE GD3 OMEGA MATRIX FROM G3 MATRIX
!
        call daread (x, msq, jgarc + ia) 
        call fhpatn (gd3, x, norbs, 2, one) 
!
! MAKE UD3 OMEGA MATRIX FROM U3 OMEGA MATRIX
!
        call daread (x, msq, juarc + ia) 
        call fhpatn (ud3, x, norbs, 2, (-one)) 
!
        yy = 0.0D00 
        imove = 1 
   30   continue 
!
        select case (imove)  
!
        case default 
          j2 = ib 
          j34 = icd 
        case (2)  
          j2 = ic 
          j34 = ibd 
        case (3)  
          j2 = id 
          j34 = ibc 
        end select 
!
!  READ IN G1,U1,GS,US,UMS,EPS
!
!  GET  UB
        call daread (u1, msq, jurec + j2) 
!  GET  GB
        call daread (g1, msq, jgrec + j2) 
!  GET  GCD
        call daread (gs, msq, jg2rec + j34) 
!  GET  UCD
        call daread (us, msq, ju2rec + j34) 
!  GET  USMD
        call daread (usmd, msq, ju2mrc + j34) 
!  GET  EPCD
        call daread (eps, msq, jeprec + j34) 
!
!
! FIRST KIND
!
        yy = yy + trsub(ud3,g1,us,nclose,norbs,norbs) 
        yy = yy - trsub(usmd,g1,ud3,nclose,norbs,norbs) 
        yy = yy - trsub(ud3,g1,us,norbs,nclose,norbs) 
        yy = yy + trsub(usmd,g1,ud3,norbs,nclose,norbs) 
!
! SECOND KIND
!
        yy = yy + trsub(ud3,gs,u1,nclose,norbs,norbs) 
        yy = yy + trsub(u1,gs,ud3,nclose,norbs,norbs) 
        yy = yy - trsub(ud3,eps,u1,norbs,nclose,norbs) 
        yy = yy - trsub(u1,eps,ud3,norbs,nclose,norbs) 
!
! THIRD KIND
!
        yy = yy + trsub(u1,gd3,us,nclose,norbs,norbs) 
        yy = yy - trsub(usmd,gd3,u1,nclose,norbs,norbs) 
        yy = yy - trsub(u1,gd3,us,norbs,nclose,norbs) 
        yy = yy + trsub(usmd,gd3,u1,norbs,nclose,norbs) 
!
        imove = imove + 1 
        if (imove <= 3) go to 30 
!
        gamma(ie) = yy 
!
! CALCULATE THE AVERAGE GAMMA VALUE
!
        gav = gav + yy 
!
! WRITE GAMMA(ABCD)
!
        write (iw,"(' GAMMA(',a1,',',a1,',',a1,',',a1,') = ',f19.5)") &
        alab(ia), alab(ib), alab(ic), alab(id), gamma(ie) 
      end do 
      gave = gav/5.0D+00 
      write (iw, "(/,/,' AVERAGE GAMMA VALUE AT ',f10.5,' = ',f19.5,a,f19.5,a,/,/)") &
      omega, gave, " a.u. =", gave/1.985425939776565, " ESU (X10-39)"
      return  
      end subroutine ngamtg 


      subroutine ngefis(x, gd3, ud3, g1, u1, gs, usmd, eps, us) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use molkst_C, only : norbs, nclose
      use chanel_C, only : iw
       USE polar_C, only : omega 
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use daread_I 
      use fhpatn_I 
      use trsub_I 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double)  :: x(norbs,norbs) 
      real(double)  :: gd3(norbs,norbs) 
      real(double)  :: ud3(norbs,norbs) 
      real(double)  :: g1(norbs,norbs) 
      real(double)  :: u1(norbs,norbs) 
      real(double)  :: gs(norbs,norbs) 
      real(double)  :: usmd(norbs,norbs) 
      real(double)  :: eps(norbs,norbs) 
      real(double)  :: us(norbs,norbs) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(15) :: ida, idb, idc, idd 
      integer , dimension(3,3) :: ip, ipair 
      integer :: msq, jgarc, juarc, jurec, jgrec, jg2rec, ju2rec, ju2mrc, &
        jeprec, ie, ia, ib, ic, id, icd, ibd, ibc, imove, j2, j34, j3u, j3g, &
        j3e, j3um 
      real(double), dimension(15) :: gamma 
      real(double) :: one, gav, yy, gave 
      character, dimension(3) :: alab 

      save alab, ida, idb, idc, idd, ip, ipair 
!-----------------------------------------------
!.....................................................................
!  CALCULATE GAMMA(DC-EFISHG) IN A NONITERATIVE FASHION
!..................................................................... 
      data alab/ 'X', 'Y', 'Z'/  
!
      data ida/ 1, 2, 3, 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3/  
      data idb/ 1, 2, 3, 2, 3, 1, 3, 1, 2, 1, 1, 2, 2, 3, 3/  
      data idc/ 1, 2, 3, 1, 1, 2, 2, 3, 3, 2, 3, 1, 3, 1, 2/  
      data idd/ 1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 3, 1, 3, 1, 2/  
      data ip/ 1, 2, 3, 2, 4, 5, 3, 5, 6/  
      data ipair/ 1, 4, 7, 2, 5, 8, 3, 6, 9/  
      one = 1.D0 
      msq = norbs*norbs 
      write (iw, 10) omega 
   10 format(/,/,' GAMMA (DC-EFISHG) AT ',f10.5,' EV.'/,/) 
!
! GET DATA FROM ALPHA  AND ITERATIVE BETA CALCULATIONS
!
!   REQUIRED RECORDS FROM POLARIZABILITY CALCULATIONS
!   -------------------------------------------------------
!        0    W    2W    3W
!
!       -02- -08-  -14-  -20-  -U- MATRIX FOR -X- DIRECTION
!       -03- -09-  -15-  -21-  -U- MATRIX FOR -Y- DIRECTION
!       -04- -10-  -16-  -22-  -U- MATRIX FOR -Z- DIRECTION
!       -05- -11-  -17-  -23-  -G- MATRIX FOR -X- DIRECTION
!       -06- -12-  -18-  -24-  -G- MATRIX FOR -Y- DIRECTION
!       -07- -13-  -19-  -25-  -G- MATRIX FOR -Z- DIRECTION
!   -------------------------------------------------------
!      (0,0)    (W,W)    (0,W)    (W,-W)
!
!      -26-     -50-     -74-     -110-   -U- MATRIX FOR -XX- DIRECTION
!      -27-     -51-     -75-     -111-   -U- MATRIX FOR -XY- DIRECTION
!      -28-     -52-     -76-     -112-   -U- MATRIX FOR -XZ- DIRECTION
!                        -77-     -113-   -U- MATRIX DOE -YX- DIRECTION
!      -29-     -53-     -78-     -114-   -U- MATRIX FOR -YY- DIRECTION
!      -30-     -54-     -79-     -115-   -U- MATRIX FOR -YZ- DIRECTION
!                        -80-     -116-   -U- MATRIX FOR -ZX- DIRECTION
!                        -81-     -117-   -U- MATRIX FOR -ZY- DIRECTION
!      -31-     -55-     -82-     -118-   -U- MATRIX FOR -ZZ- DIRECTION
!   ------------------------------------------------------------------
!
!      -32-     -56-     -83-     -119-   -G- MATRIX FOR -XX- DIRECTION
!      -33-     -57-     -84-     -120-   -G- MATRIX FOR -XY- DIRECTION
!      -34-     -58-     -85-     -121-   -G- MATRIX FOR -XZ- DIRECTION
!                        -86-     -122-   -G- MATRIX FOR -YX- DIRECTION
!      -35-     -59-     -87-     -123-   -G- MATRIX FOR -YY- DIRECTION
!      -36-     -60-     -88-     -124-   -G- MATRIX FOR -YZ- DIRECTION
!                        -89-     -125-   -G- MATRIX FOR -ZX- DIRECTION
!                        -90-     -126-   -G- MATRIX FOR -ZY- DIRECTION
!      -37-     -61-     -91-     -127-   -G- MATRIX FOR -ZZ- DIRECTION
!   ------------------------------------------------------------------
!
!      -38-     -62-     -92-     -128-   -E- MATRIX FOR -XX- DIRECTION
!      -39-     -63-     -93-     -129-   -E- MATRIX FOR -XY- DIRECTION
!      -40-     -64-     -94-     -130-   -E- MATRIX FOR -XZ- DIRECTION
!                        -95-     -131-   -E- MATRIX FOR -YX- DIRECTION
!      -41-     -65-     -96-     -132-   -E- MATRIX FOR -YY- DIRECTION
!      -42-     -66-     -97-     -133-   -E- MATRIX FOR -YZ- DIRECTION
!                        -98-     -134-   -E- MATRIX FOR -ZX- DIRECTION
!                        -99-     -135-   -E- MATRIX FOR -ZY- DIRECTION
!      -43-     -67-     -100-    -136-   -E- MATRIX FOR -ZZ- DIRECTION
!   ------------------------------------------------------------------
!
!      -44-     -68-     -101-    -137-   -UM- MATRIX FOR -XX- DIRECTION
!      -45-     -69-     -102-    -138-   -UM- MATRIX FOR -XY- DIRECTION
!      -46-     -70-     -103-    -139-   -UM- MATRIX FOR -XZ- DIRECTION
!                        -104-    -140-   -UM- MATRIX FOR -YX- DIRECTION
!      -47-     -71-     -105-    -141-   -UM- MATRIX FOR -YY- DIRECTION
!      -48-     -72-     -106-    -142-   -UM- MATRIX FOR -YZ- DIRECTION
!                        -107-    -143-   -UM- MATRIX FOR -ZX- DIRECTION
!                        -108-    -144-   -UM- MATRIX FOR -ZY- DIRECTION
!      -49-     -73-     -109-    -145-   -UM- MATRIX FOR -ZZ- DIRECTION
!   ------------------------------------------------------------------
!
! CALCULATION OF DIFFERENT GAMMA VALUES IN A NONITERATIVE METHOD.
!
! IGAM=2 (DC-ELECTIC FIELD INDUCED SECOND HARMONIC GENERATION)
!
      jgarc = 16 
      juarc = 13 
      jurec = 1 
      jgrec = 4 
      jg2rec = 55 
      ju2rec = 49 
      ju2mrc = 67 
      jeprec = 61 
! LOOP BEGINS FOR THE CALCULATION OF GAMMA(ABCD)
!
      gav = 0.0D00 
      do ie = 1, 15 
        ia = ida(ie) 
        ib = idb(ie) 
        ic = idc(ie) 
        id = idd(ie) 
        icd = ip(ic,id) 
        ibd = ipair(ib,id) 
        ibc = ipair(ib,ic) 
!
!  READ IN THE FIRST ORDER U3 OMEGA AND G3 OMEGA IN THE DIRECTION A
!  MAKE GD3 OMEGA MATRIX FROM G3 MATRIX
!
        call daread (x, msq, jgarc + ia) 
        call fhpatn (gd3, x, norbs, 2, one) 
!
! MAKE UD3 OMEGA MATRIX FROM U3 OMEGA MATRIX
!
        call daread (x, msq, juarc + ia) 
        call fhpatn (ud3, x, norbs, 2, (-one)) 
        yy = 0.0D00 
        imove = 1 
   20   continue 
!
! DC EFISHG
!
        select case (imove)  
        case default 
          j2 = ib 
          j34 = icd 
        case (2)  
          j2 = ic + 6 
!                 J34=IBD+24
          j3u = ibd + 24 
          j3g = ibd + 27 
          j3e = ibd + 30 
          j3um = ibd + 33 
        case (3)  
          j2 = id + 6 
          j3u = ibc + 24 
          j3g = ibc + 27 
          j3e = ibc + 30 
          j3um = ibc + 33 
        end select 
!
!  READ IN G1,U1,GS,US,UMS,EPS
!
!  CALL UB
!
        call daread (u1, msq, jurec + j2) 
!  CALL GB
        call daread (g1, msq, jgrec + j2) 
        if (imove == 1) then 
!  CALL GCD
          call daread (gs, msq, jg2rec + j34) 
!  CALL UCD
          call daread (us, msq, ju2rec + j34) 
!  CALL USMD
          call daread (usmd, msq, ju2mrc + j34) 
!  CALL EPCD
          call daread (eps, msq, jeprec + j34) 
!
        else 
!  CALL GCD
          call daread (gs, msq, jg2rec + j3g) 
!  CALL UCD
          call daread (us, msq, ju2rec + j3u) 
!  CALL USMD
          call daread (usmd, msq, ju2mrc + j3um) 
!  CALL EPCD
          call daread (eps, msq, jeprec + j3e) 
        endif 
!
! FIRST KIND
!
        yy = yy + trsub(ud3,g1,us,nclose,norbs,norbs) 
        yy = yy - trsub(usmd,g1,ud3,nclose,norbs,norbs) 
        yy = yy - trsub(ud3,g1,us,norbs,nclose,norbs) 
        yy = yy + trsub(usmd,g1,ud3,norbs,nclose,norbs) 
!
! SECOND KIND
!
        yy = yy + trsub(ud3,gs,u1,nclose,norbs,norbs) 
        yy = yy + trsub(u1,gs,ud3,nclose,norbs,norbs) 
        yy = yy - trsub(ud3,eps,u1,norbs,nclose,norbs) 
        yy = yy - trsub(u1,eps,ud3,norbs,nclose,norbs) 
!
! THIRD KIND
!
        yy = yy + trsub(u1,gd3,us,nclose,norbs,norbs) 
        yy = yy - trsub(usmd,gd3,u1,nclose,norbs,norbs) 
        yy = yy - trsub(u1,gd3,us,norbs,nclose,norbs) 
        yy = yy + trsub(usmd,gd3,u1,norbs,nclose,norbs) 
!
        imove = imove + 1 
        if (imove <= 3) go to 20 
!
        gamma(ie) = yy 
! CALCULATE THE AVERAGE GAMMA VALUE
        if (ie <= 3) then 
          gav = gav + 3*yy 
        else if (ie > 9) then 
          gav = gav + yy 
        else 
          gav = gav + 2*yy 
        endif 
!
! WRITE GAMMA(ABCD)
!
        write (iw, 80) alab(ia), alab(ib), alab(ic), alab(id), gamma(ie) 
   80   format(' GAMMA(',a1,',',a1,',',a1,',',a1,') = ',1p,d14.7) 
!
      end do 
      gave = gav/15.0D+00 
      write (iw, 100) omega, gave, " a.u. =", gave/1.985425939776565, " ESU (X10-39)"
  100 format(/,/,' AVERAGE GAMMA VALUE AT ',f10.5,' EV = ',1p,d14.7,a,d14.7,a,/,/) 
      return  
      end subroutine ngefis 


      subroutine ngidri(x, gd3, ud3, g1, u1, gs, usmd, eps, us) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use molkst_C, only : norbs, nclose
       USE polar_C, only : omega 
      use chanel_C, only : iw
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use daread_I 
      use fhpatn_I 
      use trsub_I 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double)  :: x(norbs,norbs) 
      real(double)  :: gd3(norbs,norbs) 
      real(double)  :: ud3(norbs,norbs) 
      real(double)  :: g1(norbs,norbs) 
      real(double)  :: u1(norbs,norbs) 
      real(double)  :: gs(norbs,norbs) 
      real(double)  :: usmd(norbs,norbs) 
      real(double)  :: eps(norbs,norbs) 
      real(double)  :: us(norbs,norbs) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(15) :: ida, idb, idc, idd 
      integer , dimension(3,3) :: ip, ipair 
      integer :: msq, jgarc, juarc, jurec, jgrec, ie, ia, ib, ic, id, icd, ibd&
        , ibc, imove, j2, j34, jg2rec, ju2rec, ju2mrc, jeprec 
      real(double), dimension(15) :: gamma 
      real(double) :: one, gav, yy, gave 
      character, dimension(3) :: alab 

      save alab, ida, idb, idc, idd, ip, ipair 
!-----------------------------------------------
!.....................................................................
!  CALCULATE GAMMA(IDRI) IN A NONITERATIVE FASHION
!.....................................................................
      data alab/ 'X', 'Y', 'Z'/  
!
      data ida/ 1, 2, 3, 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3/  
      data idb/ 1, 2, 3, 1, 1, 2, 2, 3, 3, 2, 3, 1, 3, 1, 2/  
      data idc/ 1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 3, 1, 3, 1, 2/  
      data idd/ 1, 2, 3, 2, 3, 1, 3, 1, 2, 1, 1, 2, 2, 3, 3/  
      data ip/ 1, 2, 3, 2, 4, 5, 3, 5, 6/  
      data ipair/ 1, 4, 7, 2, 5, 8, 3, 6, 9/  
      one = 1.0D00 
      msq = norbs*norbs 
!
      write (iw, 10) omega 
   10 format(/,/,' GAMMA (IDRI) AT ',f10.5,' EV.'/,/) 
!
! GET DATA FROM ALPHA  AND ITERATIVE BETA CALCULATIONS
!
!
! IGAM=3 (INTENSITY DEPENDENT REFRACTIVE INDEX OR DEGENERATED FOUR
! WAVE MIXING)
!
      jgarc = 10 
      juarc = 7 
      jurec = 7 
      jgrec = 10 
! LOOP BEGINS FOR THE CALCULATION OF GAMMA(ABCD)
!
      gav = 0.0D+00 
      do ie = 1, 15 
        ia = ida(ie) 
        ib = idb(ie) 
        ic = idc(ie) 
        id = idd(ie) 
        icd = ipair(ic,id) 
        ibd = ipair(ib,id) 
        ibc = ip(ib,ic) 
!
!  READ IN THE FIRST ORDER U3 OMEGA AND G3 OMEGA IN THE DIRECTION A
!
! MAKE GD3 OMEGA MATRIX FROM G3 MATRIX
!
        call daread (x, msq, jgarc + ia) 
        call fhpatn (gd3, x, norbs, 2, one) 
!
! MAKE UD3 OMEGA MATRIX FROM U3 OMEGA MATRIX
!
        call daread (x, msq, juarc + ia) 
        call fhpatn (ud3, x, norbs, 2, (-one)) 
!
        yy = 0.0D00 
        imove = 1 
   20   continue 
!
!
! IDRI
!
        select case (imove)  
        case default 
          j2 = ib 
          j34 = icd 
          jg2rec = 118 
          ju2rec = 109 
          ju2mrc = 136 
          jeprec = 127 
        case (2)  
          j2 = ic 
          j34 = ibd 
          jg2rec = 118 
          ju2rec = 109 
          ju2mrc = 136 
          jeprec = 127 
        case (3)  
          j2 = id 
          j34 = ibc 
          jg2rec = 55 
          ju2rec = 49 
          ju2mrc = 67 
          jeprec = 61 
        end select 
!  READ IN G1,U1,GS,US,UMS,EPS
!
!  CALL UB
        if (imove == 3) then 
          call daread (x, msq, jurec + j2) 
          call fhpatn (u1, x, norbs, 2, (-one)) 
        else 
          call daread (u1, msq, jurec + j2) 
        endif 
!  CALL GB
        if (imove == 3) then 
          call daread (x, msq, jgrec + j2) 
          call fhpatn (g1, x, norbs, 2, one) 
        else 
          call daread (g1, msq, jgrec + j2) 
        endif 
!  CALL GCD
        call daread (gs, msq, jg2rec + j34) 
!  CALL UCD
        call daread (us, msq, ju2rec + j34) 
!  CALL USMD
        call daread (usmd, msq, ju2mrc + j34) 
!  CALL EPCD
        call daread (eps, msq, jeprec + j34) 
!
! FIRST KIND
!
        yy = yy + trsub(ud3,g1,us,nclose,norbs,norbs) 
        yy = yy - trsub(usmd,g1,ud3,nclose,norbs,norbs) 
        yy = yy - trsub(ud3,g1,us,norbs,nclose,norbs) 
        yy = yy + trsub(usmd,g1,ud3,norbs,nclose,norbs) 
!
! SECOND KIND
!
        yy = yy + trsub(ud3,gs,u1,nclose,norbs,norbs) 
        yy = yy + trsub(u1,gs,ud3,nclose,norbs,norbs) 
        yy = yy - trsub(ud3,eps,u1,norbs,nclose,norbs) 
        yy = yy - trsub(u1,eps,ud3,norbs,nclose,norbs) 
!
! THIRD KIND
!
        yy = yy + trsub(u1,gd3,us,nclose,norbs,norbs) 
        yy = yy - trsub(usmd,gd3,u1,nclose,norbs,norbs) 
        yy = yy - trsub(u1,gd3,us,norbs,nclose,norbs) 
        yy = yy + trsub(usmd,gd3,u1,norbs,nclose,norbs) 
!
        imove = imove + 1 
        if (imove <= 3) go to 20 
!
        gamma(ie) = yy 
!
! CALCULATE THE AVERAGE GAMMA VALUE
!
        if (ie <= 3) then 
          gav = gav + 3.0D0*yy 
        else if (ie > 9) then 
          gav = gav + yy 
        else 
          gav = gav + 2.0D0*yy 
        endif 
!
! WRITE GAMMA(ABCD)
!
        write (iw, 70) alab(ia), alab(ib), alab(ic), alab(id), gamma(ie) 
   70   format(' GAMMA(',a1,',',a1,',',a1,',',a1,') = ',1p,d14.7) 
!
      end do 
      gave = gav/15.0D+00 
      write (iw, 90) omega, gave, " a.u. =", gave/1.985425939776565, " ESU (X10-39)"
   90 format(/,/,'  AVERAGE GAMMA VALUE AT ',f10.5,' = ',1p,d14.7,a,d14.7,a,/,/) 
      return  
      end subroutine ngidri 


      subroutine ngoke(igam, x, gd3, ud3, g1, u1, gs, usmd, eps, us) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use molkst_C, only : norbs, nclose
      use chanel_C, only : iw
       USE polar_C, only : omega 
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use daread_I 
      use fhpatn_I 
      use trsub_I 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: igam 
      real(double)  :: x(norbs,norbs) 
      real(double)  :: gd3(norbs,norbs) 
      real(double)  :: ud3(norbs,norbs) 
      real(double)  :: g1(norbs,norbs) 
      real(double)  :: u1(norbs,norbs) 
      real(double)  :: gs(norbs,norbs) 
      real(double)  :: usmd(norbs,norbs) 
      real(double)  :: eps(norbs,norbs) 
      real(double)  :: us(norbs,norbs)  
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(15) :: ida, idb, idc, idd 
      integer , dimension(3,3) :: ip, ipair 
      integer :: msq, jgarc, juarc, jurec, jgrec, ie, ia, ib, ic, id, icd, ibd&
        , ibc, imove, j2, j34, jg2rec, ju2rec, ju2mrc, jeprec 
      real(double), dimension(15) :: gamma 
      real(double) :: one, gav, yy, gave 
      character, dimension(3) :: alab 

      save alab, ida, idb, idc, idd, ip, ipair 
!-----------------------------------------------
!.....................................................................
!  CALCULATE GAMMA(OKE) IN A NONITERATIVE FASHION
!.....................................................................
      data alab/ 'X', 'Y', 'Z'/  
!
      data ida/ 1, 2, 3, 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3/  
      data idb/ 1, 2, 3, 1, 1, 2, 2, 3, 3, 2, 3, 1, 3, 1, 2/  
      data idc/ 1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 3, 1, 3, 1, 2/  
      data idd/ 1, 2, 3, 2, 3, 1, 3, 1, 2, 1, 1, 2, 2, 3, 3/  
      data ip/ 1, 2, 3, 2, 4, 5, 3, 5, 6/  
      data ipair/ 1, 4, 7, 2, 5, 8, 3, 6, 9/  
      one = 1.0D00 
      msq = norbs*norbs 
!
      if (igam == 3) then 
        write (iw, 10) omega 
   10   format(/,/,' GAMMA (IDRI) AT ',f10.5,' EV.'/,/) 
      else 
        write (iw, 20) omega 
   20   format(/,/,' GAMMA (OKE) AT ',f10.5,' EV.'/,/) 
      endif 
!
! DATA INCLUDING YX, ZY, ZX DIRECTIONS
! GET DATA FROM ALPHA  AND ITERATIVE BETA CALCULATIONS
!
!   REQUIRED RECORDS FROM POLARIZABILITY CALCULATIONS
!
!   -------------------------------------------------------
!        0    W    2W    3W
!
!       -02- -08-  -14-  -20-  -U- MATRIX FOR -X- DIRECTION
!       -03- -09-  -15-  -21-  -U- MATRIX FOR -Y- DIRECTION
!       -04- -10-  -16-  -22-  -U- MATRIX FOR -Z- DIRECTION
!       -05- -11-  -17-  -23-  -G- MATRIX FOR -X- DIRECTION
!       -06- -12-  -18-  -24-  -G- MATRIX FOR -Y- DIRECTION
!       -07- -13-  -19-  -25-  -G- MATRIX FOR -Z- DIRECTION
!   -------------------------------------------------------
!      (0,0)    (W,W)    (0,W)    (W,-W)
!
!      -26-     -50-     -74-     -110-   -U- MATRIX FOR -XX- DIRECTION
!      -27-     -51-     -75-     -111-   -U- MATRIX FOR -XY- DIRECTION
!      -28-     -52-     -76-     -112-   -U- MATRIX FOR -XZ- DIRECTION
!                        -77-     -113-   -U- MATRIX DOE -YX- DIRECTION
!      -29-     -53-     -78-     -114-   -U- MATRIX FOR -YY- DIRECTION
!      -30-     -54-     -79-     -115-   -U- MATRIX FOR -YZ- DIRECTION
!                        -80-     -116-   -U- MATRIX FOR -ZX- DIRECTION
!                        -81-     -117-   -U- MATRIX FOR -ZY- DIRECTION
!      -31-     -55-     -82-     -118-   -U- MATRIX FOR -ZZ- DIRECTION
!   ------------------------------------------------------------------
!
!      -32-     -56-     -83-     -119-   -G- MATRIX FOR -XX- DIRECTION
!      -33-     -57-     -84-     -120-   -G- MATRIX FOR -XY- DIRECTION
!      -34-     -58-     -85-     -121-   -G- MATRIX FOR -XZ- DIRECTION
!                        -86-     -122-   -G- MATRIX FOR -YX- DIRECTION
!      -35-     -59-     -87-     -123-   -G- MATRIX FOR -YY- DIRECTION
!      -36-     -60-     -88-     -124-   -G- MATRIX FOR -YZ- DIRECTION
!                        -89-     -125-   -G- MATRIX FOR -ZX- DIRECTION
!                        -90-     -126-   -G- MATRIX FOR -ZY- DIRECTION
!      -37-     -61-     -91-     -127-   -G- MATRIX FOR -ZZ- DIRECTION
!   ------------------------------------------------------------------
!
!      -38-     -62-     -92-     -128-   -E- MATRIX FOR -XX- DIRECTION
!      -39-     -63-     -93-     -129-   -E- MATRIX FOR -XY- DIRECTION
!      -40-     -64-     -94-     -130-   -E- MATRIX FOR -XZ- DIRECTION
!                        -95-     -131-   -E- MATRIX FOR -YX- DIRECTION
!      -41-     -65-     -96-     -132-   -E- MATRIX FOR -YY- DIRECTION
!      -42-     -66-     -97-     -133-   -E- MATRIX FOR -YZ- DIRECTION
!                        -98-     -134-   -E- MATRIX FOR -ZX- DIRECTION
!                        -99-     -135-   -E- MATRIX FOR -ZY- DIRECTION
!      -43-     -67-     -100-    -136-   -E- MATRIX FOR -ZZ- DIRECTION
!   ------------------------------------------------------------------
!
!      -44-     -68-     -101-    -137-   -UM- MATRIX FOR -XX- DIRECTION
!      -45-     -69-     -102-    -138-   -UM- MATRIX FOR -XY- DIRECTION
!      -46-     -70-     -103-    -139-   -UM- MATRIX FOR -XZ- DIRECTION
!                        -104-    -140-   -UM- MATRIX FOR -YX- DIRECTION
!      -47-     -71-     -105-    -141-   -UM- MATRIX FOR -YY- DIRECTION
!      -48-     -72-     -106-    -142-   -UM- MATRIX FOR -YZ- DIRECTION
!                        -107-    -143-   -UM- MATRIX FOR -ZX- DIRECTION
!                        -108-    -144-   -UM- MATRIX FOR -ZY- DIRECTION
!      -49-     -73-     -109-    -145-   -UM- MATRIX FOR -ZZ- DIRECTION
!   ------------------------------------------------------------------
!
! GET DATA FROM ALPHA  AND ITERATIVE BETA CALCULATIONS
!
!
! IGAM=4 (OPTICAL KERR EFFECT)
!
      jgarc = 10 
      juarc = 7 
      jurec = 1 
      jgrec = 4 
! LOOP BEGINS FOR THE CALCULATION OF GAMMA(ABCD)
!
      gav = 0.0D+00 
      do ie = 1, 15 
        ia = ida(ie) 
        ib = idb(ie) 
        ic = idc(ie) 
        id = idd(ie) 
        icd = ipair(ic,id) 
        ibd = ipair(ib,id) 
        ibc = ip(ib,ic) 
!
!  READ IN THE FIRST ORDER U3 OMEGA AND G3 OMEGA IN THE DIRECTION A
!
! MAKE GD3 OMEGA MATRIX FROM G3 MATRIX
!
        call daread (x, msq, jgarc + ia) 
        call fhpatn (gd3, x, norbs, 2, one) 
!
! MAKE UD3 OMEGA MATRIX FROM U3 OMEGA MATRIX
!
        call daread (x, msq, juarc + ia) 
        call fhpatn (ud3, x, norbs, 2, (-one)) 
!
        yy = 0.0D00 
        imove = 1 
   30   continue 
!
! OKE
        select case (imove)  
        case default 
          j2 = ib 
          j34 = icd 
          jg2rec = 82 
          ju2rec = 73 
          ju2mrc = 100 
          jeprec = 91 
        case (2)  
          j2 = ic 
          j34 = ibd 
          jg2rec = 82 
          ju2rec = 73 
          ju2mrc = 100 
          jeprec = 91 
        case (3)  
          j2 = id 
          j34 = ibc 
          jg2rec = 31 
          ju2rec = 25 
          ju2mrc = 43 
          jeprec = 37 
        end select 
!  READ IN G1,U1,GS,US,UMS,EPS
!
!  CALL UB
        if (imove == 3) then 
          call daread (u1, msq, juarc + j2) 
        else 
          call daread (u1, msq, jurec + j2) 
        endif 
!  CALL GB
        if (imove == 3) then 
          call daread (g1, msq, jgarc + j2) 
        else 
          call daread (g1, msq, jgrec + j2) 
        endif 
!  CALL GCD
        call daread (gs, msq, jg2rec + j34) 
!  CALL UCD
        call daread (us, msq, ju2rec + j34) 
!  CALL USMD
        call daread (usmd, msq, ju2mrc + j34) 
!  CALL EPCD
        call daread (eps, msq, jeprec + j34) 
!
! FIRST KIND
!
        yy = yy + trsub(ud3,g1,us,nclose,norbs,norbs) 
        yy = yy - trsub(usmd,g1,ud3,nclose,norbs,norbs) 
        yy = yy - trsub(ud3,g1,us,norbs,nclose,norbs) 
        yy = yy + trsub(usmd,g1,ud3,norbs,nclose,norbs) 
!
! SECOND KIND
!
        yy = yy + trsub(ud3,gs,u1,nclose,norbs,norbs) 
        yy = yy + trsub(u1,gs,ud3,nclose,norbs,norbs) 
        yy = yy - trsub(ud3,eps,u1,norbs,nclose,norbs) 
        yy = yy - trsub(u1,eps,ud3,norbs,nclose,norbs) 
!
! THIRD KIND
!
        yy = yy + trsub(u1,gd3,us,nclose,norbs,norbs) 
        yy = yy - trsub(usmd,gd3,u1,nclose,norbs,norbs) 
        yy = yy - trsub(u1,gd3,us,norbs,nclose,norbs) 
        yy = yy + trsub(usmd,gd3,u1,norbs,nclose,norbs) 
!
        imove = imove + 1 
        if (imove <= 3) go to 30 
!
        gamma(ie) = yy 
!
! CALCULATE THE AVERAGE GAMMA VALUE
!
        if (ie <= 3) then 
          gav = gav + 3.0D0*yy 
        else if (ie > 9) then 
          gav = gav + yy 
        else 
          gav = gav + 2.0D0*yy 
        endif 
!
! WRITE GAMMA(ABCD)
!
        write (iw, 80) alab(ia), alab(ib), alab(ic), alab(id), gamma(ie) 
   80   format(' GAMMA(',a1,',',a1,',',a1,',',a1,') = ',1p,d14.7) 
!
      end do 
      gave = gav/15.0D+00 
      write (iw, 100) omega, gave, " a.u. =", gave/1.985425939776565, " ESU (X10-39)"
  100 format(/,/,'  AVERAGE GAMMA VALUE AT ',f10.5,' = ',1p,d14.7,a,d14.7,a,/,/) 
      return  
      end subroutine ngoke 


      subroutine nonbet(u1x, u1y, u1z, u2x, u2y, u2z, g1x, g1y, g1z, g2x, g2y, g2z) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use molkst_C, only : norbs, nclose, line
      use chanel_C, only : iw
       USE polar_C, only : omega 
!
! THIS SUBROUTINE CALCULATES SECOND HARMONIC GENERATION IN A
! NONITERATIVE WAY.
!
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use daread_I 
      use betcom_I 
      use betall_I 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double)  :: u1x(norbs,norbs) 
      real(double)  :: u1y(norbs,norbs) 
      real(double)  :: u1z(norbs,norbs) 
      real(double)  :: u2x(norbs,norbs) 
      real(double)  :: u2y(norbs,norbs) 
      real(double)  :: u2z(norbs,norbs) 
      real(double)  :: g1x(norbs,norbs) 
      real(double)  :: g1y(norbs,norbs) 
      real(double)  :: g1z(norbs,norbs) 
      real(double)  :: g2x(norbs,norbs) 
      real(double)  :: g2y(norbs,norbs) 
      real(double)  :: g2z(norbs,norbs) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: maxsq 
      real(double) :: bavx, bavy, bavz, bxxx, byxx, bzxx, bxxy, byxy, bzxy, &
        bxxz, byxz, bzxz, bxyx, byyx, bzyx, bxyy, byyy, bzyy, bxyz, byyz, bzyz&
        , bxzx, byzx, bzzx, bxzy, byzy, bzzy, bxzz, byzz, bzzz, bvec 
!-----------------------------------------------
!.....................................................................
!  CALCULATE BETA IN A NONITERATIVE FASHION
!.....................................................................
!
! GET DATA FROM ALPHA CALCULATIONS
!
      maxsq = norbs*norbs 
      bavx = 0.0D+00 
      bavy = 0.0D+00 
      bavz = 0.0D+00 
      call daread (u1x, maxsq, 8) 
      call daread (u1y, maxsq, 9) 
      call daread (u1z, maxsq, 10) 
      call daread (g1x, maxsq, 11) 
      call daread (g1y, maxsq, 12) 
      call daread (g1z, maxsq, 13) 
      call daread (u2x, maxsq, 14) 
      call daread (u2y, maxsq, 15) 
      call daread (u2z, maxsq, 16) 
      call daread (g2x, maxsq, 17) 
      call daread (g2y, maxsq, 18) 
      call daread (g2z, maxsq, 19) 
! XXX
      call betcom (u1x, g1x, u2x, g2x, nclose, norbs, bxxx) 
      bavx = bavx + 3.0D0*bxxx 
! YXX
      call betcom (u1x, g1x, u2y, g2y, nclose, norbs, byxx) 
      bavy = bavy + byxx 
! ZXX
      call betcom (u1x, g1x, u2z, g2z, nclose, norbs, bzxx) 
      bavz = bavz + bzxx 
! XXY
      call betall (u2x, g2x, u1x, g1x, u1y, g1y, nclose, norbs, bxxy) 
      bavy = bavy + bxxy 
! YXY
      call betall (u2y, g2y, u1x, g1x, u1y, g1y, nclose, norbs, byxy) 
      bavx = bavx + byxy 
! ZXY
      call betall (u2z, g2z, u1x, g1x, u1y, g1y, nclose, norbs, bzxy) 
! XXZ
      call betall (u2x, g2x, u1x, g1x, u1z, g1z, nclose, norbs, bxxz) 
      bavz = bavz + bxxz 
! YXZ
      call betall (u2y, g2y, u1x, g1x, u1z, g1z, nclose, norbs, byxz) 
! ZXZ
      call betall (u2z, g2z, u1x, g1x, u1z, g1z, nclose, norbs, bzxz) 
      bavx = bavx + bzxz 
! XYX
      call betall (u2x, g2x, u1y, g1y, u1x, g1x, nclose, norbs, bxyx) 
      bavy = bavy + bxyx 
! YYX
      call betall (u2y, g2y, u1y, g1y, u1x, g1x, nclose, norbs, byyx) 
      bavx = bavx + byyx 
! ZYX
      call betall (u2z, g2z, u1y, g1y, u1x, g1x, nclose, norbs, bzyx) 
! XYY
      call betcom (u1y, g1y, u2x, g2x, nclose, norbs, bxyy) 
      bavx = bavx + bxyy 
! YYY
      call betcom (u1y, g1y, u2y, g2y, nclose, norbs, byyy) 
      bavy = bavy + 3.0D0*byyy 
! ZYY
      call betcom (u1y, g1y, u2z, g2z, nclose, norbs, bzyy) 
      bavz = bavz + bzyy 
! XYZ
      call betall (u2x, g2x, u1y, g1y, u1z, g1z, nclose, norbs, bxyz) 
! YYZ
      call betall (u2y, g2y, u1y, g1y, u1z, g1z, nclose, norbs, byyz) 
      bavz = bavz + byyz 
! ZYZ
      call betall (u2z, g2z, u1y, g1y, u1z, g1z, nclose, norbs, bzyz) 
      bavy = bavy + bzyz 
! XZX
      call betall (u2x, g2x, u1z, g1z, u1x, g1x, nclose, norbs, bxzx) 
      bavz = bavz + bxzx 
! YZX
      call betall (u2y, g2y, u1z, g1z, u1x, g1x, nclose, norbs, byzx) 
! ZZX
      call betall (u2z, g2z, u1z, g1z, u1x, g1x, nclose, norbs, bzzx) 
      bavx = bavx + bzzx 
! XZY
      call betall (u2x, g2x, u1z, g1z, u1y, g1y, nclose, norbs, bxzy) 
! YZY
      call betall (u2y, g2y, u1z, g1z, u1y, g1y, nclose, norbs, byzy) 
      bavz = bavz + byzy 
! ZZY
      call betall (u2z, g2z, u1z, g1z, u1y, g1y, nclose, norbs, bzzy) 
      bavy = bavy + bzzy 
! XZZ
      call betcom (u1z, g1z, u2x, g2x, nclose, norbs, bxzz) 
      bavx = bavx + bxzz 
! YZZ
      call betcom (u1z, g1z, u2y, g2y, nclose, norbs, byzz) 
      bavy = bavy + byzz 
! ZZZ
      call betcom (u1z, g1z, u2z, g2z, nclose, norbs, bzzz) 
      bavz = bavz + 3.0D0*bzzz 
!
      bavx = bavx/5.0D+00 
      bavy = bavy/5.0D+00 
      bavz = bavz/5.0D+00 
!
      bvec = (bavx*bavx + bavy*bavy + bavz*bavz)**0.5D+00 
      write (iw, 10) 
   10 format(/,/,' BETA (SECOND HARMONIC GENERATION)'/,/) 
      write (iw, 20) bxxx, byxx, bzxx, bxxy, byxy, bzxy, bxxz, byxz, bzxz, bxyx&
        , byyx, bzyx, bxyy, byyy, bzyy, bxyz, byyz, bzyz, bxzx, byzx, bzzx, &
        bxzy, byzy, bzzy, bxzz, byzz, bzzz 
   20 format(/,/,'  BXXX  ',d15.8,'  BYXX ',d15.8,'  BZXX ',d15.8,/,'  BXXY  ',&
        d15.8,'  BYXY ',d15.8,'  BZXY ',d15.8,/,'  BXXZ  ',d15.8,'  BYXZ ',&
        d15.8,'  BZXZ ',d15.8,/,'  BXYX  ',d15.8,'  BYYX ',d15.8,'  BZYX ',&
        d15.8,/,'  BXYY  ',d15.8,'  BYYY ',d15.8,'  BZYY ',d15.8,/,'  BXYZ  ',&
        d15.8,'  BYYZ ',d15.8,'  BZYZ ',d15.8,/,'  BXZX  ',d15.8,'  BYZX ',&
        d15.8,'  BZZX ',d15.8,/,'  BXZY  ',d15.8,'  BYZY ',d15.8,'  BZZY ',&
        d15.8,/,'  BXZZ  ',d15.8,'  BYZZ ',d15.8,'  BZZZ ',d15.8) 
!
      if (bvec < 1.d10) then
        line = "(' AVERAGE BETA',a,'(SHG) VALUE AT',f10.5,' EV = ',f11.4, a)"
      else if (bvec < 1.d15) then
        line = "(' AVERAGE BETA',a,'(SHG) VALUE AT',f10.5,' EV = ',f16.4, a)"
      else
        line = "(' AVERAGE BETA',a,'(SHG) VALUE AT',f10.5,' EV = ',f21.4, a)"
      end if
      write (iw,"(/)")
      write (iw, trim(line)) "X", omega, bavx, " a.u."
      write (iw, trim(line)) "Y", omega, bavy, " a.u."
      write (iw, trim(line)) "Z", omega, bavz, " a.u."
      write (iw,"(2/)")
      write (iw, trim(line)) " ", omega, bvec, " a.u."
      write (iw,"(/)")
      return  
      end subroutine nonbet 


      subroutine nonope(u0x, u1y, u1z, u1x, u0y, u0z, g0x, g1y, g1z, g1x, g0y, &
        g0z) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use molkst_C, only : norbs, nclose, line
      use chanel_C, only : iw 
      USE polar_C, only : omega 
!
! THIS SUBROUTINE CALCULATES ELECTROOPTIC POCKEL'S EFFECT
! IN A NONITERATIVE WAY.
!
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use daread_I 
      use betall_I 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double)  :: u0x(norbs,norbs) 
      real(double)  :: u1y(norbs,norbs) 
      real(double)  :: u1z(norbs,norbs) 
      real(double)  :: u1x(norbs,norbs) 
      real(double)  :: u0y(norbs,norbs) 
      real(double)  :: u0z(norbs,norbs) 
      real(double)  :: g0x(norbs,norbs) 
      real(double)  :: g1y(norbs,norbs) 
      real(double)  :: g1z(norbs,norbs) 
      real(double)  :: g1x(norbs,norbs) 
      real(double)  :: g0y(norbs,norbs) 
      real(double)  :: g0z(norbs,norbs) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: maxsq 
      real(double) :: bavx, bavy, bavz, bxxx, byxx, bzxx, bxxy, byxy, bzxy, &
        bxxz, byxz, bzxz, bxyx, byyx, bzyx, bxyy, byyy, bzyy, bxyz, byyz, bzyz, &
        bxzx, byzx, bzzx, bxzy, byzy, bzzy, bxzz, byzz, bzzz, bvec, &
        au_to_esu = 8.639418d-33
!-----------------------------------------------
!.....................................................................
!  CALCULATE BETA IN A NONITERATIVE FASHION
!..................................................................... 
!
      maxsq = norbs*norbs 
!
! READ DATA FROM ALPHA CALCULATION
!
      bavx = 0.0D+00 
      bavy = 0.0D+00 
      bavz = 0.0D+00 
!
      call daread (u0x, maxsq, 2) 
      call daread (u0y, maxsq, 3) 
      call daread (u0z, maxsq, 4) 
      call daread (g0x, maxsq, 5) 
      call daread (g0y, maxsq, 6) 
      call daread (g0z, maxsq, 7) 
      call daread (u1x, maxsq, 8) 
      call daread (u1y, maxsq, 9) 
      call daread (u1z, maxsq, 10) 
      call daread (g1x, maxsq, 11) 
      call daread (g1y, maxsq, 12) 
      call daread (g1z, maxsq, 13) 
! XXX
      call betall (u1x, g1x, u0x, g0x, u1x, g1x, nclose, norbs, bxxx) 
      bavx = bavx + 3.0D0*bxxx 
! YXX
      call betall (u1y, g1y, u0x, g0x, u1x, g1x, nclose, norbs, byxx) 
      bavy = bavy + byxx 
! ZXX
      call betall (u1z, g1z, u0x, g0x, u1x, g1x, nclose, norbs, bzxx) 
      bavz = bavz + bzxx 
! XXY
      call betall (u1x, g1x, u0x, g0x, u1y, g1y, nclose, norbs, bxxy) 
      bavy = bavy + bxxy 
! YXY
      call betall (u1y, g1y, u0x, g0x, u1y, g1y, nclose, norbs, byxy) 
      bavx = bavx + byxy 
! ZXY
      call betall (u1z, g1z, u0x, g0x, u1y, g1y, nclose, norbs, bzxy) 
! XXZ
      call betall (u1x, g1x, u0x, g0x, u1z, g1z, nclose, norbs, bxxz) 
      bavz = bavz + bxxz 
! YXZ
      call betall (u1y, g1y, u0x, g0x, u1z, g1z, nclose, norbs, byxz) 
! ZXZ
      call betall (u1z, g1z, u0x, g0x, u1z, g1z, nclose, norbs, bzxz) 
      bavx = bavx + bzxz 
! XYX
      call betall (u1x, g1x, u0y, g0y, u1x, g1x, nclose, norbs, bxyx) 
      bavy = bavy + bxyx 
! YYX
      call betall (u1y, g1y, u0y, g0y, u1x, g1x, nclose, norbs, byyx) 
      bavx = bavx + byyx 
! ZYX
      call betall (u1z, g1z, u0y, g0y, u1x, g1x, nclose, norbs, bzyx) 
! XYY
      call betall (u1x, g1x, u0y, g0y, u1y, g1y, nclose, norbs, bxyy) 
      bavx = bavx + bxyy 
! YYY
      call betall (u1y, g1y, u0y, g0y, u1y, g1y, nclose, norbs, byyy) 
      bavy = bavy + 3.0D0*byyy 
! ZYY
      call betall (u1z, g1z, u0y, g0y, u1y, g1y, nclose, norbs, bzyy) 
      bavz = bavz + bzyy 
! XYZ
      call betall (u1x, g1x, u0y, g0y, u1z, g1z, nclose, norbs, bxyz) 
! YYZ
      call betall (u1y, g1y, u0y, g0y, u1z, g1z, nclose, norbs, byyz) 
      bavz = bavz + byyz 
! ZYZ
      call betall (u1z, g1z, u0y, g0y, u1z, g1z, nclose, norbs, bzyz) 
      bavy = bavy + bzyz 
! XZX
      call betall (u1x, g1x, u0z, g0z, u1x, g1x, nclose, norbs, bxzx) 
      bavz = bavz + bxzx 
! YZX
      call betall (u1y, g1y, u0z, g0z, u1x, g1x, nclose, norbs, byzx) 
! ZZX
      call betall (u1z, g1z, u0z, g0z, u1x, g1x, nclose, norbs, bzzx) 
      bavx = bavx + bzzx 
! XZY
      call betall (u1x, g1x, u0z, g0z, u1y, g1y, nclose, norbs, bxzy) 
! YZY
      call betall (u1y, g1y, u0z, g0z, u1y, g1y, nclose, norbs, byzy) 
      bavz = bavz + byzy 
! ZZY
      call betall (u1z, g1z, u0z, g0z, u1y, g1y, nclose, norbs, bzzy) 
      bavy = bavy + bzzy 
! XZZ
      call betall (u1x, g1x, u0z, g0z, u1z, g1z, nclose, norbs, bxzz) 
      bavx = bavx + bxzz 
! YZZ
      call betall (u1y, g1y, u0z, g0z, u1z, g1z, nclose, norbs, byzz) 
      bavy = bavy + byzz 
! ZZZ
      call betall (u1z, g1z, u0z, g0z, u1z, g1z, nclose, norbs, bzzz) 
      bavz = bavz + 3.0D0*bzzz 
!
      bavx = bavx/5.0D+00 
      bavy = bavy/5.0D+00 
      bavz = bavz/5.0D+00 
!
      bvec = (bavx*bavx + bavy*bavy + bavz*bavz)**0.5D+00 
      write (iw, *) '  BETA (ELECTOPTIC POCKELS EFFECT) ' 
      write (iw, 10) bxxx, byxx, bzxx, bxxy, byxy, bzxy, bxxz, byxz, bzxz, bxyx&
        , byyx, bzyx, bxyy, byyy, bzyy, bxyz, byyz, bzyz, bxzx, byzx, bzzx, &
        bxzy, byzy, bzzy, bxzz, byzz, bzzz 
   10 format(/,/,'  BXXX  ',d15.8,'  BYXX ',d15.8,'  BZXX ',d15.8,/,'  BXXY  ',&
        d15.8,'  BYXY ',d15.8,'  BZXY ',d15.8,/,'  BXXZ  ',d15.8,'  BYXZ ',&
        d15.8,'  BZXZ ',d15.8,/,'  BXYX  ',d15.8,'  BYYX ',d15.8,'  BZYX ',&
        d15.8,/,'  BXYY  ',d15.8,'  BYYY ',d15.8,'  BZYY ',d15.8,/,'  BXYZ  ',&
        d15.8,'  BYYZ ',d15.8,'  BZYZ ',d15.8,/,'  BXZX  ',d15.8,'  BYZX ',&
        d15.8,'  BZZX ',d15.8,/,'  BXZY  ',d15.8,'  BYZY ',d15.8,'  BZZY ',&
        d15.8,/,'  BXZZ  ',d15.8,'  BYZZ ',d15.8,'  BZZZ ',d15.8) 
!
      if (bvec < 1.d10) then
        line = "(' AVERAGE BETA',a,'VALUE AT',f10.5,' EV = ',f15.4,a,g14.6,a)"
      else 
        line = "(' AVERAGE BETA',a,'VALUE AT',f10.5,' EV = ',g13.6,a,g14.6,a)"
      end if
      write (iw,"(/)")
      write (iw, trim(line)) "X      ", omega, bavx, " a.u. =", bavx*au_to_esu, " ESU"
      write (iw, trim(line)) "Y      ", omega, bavy, " a.u. =", bavy*au_to_esu, " ESU"
      write (iw, trim(line)) "Z      ", omega, bavz, " a.u. =", bavz*au_to_esu, " ESU"
      write (iw,"(2/)")
      write (iw,"(/)")


      return  
      end subroutine nonope 


      subroutine nonor(u0x, u1y, u1z, u1x, u0y, u0z, g0x, g1y, g1z, g1x, g0y, &
        g0z) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use molkst_C, only : norbs, nclose
      use chanel_C, only : iw
       USE polar_C, only : omega 
!
! THIS SUBROUTINE CALCULATES OPTICAL RECTIFICATION IN A
! NONITERATIVE WAY
!
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use daread_I 
      use betal1_I 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double)  :: u0x(norbs,norbs) 
      real(double)  :: u1y(norbs,norbs) 
      real(double)  :: u1z(norbs,norbs) 
      real(double)  :: u1x(norbs,norbs) 
      real(double)  :: u0y(norbs,norbs) 
      real(double)  :: u0z(norbs,norbs) 
      real(double)  :: g0x(norbs,norbs) 
      real(double)  :: g1y(norbs,norbs) 
      real(double)  :: g1z(norbs,norbs) 
      real(double)  :: g1x(norbs,norbs) 
      real(double)  :: g0y(norbs,norbs) 
      real(double)  :: g0z(norbs,norbs) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: maxsq 
      real(double) :: bavx, bavy, bavz, bxxx, byxx, bzxx, bxxy, byxy, bzxy, &
        bxxz, byxz, bzxz, bxyx, byyx, bzyx, bxyy, byyy, bzyy, bxyz, byyz, bzyz&
        , bxzx, byzx, bzzx, bxzy, byzy, bzzy, bxzz, byzz, bzzz, bvec 
!-----------------------------------------------
!.....................................................................
!  CALCULATE BETA IN A NONITERATIVE FASHION
!.....................................................................
      maxsq = norbs*norbs 
!
! READ DATA FROM ALPHA CALCULATION
!
      bavx = 0.0D+00 
      bavy = 0.0D+00 
      bavz = 0.0D+00 
!
      call daread (u0x, maxsq, 2) 
      call daread (u0y, maxsq, 3) 
      call daread (u0z, maxsq, 4) 
      call daread (g0x, maxsq, 5) 
      call daread (g0y, maxsq, 6) 
      call daread (g0z, maxsq, 7) 
      call daread (u1x, maxsq, 8) 
      call daread (u1y, maxsq, 9) 
      call daread (u1z, maxsq, 10) 
      call daread (g1x, maxsq, 11) 
      call daread (g1y, maxsq, 12) 
      call daread (g1z, maxsq, 13) 
!
! NONITERATIVE BETA CALCULATION
!
! XXX
      call betal1 (u0x, g0x, u1x, g1x, u1x, g1x, nclose, norbs, bxxx) 
      bavx = bavx + 3.0D0*bxxx 
! YXX
      call betal1 (u0y, g0y, u1x, g1x, u1x, g1x, nclose, norbs, byxx) 
      bavy = bavy + byxx 
! ZXX
      call betal1 (u0z, g0z, u1x, g1x, u1x, g1x, nclose, norbs, bzxx) 
      bavz = bavz + bzxx 
! XXY
      call betal1 (u0x, g0x, u1x, g1x, u1y, g1y, nclose, norbs, bxxy) 
      bavy = bavy + bxxy 
! YXY
      call betal1 (u0y, g0y, u1x, g1x, u1y, g1y, nclose, norbs, byxy) 
      bavx = bavx + byxy 
! ZXY
      call betal1 (u0z, g0z, u1x, g1x, u1y, g1y, nclose, norbs, bzxy) 
! XXZ
      call betal1 (u0x, g0x, u1x, g1x, u1z, g1z, nclose, norbs, bxxz) 
      bavz = bavz + bxxz 
! YXZ
      call betal1 (u0y, g0y, u1x, g1x, u1z, g1z, nclose, norbs, byxz) 
! ZXZ
      call betal1 (u0z, g0z, u1x, g1x, u1z, g1z, nclose, norbs, bzxz) 
      bavx = bavx + bzxz 
! XYX
      call betal1 (u0x, g0x, u1y, g1y, u1x, g1x, nclose, norbs, bxyx) 
      bavy = bavy + bxyx 
! YYX
      call betal1 (u0y, g0y, u1y, g1y, u1x, g1x, nclose, norbs, byyx) 
      bavx = bavx + byyx 
! ZYX
      call betal1 (u0z, g0z, u1y, g1y, u1x, g1x, nclose, norbs, bzyx) 
! XYY
      call betal1 (u0x, g0x, u1y, g1y, u1y, g1y, nclose, norbs, bxyy) 
      bavx = bavx + bxyy 
! YYY
      call betal1 (u0y, g0y, u1y, g1y, u1y, g1y, nclose, norbs, byyy) 
      bavy = bavy + 3.0D0*byyy 
! ZYY
      call betal1 (u0z, g0z, u1y, g1y, u1y, g1y, nclose, norbs, bzyy) 
      bavz = bavz + bzyy 
! XYZ
      call betal1 (u0x, g0x, u1y, g1y, u1z, g1z, nclose, norbs, bxyz) 
! YYZ
      call betal1 (u0y, g0y, u1y, g1y, u1z, g1z, nclose, norbs, byyz) 
      bavz = bavz + byyz 
! ZYZ
      call betal1 (u0z, g0z, u1y, g1y, u1z, g1z, nclose, norbs, bzyz) 
      bavy = bavy + bzyz 
! XZX
      call betal1 (u0x, g0x, u1z, g1z, u1x, g1x, nclose, norbs, bxzx) 
      bavz = bavz + bxzx 
! YZX
      call betal1 (u0y, g0y, u1z, g1z, u1x, g1x, nclose, norbs, byzx) 
! ZZX
      call betal1 (u0z, g0z, u1z, g1z, u1x, g1x, nclose, norbs, bzzx) 
      bavx = bavx + bzzx 
! XZY
      call betal1 (u0x, g0x, u1z, g1z, u1y, g1y, nclose, norbs, bxzy) 
! YZY
      call betal1 (u0y, g0y, u1z, g1z, u1y, g1y, nclose, norbs, byzy) 
      bavz = bavz + byzy 
! ZZY
      call betal1 (u0z, g0z, u1z, g1z, u1y, g1y, nclose, norbs, bzzy) 
      bavy = bavy + bzzy 
! XZZ
      call betal1 (u0x, g0x, u1z, g1z, u1z, g1z, nclose, norbs, bxzz) 
      bavx = bavx + bxzz 
! YZZ
      call betal1 (u0y, g0y, u1z, g1z, u1z, g1z, nclose, norbs, byzz) 
      bavy = bavy + byzz 
! ZZZ
      call betal1 (u0z, g0z, u1z, g1z, u1z, g1z, nclose, norbs, bzzz) 
      bavz = bavz + 3.0D0*bzzz 
!
      bavx = bavx/5.0D+00 
      bavy = bavy/5.0D+00 
      bavz = bavz/5.0D+00 
!
      bvec = (bavx*bavx + bavy*bavy + bavz*bavz)**0.5D+00 
      write (iw, 10) 
   10 format(/,/,' BETA (OPTICAL RECTIFICATION) ') 
      write (iw, 20) bxxx, byxx, bzxx, bxxy, byxy, bzxy, bxxz, byxz, bzxz, bxyx&
        , byyx, bzyx, bxyy, byyy, bzyy, bxyz, byyz, bzyz, bxzx, byzx, bzzx, &
        bxzy, byzy, bzzy, bxzz, byzz, bzzz 
   20 format(/,/,'  BXXX  ',d15.8,'  BYXX ',d15.8,'  BZXX ',d15.8,/,'  BXXY  ',&
        d15.8,'  BYXY ',d15.8,'  BZXY ',d15.8,/,'  BXXZ  ',d15.8,'  BYXZ ',&
        d15.8,'  BZXZ ',d15.8,/,'  BXYX  ',d15.8,'  BYYX ',d15.8,'  BZYX ',&
        d15.8,/,'  BXYY  ',d15.8,'  BYYY ',d15.8,'  BZYY ',d15.8,/,'  BXYZ  ',&
        d15.8,'  BYYZ ',d15.8,'  BZYZ ',d15.8,/,'  BXZX  ',d15.8,'  BYZX ',&
        d15.8,'  BZZX ',d15.8,/,'  BXZY  ',d15.8,'  BYZY ',d15.8,'  BZZY ',&
        d15.8,/,'  BXZZ  ',d15.8,'  BYZZ ',d15.8,'  BZZZ ',d15.8) 
!
      write (iw, 30) omega, bavx, " a.u."
   30 format(/,/,' AVERAGE BETAX VALUE AT ',f10.5,'EV = ',1f15.5, a) 
!
      write (iw, 40) omega, bavy, " a.u."
   40 format(' AVERAGE BETAY VALUE AT ',f10.5,'EV = ',1f15.5, a) 
!
      write (iw, 50) omega, bavz, " a.u."
   50 format(' AVERAGE BETAZ VALUE AT ',f10.5,'EV = ',1f15.5,a,/,/) 
!
      write (iw, 60) omega, bvec, " a.u."
   60 format(/,/,' AVERAGE BETA(OR) VALUE AT ',f10.5,'EV = ',1f15.5,a,/,/) 
!
      return  
      end subroutine nonor 


      subroutine openda(irest) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      use chanel_C, only : idaf, irecln, irecst, ioda, ifilen, &
      pol_fn
!
!     - - - - OPEN MASTER DICTIONARY FILE 10 - - - -
!
!...Translated by Pacific-Sierra Research 77to90  4.4G  21:23:38  03/14/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: irest 
!-----------------------------------------------
!
!
      idaf = 17 
      irecln = 1023 
      open(unit=idaf, file=pol_fn, status='UNKNOWN', access=&
        'DIRECT', form='UNFORMATTED', recl=8*irecln) 
!
!     ----- IS THIS A NEW OR OLD DAF FILE -----
!
      if (irest == 0) then 
!
!        ----- MARK THE NEW DAF RECORDS AS EMPTY -----
!
        irecst = 1 
        ioda = -1 
        irecst = irecst + 1 
        write (unit=idaf, rec=1) irecst, ioda, ifilen 
        return  
      endif 
!
!     ----- LOAD THE OLD DAF DIRECTORY -----
!
      read (unit=idaf, rec=1) irecst, ioda, ifilen 
      return  
      end subroutine openda 
  double precision function pol_vol(average)
  use molkst_C, only: numat, method_MNDO, method_AM1, method_PM3, method_PM6, method_PM7
  use parameters_C, only: polvol
  use common_arrays_C, only: labels
  use funcon_C, only: a0
  implicit none
  !
  !.. Formal Arguments ..
  double precision, intent (in) :: average
  double precision :: polarizability
  integer :: i
  double precision :: C_1, C_2
!
!  The corrected polarizability (alpha) is given by:
!
!          P = Sum_i C_i  +  P'*C_1 + C_2  
!
!   where C_i are constants per atom, C_1 and C_2 are constants.
!
!  If these constants need to be altered, use v_par in PARAM to get new values.
!
     if (method_MNDO) then
      C_1 = 0.5651170d0
      C_2 = 0.2859620d0
    else if (method_AM1) then
      C_1 = 0.54687d0
      C_2 = 0.321973d0
    else if (method_PM3) then
      C_1 = 0.2465020d0
      C_2 = 0.0d0
     else if (method_PM6) then
      C_1 = 0.791396d0
      C_2 = 0.373638d0
    else if (method_PM7) then
      C_1 =  0.445109d0
      C_2 =  0.351109d0
    else
      C_1 = 0.d0
      C_2 = 0.d0
    end if
    polarizability = average*a0**3
    if (C_1 > 1.d-4) then
      polarizability = polarizability * C_1 + C_2
      do i = 1, numat
        polarizability = polarizability + polvol(labels(i))
      end do
    end if
    pol_vol = polarizability/a0**3
    return
  end function pol_vol
