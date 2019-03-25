      subroutine symr 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use molkst_C, only : numat
      use common_arrays_C, only : coord 
      use chanel_C, only : iw
      use symmetry_C, only : nclass, r, ipo, elem, nsym, nent
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  11:05:03  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use mopend_I 
      use symp_I 
      implicit none
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: maxent = 6 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nvalue, i, ndflt, k, j, n 
      real(double) :: x, y, z, xa, ya, za, dist 
      logical :: prob 
!-----------------------------------------------
!
!     R(9,*):   The 9 elements of each record are a packed 3 by 3
!          array of a given symmetry operations.
!    IPO(n,*):  A vector that contains the symmetry mapping of atomic ce
!
      if (allocated(ipo)) deallocate(ipo)
      allocate(ipo(numat,120))
      prob = .FALSE. 
      nvalue = 0 
!  Get the symmetry functions: (NOTE: THE FIRST IS ALWAYS E)
      r(1,1) = 1.D0 
      r(2,1) = 0.D0 
      r(3,1) = 0.D0 
      r(4,1) = 0.D0 
      r(5,1) = 1.D0 
      r(6,1) = 0.D0 
      r(7,1) = 0.D0 
      r(8,1) = 0.D0 
      r(9,1) = 1.D0 
!
!  CENTER THE MOLECULE
!
      x = 0.D0 
      y = 0.D0 
      z = 0.D0 
      do i = 1, numat 
        x = x + coord(1,i) 
        y = y + coord(2,i) 
        z = z + coord(3,i) 
        ipo(i,1) = i 
      end do 
      xa = x/dble(numat) 
      ya = y/dble(numat) 
      za = z/dble(numat) 
      coord(1,:numat) = (-xa) + coord(1,:numat) 
      coord(2,:numat) = (-ya) + coord(2,:numat) 
      coord(3,:numat) = (-za) + coord(3,:numat) 
!
      nent = 1 
      nsym = 0 
      ndflt = 1 
   30 continue 
      nsym = nsym + 1 
      ndflt = ndflt + 1 
      if (ndflt <= nclass) then 
!
!   Copy Symmetry Operation from ELEM
!
        k = 0 
        nvalue = 1 
        do i = 1, 3 
          r(k+1:3+k,1+nent) = elem(i,:,ndflt) 
          k = 3 + k 
        end do 
!  NOW, TO CALCULATE THE IPO OF THIS FUNCTION
        nent = 1 + nent 
        n = nent 
!  Now, to initialize IPO(n) and
!  Perform R on each atomic center and determine where it maps to.
        l80: do i = 1, numat 
          x = coord(1,i)*r(1,n) + coord(2,i)*r(2,n) + coord(3,i)*r(3,n) 
          y = coord(1,i)*r(4,n) + coord(2,i)*r(5,n) + coord(3,i)*r(6,n) 
          z = coord(1,i)*r(7,n) + coord(2,i)*r(8,n) + coord(3,i)*r(9,n) 
          ipo(i,n) = 0 
          do j = 1, numat 
            dist = abs(x - coord(1,j)) + abs(y - coord(2,j)) + abs(z - coord(3,&
              j)) 
            if (dist >= 0.6D0) cycle  
            if (ipo(i,n) == 0) then 
              ipo(i,n) = j 
            else 
              write (iw, 50) 
              prob = .TRUE. 
              exit  l80 
   50         format('  ONE ATOM MAPS ONTO TWO DIFFERENT ATOMIC C','ENTERS',/,&
                '  ADD KEYWORD '' NOSYM'' AND RE-RUN') 
            endif 
          end do 
          if (ipo(i,n) /= 0) cycle  l80 
          write (iw, 70) 
   70     format('  ONE ATOM MAPS ONTO NO OTHER ATOM ',/,&
            '  ADD KEYWORD '' NOSYM'' AND RE-RUN') 
          prob = .TRUE. 
          exit  l80 
        end do l80 
!
!
      endif 
      if (nvalue/=0 .and. nsym<maxent) go to 30 
!
!  If a problem exists.  Stop the program.
!
      if (prob) then 
        write (iw, *) ' PROBLEM IN SYMR' 
        call mopend ('PROBLEM IN SYMR') 
        return  
      endif 
      nsym = nent 
!
!  NEXT, EXPAND THE EXISTING OPERATORS TO THE FULL SET
!
      call symp 
!
      return  
      end subroutine symr 
      subroutine symp 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use symmetry_C, only : nsym, r, ipo, name
      USE molkst_C, ONLY: numat 
      USE chanel_C, only : iw
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  15:34:19  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: maxfun = 120 
      real(double), parameter :: tol = 1D-2 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, n, m 
      real(double) :: res 
!-----------------------------------------------
!****************************************************************
!
!   ON INPUT   R    = SYMMETRY OPERATIONS (7 MAX)
!              IPO  = PERM OPR FOR ABOVE OPERATIONS
!              NSYM = CURRENT NUMBER OF SYMMETRY OPERATIONS
!              NENT = NUMBER OF USER SUPPLIED OPERATIONS
!
!   ON OUTPUT  R    = SYMMETRY OPERATIONS (120 MAX)
!              IPO  = PERMUTATION OPERATOR FOR SYMMETRY OPERATIONS
!              NSYM = NUMBER OF SYMMETRY OPERATIONS
!
!****************************************************************
!
!  A SUBROUTINE THAT WILL EXPAND THE SYMMETRY OPERATIONS READ IN INTO
!     THE COMPLETE SET.  NOTE: VERY FEW OPERATIONS ARE REQUIRED TO
!     GENERATE EVEN VERY LARGE GROUPS OF OPERATIONS.
!
!
!
!
!  Variables used:  (n represents the number of atomic centers)
!
!    For the next two items, the last index represents the symmetry
!        operation number.
!     R(9,*):   The 9 elements of each record are a packed 3 by 3
!          array of a given symmetry operation.
!    IPO(n,*):  A vector that contains the symmetry mapping of atomic ce
!
!  NSYM IS ALWAYS THE UPPER BOUND OF THE VALID FUNCTIONS.  QUIT IF IT
!     REACHES 120.
!  I IS THE SLOW INDEX OF FUNCTIONS TO MULTIPLY
!  J IS THE FAST INDEX OF FUNCTIONS TO MULTIPLY
!  ALWAYS DO R(I)*R(J) AND TAKE I,J FROM 2 TO NSYM
!
      i = 2 
      j = 1 
!
!  DETERMINE IF IT IS TIME TO STOP
!
   10 continue 
      j = j + 1 
      if (j > nsym) then 
        j = 2 
        i = i + 1 
        if (i > nsym) go to 50 
      endif 
      if (nsym == maxfun) go to 50 
!
!  NOW TO START THE MULTIPLICATION
!
      r(1,nsym+1) = r(1,i)*r(1,j) + r(2,i)*r(4,j) + r(3,i)*r(7,j) 
      r(2,nsym+1) = r(1,i)*r(2,j) + r(2,i)*r(5,j) + r(3,i)*r(8,j) 
      r(3,nsym+1) = r(1,i)*r(3,j) + r(2,i)*r(6,j) + r(3,i)*r(9,j) 
      r(4,nsym+1) = r(4,i)*r(1,j) + r(5,i)*r(4,j) + r(6,i)*r(7,j) 
      r(5,nsym+1) = r(4,i)*r(2,j) + r(5,i)*r(5,j) + r(6,i)*r(8,j) 
      r(6,nsym+1) = r(4,i)*r(3,j) + r(5,i)*r(6,j) + r(6,i)*r(9,j) 
      r(7,nsym+1) = r(7,i)*r(1,j) + r(8,i)*r(4,j) + r(9,i)*r(7,j) 
      r(8,nsym+1) = r(7,i)*r(2,j) + r(8,i)*r(5,j) + r(9,i)*r(8,j) 
      r(9,nsym+1) = r(7,i)*r(3,j) + r(8,i)*r(6,j) + r(9,i)*r(9,j) 
!
!  IS IT UNIQUE?
!
      do n = 1, nsym 
        res = 0.D0 
        do m = 1, 9 
          res = res + abs(r(m,n)-r(m,nsym+1)) 
        end do 
        if (res < tol) go to 10 
      end do 
!
!  YES, IT IS UNIQUE.  NOW, GENERATE THE NEW IPO(,NSYM)
!
      nsym = nsym + 1 
      do n = 1, numat 
        ipo(n,nsym) = ipo(ipo(n,j),i) 
      end do 
!
!     ALL DONE ADDING THE NEW FUNCTION.  GO TRY TO FIND A NEW ONE.
!
      go to 10 
!
!
   50 continue 
      write (iw, 60) name, nsym 
   60 format(/,'    FOR POINT-GROUP ',a4,' THERE ARE ',i3,&
        ' UNIQUE SYMMETRY FUNCTIONS.',/) 
!#       WRITE(IW,*)' Symmetry Operations'
!#       DO 14 I=1,NSYM
!#       WRITE(IW,*)
!#       WRITE(IW,'(3f12.6)')(R(J,I),J=1,9)
!#  14   continue
!
      return  
      end subroutine symp 

