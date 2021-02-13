      subroutine nuchar(line, l_line, value, nvalue) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: l_line
      integer , intent(out) :: nvalue 
      character  :: line*(*)
      double precision , intent(out) :: value(40) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(40) :: istart 
      integer :: i 
      logical :: leadsp 
      character :: tab, comma, space 
      double precision, external :: reada

      save comma, space 
!-----------------------------------------------
!***********************************************************************
!
!   NUCHAR  DETERMINS AND RETURNS THE REAL VALUES OF ALL NUMBERS
!           FOUND IN 'LINE'. ALL CONNECTED SUBSTRINGS ARE ASSUMED
!           TO CONTAIN NUMBERS
!   ON ENTRY LINE    = CHARACTER STRING
!   ON EXIT  VALUE   = ARRAY OF NVALUE REAL VALUES
!
!***********************************************************************
      data comma, space/ ',', ' '/  
      tab = char(9) 
!
! CLEAN OUT TABS AND COMMAS
!
      do i = 1, l_line 
        if (line(i:i)/=tab .and. line(i:i)/=comma) cycle  
        line(i:i) = space 
      end do 
!
! FIND INITIAL DIGIT OF ALL NUMBERS, CHECK FOR LEADING SPACES FOLLOWED
!     BY A CHARACTER
!
      leadsp = .TRUE. 
      nvalue = 0 
      do i = 1, l_line 
        if (leadsp .and. line(i:i)/=space) then 
          nvalue = nvalue + 1 
          istart(nvalue) = i 
        end if 
        leadsp = line(i:i) == space 
      end do 
!
! FILL NUMBER ARRAY
!
      do i = 1, nvalue 
        value(i) = reada(line,istart(i)) 
      end do 
      return  
      end subroutine nuchar 
