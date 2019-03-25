      subroutine upcase(keywrd, n) 
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:36:11  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      character , intent(inout) :: keywrd*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icapa, ilowa, ilowz, i, iline, j 
      character :: keybuf*2000 
!-----------------------------------------------
!
!  UPCASE WILL TAKE A CHARACTER STRING, IN KEYWRD, AND PUT IT INTO
!  UPPER CASE.  KEYWRD IS LIMITED TO 80 CHARACTERS
!
      icapa = ichar('A') 
      ilowa = ichar('a') 
      ilowz = ichar('z') 
      keybuf = keywrd 
      do i = 1, n 
        iline = ichar(keywrd(i:i)) 
        if (iline>=ilowa .and. iline<=ilowz) &
        keywrd(i:i) = char(iline + icapa - ilowa) 
        if (iline /= 9) cycle  
!
!  Change tabs to spaces.  A tab is ASCII character 9.
!
        keywrd(i:i) = ' ' 
      end do 
!
!   If the word EXTERNAL is present, do NOT change case of the following
!   character string.
!
      i = index(keywrd,'EXTERNAL=') 
      if (i /= 0) then 
        j = index(keywrd(i+1:),' ') + i 
        keywrd(i+9:j) = keybuf(i+9:j) 
      endif 
      return  
      end subroutine upcase 
