      logical function myword (keywrd, testwd) 
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:35:31  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character , intent(inout) :: keywrd*(*) 
      character , intent(in) :: testwd*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, k, len_key
      logical :: quote
!-----------------------------------------------
      myword = .FALSE. 
      len_key = len_trim(keywrd)
!
!  If the keyword containes quoted text, then ignore text between quotation marks.
!  quote starts off .FALSE.  
!  quote is turned .TRUE. if a quotation mark, '"', is found
!  quote is turned .FALSE. when another quotation mark is found.
      quote = .false.
   10 continue 
      j = index(keywrd,testwd) 
      if (j /= 0) then 
!
!  Keyword found.  Now delete keyword in its entirety
!
   20   continue 
        do while(keywrd(j:j) == ' ') 
          j = j + 1 
        end do 
        myword = .TRUE. 
        do k = j, len_key 
          if (keywrd(k:k) == '"') quote = (.not. quote)
          if (.not. quote) then
            if (keywrd(k:k)=='=' .or. keywrd(k:k)==' ') then 
!
!     CHECK FOR ATTACHED '=' SIGN
!
              j = k 
              if (keywrd(j:j) == '=') go to 50 
!
!     CHECK FOR SEPARATED '=' SIGN
!
              do j = k + 1, len_key
                if (keywrd(j:j) == '=') go to 50            
                if (keywrd(j:j) /= ' ') go to 10 
              end do 
!
!    THERE IS NO '=' SIGN ASSOCIATED WITH THIS KEYWORD
!
              go to 10 
     50       continue 
              keywrd(j:j) = ' ' 
!
!   THERE MUST BE A NUMBER AFTER THE '=' SIGN, SOMEWHERE
!
              go to 20 
            end if
          end if
          keywrd(k:k) = ' ' 
        end do 
      endif 
      return  
      end function myword 
