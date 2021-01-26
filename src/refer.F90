      subroutine refer 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE journal_references_C, only : allref
      use chanel_C, only : iw
      use common_arrays_C, only : nat
      use molkst_C, only : numat, keywrd, is_PARAM, method_rm1, &
      method_mndo, method_am1, method_pm3, method_pm6, method_mndod, &
      method_pm7, sparkle, method_pm8, method_indo
      use parameters_C, only : gss
      use parameters_for_PM6_Sparkles_C, only : gss6sp
      use parameters_for_AM1_Sparkles_C, only : gssam1sp
      use parameters_for_PM3_Sparkles_C, only : gssPM3sp
      use reimers_C, only : isok
!***********************************************************************
!DECK MOPAC
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use mopend_I  
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, mode
      logical :: allok 
      logical , dimension(107) :: elemns 
      logical :: mixok, mix, exter

      save mix 
!-----------------------------------------------
      data mix/ .FALSE./  
      elemns(:102) = .FALSE.
      elemns(nat(:numat)) = .TRUE. 
      if (method_indo) then
        allok = .true.
        do i = 1,numat
          if (.not. isok(nat(i))) allok = .false.
        end do
        if (.not. allok) then
          write (iw, 40) 'SOME ELEMENTS HAVE BEEN SPECIFIED FOR WHICH', &
            'NO PARAMETERS ARE AVAILABLE.  CALCULATION STOPPED.'
          call mopend("Parameters for some elements are missing")
        end if
        return
      end if
      if (method_pm6) then
        allok = .true.
        if (sparkle) then !  Check for Lanthanides
          write(iw,*)
          write (iw, '(/,a)') " References for PM6:"   
          do i = 57, 71
            if (elemns(i)) then
              write (iw, '(A)') allref(i,2)(:len_trim(allref(58,2)))
            end if
          end do
        else
          do i = 58, 70
            if (elemns(i)) then
              write (iw, '(A,I3)') ' DATA ARE NOT AVAILABLE FOR ELEMENT NO.', i 
              if (Abs(gss6sp(i))   > 0.1d0) write(iw,*)" (Parameters are available if PM6 SPARKLE is used)"
              allok = .FALSE. 
            end if
          end do
        end if
        write (iw, '(/,a)') " General Reference for PM6:"        
        write (iw, '(1x,a)') &
          &  """Optimization of Parameters for Semiempirical Methods V: Modification of NDDO Approximations", &
          &  "and Application to 70 Elements"", J. J. P. Stewart, J. Mol. Mod., 13, 1173-1213 (2007)", &
          &  "URL: https://link.springer.com/article/10.1007/s00894-007-0233-4"
        if (index(keywrd," PM6-DH+") /= 0) then
          write (iw, '(/,a)') " Reference for PM6-DH+:"  
          write(iw,'(1x,a)') &
            & """Third-Generation Hydrogen-Bonding Corrections for Semiempirical QM Methods and Force Fields""", &
            & "Martin Korth, J. Chem. Theory Comput., 6 (12), pp 3808-3816 (2010)", &
            & "URL: http://pubs.acs.org/doi/abs/10.1021/ct100408b"
        end if
        if (allok .or. is_PARAM) return  
        write (iw, 40) 'SOME ELEMENTS HAVE BEEN SPECIFIED FOR WHICH', &
          'NO PARAMETERS ARE AVAILABLE.  CALCULATION STOPPED.' 
        call mopend("Parameters for some elements are missing")
        return  
      end if
      if (method_pm7) then
        allok = .true.
        if (sparkle) then !  Check for Lanthanides
          write(iw,*)
          write (iw, '(/,a)') " References for PM7 Sparkles:"   
          do i = 57, 71
            if (elemns(i)) then
              write (iw, '(A)') " ""Sparkle/PM7 Lanthanide Parameters for the Modeling of Complexes and Materials""", &
                & " J. D. L. Dutra, M. A. M. Filho1, R. O. Freire, G. B. Rocha, A. M. Simas, and J. J. P. Stewart", &
                & " Chemistry of Materials (submitted)"
              exit
            end if
          end do
        else
          do i = 58, 70
            if (elemns(i)) then
              write (iw, '(A,I3)') ' DATA ARE NOT AVAILABLE FOR ELEMENT NO.', i 
              if (Abs(gss6sp(i))   > 0.1d0) write(iw,*)" (Parameters are available if SPARKLE is used)"
              allok = .FALSE. 
            end if
          end do
        end if
        write (iw, '(/,a)') " General Reference for PM7:"        
        write (iw, '(1x,a)') &
          &  """Optimization of Parameters for Semiempirical Methods VI: More Modifications to the ", &
          & "NDDO Approximations and Re-optimization of Parameters"", J. J. P. Stewart, J. Mol. Mod., 1:32, 19 (2013)", &
          & "https://link.springer.com/article/10.1007/s00894-012-1667-x"
        if (allok .or. is_PARAM) return  
        write (iw, 40) 'SOME ELEMENTS HAVE BEEN SPECIFIED FOR WHICH', &
          'NO PARAMETERS ARE AVAILABLE.  CALCULATION STOPPED.' 
        call mopend("Parameters for some elements are missing")
        return  
      end if
      mixok = index(keywrd,'PARASOK') /= 0 
      exter = index(keywrd," EXTERNAL") /= 0
       
      if (method_pm7)   mode = 7   
      if (method_rm1)   then
        if (index(keywrd, "SPARKLES") /= 0) then
          do i = 58, 71
            allref(i,6) = allref(i,8)
          end do
        end if
        mode = 6
      end if
      if (method_mndod) mode = 5 
      if (method_pm3)   mode = 4 
      if (method_am1)   mode = 3 
      if (method_pm6)   mode = 2 
      if (method_mndo)  mode = 1  
      if (method_pm8)   mode = 9  
      allref(99,mode) = &
        ' DUMMY ATOMS ARE USED; THESE DO NOT AFFECT THE CALCULATION' 
      allref(100,mode) = ' '      
      allok = .TRUE. 
       write (iw, '(A)')
      do i = 1, 102 
        if (.not.elemns(i)) cycle  
        if (i < 99 .and. .not.mix .and. mode==3) mix = index(allref(i,3),'MNDO') /= 0 
        if (allref(i,5)(1:1) /= ' ' .or. allref(i,5)(1:4) == "    " ) allref(i,5) = allref(i,1) 
        if (allref(i,mode)(1:1) /= ' ' .or. &
          allref(i,mode)(1:4) == "    " .or. Abs(gss(i)) < 0.1d0 ) then 
          if (.not. exter) then
          write (iw, '(A,I3)') ' DATA ARE NOT AVAILABLE FOR ELEMENT NO.', i 
          if (Abs(gssam1sp(i)) > 0.1d0) write(iw,*)" (Parameters are available if AM1 SPARKLE is used)"
          if (Abs(gssPM3sp(i)) > 0.1d0) write(iw,*)" (Parameters are available if PM3 SPARKLE is used)"
          if (Abs(gss6sp(i))   > 0.1d0) write(iw,*)" (Parameters are available if PM6 SPARKLE is used)"
          allok = .FALSE. 
          end if
        else 
          write (iw, '(A)') allref(i,mode)(:len_trim(allref(i,mode)))
          if (mode == 7) exit
        endif 
      end do 
      if (mix .and. .not.mixok) then 
        write (iw, 40) 'SOME ELEMENTS HAVE BEEN SPECIFIED FOR WHICH ONLY MNDO'&
          , 'PARAMETERS ARE AVAILABLE.  SUCH MIXTURES OF METHODS ARE', &
          'VERY RISKY AND HAVE NOT BEEN FULLY TESTED.  IF YOU FEEL', &
          'THE RISK IS WORTH WHILE - CHECK THE MANUAL FIRST - THEN', &
          'SPECIFY "PARASOK" IN THE KEYWORDS' 
          call mopend ('MIXED PARAMETER SETS.  USE "PARASOK" TO CONTINUE') 
        return  
      endif 
      if (allok .or. is_PARAM .or. index(keywrd, "0SCF") /= 0) return  
      write (iw, 40) 'SOME ELEMENTS HAVE BEEN SPECIFIED FOR WHICH', &
        'NO PARAMETERS ARE AVAILABLE.  CALCULATION STOPPED.' 
        call mopend("Parameters for some elements are missing")
      return  
   40 format(/,/,/,/,/,10x,a,4(/,10x,a))  
      end subroutine refer 
