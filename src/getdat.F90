      subroutine getdat(input, output)  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      use molkst_C, only : natoms, jobnam, run, ijulian, verson, &
      gui, line, ncomments, is_PARAM, keywrd, arc_hof_1, arc_hof_2
      use chanel_C, only : iw0, job_fn, input_fn, iw
      use common_arrays_C, only : all_comments
!...Translated by Pacific-Sierra Research 77to90  4.4G  22:53:42  03/15/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!----------------------------------------------- 
      use mopend_I 
      use to_screen_I
      implicit none
!-----------------------------------------------
      integer, intent (in) :: input, output
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, parameter :: from_data_set = 7
      integer :: i, j, io_stat, l, nlines
      logical :: exists, arc_file, comments = .true.
      character :: text*90, line1*1000, num1*1, num2*1
      character, allocatable :: tmp_comments(:)*120
      double precision, external :: reada
      save i 
!-----------------------------------------------
!
!***********************************************************************
!
!   GETDAT READS IN ALL THE DATA USING CHANEL "from_data_set", AND WRITES IT
!   TO SCRATCH CHANNEL "input".  THIS WAY THE ORIGINAL DATA-SET IS
!   FREED UP AS SOON AS THE JOB STARTS.
!******************************************************************** 
      natoms = 0 
      call to_screen("To_file: getdat")
      if (gui) then
        jobnam = "MOPAC input"
        natoms = 1
      else
        if (run /= 2 .or.jobnam ==" ") then
          i = iargc()
          if (i >= run) then
            call getarg (run, jobnam)
            natoms = 1
            do i = len_trim(jobnam), 1, -1   !  Remove any unprintable characters from the end of the file-name
              if (ichar(jobnam(i:i)) > 39 .and. ichar(jobnam(i:i)) < 126 .or. jobnam(i:i) =="'") exit
            end do
            jobnam(i + 1:) = " "
          else if (i == 0) then
            if (is_PARAM) then
              write(line,'(2a)')" PARAM is the parameter optimization program for use with MOPAC"
              write(0,'(//10x,a,/)')trim(line)
              write(0,'(10x,a)')" It uses a single argument, the PARAM data-set"
              write(0,'(10x,a)')" The command to run PARAM is 'PARAM.exe <data-set>.dat'"
             ! call web_message(0,"running_MOPAC.html")
              write(0,'(10x,a)')" Press (return) to continue"
              read(5,*, iostat=i) 
              return
            else
              write(0,'(//10x,a,/)')" MOPAC is a semiempirical quantum chemistry program"
              write(0,'(10x,a)')" It uses a single argument, the MOPAC data-set"
              write(0,'(10x,a)')" The command to run MOPAC is 'MOPAC2016.exe <data-set>.mop'"
              call web_message(0,"running_MOPAC.html")
              write(0,'(10x,a)')" Press (return) to continue"
              read(5,*, iostat=i) 
              return
            end if            
          end if
        else
          natoms = 1
        end if
      end if
      if (natoms == 0) return
!
! Check for the data set in the order: <file>.mop, <file>.dat, <file>
!
      line = jobnam
      exists = .false.
      call upcase(line, len_trim(line))
      i = Index(line,".MOP") + Index(line,".DAT") + Index(line, ".ARC")
      arc_file = (Index(line, ".ARC") > 0)
      if (i > 0) then  ! User has supplied a suffix - use it     
        line1 = trim(line)
        line = trim(jobnam)
        inquire(file=line, exist=exists)
        j = len_trim(line1)
        if (line1(j - 3:j) == ".TXT") then
          jobnam(j - 3:j) = " "
          if (j - i > 7) then
          line1 = "There must not be any text between """//line1(i:i+3)// &
          """ and "".TXT"" in the file-name"
           open(unit=iw, file=trim(jobnam)//'.out') 
           call mopend(trim(line1))
           write(iw,'(/10x,a)')"(End of file name: """//trim(line(i:))//""")"           
           return
          end if          
        end if
        jobnam(i:i+3) = " "
      else  ! No suffix supplied, try the file, then .mop, then .dat
        line = trim(jobnam)
        inquire(file=line, exist=exists)
        if (exists) then ! Check that it is not a folder
          open(unit=from_data_set, file=line, status='OLD', position=&
          'asis', iostat=io_stat) 
          if (io_stat == 9) then
            exists = .false.
          else
            close (from_data_set)
          end if
        end if
        if ( .not. exists) then
          line = trim(jobnam)//'.mop'
          inquire(file=line, exist=exists)
        end if
        if ( .not. exists) then
          line = trim(jobnam)//'.dat'
          inquire(file=line, exist=exists)
        end if
      end if
!
!  Check that a password has not been used
!
      j = 0
      do i = len_trim(jobnam), 1, -1
        if (ichar(jobnam(i:i)) > ichar("9") .or. ichar(jobnam(i:i)) < ichar("0")) exit
        j = j + 1
      end do
      if (j == 7 .or. j == 8) then
        if(jobnam(i:i) == "a" .or. jobnam(i:i) == "A") then
          j = 0
          do i = i - 1, 1, -1
            if (ichar(jobnam(i:i)) > ichar("9") .or. ichar(jobnam(i:i)) < ichar("0")) exit
            j = j + 1
          end do
          if (j == 7 .or. j == 8) then
            write (line, '(a)') ' Password detected, but password has already been correctly installed.'
            open(unit=iw, file=trim(jobnam)//'.out') 
            write(iw,"(//10x,a)")trim(line)
            write(*,"(//,a)")trim(line)
            write (line, '(a)') ' MOPAC is ready to run data sets.' 
            write(*,"(/,a,/)")trim(line)
            call mopend (trim(line)) 
            return  
          end if
        end if
      end if
  98  if (exists) then 
         if (iw0 > -1) then           
           if (ijulian < 365) then
             call to_screen("********************************************************************************")
             write (text,"(a,a,a,i4,a)") &
!            " <-Start of text                                                              End of text->"
             "**                          MOPAC (",verson,")        Days remaining",ijulian,"         **"
             call to_screen(text)
             call to_screen("********************************************************************************")
           end if
           call to_screen("Preparing to read the following MOPAC file: ")
           i = min(len_trim(line), 240)
                         call to_screen(line(:min(i,120)))
           if (i > 120)  call to_screen(line(121:min(i,240)))
         end if
         job_fn = line(:len(job_fn))
        open(unit=from_data_set, file=job_fn, status='OLD', position='asis', iostat=io_stat) 
        if (io_stat /= 0) then
          write(line,'(2a)')" Data file: '"//trim(job_fn)//"' exists, but it cannot be opened."
          write(0,'(//10x,a,//)')trim(line)
          call to_screen(" File '"//job_fn(:len_trim(job_fn))//"' cannot be opened")
          call mopend("File '"//job_fn(:len_trim(job_fn))//"' cannot be opened")
          return
        end if
!
!  Now that the name of the data-set is known, set up all the other file-names
!
        call init_filenames 
      endif 
!
!  CLOSE UNIT IFILES(5) IN CASE IT WAS ALREADY PRE-ASSIGNED.
!INPUT FILE MISSING 
      close(input) 
      open(unit=input, status='SCRATCH', iostat = io_stat) 
      if (io_stat /= 0) then
        if (io_stat == 30) then
          call to_screen(" The file'"//input_fn(:len_trim(input_fn))//"' is busy")
          call to_screen(" Correct this problem, and re-submit job")
          call mopend("Temporary file '"//input_fn(:len_trim(input_fn))//"' is unavailable for use")
          return
        end if
      end if
      rewind input 
      rewind from_data_set 
      arc_hof_1 = 0.d0
      arc_hof_2 = 0.d0
      if (arc_file) then
        do 
          read (from_data_set, '(A1000)', iostat = i) line
          if (i /= 0) then
            i = 0
            ncomments = 0
            rewind(from_data_set)
            exit
          end if 
          if (index(line, "HEAT OF FORMATION") > 0) arc_hof_1 = reada(line,20)
          if (index(line, "FINAL GEOMETRY OBTAINED") > 0) exit   
          if (index(line, "GEOMETRY IN CARTESIAN COORDINATE") > 0) exit 
          if (index(line, "GEOMETRY IN MOPAC Z-MATRIX") > 0) exit            
        end do    
      end if
      nlines = 0 
      ncomments = 0
      if (.not. allocated(tmp_comments)) allocate(tmp_comments(10000))
      keywrd = " "
!
!  The size of the comments is not known, so set up a temporary array of size 10000 lines
!
      do        
        read (from_data_set, '(A1000)', end=30, err=30) line 
        nlines = nlines + 1 
        if (nlines == 1) then
          j = len_trim(line1)
          line1(:j) = line(:j)
          call upcase(line1, j)
          j = Index(line1,"DATA=")
          if (j > 0) exit
        end if
        if (line(1:1) /= '*') then
            line1 = trim(line)    
            i = 0
            do j = 1, len_trim(line1)
              if (ichar(line1(j:j)) == 9) then
      !
      ! convert tab to space(s).  Align with 8 character boundary
      !
                i = i + 1
                l = mod(i,8)
                line(i:) = " "
                if (l /= 0) i = i + 8 - l                    
              else
                i = i + 1
                line(i:i) = line1(j:j)
              end if
            end do
          write (input, '(A)', iostat=io_stat) trim(line)
          if (io_stat /= 0) then
            write (line, '(a)') ' The run-time temporary file "'//trim(jobnam)//'.temp" cannot be written to.'
            open(unit=iw, file=trim(jobnam)//'.out') 
            write(iw,"(a)")line
            if( .not. gui) write(0,"(///10x,a)")line
            call to_screen(line)
            call mopend (trim(line)) 
            return
          end if
          comments = .false.
          if (keywrd == " ") keywrd = line      
        else if (comments) then
          ncomments = ncomments + 1
          tmp_comments(ncomments) = trim(line)
          i = len_trim(line)
          if (i > 81) then
            line = trim(line(2:))
            if (index(line(:4),"ATOM") + index(line(:6),"HETATM") + index(line,"TITLE") + &
            index(line,"HEADER") + index(line,"ANISOU") + index(line,"COMPND") + &
            index(line,"SOURCE") + index(line,"KEYWDS") + index(line,"USER ") + &
            index(line,"HELIX") + index(line,"SHEET") + index(line,"REMARK") + &
            index(line, "SEQRES") /= 0) then
              open(unit=iw, file=trim(jobnam)//'.out') 
              num1 = char(ichar("2") + int(log10(ncomments*1.01)))
              num2 = char(ichar("2") + int(log10(i*1.01)))
              write(iw,'(/10x, a, i'//num1//', a,i'//num2//', a)') "Comment number", ncomments, &
                " is", i -1, " characters long. "
              write(iw,'(10x,a)')"This comment is recognized as being in Protein Data Bank format.", &
                "PDB lines are limited to 80 characters."
              write(iw,'(84x,a)')"Limit => |"
              write(iw,'(a)')"(Comment is:'*"//line(:i - 1)//"')"
              call mopend("Either make the line non-PDB format or shorten it. ")
              return
            end if
          end if
        end if
      end do
30    continue 
      call upcase(keywrd, len_trim(keywrd))
      if (keywrd(1:1) /= " ") keywrd = " "//trim(keywrd)
      i = index(keywrd, "GEO-DAT")
      if (i /= 0) keywrd(i:i+6) = "GEO_DAT"
      i = index(keywrd, "GEO-REF")
      if (i /= 0) keywrd(i:i+6) = "GEO_REF"
      if (index(keywrd, " GEO_DAT") /= 0) then
        nlines = nlines + 3
      else if (.not. is_PARAM .and. nlines < 4) then
        inquire(unit=iw, opened=exists) 
        if (.not. exists) open(unit=iw, file=trim(jobnam)//'.out') 
        if (keywrd /= " ") write(iw,'(3/10x,a,/)')" Data set does not contain "//&
           "any atoms and GEO_DAT is not present on keyword line"
      end if
      keywrd = "  "
      if (nlines == 1 .and. Len_trim(line1) > 0 .and. .not. is_PARAM) then  
!
!  Data set points to a MOPAC data set.  
!
       i = Len_trim(line1)
        job_fn = line(j+5:i)
        if(line1(i:i) == "+") then
!
!  File name runs onto two lines
!
           read (input, "(A120)", end=1000, err=1000) line
           i = i - j - 6
           job_fn(i:) = line
        end if
        i = Len_trim(line)
        if(line(i:i) == "+") then
!
!  File name runs onto three lines
!
           read (input, "(A120)", end=1000, err=1000) line
           i = Len_trim(job_fn) - 1
           job_fn(i:) = line
        end if
!
!  Make sure that the file ends in ".mop", ".dat", or "arc"
!
        i = Len_trim(job_fn) - 1
        line = job_fn
        call upcase(line,i)
        i = index(line,".MOP") + index(line,".DAT") + index(line,".ARC")
        if (i /= 0) job_fn(i + 4:) = " " 
        if (job_fn(1:1) == '"') job_fn = trim(job_fn(2:))
        do j = 1,Len_trim(job_fn) 
          if (job_fn(j:j) == char(92)) job_fn(j:j) = "/"
        end do
!
!  Clean up path before continuing
!
        do j = 1, Len_trim(jobnam)
          if (jobnam(j:j) == char(92)) jobnam(j:j) = "/"
        end do
        l = Index(job_fn,"/")
        if (l == 0) then
          open(unit=iw, file=trim(jobnam)//'.out') 
          call mopend ('INPUT FILE PATH MUST BE DEFINED IN ''DATA="<file plus path>"''')
          call web_message(iw,"DATA.html")
          return
        end if                       
        call to_screen ("Preparing to read the following MOPAC file: ")
        j = Len_trim (job_fn)
        if (j > 160) then
          call to_screen ("'"//job_fn(1:80))
          call to_screen (job_fn(81:160))
          call to_screen (job_fn(161:j)//"'")
        else if (j> 80) then
          call to_screen ("'"//job_fn(1:80))
          call to_screen (job_fn(81:j)//"'")
        else
          call to_screen ("'"//job_fn(1:j)//"'")
        end if
        line = job_fn
        goto 98
    end if   ! Line was one of first 5 lines
!
! Now the size of the comments list is known.  Allocate all_comments
!
      if ( .not. allocated(all_comments)) allocate(all_comments(ncomments + 100))
      do i = 1, ncomments
        all_comments(i) = trim(tmp_comments(i))
      end do
      deallocate(tmp_comments)
      line = ' ' 
      write (input, '(A241)') line 
      rewind input 
1000  if (nlines < 3 .and. .not. is_PARAM) then 
        inquire(unit=output, opened=exists)
        if (.not. exists) open(unit=output, file=trim(jobnam)//'.out') 
        call getarg (run, jobnam)
        write (0, '(A)') ' INPUT FILE "'//trim(jobnam)//'" MISSING OR EMPTY'  
        call mopend ( ' INPUT FILE "'//trim(jobnam)//'" MISSING OR EMPTY') 
        return  
      endif 
      natoms = nlines
      close(from_data_set) 
      return  
      end subroutine getdat 
 
