! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
 
 subroutine parkey (keywrd)
    use param_global_C, only : ifiles_8
    use molkst_C, only : method_PM7, method_PM6_ORG
    implicit none
    character (len=*), intent (in) :: keywrd
    integer :: i, j, k
    character :: num*1, line*300
    character, external :: get_a_name*300
    integer, external :: end_of_keyword
    double precision, external :: reada
     if (Index (keywrd, " MNDO ") /= 0) &
       write (ifiles_8, "(' *  MNDO        -  Default parameter set: MNDO')")
     if (Index (keywrd, " AM1") /= 0) &
       write (ifiles_8, "(' *  AM1         -  Default parameter set: AM1')")
     if (Index (keywrd, " PM3") /= 0) &
       write (ifiles_8, "(' *  PM3         -  Default parameter set: PM3')")
     if (Index (keywrd, " RM1") /= 0) &
       write (ifiles_8, "(' *  RM1         -  Default parameter set: RM1')")
     if (Index (keywrd, " PM6-DH") /= 0)  then
       write (ifiles_8, "(' *  PM6-DH2     -  Default parameter set: PM6-DH2')")
! Gallo
! These lines only print in the output file the chosen method. That way,
! filename.out will also acknowledge the use of PM6-D3H4 and PM6-D3(H4)
     else if (Index (keywrd, " PM6-D3H4") /= 0)  then
       write (ifiles_8, "(' *  PM6-D3H4    -  Default parameter set: PM6-D3H4')")
     else if (Index (keywrd, " PM6-D3(H4)") /= 0)  then
       write (ifiles_8, "(' *  PM6-D3(H4)  -  Default parameter set: PM6-D3(H4)')")
! Gallo
     else  if (Index (keywrd, " PM6 ") /= 0)  then
       write (ifiles_8, "(' *  PM6         -  Default parameter set: PM6')")
     end if
     if (Index (keywrd, " MNDOD") /= 0) &
       write (ifiles_8, "(' *  MNDOD       -  Default parameter set: MNDOD')")
     if (method_PM7) write (ifiles_8,     "(' *  PM7         -  Default parameter set: PM7')")
     if (method_PM6_ORG) write (ifiles_8, "(' *  PM6-ORG     -  Default parameter set: PM6-ORG')")
     if (Index (keywrd, " PM5D ") /= 0) &
       write (ifiles_8, "(' *  PM5D        -  Default parameter set: PM5D')")
      if (Index (keywrd, "LARGE") /= 0) then
10000 format (" *  LARGE       -  GENERATE STANDARD MOPAC OUTPUT")
      write (ifiles_8, 10000)
    end if
    if (Index (keywrd, "SAVEG") /= 0) then
10010 format (" *  SAVEG      -  SAVE OPTIMIZED GEOMETRIES")
      write (ifiles_8, 10010)
    end if
    if (Index (keywrd, " NOSAVEP") /= 0) then
10011 format (" *  NOSAVEP    -  DO NOT SAVE OPTIMIZED PARAMETERS")
      write (ifiles_8, 10011)
    end if
    if (Index (keywrd, "FULL") /= 0) then
10021 format (" *  FULL       -  PRINT ALL DERIVATIVES")
      write (ifiles_8, 10021)
    end if
    if (Index (keywrd, "PARALLEL") /= 0) then
10022 format (" *  PARALLEL   -  COMBINE PARAMETERS FROM OTHER RUNNING PARAM JOBS")
      write (ifiles_8, 10022)
    end if
     if (Index (keywrd, "EXPORT") /= 0) &
     & write (ifiles_8, '(" *  EXPORT     -  REMOVE REFERENCE INFO FROM NEW DATA FILE")')
    if (Index (keywrd, "LET") /= 0) then
10040 format (" *  LET        -  CARRY OUT CALCULATION EVEN IF ", "SOME SYSTEMS &
     &ARE FAULTY")
      write (ifiles_8, 10040)
    end if
    if (Index(keywrd, " GNORM") /= 0) then
10041 format (" *  GNORM=     -  EXIT WHEN GRADIENT NORM DROPS BELOW ", g10.3)
      write (ifiles_8, 10041) reada (keywrd, Index (keywrd, " GNORM"))
    end if
    if (Index (keywrd, "CLEAN") /= 0) write(ifiles_8,'(a)') &
    & " * CLEAN       -  REMOVE DIAGNOSTIC KEYWORDS FROM REFERENCE DATA SETS"
    if (Index (keywrd, "ALL") /= 0) then
10042 format (" *  ALL         -  USE ALL REFERENCE DATA, INCLUDING SYSTEMS NOT BEING OPTIMIZED")
      write (ifiles_8, 10042)
    end if
    if (Index (keywrd, "ONLY") /= 0) then
10043 format (" *  ONLY        -  Use only compounds that contain elements being parameterized.")
      write (ifiles_8, 10043)
    end if
    if (Index (keywrd, "1SCF") /= 0) then
10050 format (" *  1SCF       -  GEOMETRIES ARE NOT TO BE OPTIMIZED ")
      write (ifiles_8, 10050)
    end if
    if (Index (keywrd, " PRECISE") /= 0) then
10055 format (" *  PRECISE    -  USE ""PRECISE"" IN MOPAC")
      write (ifiles_8, 10055)
    end if
    i = Index (keywrd, " PARAB=")
    if (i /= 0) then
10061 format (" *  PARAB=n     -  DAMPENING FACTOR OF",f13.2," USED")
      write (ifiles_8, 10061) Reada(keywrd,i)
    end if
     i = Index (keywrd, " POWER=")
    if (i /= 0) then
10063 format (" *  POWER=n    -  Power to be minimized",f5.2," USED")
      write (ifiles_8, 10063) Reada(keywrd,i)
    end if
    i = Index (keywrd, " RELH=")
    if (i /= 0) then
10051 format (" *  RELH=n     -  Relative weight for Heats of Formation:",f4.1)
      write (ifiles_8, 10051) Reada(keywrd,i)
    end if
    i = Index (keywrd, " RELD=")
    if (i /= 0) then
10052 format (" *  RELD=n     -  Relative weight for Dipole Moments:",f4.1)
      write (ifiles_8, 10052) Reada(keywrd,i)
    end if
     i = Index (keywrd, " RELI=")
    if (i /= 0) then
10053 format (" *  RELI=n     -  Relative weight for Ionization Potentials:",f4.1)
      write (ifiles_8, 10053) Reada(keywrd,i)
    end if
     i = Index (keywrd, " RELG=")
    if (i /= 0) then
10054 format (" *  RELG=n     -  Relative weight for Geometries:",f4.1)
      write (ifiles_8, 10054) Reada(keywrd,i)
    end if
    if (Index (keywrd, " NOSCALE") /= 0)  write (ifiles_8, '(a)') &
        " *  NOSCALE     -  DO NOT REDUCE WEIGHTING WHEN PREDICTED VALUES ARE IN ERROR"
    i = Index (keywrd, " CYCLES=")
    if (i /= 0) then
      j = Nint(Reada(keywrd,i))
      num = char(ichar("2") +int(log10(j + 0.05)))
      write (ifiles_8, '(" *  CYCLES=n    -  RUN A MAXIMUM OF",I'//num//'," CYCLES")') j
    end if
     i = Index (keywrd, " NUMOPT=")
    if (i /= 0) then
      j = Nint(Reada(keywrd,i))
      num = char(ichar("2") +int(log10(j + 0.05)))
      write (ifiles_8, '(" *  NUMOPT=n    -  USE THE",I'//num//'," PARAMETERS WITH THE LARGEST WEIGHTED GRADIENTS")') j
    end if
    i = Index (keywrd, " OUT=")
    if (i /= 0) then
      j = Index (keywrd(i+1:), " ") + i
10060 format (" *  OUT=n       -  WRITE OUTPUT TO DIRECTORY ",/19x, a)
      write (ifiles_8, 10060) keywrd(i+5:j)
    end if
    i = Index (keywrd, " PARAMS=")
    if (i /= 0) then
      j = Index (keywrd(i+1:), " ") + i
10071 format (" *  PARAMS=n    -  REFERENCE PARAMETERS TO BE READ FROM DIRECTORY ",/16 x, a)
      write (ifiles_8, 10071) keywrd(i+8:j)
    end if
    i = Index (keywrd, " NEW_REF=")
    if (i /= 0) then
      line = get_a_name(keywrd(i + 9:), len_trim(keywrd(i + 9:)))
      write (ifiles_8, '(" *  NEW_REF=n   -  OPTIMIZED REFERENCE DATA ARE WRITTEN TO DIRECTORY: ",/16 x, a)') trim(line)
      call l_control("CYCLES=0", len("CYCLES=0"), 1)
    end if
    i = Index (keywrd, " NEW_REFP")
    if (i /= 0) then
      j = Index (keywrd(i+1:), " ") + i
10081 format (" *  NEW_RP=n    -  OPTIMIZED PARAMETERS ARE WRITTEN TO DIRECTORY ",/16 x, a)
      write (ifiles_8, 10081) keywrd(i+14:j)
    end if
    if (Index (keywrd, "CHKPAR") /= 0) then
10090 format (" *  CHKPAR      -  CHECK THAT PARAMETER SET IS VALID")
      write (ifiles_8, 10090)
    end if
    if (Index (keywrd, " SURVEY") /= 0) then
10091 format (" *  SURVEY      -  GENERATE TABLES FOR STATISTICS, ETC")
      write (ifiles_8, 10091)
    end if
    if (Index (keywrd, "NOLIM") /= 0) write(ifiles_8,'(" *  NOLIM       -  REMOVE LIMIT CONSTRAINTS ON PARAMETER VALUES")')
    if (index(keywrd, " SET") + index(keywrd, " REF") /= 0) &
      write (ifiles_8, "(' *',/1X,15('*****'))")
    i = Index (keywrd, " SET")
    if (i /= 0) then
      i = i + 5
      line = get_a_name(keywrd(i:), len_trim(keywrd(i:)))
      j = i + len_trim(line)
      if (keywrd(j + 1: j + 1) == '"') j = j + 2  
      if (keywrd(j:j) == ";") then
        write (ifiles_8, '(" *",/," *  SET=n       -  THE LISTS OF COMPOUNDS TO BE USED ARE IN FILES: ")')
      else
        write (ifiles_8, '(" *",/," *  SET=n       -  THE LIST OF COMPOUNDS TO BE USED IS IN FILE: ")')
      end if
      write (ifiles_8, '(" *",17x,a)')'"'//trim(line)//'"'
      if (keywrd(i + 1: i + 1) == '"') i = i + 2     
      do
        i = i + len_trim(line)
        if (keywrd(i + 1: i + 1) == '"') i = i + 2  
        if (keywrd(i:i) == ";") then
          i = i + 1 
          line = get_a_name(keywrd(i:), len_trim(keywrd(i:)))
          write (ifiles_8, '(" *", 10x, a)')'   and "'//trim(line)//'"'
        else
          exit
        end if
      end do      
    end if
    i = Index (keywrd, " REF=")
    if (i /= 0) then
      i = i + 5
      k = end_of_keyword(keywrd, len_trim(keywrd), i)
      line = get_a_name(keywrd(i:k), len_trim(keywrd(i:k)))
      write (ifiles_8, '(" *",/," *  REF=n       -  REFERENCE DATA ARE IN DIRECTORY: ",/," *",17x, a)') '"'//trim(line)//'"'
      do
        j = index(keywrd(i:k), ";")
        if (j /= 0) then
          i = i + j 
          line = get_a_name(keywrd(i:), len_trim(keywrd(i:)))
          write (ifiles_8, '(" *", 10x, a)')'   and "'//trim(line)//'"'
        else
          exit
        end if
      end do      
    end if
    i = Index(keywrd, "EXTERNAL=") + Index(keywrd, "PARAMS=")
    if (i /= 0) then
      i = index(keywrd(i:), "=") + i
      k = end_of_keyword(keywrd, len_trim(keywrd), i)
      line = get_a_name(keywrd(i:k), len_trim(keywrd(i:k)))
      write (ifiles_8, '(" *",/," *  EXTERNAL=n  -  DEFAULT PARAMETERS RESET USING DATA IN FILES: ",/," *",17x, a)') '"'&
      &//trim(line)//'"'
      do
        j = index(keywrd(i:k), ";")
        if (j /= 0) then
          i = i + j 
          line = get_a_name(keywrd(i:), len_trim(keywrd(i:)))
          write (ifiles_8, '(" *", 10x, a)')'   and "'//trim(line)//'"'
        else
          exit
        end if
      end do      
    end if
    return
  end subroutine parkey
