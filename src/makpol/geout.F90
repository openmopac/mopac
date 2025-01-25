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

subroutine geout (mode1, geo, na, nb, nc, labels, loc)
!*********************************************************************
!
!   GEOUT PRINTS THE CURRENT GEOMETRY.  IT CAN BE CALLED ANY TIME,
!         FROM ANY POINT IN THE PROGRAM AND DOES NOT AFFECT ANYTHING.
!         IF MODE1 .EQ.1 THEN GEOMETRY IS PRINTED IN USUAL .OUT FORMAT
!                  .GT.1 THEN GEOMETRY IS PRINTED IN .DAT FORMAT to
!                  CHANNEL MODE1
!
!*********************************************************************
    use common_systm, only : natoms, iw, maxtxt, ndep, numat, line
    use common_symult, only : depmul
    use common_keywrd, only : keywrd
    use common_geosym, only : locpar, idepfn, locdep
    use common_elemts, only : elemnt
    implicit none
    integer, intent (in) :: mode1, labels(natoms)
    integer, intent (in) :: na(natoms), nb(natoms), nc(natoms)
    integer, intent (inout) :: loc(2,3*natoms + 6)
    double precision, dimension (3, natoms), intent (in) :: geo
    integer :: mode, iprt, i, j, n, ia, ii, k, igui, store_maxtxt
    double precision :: degree, w, x 
    logical :: cart, lxyz, charge, gui = .false.
    character , dimension(3) :: q*2 
    character :: flag1*2, flag0*2,  blank*80, fmt1*4, fmt23*4
!*********************************************************************
!
!   GEOUT PRINTS THE CURRENT GEOMETRY.  IT CAN BE CALLED ANY TIME,
!         FROM ANY POINT IN THE PROGRAM AND DOES NOT AFFECT ANYTHING.
!
!   mode1:   1 write geometry in normal MOPAC *.out format
!           -n write geometry in normal MOPAC *.arc format, but do not print keywords, title
!              symmetry data, etc.
!            n write geometry in normal MOPAC *.arc format
!
!*********************************************************************
      mode = mode1 
      igui = -10 ! Set to impossible value
      store_maxtxt = maxtxt
      charge = (index(keywrd, " PRTCHAR") /= 0)
      lxyz = index(keywrd,' COORD') + index(keywrd,'VELO') /= 0 
      if (index(keywrd,' 0SCF')/=0 .and. lxyz) then 
!
!  If 0SCF and coord and a polymer, then get rid of TV.
!
        natoms = numat 
        lxyz = .FALSE. 
      endif 
      if (mode == 1) then 
        flag1 = ' *' 
        flag0 = '  ' 
        iprt = iw 
      else 
        if (gui) then
          flag1 = ' 1' 
          flag0 = ' 0' 
        else
          flag1 = '+1' 
          flag0 = '+0' 
        end if
        iprt = abs(mode) 
      endif 
      degree = 57.29577951308232D0 
      if (lxyz) degree = 1.D0 
      blank = ' '  
      cart = .true.
      do i = 1, natoms
        if (na(i) > 0) cart = .false.
      end do     
      fmt1  = "13.8"
      fmt23 = "13.7"
      if (maxtxt /= 0)  maxtxt = maxtxt + 2
      if (cart) then 
        x = 0.d0
        do i = 1, natoms
          x = min(x, geo(1,i), geo(2,i), geo(3,i))
        end do
        if (x < -99.99d0) then
          fmt1  = "14.8"
          fmt23 = "14.8"
        else
          fmt23 = "13.8"
        end if
        if (mode == 1) then
          if (maxtxt == 0) then
            i = 6
          else
            i = maxtxt/2 + 1
          end if
          write (iprt,'(4a)') &
             & "   ATOM "//blank(:maxtxt/2 + 1)//" CHEMICAL  "//blank(:i)//"  X               Y               Z"
          write (iprt,'(4a)') "  NUMBER "//blank(:maxtxt/2 + 1)//" SYMBOL" & 
          &//blank(:i)//"(ANGSTROMS)     (ANGSTROMS)     (ANGSTROMS)"
          write (iprt,*)
        else if (mode > 0) then         
          call wrttxt (iprt)
        end if
      else if (mode == 1) then
        j = max(9,maxtxt + 2) 
        if (maxtxt == 0) j = 8
        write (iprt,'(4a)') &
           & "  ATOM"//blank(:j/2)//"CHEMICAL "//blank(:(j + 1)/2)//" BOND LENGTH      BOND ANGLE     TWIST ANGLE "
        write (iprt,'(4a)') &
           & " NUMBER"//blank(:j/2)//"SYMBOL   "//blank(:(j + 1)/2)//"(ANGSTROMS)      (DEGREES)       (DEGREES) "
        write (iprt, "(A)") &
           & "   (I)       "//blank(:j)//"      NA:I           NB:NA:I    " // &
           & "   NC:NB:NA:I " // "      NA    NB    NC "
      else if (mode > 0) then
        call wrttxt (iprt)
      end if 
      if (.false. .and. cart) then
        n = 0
        do i = 1, natoms
          do j = 1, 3
            n = n + 1
            loc(1,n) = i
            loc(2,n) = j
          end do
        end do
      end if
      n = 1 
      ia = loc(1,1) 
      ii = 0 
      blank = " "
      do i = 1, natoms 
        do j = 1, 3 
          q(j) = flag0 
          if (ia /= i) cycle  
          if (j/=loc(2,n)) cycle  
          q(j) = flag1 
          n = n + 1 
          ia = loc(1,n) 
        end do         
        if (na(i) > 0) then
          w = geo(2,i)*degree 
          x = geo(3,i)*degree  
!
!  CONSTRAIN ANGLE TO DOMAIN 0 - 180 DEGREES
!
          w = w - aint(w/360.D0)*360.D0 
          if (w < -1.d-6) w = w + 360.D0 
          if (w > 180.000001D0) then 
            x = x + 180.D0 
            w = 360.D0 - w 
          endif 
!
!  CONSTRAIN DIHEDRAL TO DOMAIN -180 - 180 DEGREES
!
          x = x - aint(x/360.D0 + sign(0.5D0 - 1.D-9,x) - 1.D-9)*360.D0 
        else
          w = geo(2,i) 
          x = geo(3,i)  
        endif 
       
        blank = elemnt(labels(i))
        if (labels(i) /= 99 .and. labels(i) /= 107) ii = ii + 1 
      
        k = 4
        
        k = 0
      
        if (index(blank(:j),"(") /= 0) then
          if (index(blank(:j),")") == 0) then
            if (index(blank(j + 1:j + 1), ")") /= 0) j = j + 1
          end if
        end if
        if (labels(i) == 0) cycle  
    
        if (na(i) == igui .or. cart) then 
          if (mode /= 1) then !  Print suitable for reading as a data-set
            if (labels(i)/=99 .and. labels(i)/=107) then 
              if (maxtxt == 0 ) j = 4    
              if (charge) then
                write(line,'(a,f8.4)')blank(41:59 + k)
              else
                write(line,'(a,f8.4)')blank(41:59 + k)
              end if
              write (iprt, '(1X,A,F'//fmt1//',1X,A2,F'//fmt23//',1X,A2,F'//fmt23//',1X, A2, A)') &
                 blank(:j), geo(1,i), q(1), w, q(2), x, q(3),trim(line)
            else 
              write (iprt, '(1X,A,F'//fmt1//',1X,A2,F'//fmt23//',1X,A2,F'//fmt23//',1X,A2,a)') &
              blank(:j), geo(1,i), q(1), w, q(2), x, q(3), " "
            endif 
          else !  Print in output style
            if (maxtxt == 0) j = 9
            write (iprt, '(I6,6X,A,F'//fmt1//',1X,A2,F'//fmt23//',1X,A2,F'//fmt23//',1X,A2)') i&
              , blank(:j), geo(1,i), q(1), w, q(2), x, q(3) 
          endif 
        else 
          if (mode /= 1) then  !  Print suitable for reading as a data-set
            if (maxtxt == 0 ) j = 4
            if (labels(i) /= 99 .and. labels(i) /= 107) then 
              if (charge) then
                write(line,'(a,f8.4)')blank(41:41 + k)
              else
                write(line,'(a,f8.4)')blank(41:41 + k)
              end if
              write (iprt, &
                '(1X,A,F'//fmt1//',1X,A2,F'//fmt23//',1X,A2,F'//fmt23//',1X,A2,3I6,A)') &
                & blank(:j), geo(1,i), q(1), w, q(2), x, q(3), na(i), nb(i), nc(i), &
                trim(line) 
            else 
              write (iprt, '(1X,A,F'//fmt1//',1X,A2,F'//fmt23//',1X,A2,F'//fmt23//',1X,A2,3I6)') &
                blank(:j), geo(1,i), q(1), w, q(2), x, q(3), na(i), nb(i), nc(i) 
            endif 
          else !  Print in output style
            if (maxtxt == 0) j = 9
            write (iprt, &
              '(I6,6X,A,F'//fmt1//',1X,A2,F'//fmt23//',1X,A2,F'//fmt23//',1X,A2,I6,2I6)') i, &
              blank(:j), geo(1,i), q(1), w, q(2), x, q(3), na(i), nb(i), nc(i) 
          endif 
        endif 
      end do 
      maxtxt = store_maxtxt
      if (mode == 1) return  
      write (iprt, *) 
      if (ndep /= 0) then  
!
!   OUTPUT SYMMETRY DATA.
!
        n = 1
        i = 1
        outer_loop: do
          j = i
          do
            if (j == ndep) exit outer_loop
             !
             !  Group together symmetry functions of the same type
             !  (same reference atom, same reference function, same multiplier,
             !   if function 18 or 19)
             !  (Maximum number of dependent atoms on a line: 9)
             !
            if (locpar(j) /= locpar(j+1) .or. idepfn(j) /= idepfn(j+1) .or. &
                 & j-i >= 9) exit
            if (idepfn(i) == 18 .or. idepfn(i) == 19) then
              if (Abs(depmul(n) - depmul(n+1)) > 1.d-10) exit
              n = n + 1
            end if
            j = j + 1
          end do
          if (idepfn(i) == 18 .or. idepfn(i) == 19) then
            write (iprt, "(I4,I3,F13.9,10I5)") locpar (i), idepfn (i), &
                 & depmul(n), (locdep(k), k=i, j)
            n = n + 1
          else
            write (iprt, "(I4,I3,10I5)") locpar (i), idepfn (i), &
                 & (locdep(k), k=i, j)
          end if
          i = j + 1
        end do outer_loop
        if (idepfn(i) == 19 .or. idepfn(i) == 18) then
          write (iprt, "(I4,I3,F13.9,10I5)") locpar (i), idepfn (i), &
               & depmul(n), (locdep(k), k=i, j)
        else
          write (iprt, "(I4,I3,10I5)") locpar (i), idepfn (i), &
               & (locdep(k), k=i, j)
        end if
        write (iprt,*)
      end if
     
      end subroutine geout 
