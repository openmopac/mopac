! Molecular Orbital PACkage (MOPAC)
! Copyright (C) 2021, Virginia Polytechnic Institute and State University
!
! MOPAC is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! MOPAC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

    subroutine grid
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use chanel_C, only : iw0, iw, iarc, ires, iump, archive_fn, ump_fn, &
      restart_fn
!
      use common_arrays_C, only : geo, xparam, geoa, labels, coord, &
      na, nb, nc, pa, pb, p
!      
      use elemts_C, only : elemnt
!
      use molkst_C, only : nvar, keywrd, tleft, line, norbs, &
      natoms, moperr, uhf, numat, mpack, gui
!
      use maps_C, only : rxn_coord1, rxn_coord2, ione, ijlp, ilp, jlp, jlp1, surf, &
      lpara1, latom1, lpara2, latom2
!
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: npts1, npts2, maxcyc, i, iloop, jloop, j, k, ij, l, &
      iw00, percent = 0, max_count, big_loop, loop, io_stat
      double precision :: step1, step2, degree, c1, c2, cputot, &
        escf, cpu1, cpu2, cpu3, geo11, geo22, sum
      character :: formt*4, num*1
      logical :: restrt, useef, opend, minimize_energy_in_grid = .false., use_p
      double precision, dimension (:), allocatable :: all_points1, all_points2
      double precision, dimension (:,:), allocatable ::  xy, surfac, all_pa, all_pb
      double precision, dimension (:,:,:), allocatable :: all_geo
      integer, dimension (:,:,:), allocatable :: all_nabc
      character :: txt*1
      double precision, external :: reada, seconds
!-----------------------------------------------
      if (allocated(geoa)) deallocate (geoa)
      allocate(geoa(3, natoms))
      geoa = geo
      iw00 = iw0
      iw0 = -1  !  While in "grid", do not print working in other subroutines
      useef = (Index (keywrd, " BFGS") == 0 .or. nvar < 2)
      i = Index (keywrd, "STEP1") + 6
      step1 = reada (keywrd, i)
      i = Index (keywrd, "STEP2") + 6
      step2 = reada (keywrd, i)
      npts1 = 11
      npts2 = 11
      maxcyc = 100000
      minimize_energy_in_grid = (index(keywrd," SMOOTH") /= 0)
      if (Index (keywrd, " BIGCYCLES") /= 0) then
        maxcyc = Nint (reada (keywrd, Index (keywrd, " BIGCYCLES")))
      end if
      if (Index (keywrd, "POINT1") /= 0) then
        npts1 = Nint (Abs(reada (keywrd, Index (keywrd, "POINT1")+7)))
      end if
      if (Index (keywrd, "POINT2") /= 0) then
        npts2 = Nint (Abs(reada (keywrd, Index (keywrd, "POINT2")+7)))
      end if
      restrt = (Index (keywrd, " RESTART") /= 0)
      if (allocated(surf))   deallocate(surf)
      allocate (surf(npts1*npts2), surfac(npts2, npts1), xy(2,npts1*npts2*4))
  !
  !  Basic "SMOOTH" grid - a four-fold ploughed field.
  !
  !   Left to right, down one, right to left, down one, repeat to end of grid.
  !   Reverse this plot
  !   Top to bottom, left one, bottom to top, left one, repeat to end of grid.
  !   Reverse this plot
  !
  !  Set up grid of points for the first pass
  !
      k = 0
      do i = 1,npts1
        do j = 1,npts2
          k = k + 1
          xy(1,k) = i - 1
          if (mod(i,2) == 1) then
            xy(2,k) = j - 1
          else
            xy(2,k) = npts2 - j
          end if
        end do
      end do
      ij = npts1*npts2
  !
  !  Set up grid of points for the second pass
  !
      do i = 1, ij
        k = k + 1
        xy(1,k) = xy(1,ij - i + 1)
        xy(2,k) = xy(2,ij - i + 1)
      end do
  !
  !  Set up grid of points for the third pass
  !
      do j = 1,npts2
        do i = 1,npts1
          k = k + 1
          xy(2,k) = j - 1
          if (mod(j,2) == 1) then
            xy(1,k) = i - 1
          else
            xy(1,k) = npts1 - i
          end if
        end do
      end do
  !
  !  Set up grid of points for the fourth and final pass
  !
      ij = npts1*npts2
      l = 3*ij + 1
      do i = 1, ij
        k = k + 1
        xy(1,k) = xy(1,l - i)
        xy(2,k) = xy(2,l - i)
      end do
     !
     !  THE TOP-LEFT VALUE OF THE FIRST AND SECOND DIMENSIONS ARE
     !      GEO(LPARA1,LATOM1) AND GEO(LPARA2,LATOM2)
     !
      degree = 57.29577951308232d0
      if (lpara1 /= 1) then
        step1 = step1 / degree
      end if
      if (lpara2 /= 1) then
        step2 = step2 / degree
      end if
      geo11 = geo(lpara1, latom1)
      geo22 = geo(lpara2, latom2)
      max_count = 0
      i = npts1*npts2
      if (uhf) then
        j = mpack
      else
        j = 1
      end if
      allocate (all_nabc(i,3,natoms), all_points1(i), all_points2(i), &
              & all_geo(i,3,natoms))
      allocate (all_pa(i,mpack), all_pb(i,j), stat = k)
      use_p = (i == 0)
      surf = 1.d9
      if (lpara1 /= 1 .and. na(latom1) > 0) then
        c1 = degree
      else
        c1 = 1.d0
      end if
      if (lpara2 /= 1 .and. na(latom2) > 0) then
        c2 = degree
      else
        c2 = 1.d0
      end if
      cputot = 0.d0
      if (restrt) then
          open(unit=ires, file=restart_fn, status='UNKNOWN', form='UNFORMATTED', &
          position='asis', iostat = io_stat)
          if (io_stat /= 0) then
            call mopend ("Restart file either does not exist or is not available for reading")
            return
          end if
          rewind ires
          if (nvar > 0) then
            read (ires, iostat = io_stat)i,j
            if ((norbs /= j .or. numat /= i) .and. (norbs /= i .or. numat /= j) ) then
              call mopend("Restart file read in does not match current data set")
              return
            end if
          end if
        if (useef .and. nvar > 1) then
!
!  Dummy read over the EF restart file
!
          do i = 1, 7
            read(ires, iostat = io_stat)txt
            if (txt == "x") j = 1
          end do
        else
          if (nvar > 0) then
            do i = 1, 3
              read(ires, iostat = io_stat)txt
            end do
          end if
        end if
        read(ires, iostat = io_stat) max_count, ijlp, ilp, jlp, jlp1, ione
        read(ires, iostat = io_stat) rxn_coord1, rxn_coord2
        read(ires, iostat = io_stat)all_geo(:max_count,:,:natoms)
        read(ires, iostat = io_stat)all_nabc(:max_count,:,:natoms)
        read(ires, iostat = io_stat)surf(:max_count), all_points1(:max_count), all_points2(:max_count)
      else
        ijlp = 1
      end if
      if ( minimize_energy_in_grid )  then
        big_loop = 4*npts1*npts2
      else
        big_loop = npts1*npts2
      end if
      if (iw00 > -1) &
      call to_screen('       FIRST VARIABLE   SECOND VARIABLE        FUNCTION      DONE LEFT')
 Main_Loop:  do loop = ijlp, big_loop
        if (loop >= maxcyc) then
          tleft = -100.d0
        end if
        geo (lpara1, latom1) = geo11 + xy(1,loop)*step1
        geo (lpara2, latom2) = geo22 + xy(2,loop)*step2
        cpu1 = seconds (2)
        if (useef .and. nvar > 1) then
          call ef (xparam, escf)
        else
          call flepo (xparam, nvar, escf)
        end if
        cpu2 = seconds (2)
        cpu3 = cpu2 - cpu1
        cputot = cputot + cpu3
        jlp = jlp + 1
        ijlp = ijlp + 1
 !
 !  Find the point in all_points1(1:max_count) and all_points2(1:max_count)
 !  defined by geo(lpara1, latom1) and geo(lpara2, latom2).
 !  If it does not exist, then create it.
 !
        k = 0
        do i = 1, max_count
          if( Abs(all_points1(i) - geo (lpara1, latom1)) < 1.d-6) then
            if( Abs(all_points2(i) - geo (lpara2, latom2)) < 1.d-6) then
 !
 !  Point located in "all" arrays
 !
              k = i
              exit
            end if
          end if
          if (k /= 0) exit
        end do
        if (k == 0) then
          max_count = max_count + 1
          k = max_count
          all_points1(k) = geo (lpara1, latom1)
          all_points2(k) = geo (lpara2, latom2)
        end if
 !
 !  If the current escf is lower than the stored value (surf(k)) then store
 !  the current value.
 !
        if (escf + 1.d-10 < surf(k)) then

 !
 !  Store all data for this point
 !
          surf(k) = escf
          do i = 1,natoms
            do j = 1,3
              all_geo(k,j,i) = geo(j,i)
            end do
            all_nabc(k,1,i) = na(i)
            all_nabc(k,2,i) = nb(i)
            all_nabc(k,3,i) = nc(i)
          end do
          if (gui) then
            if (use_p) then
              all_pa(k,:mpack) = pa(:mpack)
              if (uhf) then
                all_pb(k,:mpack) = pb(:mpack)
              end if
            end if
          end if
          write (line, '('' :'',F16.5,F16.5,F21.6, i10,i5)') &
            geo(lpara1,latom1)*c1, geo(lpara2,latom2)*c2, escf, &
            loop , big_loop - loop

          if (iw00 > -1) then
            i = nint((100.0*loop)/big_loop)
            if (i /= percent) then
              percent = i
              write(line,"(i4,a)")percent, "% of Grid Surface done"
              call to_screen(line)
            end if
            call to_screen(line)
          end if
          write(iw,'(a)')trim(line)
          if (index(keywrd, " DEBUG") /= 0) then
            if (index(keywrd, " PDBOUT") /= 0) then
              call pdbout(iw)
            else
              call geout(iw)
            end if
          end if
        end if
        if (tleft < 0.d0 .or. moperr) then
          if (tleft < 0.d0) then
            inquire(unit=ires, opened=opend)
            if (.not. opend) &
              open(unit=ires, file=restart_fn, status='UNKNOWN', form='UNFORMATTED', &
              position='asis', iostat = io_stat)
!
!  Save everything that would be needed by GRID when the job is restarted.
!
            write(ires, iostat = io_stat) max_count, ijlp, ilp, jlp, jlp1, ione
            write(ires, iostat = io_stat) rxn_coord1, rxn_coord2
            write(ires, iostat = io_stat)all_geo(:max_count,:,:natoms)
            write(ires, iostat = io_stat)all_nabc(:max_count,:,:natoms)
            write(ires, iostat = io_stat)surf(:max_count), all_points1(:max_count), all_points2(:max_count)
          end if
          deallocate (surf, surfac)
          return
        end if
      end do   Main_Loop

      if (minimize_energy_in_grid) write(iw,*)" Survey complete.  About to print points"
      do iloop = 1,npts1
        do jloop = 1, npts2
          geo (lpara1, latom1) = geo11 + (iloop - 1)*step1
          geo (lpara2, latom2) = geo22 + (jloop - 1)*step2
!
!  Find the correct point to print
!
          k = 0
          do i = 1, max_count
            if( Abs(all_points1(i) - geo (lpara1, latom1)) < 1.d-6) then
              if( Abs(all_points2(i) - geo (lpara2, latom2)) < 1.d-6) then
!
!  Point has been located in "all" arrays
!
                k = i
                exit
              end if
            end if
          if (k /= 0) exit
          end do
!
!  Restore all data for this point
!
          escf = surf(k)
          do i = 1,natoms
            do j = 1,3
              geo(j,i) = all_geo(k,j,i)
            end do
            na(i) = all_nabc(k,1,i)
            nb(i) = all_nabc(k,2,i)
            nc(i) = all_nabc(k,3,i)
          end do
          surfac(jloop, iloop) = escf
          write (iw, "(/'       FIRST VARIABLE   ',  'SECOND VARIABLE FUNCTION')")
          write (iw, "(' :',F16.5,F16.5,F16.6)") geo(lpara1, latom1)*c1, geo(lpara2, latom2)*c2, escf
          call geout (iw)
          if (index(keywrd, " PRTXYZ") /= 0) then
            write (iw, '(29X,''CARTESIAN COORDINATES '',/)')
            call gmetry(geo, coord)
            l = 0 
            do i = 1, natoms 
              if (labels(i) == 99 .or. labels(i) == 107) cycle
              l = l + 1
              write (iw, '(I4,3X,A2,3x, 3F16.9)') l, elemnt(labels(i)), (coord(k,l),k=1,3) 
            end do
          end if
          if (gui) then
            if (use_p) then
              ij = 0
              do i = 1, norbs
                do j = 1, i
                ij = ij + 1
                pa(ij) = all_pa(k,ij)
                end do
              end do
              if (uhf) then
                ij = 0
                do i = 1, norbs
                  do j = 1, i
                  ij = ij + 1
                  pb(ij) = all_pb(k,ij)
                  end do
                end do
              else
                pa = 0.d0
                pb = 0.d0
              end if
              p = pa + pb
            else
              p = 2*pa
            end if
            call bonds()
          end if
        end do
      end do
      write (iw, "(/10x,'HORIZONTAL: VARYING SECOND PARAMETER,',/10x, &
   & 'VERTICAL:   VARYING FIRST PARAMETER')")
    write (iw, "(/10X,'WHOLE OF GRID, SUITABLE FOR PLOTTING',//)")
   !
   !  ARCHIVE
      inquire (unit=iarc, opened=opend)
      if (opend) then
        close (unit=iarc, status="KEEP")
      end if
      open (unit=iarc, file=archive_fn, status="UNKNOWN")
10000 format (" ARCHIVE FILE FOR GRID CALCULATION"/"GRID OF HEATS" /)
      write (iarc, 10000)
      call wrttxt (iarc)
      write (iarc, "(/' TOTAL JOB TIME IN FLEPO : ',F10.3/)") cputot
   !
   !  WRITE OUT THE GRIDS
   !
      sum = 0.d0
      do i = 1, npts1
        do j = 1, npts2
          if (abs(surfac(j,i)) > sum) sum = abs(surfac(j,i))
        end do
      end do
      if (sum > 0.9999999d3) then
        num = char(ichar("5") + int(log10(sum)))
        formt = " "//num//".1"
      else
        num = char(ichar("7") + int(log10(sum)))
        formt = " "//num//".3"
      end if
      write(iw  ,'('//num//'x,1000f'//formt//')') ((geo22 + (j - 1)*step2)*c2,j=1,npts2)
      write(iarc,'('//num//'x,1000f'//formt//')') ((geo22 + (j - 1)*step2)*c2,j=1,npts2)
      open(unit=iump, file=ump_fn, status='UNKNOWN', position='asis')
      do i = 1, npts1
        write (iw,   '(1000f'//formt//')') (geo11 + (i - 1)*step1)*c1, (surfac(j,i),j=1,npts2)
        write (iarc, '(1000f'//formt//')') (geo11 + (i - 1)*step1)*c1, (surfac(j,i),j=1,npts2)
        do j = 1, npts2
          write (iump, '(2f9.3,F'//formt//')')(rxn_coord1 + step1*(i - 1))*c1, &
          (rxn_coord2 + step2*(j - 1))*c2, surfac(j,i)
        end do
      end do
      close(iump)
    if (allocated(surf)) then
      deallocate (surf, surfac)
    end if
    if (allocated(all_nabc)) &
    deallocate (all_nabc, all_points1, all_points2, all_geo)
    iw0 = iw00
end subroutine grid
