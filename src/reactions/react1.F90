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

      subroutine react1
      USE molkst_C, only : keywrd, natoms, moperr, numat, nvar, keywrd_txt, &
      tleft,  gnorm, escf, cosine, mpack, step, time0, iflepo, line, density
      use chanel_C, only : iw, ir
      use elemts_C, only : elemnt
      use common_arrays_C, only : p, pa, pb, geo, geoa, coord, &
      & labels, na, nb, nc, xparam, grad, nat, loc
      implicit none
      integer :: linear, iflag, i, maxstp, maxcyc, j, k, iloop, lopt(3,natoms)
      double precision, dimension(3*numat) :: xold, grold
      double precision :: stepmx,  x, sumx, sumy, sumz, sum, &
        step0, one, dell, eold, time1, swap, funct1, time2, gold, c1, c2, &
        dist
      double precision, dimension(:,:), allocatable :: xyz
      double precision, dimension(:), allocatable ::  pastor, pbstor
      logical :: gradnt, finish, intl
      logical , dimension(2) :: gok
      logical :: lef
      double precision, external :: ddot, reada, seconds
!-----------------------------------------------
!***********************************************************************
!
!  REACT1 DETERMINES THE TRANSITION STATE OF A CHEMICAL REACTION.
!
!   REACT WORKS BY USING TWO SYSTEMS SIMULTANEOUSLY, THE HEATS OF
!   FORMATION OF BOTH ARE CALCULATED, THEN THE MORE STABLE ONE
!   IS MOVED IN THE DIRECTION OF THE OTHER. AFTER A STEP THE
!   ENERGIES ARE COMPARED, AND THE NOW LOWER-ENERGY FORM IS MOVED
!   IN THE DIRECTION OF THE HIGHER-ENERGY FORM. THIS IS REPEATED
!   UNTIL THE SADDLE POINT IS REACHED.
!
!   IF ONE FORM IS MOVED 3 TIMES IN SUCCESSION, THEN THE HIGHER ENERGY
!   FORM IS RE-OPTIMIZED WITHOUT SHORTENING THE DISTANCE BETWEEN THE TWO
!   FORMS. THIS REDUCES THE CHANCE OF BEING CAUGHT ON THE SIDE OF A
!   TRANSITION STATE.
!
!***********************************************************************
      linear = 0
      density = 0.d0
      iflag = 1
      gok(1) = .FALSE.
      gok(2) = .FALSE.
      lef = index(keywrd,' EF') /= 0
      gradnt = index(keywrd,'GRAD') /= 0
      i = index(keywrd,' BAR')
      stepmx = 0.01D0
      if (i /= 0) stepmx = reada(keywrd,i)
      maxstp = 10000
      allocate(pastor(mpack), pbstor(mpack))
      if (index(keywrd_txt, "GEO_REF") == 0) then
!
!    READ IN THE SECOND GEOMETRY.
!
        if (allocated (geoa)) deallocate(geoa)
        allocate(geoa(3,natoms))
        i = numat
        call getgeo (ir, labels, geoa, coord, lopt, na, nb, nc, intl)
        if (numat /= i) then
          rewind (ir)
          write(iw,"(a,i4)")" Number of atoms in the first geometry: ",i, &
                            " Number of atoms in the second geometry:", numat
          write (iw,"(a)") " Second geometry has different number of atoms from first geometry", &
        " Data set suppled:"
          do i = 1, 100000
            read(ir,"(a)", iostat = j) line
            if (j < 0) exit
            write(iw,"(a)") Trim (line)
          end do
          call mopend("Geometry supplied is faulty - see output file")
  ! Dummy use of intl
          if (intl .or. lopt(1,1) /= 0) return
          return
        end if
      else
        call gmetry (geoa, coord)
      end if
      call delete_ref_key("GEO_REF=", len_trim("GEO_REF="), '" ', 2)
      call delete_ref_key("SADDLE", len_trim("SADDLE"), ' ', 1)
      call delete_ref_key("BAR", len_trim("BAR"), ' ', 1)
      intl = (index(keywrd, " XYZ") == 0)
      maxcyc = 100000
      if (Index (keywrd, " BIGCYCLES") /= 0) &
        maxcyc = Nint (reada (keywrd, Index (keywrd, " BIGCYCLES")))
   !
   !  SWAP FIRST AND SECOND GEOMETRIES AROUND
   !  SO THAT GEOUT CAN OUTPUT DATA ON SECOND GEOMETRY.
   !
      allocate(xyz(3,numat))
      xyz = geo(:,:numat)
      geo(:,:numat) = geoa(:,:numat)
      geoa(:,:numat) = xyz
      sumx = 0.d0
      sumy = 0.d0
      sumz = 0.d0
      do j = 1, numat
        sumx = sumx + xyz(1, j)
        sumy = sumy + xyz(2, j)
        sumz = sumz + xyz(3, j)
      end do
      sumx = sumx / numat
      sumy = sumy / numat
      sumz = sumz / numat
      do j = 1, numat
        geoa(1, j) = xyz(1, j) - sumx
        geoa(2, j) = xyz(2, j) - sumy
        geoa(3, j) = xyz(3, j) - sumz
      end do
      write (iw, "(//,'  CARTESIAN GEOMETRY OF FIRST SYSTEM',//)")
      write (iw, "(i4,3x,a2,3x,3F14.5)") (i,elemnt(nat(i)),(geo(j, i), j=1, 3), i=1, numat)
      call dock (geoa, geo, dist)
      write (iw, "(//,'  CARTESIAN GEOMETRY OF SECOND SYSTEM',//)")
        write (iw, "(i4,3x,a2,3x,3F14.5)") (i,elemnt(nat(i)),(geoa(j, i), j=1, 3), i=1, numat)
      write (iw, "(//,'   ""DISTANCE"":',F13.6)") dist
      write (iw, "(//,'  REACTION COORDINATE VECTOR',//)")
      write (iw, "(i4,3x,a2,3x,3F14.5)") (i,elemnt(nat(i)), &
      (geo(j, i) - geoa(j, i), j=1, 3), i=1, numat)
      !
      !   XPARAM HOLDS THE VARIABLE PARAMETERS FOR GEOMETRY IN GEO
      !   XOLD   HOLDS THE VARIABLE PARAMETERS FOR GEOMETRY IN GEOA
      !
      sum = 0.d0
      i = 0
      do j = 1, numat
        do k = 1,3
          i = i + 1
          grold(i) = 1.d0
          xparam(i) = geo(k,j)
          xold(i) = geoa(k,j)
          sum = sum + (xparam(i)-xold(i)) ** 2
        end do
      end do
      step0 = Sqrt (sum)
      if (step0 < 1.d-2) then
        write (iw, "(//,3(5X,A,/))") " THE TWO GEOMETRIES ARE IDENTICAL OR ALMOST IDENTICAL", &
       & " A SADDLE CALCULATION INVOLVES A REACTANT AND A PRODUCT", &
       & " THESE MUST BE DIFFERENT GEOMETRIES"
        call mopend ("THE TWO GEOMETRIES IN SADDLE ARE IDENTICAL.")
        return
      end if
      time0 = seconds(2)
      maxcyc = 100000
      if (index(keywrd,' BIGCYCLES') /= 0) maxcyc = nint(reada(keywrd,index(&
        keywrd,' BIGCYCLES')))

      one = 1.D0
      dell = 0.1D0
      eold = -2000.D0
      time1 = seconds(2)
      swap = 0.D0
      do iloop = 1, maxstp
        if (iloop >= maxcyc) tleft = -100.D0
        write (iw, '('' '',40(''*+''))')
!
!   THIS METHOD OF CALCULATING 'STEP' IS QUITE ARBITARY, AND NEEDS
!   TO BE IMPROVED BY INTELLIGENT GUESSWORK!
!
        gnorm = dmax1(1.D-3,gnorm)
        step = Min (0.2d0, Min (swap, 0.5d0, 6.d0/gnorm, dell, stepmx*step0)/step0) * step0
        swap = swap + 1.D0
        dell = dell + 0.1D0
        write (iw, '(''            BAR SHORTENED BY'',F12.3,'' PERCENT'')') step/step0*100.D0
        step0 = step0 - step
        if (step0 < 0.1D0) exit
        step = step0
        if (lef) then
          call ef (xparam, escf)
          if (moperr) return
        else
          call flepo (xparam, nvar, escf)
          if (moperr) return
        end if
        if (linear == 0 .and. allocated(pa)) then
          linear = mpack
          pastor = pa
          pbstor = pb
        end if
        i = 0
        do j = 1, numat
          do k = 1,3
            i = i + 1
            xparam(i) = geo(k,j)
          end do
        end do
        if (iflag == 1) then
          write (iw, '(2/10X,''FOR POINT'',I5,'' SECOND STRUCTURE'')') iloop
        else
          write (iw, '(2/10X,''FOR POINT'',I5,'' FIRST  STRUCTURE'')') iloop
        end if
        write (iw, '('' DISTANCE A - B  '',F12.6)') step
        sum = step
!
!   NOW TO CALCULATE THE "CORRECT" GRADIENTS, SWITCH OFF 'STEP'.
!
        step = 0.D0
        grad(:nvar) = grold(:nvar)
        call compfg (xparam, .TRUE., funct1, .TRUE., grad, .TRUE.)
        if (moperr) return
        grold(:nvar) = grad(:nvar)
        if (gradnt) then
          write (iw, '(''  ACTUAL GRADIENTS OF THIS POINT'')')
          write (iw, '(8F10.4)') (grad(i),i=1,nvar)
        end if
        write (iw, '('' HEAT            '',F12.6)') funct1
        gnorm = dsqrt(ddot(nvar,grad,1,grad,1))
        write (iw, '('' GRADIENT NORM   '',F12.6)') gnorm
        cosine = cosine*one
        write (iw, '('' DIRECTION COSINE'',F12.6)') cosine
        if (mod(iloop - 1, 50) == 0) &
        call to_screen(" Structure  No.  Distance from A - B  Heat of Formation Gradient Norm  Cosine of angle")
        if (iflag == 1) then
          write (line, '( "  SECOND ",i5,2f18.6,3f16.6)') iloop, sum,  funct1, gnorm, cosine
        else
          write (line, '( "  FIRST  ",i5,2f18.6,3f16.6)') iloop, sum, funct1, gnorm, cosine
        end if
        call to_screen(trim(line))
        if (intl) then
          xyz = geo(:,:numat)
          call xyzint (xyz, numat, na, nb, nc, 1.d0, geo)
          k = 0
          do i = 1, numat
            do j = 1, min(i - 1,3)
              k = k + 1
              loc(1,k) = i
              loc(2,k) = j
            end do
          end do
        end if
        call geout (iw)
        if (intl) then
          na = 0
          geo(:,:numat) = xyz
          k = 0
          do i = 1, numat
            do j = 1, 3
              k = k + 1
              loc(1,k) = i
              loc(2,k) = j
            end do
          end do
        end if
        if (cosine < 0.d0) then
         i = 0
        end if
        if (swap > 2.9D0 .or. iloop > 3 .and. cosine < 0.D0 .or. escf > eold) then
          if (swap > 2.9D0) then
            swap = 0.D0
          else
            swap = 0.5D0
          end if
!
!   SWAP REACTANT AND PRODUCT AROUND
!
          finish = (gok(1) .and. gok(2) .and. cosine < -0.8D0)
          if (finish) then
            call mopend("BOTH SYSTEMS ARE ON THE SAME SIDE OF THE TRANSITION STATE.")
            return
          end if
          time2 = seconds(2)
          write (iw, '('' TIME='',F9.2)') time2 - time1
          time1 = time2
          write (iw, '(''  REACTANTS AND PRODUCTS SWAPPED AROUND'')')
          iflag = 1 - iflag
          one = -1.D0
          eold = escf
          if (gnorm > 10.D0) gok(iflag + 1) = .TRUE.
          gnorm = sum
          do i = 1, natoms
            do j = 1, 3
              x = geo(j,i)
              geo(j,i) = geoa(j,i)
              geoa(j,i) = x
            end do
          end do
          do i = 1, nvar
            x = xold(i)
            xold(i) = xparam(i)
            xparam(i) = x
          end do
!
!
!    SWAP AROUND THE DENSITY MATRICES.
!
          do i = 1, linear
            x = pastor(i)
            pastor(i) = pa(i)
            pa(i) = x
            x = pbstor(i)
            pbstor(i) = pb(i)
            pb(i) = x
            p(i) = pa(i) + pb(i)
          end do
          if (finish) exit
        else
          one = 1.D0
        end if
      end do
      iflepo = 18
      write (iw, '('' AT END OF STEPWISE ASCENT OF REACTION BARRIER'')')
      gold = dsqrt(ddot(nvar,grad,1,grad,1))
      call compfg (xparam, .TRUE., funct1, .TRUE., grad, .TRUE.)
      if (moperr) return
      gnorm = dsqrt(ddot(nvar,grad,1,grad,1))
      grold(:nvar) = xparam(:nvar)
      call writmo ()
      if (moperr) return
!
! THE GEOMETRIES HAVE (A) BEEN OPTIMIZED CORRECTLY, OR
!                     (B) BOTH ENDED UP ON THE SAME SIDE OF THE T.S.
!
!  TRANSITION STATE LIES BETWEEN THE TWO GEOMETRIES
!
      c1 = gold/(gold + gnorm)
      c2 = 1.D0 - c1

      xparam(:nvar) = c1*grold(:nvar) + c2*xold(:nvar)
      step = 0.D0
      sum = gnorm
      grold(:nvar) = grad(:nvar)
      call compfg (xparam, .TRUE., funct1, .TRUE., grad, .TRUE.)
      if (moperr) return
      gnorm = dsqrt(ddot(nvar,grad,1,grad,1))
      if (gnorm < sum) then
!
!  Gradient has been minimized - write out "better" ts
!
        write (iw, '('' BEST ESTIMATE GEOMETRY OF THE TRANSITION STATE'')')
        write (iw, '(2/10X,'' C1='',F8.3,5X,''C2='',F8.3)') c1, c2
        escf = funct1
        call writmo ()
      else
!
!  Gradient is worse - restore last point
!
        grad = gold
        xparam(:nvar) = grold
      end if
      iflepo = -1 ! Prevent any more calls to writmo
      if (moperr) return
      end subroutine react1
      subroutine dock (geoa, geo, dist)
   !***********************************************************************
   !
   !  DOCK rotates and translates the geometry in GEOA so that the root
   !       mean square difference between GEOA and GEO is a minimum.
   !
   !  On input  GEOA = Cartesian coordinate geometry of first system.
   !            GEO  = Cartesian coordinate geometry of second system.
   !            NUMAT= Number of atoms.
   !
   ! On output  GEOA = "best fit" Cartesian geometry of first system.
   !            GEO  = Unchanged
   !            DIST = Scalar of distance between the two systems.
   !
   !***********************************************************************
    use molkst_C, only : numat, id
    implicit none
    double precision, dimension (3, numat), intent (inout) :: geoa
    double precision, dimension (3, numat), intent (inout) :: geo
    double precision, intent (out) :: dist
    integer :: i, ir, j, jr, k, l
    double precision :: ca, sa, sum, summ, sum_b(3), x
    double precision, save :: sum_a(3) = 0.d0
    integer, dimension (2, 3), save :: irot
    logical, save :: first = .true.
    intrinsic Abs, Cos, Sqrt
    data irot / 1, 2, 1, 3, 2, 3 /
    if (first) then
      first = .false.
      do j = 1, numat - id
        sum_a(:) = sum_a(:) + geo(:,j)
      end do
      sum_a = sum_a/numat
    end if
    do j = 1, numat - id
        geo(:,j) = geo(:,j) - sum_a(:)
    end do
    sum_b = 0.d0
    do j = 1, numat - id
      sum_b(:) = sum_b(:) + geoa(:, j)
    end do
    sum = 0.d0
    sum_b = sum_b / numat
    do j = 1, numat - id
      geoa(:, j) = geoa(:, j) - sum_b(:)
      sum = sum + (geo(1, j)-geoa(1, j)) ** 2 + (geo(2, j)-geoa(2, j)) ** 2 + &
     & (geo(3, j)-geoa(3, j)) ** 2
    end do
    do l = 3, -3, -1
      !
      !     DOCKING IS DONE IN STEPS OF 16, 4, AND 1 DEGREES AT A TIME.
      !
      ca = Cos (4.d0**(l-1)*0.0174532925199432957d0)
      sa = Sqrt (Abs (1.d0-ca**2))
      do j = 1, 3
        ir = irot(1, j)
        jr = irot(2, j)
        do i = 1, 40
          summ = 0.d0
          do k = 1, numat
            x = ca * geoa(ir, k) + sa * geoa(jr, k)
            geoa(jr, k) = -sa * geoa(ir, k) + ca * geoa(jr, k)
            geoa(ir, k) = x
            summ = summ + (geo(1, k)-geoa(1, k)) ** 2 + &
                          (geo(2, k)-geoa(2, k)) ** 2 + &
                          (geo(3, k)-geoa(3, k)) ** 2
          end do
          if (summ > sum) then
            if (i > 1) go to 1000
            sa = -sa
          else
            sum = summ
          end if
          if (i == 1) sum = summ
        end do
        cycle
1000    continue
        sa = -sa
        do k = 1, numat
          x = ca * geoa(ir, k) + sa * geoa(jr, k)
          geoa(jr, k) = -sa * geoa(ir, k) + ca * geoa(jr, k)
          geoa(ir, k) = x
        end do
      end do
    end do
    do j = 1, numat - id
      geo(:,j) = geo(:,j) + sum_a(:)
      geoa(:,j) = geoa(:,j) + sum_a(:)
    end do
    dist = sqrt(sum)
    return
end subroutine dock
