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

subroutine solir (eigs, cr, nvecs, xyzk, itriv, xyzxyz)
  use common_common, only : numat, iop, xop, id, phonon, iw, iw_new, &
    line, nfirst, per_atom, title, allt1, mop, keywrd, data_set_name, jobnam
  implicit none
  integer, intent(in) :: itriv, nvecs
  complex, dimension(nvecs, nvecs), intent(in) :: cr
  double precision, dimension(nvecs+1), intent(inout) :: eigs
  double precision, dimension(3), intent(in) :: xyzk, xyzxyz
  ! 
  !.. Local Scalars .. 
  !
  logical :: first = .true., prtvec = .false., opend
  integer :: i, i1, iatom, il, iloop, iold, ir, j, jop, k, kk, klim, &
       & l, m, nnnew, nvecs1, nreps, alloc_stat
  complex :: angle, xc
  double precision :: deltak, eig1, eig2, elim, optlim, phalim, rdotk, repp, &
       & twopi, eigs2(120), t1(9,9)
  character :: little(149)*12, num*1, num2*1
  integer :: nclass(iop), itype(nvecs)
  complex :: rotc(iop), vecc(nvecs), phaser(iop, nvecs), refchi(iop, nvecs)
  real :: sum
  complex, dimension(:,:), allocatable :: wfnr
  double precision, external :: reada
  save
  allocate (wfnr(nvecs, nvecs), stat=alloc_stat)
  if (alloc_stat /= 0) then
    write(6,*)"Failed to allocate array wfnr in solir"
    stop
  end if
  inquire(unit=iw_new, opened=opend) 
  if (.not. opend) open (unit = iw_new, file = trim(data_set_name)//".txt", form = "FORMATTED", status = "UNKNOWN", iostat = i)
  phaser = (0.0, 0.0)
  if (first) then
    first = .false.
    twopi = 2.d0 * Acos (-1.d0)
    elim = 0.1d0
    if ( phonon ) elim = 1.0d0
    prtvec = (index(keywrd," PRTVEC") /= 0)
    if (Index (keywrd, "DEGEN") /= 0) then
      elim = reada (keywrd, Index (keywrd, "DEGEN"))
    end if
    if (Index (keywrd, " LET") /= 0) then
      optlim = 10.d0
      phalim = 1.d0
      klim = -1
    else
      if (phonon) then
        optlim = 0.01d0
        phalim = 0.3d0
        klim = -1
      else
        optlim = 0.01d0
        phalim = 0.3d0
        klim = 0
      end if      
    end if
    write(iw_new,'(/10x,a,/)')"BZ Data-set: """//trim(jobnam)//""""
  end if
  !
  !  Evaluate the characters for each operation.
  !
  do jop = 1, iop
    do i = 1, 9
      do j = 1, 9
        t1(i, j) = allt1(i, j, jop)
      end do
    end do
    !
    !  T1 now holds the transform of the atomic orbitals under
    !     operation JOP
    !
    do iatom = 1, numat
      k = iatom
      ir = per_atom(iatom)
      rdotk = mop(2, k, jop)*xyzk(1) + mop(3, k, jop)*xyzk(2) + mop(4, k, jop)*xyzk(3)
      deltak = xyzk(1)*xop(2, jop) + xyzk(2)*xop(3, jop) + xyzk(3)*xop(4, jop)
      k = mop(1, k, jop)
      !
      !  ROTC is the complex phase which arises due to non-trivial
      !       translations in an operation.
      !
      rotc(jop) = Cmplx (Cos(deltak*twopi), -Sin(deltak*twopi))
      il = nfirst(k) - 1
      l = nfirst(iatom) - 1
      !
      !  ANGLE is the complex phase resulting from an atom being moved
      !        into a different unit cell
      !
      angle = Cmplx (Cos(rdotk*twopi), Sin(rdotk*twopi))
      do iloop = 1, nvecs
        do i = 1, nvecs
          vecc(i) = cr(i, iloop)
        end do
        do i = 1, ir
          k = i + il
          xc = 0.0
          do j = 1, ir
            i1 = j + l
            xc = xc + vecc(i1)*real(t1(j, i))
          end do
          wfnr(k, iloop) = xc * angle
        end do
      end do
    end do
    if (prtvec) then
      m = Min (40, nvecs)
      write(iw, "(//,A)")"Complex Eigenvectors, before Operation "//title(jop)// "  Real part"
      write(iw_new, "(//,A)")"Complex Eigenvectors, before Operation "//title(jop)// "  Real part"
      do i = 1, nvecs
        write(iw, "(40f10.4)") (Real (cr(i, j)), j = 1, m)
        write(iw_new, "(40f10.4)") (Real (cr(i, j)), j = 1, m)
      end do
      write(iw, "(/,A)") "Complex Eigenvectors, before Operation "//title(jop)//"  Imaginary part"
      write(iw_new, "(/,A)") "Complex Eigenvectors, before Operation "//title(jop)//"  Imaginary part"
      do i = 1, nvecs
        write(iw, "(40f10.4)") (Aimag (cr(i, j)), j = 1, m)
        write(iw_new, "(40f10.4)") (Aimag (cr(i, j)), j = 1, m)
      end do
      write(iw, "(//,A)")"Complex Eigenvectors, After Operation "//title(jop)//"  Real part"
      write(iw_new, "(//,A)")"Complex Eigenvectors, After Operation "//title(jop)//"  Real part"
      do i = 1, nvecs
        write(iw, "(40f10.4)") (Real (wfnr(i, j)), j = 1, m)
        write(iw_new, "(40f10.4)") (Real (wfnr(i, j)), j = 1, m)
      end do
      write(iw, "(/,A)")"Complex Eigenvectors, After Operation "//title(jop)//"  Imaginary part"
      write(iw_new, "(/,A)")"Complex Eigenvectors, After Operation "//title(jop)//"  Imaginary part"
      do i = 1, nvecs
        write(iw, "(40f10.4)") (Aimag (wfnr(i, j)), j = 1, m)
        write(iw_new, "(40f10.4)") (Aimag (wfnr(i, j)), j = 1, m)
      end do
    end if
    do iloop = 1, nvecs
      xc = 0.0
      do i = 1, nvecs
        xc = xc + cr(i, iloop) * Conjg (wfnr(i, iloop))
      end do
      phaser(jop, iloop) = xc
    end do
  end do

  deallocate (wfnr, stat=alloc_stat)
  if (alloc_stat /= 0) then
    write(6,*)"Array wfnr in solir corrupted"
    stop
  end if

  if (itriv == 3) then
    do j = 1, iop
      if (Abs (xop(2, j)) + Abs (xop(3, j)) + Abs (xop(4, j)) > 0.1d0) then
        xc = phaser(j, 1)
        do i = 1, nvecs
          phaser(j, i) = phaser(j, i) * Conjg (rotc(j))
        end do
      end if
    end do
  elseif (itriv == 2) then
    do j = 1, iop
      if (Abs (xop(2, j))+Abs (xop(3, j))+Abs (xop(4, j)) > 0.1d0) then
        xc = phaser(j, 1)
        do i = 1, nvecs
          phaser(j, i) = phaser(j, i) * Conjg (xc)
        end do
      end if
    end do
  end if
  !***********************************************************************
  !
  !    Start of Tidy-Up
  !
  !   Compress levels so that degenerate levels are printed once only.
  !
  eig1 = eigs(1)
  nvecs1 = nvecs + 1
  eigs(nvecs1) = 1.d10
  do i = 1, nvecs
    do j = 1, iop
      refchi(j, i) = 0.0
    end do
  end do
  j = 0
  iold = 1
  do i = 2, nvecs1
    eig2 = eigs(i)
    if (eig2-eig1 >= elim) then
      j = j + 1
      eigs2(j) = eigs(iold)
      do l = iold, i-1
        do k = 1, iop
          refchi(k, j) = refchi(k, j) + phaser(k, l)
        end do
      end do
      iold = i
    end if
    eig1 = eig2
  end do
  nnnew = j
  do i = 1, nnnew
    do j = 1, iop
      phaser(j, i) = refchi(j, i)
    end do
  end do
  !
  !    Eliminate all operations which are not operations of the group
  !
  l = 0
  do i = 1, iop
    repp = 0
    xc = 0.0
    k = 0
    do j = 1, nnnew
      xc = phaser(i, j) * Conjg (phaser(i, j))
      kk = Nint (Real (xc))
      k = k + kk
      if (Abs (kk - xc) < optlim .and. kk /= 0 .or. Abs (kk - xc) < 0.1d0*optlim) then
        repp = repp + 1
      end if
    end do
    if (repp/nnnew > 0.6d0 .and. k /= klim) then
      l = l + 1
      little(l) = title(i)
      rotc(l) = rotc(i)
      do j = 1, nnnew
        phaser(l, j) = phaser(i, j)
      end do
    end if
  end do
  !
  !  Eliminate Duplicate Operations
  !
  do i = 1, l
    nclass(i) = 1
  end do
  jop = 0
  do i = 1, l
    if (i /= 1) then
      loop1: do j = 1, jop
        do k = 1, nnnew
          if (Abs (phaser(i, k) - phaser(j, k)) > phalim) cycle loop1
        end do
        if (little(j) == little(i)) goto 1000
      end do loop1
      goto 1100
1000  continue
      nclass(j) = nclass(j) + 1
      cycle
    end if
    !
    !   Characters for I are different to those of all known operations.
    !
1100 continue
    jop = jop + 1
    little(jop) = little(i)
    rotc(jop) = rotc(i)
    do k = 1, nnnew
      phaser(jop, k) = phaser(i, k)
    end do
  end do
  !
  !   Label Irreducible Representations
  !
  nreps = 0
  do i = 1, nnnew
    loop2: do j = 1, nreps
      do k = 1, jop
        if (Abs (phaser(k, i) - refchi(k, j)) > 1.d-1) cycle loop2
      end do
      goto 1200
    end do loop2
    nreps = nreps + 1
    l = nreps
    do k = 1, jop
      refchi(k, nreps) = phaser(k, i)
    end do
    goto 1300
    !
    !   Level I has same I.R. as known representation J
    !
1200 continue
    l = j
1300 continue
    itype(i) = l
  end do
  if (line /= " ") then
    line = " "//trim(line)//" = (" 
  else
    line = " = ("
  end if
  num = char(id + ichar("0"))
  sum = 0.0
  do i = 1, id
    sum = max(sum, real(abs(xyzxyz(i))))
  end do
  if (id == 1) then
    num2 =char(int(log10(sum)) + ichar("7"))
    write(iw, '(///10x, "EXPECTATION VALUE OF OPERATORS FOR POINT", a, '//num//'f'//num2//'.4,a)') &
    trim(line), xyzxyz(1),")"
  else
    num2 =char(int(log10(sum)) + ichar("8"))
    write(iw, '(///10x, "EXPECTATION VALUE OF OPERATORS FOR POINT", a, '//num//'f'//num2//'.4,a)') &
    trim(line), (xyzxyz(i), i = 1, id),")"
  end if
  write(iw_new, '(///10x, "EXPECTATION VALUE OF OPERATORS FOR POINT",a, '//num//'f'//num2//'.4,a)') &
  trim(line), (xyzxyz(i), i = 1, id),")"
  write(iw, *)
  write(iw, "(18X,50(I4,1X,A12))") (nclass(i), little(i), i = 1, jop)
  write(iw_new, '(a)')
  write(iw_new, "(18X,50(I4,1X,A12))") (nclass(i), little(i), i = 1, jop)
  k = 1
  do i = 1, nnnew
    do j = 1, jop
      if (real(phaser(j, i)) < 0.0000 .and. real(phaser(j, i)) > -0.0005) phaser(j,i) = phaser(j,i) + (0.0005, 0.0)
      if (imag(phaser(j, i)) < 0.0000 .and. imag(phaser(j, i)) > -0.0005) phaser(j,i) = phaser(j,i) + (0.0, 0.00049)
    end do
    write(iw, "(i3,f10.4,i3,50(a2,f6.3,f7.3,a2))") &
            & k, eigs2(i), itype(i), (" (", phaser(j, i), ") ", j = 1, jop)
    write(iw_new, "(I3,F10.4,I3,50(A2,F6.3,F7.3,A2))") &
            & k, eigs2(i), itype(i), (" (", phaser(j, i), ") ", j = 1, jop)
            k = k + nint(real(phaser(1, i)))
  end do
end subroutine solir
