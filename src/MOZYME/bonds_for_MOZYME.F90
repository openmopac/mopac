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

   subroutine bonds_for_MOZYME ()
    use molkst_C, only: numat, keywrd, nl_atoms, maxtxt
    use chanel_C, only: iw
    use common_arrays_C, only : p, nat, l_atom
    use MOZYME_C, only : iorbs
    implicit none
    logical :: lall, mozyme_style
    integer :: i, io, j, jo, k, kk, kl, ku, l, newl, j_lim
    double precision :: sum, sumlim, valenc
    integer, allocatable, dimension (:) :: ibab
    logical, allocatable, dimension (:) :: l_atom_store
    double precision, allocatable, dimension (:) :: bab
    intrinsic Index
    integer, external :: ijbo
  !
    lall = (Index (keywrd, " ALLBOND") /= 0)
    sumlim = 0.01d0                  !  Report all bonds greater than 0.01 not involving hydrogen
    if (lall) sumlim = sumlim*0.1d0  !  Report all bonds greater than 0.001
    mozyme_style = .true.
    if (index(keywrd, " IRC") + index(keywrd, " DRC") + index(keywrd, " MINI") /= 0) mozyme_style = .false.
    if (mozyme_style) then
      if (maxtxt > 25) then
        write (iw, '(1X,2/37X,"(VALENCIES)",1x,"BOND ORDER       TO",/)')
      else
        write (iw, '(1X,2/12X,"(VALENCIES)   BOND ORDERS",/)')
      end if
      allocate (bab(100), ibab(100))
    else
      write (iw, '(1X,2/20X,''BOND ORDERS AND VALENCIES'')')
      sumlim = -1.d0
      i = (nl_atoms*(nl_atoms + 1))/2
      allocate (bab(i), l_atom_store(nl_atoms))
    end if
    newl = 0
    l = 0
    do i = 1, numat
    if ( .not. mozyme_style .and.  .not. l_atom(i)) cycle
      newl = newl + 1
      valenc = 0.d0
      io = iorbs(i)
      if (mozyme_style) then
          j_lim = numat
          l = 0
        else
          j_lim = i - 1
        end if
      if (nat(i) /= 1 .or. lall) then
        if (io == 1) then
          kk = ijbo (i, i) + 1
          valenc = 2.d0 * p(kk) - p(kk) ** 2
        else
          kk = ijbo (i, i)
          do j = 1, iorbs(i)
            do k = 1, j
              kk = kk + 1
              valenc = valenc - 2.d0*p(kk) * p(kk)
            end do
            valenc = valenc + p(kk) * p(kk)
            valenc = valenc + 2.d0 * p(kk)
          end do
        end if
        do j = 1, j_lim
          if ( .not. mozyme_style .and.  .not. l_atom(j)) cycle
          if (nat(j) /= 1 .or. lall) then
            jo = iorbs(j)
            if (i /= j .and. ijbo (i, j) >= 0) then
              kl = ijbo (i, j) + 1
              ku = kl + io * jo - 1
              sum = 0.d0
              do k = kl, ku
                sum = sum + p(k) ** 2
              end do
              if (sum > sumlim) then
                l = l + 1
                bab(l) = sum
                if (mozyme_style) ibab(l) = j
              end if
            else
              if ( .not. mozyme_style) then
                l = l + 1
                bab(l) = 0.d0
              end if
            end if
          end if
        end do
        if (mozyme_style) then
          call print_bonds_compact(i,l, valenc, bab, ibab)
        else
          l = l + 1
          bab(l) = valenc
        end if
      end if
    end do
    if ( .not. mozyme_style) then
      l_atom_store = l_atom(:nl_atoms)
      l_atom(:nl_atoms) = .true.
      call vecprt (bab, newl)
      l_atom(:nl_atoms) = l_atom_store
      deallocate (l_atom_store, bab)
    end if
  end subroutine bonds_for_MOZYME
  subroutine print_bonds_compact(i_atom,n_non_zero, valenc, bab, ibab)
!
! Print bonds in a compact form.
! On input: i_atom = atom number
!           n_non_zero = number of non-zero bonds
!
  use molkst_C, only: maxtxt
  use common_arrays_C, only : txtatm,nat
  use elemts_C, only: elemnt
  use chanel_C, only: iw
  implicit none
  double precision, intent (inout) :: bab(100)
  double precision, intent (in) :: valenc
  integer, intent (inout) :: ibab(100)
  integer, intent (in) :: i_atom, n_non_zero
!
  integer :: j, k, m, mtx, nper, ju
  double precision :: sum
  character :: nc*2, blank*60,  f5p3*4 = "f5.3"
  save :: f5p3
    if (n_non_zero == 0) return
    if (maxtxt /= 0) then
      write(nc,'(i2)')maxtxt + 20
      mtx = maxtxt + 23
      nper = 33
    else
      mtx = 21
      nper = 6
    end if
    blank = " "
    do j = 1, n_non_zero
      m = 0
      sum = 0.d0
      do k = j, n_non_zero
        if (bab(k) > sum) then
          m = k
          sum = bab(k) * (1.d0+1.d-4)
        end if
      end do
      k = ibab(j)
      ibab(j) = ibab(m)
      ibab(m) = k
      sum = bab(j)
      bab(j) = bab(m)
      bab(m) = sum
    end do
    ju = Min (nper, n_non_zero)

    if (maxtxt /= 0) then
      write (iw, "(I6, 1X, A, A, "//f5p3//", A1, F8.3, 2x, a)") i_atom, elemnt (nat(i_atom)) &
      & //"("// txtatm (i_atom)(:maxtxt)//")", "  (", valenc, ")      ", &
      & bab(1), elemnt(nat(ibab(1)))//"("//txtatm(ibab(1))(:maxtxt)//")"
      do j = 2, ju
        write(iw,'('//nc//'x, f8.3, 2x, a)') bab(j), elemnt(nat(ibab(j)))//"("//txtatm(ibab(j))(:maxtxt)//")"
      end do
    else
      write (iw, "(I6,1X,A,A,"//f5p3//",')',6(I6,1x,A2,F6.3))") i_atom, elemnt(nat(i_atom)), &
      "     (", valenc, (ibab(j), elemnt(nat(ibab(j))), &
      & bab(j), j=1, ju)
      do m = nper + 1, n_non_zero, nper
        ju = Min (m+nper-1, n_non_zero)
        write (iw, "(A,6(I6,1x,A2,F6.3))") blank (1:mtx), (ibab(j), &
        & elemnt(nat(ibab(j))), bab(j), j=m, ju)
      end do
    end if
    write (iw,*)
  end subroutine print_bonds_compact
