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

subroutine prtlmo ()
    use molkst_C, only: norbs, keywrd, numat
!
    use MOZYME_C, only: isort, nncf, icocc, ncocc, &
       & cocc, ncf, nnce, icvir, ncvir, cvir, nce, noccupied, nvirtual
!
    use common_arrays_C, only : eigs
!
    use chanel_C, only: iw
    implicit none
    integer :: i
    character :: num_1*1
    if (.not. Allocated (isort)) then
      allocate (isort(norbs), stat=i)
      if (i /= 0) then
        call memory_error ("prtlmo")
        return
      end if
      isort = 0
    end if
    num_1 = char(ichar("1") +int(log10(numat + 0.05)))
    write (iw,'(10x,a)') "LOCALIZED MOLECULAR ORBITALS"
   !
   !                Print the occupied set of LMOs
   !
    write (iw,'(/10x,a,/)') " Occupied Set"
    write (iw, "("//num_1//"x,'NUMBER OF CENTERS  LMO ENERGY     COMPOSITION OF ORBITALS ')")
    write (iw, "("//num_1//"x,34x,'(AS PERCENT OF THE LMO)',/)")
    call prtlmn (nncf, icocc, ncocc, cocc, ncf, isort, eigs, noccupied, 0)
    if (index(keywrd, " ALLVEC") /= 0) then
   !
   !                Print the virtual set of LMOs
   !
      write (iw,'(/10x,a,/)') " Virtual Set (NOT re-localized!)"
      write (iw, "("//num_1//"x,'NUMBER OF CENTERS  LMO ENERGY     COMPOSITION OF ORBITALS ')")
      write (iw, "("//num_1//"x,34x,'(AS PERCENT OF THE LMO)',/)")
      call prtlmn (nnce, icvir, ncvir, cvir, nce, isort(noccupied + 1:), &
           & eigs(noccupied + 1:), nvirtual, noccupied)
    else
      write (iw,'(/10x,a,/)') "To print Virtual Set as well, use keyword ""ALLVEC"""
    end if
    return
end subroutine prtlmo
subroutine prtlmn (nncx, icxxx, ncxxx, cxxx, ncx, isort, eigs, mmos, i_offset)
   !***********************************************************************
   !
   !  PRTLMN (PRTLMO) prints the Localized Molecular Orbitals as the
   !         scalar of the atomic contribution for the top few atoms.
   !  The value is multipled by 100, and printed as a percentage.
   !
   !***********************************************************************
    use MOZYME_C, only: iorbs, noccupied
!
    use molkst_C, only: norbs, numat
!
    use chanel_C, only: iw
!
    use elemts_C, only: elemnt
!
    use common_arrays_C, only : nat
!
    implicit none
    integer, intent (in) :: mmos, i_offset
    integer, dimension (norbs), intent (inout) :: isort
    integer, dimension (*), intent (in) :: ncx, nncx, ncxxx
    integer, dimension (*), intent (in) :: icxxx
    double precision, dimension (norbs), intent (in) :: eigs
    double precision, dimension (*), intent (in) :: cxxx
!
    integer :: i, ii, iunsrt, j, jl, ju, k, kk, l, lj, loop, m, n, nj, ibig(100), jbig(100)
    double precision :: const, sum, eig_min, xbig(100)
    double precision, dimension(:), allocatable :: w, eigs_temp
    integer, dimension(:), allocatable :: iscrch, jat
    character :: num_1*1, num_2*1

    allocate (w(Max(1,noccupied*20)), iscrch(Max(1,noccupied*20)), jat(norbs), &
            stat=i)
    if (i /= 0) then
      call memory_error ("prtlmo")
      return
    end if
!
! Construct map of LMO energy levels
!
    allocate (eigs_temp(mmos))
    eigs_temp(:mmos) = eigs(:mmos)
    do i = 1, mmos
      eig_min = 1.d7
      ii = 0
      do j = 1, mmos
        if (eigs_temp(j) < eig_min) then
          eig_min = eigs_temp(j)
          ii = j
        end if
      end do
      isort(i) = ii
      eigs_temp(ii) = 1.d8
    end do
   !
   !   Put LMO contributions per atom into an atom-based array
   !   Store number of atoms in each LMO in JAT
   !   Store number of times an atom is used in the LMOs in IAT
   !   Store the address of each atom in the square array IM
   !
    j = 0
    const = 2.d-4
    do iunsrt = 1, mmos
      i = isort(iunsrt)
      n = ncxxx(i)
      k = 0
      do ii = nncx(i) + 1, nncx(i) + ncx(i)
        l = icxxx(ii)
        nj = iorbs(l)
        sum = 0.d0
        do m = 1, nj
          sum = sum + cxxx(m+n) ** 2
        end do
        if (sum > const) then
          j = j + 1
          k = k + 1
          iscrch(j) = l
          w(j) = sum * 2.d0
        end if
        n = n + nj
      end do
      !
      !   JAT holds the number of atoms that are important in the LMO(I)
      !
      jat(i) = k
    end do
    ju = 0
    num_1 = char(ichar("3") +int(log10(numat + 0.05)))
    num_2 = char(ichar("2") +int(log10(norbs + 0.05)))
    do iunsrt = 1, mmos
      i = isort(iunsrt)
      jl = ju + 1
      ju = jl + jat(i) - 1
      k = 0
      do j = jl, ju
        k = k + 1
        xbig(k) = w(j)
        ibig(k) = iscrch(j)
      end do
      kk = Min (20, k)
      do loop = 1, kk
        sum = -1.d0
        l = 0
        lj = 0
        do j = loop, k
          if (xbig(j) > sum) then
            sum = xbig(j)
            l = ibig(j)
            lj = j
          end if
        end do
        if (sum < const*100.d0) exit
        xbig(lj) = xbig(loop)
        ibig(lj) = ibig(loop)
        jbig(loop) = Int (sum/const)
        ibig(loop) = l
      end do
      loop = loop - 1
      sum = 0.d0
      do j = 1, loop
        sum = sum + jbig(j)**2
      end do
      k = iunsrt + i_offset
      sum = 1.d8/sum
      if (loop == 1) then
        if (jbig(1)*0.01d0> 99.949d0) then
          write (iw, '(i'//num_1//',f10.4,f17.5, 3x,a2,i'//num_2//',f6.1)') &
                & k, sum,eigs(i), elemnt(nat(ibig(1))),ibig(1),jbig(1)*0.01d0
        else
          write (iw, '(i'//num_1//',f10.4,f17.5, 3x,a2,i'//num_2//',f6.2)') &
                & k, sum,eigs(i), elemnt(nat(ibig(1))),ibig(1),jbig(1)*0.01d0
        end if
      else
        if (loop < 6) then
          write (iw, '(i'//num_1//',f10.4,f17.5, 4(5(3x,a2,i'//num_2//',f6.2)))') &
          k, sum,eigs(i), (elemnt(nat(ibig(j))),ibig(j),jbig(j)*0.01d0, j=1, loop)
        else
          write (iw, '(i'//num_1//',f10.4,f17.5, 4(5(3x,a2,i'//num_2//',f6.2),/,'//num_1//'x,27x))') &
          k, sum,eigs(i), (elemnt(nat(ibig(j))),ibig(j),jbig(j)*0.01d0, j=1, loop)
        end if
      end if
    end do
    deallocate (w, iscrch, jat)
    return
  end subroutine prtlmn
