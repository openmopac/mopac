subroutine prtlmo ()
    use molkst_C, only: norbs
!
    use MOZYME_C, only: isort, nncf, icocc, ncocc, &
       & cocc, ncf, nnce, icvir, ncvir, cvir, nce, noccupied, nvirtual
!
    use common_arrays_C, only : eigs
!
    use chanel_C, only: iw
    implicit none
    integer :: i
    if (.not. Allocated (isort)) then
      allocate (isort(norbs), stat=i)
      if (i /= 0) then
        call memory_error ("prtlmo")
        return
      end if
      isort = 0
    end if
     write (iw,*) " LOCALIZED MOLECULAR ORBITALS"
     write (iw, "(/,2A,/)") " (PSI**2 IS 10000 TIMES THE ", &
   & "SQUARE OF THE WAVE-FUNCTION ON AN ATOM)"
   !
   !                Print the occupied set of LMOs
   !
    write (iw,'(/10x,a,/)') " Occupied Set"
    write (iw, "(2A)") "  LMO   ENERGY    ATOM PSI**2  ATOM", &
   & " PSI**2  ATOM PSI**2  ATOM PSI**2  ATOM PSI**2"
    call prtlmn (nncf, icocc, ncocc, cocc, ncf, isort, eigs, noccupied, 0)
   !
   !                Print the virtual set of LMOs
   !
    write (iw,'(/10x,a,/)') " Virtual Set"
    write (iw, "(2A)") "  LMO   ENERGY    ATOM PSI**2  ATOM", &
   & " PSI**2  ATOM PSI**2  ATOM PSI**2  ATOM PSI**2"
    call prtlmn (nnce, icvir, ncvir, cvir, nce, isort(noccupied + 1:), &
           & eigs(noccupied + 1:), nvirtual, noccupied)   
    return
end subroutine prtlmo
subroutine prtlmn (nncx, icxxx, ncxxx, cxxx, ncx, isort, eigs, mmos, i_offset)
   !***********************************************************************
   !
   !  PRTLMN (PRTLMO) prints the Localized Molecular Orbitals as the
   !         scalar of the atomic contribution for the top few atoms.
   !  The value is multipled by 10,000, and printed as an integer.
   !
   !***********************************************************************
    use MOZYME_C, only: iorbs, noccupied
!
    use molkst_C, only: norbs
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
    integer :: i, ii, iunsrt, j, jl, ju, k, kk, l, lj, loop, m, n, nj
    double precision :: const, sum, eig_min
    integer, dimension (100) :: ibig, jbig
    double precision, dimension (100) :: xbig
    double precision, dimension(:), allocatable :: w, eigs_temp
    integer, dimension(:), allocatable :: iscrch, jat
    
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
        do j = loop, k
          if (xbig(j) > sum) then
            sum = xbig(j)
            l = ibig(j)
            lj = j
          end if
        end do
        if (sum < const) exit
        xbig(lj) = xbig(loop)
        ibig(lj) = ibig(loop)
        jbig(loop) = Int (sum/const)
        ibig(loop) = l
      end do
      write (iw, "(I5,F9.4,50(5(I5,1X,A2,I5),/14X))") iunsrt + i_offset, eigs (i), &
     & (ibig(j), elemnt(nat(ibig(j))), jbig(j), j=1, loop-1)
    end do
    deallocate (w, iscrch, jat)
    return
end subroutine prtlmn
