subroutine pinout (mode)
   !**********************************************************************
   !
   !   POUT WRITES ALL INFORMATION REGARDING THE LMO'S TO DISK
   !
   !**********************************************************************
    use molkst_C, only: numat, keywrd, nelecs, norbs, line
    use MOZYME_C, only : icocc_dim, iorbs, &
    nce, ncf, ncvir, nncf, nnce, ncocc, &
       & icvir_dim, cocc_dim, cvir_dim, icocc, icvir, cocc, cvir
    use chanel_C, only: iw, iden, density_fn
    use common_arrays_C, only : nbonds, ibonds

    implicit none
    integer :: mode
!
    logical :: opend, exists
    integer :: i, j, k, l, nocc, nvir
!
    nocc = nelecs / 2
    nvir = norbs - nocc
    inquire (unit=iden, opened=opend)
    if (opend) then
      close (unit=iden, status="KEEP")
    end if
    if (mode == 0) then
      inquire (file=density_fn, exist = exists)
      if ( .not. exists) then
        write (line,"(10x,a)") " FILE " // trim(density_fn) // " IS MISSING"
        write(iw,'(//,a)')trim(line)
        call to_screen(line)
        call mopend ("DENSITY file is missing")
        return
      end if
    end if
    open (unit=iden, file=density_fn, status="UNKNOWN", form="UNFORMATTED")
    rewind (iden)
    if (mode == 1) then
      !
      !   OUTPUT ALL DATA FOR LMO'S
      !
      write (iden, err=1030) (ncf(i), i=1, nocc)
      write (iden) (nce(i), i=1, nvir)
      do i = 1, nocc
        write (iden) (icocc(j), j=nncf(i)+1, nncf(i)+ncf(i))
      end do
      do i = 1, nvir
        write (iden) (icvir(j), j=nnce(i)+1, nnce(i)+nce(i))
      end do
      write (iden) (iorbs(i), i=1, numat)
      write (iden) (nbonds(i), i=1, numat)
      write (iden) ((ibonds(j, i), j=1, 9), i=1, numat)
      !
      !  NOW TO WRITE THE LMOS
      !
      do i = 1, nocc
        l = 0
          do j = nncf(i) + 1, nncf(i) + ncf(i)
            l = l + iorbs(icocc(j))
        end do
        write (iden) (cocc(j), j=ncocc(i)+1, ncocc(i)+l)
      end do
      do i = 1, nvir
        l = 0
        do j = nnce(i) + 1, nnce(i) + nce(i)
          l = l + iorbs(icvir(j))
        end do
        write (iden) (cvir(j), j=ncvir(i)+1, ncvir(i)+l)
      end do
      close (iden)
      go to 1000
    else
    !
    !   READ IN ALL DATA FOR LMO'S
    !
      read (iden, end=1010, err=1010) (ncf(i), i=1, nocc)
      read (iden, end=1020, err=1020) (nce(i), i=1, nvir)
      l = 0
      do i = 1,nocc
        l = l + ncf(i)
      end do
!
!  Check to see if there is enough storage for incoming LMOs.
!  If not, then delete old storage, and create new, larger, storage.
!
!  Use a safety marging of 60%
!
        if ((l*8)/5 > icocc_dim) then
!
!  Delete old memory
!
          deallocate (icocc, cocc)
!
!  Re-allocate more memory
!
          cocc_dim = Nint(cocc_dim*(l*1.6/icocc_dim))
          icocc_dim = (l*8)/5
          allocate (icocc(icocc_dim), cocc(cocc_dim))
        end if
        l = 0
        do i = 1,nvir
          l = l + nce(i)
        end do
        if ((l*8)/5 > icvir_dim) then
!
!  Delete old memory
!
          deallocate (icvir, cvir)
!
!  Re-allocate more memory
!
          cvir_dim = Nint(cvir_dim*(l*1.6/icvir_dim))
          icvir_dim = (l*8)/5
          allocate (icvir(icvir_dim), cvir(cvir_dim))
        end if
      !
      !   COMPRESS INCOMING DATA.  
      !
        j = 0
        do i = 1, nocc
          nncf(i) = j
          j = j + ncf(i)
        end do
        j = 0
        do i = 1, nvir
          nnce(i) = j
          j = j + nce(i)
        end do
        do i = 1, nocc
          read (iden, end=1020, err=1020) (icocc(j), j=nncf(i)+1, &
         & nncf(i)+ncf(i))
        end do
        do i = 1, nvir
          read (iden, end=1020, err=1020) (icvir(j), j=nnce(i)+1, &
         & nnce(i)+nce(i))
        end do
        read (iden, end=1020, err=1020) (iorbs(i), i=1, numat)
        read (iden, end=1020, err=1020) (nbonds(i), i=1, numat)
        read (iden, end=1020, err=1020) ((ibonds(j, i), j=1, 9), i=1, numat)
      !
      !  NOW TO READ THE LMOS
      !
        k = 0
        cocc = 0.d0
        do i = 1, nocc
          ncocc(i) = k
          l = 0
          do j = nncf(i) + 1, nncf(i) + ncf(i)
            l = l + iorbs(icocc(j))
          end do
          k = k + l
          read (iden, end=1020, err=1020) (cocc(j), j=ncocc(i)+1, ncocc(i)+l)
        end do
        k = 0
        cvir = 0.d0
        do i = 1, nvir
          ncvir(i) = k
          l = 0
          do j = nnce(i) + 1, nnce(i) + nce(i)
            l = l + iorbs(icvir(j))
          end do
          k = k + l
          read (iden, end=1020, err=1020) (cvir(j), j=ncvir(i)+1, ncvir(i)+l)
        end do
        close (iden)
        go to 1000
      !
1010    continue
      !
      !
        call mopend (" FILE " // trim(density_fn) // " EXISTS, BUT IS FAULTY")
        return
      !
      !
1020    continue
      !
          call mopend (" OLDENS FILE FOR " // trim(density_fn) // " IS CORRUPT")
          return
      end if
     !
1030  call mopend(" Cannot write density matrix to """ // density_fn//"""")
      return
   !
1000 if (Index (keywrd, " PINOUT") /= 0) then
      call prtlmo ()
    end if
  return
   !
end subroutine pinout
