subroutine resetp (mode, loop)
    use molkst_C, only : mpack, keywrd
    use param_global_C, only : pas, pbs
    use common_arrays_C, only : p, pa, pb
!--------------------------------------------------------------------
    implicit none
    integer, intent (in) :: loop, mode
!--------------------------------------------------------------------
    integer :: i
    integer, save :: ilin
    if (loop == 1) then
      ilin = 0
    end if
    if (mode == 0) then
    !
    !  Restore density matrix from store (PAS and PBS) to current array (P)
    !
      do i = 1, mpack
        pa(i) = pas(i+ilin)
        pb(i) = pbs(i+ilin)
        p(i) = pa(i) + pb(i)
      end do
    else
    !
    !  Store density matrix
    !
      if(index(keywrd," UHF") /=0)then
        do i = 1, mpack
          pas(i+ilin) = pa(i)
          pbs(i+ilin) = pb(i)
        end do
      else
       do i = 1, mpack
          pas(i+ilin) = pa(i)
          pbs(i+ilin) = pa(i)
        end do
      end if
    end if
    ilin = ilin + mpack
end subroutine resetp
