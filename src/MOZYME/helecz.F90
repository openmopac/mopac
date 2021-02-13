double precision function helecz ()
!
!   Return the total electronic energy,
!   from e = 0.5P(H + F)
!
    use molkst_C, only: numat
    use MOZYME_C, only : iorbs
    use common_arrays_C, only : p, h, f
    implicit none
    integer :: i, j, k, l, l1, l2, m
    double precision :: ed, ee
    integer, external :: ijbo
    ed = 0.0d00
    ee = 0.0d00
    do i = 1, numat
      do j = 1, i - 1
          k = ijbo (i, j)
          if (k >= 0) then
            l = k + iorbs(i) * iorbs(j)
            do m = k + 1, l
              ee = ee + p(m) * (h(m)+f(m))
            end do
          end if
      end do
        k = ijbo (i, i)
        do l1 = 1, iorbs(i)
          do l2 = 1, l1 - 1
            k = k + 1
            ee = ee + p(k) * (h(k)+f(k))
          end do
          k = k + 1
          ed = ed + p(k) * (h(k)+f(k))
        end do
    end do
    ee = ee + .5d00 * ed
    helecz = ee
end function helecz
