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

subroutine mbonds (locc, lvir, f, catom, nfirst, ii, jj, ui, uj, lok, iorbs, &
     & cocc, cvir, cocc_dim, cvir_dim, numat, norbs, morb, mpack1)
    use chanel_C, only: iw
    use molkst_C, only: nelecs
    implicit none
    logical, intent (out) :: lok
    integer, intent (in) :: ii, jj, locc, lvir, morb, mpack1, cocc_dim, &
         & cvir_dim, norbs, numat
    logical, dimension (*), intent (inout) :: ui, uj
    integer, dimension (numat), intent (in) :: iorbs, nfirst
    double precision, dimension (mpack1), intent (in) :: f
    double precision, dimension (cocc_dim), intent (inout) :: cocc
    double precision, dimension (cvir_dim), intent (inout) :: cvir
    double precision, dimension (morb, norbs), intent (in) :: catom
!
    logical :: ok
    integer :: nocc, nvir
    integer :: i, i1 = 0, i2 = 0, iorbsi, iorbsj, j, j1 = 0, j2 = 0, k, l, m, ni, nj
    double precision :: alpha, beta, d, e, e11, e111, e12, e22, e221, one, &
   & sum, summax, summin
    integer, external :: ijbo
    lok = .false.
   !
   !   Four cases: Both light, one light and one heavy,
   !               one heavy and one light, and both heavy
   !
    iorbsi = iorbs(ii)
    iorbsj = iorbs(jj)
   !
    ni = nfirst(ii) - 1
    nj = nfirst(jj) - 1
    summin = 0.d0
    summax = -1.d5
    ok = .false.
    do i = 1, iorbsi
      if ( .not. ui(i)) then
        do j = 1, iorbsj
          if ( .not. uj(j)) then
            ok = .true.
            sum = 0.d0
            l = ijbo (ii, jj)
            do m = 1, iorbsj
              do k = 1, iorbsi
                l = l + 1
                sum = sum + catom(k, i+ni) * f(l) * catom(m, j+nj)
              end do
            end do
            if (sum < summin) then
              i1 = i
              j1 = j
              summin = sum
            end if
            if (sum > summax) then
              i2 = i
              j2 = j
              summax = sum
            end if
          end if
        end do
      end if
    end do
    if ( .not. ok) return
    lok = .true.
    k = 0
    do i = 1, iorbsi
      if (ui(i)) k = k + 1
    end do
    do i = 1, iorbsj
      if (uj(i)) k = k + 1
    end do
   !
    if ((locc + iorbsi + iorbsj - k) > cocc_dim .or. &
      & (lvir + iorbsi + iorbsj - k) > cvir_dim) then
      nocc = Max (1, nelecs/2)
      nvir = Max (1, norbs - nocc)
      write (iw, "(A,I3,A)") " THIS FAULT CAN PROBABLY BE " // &
     & "CORRECTED BY USE OF KEYWORD 'NLMO=", &
     & Max (Nint (Real (locc + iorbsi + iorbsj - k) / Real (nocc) &
          &     * Real (numat) / Real (norbs)), &
          & Nint (Real (lvir + iorbsi + iorbsj - k) / Real (nvir) &
          &     * Real (numat) / Real (norbs))), "'"
      write (iw,*)
      call mopend ("VALUE OF NLMO IS TOO SMALL")
      return
    end if
   !
   !   'best' hybrid pair is I1 and J1
   !
   !  Now form LCAO for the two hybrid orbitals
   !
    if (summax >-summin) then
      i1 = i2
      j1 = j2
      one = -1.d0
      e12 = -summax
    else
      one = 1.d0
      e12 = summin
    end if
    e11 = 0.d0
    e111 = 0.d0
    l = ijbo (ii, ii)
    do i = 1, iorbsi
      do j = 1, i - 1
        l = l + 1
        e111 = e111 + catom(i, ni+i1) * f(l) * catom(j, ni+i1)
      end do
      l = l + 1
      e11 = e11 + catom(i, ni+i1) * f(l) * catom(i, ni+i1)
    end do
    e11 = e11 + e111 * 2.d0
    l = ijbo (jj, jj)
    e22 = 0.d0
    e221 = 0.d0
    do i = 1, iorbsj
      do j = 1, i - 1
        l = l + 1
        e221 = e221 + catom(i, nj+j1) * f(l) * catom(j, nj+j1)
      end do
      l = l + 1
      e22 = e22 + catom(i, nj+j1) * f(l) * catom(i, nj+j1)
    end do
    e22 = e22 + e221 * 2.d0
    d = e11 - e22
    e = Sign (Sqrt(4.d0*e12*e12+d*d), d)
    alpha = Min (0.866d0, Max (Sqrt(0.5d0*(1.d0+d/e)), 0.5d0))
    beta = -Sign (Sqrt(1.d0-alpha*alpha), e12)
    do i = 1, iorbsi
    if (locc + i == 1225)then
      alpha = alpha
    end if
      cocc(locc+i) = catom(i, ni+i1) * alpha
      cvir(lvir+i) = catom(i, ni+i1) * beta
    end do
    do i = 1, Min(iorbsj, Min(cocc_dim -locc -iorbsi, cvir_dim -lvir -iorbsi))
      cocc(locc+iorbsi+i) = catom(i, nj+j1) * beta * one
      cvir(lvir+iorbsi+i) = -catom(i, nj+j1) * alpha * one
    end do
    ui(i1) = .true.
    uj(j1) = .true.
end subroutine mbonds
