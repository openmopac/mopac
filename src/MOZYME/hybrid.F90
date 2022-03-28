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

subroutine hybrid (catom)
!
!   Construct hybrid atomic orbitals for each atom.
!
!   The hybrid orbitals point to the nearest atoms.
!
    use molkst_C, only: numat, norbs, keywrd
    use chanel_C, only: iw
    use common_arrays_C, only : nbonds, ibonds, f, nat
    use MOZYME_C, only : morb, iorbs
    use elemts_C, only: atom_names
    implicit none
    double precision, dimension (morb, norbs), intent (out) :: catom
    logical :: debug
    integer :: i, ii, j, jj, jl, k, l, loop, m, nb
    double precision :: root2, rt2
    integer, dimension (16) :: nf_loc, nl_loc
    double precision, dimension (10) :: eig
    double precision, dimension (190) :: a
    double precision, dimension (145) :: c
    character :: element*13
    integer, external :: ijbo
    data nf_loc / 1, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 /
    data nl_loc / 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 /
   !
   ! ... Executable Statements ...
   !
    debug = (Index (keywrd, " HYBRID") /= 0)
    root2 = Sqrt (2.d0)
    loop = 0
    do ii = 1, numat
      nb = nbonds(ii)
      if (nb > 15) then
        go to 1000
      else if (iorbs(ii) /= 0) then
        if (iorbs(ii) == 1) then
            !
            !  Hydrogen: single orbital = hybrid orbital
            !
          loop = loop + 1
          catom(1, loop) = 1.d0
        else
            !
            !  Heavy atom, therefore make hydrids.
            !
          k = ((nb+4)*(nb+5)) / 2
          do j = 1, k
            a(j) = 0.d0
          end do
          k = 10
          do j = 1, nb
            jj = ibonds(j, ii)
            jl = ijbo (jj, ii)
            m = iorbs(jj)
            if (ibonds(j, ii) < ii) then
              do l = 1, 4
                a(k+l) = f(jl-m+1+m*l)
              end do
            else
              do l = 1, 4
                a(k+l) = f(jl+l)
              end do
            end if
            k = k + 4 + j
          end do
          k = 4 + nb
            !
            !   Perturb A so that any symmetry that exists is destroyed.
            !
          do i = 1, 4
            do l = 1, 4
              j = (i*(i-1)) / 2 + l
              a(j) = a(j) + j * 1.d-8
            end do
          end do
          do i = 5, k
            do l = 1, 4
              j = (i*(i-1)) / 2 + l
              a(j) = a(j) + j * 2.d-4
            end do
          end do
          if (debug) then
            write (iw, "(A,I4)") " ATOM NUMBER", ii
            do i = 1, k
              write (iw, "(8F10.4)") (a(j), j=(i*(i-1))/2+1, (i*(i+1))/2)
            end do
            write (iw,*)
          end if
          if (k == 1) then
            eig(1) = a(1)
            c(1) = 1.d0
          else
            call rsp (a, k, eig, c)
               !
               !  Set phase of coefficient 1 in each M.O. to +1
               !
            do i = 1, k
              if (c((i-1)*k+1) < 1.d-14) then
                do j = 1, k
                  c((i-1)*k+j) = -c((i-1)*k+j)
                end do
              end if
            end do
          end if
          if (debug) then
            write (iw, "(8F10.4)") (eig(j), j=1, k)
            write (iw,*)
            do i = 1, 4
              write (iw, "(8F10.4)") (c(j), j=(i-1)*k+1, i*k)
            end do
          end if
          j = 8 - k
          if (j > 1) then
            call minloc (c, k, j)
          end if
          call local2 (c, k, 4, nf_loc, nl_loc, nb+1)
          if (debug) then
            write (iw,*) "Hybrid orbitals"
            write (iw,*)
            do i = 1, 4
              write (iw, "(8F10.4)") (c(j), j=(i-1)*k+1, i*k)
            end do
          end if
          rt2 = root2
          do i = 1, 4
            rt2 = 0.d0
            do j = 1, 4
              rt2 = rt2 + c((i-1)*k+j) ** 2
            end do
            rt2 = 1.d0 / Sqrt (rt2)
            do j = 1, 4
              catom(j, i+loop) = c((i-1)*k+j) * rt2
            end do
          end do
            !
            !  Make "d" orbitals the simple orthogonal set
            !
          k = Min (4, iorbs(ii))
          do i = 1, k
            do j = 5, iorbs(ii)
              catom(j, i+loop) = 0.d0
            end do
          end do
          do i = 5, iorbs(ii)
            do j = 1, iorbs(ii)
              catom(j, i+loop) = 0.d0
            end do
            catom(i, i+loop) = 1.d0
          end do
          if (debug) then
            write (iw,*)
            do i = 1, 4
              write (iw, "(20X,4F12.5)") (catom(j, i+loop), j=1, 4)
            end do
          end if
          loop = loop + iorbs(ii)
        end if
      end if
    end do
    return
1000 element = atom_names(nat(ii))
    do
      if(element(1:1) /= " ") exit
      element(:12) = element(2:)
    end do
    write (iw, "(A,I4,A,I8)") " For atom", ii, ", a "//trim(element)//", number of bonds:", nb
    stop
end subroutine hybrid
subroutine minloc (vecs, nvec, n)
   !***********************************************************************
   !                                                                      *
   ! MINLOC puts the hybrid orbitals that are ill-defined into a          *
   ! definite form.  If a heavy atom bonds to one other atom, then only   *
   ! the orientation of the sigma bond is defined, the orientations of    *
   ! the two or three other hybrids is not well defined (is ill-defined). *
   !                                                                      *
   ! By using a simple unitary transform, the orientations of the two or  *
   ! ill-defined hybrids becomes well-defined.  Note:  Although they are  *
   ! now well-defined, they do not need to make chemical sense.  This     *
   ! operation is performed to allow MOZYME to be more easily ported.     *
   ! As with the original hybrids, the new hybrids, although              *
   ! well-defined, are not rotationally invariant.  Rotational            *
   ! invariance is only achieved at self-consistency.                     *
   !***********************************************************************
   !
   !.. Implicit Declarations ..
    implicit none
   !
   !.. Formal Arguments ..
    integer, intent (in) :: nvec
    double precision, dimension (nvec, 4), intent (inout) :: vecs
    integer, intent (in) :: n
   !
   !.. Local Scalars ..
    integer :: i, j, l1, l2
    double precision :: beta
    double precision :: rot = 0.999d0
    double precision :: alpha, sum
   !
   !.. Intrinsic Functions ..
    intrinsic Sqrt
   !
   ! ... Executable Statements ...
   !
    if (n /= 2) then
      !
      !                             x   x   x   x  x
      !   Pattern of Coefficients:  x   x   x   x  x
      !                             x   x   x   x  x
      !                             x   x   x   x  x
      !
      !   Re-hybridize three orbitals.  First, annihilate element VECS(2,3)
      !
      do i = 2, 4
        sum = vecs(i, 2) ** 2 + vecs(i, 3) ** 2
        if (sum > 0.1d0) go to 1000
      end do
1010  sum = vecs(i, 4) ** 2 + vecs(i, 2) ** 2
      !
      !                             x   x   x   x  x
      !   Pattern of Coefficients:  x   x   x   x  x
      !                             x   0   x   x  x
      !                             x   x   x   x  x
      !
      !                                  Next, annihilate element VECS(2,4)
      !
      sum = 1.d0 / Sqrt (sum)
      alpha = vecs(i, 4) * sum
      beta = vecs(i, 2) * sum
      do j = 1, nvec
        sum = alpha * vecs(j, 4) + beta * vecs(j, 2)
        vecs(j, 4) = -beta * vecs(j, 4) + alpha * vecs(j, 2)
        vecs(j, 2) = sum
      end do
      go to 1020
1000  sum = 1.d0 / Sqrt (sum)
      alpha = vecs(i, 2) * sum
      beta = vecs(i, 3) * sum
      do j = 1, nvec
        sum = alpha * vecs(j, 2) + beta * vecs(j, 3)
        vecs(j, 3) = -beta * vecs(j, 2) + alpha * vecs(j, 3)
        vecs(j, 2) = sum
      end do
      go to 1010
    end if
   !
   !                             x   x   x   x  x
   !   Pattern of Coefficients:  x   x   x   x  x
   !                             x   0   x   x  x
   !                             x   0   x   x  x
   !
1020 do i = 2, 4
           !
           !                                  Annihilate element VECS(i,4)
           !
      sum = vecs(i, 4) ** 2 + vecs(i, 3) ** 2
      if (sum > 0.1d0) go to 1030
    end do
    return
1030 sum = 1.d0 / Sqrt (sum)
    alpha = vecs(i, 4) * sum
    beta = vecs(i, 3) * sum
    do j = 1, nvec
      sum = alpha * vecs(j, 4) + beta * vecs(j, 3)
      vecs(j, 4) = -beta * vecs(j, 4) + alpha * vecs(j, 3)
      vecs(j, 3) = sum
    end do
   !
   !   Finally, make a small perturbation to destroy any accidental
   !   symmetry that might exist.  This symmetry can prevent an SCF
   !   from forming.
   !
    do l1 = 1, 4
      do l2 = l1 + 1, 4
        alpha = rot
        beta = Sqrt (1.d0-alpha**2)
        do j = 1, nvec
          sum = alpha * vecs(j, l1) + beta * vecs(j, l2)
          vecs(j, l1) = -beta * vecs(j, l1) + alpha * vecs(j, l2)
          vecs(j, l2) = sum
        end do
      end do
    end do
!
!                             x   x   x   x  x
!   Pattern of Coefficients:  x   x   x   x  x
!                             x   0   x   x  x
!                             x   0   0   x  x
!
end subroutine minloc
