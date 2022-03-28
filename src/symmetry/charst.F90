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

module get_det
  implicit none
contains
  recursive function determinant(matrix, n) result(det)
!
! Evaluate the determinant of the matrix "matrix"
!
! On input,  matrix      = square matrix of size "n"
! On output, det         = value of the determinant of matrix "matrix"
!
    implicit none
    integer, intent (in) :: n
    double precision, intent (in) :: matrix(n,n)
!
! Local variables
!
    double precision :: det
    integer :: i
    double precision :: sum
    double precision, allocatable :: cf(:,:)
    if (n == 0) then
      sum = 1.d0
    else if (n == 1) then
      sum = matrix(1,1)
    else
      sum = 0
      do i = 1, n
        allocate(cf(n - 1, n - 1))
        cf = cofactor(matrix, n, i)
        sum = sum + ((-1)**(i + 1))*matrix(i, 1)*determinant(cf, n - 1)
        deallocate(cf)
      end do
    end if
    det = sum
  end function determinant
!
  function cofactor(matrix, n, mI)
    implicit none
    integer, intent (in) :: n, mI
    double precision, intent (in) :: matrix(n, n)
!
! Local variables
!
    integer ::  i, j, k, l
    double precision, dimension(:,:), allocatable :: cofactor
    allocate(cofactor(n - 1, n - 1))
    l = 0
    k = 1
    do i = 1, n
      if (i /= mI) then
        l = 1
        do j = 2, n
          cofactor(k, l) = matrix(i, j)
          l = l + 1
        end do
        k = k + 1
      end if
    end do
    return
  end function cofactor
end module

  double precision function charst (vects, ntype, istate, ioper, r,nvecs, first)
!-----------------------------------------------
!   M o d u l elem s
!-----------------------------------------------
      USE meci_C, only : nalmat, microa, microb, lab, conf, nmos
      use chanel_C, only : iw
      use symmetry_C, only : elem, jelem
      use molkst_C, only : keywrd, numat, norbs
      use get_det
      implicit none
      integer , intent(in) :: istate
      integer  :: ioper
      integer , intent(in) :: nvecs
      logical  :: first
      integer , intent(in) :: ntype(norbs)
      double precision , intent(in) :: vects(nvecs,nmos)
      double precision  :: r(3,3)
!-----------------------------------------------
!   L o c a l   P a r a m elem t elem r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l elem s
!-----------------------------------------------
      integer , dimension(2,3) :: ip
      integer , dimension(2,5) :: id
      integer , dimension(2,9) :: loc
      integer , dimension(:), allocatable :: iphase
      integer , dimension(:,:), allocatable :: iperma, ipermb
      integer :: nstate, j, i, iloop, iatom, jatom, ibase, kj, icheck, jcheck, &
        ii, jj, k, l, ne, nai, nbi
      double precision, dimension(:,:), allocatable :: vect1, vect2
      double precision :: h(5), p(3), d(5), matrix(36)
      double precision, dimension(:,:), allocatable :: t2, t4
      double precision, dimension(:), allocatable :: work
      double precision :: sum, det, suma, sumb
      logical :: posita, positb, debug
      double precision, external :: ddot
!***********************************************************************
!
!    CHARST evaluates the character of the State ISTATE under the
!           operation IOPER.
!
!   Info:   VECTS  = Molecular Orbitals
!           CONF   = State Eigenvectors
!           NVECS  = Number of atomic bases
!           NMOS   = Number of M.O.s in active space
!           NSTATE    = Number of Microstates in each State
!
!***********************************************************************
      data nstate/ 0/
      data posita, positb/ 2*.FALSE./
      save vect1, vect2, t2, t4, iphase, iperma, ipermb, work
      charst = 1.D0
      if (istate < 0) then
!
!   Reset MICROA and MICROB, if necessary.
!
        if (posita) then
          nalmat(:nstate) = nmos - nalmat(:nstate)
          microa(:nmos,:nstate) = 1 - microa(:nmos,:nstate)
        end if
        if (positb) then
          microb(:nmos,:nstate) = 1 - microb(:nmos,:nstate)
        end if
        charst = 0.D0
        if (allocated(vect1))  deallocate(vect1)
        if (allocated(vect2))  deallocate(vect2)
        if (allocated(t2))     deallocate(t2)
        if (allocated(t4))     deallocate(t4)
        if (allocated(iperma)) deallocate(iperma)
        if (allocated(ipermb)) deallocate(ipermb)
        if (allocated(work))   deallocate(work)
        if (allocated(iphase)) deallocate(iphase)
        return
      end if
      debug = index(keywrd,'CHARST')/=0 .and. index(keywrd,'DEBUG')/=0
!
!  Trivial case:  Operation is 'elem', the identity.
!
      charst = 1.D0
      if (ioper == 1) return
!
!   Non-trivial case
!
      if (istate == 1) then
        if (ioper == 2) then
         i = max(lab, nmos**2)
        allocate(vect1(nvecs,nmos), vect2(nvecs, nmos), t2(nmos, nmos), &
      & t4(lab,lab),iphase(lab), iperma(nmos + 3,lab), ipermb(nmos + 3,lab), &
      & work(i))
      end if
        nstate = lab
!
!    Set up the M.O. unitary matrix, <psi | IOPER | psi>
!
        do iloop = 1, nmos
          do iatom = 1, numat
            jatom = jelem(ioper,iatom)
            ibase = 0
            kj = 0
            do i = 1, nvecs
              icheck = ntype(i)/100
              if (icheck == iatom) then
                ibase = ibase + 1
                loc(1,ibase) = i
              end if
              if (icheck /= jatom) cycle
              kj = kj + 1
              loc(2,kj) = i
            end do
            if (ibase == 0) cycle
!
!   's'-type basis function
!
            icheck = loc(1,1)
            jcheck = loc(2,1)
            vect1(icheck,iloop) = vects(icheck,iloop)
            vect2(jcheck,iloop) = vects(icheck,iloop)
            if (ibase < 4) cycle
!
!    Atom I had a 'p' shell
!
            ip(1,:) = 0
            id(1,:3) = 0
            id(1,4) = 0
            id(1,5) = 0
            do i = 2, ibase
              icheck = loc(1,i)
              if (i <= 4) then
                p(i-1) = vects(icheck,iloop)
                ip(1,i-1) = loc(1,i)
                ip(2,i-1) = loc(2,i)
              else
                d(i-4) = vects(icheck,iloop)
                id(1,i-4) = loc(1,i)
                id(2,i-4) = loc(2,i)
              end if
            end do
            if (ibase /= 1) then
!
!    'p' transform
!
              h(1) = r(1,1)*p(1) + r(2,1)*p(2) + r(3,1)*p(3)
              h(2) = r(1,2)*p(1) + r(2,2)*p(2) + r(3,2)*p(3)
              h(3) = r(1,3)*p(1) + r(2,3)*p(2) + r(3,3)*p(3)
              p(1) = elem(1,1,ioper)*h(1) + elem(1,2,ioper)*h(2) + elem(1,3,ioper)*h(3)
              p(2) = elem(2,1,ioper)*h(1) + elem(2,2,ioper)*h(2) + elem(2,3,ioper)*h(3)
              p(3) = elem(3,1,ioper)*h(1) + elem(3,2,ioper)*h(2) + elem(3,3,ioper)*h(3)
              do i = 1, 3
                if (ip(1,i) < 1) return
                ii = ip(1,i)
                jj = ip(2,i)
                vect1(ii,iloop) = h(i)
                vect2(jj,iloop) = p(i)
              end do
            end if
            if (ibase /= 9) cycle
!
!   'd' transform
!
            h = d
            call dtrans (d, h, ioper, first, r)
            do i = 1, 5
              if (id(1,i) < 1) return
              ii = id(1,i)
              jj = id(2,i)
              vect1(ii,iloop) = h(i)
              vect2(jj,iloop) = d(i)
            end do
          end do
        end do
!
!    VECT1 holds the molecular orbitals in the coordinate system of
!          SYMTRZ.
!    VECT2 holds the same orbitals, after being operated on by the
!          symmetry operation IOPER.
!
!     T2 hold the unitary transform for the set of M.O.s
!        <Vect1|IOPER|VECT2>
!
        do i = 1, nmos
          do j = 1, nmos
            sum = 0.D0
            do k = 1, nvecs
              sum = sum + vect1(k,i)*vect2(k,j)
            end do
            t2(i,j) = sum
          end do
        end do
        if (debug) then
          write (iw, '(A)') ' Symmetry Operation in CHARST'
          write (iw, '(3F12.6)') ((elem(i,j,ioper),i=1,3),j=1,3)
          write (iw, '(A)') ' Transform of M.O.s'
          do i = 1, nmos
            write (iw, '(8F12.6)') (t2(j,i),j=1,nmos)
          end do
        end if
!***********************************************************************
!
!   For the microstates, take the positron equivalent if more than
!   half filled.
!
        if (ioper==2 .and. istate==1) then
          do i = 1, nstate
            l = 0
            do j = 1, nmos
              if (microa(j,i) == 0) cycle
              do k = j, nmos
                l = l + microb(k,i)
              end do
            end do
            iphase(i) = 1 - 2*mod(l,2)
          end do
          k = 0
          do i = 1, nmos
            k = k + microa(i,1)
          end do
!            NE=NE+K
          posita = k > nmos/2
          if (posita) then
            do j = 1, nstate
              nalmat(j) = nmos - nalmat(j)
              l = (2 - iphase(j))/2
              do i = 1, nmos
                if (microa(i,j) == 0) then
                  do k = i + 1, nmos
                    l = l + microa(k,j)
                  end do
                end if
                microa(i,j) = 1 - microa(i,j)
              end do
              iphase(j) = 1 - 2*mod(l,2)
            end do
          end if
          k = 0
          do i = 1, nmos
            k = k + microb(i,1)
          end do
          positb = k > nmos/2
          if (positb) then
            do j = 1, nstate
              l = (2 - iphase(j))/2
              do i = 1, nmos
                if (microb(i,j) == 0) then
                  do k = i + 1, nmos
                    l = l + microb(k,j)
                  end do
                end if
                microb(i,j) = 1 - microb(i,j)
              end do
              iphase(j) = 1 - 2*mod(l,2)
            end do
          end if
        end if
        ne = 0
        do i = 1, nmos
          ne = ne + microb(i,1) + microa(i,1)
        end do
!***********************************************************************
!
!   Now to work out the phase of the microstates.  The defined order
!   of M.O. occupancy is (alpha-1)(beta-1)(alpha-2)(beta-2) ...
!   If this order is permuted an odd number of times, the phase
!   will be negative.
!
        do i = 1, nstate
!
!    Load into IPERMA(1-2,I) and IPERMB(1-2,I) the locations of the
!    electrons.
!
          k = 0
          do j = 1, nmos
            if (microa(j,i) /= 1) cycle
            k = k + 1
            iperma(k,i) = j
          end do
          k = 0
          do j = 1, nmos
            if (microb(j,i) /= 1) cycle
            k = k + 1
            ipermb(k,i) = j
          end do
        end do
        if (posita .neqv. positb) then
!
!  Calculate determinant of M.O. transform in order to
!  define character of half-filled shell
!
          k = 0
          do i = 1, nmos
            work(k+1:nmos+k) = t2(i,:nmos)
            k = nmos + k
          end do
          call minv (work, nmos, det)
        else
          det = 1.D0
        end if
!
!   The big loop to fill T4
!
        do i = 1, nstate
          nai = nalmat(i)
          do j = 1, nstate
            if (nalmat(j) /= nai) cycle
             nbi = ne - nai
!
!    NAI = Number of alpha electrons
!    NBI = Number of beta electrons
!
! Fill the matrix for alpha electrons
!
            do ii = 1, nai
              do jj = 1,nai
                k = (ii - 1)*nai + jj
                matrix(k) = t2(iperma(jj,j),iperma(ii,i))
              end do
            end do
            suma = determinant(matrix, nai )
!
! Fill the matrix for beta electrons
!
            do ii = 1, nbi
              do jj = 1,nbi
                k = (ii - 1)*nbi + jj
                matrix(k) = t2(ipermb(jj,j),ipermb(ii,i))
              end do
            end do
            sumb = determinant(matrix, nbi )
            t4(i,j) = suma*sumb*iphase(i)*iphase(j)*det
          end do
        end do
        if (debug) then
          write (iw, *) ' State Transform for State', istate, &
            ' under Operation', ioper
          call matout (t4, t4, nstate, nstate, lab)
        end if
      end if
!
!    Now to perform <State(ISTATE) | Transform(IOPER) | State(ISTATE)>
!
      do j = 1, nstate
        sum = 0.D0
        do k = 1, nstate
          sum = sum + conf(k+(istate-1)*nstate)*t4(j,k)
        end do
        work(j) = sum
      end do
      sum = ddot(nstate,work(:nstate),1,conf(1+(istate-1)*nstate:istate*nstate),1)
      charst = sum
      return
      end function charst
