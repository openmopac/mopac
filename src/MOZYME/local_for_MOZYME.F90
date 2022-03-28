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

subroutine local_for_MOZYME (type)
    use molkst_C, only: norbs, numat, nelecs
    use MOZYME_C, only : cocc_dim, icocc_dim, &
       & cvir_dim, icvir_dim, cocc, icocc, ncf, ncocc, iorbs, nncf, &
       & cvir, icvir, nce, ncvir, nnce
    use chanel_C, only: iw
    implicit none
    character (len=*), intent (in) :: type
    integer :: i, nocc, nvir, alloc_stat
    double precision :: totij, sum
    double precision, dimension(:), allocatable :: psi1, psi2, axiiii
    integer, dimension(:), allocatable :: nf, nl
    integer, dimension(:,:), allocatable :: ioc
    allocate (psi1(norbs), psi2(norbs), axiiii(norbs), nf(numat), &
         & nl(numat), ioc(2,numat), stat=alloc_stat)
    if (alloc_stat /= 0) then
      call memory_error ("localz")
      goto 1100
    end if
    if (type == "OCCUPIED") then
      nocc = nelecs / 2
      !
      do i = 1, 100
        call localize_for_MOZYME (cocc, cocc_dim, icocc, icocc_dim, ncf, ncocc, &
             & nocc, iorbs, psi1, psi2, axiiii, nf, nl, ioc, nncf, totij, sum)
   !     write (iw, "(10x,'NUMBER OF ITERATIONS =',i4,/,10x,'LOCALIZATION VALUE =',f14.9)") i, sum
        if (totij < 1.d-5) exit
      end do
      write (iw, "(10x,'NUMBER OF ITERATIONS =',i4,/,10x,'LOCALIZATION VALUE =',f14.9,/)") i, sum
!
      call MOZYME_eigs(nocc)
!
    else if (type == "VIRTUAL") then
      nvir = norbs - nelecs / 2
      !
      do i = 1, 100
        call localize_for_MOZYME (cvir, cvir_dim, icvir, icvir_dim, nce, ncvir, &
             & nvir, iorbs, psi1, psi2, axiiii, nf, nl, ioc, &
             & nnce, totij, sum)
        if (totij < 1.d-5) exit
      end do
    else
      write (iw,*) " Error"
      call mopend ("Error in LOCAL")
    end if
    deallocate (psi1, psi2, axiiii, nf, nl, ioc)
1100 continue
end subroutine local_for_MOZYME
subroutine localize_for_MOZYME (c, n289, ic, n267, nc, ncstrt, nmos_loc, iorbs, psi1, &
     & psi2, axiiii, nf, nl, ioc, nnc_loc, totij, total)
    use molkst_C, only: numat, norbs, natoms
    implicit none
    integer, intent (in) :: n267, n289, nmos_loc
    double precision, intent (out) :: totij, total
    integer, dimension (n267), intent (in) :: ic
    integer, dimension (nmos_loc), intent (in) :: nc, ncstrt, nnc_loc
    integer, dimension (numat), intent (inout) :: nf, nl
    integer, dimension (numat), intent (in) :: iorbs
    integer, dimension (2, numat), intent (inout) :: ioc
    double precision, dimension (n289), intent (inout) :: c
    double precision, dimension (nmos_loc), intent (inout) :: axiiii
    double precision, dimension (norbs), intent (inout) :: psi1, psi2
    integer :: i, i1, i5, i8, ii, ij, ijorb, il, il1, iloop, j, j1, j5, j8, &
   & jl, jl1, jloop, k, k1, l, m
    double precision :: aij, bij, ca, dii, dij, djj, sa, xiiii, xiijj, &
   & xijij, xijjj, xjiii, xjjjj
   !
   !**********************************************************************
   !
   !   LOCALISATION SUBROUTINE
   ! ON INPUT
   !        C      = EIGENVECTORS IN AN MDIM*MDIM MATRIX
   !        NCSTRT = STARTING ADDRESS FOR ATOMIC ORBITALS IN LMO'S
   !        NMOS = NUMBER OF LMO'S TO BE USED
   !        NNC    = STARTING ADDRESS OF ATOMS IN LMO'S
   !        NLAST   = INTEGER ARRAY OF ATOM ORBITAL COUNTERS
   !        NFIRST   = INTEGER ARRAY OF ATOM ORBITAL COUNTERS
   !
   !       SUBROUTINE MAXIMIZES (PSI)**4
   !
   !**********************************************************************
   !  DUMMY STATEMENT TO USE COMMON BLOCK
    i = natoms
   !
   !   Build all the DII
   !
    l = 0
    xiiii = 0.d0
    do k = 1, nmos_loc
      l = nnc_loc(k)
      m = ncstrt(k)
      axiiii(k) = 0.d0
      do i = 1, nc(k)
        l = l + 1
        i1 = ic(l)
        dii = 0.d0
        do j = 1, iorbs(i1)
          m = m + 1
          dii = dii + c(m) ** 2
        end do
        axiiii(k) = axiiii(k) + dii * dii
      end do
    end do
    total = 0.d0
    totij = 0.d0
    do iloop = 1, nmos_loc
      i5 = ncstrt(iloop)
      i8 = nnc_loc(iloop)
      do jloop = 1, nmos_loc
        j8 = nnc_loc(jloop)
        j5 = ncstrt(jloop)
        if (jloop /= iloop) then
          do i = 1, 2
            do j = 1, 2
              if (ic(i+i8) == ic(j+j8)) go to 1000
            end do
          end do
          cycle
1000      ij = 0
          ijorb = 0
          il = 0
          do i = 1, nc(iloop)
            i1 = ic(i+i8)
            jl = 0
            do j = 1, nc(jloop)
              j1 = ic(j+j8)
              if (i1 == j1) then
                ij = ij + 1
                nf(ij) = ijorb + 1
                nl(ij) = ijorb + iorbs(i1)
                ioc(1, ij) = il
                ioc(2, ij) = jl
                jl1 = jl
                il1 = il
                do k = 1, iorbs(i1)
                  ijorb = ijorb + 1
                  il1 = il1 + 1
                  jl1 = jl1 + 1
                  psi1(ijorb) = c(il1+i5)
                  psi2(ijorb) = c(jl1+j5)
                end do
              end if
              jl = jl + iorbs(j1)
            end do
            il = il + iorbs(i1)
          end do
          xijjj = 0.d0
          xjiii = 0.d0
          xijij = 0.d0
          xiijj = 0.d0
          do k1 = 1, ij
            dij = 0.d0
            dii = 0.d0
            djj = 0.d0
            do k = nf(k1), nl(k1)
              dij = dij + psi1(k) * psi2(k)
              dii = dii + psi1(k) * psi1(k)
              djj = djj + psi2(k) * psi2(k)
            end do
            xijjj = xijjj + dij * djj
            xjiii = xjiii + dij * dii
            xijij = xijij + dij * dij
            xiijj = xiijj + dii * djj
          end do
          if (xiijj >= 0.001d0) then
            xiiii = axiiii(iloop)
            xjjjj = axiiii(jloop)
            aij = xijij - (xiiii+xjjjj-2.0d0*xiijj) / 4.0d0
            bij = xjiii - xijjj
            ca = Sqrt (aij*aij+bij*bij)
            sa = aij + ca
            if (sa > 1.d-14) then
              ca = (1.0d0+Sqrt((1.0d0-aij/ca)/2.0d0)) / 2.0d0
              sa = Sqrt (1.0d0-ca)
              ca = Sqrt (ca)
              totij = totij + sa
              ii = 0
              do k = 1, ij
                il = 0
                do i = nf(k), nl(k)
                  il = il + 1
                  ii = ii + 1
                  c(ioc(1, k)+il+i5) = ca * psi1(ii) + sa * psi2(ii)
                  c(ioc(2, k)+il+j5) = -sa * psi1(ii) + ca * psi2(ii)
                end do
              end do
            end if
          end if
        end if
      end do
      total = total + xiiii
    end do
    return
end subroutine localize_for_MOZYME
