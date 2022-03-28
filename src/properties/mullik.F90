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

      subroutine mullik()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!
      use molkst_C, only : numat, nelecs, nclose, nopen, fract,  &
        keywrd, norbs, id, verson, method_pm6, uhf, nalpha, nbeta, &
        numcal, escf, line
!
      use symmetry_C, only : jndex, namo
!
      use maps_C, only : rxn_coord
!
      use common_arrays_C, only : nfirst, nlast, nat, coord, &
      & c, h, pb, tvec, eigs, q, eigb, p, ifact, cb
!
      use parameters_C, only : zs, zp, zd, betas, betap, tore
!
      use chanel_C, only : igpt, gpt_fn, iw
!
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, if, il, im1, k, ii, j, jf, jl, ij, icalcn=0, mo_l, mo_u
      character, dimension(:), allocatable :: namo_tmp*4
      character :: paras*1
      double precision, dimension(norbs) :: eig
      double precision, dimension(:), allocatable :: store, vecs, store_h
      double precision, dimension(:, :), allocatable :: work
      double precision :: bi, bj, sum, summ, q2(numat)
      logical :: graph, graph_formatted, namo_ok
! GBR_new_addition
      integer :: nlower
!
      save :: graph, graph_formatted, icalcn, mo_l, mo_u
!

!-----------------------------------------------
!*********************************************************************
!
!  MULLIK DOES A MULLIKEN POPULATION ANALYSIS
! ON INPUT     C      =  SQUARE ARRAY OF EIGENVECTORS.
!              H      =  PACKED ARRAY OF ONE-ELECTRON MATRIX
!              VECS   =  WORKSTORE OF SIZE AT LEAST NORBS*NORBS
!              STORE  =  WORKSTORE OF SIZE AT LEAST (NORBS*(NORBS+1))/2
!
!*********************************************************************


     allocate(store((norbs*(norbs + 1))/2), store_h((norbs*(norbs + 1))/2), &
       vecs(norbs**2), stat = i)
     if (i /= 0) then
        write(iw,*)" Unable to allocate memory for eigenvector matrices in MULLIK"
        call mopend("Unable to allocate memory for eigenvector matrices in MULLIK")
        return
      end if
!*********************************************************************
!
!  FIRST, RE-CALCULATE THE OVERLAP MATRIX
!
!*********************************************************************
      if (icalcn /= numcal) then
        icalcn = numcal
        graph = index(keywrd,'GRAPH') /= 0
        graph_formatted = index(keywrd,'GRAPHF') /= 0
        if (allocated(ifact)) deallocate(ifact)
        allocate (ifact(norbs + 1))
        do i = 1, norbs
          ifact(i) = (i*(i - 1))/2
        end do
        ifact(norbs+1) = (norbs*(norbs + 1))/2
        if (graph) then
          close (unit=igpt, iostat=i)
          if (graph_formatted) then
            open(unit=igpt, file=gpt_fn(:len_trim(gpt_fn) - 3)//"mgf", form='FORMATTED', &
              status="UNKNOWN", iostat = i)
          else
            open(unit=igpt, file=gpt_fn, form='UNFORMATTED', status="UNKNOWN", iostat = i)
          end if
          if (i /= 0) then
            write(iw,*)" File '"//gpt_fn(:len_trim(gpt_fn))//"' is unavailable for use"
            call mopend("File '"//gpt_fn(:len_trim(gpt_fn))//"' is unavailable for use")
            return
          end if
          rewind igpt
        end if
        mo_l = 1
        mo_u = norbs
      end if
      call chrge (p, q2)
      do i = 1, numat
        q(i) = tore(nat(i) ) - q2(i)
      end do
      store_h = h
      do i = 1, numat
        if = nfirst(i)
        il = nlast(i)
        im1 = i - 1
        bi = betas(nat(i))
        do k = if, il
          ii = (k*(k - 1))/2
          do j = 1, im1
            jf = nfirst(j)
            jl = nlast(j)
            bj = betas(nat(j))
!  THE  +1.D-14 IS TO PREVENT POSSIBLE ERRORS IN THE DIAGONALIZATION.
            ij = ii + jf
            h(ij) = 2.D0*h(ij)/(bi + bj) + 1.D-14
            store(ij) = h(ij)
            bj = betap(nat(j))
            bj = betap(nat(j))
            h(ii+jf+1:jl+ii) = 2.D0*h(ii+jf+1:jl+ii)/(bi + bj) + 1.D-14
!  THE  +1.D-14 IS TO PREVENT POSSIBLE ERRORS IN THE DIAGONALIZATION.
            store(ii+jf+1:jl+ii) = h(ii+jf+1:jl+ii)
          end do
          store(ii+if:k+ii) = 0.D0
          h(ii+if:k+ii) = 0.D0
          bi = betap(nat(i))
        end do
      end do
      store(ifact(2:norbs+1)) = 1.D0
      h(ifact(2:norbs+1)) = 1.D0
      call rsp (h, norbs, eig, vecs)
      do i = 1, norbs
        eig(i) = 1.D0/sqrt(abs(eig(i)))
      end do
      if (.not. allocated(work)) allocate(work(norbs,norbs))
      ij = 0
      do i = 1, norbs
        do j = 1, i
          ij = ij + 1
          sum = 0.D0
          do k = 1, norbs
            sum = sum + vecs(i+(k-1)*norbs)*eig(k)*vecs(j+(k-1)*norbs)
          end do
          work(i,j) = sum
          work(j,i) = sum
        end do
      end do
      if (graph) then
        if (graph_formatted) then
          write(igpt,"(i5, a)") numat," MOPAC-Graphical data Version "//verson
          do i = 1, numat
            write(igpt,"(i4,f13.7,2f12.7,f9.4)")nat(i),(coord(j,i),j=1,3), q(i)
          end do
          do i = 1, numat
            write(igpt,"(3f11.7)")zs(nat(i)), zp(nat(i)), zd(nat(i))
          end do
          if (index(keywrd," ALLV") == 0) then
            mo_l = max(1, max(nalpha,nclose) - 8)
            mo_u = min(norbs, max(nalpha,nclose) + 7)
          end if
          paras = char(ichar("1") +int(log10(norbs + 0.05)))
          namo_ok =allocated(namo)
          do i = mo_l, mo_u
            j = 0
            if (uhf .and. i <= nalpha) j = 1
            if ( .not. uhf) then
              if (i <= nclose) then
                j = 2
              else if ( i <= nopen .and. abs(fract - 1.d0) < 1.d-4) then
                j = 1
              end if
            end if
            if (namo_ok) then
              write(igpt,"(a,i2,1x,i"//paras//",a,f10.4)")' ORBITAL', j, jndex(i),namo(i), eigs(i)
            else
              write(igpt,"(a,i2,1x,'NONE',f10.4)")' ORBITAL', j, eigs(i)
            end if
            write(igpt,"(5d15.8)")(c(j,i), j=1, norbs)  !   Write out the spin-free M.O.s, or alpha M.O.s, if UHF
          end do
          write(igpt,"(a,i"//paras//",a,i"//paras//",a)")' INVERSE_MATRIX[',norbs,'x',norbs,']='
          do i = 1, norbs
            write(igpt,"(5d15.8)")(work(j,i), j=1, i) !   This is the inverse-square-root of the overlap matrix, lower half triangle.
          end do
          if (rxn_coord < 1.d8) then
            write(igpt,'(a,f12.4)')" Reaction coordinate:", rxn_coord
            write(igpt,'(a,f12.4)')" Heat_of_formation:  ", escf
          end if
          line = trim(keywrd)
          if (method_pm6 .and. index(keywrd," PM6") == 0) line = " PM6"//trim(line)
          if (uhf .and. index(keywrd," UHF") == 0) line = " UHF"//trim(line)
          write(igpt,"(a)")" Keywords:"//trim(line)
          if (uhf) then
            if (index(keywrd," ALLV") == 0) then
              mo_l = max(1, max(nbeta,nclose) - 8)
              mo_u = min(norbs, max(nbeta,nclose) + 7)
            end if
            if (namo_ok) then
              allocate(namo_tmp(norbs))
              namo_tmp(:norbs) = namo(:norbs)
              call symtrz (cb, eigb, 1, .TRUE.)
            end if
            do i = mo_l, mo_u
              j = 0
              if (uhf .and. i <= nbeta) j = 1
              if (namo_ok) then
                write(igpt,"(a,i2,1x,i"//paras//",a,f10.4)")' ORBITAL', j, jndex(i),namo(i), eigb(i)
              else
                write(igpt,"(a,i2,1x,'NONE',f10.4)")' ORBITAL', j, eigb(i)
              end if
              write(igpt,"(5d15.8)")(cb(j,i), j=1, norbs)  !   Write out the beta M.O.s, if UHF
            end do
            if (namo_ok) then
              namo(:norbs) = namo_tmp(:norbs)
              deallocate(namo_tmp)
            end if
          end if
        else
!
! WRITE TO DISK THE FOLLOWING DATA FOR GRAPHICS CALCULATION, IN ORDER:
!
!      NUMBER OF ATOMS, ORBITAL, ELECTRONS
!      ALL ATOMIC COORDINATES
!      ORBITAL COUNTERS
!      ORBITAL EXPONENTS, S, P, AND D, AND ATOMIC NUMBERS
!      EIGENVECTORS (M.O.S NOT RE-NORMALIZED)
!      INVERSE-SQUARE ROOT OF THE OVERLAP MATRIX.
!

          write (igpt) numat, norbs, nelecs, ((coord(i,j),j=1,numat),i=1,3)
          write (igpt) (nlast(i),nfirst(i),i=1,numat)
          write (igpt) (zs(nat(i)),i=1,numat), (zp(nat(i)),i=1,numat), &
          (zd(nat(i)),i=1,numat), (nat(i),i=1,numat)
          write (igpt) c
          write (igpt) work
          write (igpt) id, tvec
        end if
        h = store_h
        if (index(keywrd,'MULLIK') == 0) return
      end if
      call mult (c, work, vecs, norbs)
      i = -1
      nlower = (norbs*(norbs + 1))/2
      call density_for_GPU (vecs, fract, nclose, nopen, 2.d0, nlower, norbs, 2, pb, 3)
!
      pb = pb*store
      summ = 0.D0
      do i = 1, norbs
        sum = 0
        do j = 1, i
          sum = sum + pb(ifact(i)+j)
        end do
        do j = i + 1, norbs
          sum = sum + pb(ifact(j)+i)
        end do
        summ = summ + sum
        pb(ifact(i+1)) = sum
      end do
      h = store_h
      return
      end subroutine mullik
