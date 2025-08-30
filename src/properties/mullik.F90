! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

      subroutine mullik(popmat)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!
      use molkst_C, only : numat, nelecs, nclose, nopen, fract,  &
        keywrd, norbs, id, verson, method_pm6, uhf, nalpha, nbeta, &
        numcal, escf, line, mpack, nalpha_open, nbeta_open
!
      use symmetry_C, only : jndex, namo
!
      use maps_C, only : rxn_coord
!
      use common_arrays_C, only : nfirst, nlast, nat, coord, &
      & c, h, pb, tvec, eigs, q, eigb, p, ifact, cb
!
      use parameters_C, only : zs, zp, zd, betas, betap, betad, tore, natorb
!
      use chanel_C, only : igpt, gpt_fn, iw
!
      implicit none
      double precision, intent(out) :: popmat(mpack)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, if, il, im1, k, ii, j, jf, jl, ij, icalcn=0, mo_l, mo_u, norbi, norbj
      character, dimension(:), allocatable :: namo_tmp*4
      character :: paras*1
      double precision, dimension(norbs) :: eig
      double precision, dimension(:), allocatable :: store, vecs, overlap
      double precision, dimension(:, :), allocatable :: work
      double precision :: bi(9), bj(9), sum, summ, q2(numat)
      logical :: graph, graph_formatted, namo_ok
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
! ON OUTPUT    POPMAT =  PACKED ARRAY OF MULLIKEN POPULATION MATRIX
!
!*********************************************************************


     allocate(p_mullik(mpack), overlap(mpack), vecs(norbs**2), stat = i)
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
              iostat = i)
          else
            open(unit=igpt, file=gpt_fn, form='UNFORMATTED', iostat = i)
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
      ! overlap matrix construction from to_screen, including d orbitals
      overlap = 0.d0
      do i = 1, numat
        if = nfirst(i)
        im1 = i - 1
        ni = nat(i)
        bi = betas(ni)
        bi(1) = betas(ni)*0.5D0
        bi(2) = betap(ni)*0.5D0
        bi(3) = bi(2)
        bi(4) = bi(2)
        bi(5) = betad(ni)*0.5D0
        bi(6) = bi(5)
        bi(7) = bi(5)
        bi(8) = bi(5)
        bi(9) = bi(5)
        norbi = natorb(ni)
        do j = 1, im1
          nj = nat(j)
          bj(1) = betas(nj)*0.5D0
          bj(2) = betap(nj)*0.5D0
          bj(3) = bj(2)
          bj(4) = bj(2)
          bj(5) = betad(nj)*0.5D0
          bj(6) = bj(5)
          bj(7) = bj(5)
          bj(8) = bj(5)
          bj(9) = bj(5)
          norbj = natorb(nj)
          jf = nfirst(j)
          do ii = 1, norbi
            do jj = 1, norbj
              ij = ((if + ii - 1)*(if + ii - 2))/2 + jf + jj - 1
              overlap(ij) = h(ij)/(bi(ii) + bj(jj))
            end do
          end do
        end do
      end do
      overlap(ifact(2:norbs+1)) = 1.D0
      call rsp (overlap, norbs, eig, vecs)
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
        if (index(keywrd,'MULLIK') == 0) return
      end if
      if (uhf) then
        call mult (c, work, vecs, norbs)
        call density_for_GPU (vecs, fract, nalpha, nalpha_open, 1.d0, mpack, norbs, 2, p_mullik, 3)
        popmat(:mpack) = p_mullik(:mpack)*overlap(:mpack)
        call mult (cb, work, vecs, norbs)
        call density_for_GPU (vecs, fract, nbeta, nbeta_open, 1.d0, mpack, norbs, 2, p_mullik, 3)
        popmat(:mpack) = popmat(:mpack) + p_mullik(:mpack)*overlap(:mpack)
      else
        call mult (c, work, vecs, norbs)
        call density_for_GPU (vecs, fract, nclose, nopen, 2.d0, mpack, norbs, 2, p_mullik, 3)
        popmat(:mpack) = p_mullik(:mpack)*overlap(:mpack)
      end if
      summ = 0.D0
      do i = 1, norbs
        sum = 0
        do j = 1, i
          sum = sum + popmat(ifact(i)+j)
        end do
        do j = i + 1, norbs
          sum = sum + popmat(ifact(j)+i)
        end do
        summ = summ + sum
        popmat(ifact(i+1)) = sum
      end do
      return
      end subroutine mullik
