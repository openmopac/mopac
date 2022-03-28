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

subroutine modchg ()
   !***********************************************************************
   !
   !   MODGRA prints the charge due to (a) backbone residue atoms, and
   !   (b) all side-chain atoms in each residue.
   !
   !***********************************************************************
    use molkst_C, only: numat, id, maxtxt
    use MOZYME_C, only : nres, at_res, res_start
    use common_arrays_C, only : q, txtatm
    use chanel_C, only: iw
    implicit none
    double precision, allocatable :: work(:)
    logical, allocatable :: l_used(:)
    logical :: first = .true.
    character :: nam_het(600)*3, het_num(600)*4
    double precision :: total_charge = 0.d0, res_charge = 0.d0
    integer :: i, j, n_cat, n_ani, cations(50), anions(50), n_het = 0
   !
   ! ... Executable Statements ...
   !
    allocate (work(0:numat + id))
    call build_res_start_etc()
    work = 0.d0
    if (.not. allocated(q)) then
      write(iw,"(a)") "Density matrix is not available.  Run a 1SCF."
      return
    end if
    allocate(l_used(numat + id))
   !
   !    WORK collects charge contributions from atoms.  The negative
   !    addresses refer to the residue backbone, the positive addresses
   !    refer to the side-chain.
   !
   !    (Backbone atoms are -NH-CH-CO-)
   !
    work = 0.d0
    do i = 1, numat
      j = at_res(i)
      work(j) = work(j) + q(i)
      total_charge = total_charge + q(i)
      l_used(i) = (j > 0)
    end do
    write (iw,"(/9x,a,/)") " NET CHARGE ON RESIDUES"
    write (iw,"(a)")"      Residue         Charge  Anion or"
    write (iw,"(a)")"                              Cation?"
    n_cat = 0
    n_ani = 0
    do i = 1, nres
      if (res_start(i) < 1) cycle
      if (work(i) > 0.5d0) then
         write (iw, "(6x,a,SP,F13.3, a)") txtatm(res_start(i))(18:maxtxt), work(i), "  CATION"
         n_cat = min(50, n_cat + 1)
         cations(n_cat) = res_start(i)
      else if (work(i) < -0.5d0) then
         write (iw, "(6x,a,F13.3, a)") txtatm(res_start(i))(18:maxtxt), work(i), "  ANION"
         n_ani = min(50, n_ani + 1)
         anions(n_ani) = res_start(i)
      else
        if (work(i) < 0.d0 .and. work(i) > -0.0005d0) work(i) = 0.d0
        if (txtatm(res_start(i))(18:20) /= "HOH" .or. abs(work(i)) > 0.1d0) &
         write (iw, "(6x,a,SP,F13.3)") txtatm(res_start(i))(18:maxtxt), work(i)
      end if
      res_charge = res_charge + work(i)
    end do
    work = 0.d0
    do i = 1, numat
      if (.not. l_used(i)) then
        do j = 1, n_het
          if (txtatm(i)(18:20) == nam_het(j) .and. txtatm(i)(23:26) == het_num(j)) exit
        end do
        if (j > n_het) then
          n_het = n_het + 1
          if (n_het > 600) exit
          nam_het(j) = txtatm(i)(18:20)
          het_num(j) = txtatm(i)(23:26)
        end if
        work(j) = work(j) + q(i)
      end if
    end do
    if (n_het > 0) then
      do i = 1, n_het
        if (abs(work(i)) < 0.1d0) cycle
        if (first) then
          first = .false.
          write (iw,"(/6x,a,/2x,a,/)") "NET CHARGES ON HETERO-GROUPS  (|Charges| less than 0.1 not printed)"
          write (iw,"(a)")"    HETERO-GROUP      Charge "
          write (iw,"(a)")" "
        end if
        write (iw, "(6x,a,F14.3)") txtatm(res_start(i))(18:maxtxt), work(i)
      end do
    end if
    write(iw,*)
!
!  Print out salt bridges
!
    if (n_cat > 0 .and. n_ani > 0) then
      call set_up_dentate()
      call find_salt_bridges(cations, anions, n_cat, n_ani)
    end if
    write(iw,'(/)')
    return
  end subroutine modchg
  subroutine build_res_start_etc()
    use MOZYME_C, only : nres, at_res, res_start
    use molkst_C, only: numat, id
    use common_arrays_C, only : txtatm, nat, nbonds, ibonds
    implicit none
    integer :: i, j
    character :: ren_name(numat)*9, residue*9
!
! Work out which residue each atom belongs to, and the location of the first atom
! in each residue.
!
    nres = 0
    ren_name = " "
    if (.not. allocated(at_res)) allocate(at_res(numat + id))
!
!  Do all non-hydrogen atoms
!
     do i = 1, numat
      if (nat(i) == 1 .and. nbonds(i) == 1) cycle
      residue = txtatm(i)(18:26)
      do j = 1, nres
        if (ren_name(j) == residue) exit
      end do
      if (j > nres) then
        nres = nres + 1
        ren_name(j) = residue
        res_start(j) = i
      end if
      at_res(i) = j
     end do
!
! Do all hydrogen atoms
!
    do i = 1, numat
      if (nat(i) /= 1 .or. nbonds(i) /= 1) cycle
        j = ibonds(1,i)
        residue = txtatm(j)(18:26)
      do j = 1, nres
        if (ren_name(j) == residue) exit
      end do
      at_res(i) = j
    end do
    return
  end subroutine build_res_start_etc
