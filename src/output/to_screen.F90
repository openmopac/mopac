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

  subroutine to_screen(text0)
  use chanel_C, only : iw0
  use molkst_C, only : keywrd, keywrd_txt
  implicit none
  character (len=*) :: text0
  character (len=200) :: text
  integer :: i
  i = len_trim(text0)
  if (i == 0) then
!
! Length not passed, so work it out
!
    do i = 1, 200
      if (ichar(text0(i:i)) == 0) exit
    end do
    i = i - 1
  end if
  text = text0(:i)
  if (text(:min(len_trim(text), 8)) /= "To_file:") then
    if (iw0 > -1) then
      write(iw0,"(a)") trim(text)
      call flush (iw0)
    end if
  else
    if (index(keywrd, " AUX") == 0 .and. index(keywrd_txt, " AUX") == 0) return
    call current_version(text)
  end if
  end subroutine to_screen

  subroutine current_version (text)
  use chanel_C, only : output_fn
!
  use to_screen_C, only : rot, xyzmom, dip, dipt, travel, freq, &
  redmas, fcint, cnorml
!
  use common_arrays_C, only : nat, tvec, coord, nfirst, nlast, c, grad, &
  h, eigs, q, geo, na, nb, nc, p, pa, pb, eigb, cb, fmatrx, atmass, &
  txtatm, chrg, ipKa_sorted, pKa_sorted, f, fb, bondab, T_range, HOF_tot, H_tot, Cp_tot, S_tot
!
  use parameters_C, only : zs, zp, zd, npq, betas, betap, betad, tore, natorb
!
  use symmetry_C, only : name, jndex, namo, state_spin, state_Irred_Rep, state_QN
!
  use polar_C, only : omega, alpavg
!
  use molkst_C, only : numat, norbs, escf, nelecs, nclose, nopen, verson, &
  method_am1, method_mndo, method_pm3, method_rm1, method_mndod, method_pm6, &
  method_pm7, nvar, koment, keywrd, zpe, id, density, natoms, formula, press, voigt, &
  uhf, nalpha, nbeta,  gnorm, mozyme, mol_weight, ilim, &
  line, nscf, time0, sz, ss2, no_pKa, title, jobnam, job_no, fract
!
  use MOZYME_C, only : ncf, ncocc, noccupied, icocc_dim, cocc_dim, nvirtual, icvir_dim, &
  nncf, iorbs, cocc, icocc, ncvir, nnce, nce, icvir, cvir, tyres, size_mres, &
  cvir_dim, idiag
!
  use elemts_C, only : elemnt
!
  use chanel_c, only : iw
!
  use funcon_C, only : a0
!
  use cosmo_C, only : area, solv_energy, cosvol
!
  use meci_C, only : deltap, nmos, occa, microa, microb, lab, nstate, vectci, &
    root_requested, eig, nelec, rjkab
!
  use maps_C, only : rxn_coord, rc_escf, ekin, lparam, latom
!
  use drc_C, only: time
!
#if MOPAC_F2003
  USE, INTRINSIC :: IEEE_ARITHMETIC
#endif
!
  implicit none
  character (len=*) :: text

  integer :: i, j, k, l, hook = 50, if, im1, jj, ii, jf, ij, opt_hook, ni, nj, &
  norbi, norbj, kl, ku, i1, j1, ic, mos, moa_lower, moa_upper, jloop = -1, &
  mob_lower, mob_upper, i_map
  logical :: normal, force, geom_opt, int_force_const, error, finished, &
  first = .true., LMO, LMO_alpha = .true., end_of_job, opend, compressed, &
  eigen, mullik, esp, irc_drc, irc, lnmoas, lnmobs, loc_mos, new_calcn = .true., &
  polar, reaction_path, l_RC, L_overlap, L_MO_s, L_Density
  character :: atorbs(9)*2, idate*24, atoms*3, paras*3, orbs*3, orbs2*5, fmt*3, &
  fmt1*5
  double precision, dimension(3) :: convert = (/1.d0,57.29577951d0,57.29577951d0/)
  double precision :: work(norbs),  bi(9), bj(9), sum, limit, bk, eig_min
  double precision, allocatable :: overlap(:), c_lmo(:), comp(:), store_eigs(:)
  double precision, allocatable :: overlap2(:,:)
  double precision, external :: seconds, reada
  character, dimension(:), allocatable :: namo_tmp*4, letters*1
  character :: fmt8p4*5, fmt9p3*5, fmt10p1*5, fmt10p2*5, fmt10p3*5, fmt9p4*5, fmt10p4*5, fmtnnp4*5, fmt9p5*5, &
  fmt13p5*5, fmt13p6*5, fmt7p4*5, irc_or_drc*1, num*1
  integer :: iwork(norbs), ifact(norbs + 1)
  integer, allocatable :: icomp(:), eigs_map(:)
  integer, external :: ijbo
  character, dimension (23) :: tyr
  save :: first, opt_hook, atoms, paras, orbs, compressed, fmt8p4, fmt9p3, fmt10p1, fmt10p2, fmt10p3, fmt9p4, fmt10p4, &
  fmt9p5, fmt13p5, fmt13p6, fmt7p4, mos, irc, irc_or_drc, moa_lower, moa_upper, jloop, L_overlap, L_MO_s, L_Density, &
  mob_lower, mob_upper, lnmoas, lnmobs, loc_mos, new_calcn, orbs2, l_RC
  data atorbs/ ' S', 'PX', 'PY', 'PZ', 'X2', 'XZ', 'Z2', 'YZ', 'XY'/
  data tyr / "G", "A", "V", "L", "I", "S", "T", "D", "N", "K", "E", "Q", &
     & "R", "H", "F", "C", "W", "Y", "M", "P", "P", "P", "?" /

!
!   Take special action to output essential data in a compact ASCII form
!
  inquire(unit=hook, opened=opend)
  if (.not. opend .or. new_calcn) then
    new_calcn = .false.
    first = .true.
    if (index(text, "Leaving MOPAC") + index(text, "END_OF_JOB") /= 0) then
      new_calcn = .true.
      return
    end if
    l_RC = (index(keywrd, " IRC") + index(keywrd," DRC") > 0)
    open(unit=hook,file=output_fn(:len_trim(output_fn) - 4)//".aux")
    if (.not. opend) then
      write(hook,"(a)")" START OF MOPAC PROGRAM"
    end if
    write(hook,"(a)")" START OF MOPAC FILE"
    write(hook,"(a)")" ####################################"
    write(hook,"(a)")" #                                  #"
    write(hook,"(a)")" #       Start of Input data        #"
    write(hook,"(a)")" #                                  #"
    write(hook,"(a)")" ####################################"
    write(hook,"(a,a)")" MOPAC_VERSION=",verson
    call fdate (idate)
    write(hook,"(3a)")" DATE=""",idate,""""
    if (method_mndo) then
      write(hook,"(a)")" METHOD=MNDO"
    else if (method_am1) then
      write(hook,"(a)")" METHOD=AM1"
    else if (method_pm3) then
      write(hook,"(a)")" METHOD=PM3"
    else if (method_pm6) then
      i = index(keywrd, " PM6") + 1
      j = index(keywrd(i:), " ") + i - 2
      write(hook,"(a)")" METHOD="//keywrd(i:j)
    else if (method_pm7) then
      write(hook,"(a)")" METHOD=PM7"
    else if (method_rm1) then
      write(hook,"(a)")" METHOD=RM1"
    else if (method_mndod) then
      write(hook,"(a)")" METHOD=MNDOD"
    end if
    line = koment
    if (len_trim(line) > 1) then
      do
        if (line(1:1) /= " ") exit
        line = line(2:)
      end do
    end if
    write(hook,"(a,a)")" TITLE=""",trim(line)//'"'
    line = trim(keywrd)
    do i = 1, len_trim(line)
      do
        if (line(i:i + 1) /= "  ") exit
        if (line(i:) == " ") exit
        line(i:) = line(i + 1:)
      end do
    end do
    write(hook,"(a,a)")" KEYWORDS=""",trim(line)//'"'
    line = title
    if (len_trim(line) > 0) then
      do
        if (line(1:1) /= " ") exit
        line = line(2:)
      end do
    end if
    write(hook,"(a,a)")" COMMENTS=""",trim(line)//'"'
    j = int(log10(numat*1.0001)) + 2
    write(atoms,'(i1,a,i1)')j, ".", j
    j = int(log10(numat*3.0001)) + 2
    write(paras,'(i1,a,i1)')j, ".", j
    j = int(log10(norbs*1.0001)) + 2
    write(orbs,'(i1,a,i1)')j, ".", j
    j = int(log10(norbs**2*1.0001)) + 2
    if (j > 9) then
      write(orbs2,'(i2,a,i2)')j, ".", j
    else
      write(orbs2,'(i1,a,i1)')j, ".", j
    end if
    mos = 10
    i = index(keywrd, " AUX(")
    if (i /= 0) then
      j = index(keywrd(i:),") ") + i
      k = i + 5
      if (ichar(keywrd(k:k)) >= ichar("0") .and. ichar(keywrd(k:k)) <= ichar("9")) then
        opt_hook = nint(reada(keywrd(i:j),1))
        opt_hook = max(0, min(opt_hook,100))
      else
        opt_hook = hook
      end if
      compressed = (index(keywrd(i:j), "COMP") /= 0)
      L_overlap = (index(keywrd(i:j), "XS") == 0)
      L_MO_s    = (index(keywrd(i:j), "XW") == 0)
      L_Density = (index(keywrd(i:j), "XP") == 0)
      k = index(keywrd(i:j), "MOS")
      if (k > 0) mos = nint(reada(keywrd(i:), k))
    else
      compressed = .false.
      opt_hook = hook
      L_overlap = .true.
      L_MO_s    = .true.
      L_Density = .true.
    end if
    if(index(keywrd," LARGE") == 0) then
        moa_lower = max(1, nclose - mos + 1, nalpha - mos + 1)
        moa_upper = min(norbs, max(nclose + mos, nalpha + mos))
        lnmoas = (moa_upper - moa_lower > -1)
        mob_lower = max(1, nclose - mos + 1, nalpha - mos + 1)
        mob_upper = min(norbs, max(nclose + mos, nalpha + mos))
        lnmobs = (mob_upper - mob_lower > -1)
      else
        moa_lower = 1
        moa_upper = norbs
        mob_lower = 1
        mob_upper = norbs
        lnmobs = .true.
        lnmoas = .true.
      end if
    k = 0
    if (i > 0) k = index(keywrd(i:j), "PRECISION")
    if (k > 0) k = nint(reada(keywrd(i:), k))
    write(fmt7p4,"(i2.2,'.',i2.2)")7 + k, 4 + k
    write(fmt8p4,"(i2.2,'.',i2.2)")8 + k, 4 + k
    write(fmt9p3,"(i2.2,'.',i2.2)")9 + k, 3 + k
    write(fmt9p4,"(i2.2,'.',i2.2)")9 + k, 4 + k
    write(fmt9p5,"(i2.2,'.',i2.2)")9 + k, 5 + k
    write(fmt10p1,"(i2.2,'.',i2.2)")10 + k, 1 + k
    write(fmt10p2,"(i2.2,'.',i2.2)")10 + k, 2 + k
    write(fmt10p3,"(i2.2,'.',i2.2)")10 + k, 3 + k
    write(fmt10p4,"(i2.2,'.',i2.2)")10 + k, 4 + k
    write(fmt13p5,"(i2.2,'.',i2.2)")13 + k, 5 + k
    write(fmt13p6,"(i2.2,'.',i2.2)")13 + k, 6 + k
!
!  All the data used in defining the starting system
!
    write(hook,"(a,i"//atoms//",a)")" ATOM_EL[",numat,"]="
    write(hook,"(40(' ',a2))")(elemnt(nat(i)), i=1,numat)
    write(hook,"(a,i"//atoms//",a)")" ATOM_CORE[",numat,"]="
    write(hook,"(40(' ',i2))")(nint(tore(nat(i))), i=1,numat)
    write(hook,"(a,i"//paras//",a)")" ATOM_X:ANGSTROMS[",3*numat, "]="
    write(hook,"(3f"//fmt10p4//")") ((coord(j,i),j=1,3), i=1,numat)
    if (norbs > 0) then
      write(hook,"(a,i"//orbs//",a)")" AO_ATOMINDEX[",norbs,"]="
      j = int(log10(norbs*1.0001)) + 2
      write(fmt1,'(i2,a,i2)')120/j, "i", j
      write(hook,"("//fmt1//")")((i,j = nfirst(i),nlast(i)), i = 1, numat)
      write(hook,"(a,i"//orbs//",a)")" ATOM_SYMTYPE[", norbs,"]="
      write(hook,"(40(' ',a2))")&
       ((atorbs(j - nfirst(i) + 1),j = nfirst(i),nlast(i)), i = 1, numat)
      write(hook,"(a,i"//orbs//",a)")" AO_ZETA[",norbs,"]="
!
!   Set up an array to hold all the atomic orbital exponents and principal quantum numbers
!
      j = 1
      do i = 1, numat
        if (nlast(i) - nfirst(i) > -1) then
          work(j) = zs(nat(i))
          iwork(j) = npq(nat(i),1)
          j = j + 1
        end if
        if (nlast(i) - nfirst(i) > 2) then
          work(j:j + 2) = zp(nat(i))
          iwork(j:j + 2) = npq(nat(i),2)
          j = j + 3
        end if
        if (nlast(i) - nfirst(i) > 7) then
          work(j:j + 4) = zd(nat(i))
          iwork(j:j + 4) = npq(nat(i),3)
          j = j + 5
        end if
      end do
      write(hook,"(10f"//fmt8p4//")") (work(i), i=1,norbs)
      write(hook,"(a,i"//orbs//",a)")" ATOM_PQN[",norbs,"]="
      write(hook,"(40i2)")(iwork(i), i=1,norbs)
    end if
    if (norbs > 0) write(hook,"(a,i"//orbs//")")" NUM_ELECTRONS=",nelecs
    if (id > 0) then
      write(hook,"(a,i1,a)")" INITIAL_TRANS_VECTS:ANGSTROMS[",id*3,"]="
      write(hook,"(3f"//fmt9p4//")")((tvec(j,i),j=1,3),i=1,id)
      if (density > 1.d-1) write(hook,"(a,d"//fmt13p6//",a)")" DENSITY:G/CM^3=",density
    end if
    write(hook,"(a)")" EMPIRICAL_FORMULA="""//formula(31:len_trim(formula))//""""
!
!  End of definition of starting system
!
  end if
  if (L_MO_s) L_MO_s = (moa_upper - moa_lower > -1)
  if (nelecs == 0) then
    L_MO_s = .false.
    L_overlap = .false.
    L_Density = .false.
    eigen = .false.
    lnmobs = .false.
    lnmoas = .false.
  end if
!
!  Set up options for various types of calculation
!
  geom_opt         = (index(text,"Geometry optimizing") /= 0)
  normal           = (index(text,"Normal output") /= 0)
  force            = (index(text,"Force output") /= 0)
  LMO              = (index(text,"LMO") /= 0)
  int_force_const  = (index(text,"Internal Force Constants") /= 0)
  error            = (index(text,": ERROR:") /= 0)
  finished         = (index(text,"Leaving MOPAC") /= 0)
  end_of_job       = (index(text,"END_OF_JOB") /= 0)
  mullik           = (index(text,"Mullik") /= 0)
  esp              = (index(text,"Esp") /= 0)
  irc_drc          = (index(text,"IRC-DRC") /= 0)
  loc_mos          = (index(text,"Localized") /= 0)
  polar            = (index(text,"POLAR") /= 0)
  reaction_path    = (index(text,"Reaction path") /= 0)
  if (geom_opt) geom_opt = (lparam == 0 .and. latom == 0)
  if (.not. l_RC .and. (geom_opt .or. first)) then
    if (first) then
      write(opt_hook,"(a)",iostat=i)" ####################################"
      if (i /= 0) then
        write(iw,'(/5x,a,i3,a)') "WARNING: Channel", opt_hook, " cannot be opened"
        open(unit=opt_hook, file=trim(jobnam)//'_opt.aux')
        write(iw,'(3x,a)') "- intermediate results will be sent to '"//trim(jobnam)//"_opt.aux'"
        if (opt_hook /= 0) then
          write(iw,'(5x,a,/)') "(Try using chanel 0 if you want to avoid having intermediate results written out.)"
        else
          write(iw,'(5x,a)') " "
        end if
        write(opt_hook,"(a)")" ####################################"
      end if
      write(opt_hook,"(a)")" #                                  #"
      if (lparam /= 0 .and. latom /= 0 .or. l_RC) then
        write(opt_hook,"(a)")" #          Reaction path           #"
      else
        write(opt_hook,"(a)")" #      Geometry optimization       #"
      end if
      write(opt_hook,"(a)")" #                                  #"
      write(opt_hook,"(a)")" ####################################"
      first = .false.
    end if
!
!  Geometry updated.
!
    if (abs(escf) > 1.d-30) then
      write(opt_hook,"(a,sp, d"//fmt13p6//",a)")" HEAT_OF_FORM_UPDATED:KCAL/MOL=",escf
      write(opt_hook,"(a,sp, d"//fmt13p6//",a)")" GRADIENT_NORM_UPDATED:KCAL/MOL/ANG=",gnorm

      write(opt_hook,"(a,i"//paras//",a)")" ATOM_X_UPDATED:ANGSTROMS[",3*numat, "]="
      write(opt_hook,"(3f"//fmt10p4//")") ((coord(j,i),j=1,3), i=1,numat)
      if (id > 0) then
        write(opt_hook,"(a,i1,a)")" TRANS_VECTS_UPDATED:ANGSTROMS[",id*3,"]="
        write(opt_hook,"(3f"//fmt9p4//")")((tvec(j,i),j=1,3),i=1,id)
        if (density > 1.d-1) write(opt_hook,"(a,d"//fmt13p6//",a)")" DENSITY:G/CM^3=",density
      end if
      if (id == 3 .and. nvar == 3*natoms) then
        sum = 0.d0
        do i = 1, 6
          if (abs(voigt(i)) > sum) sum = abs(voigt(i))
        end do
        if (sum /= 0.d0) then
          write(hook,"(a,i1,a)")" VOIGT_STRESS_UPDATED:GIGAPASCALS[6]="
          write(hook,"(6f"//fmt13p5//")")(voigt(i),i=1,6)
        end if
      end if
      if (nvar > 0) then
        sum = 0.d0
        do i = 1, nvar
          if (abs(grad(i)) > sum) sum = abs(grad(i))
        end do
        if (sum > 1.d-4) then
          write(hook,"(a,i"//paras//",a)")" GRADIENTS_UPDATED:KCAL/MOL/ANGSTROM[",nvar, "]="
          read(fmt9p4,'(3x,i2)')i
          j = max(int(log10(sum)),1)
          write(fmtnnp4,"(i2.2,'.',i2.2)")4 + j + i, i
          write(hook,"(10f"//fmtnnp4//")") (grad(i), i=1,nvar)
        end if
      end if
    end if
    if (opt_hook == 0) call flush (0)
    return
  else if (irc_drc) then
    if (first) then
      if (index(keywrd, " DRC") == 0) then
        irc_or_drc = "I"
        irc = .true.
      else
        irc_or_drc = "D"
        irc = .false.
      end if
      write(opt_hook,"(a)")" ####################################"
      write(opt_hook,"(a)")" #                                  #"
      write(opt_hook,"(a)")" #               "//irc_or_drc//"RC                #"
      write(opt_hook,"(a)")" #                                  #"
      write(opt_hook,"(a)")" ####################################"
      first = .false.
    end if
!
!  Geometry updated.
!
    jloop = jloop + 1
    write(opt_hook,"(a,i5.5)")" REF_POINT=",jloop
    write(opt_hook,"(a,sp, d"//fmt13p6//",a)")" MOVEMENT:ANGSTROMS=",rxn_coord
    write(opt_hook,"(a,sp, d"//fmt13p6//",a)")" POTENTIAL:KCAL/MOL=",rc_escf
    if (irc) then
      write(opt_hook,"(a,sp, d"//fmt13p6//",a)")" ENERGY_LOST:KCAL/MOL=",ekin
    else
      write(opt_hook,"(a,sp, d"//fmt13p6//",a)")" KINETIC_ENERGY:KCAL/MOL=",ekin
      write(opt_hook,"(a,sp, d"//fmt13p6//",a)")" ELAPSED_TIME:FS=",time
    end if
    write(opt_hook,"(a,i"//paras//",a)")" ATOM_X_UPDATED:ANGSTROMS[",3*numat, "]="
    write(opt_hook,"(3f"//fmt10p4//")") ((coord(j,i),j=1,3), i=1,numat)
    if (id > 0) then
      write(opt_hook,"(a,i1,a)")" TRANS_VECTS_UPDATED:ANGSTROMS[",id*3,"]="
      write(opt_hook,"(3f"//fmt9p4//")")((tvec(j,i),j=1,3),i=1,id)
      if (density > 1.d-1) write(opt_hook,"(a,d"//fmt13p6//",a)")" DENSITY:G/CM^3=",density
    end if
    if (id == 3 .and. nvar == 3*natoms) then
      sum = 0.d0
      do i = 1, 6
        if (abs(voigt(i)) > sum) sum = abs(voigt(i))
      end do
      if (sum /= 0.d0) then
        write(hook,"(a,i1,a)")" VOIGT_STRESS_UPDATED:GIGAPASCALS[6]="
        write(hook,"(6f"//fmt13p5//")")(voigt(i),i=1,6)
      end if
    end if
    if (nvar > 0) then
      sum = 0.d0
      do i = 1, nvar
        if (abs(grad(i)) > sum) sum = abs(grad(i))
      end do
      if (sum > 1.d-4) then
        write(hook,"(a,i"//paras//",a)")" GRADIENTS_UPDATED:KCAL/MOL/ANGSTROM[",nvar, "]="
        read(fmt9p4,'(3x,i2)')i
        j = max(int(log10(sum)),1)
        write(fmtnnp4,"(i2.2,'.',i2.2)")4 + j + i, i
        write(hook,"(10f"//fmtnnp4//")") (grad(i), i=1,nvar)
      end if
    end if
    if (L_MO_s) &
    call print_conventional_M_O_s(opt_hook, compressed, moa_lower, moa_upper, mob_lower, mob_upper, fmt9p4, orbs2)
    if (L_MO_s) &
    call print_M_O_data(hook, moa_lower, moa_upper, lnmoas, mob_lower, mob_upper, lnmobs, orbs, fmt9p3)
    call chrge (p, q)
    q(:numat) = tore(nat(:numat)) - q(:numat)
    write(opt_hook,"(a,i"//atoms//",a)")" ATOM_CHARGES[",numat,"]="
    write(opt_hook,"(sp,10f"//fmt9p5//")") (q(i), i=1,numat)
    if (index(keywrd, " BOND") /= 0)  then
      if (allocated(bondab)) deallocate(bondab)
      call write_screen_bonds(compressed, orbs2, hook, fmt9p4)
    end if
    if (opt_hook == 0) call flush (0)
    return
  else if (normal) then!  Calculated results common to most types of calculation
    write(hook,"(a)")" ####################################"
    write(hook,"(a)")" #                                  #"
    write(hook,"(a)")" #        Final SCF results         #"
    write(hook,"(a)")" #                                  #"
    write(hook,"(a)")" ####################################"

    write(hook,"(a,sp, d"//fmt13p6//",a)")" HEAT_OF_FORMATION:KCAL/MOL=",escf
    write(hook,"(a,sp, d"//fmt13p6//",a)")" GRADIENT_NORM:KCAL/MOL/ANGSTROM=",gnorm
    write(hook,"(a,a)")" POINT_GROUP=",name
    if (id > 0) then
      write(hook,"(a,i1,a)")" TRANS_VECTS:ANGSTROMS[",id*3,"]="
      write(hook,"(3f"//fmt9p4//")")((tvec(j,i),j=1,3),i=1,id)
      if (density > 1.d-1) write(hook,"(a,d"//fmt13p6//",a)")" DENSITY:G/CM^3=",density
      if (Abs(press(1)) > 1.d-20) then
        write(hook,"(a,i1,a)")" RESTRAINING_PRESSURE:GIGAPASCALS[3]="
        write(hook,"(3f"//fmt13p5//")")(press(i),i=1,3)
      end if
    end if
    if (id == 3 .and. nvar == 3*natoms) then
      sum = 0.d0
      do i = 1, 6
        if (abs(voigt(i)) > sum) sum = abs(voigt(i))
      end do
      if (sum /= 0.d0) then
        write(hook,"(a,i1,a)")" VOIGT_STRESS:GIGAPASCALS[6]="
        write(hook,"(6f"//fmt13p5//")")(voigt(i),i=1,6)
      end if
    end if
#ifdef MOPAC_F2003
    if (.not. ieee_is_nan(dip(4,3))) then
#else
    if (.not. isnan(dip(4,3))) then
#endif
      if (Abs(dip(4,3)) > 1.d-20) then
        write(hook,"(a,sp, d"//fmt13p6//", a)")" DIPOLE:DEBYE=",dip(4,3)
        write(hook,"(a,sp, 3d"//fmt13p5//", a)")" DIP_VEC:DEBYE[3]=",(dip(i,3), i = 1, 3)
      end if
    end if
    if (no_pKa > 0) then
      num = char(ichar("1") +int(log10(ipKa_sorted(1) + 0.05)))
      write(hook,"(a,i"//num//")")" HYDROGEN_ATOM_NO_FOR_PKA=",ipKa_sorted(1)
      write(hook,"(a,sp, d"//fmt13p6//",a)")" PKA_OF_MOST_IONIZABLE_HYDROGEN=",pKa_sorted(1)
    end if
    if (area > 1.d-6) then
      write(hook,"(a,sp, d"//fmt13p6//",a)")" AREA:SQUARE ANGSTROMS=",area
      write(hook,"(a,sp, d"//fmt13p6//",a)")" VOLUME:CUBIC ANGSTROMS=",cosvol
      if (Abs(solv_energy) > 1.d-6) write(hook,"(a,sp, d"//fmt13p6//",a)")" DIEL_ENER:EV=",solv_energy
    end if
    if (.not. mozyme) then
      if (nalpha > 0) then
        eig_min = -eigs(nalpha)
        if (nbeta > 0) eig_min = min(eig_min,-eigb(nbeta))
      else if (nelecs == 1) then
        eig_min = -eigs(1)
      else if (nelecs > 1) then
        if (nopen > 0) eig_min = -eigs(nopen)
        if (nclose > 0) eig_min = min(eig_min,-eigs(nclose))
!   CORRECTION TO I.P. OF DOUBLETS
        if (nopen - nclose == 1 .and. fract <= 1.99D0) eig_min = eig_min + 0.5D0*rjkab(1,1)
      end if
      if (nelecs > 0) write(hook,"(a,sp, d"//fmt13p6//",a)")" IONIZATION_POTENTIAL:EV=",eig_min
    end if
    write(hook,"(a,sp, d"//fmt13p6//",a)")" SPIN_COMPONENT=",sz
    write(hook,"(a,sp, d"//fmt13p6//",a)")" TOTAL_SPIN=",ss2
    num = char(ichar("1") +int(log10(nscf + 0.05)))
    write(hook,"(a,i"//num//")")" NUMBER_SCF_CYCLES=",nscf
    if (uhf) then
      num = char(ichar("1") +int(log10(nalpha + 0.05)))
      write(hook,"(a,i"//num//")")" NUM_ALPHA_ELECTRONS=",nalpha
      num = char(ichar("1") +int(log10(nbeta + 0.05)))
      write(hook,"(a,i"//num//")")" NUM_BETA_ELECTRONS=",nbeta
    end if
      sum = seconds(1) - time0
      i = int(sum*0.000001D0)
      sum = sum - i*1000000
    write(hook,"(a,sp, d"//fmt13p6//",a)")" CPU_TIME:SEC=",sum
    write(hook,"(a,sp, d"//fmt13p6//",a)")" MOLECULAR_WEIGHT:AMU=",mol_weight
    write(hook,"(a,i"//paras//",a)")" ATOM_X_OPT:ANGSTROMS[",3*numat, "]="
    write(hook,"(3f"//fmt10p4//")") ((coord(j,i),j=1,3), i=1,numat)
    if(nelecs == 0)then
      call chrge (p, q)
      q(:numat) = tore(nat(:numat)) - q(:numat)
    end if
    write(hook,"(a,i"//atoms//",a)")" ATOM_CHARGES[",numat,"]="
    write(hook,"(sp,10f"//fmt9p5//")") (q(i), i=1,numat)
    write(hook,"(a,i"//orbs//",a)")" AO_CHARGES[",norbs,"]="
    if (mozyme) then
      write(hook,"(10f"//fmt9p5//")") (p(idiag(i)), i=1,norbs)
    else
      write(hook,"(10f"//fmt9p5//")") (p((i*(i+1))/2), i=1,norbs)
      if (uhf) then
        write(hook,"(a,i"//orbs//",a)")" AO_SPINS[",norbs,"]="
        write(hook,"(sp,10f"//fmt9p5//")") (pa((i*(i+1))/2)-pb((i*(i+1))/2), i=1,norbs)
      end if
    end if
    if (nvar > 0) then
      sum = 0.d0
      do i = 1, nvar
        if (abs(grad(i)) > sum) sum = abs(grad(i))
      end do
      if (sum > 1.d-4) then
        write(hook,"(a,i"//paras//",a)")" GRADIENTS:KCAL/MOL/ANGSTROM[",nvar, "]="
        read(fmt9p4,'(3x,i2)')i
        j = max(int(log10(sum)),1)
        write(fmtnnp4,"(i2.2,'.',i2.2)")4 + j + i, i
        write(hook,"(10f"//fmtnnp4//")") (grad(i), i=1,nvar)
      end if
    end if
!
! Write out the residue letter, if it exists.
!
    do i = 1, numat
      do j = 1, size_mres
        if (index(txtatm(i)(8:11), tyres(j)) > 0) exit
      end do
      if (j <= size_mres) exit
    end do
    if (i <= numat) then
      allocate(letters(numat))
      do i = 1, numat
        do j = 1, size_mres
           if (index(txtatm(i)(8:11), tyres(j)) > 0) exit
        end do
        if (j <= size_mres) then
          letters(i) = tyr(j)
        else
          letters(i) = "?"
        end if
      end do
      write(hook,"(a,i"//atoms//",a)")" RESIDUE_LETTER[",numat,"]="
      write(hook,"(100a1)")(letters(i), i = 1, numat)
      deallocate(letters)
    end if
    if (mozyme) then
!
!   Start of MOZYME
!
!
!  Set up an array to hold the overlap matrix
!
      if (L_overlap) then
        if (compressed) then
          write(hook,"(a)")" ####################################"
          write(hook,"(a)")" #                                  #"
          write(hook,"(a)")" #    Compressed Overlap Matrix     #"
          write(hook,"(a)")" #                                  #"
          write(hook,"(a)")" ####################################"
          write(hook,"(a)")" #  Lower half triangle only"
          limit = 0.05d0
          i = min(1000000, max(10000, 50*norbs))
          allocate(comp(i), icomp(i))
          ic = 0
        else
          write(hook,"(a,i"//orbs2//",a)")" OVERLAP_MATRIX[",(norbs*(norbs + 1))/2,"]="
          write(hook,"(a)")" #  Lower half triangle only"
          allocate(overlap2(norbs,9))
        end if
        do i = 1, numat
          if (.not. compressed) overlap2 = 0.d0
          ii = iorbs(i)
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
          do j = 1, i - 1
            jj = iorbs(j)
            if (i /= j .and. ijbo (i, j) >= 0) then
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
              kl = ijbo (i, j)! starting address in H matrix  minus 1
              do i1 = 1, ii
                do j1 = 1, jj
                  sum =h(kl + (i1 - 1)* jj + j1)/(bi(i1) + bj(j1))
                  if (compressed) then
                    if (abs(sum) > limit) then
                      ic = ic + 1
                      k = nfirst(i) + i1 - 1
                      l = nfirst(j) + j1 - 1
                      icomp(ic) = (k*(k - 1))/2 + l
                      comp(ic) = sum
                    end if
                  else
                    overlap2(j1 + nfirst(j) - 1, i1) = sum
                  end if
                end do
              end do
            end if
          end do
!
! Diagonal term
!
          if (compressed) then
            do i1 = 1, ii
              k = nfirst(i) + i1 - 1
              ic = ic + 1
              icomp(ic) = (k*(k + 1))/2
              comp(ic) = 1.d0
            end do
          else
            do i1 = 1, ii
              overlap2(i1 + nfirst(j) - 1, i1) = 1.d0
            end do
          end if
!
! At this point, overlap2 contains the entire overlap matrix for atom i up to atom i
!
          if (.not. compressed) then
            do i1 = 1,ii
              write(hook,"(10f"//fmt9p4//")") (overlap2(j1,i1), j1 = 1, nfirst(i) + i1 -1)
            end do
          end if
        end do
        if (compressed) then
          if (ic < 999) then
            fmt = "3.3"
          else
            fmt = "6.6"
          end if
          j = int(log10((norbs*(norbs + 1))/2*1.0001)) + 2
          write(fmt1,'(i2,a,i2)')120/j, "i", j
          write(hook,"(1x,a,i"//fmt//",a)")" COMPRESSED_OVERLAP_INDICES[",ic,"]="
          write(hook,"("//fmt1//")") (icomp(ii),ii=1,ic)
          write(hook,"(1x,a,i"//fmt//",a)")" COMPRESSED_OVERLAP_VALUES[",ic,"]="
          write(hook,"(10f"//fmt9p4//")") (comp(ii),ii=1,ic)
        end if
      end if
      if (index(keywrd, " BOND") /= 0) then
        if (compressed) then
          write(hook,"(a)")" ####################################"
          write(hook,"(a)")" #                                  #"
          write(hook,"(a)")" #     Compressed Bond Orders       #"
          write(hook,"(a)")" #                                  #"
          write(hook,"(a)")" ####################################"
          write(hook,"(a)")" #  Lower half triangle only"
          limit = 0.05d0
          if (allocated(comp)) deallocate (comp,icomp)
          i =  50*numat
          allocate(comp(i), icomp(i))
        else
          write(hook,"(a,i"//orbs2//",a)")" BOND_ORDERS[",(numat*(numat + 1))/2,"]="
          write(hook,"(a)")" #  Lower half triangle only"
        end if
        ic = 0
        do i = 1, numat
          ii = iorbs(i)
          do j = 1, i
            jj = iorbs(j)
            if (i /= j .and. ijbo (i, j) >= 0) then
              kl = ijbo (i, j) + 1
              ku = kl + ii * jj - 1
              sum = 0.d0
              do k = kl, ku
                sum = sum + p(k) ** 2
              end do
              if (compressed) then
                if (sum > limit) then
                  ic = ic + 1
                  comp(ic) = sum
                  icomp(ic) = (i*(i - 1))/2 + j
                end if
              else
                work(j) = sum
              end if
            end if
          end do
          if (.not. compressed) then
            work(i) = 0.d0
            write(hook,"(10f"//fmt9p4//")") (work(j), j=1, i)
          end if
        end do
        if (compressed) then
          if (ic < 999) then
            fmt = "3.3"
          else
            fmt = "6.6"
          end if
          j = int(log10((numat*(numat + 1))/2*1.0001)) + 2
          write(fmt1,'(i2,a,i2)')120/j, "i", j
          write(hook,"(1x,a,i"//fmt//",a)")" COMPRESSED_BOND_ORDERS_INDICES[",ic,"]="
          write(hook,"("//fmt1//")") (icomp(ii),ii=1,ic)
          write(hook,"(1x,a,i"//fmt//",a)")" COMPRESSED_BOND_ORDERS_VALUES[",ic,"]="
          write(hook,"(10f"//fmt9p4//")") (comp(ii),ii=1,ic)
        end if
      end if! Bond orders
!
! Construct map of LMO energy levels
!
      allocate (store_eigs(norbs), eigs_map(norbs))
      store_eigs(:norbs) = eigs(:norbs)
      ii = 1
      do i = 1, noccupied
        eig_min = 1.d7
        do j = 1, noccupied
          if (eigs(j) < eig_min) then
            eig_min = eigs(j)
            ii = j
          end if
        end do
        eigs_map(i) = ii
        eigs(ii) = 1.d8
      end do
      do i = noccupied + 1, norbs
        eig_min = 1.d7
        do j = noccupied + 1, norbs
          if (eigs(j) < eig_min) then
            eig_min = eigs(j)
            ii = j
          end if
        end do
        eigs_map(i) = ii
        eigs(ii) = 1.d8
      end do
      eigs(:norbs) = store_eigs(:norbs)
      if (L_MO_s) then
        write(hook,"(a,2i8)")" SET_OF_MOS=",moa_lower, moa_upper
        if (moa_lower > 1 .or. moa_upper < norbs) write(hook,"(a)") &
        " # To print all the molecular orbitals, add keyword LARGE or use MOS=99999 in AUX  #"
        i = max(noccupied, nvirtual)
        allocate (c_lmo(norbs*i))
        eigen = (index(keywrd, " EIGEN") /= 0)
        if (eigen) then
          call lmo_to_eigenvectors(noccupied, ncf, nncf, ncocc, noccupied, &
           & icocc, icocc_dim, cocc, cocc_dim, eigs, c_lmo)
          if (compressed) then
            do i = moa_lower, noccupied
              call write_comp_vect(hook, c_lmo((i - 1)*norbs + 1), norbs, 0.002d0, .true., &
              "EIGENVECTOR_INDICES","EIGENVECTOR_COEFFICIENTS", fmt9p4)
            end do
          else
            write(hook,"(a,i"//orbs2//",a)")" EIGENVECTORS[",norbs*(moa_upper - moa_lower + 1),"]="
            write(hook,"(10f"//fmt9p4//")") ((c_lmo(j + (i - 1)*norbs), j=1,norbs), i=moa_lower,noccupied)
          end if
          call lmo_to_eigenvectors(nvirtual, nce, nnce, ncvir, nvirtual, icvir, &
           & icvir_dim, cvir, cvir_dim, eigs(noccupied + 1), c_lmo)
          if (compressed) then
            do i = 1, moa_upper - noccupied
              call write_comp_vect(hook, c_lmo((i - 1)*norbs + 1), norbs, 0.002d0, .true., &
              "EIGENVECTOR_INDICES","EIGENVECTOR_COEFFICIENTS", fmt9p4)
            end do
          else
            write(hook,"(10f"//fmt9p4//")") ((c_lmo(j + (i - 1)*norbs), j=1,norbs), i= 1, moa_upper - noccupied)
          end if
        else
          if (compressed) then
            write(hook,"(a)")" ####################################"
            write(hook,"(a)")" #                                  #"
            write(hook,"(a)")" #     Compressed LMO vectors       #"
            write(hook,"(a)")" #                                  #"
            write(hook,"(a)")" ####################################"
          else
            write(hook,"(a,i"//orbs2//",a)")" LMO_VECTORS[",norbs*(moa_upper - moa_lower +1),"]="
          end if
!
          do i_map = moa_lower, noccupied
            i = eigs_map(i_map)
!
! Expand i'th occupied LMO into uncompressed format
!
            c_lmo = 0.d0
            ii = ncocc(i)
            do j = nncf(i) + 1, nncf(i) + ncf(i)
              jj = icocc(j)
              do k = nfirst(jj), nlast(jj)
                ii = ii + 1
                c_lmo(k) = cocc(ii)
              end do
            end do
            if (compressed) then
              call write_comp_vect(hook, c_lmo, norbs, 0.002d0, .true., &
                "LMO_INDICES","LMO_COEFFICIENTS", fmt9p4)
            else
              write(hook,"(10f"//fmt9p4//")") (c_lmo(j), j=1,norbs)
            end if
          end do
          do i_map = 1, moa_upper - noccupied
             i = eigs_map(i_map + nclose) - nclose
!
! Expand i'th virtual LMO into uncompressed format
!
            c_lmo = 0.d0
            ii = ncvir(i)
            do j = nnce(i) + 1, nnce(i) + nce(i)
              jj = icvir(j)
              do k = nfirst(jj), nlast(jj)
                ii = ii + 1
                c_lmo(k) = cvir(ii)
              end do
            end do
            if (compressed) then
              call write_comp_vect(hook, c_lmo, norbs, 0.002d0, .true., &
              "LMO_INDICES","LMO_COEFFICIENTS", fmt9p4)
            else
              write(hook,"(10f"//fmt9p4//")") (c_lmo(j), j=1,norbs)
            end if
          end do
        end if
      end if

      if (L_Density) then
        if (compressed) then
          write(hook,"(a)")" ####################################"
          write(hook,"(a)")" #                                  #"
          write(hook,"(a)")" #    Compressed Density Matrix     #"
          write(hook,"(a)")" #                                  #"
          write(hook,"(a)")" ####################################"
          write(hook,"(a)")" #  Lower half triangle only"
          limit = 0.05d0
          ic = 0
        else
          write(hook,"(a,i"//orbs2//",a)")" DENSITY_MATRIX[",(norbs*(norbs + 1))/2,"]="
          write(hook,"(a)")" #  Lower half triangle only."
        end if
        if (.not. allocated(overlap2)) allocate(overlap2(norbs,9))
        if (allocated(comp)) deallocate(comp, icomp)
        i = min(1000000, max(10000, 50*norbs))
        allocate(comp(i), icomp(i))
        do i = 1, numat
          if(.not. compressed) overlap2 = 0.d0
          ii = iorbs(i)
          do j = 1, i - 1
            jj = iorbs(j)
            if (i /= j .and. ijbo (i, j) >= 0) then
              ij = ijbo (i, j)! starting address in P matrix  minus 1
              do i1 = 1, ii
                do j1 = 1, jj
                  sum = p(ij + ((i1 - 1)*i1)/2 + j1)
                  if (compressed) then
                    if (abs(sum) > limit) then
                      ic = ic + 1
                      k = nfirst(i) + i1 - 1
                      l = nfirst(j) + j1 - 1
                      icomp(ic) = (k*(k - 1))/2 + l
                      comp(ic) = sum
                    end if
                  else
                    overlap2(j1 + nfirst(j) - 1, i1) = sum
                  end if
                end do
              end do
            end if
          end do
!
! Diagonal term
!
          ij = ijbo (i, i)
          do i1 = 1, ii
            do j1 = 1, i1
              sum = p(ij + ((i1 - 1)*i1)/2 + j1)
              if (compressed) then
                if (abs(sum) > limit) then
                  ic = ic + 1
                  k = nfirst(i) + i1 - 1
                  l = nfirst(j) + j1 - 1
                  icomp(ic) = (k*(k - 1))/2 + l
                  comp(ic) = sum
                end if
              else
                overlap2(j1 + nfirst(j) - 1, i1) = sum
              end if
            end do
          end do
!
! At this point, overlap2 contains the entire density matrix for atom i up to atom i
!
          if ( .not. compressed) then
            do i1 = 1,ii
              write(hook,"(10f"//fmt9p4//")") (overlap2(j1,i1), j1 = 1, nfirst(i) + i1 -1)
            end do
          end if
        end do
        if (compressed) then
          if (ic < 999) then
            fmt = "3.3"
          else
            fmt = "6.6"
          end if
          j = int(log10((norbs*(norbs + 1))/2*1.0001)) + 2
          write(fmt1,'(i2,a,i2)')120/j, "i", j
          write(hook,"(1x,a,i"//fmt//",a)")" COMPRESSED_DENSITY_INDICES[",ic,"]="
          write(hook,"("//fmt1//")") (icomp(ii),ii=1,ic)
          write(hook,"(1x,a,i"//fmt//",a)")" COMPRESSED_DENSITY_VALUES[",ic,"]="
          write(hook,"(10f"//fmt9p4//")") (comp(ii),ii=1,ic)
        end if
      end if
      if (eigen) then
        write(hook,"(a,i"//orbs//",a)")" EIGENVALUES[",moa_upper - moa_lower + 1,"]="
      else
        if (lnmoas) &
        write(hook,"(a,i"//orbs//",a)")" LMO_ENERGY_LEVELS[",moa_upper - moa_lower + 1,"]="
      end if
      write(hook,"(10f"//fmt9p3//")") (eigs(eigs_map(i)), i=moa_lower, moa_upper)
      if (compressed) then
        if (allocated(icomp)) deallocate (icomp)
        if (allocated(comp)) deallocate (comp)
      else
        if (allocated(overlap2)) deallocate (overlap2)
      end if
      if (allocated(c_lmo)) deallocate (c_lmo)
      if (allocated(eigs_map)) deallocate (eigs_map)
      if (allocated(store_eigs)) deallocate (store_eigs)
!
!   End of MOZYME part
!
    else
!
!   Conventional matrix-algebra method
!
!
!  Set up an array to hold the overlap matrix
!
      do i = 1, norbs
        ifact(i) = (i*(i - 1))/2
      end do
      ifact(norbs+1) = (norbs*(norbs + 1))/2
      allocate(overlap((norbs*(norbs + 1))/2))
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
      if (compressed .and. L_Overlap) then
        call write_comp_vect(hook, overlap, ((norbs*(norbs + 1))/2), 0.0001d0, .false., &
        "OVERLAP_INDICES","OVERLAP_COEFFICIENTS", fmt9p4)
      else if (L_Overlap) then
        write(hook,"(a,i"//orbs2//",a)")" OVERLAP_MATRIX[",(norbs*(norbs + 1))/2,"]="
        write(hook,"(a)")" #  Lower half triangle only"
        write(hook,"(10f"//fmt9p4//")") (overlap(i), i=1,(norbs*(norbs + 1))/2)
      end if
      deallocate (overlap)

      if (L_MO_s) &
      call print_conventional_M_O_s(hook, compressed, moa_lower, moa_upper, mob_lower, mob_upper, fmt9p4, orbs2)

      if (L_Density) then
        if (uhf) then
          if (compressed) then
            call write_comp_vect(hook, pa, ((norbs*(norbs + 1))/2), 0.0005d0, .false., &
              "ALPHA_DENSITY_MATRIX_INDICES","ALPHA_DENSITY_MATRIX_COEFFICIENTS", fmt9p4)
            call write_comp_vect(hook, pb, ((norbs*(norbs + 1))/2), 0.0005d0, .false., &
              "BETA_DENSITY_MATRIX_INDICES","BETA_DENSITY_MATRIX_COEFFICIENTS", fmt9p4)
          else
            write(hook,"(a,i"//orbs2//",a)")" ALPHA_DENSITY_MATRIX[",(norbs*(norbs + 1))/2,"]="
            write(hook,"(a)")" #  Lower half triangle only"
            write(hook,"(10f"//fmt9p4//")") (pa(i), i=1,(norbs*(norbs + 1))/2)
            write(hook,"(a,i"//orbs2//",a)")" BETA_DENSITY_MATRIX[",(norbs*(norbs + 1))/2,"]="
            write(hook,"(a)")" #  Lower half triangle only"
            write(hook,"(10f"//fmt9p4//")") (pb(i), i=1,(norbs*(norbs + 1))/2)
          end if
        else
          if (compressed) then
            call write_comp_vect(hook, p, ((norbs*(norbs + 1))/2), 0.0005d0, .false., &
              "DENSITY_MATRIX_INDICES","DENSITY_MATRIX_COEFFICIENTS", fmt9p4)
          else if (L_Density) then
            write(hook,"(a,i"//orbs2//",a)")" TOTAL_DENSITY_MATRIX[",(norbs*(norbs + 1))/2,"]="
            write(hook,"(a)")" #  Lower half triangle only"
            write(hook,"(10f"//fmt9p4//")") (p(i), i=1,(norbs*(norbs + 1))/2)
          end if
        end if
      end if
      if (index(keywrd, " FOCK") /= 0) then
        if (uhf) then
          write(hook,"(a,i"//orbs2//",a)")" ALPHA_FOCK_MATRIX[",(norbs*(norbs + 1))/2,"]="
        else
          write(hook,"(a,i"//orbs2//",a)")" FOCK_MATRIX[",(norbs*(norbs + 1))/2,"]="
        end if
        write(hook,"(a)")" #  Lower half triangle only"
        write(hook,"(10f"//fmt9p4//")") (f(i), i=1,(norbs*(norbs + 1))/2)
        if (uhf) then
          write(hook,"(a,i"//orbs2//",a)")" BETA_FOCK_MATRIX[",(norbs*(norbs + 1))/2,"]="
          write(hook,"(a)")" #  Lower half triangle only"
          write(hook,"(10f"//fmt9p4//")") (fb(i), i=1,(norbs*(norbs + 1))/2)
        end if
      end if
      if (index(keywrd, " BOND") /= 0)  &
        call write_screen_bonds(compressed, orbs2, hook, fmt9p4)
      if (L_MO_s) &
      call print_M_O_data(hook, moa_lower, moa_upper, lnmoas, mob_lower, mob_upper, lnmobs, orbs, fmt9p3)
    end if
!
!   Molecular orbital occupancies
!
    if (uhf) then
      if (lnmoas) then
        write(hook,"(a,i5.5,a)")" ALPHA_MOLECULAR_ORBITAL_OCCUPANCIES[",moa_upper - moa_lower + 1,"]="
        write(hook,"(40i2)")(1, i = moa_lower, nalpha),(0, i = nalpha + 1, moa_upper)
      end if
      if (lnmobs) then
        write(hook,"(a,i5.5,a)")" BETA_MOLECULAR_ORBITAL_OCCUPANCIES[",mob_upper - mob_lower + 1,"]="
        write(hook,"(40i2)")(1, i = mob_lower, nbeta),(0, i = nbeta + 1, mob_upper)
      end if
    else
      if (nmos == 0) then
        if (lnmoas) then
!
!   Simple closed shell RHF
!
          write(hook,"(a,i5.5,a)")" MOLECULAR_ORBITAL_OCCUPANCIES[",moa_upper - moa_lower + 1,"]="
          write(hook,"(10f"//fmt7p4//")")(2.d0, i = moa_lower, nclose),(0.d0, i = nclose + 1, moa_upper)
        end if
      else
        if (lnmoas) then
!
!   RHF-C.I.: Recalculate M.O. occupancies of the active space
!
          write(hook,"(a,i5.5,a)")" MOLECULAR_ORBITAL_OCCUPANCIES[",moa_upper - moa_lower + 1,"]="
          if (nelec > 0) then
            write(hook,"(a)")" #  Below active space"
            write(hook,"(10f"//fmt7p4//")")(2.d0, i = moa_lower, nelec)
            write(hook,"(a)")" #  Active space"
            write(hook,"(10f"//fmt7p4//")")(2.d0*occa(i) + deltap(i,i), i = 1, nmos)
            if (nelec + nmos + 1 <= moa_upper) then
              write(hook,"(a)")" #  Above active space"
              write(hook,"(10f"//fmt7p4//")")(0.d0, i = nelec + nmos + 1, moa_upper)
            end if
          end if
!
!  Microstates
!
          write(hook,"(a,i"//orbs//")")" SIZE_OF_ACTIVE_SPACE=",nmos
          write(hook,"(a,i"//orbs//")")" NUMBER_MICROSTATES=",lab
          write(hook,"(a)")" #  microstate configurations are alpha first, then beta. One configuration per line"
          write(hook,"(a,i5.5,a)")" MICROSTATE_CONFIGURATIONS[",2*lab*nmos,"]="
          do i = 1, lab
            write(hook,"(25i5)")(microa(j,i),j = 1, nmos), (microb(j,i),j = 1, nmos)
          end do
          write(hook,"(a,i5)")" STATE_REQUESTED=", root_requested
          write(hook,"(a,i"//orbs//")")" STATE_DEGENERACY=",nstate
          write(hook,"(a,i5.5,a)")" STATE_VECTOR[",lab*nstate,"]="
          do i = 1, nstate
            write(hook,"(10f"//fmt9p4//")") (vectci(j + (i - 1)*lab), j=1,lab)
          end do
          write(hook,'('' STATE="'',i2,1x,3A)') state_QN, state_spin, trim(state_Irred_Rep)//'"'
          write(hook,"(a,sp, d"//fmt13p6//",a)")" STATE_ENERGY_ABSOLUTE:EV=",eig(root_requested)
          write(hook,"(a,sp, d"//fmt13p6//",a)")" STATE_ENERGY_RELATIVE:EV=",eig(root_requested) - eig(1)
        end if
      end if
    end if

    if (LMO) then
      if (uhf) then
        if (LMO_alpha) then
          LMO_alpha = .false.
          write(hook,"(a,i"//orbs2//",a)")" ALPHA_LMO_MO[",norbs*nalpha,"]="
          write(hook,"(10f"//fmt9p4//")") ((c(j,i), j=1,norbs), i=1,nalpha)
          write(hook,"(a,i"//orbs//",a)")" ALPHA_LMO_E[",nalpha,"]="
        write(hook,"(10f"//fmt9p3//")") (eigs(i), i=1,nalpha)
        else
          write(hook,"(a,i"//orbs2//",a)")" BETA_LMO_MO[",norbs*nbeta,"]="
          write(hook,"(10f"//fmt9p4//")") ((cb(j,i), j=1,norbs), i=1,nbeta)
          write(hook,"(a,i"//orbs//",a)")" BETA_LMO_E[",nbeta,"]="
        write(hook,"(10f"//fmt9p3//")") (eigb(i), i=1,nbeta)
        end if
      else
        write(hook,"(a,i"//orbs2//",a)")" LMO_MO[",norbs*nclose,"]="
        write(hook,"(10f"//fmt9p4//")") ((c(j,i), j=1,norbs), i=1,nclose)
        write(hook,"(a,i"//orbs//",a)")" LMO_E[",nclose,"]="
        write(hook,"(10f"//fmt9p3//")") (eigs(i), i=1,nclose)
      end if
    end if
  else if (reaction_path) then
    write(hook,"(a,sp, d"//fmt13p6//",a)")" HEAT_OF_FORMATION:KCAL/MOL=",escf
    write(hook,"(a,i"//paras//",a)")" ATOM_X_OPT:ANGSTROMS[",3*numat, "]="
    write(hook,"(3f"//fmt10p4//")") ((coord(j,i),j=1,3), i=1,numat)
    if (id > 0) then
      write(opt_hook,"(a,i1,a)")" TRANS_VECTS_UPDATED:ANGSTROMS[",id*3,"]="
      write(opt_hook,"(3f"//fmt9p4//")")((tvec(j,i),j=1,3),i=1,id)
      if (density > 1.d-1) write(opt_hook,"(a,d"//fmt13p6//",a)")" DENSITY:G/CM^3=",density
    end if
    if (id == 3 .and. nvar == 3*natoms) then
      sum = 0.d0
      do i = 1, 6
        if (abs(voigt(i)) > sum) sum = abs(voigt(i))
      end do
      if (sum /= 0.d0) then
        write(hook,"(a,i1,a)")" VOIGT_STRESS_UPDATED:GIGAPASCALS[6]="
        write(hook,"(6f"//fmt13p5//")")(voigt(i),i=1,6)
      end if
    end if
    if (nvar > 0) then
      sum = 0.d0
      do i = 1, nvar
        if (abs(grad(i)) > sum) sum = abs(grad(i))
      end do
      if (sum > 1.d-4) then
        write(hook,"(a,i"//paras//",a)")" GRADIENTS_UPDATED:KCAL/MOL/ANGSTROM[",nvar, "]="
        read(fmt9p4,'(3x,i2)')i
        j = max(int(log10(sum)),1)
        write(fmtnnp4,"(i2.2,'.',i2.2)")4 + j + i, i
        write(hook,"(10f"//fmtnnp4//")") (grad(i), i=1,nvar)
      end if
    end if
    if (index(keywrd, " BOND") /= 0)  &
      call write_screen_bonds(compressed, orbs2, hook, fmt9p4)
  else if (mullik) then
    write(hook,"(a,i"//atoms//",a)")" MULLIKEN_ATOM_CHARGES[",numat,"]="
    write(hook,"(sp,10f9.5)") (chrg(i), i=1,numat)
  else if (esp) then
    write(hook,"(a,i"//atoms//",a)")" ELECTROSTATIC_POTENTIAL_CHARGES[",numat,"]="
    write(hook,"(sp,10f9.5)") (q(i), i=1,numat)
  else if (loc_mos) then
    if (uhf) then
      if ((index(text,"alpha") /= 0)) then
        write(hook,"(a,i"//orbs2//",a)")" ALPHA_LMO_MO[",norbs*nalpha,"]="
        write(hook,"(10f"//fmt9p4//")") ((c(j,i), j=1,norbs), i=1,nalpha)
        write(hook,"(a,i"//orbs//",a)")" ALPHA_LMO_E[",nalpha,"]="
      write(hook,"(10f"//fmt9p3//")") (eigs(i), i=1,nalpha)
      else
        write(hook,"(a,i"//orbs2//",a)")" BETA_LMO_MO[",norbs*nbeta,"]="
        write(hook,"(10f"//fmt9p4//")") ((cb(j,i), j=1,norbs), i=1,nbeta)
        write(hook,"(a,i"//orbs//",a)")" BETA_LMO_E[",nbeta,"]="
      write(hook,"(10f"//fmt9p3//")") (eigb(i), i=1,nbeta)
      end if
    else
      write(hook,"(a,i"//orbs2//",a)")" LMO_MO[",norbs*nclose,"]="
      write(hook,"(10f"//fmt9p4//")") ((c(j,i), j=1,norbs), i=1,nclose)
      write(hook,"(a,i"//orbs//",a)")" LMO_E[",nclose,"]="
      write(hook,"(10f"//fmt9p3//")") (eigs(i), i=1,nclose)
    end if
  else if(force) then
    write(hook,"(a)")" ####################################"
    write(hook,"(a)")" #                                  #"
    write(hook,"(a)")" #      Normal mode analysis        #"
    write(hook,"(a)")" #                                  #"
    write(hook,"(a)")" ####################################"
    write(hook,"(a,i"//paras//",a)")" ORIENTATION_ATOM_X:ANGSTROMS[",3*numat, "]="
    write(hook,"(3f"//fmt10p4//")") ((coord(j,i),j=1,3), i=1,numat)
    write(hook,"(a,i"//paras//",a)")" ISOTOPIC_MASSES[",numat, "]="
    write(hook,"(10f9.4)") (atmass(i), i=1,numat)
    write(hook,"(a,3f"//fmt13p6//",a)")" ROTAT_CONSTS:CM(-1)[3]=",rot
    write(hook,"(a,3f"//fmt13p6//",a)")" PRI_MOM_OF_I:10**(-40)*GRAM-CM**2[3]=",xyzmom
    write(hook,"(a,f"//fmt13p6//",a)")" ZERO_POINT_ENERGY:KCAL/MOL=",zpe
    write(hook,"(a,i"//paras//",a)")" VIB._FREQ:CM(-1)[", nvar, "]="
    write(hook,"(10f9.2)") (freq(i), i=1,nvar)
    if (compressed) then
      do i = 1, nvar
         call write_comp_vect(hook, cnorml((i - 1)*3*numat + 1), 3*numat, 0.0005d0, .false., &
            "NORMAL_MODE_INDICES","NORMAL_MODE_COEFFICIENTS", fmt9p4)
      end do
    else
      write(hook,"(a,i"//orbs2//",a)")" NORMAL_MODES[",3*nvar*numat,"]="
      write(hook,"(10f"//fmt9p4//")") (cnorml(i),i = 1,3*nvar*numat)
    end if
    write(hook,"(a)")" #Warning: Hessian matrix is mass-weighted."
    j = (3*numat*(3*numat + 1))/2
    if (compressed) then
      write(hook,"(a)")" #COMPRESSED_HESSIAN_MATRIX:MILLIDYNES/ANGSTROM"
      call write_comp_vect(hook, fmatrx, j, 0.0001d0, .false., &
        "HESSIAN_INDICES","HESSIAN_COEFFICIENTS", fmt9p4)
    else
      write(hook,"(a,i"//orbs2//",a)")" HESSIAN_MATRIX:MILLIDYNES/ANGSTROM/SQRT(MASS(I)*MASS(J))[",j,"]="
      write(hook,"(a)")" #  Lower half triangle only"
      write(hook,"(10f"//fmt9p4//")") (fmatrx(i), i=1,j)
    end if
    if (jndex(1) == 1) then
      allocate(namo_tmp(nvar))
      namo_tmp = namo
      do i = 1,nvar
        j = index(namo_tmp(i),'"')
        if (j /= 0) namo_tmp(i) = namo(i)(:j - 1)//"''"
      end do
      write(hook,"(a,i"//paras//",a)")" NORMAL_MODE_SYMMETRY_LABELS[",nvar,"]="
      write(hook,'(10(i4,a4))') (jndex(i),namo_tmp(i), i=1,nvar)
      deallocate(namo_tmp)
    end if
    write(hook,"(a,i"//paras//",a)")" VIB._T_DIP:ELECTRONS[", nvar, "]="
    write(hook,"(10f"//fmt9p3//")") (dipt(i), i=1,nvar)
    write(hook,"(a,i"//paras//",a)")" VIB._TRAVEL:ANGSTROMS[", nvar, "]= "
    write(hook,"(10f"//fmt9p3//")") (travel(i), i=1,nvar)
    write(hook,"(a,i"//paras//",a)")" VIB._RED_MASS:AMU[", nvar, "]= "
    write(hook,"(10f"//fmt9p3//")") (redmas(i, 1), i=1,nvar)
    write(hook,"(a,i"//paras//",a)")" VIB._EFF_MASS:AMU[", nvar, "]= "
    write(hook,"(10f"//fmt9p3//")") (min(9999.999d0, max(-999.999d0, redmas(i, 2))), i=1,nvar)
    if (ilim > 0) then
      write(hook,"(a,i"//paras//",a)")" THERMODYNAMIC_PROPERTIES_TEMPS:K[", ilim, "]= "
      write(hook,"(10f"//fmt9p3//")") (T_range(i), i=1,ilim)
      write(hook,"(a,i"//paras//",a)")" ENTHALPY_TOT:CAL/MOL[", ilim, "]= "
      write(hook,"(10f"//fmt10p1//")") (H_tot(i), i=1,ilim)
      write(hook,"(a,i"//paras//",a)")" HEAT_CAPACITY_TOT:CAL/K/MOL[", ilim, "]= "
      write(hook,"(10f"//fmt10p3//")") (Cp_tot(i), i=1,ilim)
      write(hook,"(a,i"//paras//",a)")" ENTROPY_TOT:CAL/K/MOL[", ilim, "]= "
      write(hook,"(10f"//fmt10p3//")") (S_tot(i), i=1,ilim)
      write(hook,"(a,i"//paras//",a)")" H_O_F(T):KCAL/MOL[", ilim, "]= "
      write(hook,"(10f"//fmt10p2//")")   (HOF_tot(i), i=1,ilim)
    end if
    write(hook,"(a,i"//paras//",a)")" ATOM_X_FORCE:ANGSTROMS[",3*numat, "]="
    write(hook,"(3f"//fmt10p4//")") ((coord(j,i),j=1,3), i=1,numat)
  else if (int_force_const) then
    write(hook,"(a,i"//paras//",a)")" INT_FORCE_CONSTS:MILLIDYNES/ANGSTROM[",3*numat, "]= "
    write(hook,"(3f"//fmt10p4//")") ((fcint(j,i),j=1,3), i=1,numat)
    j = 0
    do i = 2, natoms
      if (na(i) /= 0) j = j + 1
    end do
    if (j /= 0) then
      write(hook,"(a,i"//paras//",a)")" INT_COORDS:ANGSTROMS AND DEGREES[",3*numat, "]="
      do i = 1, natoms
        if (na(i) == 0) then
          write(hook,"(3f"//fmt10p4//")") (geo(j,i),j=1,3)
        else
          write(hook,"(3f"//fmt10p4//")") (geo(j,i)*convert(j),j=1,3)
        end if
      end do
      write(hook,"(a,i"//paras//",a)")" CONNECTIVITY[",3*numat, "]="
      write(hook,"(3i5)")(na(i), nb(i), nc(i), i = 1, natoms)
    end if
  else if (error) then
    write(hook,"(a,a)")" ERROR=""",text(17:len_trim(text))//""""
  else if (polar) then
    write(hook,"(a,sp, d"//fmt13p6//", a)")" POLAR_FREQUENCY:EV=",omega
    write(hook,"(a,sp, d"//fmt13p6//", a)")" POLAR_AVE_ALPHA:CUBIC ANGSTROMS=", alpavg*a0**3

  else if (end_of_job) then
    if (index(text,"JOB ENDED NORM") /= 0) then
      write(hook,"(a,a)")" TERMINATION_MESSAGE=""",text(19:len_trim(text))//'"'
    else
      write(hook,"(a,a)")" ERROR_MESSAGE=""",text(19:len_trim(text))//'"'
    end if
  else if (finished) then
    write(hook,"(a,f12.2)")" CPU_TIME:SECONDS[1]=",time0
!
!  Don't print processor-independent CPU times for quick jobs - that would waste too much time.
!
    if (time0 > 1.d0 .and. job_no < 4) then
!
!  Deliberately run a time-consuming calculation to work out CPU speed
!  The "j" index is set to use up 1.0 seconds on the development computer.
!  Do NOT change this quantity!
!
      bk = seconds(2)
      bi = 0.d0
      do j = 1,142083
        do i = 1,1000
          bi = bi + 1.d0/i
        end do
      end do
      bk = seconds(2) - bk
      if(natoms == 0) write(hook,"(a,f12.2)") " DUMMY[1]=",bi ! dummy code to force bi evaluation
      write(hook,"(a,f12.2)")" CPU_TIME:ARBITRARY_UNITS[1]=",time0/bk
    end if
    write(hook,"(a)")" END OF MOPAC FILE"
    new_calcn = .true.
  end if
  return
end subroutine current_version
subroutine write_comp_vect(output, c, size_c, cutoff, norm, text1, text2, fmt9p4)
!
!  Write out the vector "c" in compressed form
!
  implicit none
  integer, intent(in) :: output, size_c
  logical, intent (in) :: norm
  double precision, intent(in) :: c(size_c), cutoff
  character (len=*), intent (in) :: text1, text2
  character, intent (in) :: fmt9p4*5
!
!  Local
!
  integer :: n, ii, j, k
  double precision :: sum, phase
  double precision, allocatable :: c_loc(:), c1(:)
  integer, allocatable :: ivec(:), ic(:)
  character :: fmt*3, fmt1*5
  allocate(c_loc(size_c), c1(size_c), ivec(size_c), ic(size_c))
  n = 0
  do j = 1, size_c
  if (Abs(c(j)) > cutoff) then
    n = n + 1
    ivec(n) = j
    c_loc(n) = c(j)
  end if
  end do
!
!  Re-sequence the coefficients so that the largest is first
!
  do ii = 1, n
  sum = 0.d0
  k = 1
  do j = 1, n
    if (sum < Abs(c_loc(j))) then
      sum = Abs(c_loc(j))
      k = j
    end if
  end do
  ic(ii) = ivec(k)
  c1(ii) = c_loc(k)
  c_loc(k) = 0.d0
  end do
  if (norm) then
    sum = 0.d0
    do j = 1, n
      sum = sum + c1(j)**2
    end do
    phase = 1.d0/sqrt(sum)
    if (c1(1) < 0.d0) phase = -phase
  else
    phase = 1.d0
  end if
!
!  Write out all finite indices.
!
  if (n < 999) then
  fmt = "3.3"
  else
  fmt = "6.6"
  end if
  j = int(log10(size_c*1.0001)) + 2
  write(fmt1,'(i2,a,i1)')120/j, "i", j
  write(output,"(1x,a,i"//fmt//",a)")trim(text1)//"[",n,"]="
  write(output,"("//fmt1//")") (ic(ii),ii=1,n)
  write(output,"(1x,a,i"//fmt//",a)")trim(text2)//"[",n,"]="
  write(output,"(10f"//fmt9p4//")") (phase*c1(ii),ii=1,n)
  return
end subroutine write_comp_vect
!
subroutine print_conventional_M_O_s(iwrite, compressed, moa_lower, moa_upper, mob_lower, mob_upper, fmt9p4, orbs2)
  use common_arrays_C, only : c, cb
  use molkst_C, only : norbs, uhf
  implicit none
  logical :: compressed
  integer :: iwrite, moa_lower, moa_upper, mob_lower, mob_upper
  character :: fmt9p4*5, orbs2*3
!
!
  integer :: i, j
  logical :: lnmos, first = .true.
  save :: first
  if (uhf) then
    if (compressed) then
      lnmos = (moa_upper - moa_lower > -1)
      if (lnmos) then
        write(iwrite,"(a,2i8)")" SET_OF_ALPHA_MOS=",moa_lower, moa_upper
        if (first .and. (moa_lower > 1 .or. moa_upper < norbs)) write(iwrite,"(a)") &
        " # To print all the molecular orbitals, add keyword LARGE or use MOS=99999 in AUX  #"
        do i = moa_lower, moa_upper
          call write_comp_vect(iwrite, c(1,i), norbs, 0.002d0, .false., &
          "ALPHA_MO_INDICES","ALPHA_MO_COEFFICIENTS", fmt9p4)
        end do
      end if
      lnmos = (mob_upper - mob_lower > -1)
      if (lnmos) then
        write(iwrite,"(a,2i8)")" SET_OF_BETA_MOS=",mob_lower, mob_upper
        if (first .and. (mob_lower > 1 .or. mob_upper < norbs)) write(iwrite,"(a)") &
        " # To print all the molecular orbitals, add keyword LARGE or use MOS=99999 in AUX  #"
        do i = mob_lower, mob_upper
          call write_comp_vect(iwrite, cb(1,i), norbs, 0.002d0, .false., &
          "BETA_MO_INDICES","BETA_MO_COEFFICIENTS", fmt9p4)
        end do
      end if
    else
      lnmos = (moa_upper - moa_lower > -1)
      if (lnmos) then
        write(iwrite,"(a,2i8)")" SET_OF_ALPHA_MOS=",moa_lower, moa_upper
        if (first .and. (moa_lower > 1 .or. moa_upper < norbs)) write(iwrite,"(a)") &
        " # To print all the molecular orbitals, add keyword LARGE or use MOS=99999 in AUX  #"
        write(iwrite,"(a,i"//orbs2//",a)")" ALPHA_EIGENVECTORS[",norbs*(moa_upper - moa_lower + 1),"]="
        write(iwrite,"(10f"//fmt9p4//")") ((c(j,i), j=1,norbs), i=moa_lower, moa_upper)
      end if
      lnmos = (mob_upper - mob_lower > -1)
      if (lnmos) then
        write(iwrite,"(a,2i8)")" SET_OF_BETA_MOS=",mob_lower, mob_upper
        if (first .and. (mob_lower > 1 .or. mob_upper < norbs)) write(iwrite,"(a)") &
        " # To print all the molecular orbitals, add keyword LARGE or use MOS=99999 in AUX  #"
        write(iwrite,"(a,i"//orbs2//",a)")" BETA_EIGENVECTORS[",norbs*(mob_upper - mob_lower + 1),"]="
        write(iwrite,"(10f"//fmt9p4//")") ((cb(j,i), j=1,norbs), i=mob_lower, mob_upper)
      end if
    end if
  else
    lnmos = (moa_upper - moa_lower > -1)
    if (.not. lnmos) return
    write(iwrite,"(a,2i8)")" SET_OF_MOS=",moa_lower, moa_upper
    if (first .and. (moa_lower > 1 .or. moa_upper < norbs)) write(iwrite,"(a)") &
    " # To print all the molecular orbitals, add keyword LARGE or use MOS=99999 in AUX  #"
    if (compressed) then
      do i = moa_lower, moa_upper
        call write_comp_vect(iwrite, c(1,i), norbs, 0.002d0, .true., &
        "MO_INDICES","MO_COEFFICIENTS", fmt9p4)
      end do
    else
      write(iwrite,"(a,i"//orbs2//",a)")" EIGENVECTORS[",norbs*(moa_upper - moa_lower + 1),"]="
      write(iwrite,"(10f"//fmt9p4//")") ((c(j,i), j=1,norbs), i=moa_lower, moa_upper)
    end if
  end if
  first = ( .not. (moa_lower > 1 .or. moa_upper < norbs .or. mob_lower > 1 .or. mob_upper < norbs))
end subroutine print_conventional_M_O_s

Subroutine print_M_O_data(hook, moa_lower, moa_upper, lnmoas, mob_lower, mob_upper, lnmobs, orbs, fmt9p3)
!
! Write out symmetry labels and eigenvalues
!
  use molkst_C, only : norbs, uhf
!
  use symmetry_C, only : jndex, namo
!
  use common_arrays_C, only : c, eigs, eigb
!
  implicit none
  integer, intent (in) :: hook, moa_lower, moa_upper, mob_lower, mob_upper
  logical, intent (in) :: lnmoas, lnmobs

  character, intent (in) :: fmt9p3*5, orbs*3
!
!  Local quantities
!
  integer :: i, j
  character, dimension(:), allocatable :: namo_tmp*4
  logical :: lnmos
  if (allocated(jndex)) then
    if (jndex(1) == 1) then
      allocate(namo_tmp(norbs))
      namo_tmp = namo
      do i = 1,norbs
        j = index(namo_tmp(i),'"')
        if (j /= 0) namo_tmp(i) = namo(i)(:j - 1)//"''"
      end do
      if (uhf) then
        lnmos = (mob_upper - mob_lower > -1)
        if (lnmos) then
          write(hook,"(a,i"//orbs//",a)")" BETA_M.O.SYMMETRY_LABELS[",mob_upper - mob_lower + 1,"]="
          write(hook,'(10(i4,a4))') (jndex(i),namo_tmp(i), i = mob_lower, mob_upper)
        end if
        lnmos = (moa_upper - moa_lower > -1)
        if (lnmos) then
        call symtrz (c, eigs, 1, .TRUE.)
        end if
        namo_tmp = namo
        do i = 1,norbs
          j = index(namo_tmp(i),'"')
          if (j /= 0) namo_tmp(i) = namo(i)(:j - 1)//"''"
        end do
        write(hook,"(a,i"//orbs//",a)")" ALPHA_M.O.SYMMETRY_LABELS[",moa_upper - moa_lower + 1,"]="
        write(hook,'(10(i4,a4))') (jndex(i),namo_tmp(i), i = moa_lower, moa_upper)
      else
        if (lnmoas) then
          write(hook,"(a,i"//orbs//",a)")" M.O.SYMMETRY_LABELS[",moa_upper - moa_lower + 1,"]="
          write(hook,'(10(i4,a4))') (jndex(i),namo_tmp(i), i = moa_lower, moa_upper)
        end if
      end if
      deallocate(namo_tmp)
    end if
  end if
  if (uhf) then
    if (lnmoas) then
      write(hook,"(a,i"//orbs//",a)")" ALPHA_EIGENVALUES[",moa_upper - moa_lower + 1,"]="
      write(hook,"(10f"//fmt9p3//")") (eigs(i), i = moa_lower, moa_upper)
    end if
    if (lnmobs) then
      write(hook,"(a,i"//orbs//",a)")" BETA_EIGENVALUES[",mob_upper - mob_lower + 1,"]="
      write(hook,"(10f"//fmt9p3//")") (eigb(i), i = mob_lower, mob_upper)
    end if
  else
    if (lnmoas) then
      write(hook,"(a,i"//orbs//",a)")" EIGENVALUES[",moa_upper - moa_lower + 1,"]="
      write(hook,"(10f"//fmt9p3//")") (eigs(i), i = moa_lower, moa_upper)
    end if
  end if
end subroutine print_M_O_data
subroutine fill_overlap_matrix(overlap)
!
!  Set up an array to hold the overlap matrix
!
  use parameters_C, only : betas, betap, betad, natorb
  use molkst_C, only : norbs, numat
  use common_arrays_C, only : nat, nfirst, h
  implicit none
  double precision, intent (out) :: overlap((norbs*(norbs + 1))/2)
  integer :: i, im1, j, ii, jj, ij, if, jf, ni, nj, norbi, norbj, ifact(norbs + 1)
  double precision :: bi(9), bj(9)
    do i = 1, norbs
      ifact(i) = (i*(i - 1))/2
    end do
    ifact(norbs+1) = (norbs*(norbs + 1))/2
    overlap = 0.d0
    do i = 1, numat
      if = nfirst(i)
      im1 = i - 1
      bi = betas(nat(i))
      ni = nat(i)
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
  end subroutine fill_overlap_matrix

  subroutine write_screen_bonds(compressed, orbs2, hook, fmt9p4)
  use common_arrays_C, only : bondab
  use molkst_C, only : numat
  use chanel_C, only : iw
  implicit none
  integer, intent (in) :: hook
  logical, intent (in) :: compressed
  character, intent (in) :: orbs2*3, fmt9p4*5
  character :: fmt*3, fmt1*5
  integer :: i, j, ic
  double precision :: sum
  double precision, allocatable :: comp(:)
  integer, allocatable :: icomp(:)
    if (.not. allocated(bondab))  then
      i = iw
      iw = 88
      open(unit=iw, status='SCRATCH')
      call bonds ()
      close(iw)
      iw = i
    end if
    if (compressed) then
      i = (numat*(numat + 1))/2
      allocate(comp(i), icomp(i))
      ic = 0
      do i = 1, numat
!
!  Work out all bond orders involving atom i
!
        do j = 1, i - 1
          sum = bondab((i*(i - 1))/2 + j)
          if (sum > 0.05d0) then
            ic = ic + 1
            comp(ic) = sum
            icomp(ic) = (i*(i - 1))/2 + j
          end if
        end do
      end do
      if (ic < 999) then
        fmt = "3.3"
      else
        fmt = "6.6"
      end if
      j = int(log10((numat*(numat + 1))/2*1.0001)) + 2
      write(fmt1,'(i2,a,i1)')120/j, "i", j
      write(hook,"(1x,a,i"//fmt//",a)")" COMPRESSED_BOND_ORDERS_INDICES[",ic,"]="
      write(hook,"("//fmt1//")") (icomp(i),i=1,ic)
      write(hook,"(1x,a,i"//fmt//",a)")" COMPRESSED_BOND_ORDERS_VALUES[",ic,"]="
      write(hook,"(10f"//fmt9p4//")") (comp(i),i=1,ic)
    else
      write(hook,"(a,i"//orbs2//",a)")" BOND_ORDERS[",(numat*(numat + 1))/2,"]="
      write(hook,"(a)")" #  Lower half triangle only"
      write(hook,"(10f"//fmt9p4//")") (bondab(i), i=1,(numat*(numat + 1))/2)
    end if
  end subroutine write_screen_bonds
