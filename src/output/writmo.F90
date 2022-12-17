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

      subroutine  writmo
      use cosmo_C, only : iseps, area, fepsi, cosvol, ediel, solv_energy
!
      use molkst_C, only : numat, nclose, nopen, fract, nalpha, nelecs, nbeta, &
      & norbs, nvar, gnorm, iflepo, enuclr,elect, ndep, nscf, numcal, escf, &
      & keywrd, os, verson, time0, moperr, last, iscf, id, pressure, mol_weight, &
      jobnam, line, mers, uhf, method_indo, &
      density, formula, mozyme, mpack, stress, &
      sz, ss2, maxtxt, E_disp, E_hb, E_hh, l_normal_html, &
      no_pKa, nalpha_open, nbeta_open, use_ref_geo, N_Hbonds, caltyp, &
      hpress, nsp2_corr, Si_O_H_corr, sum_dihed, atheat, &
      prt_gradients, prt_coords, prt_cart, prt_pops, prt_charges, pdb_label, backslash, gui
!
      use MOZYME_C, only : icocc, icvir, ncocc, ncvir, nvirtual, noccupied, &
      & nnce, nncf, cocc, cvir, ncf, nce,  cocc_dim, &
      & cvir_dim, icocc_dim, icvir_dim, iorbs
!
      use common_arrays_C, only: loc, geo, pa, pb, na, eigs, eigb, coord, &
      & dxyz, p, grad, labels, xparam, nat, nfirst, nlast, cb, q, geoa, &
      & h, c, f, tvec, pdiag, txtatm, hesinv, l_atom, ipKa_sorted, ipKa_unsorted, &
      & pKa_sorted, time_end
!
      use elemts_C, only : elemnt
!
      use maps_C, only : latom, lparam
!
      use meci_C, only : rjkab
!
      use parameters_C, only : tore, natorb
!
      use symmetry_C, only : name, state_spin, state_Irred_Rep, state_QN
!
      use funcon_C, only : fpc_10, fpc_9
!
      use chanel_C, only : iw0, iw, iarc, ibrz, brillouin_fn, archive_fn, log
!
#if MOPAC_F2003
      USE, INTRINSIC :: IEEE_ARITHMETIC
#endif
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      double precision ::  gcoord(3,numat)
      double precision, dimension (:), allocatable :: rxyz
      integer :: icalcn, i, loc11, loc21, nopn, j, k, l, m, kchrge, iwrite, mvar
      double precision :: q2(numat), degree, xreact, eionis, vol, tim, xi, sum, &
      dip, dumy(3), pKa_unsorted(numat), distortion, rms, gnorm_norm, escf_min
      logical :: ci, lprtgra, still, bcc, opend, bigcycles
      character  :: type(3)*11, idate*24, gtype*13, grtype*14, &
      flepo(19)*58, iter(2)*58, namfil*241, num*2
      character, allocatable :: old_arc_file(:)*1000
      double precision, external :: dipole, dipole_for_MOZYME, dot, meci, seconds, volume
      integer, external :: ijbo
      save type, flepo, iter, namfil, icalcn, i, bigcycles, escf_min
!***********************************************************************
!
!   WRITMO PRINTS OUT MOST OF THE RESULTS.
!         IT SHOULD NOT ALTER ANY PARAMETERS, SO THAT IT CAN BE CALLED
!         AT ANY CONVENIENT TIME.
!
!***********************************************************************
      data icalcn/ 0/
      data flepo(1), flepo(2), flepo(3)/ &
        ' 1SCF WAS SPECIFIED, SO BFGS WAS NOT USED                 ', &
        ' GRADIENTS WERE INITIALLY ACCEPTABLY SMALL                ', &
        ' HERBERTS TEST WAS SATISFIED IN BFGS                      '/
      data flepo(4), flepo(5), flepo(6)/ &
        ' THE LINE MINIMIZATION FAILED TWICE IN A ROW.   TAKE CARE!', &
        ' BFGS FAILED DUE TO COUNTS EXCEEDED. TAKE CARE!           ', &
        ' PETERS TEST WAS SATISFIED IN BFGS OPTIMIZATION           '/
      data flepo(7), flepo(8), flepo(9)/ &
        ' THIS MESSAGE SHOULD NEVER APPEAR, THERE IS A BUG IN MOPAC', &
        ' GRADIENT TEST NOT PASSED, BUT FURTHER WORK NOT JUSTIFIED ', &
        ' A FAILURE HAS OCCURRED, TREAT RESULTS WITH CAUTION!!     '/
      data flepo(10), flepo(11), flepo(12)/ &
        ' GEOMETRY OPTIMIZED USING NLLSQ. GRADIENT NORM MINIMIZED  ', &
        ' GEOMETRY OPTIMIZED USING POWSQ. GRADIENT NORM MINIMIZED  ', &
        ' CYCLES EXCEEDED, GRADIENT NOT FULLY MINIMIZED IN NLLSQ   '/
      data flepo(13), flepo(14), flepo(15)/ &
        ' 1SCF RUN AFTER RESTART.  GEOMETRY MIGHT NOT BE OPTIMIZED ', &
        ' HEAT OF FORMATION MINIMIZED IN ONE LINE SEARCH           ', &
        ' GEOMETRY OPTIMISED USING EIGENVECTOR FOLLOWING (EF).     '/
      data flepo(16), flepo(17), flepo(18)/ &
        ' NO PARAMETERS MARKED FOR OPTIMIZATION, SO 1SCF WAS USED  ', &
        ' GEOMETRY OPTIMISED USING EIGENVECTOR FOLLOWING (TS).     ', &
        ' TRANSITION STATE LOCATED USING THE "SADDLE" TECHNIQUE    '/
      data flepo(19)/ &
        ' TRANSITION GEOMETRY LOCATED USING LOCATE-TS              '/
      data iter/ ' SCF FIELD WAS ACHIEVED                                   ', &
        '  ++++----**** FAILED TO ACHIEVE SCF. ****----++++        '/
!
! SUMMARY OF RESULTS (NOTE: THIS IS IN A SUBROUTINE SO IT
!          CAN BE USED BY THE PATH OPTION)
      if (icalcn == 0) then
        namfil = '**NULL**'
        escf_min = escf
        bigcycles = (index(keywrd, " BIGCYC") /= 0)
      end if
      if (bigcycles) then
        if (escf > escf_min) return
        escf_min = escf
        rewind (iarc)
      end if
      idate = ' '
      lprtgra = (index(keywrd,' GRAD') /= 0 .and. nvar > 0)
      ci = index(keywrd,' C.I.') /= 0
      call fdate (idate)
      degree = 57.29577951308232D0
      gnorm = 0.D0
      if (nvar /= 0) gnorm = dsqrt(dot(grad,grad,nvar))
      write (iw, '(/,'' ----'',15(''-----''))')
      call wrttxt (iw)
      if (iflepo == 15 .and. index(keywrd,' TS') /=0 ) iflepo = 17
      if (iflepo == 0)                             iflepo =  7
      if (iflepo == 1 .and. nvar == 0 .and. index(keywrd,"1SCF") == 0) iflepo = 16
      write (iw, '(2/4X,A58)') flepo(iflepo)
      iscf = max(1,iscf)
      write (iw, '(4X,A58)') iter(iscf)
      write (iw, "(2/29X,A,' CALCULATION')") trim(caltyp)
      write (iw, '(55X,''MOPAC v'',a,'' '',a)') trim(verson), trim(os)
      write (iw, '(55X,A24)') idate
      if (iscf == 2) then
!
!   RESULTS ARE MEANINGLESS. DON'T PRINT ANYTHING!
!
        write (iw, '(A)') ' '
        write (iw, '(A)') ' '
        write (iw, '(A)') ' FOR SOME REASON THE SCF CALCULATION FAILED.'
        write (iw, '(A)') ' '
        write (iw, '(A)') &
          ' THE RESULTS WOULD BE MEANINGLESS, SO WILL NOT BE PRINTED.'
        write (iw, '(A)') &
          ' TRY TO FIND THE REASON FOR THE FAILURE BY USING ''PL''.'
        write (iw, '(A)') ' '
        write (iw, '(A)') &
          ' CHECK YOUR GEOMETRY AND ALSO TRY USING SHIFT OR PULAY. '
        call geout (1)
        call mopend ('THE SCF CALCULATION FAILED.')
        return
      end if
      if ( Abs (pressure) > 1.d-4) then
      !
      ! Remove energy term arising from the external pressure,
      !
        if( id == 1) then
          escf = escf + pressure * dsqrt(dot(tvec(1,1), tvec(1,1),3))
        else if (id == 3) then
          escf = escf + pressure * volume(tvec,3)
        end if
      end if
      if (use_ref_geo) then
        write (iw, &
      '(4/10X,''FINAL H.O.F PLUS STRESS ='',F17.5,'' KCAL/MOL''  ,'' ='',F14.5,'' KJ/MOL'')') escf, escf*4.184D0
        write (iw, &
      '(10X,''FINAL STRESS            ='',F17.5,'' KCAL/MOL''  ,'' ='',F14.5,'' KJ/MOL'')') stress, stress*4.184D0
        write (iw, &
      '(10X,''FINAL HEAT OF FORMATION ='',F17.5,'' KCAL/MOL''  ,'' ='',F14.5,'' KJ/MOL'')') escf - stress, (escf - stress)*4.184D0
        call  geo_diff(sum, rms, .false.)
        sum = 0.d0
        rms = 0.d0
        do i = 1, numat
          sum = sum + sqrt((geo(1,i) - geoa(1,i))**2 + &
                           (geo(2,i) - geoa(2,i))**2 + &
                           (geo(3,i) - geoa(3,i))**2)
          rms = rms +      (geo(1,i) - geoa(1,i))**2 + &
                           (geo(2,i) - geoa(2,i))**2 + &
                           (geo(3,i) - geoa(3,i))**2
        end do
        distortion = sum/numat
        write (iw, &
    '(10X,''TOTAL DISTORTION        ='',F17.5,'' Angstroms'' )') sum
        write (iw, &
    '(10X,''AVERAGE DISTORTION      ='',F17.5,'' Angstroms per atom (all atoms)'' )') distortion
        write (iw, &
    '(10X,''RMS DISTORTION          ='',F17.5,'' Angstroms per atom (all atoms)'' )') sqrt(rms/numat)
      else
        distortion = 0.d0
        rms = 0.d0
        if (index(keywrd," PM7-TS") /= 0) then
          call PM7_TS
          return
        end if
        stress = -1.1d-6
        write (iw, &
      '(4/10X,''FINAL HEAT OF FORMATION ='',F17.5,'' KCAL/MOL''  ,'' ='',F14.5,'' KJ/MOL'')') escf, escf*4.184D0
      end if
      elect = elect + solv_energy
      solv_energy = 0.d0
      if (index(keywrd," DISP") /= 0) then
        if (index(keywrd,' EPS') /= 0) then
          write(iw,'(/10x,"TOTAL ENERGY            =",f17.5,a)') &
          (elect + enuclr)*fpc_9, " KCAL/MOL = ELECTRONIC ENERGY + CORE-CORE REPULSION + SOLVATION ENERGY"
        else
          write(iw,'(/10x,"TOTAL ENERGY            =",f17.5,a)') &
          (elect + enuclr)*fpc_9, " KCAL/MOL = ELECTRONIC ENERGY + CORE-CORE REPULSION"
        end if
        write(iw,'(10x,"ENERGY OF ATOMS         =",f17.5,a)') atheat, " KCAL/MOL"
        write(iw,'(10x,"                    SUM =",f17.5,a)') &
        (elect + enuclr)*fpc_9 + atheat + solv_energy*fpc_9, " KCAL/MOL"
        if (abs(hpress) > 1.d-5)      write(iw,'(10x,"ENERGY DUE TO PRESSURE  =",f17.5,a)') hpress, " KCAL/MOL"

        write(iw,'(10x,"DISPERSION ENERGY       =", f17.5, a)') e_disp, " KCAL/MOL"
        if (E_hb < -1.d-5) write(iw,'(10x,"H-BOND ENERGY           =", f17.5, a)') e_hb, " KCAL/MOL"
        if (E_hh > 1.d-5)  write(iw,'(10x,"H - H CORRECTION ENERGY =", f17.5, a)') e_hh, " KCAL/MOL"
        if (abs(nsp2_corr) > 1.d-5)   write(iw,'(10x,"MM CORRECTION FOR >N-   =",f17.5,a)') nsp2_corr, " KCAL/MOL"
        if (abs(Si_O_H_corr) > 1.d-5) write(iw,'(10x,"MM CORR. FOR Si-O-H     =",f17.5,a)') Si_O_H_corr, " KCAL/MOL"
        if (abs(sum_dihed) > 1.d-5)   write(iw,'(10x,"MM CORR. FOR -CO-NH-    =",f17.5,a)') sum_dihed, " KCAL/MOL"
        sum = (elect + enuclr + solv_energy)*fpc_9 + atheat + hpress + nsp2_corr + Si_O_H_corr + sum_dihed + &
          e_disp + e_hb + e_hh
        write(iw,'(30x,"SUM =",f17.5,a,/)') sum, " KCAL/MOL = FINAL HEAT OF FORMATION"
        if (abs(sum - escf + stress) > 1.d-3*numat) then
          write(iw,'(5x,"*",4x,a,/)')"WARNING - An energy term is incorrect or missing!"
          write(iw, '(5x,"*",4x,''FINAL HEAT OF FORMATION     ='',F13.5,'' KCAL/MOL'')') escf - stress
          write(iw,'(5x,"*",4x,"SUM OF CONTRIBUTIONS TO HoF =",f13.5,a)') sum, " KCAL/MOL"
          write(iw,'(5x,"*",21x,"DIFFERENCE =",f13.5,a,/)') escf - stress - sum, " KCAL/MOL"
        end if
        if (N_Hbonds > 0)  write(iw,'(10x,"No. OF HYDROGEN BONDS   =", i11,7x,a)') N_Hbonds, "(H-bond Energy < -1.0 Kcal/mol)"
        if (index(keywrd, " DISP(") /= 0) then
          call l_control("PRT", len_trim("PRT"), 1)
          call post_scf_corrections(sum, .false.)
        end if
      end if
      if (numat > 1 .and. iscf == 1 .and. escf > 1.d4 .and. index(keywrd, " CHECK") == 0) &
        write(iw,'(//10x,a,//)') "Calculated Heat of Formation is very large, re-run using keyword 'CHECK'"
      if (id == 3 .or. id == 1) call write_unit_cell_HOF(iw)
      call to_screen(" Job: "//jobnam(1:len_trim(jobnam)))
      write (line,'(10x,a,f16.5,a)') "Final heat of formation = ",escf," kcal/mol"
      call to_screen(line)
      gnorm_norm =  gnorm/sqrt(1.0*numat) + 1.d-8
      write(num, '(i2)')  max(0, Int(Log10( gnorm_norm))) + 7
      write (line, '(10X,''GRADIENT NORM           ='',F17.5, '' = '',f'//num//'.5, '' PER ATOM'')') gnorm, gnorm_norm
      call to_screen(line)
      if ((id == 3 .or. id == 1) .and. iw0 > -1) call write_unit_cell_HOF(iw0)
      if (index(keywrd,' EPS') /= 0) write (iw, '(10X,A,F14.2,A)') &
        'VAN DER WAALS AREA      =', area, ' SQUARE ANGSTROMS'
      if (latom == 0) write (iw, '(/)')
      if (index(keywrd," DISP") /= 0) &
        write (iw, '(    10X,''TOTAL ENERGY            ='',F17.5,'' EV'' )') elect + enuclr + solv_energy
!
!  There is a bug in the call to symtrz.  This appears to be a compiler bug!
!
      if (norbs > 0 .and. .not. mozyme .and. numat < 200) call symtrz (c, eigs, 1, .true.)
      if (moperr) return
      if (index(keywrd," DISP") /= 0) then
        write (iw, '(10X,"ELECTRONIC ENERGY       =",F17.5," EV ")') elect
        write (iw, '(    10X,''CORE-CORE REPULSION     ='',F17.5,'' EV''    )') enuclr
        if (iseps) then
          if (abs(solv_energy) > 1.d-1) &
          write (iw, '(    10X,''SOLVATION ENERGY        ='',F17.5,'' EV''   )') solv_energy
          write (iw, '(/   10X,''DIELECTRIC ENERGY       ='',F17.5,'' EV''   )') ediel
        end if
      else
        if (iseps) then
          if (abs(solv_energy) > 1.d-1) &
          write (iw, '(    10X,''SOLVATION ENERGY        ='',F17.5,'' EV''   )') solv_energy
          write (iw, '(/   10X,''DIELECTRIC ENERGY       ='',F17.5,'' EV''   )') ediel
        end if
      end if
      if (iseps .and. Index (keywrd, "COSWRT") /= 0) call coswrt()
      hpress = 0.d0
      if (Abs (pressure) > 1.d-4) then
        if (id == 1) then
          hpress = -pressure * dsqrt (dot(tvec(1, 1), tvec(1, 1), 3))
        else if (id == 3) then
          sum =  volume (tvec, 3)
          hpress = -pressure * sum
          write (iw, '(    10X,''ENERGY DUE TO PRESSURE  ='',F17.5,'' KCAL/MOL''    )') &
          hpress
          sum = -(4184.d0*10.d0**30)*pressure/fpc_10
          if (abs(sum) > 1.d9) then
            write (iw, '(10X,''PRESSURE                ='',F17.5,'' Gp''    )') sum*1.d-9
          else
            write (iw, '(10X,''PRESSURE                ='',F17.5,'' Pascals''    )') sum
          end if
        end if
      end if
      if (Index(keywrd, " EPS") == 0 .and. id == 0 .and. (.not. mozyme .and. numat < 700)) then
        sum = 78.4d0
        fepsi = (sum-1.d0) / (sum+0.5d0)
        call cosini(.false.)
        if (moperr) then
          moperr = .false.
          area = 0.d0
          cosvol = 0.d0
        else
          call coscav
        end if
      end if
      if (moperr) then
        moperr = .false.
      else
        if (fepsi > 1.d-10 .and. area > 1.d-3) then
          write (iw, "(10X,A,F14.2,A)") "COSMO AREA              =", area, &
             & " SQUARE ANGSTROMS"
          write (iw, "(10X,A,F14.2,A)") "COSMO VOLUME            =", cosvol, &
             & " CUBIC ANGSTROMS"
        end if
      end if
      if (index(keywrd," PKA") /= 0) then
      if (allocated(ipKa_sorted))  deallocate (ipKa_sorted)
      if (allocated(ipKa_unsorted))  deallocate (ipKa_unsorted)
      if (allocated(pKa_sorted))  deallocate (pKa_sorted)
      allocate(ipKa_sorted(numat), ipKa_unsorted(numat), pKa_sorted(numat))
        call prtpka(ipKa_sorted, pKa_sorted, ipKa_unsorted, pKa_unsorted, no_pKa)
        if (moperr) return
        if (no_pKa > 0) then
          write(iw,'(/10x,a,/)')"                pKa for hydroxyl hydrogens"
          j = index(txtatm(ipKa_sorted(1)), "(")
          if (j /= 0) then
             write(iw, '(10x,a)') &
          &"        Sorted by pKa                    Sorted by atom", &
          &"Atom        Label         pKa      Atom        Label         pKa"
            do i = 1, no_pKa
              j = index(txtatm(ipKa_sorted(i)), "(")
              k = index(txtatm(ipKa_sorted(i)), ")")
              l = index(txtatm(ipKa_unsorted(i)), "(")
              m = index(txtatm(ipKa_unsorted(i)), ")")
              write(iw,"(9x,i5, a, f8.3,i9,a,f8.3)")ipKa_sorted(i), "  "//txtatm(ipKa_sorted(i))(j:k), &
                 pKa_sorted(i), ipKa_unsorted(i), "  "//txtatm(ipKa_unsorted(i))(l:m), pKa_unsorted(i)
            end do
          else
            write(iw, '(22x,a)') &
            &"Sorted by pKa       Sorted by atom", &
            &"Atom      pKa       Atom      pKa"
            do i = 1, no_pKa
              write(iw,"(21x,i5, f11.3,i9,f11.3)")ipKa_sorted(i), &
                 pKa_sorted(i), ipKa_unsorted(i), pKa_unsorted(i)
            end do
          end if
          write(iw,*)
        end if
      end if
      if (latom == 0) write (iw, '(1X)')
      if (lprtgra .or. gnorm /= 0.D0) &
        write (iw, '(10X,"GRADIENT NORM           =",F17.5, 10x, "=", 7x, f'//num//'.5, '' PER ATOM'')') &
        gnorm, gnorm/sqrt(1.0*numat)
      if (gnorm > 2.D0 .and. fract > 0.05D0 .and. fract < 1.95D0 .and. &
        index(keywrd,' 1SCF') == 0 .and. index(keywrd,' UHF') == 0 .and. &
        index(keywrd,' NOANCI') == 0) write (iw, *) &
      ' TO REDUCE GNORM FURTHER, TRY ADDING KEYWORD ''NOANCI'' AND RE-RUN THE JOB'
      still = .TRUE.
      if (latom == 0) then
        if (index(keywrd,' AIDER') == 0 .and. nvar > 0) then
#ifdef MOPAC_F2003
          if (.not. ieee_is_nan(dxyz(1)) .and. (index(keywrd,' 1SCF') == 0 .or. index(keywrd,' GRAD') /= 0)) then
#else
          if (.not. isnan(dxyz(1)) .and. (index(keywrd,' 1SCF') == 0 .or. index(keywrd,' GRAD') /= 0)) then
#endif
!
!   CHECK THAT THE CARTESIAN COORDINATE GRADIENT IS ALSO SMALL
!
            sum = dsqrt(dot(dxyz, dxyz, 3*numat))
!
! Not at a stationary point if:
! (A) "sum" is more than 2 times "gnorm" scaled from (nvar/3) to numat and
! (B) "sum" is less than 0.1 times the square-root of numat
!
! where "sum" is the gradient norm of the entire system.
!
            if (nvar > 0 .and. nvar < numat*3 - 5 .and. sum > max(2.d0, 2.D0*gnorm*sqrt((numat*3.d0)/nvar)) .and. &
              sum < sqrt(0.1d0*numat) .and. nclose == nopen .and. id == 0 .and. index(keywrd, "NOANCI") == 0) then
              if (nvar /= 1 .or. index(keywrd," GRAD") + index(keywrd, "DERIV") > 0) then
                write (iw, '(9x,A)') &
                ' WARNING -- GEOMETRY IS NOT AT A STATIONARY POINT'
                call web_message(iw,"Warning_geometry.html")
                still = .FALSE.
              end if
            end if
          end if
        end if
      else
!
!   WE NEED TO CALCULATE THE REACTION COORDINATE GRADIENT.
!
        mvar = nvar
        loc11 = loc(1,1)
        loc21 = loc(2,1)
        nvar = 1
        loc(1,1) = latom
        loc(2,1) = lparam
        xreact = geo(lparam,latom)
        gcoord = 0.d0
        call deriv (geo, gcoord)
        nvar = mvar
        loc(1,1) = loc11
        loc(2,1) = loc21
        grtype = ' KCAL/ANGSTROM'
        if (lparam == 1 .or. na(latom) == 0) then
          write (iw, &
      '(10X,''FOR REACTION COORDINATE ='',F17.5,'' ANGSTROMS'')') xreact
        else
          if (na(latom) /= 0) grtype = ' KCAL/RADIAN  '
          write (iw, &
      '(10X,''FOR REACTION COORDINATE ='',F17.5,'' DEGREES'')') xreact*degree
        end if
        write (iw, '(10X,''REACTION GRADIENT       ='',F17.5,A14)') &
        gcoord(1,1), grtype
      end if
      eionis = 0.D0
      if (nalpha > 0 .and. nbeta > 0) then
        eionis = -max(eigs(nalpha),eigb(nbeta))
      else if (nelecs == 1) then
        eionis = -eigs(1)
      else if (nelecs > 1) then
        if (nopen > 0) eionis = -eigs(nopen)
        if (nclose > 0) eionis = min(eionis,-eigs(nclose))
      end if
      i = nclose
      if (fract > 1.99D0) i = nopen
      nopn = nopen - i
      if ( .not. mozyme ) then
!   CORRECTION TO I.P. OF DOUBLETS
        if (nopn == 1) eionis = eionis + 0.5D0*rjkab(1,1)
        if (abs(eionis) > 1.D-5 .and. nopn < 2) then
          if (state_Irred_Rep /= " ") then
            write(line, '(10x, "STATE:  ",i2,1x,3A)') state_QN, state_spin, state_Irred_Rep
          else
            line = " "
          end if
          write (iw, '(  10X,"IONIZATION POTENTIAL    =  ",F16.6," EV", a)') eionis, trim(line)
        end if
        if (uhf) then
          if (nalpha >= 1) write (iw, &
          '(  10X,''ALPHA SOMO LUMO (EV)    =  '',F13.3,F7.3)') eigs(nalpha), &
          (eigs(i),i=nalpha + 1,norbs,10000)
          if (nbeta >= 1) write (iw, &
          '(  10X,''BETA  SOMO LUMO (EV)    =  '',F13.3,F7.3)') eigb(nbeta), &
          (eigb(i),i=nbeta + 1,norbs,10000)
        else if (nopen == nclose) then
          if (nopen >= 1) write (iw, &
          '(  10X,''HOMO LUMO ENERGIES (EV) =  '',F13.3,F7.3)') eigs(nopen), &
          (eigs(i),i=nopen + 1,norbs,10000)
        else if (nopn == 1) then
            if (nclose >= 1) then
            write (iw, &
      '(  10X,''HOMO (SOMO) LUMO (EV)   ='',F13.3,   '' ('',F7.3,'')'',F7.3)') &
            eigs(nclose), (eigs(i),i=nclose + 1,norbs,10000), &
            (eigs(i),i=nclose + 2,norbs,10000)
          else
            write (iw, &
      '(  10X,''     (SOMO) LUMO (EV)   ='',         ''      ('',F6.3,'')'',F7.3)') &
      eigs(nclose+1), (eigs(i),i=nclose + 2,norbs,10000)
          end if
        end if
      end if
      if (uhf) then
        write (iw, '(      10X,''NO. OF ALPHA ELECTRONS  ='',I11)') nalpha + nint((nalpha_open - nalpha)*fract)
        write (iw, '(      10X,''NO. OF BETA  ELECTRONS  ='',I11)') nbeta + nint((nbeta_open - nbeta)*fract)
      else
        write (iw, '(      10X,''NO. OF FILLED LEVELS    ='',I11)') nopen - nopn
        if (nopn /= 0) write (iw, '(   10X,''AND NO. OF OPEN LEVELS  ='',I11)') nopn
      end if
      if (mol_weight > 0.1D0) then
        if (name /= " ") then
          write (iw, '(    10X,"MOLECULAR WEIGHT        =",F16.4, 9x, "POINT GROUP:", 2x, a)') mol_weight, name
        else
          write (iw, '(    10X,"MOLECULAR WEIGHT        =",F16.4 )') mol_weight
        end if
      end if
      call gmetry (geo, coord)
      if (id == 0) call dimens (coord, iw)
      if (id == 3) then
        vol = volume(tvec,3)
        call l_control("PRT", 3, 1)
        density = mol_weight*1.D24/fpc_10/vol
        call write_cell(iw)
        if (mers(1) > 1 .or. mers(2) > 1 .or. mers(3) > 1) &
          write (iw, '(/10X,''VOLUME OF CLUSTER       ='',&
          & F17.5,'' ANGSTROMS**3 ='',f9.3, '' CM**3/MOLE'')') &
        vol,  vol*fpc_10*1.d-24
        write(iw,*)
        call write_pressure(iw)
        call write_pressure(iw0)
      end if
      if (latom == 0) write (iw, '(/)')
      write (iw, '(10X,''SCF CALCULATIONS        =   '',I8 )') nscf
      tim = seconds(1) - time0
      i = int(tim*0.000001D0)
      tim = tim - i*1000000
      call date_and_time(VALUES=time_end)
      write(iw,*)
      call timout (iw)
      if (ndep /= 0) call symtry
      do i = 1, nvar
        xparam(i) = geo(loc(2,i),loc(1,i))
      end do
      call gmetry (geo, coord)
      if (lprtgra .and. gnorm > 0.2d0) call prt_sorted_gradients
      if (prt_gradients .and. (lprtgra .or. (prt_gradients .and. gnorm > 2.d0))) then
        write (iw, '(3/7X,''FINAL  POINT  AND  DERIVATIVES'',/)')
        if (mozyme) then
          call prtgra ()
          if (pdb_label) call modgra ()
        else
          write (iw, &
      '(''   PARAMETER     ATOM    TYPE            VALUE        GRADIENT'')')
          do i = 1, nvar
            j = loc(2,i)
            k = loc(1,i)
            l = labels(k)
            xi = xparam(i)
            if (j /= 1 .and. na(k) > 0) xi = xi*degree
            if (j==1 .or. na(k) == 0) then
              gtype = 'KCAL/ANGSTROM'
            else
              gtype = 'KCAL/RADIAN  '
            end if
          if (na(k) == 0) then
            type(1) = 'CARTESIAN X'
            type(2) = 'CARTESIAN Y'
            type(3) = 'CARTESIAN Z'
          else
            type(1) = 'BOND       '
            type(2) = 'ANGLE      '
            type(3) = 'DIHEDRAL   '
          end if
            write (iw, '(I7,I11,1X,A2,4X,A11,F13.6,F14.6,2X,A13)') &
            i, k, elemnt(l), type(j), xi, grad(i), gtype
          end do
        end if
      end if
      if (prt_coords) then
!
!     WRITE OUT THE GEOMETRY
!
        write (iw, '(3/)')
        call geout (1)
      end if
      if (prt_cart) then
        if (numat > 9999) then
          write (iw, '(/28x,'' CARTESIAN COORDINATES'',/10000(/,i6,3x,a2,3x,3F16.9))') &
            (i,elemnt(nat(i)),(coord(j,i),j=1,3),i=1,numat)
        else
          write (iw, '(/28x,'' CARTESIAN COORDINATES'',/10000(/,i4,3x,a2,3x,3F16.9))') &
            (i,elemnt(nat(i)),(coord(j,i),j=1,3),i=1,numat)
        end if
      end if
      if (mozyme) then
        i = 3*numat
        if (.not. log) call bridge_H()
      else
        i = max((numat*(numat + 1))/2, 3*numat)
      end if
      allocate(rxyz(i))
      if (.not. mozyme .and. index(keywrd,'PRTINT') /= 0) then
!
!   WRITE OUT THE INTERATOMIC DISTANCES
!
        l = 0
        do i = 1, numat
          do j = 1, i
            l = l + 1
            rxyz(l) = sqrt((coord(1,i)-coord(1,j))**2+(coord(2,i)-coord(2,j))**&
              2+(coord(3,i)-coord(3,j))**2)
          end do
        end do
        write (iw, '(2/10X,''  INTERATOMIC DISTANCES'')')
        call vecprt (rxyz, numat)
      end if
      if (.not. mozyme) then
        where (.not.eigs(:norbs)>=(-999.D0) .or. .not.eigs(:norbs)<=1000.D0)
          eigs(:norbs) = 1.D-12
        end where
        where (eigb(:norbs)<(-999.D0) .or. eigb(:norbs)>1000.D0)
          eigb(:norbs) = 1.D-12
        end where
      end if
      write(iw,"(/,a,/)")formula(:len_trim(formula))
      if (index(keywrd," HESSIAN") /= 0) then
        write(iw,"(a)") "    Hessian matrix from geometry optimization"
        do i = 1, nvar
          write(iw,"(10f10.4)")(hesinv((i - 1)*nvar + j), j = 1, i)
        end do
      end if
      if (norbs*nelecs > 0 ) then
        if (id == 0 .and. .not. mozyme) &
        write (iw, '(2/,''      MOLECULAR POINT GROUP   :   '',A4)') name
        if (index(keywrd, " RAMA") /= 0) then
          call output_rama()
        end if
        if (mozyme) then
          if (index(keywrd," RE-LOC") /= 0) then
            write(iw,"(10x,a,/)")"  LMOs being Re-Localized"
            call local_for_MOZYME("OCCUPIED")
!
!  Suppress re-localization of the virtual set.  Not of interest to users.
!
!            call local_for_MOZYME("VIRTUAL")
          end if
          if ((index(keywrd,' VEC') + index(keywrd,' ALLVEC'))*nelecs /= 0) then
            if (index(keywrd, " EIGEN") /= 0) then
              call eigen(.false., .true.)
            else
              call prtlmo()
            end if
          end if
        else
          if ((index(keywrd,' VEC') + index(keywrd,' ALLVEC'))*nelecs /= 0) then
            if (uhf) then
              write (iw, '(2/10X,''ALPHA EIGENVECTORS  '')')
            else
              write (iw, '(2/16X,''EIGENVECTORS  '')')
            end if

            call matou1 (c, eigs, norbs, norbs, norbs, 2)
            if (uhf) then
              write (iw, '(2/10X,'' BETA EIGENVECTORS  '')')
              call symtrz (cb, eigb, 1, .TRUE.)
              call matou1 (cb, eigb, norbs, norbs, norbs, 2)
            end if
          else
            if (uhf) then
              write (iw, '(2/10X,''ALPHA EIGENVALUES  '')')
            else
              write (iw, '(2/18X,''EIGENVALUES  '')')
            end if
            write (iw, '(8F10.5)') (eigs(i),i=1,norbs)
            if (uhf) then
              write (iw, '(2/10X,'' BETA EIGENVALUES '')')
             write (iw, '(8F10.5)') (eigb(i),i=1,norbs)
            end if
          end if
          if (index(keywrd,' SUPER') /= 0) then
            write (iw, '(/,10X,A,/)') ' SUPERDELOCALIZABILITIES'
            call superd (c, eigs, norbs, nelecs, numat, nat)
          end if
        end if
      end if
      if (nelecs /= 0) then
!
!   Correct density matrix, if necessary
!
        if (nclose/=nopen .and. abs(fract - 2.d0) > 1.d-20 .and. &
         fract > 1.d-20 .or. index(keywrd,' C.I.') /= 0 .and. .not. method_indo) call mecip ()
        if (prt_charges) then
          write (iw, '(2/13X,'' NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS'',/)')
          i = 0
          do j = 1, numat
            i = Max(i, natorb(nat(j)))
          end do
          i = Nint(Sqrt(Real(i)))*12
          line = "s-Pop       p-Pop       d-Pop"
          if  (pdb_label) then
            write (iw, "(1x,' ATOM NO.                 TYPE                        CHARGE   ',&
               &'   No. of ELECS.   ',A)") line(1:i)
          else
            write (iw, "(1x,' ATOM NO.   TYPE          CHARGE   ',&
               &'   No. of ELECS.   ',A)") line(1:i)
          end if
        end if
        call chrge (p, q2)
        sum = 0.D0
        m = 0
        do i = 1, numat
          l = nat(i)
          q(i) = tore(l) - q2(i)
          sum = sum + q(i)
          if (.not. l_atom(i)) cycle
          j = Nint(Sqrt(Real(natorb(nat(i)) + 1.d-4)))
          k = nfirst(i)
          if (j > 0) then
            if (mozyme) then
              k = ijbo(i,i)
              if(j > 0) then
                rxyz(1) = p(k+1)
                if(j > 1) then
                  rxyz(2) = p(k+3) + p(k+6) + p(k+10)
                  if(j > 2) then
                    rxyz(3) = p(k+15) + p(k+21) + p(k+28) + p(k+36) + p(k+45)
                  end if
                end if
              end if
            else
              rxyz(1) = p((k*(k+1))/2)
              if(j > 1) then
                rxyz(2) = p(((k+1)*(k+2))/2) + p(((k+2)*(k+3))/2) &
                     & + p(((k+3)*(k+4))/2)
                if(j > 2) then
                  rxyz(3) = p(((k+4)*(k+5))/2) + p(((k+5)*(k+6))/2) &
                     & + p(((k+6)*(k+7))/2) + p(((k+7)*(k+8))/2) &
                     & + p(((k+8)*(k+9))/2)
                end if
              end if
            end if
          end if
          if (prt_charges) then
            if (pdb_label) then
              m = m + 1
              if (txtatm(m)(14:14) == "X") m = m + 1
              write (iw, "(I5,9X,A,4X,F15.6,F14.4,3F12.5)") i, elemnt (l)//"("//txtatm(m)(:maxtxt)//")", &
                 & q(i), q2(i), rxyz(1:j)
            else
              write (iw, "(I5,9X,A2,4X,F15.6,F14.4,3F12.5)") i, elemnt (l), &
                 & q(i), q2(i), rxyz(1:j)
            end if
          end if
        end do
        if (pdb_label .and. index(keywrd, "CONTROL_no_MOZYME") == 0) call modchg ()
        kchrge = nint(sum)
        if (id == 0) then
          if ( mozyme ) then
            dip = dipole_for_mozyme(dumy, 1)
          else
            dip = dipole(p,coord,dumy,1)
          end if
          sum = dumy(1) ! Dummy operation - to use dumy
        end if
      end if
      if (norbs > 0) then
        if (index(keywrd,' FOCK') /= 0) then
          write (iw, '('' FOCK MATRIX '')')
          call vecprt (f, norbs)
        end if
        if (mers(1) /= 0 .and. .not. mozyme .and. index(keywrd, " BRZ") /= 0) then
          bcc = index(keywrd,' BCC') /= 0
          open(unit=ibrz, file=brillouin_fn, status='UNKNOWN')
          write (ibrz,*) norbs, (max(mers(i),1), i = 1,3), bcc
          write (ibrz,*) (f(i),i=1,(norbs*(norbs + 1))/2)
          write (ibrz,*) tvec, id, numat, ((coord(j,i),j=1,3),i=1,numat)
          write (ibrz,*) (nfirst(i),nlast(i),i=1,numat)
        end if
        if (nelecs /= 0) then
          if (index(keywrd,' DENS') /= 0) then
            write (iw, '(2/,20X,'' DENSITY MATRIX IS '')')
            call vecprt (p, norbs)
          end if
          if (prt_pops) then
            write (iw, '(2/10X,''ATOMIC ORBITAL ELECTRON POPULATIONS'' ,/)')
            j = 0
            do i = 1, numat
              j = max(nlast(i) - nfirst(i), j)
            end do
            line = " "
            if (maxtxt == 0) then
              i = 2
            else if (pdb_label) then
              i = 16
            else
              i = maxtxt/2 + 2
            end if
            if (j == 8) then
                                            !12345678901234567890123456789012345678901234567890123456789012345678901234567890
              write(iw,'(a)')line(:i)//"   Atom"//line(:i)//"  s        px        py        pz      x^2-y^2"// &
              &"     xz        z^2       yz        xy"
            else if (j == 3) then
              write(iw,'(a)')line(:i)//"   Atom"//line(:i)//"  s        px        py        pz   "
            else
              write(iw,'(a)')line(:i)//"   Atom"//line(:i)//"  s   "
            end if
          end if
          if (prt_pops) then
            if (mozyme) then
              l = 0
              do i = 1, numat
                j = ijbo (i, i)
                do k = 1, iorbs(i)
                  l = l + 1
                  j = j + k
                  pdiag(l) = p(j)
                end do
              end do
              if (maxtxt > 1) then
                m = 0
                do i = 1, numat
                  m = m + 1
                  if (txtatm(m)(14:14) == "X") m = m + 1
                  write (iw, '(i5,1x,a, 9f10.5)') i, elemnt(nat(i))//"("//txtatm(m)(:maxtxt)//")", &
                  (pdiag(j), j = nfirst(i), nlast(i))
                end do
              else
                do i = 1, numat
                  write (iw, '(i5," ",a2, 9f10.5)') i, elemnt(nat(i)), (pdiag(j), j = nfirst(i), nlast(i))
                end do
              end if
            else
              if (maxtxt > 1) then
                m = 0
                do i = 1, numat
                  m = m + 1
                  if (txtatm(m)(14:14) == "X") m = m + 1
                  write (iw, '(2x, a, 9f10.5)') elemnt(nat(i))//"("//txtatm(m)(:maxtxt)//")", &
                  (p((j*(j+1))/2), j = nfirst(i), nlast(i))
                end do
              else
                do i = 1, numat
                  write (iw, '(i5," ",a2, 9f10.5)') i, elemnt(nat(i)), (p((j*(j+1))/2), j = nfirst(i), nlast(i))
                end do
              end if
            end if
          end if
          if (index(keywrd,' PI') /= 0) then
            write (iw, '(2/10X,''SIGMA-PI BOND-ORDER MATRIX'')')
            call denrot ()
          end if
        end if
        if (uhf) then
          sz = (nalpha + (nalpha_open - nalpha)*fract - (nbeta + (nbeta_open - nbeta)*fract))*0.5D0
          ss2 = sz*sz
          l = 0
          do i = 1, norbs
            do j = 1, i
              l = l + 1
              pa(l) = pa(l) - pb(l)
              ss2 = ss2 + pa(l)**2
            end do
            ss2 = ss2 - 0.5D0*pa(l)**2
          end do
          write (iw, '(2/20X,''(SZ)    ='',F12.6)') sz
          if (fract  < 1.d-5 .or. fract > 0.9999d0) then
            write (iw, '(  20X,''(S**2)  ='',F12.6)') ss2
          else
            write(iw,'(10x,a)')" Average over configurations used, so (S**2) is not meaningful"
          end if
          if (index(keywrd,' SPIN') /= 0) then
            write (iw, '(2/10X,''SPIN DENSITY MATRIX'')')
            call vecprt (pa, norbs)
          end if
          write (iw, '(2/10X,''ATOMIC ORBITAL SPIN POPULATIONS'',/)')
          j = 0
          do i = 1, numat
            j = max(nlast(i) - nfirst(i), j)
          end do
          if (j == 8) then
                                           !123456712345671234567123456712345671234567123456712345671234567
            write(iw,'(a)')"    Atom   Total    s     px     py     pz   x^2-y^2  xz     z^2    yz     xy"
          else if (j == 3) then
            write(iw,'(a)')"    Atom   Total    s     px     py     pz   "
          else
            write(iw,'(a)')"    Atom   Total    s   "
          end if
          write(iw,'(8x,a)')" Atomic"
          do i = 1, numat
            sum = 0.d0
            do j = nfirst(i), nlast(i)
              sum = sum + pa((j*(j+1))/2)
            end do
            write (iw, '(i5," ",a2, f8.4, 9f7.3)') i, elemnt(nat(i)), sum, (pa((j*(j+1))/2), j = nfirst(i), nlast(i))
          end do
          if (index(keywrd,' HYPERFINE') /= 0) then
!
!  WORK OUT THE HYPERFINE COUPLING CONSTANTS.
!
            write (iw, '(2/10X,''    HYPERFINE COUPLING COEFFICIENTS'' ,/)')
            j = (nalpha - 1)
            q(:numat) = pa(nfirst(:numat)*(nfirst(:numat)+1)/2)*0.3333333D0 + &
             c(nfirst(:numat), j)**2*0.66666666D0
            write (iw, '(5(2X,A2,I2,F9.5,1X))') (elemnt(nat(i)),i,q(i),i=1,&
              numat)
          end if
          pa = p - pb
        end if
        if (gui .or. index(keywrd,' BONDS') + index(keywrd,' ALLBO') /= 0) then
          if ( .not. mozyme) then
            if (nbeta == 0) then
              write (iw, '(/10X,''BONDING CONTRIBUTION OF EACH M.O.'',/)')
              call molval (c, p, 2.D0)
            else
              write (iw, '(/10X,''BONDING CONTRIBUTION OF EACH ALPHA M.O.'',/)')
              call molval (c, pa, 2.D0)
              write (iw, '(/10X,''BONDING CONTRIBUTION OF EACH BETA  M.O.'',/)')
              call molval (cb, pb, 2.D0)
            end if
          end if
          call bonds ()
        end if
        call to_screen("To_file: Normal output")
        i = nclose + nalpha
        if (index(keywrd,' LOCAL') + index(keywrd,' RABBIT') + index(keywrd,' BANANA') /= 0) then
          call local (c, i, eigs, 1, "c ")
          if (nbeta /= 0) then
            write (iw, '(2/10X,'' LOCALIZED BETA MOLECULAR ORBITALS'')')
            call local (cb, nbeta, eigb, 1, "cb")
          end if
        end if
        if (index(keywrd,' 1ELE') /= 0) then
          write (iw, '('' FINAL ONE-ELECTRON MATRIX '')')
          call vecprt (h, norbs)
        end if
        if (index(keywrd,' ENPART') /= 0) call enpart ()
      end if
      if (seconds(2) - time0 > 1.d7 .or. index(keywrd,' DENOUT') /= 0) call den_in_out(1)
      if ((ci .or. nopen /= nclose .and. Abs(fract - 2.d0) > 1.d-20 .and. fract > 1.d-20 .or. &
        index(keywrd,' SIZE') /= 0) .and. index(keywrd,' MECI') + index(keywrd,' ESR')/=0) then
        write (iw, &
      '(2/10X,''MULTI-ELECTRON CONFIGURATION INTERACTION CALCULATION'',2/)')
        last = 3
        sum = meci()
        if (moperr) return
      end if
      if (index(keywrd,' MULLIK') + index(keywrd,' GRAPH') /= 0 .and. .not. gui) then
        if (index(keywrd,' MULLIK') /= 0) write (iw, &
          '(/10X,'' MULLIKEN POPULATION ANALYSIS'')')
        if (mozyme) then
          if (allocated(c)) deallocate(c)
          allocate (c(norbs, norbs), stat=i)
          if (i /= 0) then
            call memory_error ("eigen")
            return
          end if
  !
  ! Convert LMOs to Eigenvectors
  !
          call lmo_to_eigenvectors(noccupied, ncf, nncf, ncocc, noccupied, &
           & icocc, icocc_dim, cocc, cocc_dim, eigs, c)
          call lmo_to_eigenvectors(nvirtual, nce, nnce, ncvir, nvirtual, icvir, &
           & icvir_dim, cvir, cvir_dim, eigs(noccupied + 1), c(1, noccupied + 1))
!
!              Convert "h" from MOZYME form to lower-half-triangle
!
          mpack = (norbs*(norbs + 1))/2
          allocate(pb(mpack))
          call convert_mat_packed_to_triangle(h, pb)
          deallocate(h)        !   "h" must be deallocated because it might be smaller than mpack
          allocate(h(mpack))   !   Re-allocate "h"
          h(:mpack) = pb(:mpack)
        end if
        call mullik ()
        if (index(keywrd,' GRAPH') /= 0) &
          write (iw,'(/10X,'' DATA FOR GRAPH WRITTEN TO DISK'')')
      end if
!
!  NOTE: ON EXIT FROM MULLIK, PB HOLDS THE MULLIKEN ANALYSIS.
!
      if (index(keywrd," SYB") /= 0) then
        call mpcsyb(q, kchrge, eionis, dip)
      else
        if (index(keywrd,' MULLIK') /= 0)  call mpcpop(-1)
      end if
      if (icalcn /= numcal) then
        inquire(unit=iarc, opened=opend)
        if ( .not. opend) then
          if (namfil == '**NULL**') namfil = archive_fn
          open(unit=iarc, file=archive_fn, status='UNKNOWN', position='asis', iostat = i)
          if (i /= 0) then
            write(iw,*) "Could not open archive file, run continuing."
            return
          end if
          rewind iarc
        end if
        if ( index(keywrd," HTML") + index(keywrd," PDBOUT") /= 0 .and. latom == 0 ) then
          line = archive_fn(:len_trim(archive_fn) - 3)
          line = archive_fn(:len_trim(archive_fn) - 3)
          do i = len_trim(line), 1, -1
            if (line(i:i) == "/" .or. line(i:i) == backslash) exit
          end do
          line = trim(line)//"pdb"
          open(unit=31, file=trim(line), status='UNKNOWN', position='asis')
          l_normal_html = .true.
          call l_control("Write_Escf", len("Write_Escf"), 1)
          call pdbout(31)
          close (31)
        end if
        if (numcal == 2) then
          if (index(keywrd, "OLDGEO") /= 0) then
!
! Write a warning that OLDGEO has been used, so user is aware that multiple ARC files are present
!
            rewind (iarc)
            do i = 1, 100000
              read(iarc,*,iostat = j)
              if (j /= 0) exit
            end do
            i = i - 1
            allocate( old_arc_file(i))
            rewind (iarc)
            do j = 1, i
              read(iarc,'(a)') old_arc_file(j)
            end do
            rewind (iarc)
            write(iarc,'(/10x,a)')"WARNING - This file contains multiple individual ARC files."
            write(iarc,'(18x,a)')"- Be careful to select the correct ARC file."
            write(iarc,'(a)')(trim(old_arc_file(j)), j = 1, i)
          end if
        end if
        icalcn = numcal
      end if
      if (moperr) return
      iwrite = iarc
      write (iwrite, "(2/20X,' SUMMARY OF ',A,' CALCULATION',/)", iostat = j) trim(caltyp)
      if (j /= 0) then
            write(iw,'(/10x,a)') "Could not write to archive file, run continuing."
            return
          end if
      write (iwrite, '(55X,''MOPAC v'',a,'' '',a)') trim(verson), trim(os)
      write (iwrite, '(55X,A24)') idate
      write(iwrite,"(/,a,/)")formula(:len_trim(formula))
      call wrttxt (iwrite)
      write (iwrite, '(2/4X,A58)') flepo(iflepo)
      write (iwrite, '(4X,A58)') iter(iscf)
      if (stress > -1.0d-6) then
        write (iwrite, &
      '(3/10X,''FINAL H.O.F PLUS STRESS ='',F17.5,'' KCAL/MOL''  ,'' ='',F14.5,'' KJ/MOL'')') escf, escf*4.184D0
        write (iwrite, &
      '(10X,''FINAL STRESS            ='',F17.5,'' KCAL/MOL''  ,'' ='',F14.5,'' KJ/MOL'')') stress, stress*4.184D0
        write (iwrite, &
      '(10X,''HEAT OF FORMATION       ='',F17.5,'' KCAL/MOL''  ,'' ='',F14.5,'' KJ/MOL'')') escf - stress, (escf - stress)*4.184D0
        write (iwrite, &
    '(10X,''TOTAL DISTORTION        ='',F17.5,'' Angstroms'' )') distortion*numat
        write (iwrite, &
    '(10X,''AVERAGE DISTORTION      ='',F17.5,'' Angstroms per atom (all atoms)'' )') distortion
        write (iwrite, &
    '(10X,''RMS DISTORTION          ='',F17.5,'' Angstroms per atom (all atoms)'' )') sqrt(rms/numat)
      else
        write (iwrite, &
      '(/10X,''HEAT OF FORMATION       ='',F17.5,'' KCAL/MOL''  ,'' ='',F14.5,'' KJ/MOL'')') escf, escf*4.184D0
      end if
      if (index(keywrd," DISP") /= 0) then
        if (index(keywrd,' EPS') /= 0) then
          write(iwrite,'(/10x,"TOTAL ENERGY            =",f17.5,a)') &
          (elect + enuclr)*fpc_9, " KCAL/MOL = ELECTRONIC ENERGY + CORE-CORE REPULSION + SOLVATION ENERGY"
        else
          write(iwrite,'(/10x,"TOTAL ENERGY            =",f17.5,a)') &
          (elect + enuclr)*fpc_9, " KCAL/MOL = ELECTRONIC ENERGY + CORE-CORE REPULSION"
        end if
        write(iwrite,'(10x,"ENERGY OF ATOMS         =",f17.5,a)') atheat, " KCAL/MOL"
        if (index(keywrd,' EPS') /= 0) then
          if (abs(solv_energy) > 1.d-1) then
            write (iwrite, '(    10X,''SOLVATION ENERGY        ='',F17.5,'' KCAL/MOL''   )') solv_energy*fpc_9
          else if (.not. mozyme) then
            write (iwrite, '(/10X,"SOLVATION ENERGY CAN ONLY BE PRINTED WHEN MOZYME IS USED",/)')
          end if
        end if
        write(iwrite,'(10x,"                    SUM =",f17.5,a)') &
        (elect + enuclr)*fpc_9 + atheat + solv_energy*fpc_9, " KCAL/MOL"
        if (abs(hpress) > 1.d-5)      write(iwrite,'(10x,"ENERGY DUE TO PRESSURE  =",f17.5,a)') hpress, " KCAL/MOL"
        write(iwrite,'(10x,"DISPERSION ENERGY       =", f17.5, a)') e_disp, " KCAL/MOL"
        if (E_hb < -1.d-5) write(iwrite,'(10x,"H-BOND ENERGY           =", f17.5, a)') e_hb, " KCAL/MOL"
        if (E_hh > 1.d-5)  write(iwrite,'(10x,"H - H CORRECTION ENERGY =", f17.5, a)') e_hh, " KCAL/MOL"
        if (abs(nsp2_corr) > 1.d-5)   write(iwrite,'(10x,"MM CORRECTION FOR >N-   =",f17.5,a)') nsp2_corr, " KCAL/MOL"
        if (abs(Si_O_H_corr) > 1.d-5) write(iwrite,'(10x,"MM CORR. FOR Si-O-H     =",f17.5,a)') Si_O_H_corr, " KCAL/MOL"
        if (abs(sum_dihed) > 1.d-5)   write(iwrite,'(10x,"MM CORR. FOR -CO-NH-    =",f17.5,a)') sum_dihed, " KCAL/MOL"
        sum = (elect + enuclr)*fpc_9 + atheat + hpress + solv_energy*fpc_9 + nsp2_corr + Si_O_H_corr + sum_dihed + &
          e_disp + e_hb + e_hh
        write(iwrite,'(30x,"SUM =",f17.5,a,/)') sum, " KCAL/MOL = FINAL HEAT OF FORMATION"
        if (abs(sum - escf + stress) > 1.d-3*numat) then
          write(iwrite,'(5x,"*",4x,a,/)')"WARNING - An energy term is incorrect or missing!"
          write(iwrite,'(5x,"*",4x,''FINAL HEAT OF FORMATION     ='',F13.5,'' KCAL/MOL'')') escf - stress
          write(iwrite,'(5x,"*",4x,"SUM OF CONTRIBUTIONS TO HoF =",f13.5,a)') sum, " KCAL/MOL"
          write(iwrite,'(5x,"*",21x,"DIFFERENCE =",f13.5,a,/)') escf - stress - sum, " KCAL/MOL"
        end if
        if (N_Hbonds > 0)  write(iwrite,'(10x,"No. OF HYDROGEN BONDS   =", i11,7x,a, /)') &
          N_Hbonds, "(H-bond Energy < -1.0 Kcal/mol)"
      end if
      if (numat > 1 .and. iscf == 1 .and. escf > 1.d4 .and. index(keywrd, " CHECK") == 0) &
        write(iwrite,'(//10x,a,//)') "Calculated Heat of Formation is very large, re-run using keyword 'CHECK'"
      if (id == 3 .or. id == 1) call write_unit_cell_HOF(iwrite)
      if (index(keywrd," DISP") /= 0) then
        write (iwrite, '(    10X,''TOTAL ENERGY            ='',F17.5,'' EV'' )') elect + enuclr + solv_energy
        write (iwrite, '(  10X,''ELECTRONIC ENERGY       ='',F17.5,'' EV'')')  elect
        write (iwrite, '(  10X,''CORE-CORE REPULSION     ='',F17.5,'' EV'')')  enuclr
        if (abs(solv_energy) > 1.d-1) &
          write (iwrite, '(    10X,''SOLVATION ENERGY        ='',F17.5,'' EV''   )') solv_energy
      end if
      if (abs(ediel) > 1.d-5) then
        write (iwrite, '(    10X,''DIELECTRIC ENERGY       ='',F17.5,'' EV''   )') ediel
      end if
      if (Abs (pressure) > 1.d-4) then
        if (id == 1) then
          else if (id == 3) then
          sum =  volume (tvec, 3)
          write (iwrite, '(    10X,''ENERGY DUE TO PRESSURE  ='',F17.5,'' KCAL/MOL''    )') &
          hpress
          sum = -(4184.d0*10.d0**30)*pressure/fpc_10
          if (abs(sum) > 1.d9) then
            write (iwrite, '(10X,''PRESSURE                ='',F17.5,'' Gp''    )') sum*1.d-9
          else
            write (iwrite, '(10X,''PRESSURE                ='',F17.5,'' Pascals''    )') sum
          end if
        end if
      end if
      if (id == 3) then
        call l_control("PRT", 3, 1)
        call write_cell(iwrite)
        write(iwrite,'(" ")')
        if (mers(1) > 1 .or. mers(2) > 1 .or. mers(3) > 1) &
          write (iwrite, '(/10X,''VOLUME OF CLUSTER       ='',&
          & F17.5,'' ANGSTROMS**3 ='',f9.3, '' CM**3/MOLE'')') &
          vol,  vol*fpc_10*1.d-24
        write(iwrite,*)
      end if
      if (lprtgra .or. gnorm /= 0.D0) write (iwrite, &
        '(  10X, "GRADIENT NORM           =",F17.5, 10x, "=", 7x, f'//num//'.5, " PER ATOM")') gnorm, gnorm/sqrt(1.0*numat)
      if (latom == 0) then
        if (.not.still) then
          write (iwrite, '(9x,A)') &
          ' WARNING -- GEOMETRY IS NOT AT A STATIONARY POINT'
          call web_message(iwrite,"Warning_geometry.html")
        end if
      else
        grtype = ' KCAL/ANGSTROM'
        if (lparam == 1) then
          write (iwrite, &
      '(    10X,''FOR REACTION COORDINATE =  '',F16.4,'' ANGSTROMS'')')&
             xreact
        else
          if (na(latom) /= 0) grtype = ' KCAL/RADIAN  '
          write (iwrite, '(    10X,''FOR REACTION COORDINATE ='',F16.4,'' DEGREES'')') &
            xreact*degree
        end if
        write (iwrite, '(    10X,''REACTION GRADIENT       ='',F18.6,A14     )'&
          ) gcoord(1,1), grtype
      end if
      if(id == 0) then
        if(name /= " ") then
          write (iwrite, '(  10X,''DIPOLE                  ='',F17.5, " DEBYE   POINT GROUP:", 2x,A)') dip, name
        else
          write (iwrite, '(  10X,''DIPOLE                  ='',F17.5,  '' DEBYE   '')') dip
        end if
      end  if
      if (uhf) then
        write (iwrite, '(  10X,''(SZ)                    =  '',F16.6)') sz
        if (fract  < 1.d-5 .or. fract > 0.9999d0) then
            write (iwrite, '(  10X,''(S**2)                  =  '',F16.6)') ss2
          else
            write(iwrite,'(10x,a)')" Average over configurations used, so (S**2) is not meaningful"
          end if
        write (iwrite, '(  10X,''NO. OF ALPHA ELECTRONS  =  '',I9)') nalpha + nint((nalpha_open - nalpha)*fract)
        write (iwrite, '(  10X,''NO. OF BETA  ELECTRONS  =  '',I9)') nbeta + nint((nbeta_open - nbeta)*fract)
      else
        write (iwrite, '(      10X,''NO. OF FILLED LEVELS    =  '',I9, 3a)') nopen - nopn
        if (nopn /= 0) write (iwrite, &
          '(  10X,''AND NO. OF OPEN LEVELS  =  '',I9)') nopn
      end if
      if (ci .or. nopen /= nclose .and. Abs(fract - 2.d0) > 1.d-20 .and. fract > 1.d-20) &
         write (iwrite, '(  10X,''CONFIGURATION INTERACTION WAS USED'')')
      if (kchrge /= 0) write (iwrite, &
        '(  SP, 10X,''CHARGE ON SYSTEM        =  '',I9)') kchrge
      if ( .not. mozyme ) then
        if (state_Irred_Rep /= " ") then
          write(line, '(11x, "STATE:  ",i2,1x,3A)') state_QN, state_spin, state_Irred_Rep
        else
          line = " "
        end if
        write (iwrite, &
        '(  10X,"IONIZATION POTENTIAL    =  ",F16.6," EV", a)') eionis, trim(line)
        if (uhf) then
          if (nalpha >= 1) write (iwrite, &
          '(  10X,''ALPHA SOMO LUMO (EV)    =  '',F13.3,F7.3)') eigs(nalpha), &
          (eigs(i),i=nalpha + 1,norbs,10000)
          if (nbeta >= 1) write (iwrite, &
          '(  10X,''BETA  SOMO LUMO (EV)    =  '',F13.3,F7.3)') eigb(nbeta), &
          (eigb(i),i=nbeta + 1,norbs,10000)
        else if (nopen == nclose) then
          if (nopen >= 1) write (iwrite, &
          '(  10X,''HOMO LUMO ENERGIES (EV) =  '',F13.3,F7.3)') eigs(nopen), &
          (eigs(i),i=nopen + 1,norbs,10000)
          else if (nopn == 1) then
            if (nclose >= 1) then
            write (iwrite, &
      '(  10X,''HOMO (SOMO) LUMO (EV)   ='',F13.3,   '' ('',F7.3,'')'',F7.3)') &
            eigs(nclose), (eigs(i),i=nclose + 1,norbs,10000), &
            (eigs(i),i=nclose + 2,norbs,10000)
          else
            write (iwrite, &
      '(  10X,''     (SOMO) LUMO (EV)   ='',         ''      ('',F6.3,'')'',F7.3)') &
      eigs(nclose+1), (eigs(i),i=nclose + 2,norbs,10000)
          end if
        end if
      end if
      if (index(keywrd," PKA") /= 0 .and. no_pKa > 0 ) then
        write(iwrite,'(/24x,a,/)')"    pKa for hydroxyl hydrogens"
        j = index(txtatm(ipKa_sorted(1)), "(")
        if (j /= 0) then
          write(iwrite, '(10x,a)') &
          &"        Sorted by pKa                    Sorted by atom", &
          &"Atom        Label         pKa      Atom        Label         pKa"
          do i = 1, no_pKa
            j = index(txtatm(ipKa_sorted(i)), "(")
            k = index(txtatm(ipKa_sorted(i)), ")")
            l = index(txtatm(ipKa_unsorted(i)), "(")
            m = index(txtatm(ipKa_unsorted(i)), ")")
            write(iwrite,"(9x,i5, a, f8.3,i9,a,f8.3)")ipKa_sorted(i), "  "//txtatm(ipKa_sorted(i))(j:k), &
               pKa_sorted(i), ipKa_unsorted(i), "  "//txtatm(ipKa_unsorted(i))(l:m), pKa_unsorted(i)
          end do
        else
         write(iwrite, '(22x,a)') &
            &"Sorted by pKa       Sorted by atom", &
            &"Atom      pKa       Atom      pKa"
            do i = 1, no_pKa
              write(iwrite,"(21x,i5, f11.3,i9,f11.3)")ipKa_sorted(i), &
                 pKa_sorted(i), ipKa_unsorted(i), pKa_unsorted(i)
            end do
        end if
        write(iwrite,*)
      end if
      write (iwrite, '(  10X,''MOLECULAR WEIGHT        =  '',F14.4)') mol_weight
      if (id == 3) then
    !    call write_cell(iwrite)
        call write_pressure(iwrite)
      end if
      if (fepsi > 1.d-10 .and. area > 1.d-3) then
        write (iwrite, "(10X,A,F14.2,A)") "COSMO AREA              =", area, &
             & " SQUARE ANGSTROMS"
        write (iwrite, "(10X,A,F14.2,A)") "COSMO VOLUME            =", cosvol, &
             & " CUBIC ANGSTROMS"
      end if
      if (id == 0) call dimens (coord, iwrite)
      write (iwrite, '(  10X,''SCF CALCULATIONS        =  '',I9)') nscf
      call timout (iwrite)
      if (index(keywrd," PRTCHAR") /= 0) then
        write (iwrite, '(2/10X,''FINAL GEOMETRY OBTAINED'',36X,''CHARGE'')')
      else
        write (iwrite, '(2/10X,''FINAL GEOMETRY OBTAINED'')')
      end if
      call geout (iwrite)
      if (index(keywrd,' AIGOUT') /= 0) then
        write (iwrite, '(2/,A)') '  GEOMETRY IN GAUSSIAN Z-MATRIX STYLE'
        call wrttxt (iwrite)
        call geoutg (iwrite)
      end if
   !   if (in_house_only) call write_paras_used(iwrite)
      nscf = 0
      return
      end subroutine writmo
!
      subroutine empiri ()
      use molkst_C, only : numat, formula
      use Common_arrays_C, only : nat
      use chanel_C, only : log, ilog
      implicit  none
      integer :: i, j, k, nlab
      character  :: elemnt(107)*2, line*120
      integer, dimension (50) :: llab, mlab
      intrinsic Char, Ichar, Int, Log10
      data elemnt / "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", &
     & "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", "Sc", "Ti", &
     & "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", &
     & "Kr", "Rb", "Sr", "Y ", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", &
     & "Cd", "In", "Sn", "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", "Pr", &
     & "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", &
     & "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", &
     & "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U ", "Np", "Pu", "Am", "Cm", &
     & "Bk", "Mi", "XX", "+3", "-3", "Cb", "++", "+ ", "--", "- ", "Tv" /
     !
     !    The order of elements in the empirical formula is C H N O, then
     !    whatever else there is, in no particular order.  (The order of the
     !    remaining elements is the order in which they occur in the
     !    molecule.
     !
      nlab = 4
      llab = 0
      llab(1) = 6
      llab(2) = 1
      llab(3) = 7
      llab(4) = 8
      do i = 1, 4
        mlab(i) = 0
      end do
      do i = 1, numat
        j = nat(i)
        do k = 1, nlab
          if (j == llab(k)) go to 1000
        end do
        nlab = nlab + 1
        llab(nlab) = j
        mlab(nlab) = 1
        cycle
  1000  mlab(k) = mlab(k) + 1
      end do
      if (llab(nlab) == 0) nlab = nlab - 1
      j = 0
      do i = 1, nlab
        if (mlab(i) /= 0) then
          j = j + 1
          mlab(j) = mlab(i)
          llab(j) = llab(i)
        end if
      end do
      nlab = j
     !
     !   BUILD THE FORMAT STATEMENT
      line = "(10X,A,1X,"
      j = 11
      do i = 1, nlab
        if (i /= 1) then
          j = j + 1
          line (j:j) = ","
        end if
        j = j + 1
        !                       CHARACTER FORMAT
        line (j:j) = "A"
        k = 2
        if (elemnt(llab(i)) (2:2) == " ") then
          k = 1
        end if
        j = j + 1
        line (j:j) = Char (k+Ichar ("0"))
        j = j + 1
        !                       NUMBER FORMAT
        if (mlab(i) > 1) then
           !***********************************************************************
          line (j:j+5) = ",I" // Char (Int(Log10(1.0*mlab(i)))+Ichar ("1")) // ",1X"
          j = j + 5
        else
          mlab(i) = 0
  !
  !  If formula goes over nline characters, exit.  The system is obviously artificial
  !
        if ( j + 4 > 400 ) return
          line (j:j+4) = ",I1.0"
          j = j + 4
        end if
      end do
      j = j + 1
      line (j:j+2) = ")"
      write (formula, line) " Empirical Formula:", (elemnt(llab(i)), mlab(i), i=1, nlab)
      write(formula(len_trim(formula) + 1:),"(a,i6,a)")"  =",numat," atoms"
      if (log) write(ilog, line(1:len_trim(line) - 1)//",a,i6,a)") " Empirical Formula:", &
        (elemnt(llab(i)), mlab(i), i=1, nlab),"  =",numat," atoms"
  end subroutine empiri
!
  subroutine coswrt()
!
  use cosmo_C, only : area, ediel, cosvol, fepsi, amat, nsetf, &
    nps, qscnet, srad, qscat, arat, cosurf, phinet, iatsp, solv_energy
!
  use chanel_C, only : iwc, cosmo_fn, ir
!
  use molkst_C, only : keywrd, numat, escf, elect, enuclr, line, mozyme
!
  use common_arrays_C, only : nat, coord
!
  use funcon_C, only : a0, ev
!
  use linear_cosmo, only : am1dft_solve
!
  implicit none
  double precision :: sum
  integer :: i, j, k, i1, i2, ips, ix
  double precision :: fcon, dqam1dft, dd, enew
    open (unit=iwc, file=cosmo_fn, status='unknown')
    if (index(keywrd,'COSCCH') > 0) then
      fcon=a0*ev
10    read(ir,'(a)',err=30,end=30) line
        if(index(line,'dq:') == 0) goto 20
        read(line(5:),*) i1,i2,dqam1dft
        do ips=1,nps
          dd=0.d0
          do ix=1,3
           dd=dd+(cosurf(ix,ips)-(coord(ix,i1)+coord(ix,i2))/2)**2
          end do
          phinet(ips,3)=phinet(ips,3)-dqam1dft/sqrt(dd)
        end do
        goto 10
20      do i=1,numat
          read(line,*,err=30) dqam1dft
          do ips=1,nps
            dd=0.d0
            do ix=1,3
             dd=dd+(cosurf(ix,ips)-coord(ix,i))**2
            end do
            phinet(ips,3)=phinet(ips,3)-dqam1dft/sqrt(dd)
          end do
          read(ir,'(a)',err=30,end=30) line
        end do
30    continue
      if (mozyme) then
        qscnet(1:nps, 3) = qscnet(1:nps, 3) / (-fepsi)
        call am1dft_solve (qscnet(:, 3), phinet(:, 3), nps)
      else
        call coscl2(amat, nsetf, qscnet(1,3), phinet(1,3), nps)
      end if
!
! SCALE QSCNUC AND CALCULATE INTERACTION ENERGY
!
      enew = 0.d0
      do i = 1, nps
        qscnet(i,3) =  -fepsi*qscnet(i,3)
        enew = enew + qscnet(i,3)*phinet(i,3)
      end do
      ediel = fcon*enew*0.5d0
    end if
    write (iwc, "(a)") keywrd (1:80)
    write (iwc, "(a)") keywrd (81:160)
    write (iwc, "(a)") keywrd (161:240)
    sum = 0.d0
    do i = 1, nps
      sum = sum + qscnet(i, 3)
    end do
    write (iwc, "(/10x,'FINAL HEAT OF FORMATION =',f17.5,' KCAL/MOL',' = ', f12.5,' KJ/MOL')") &
       escf, escf * 4.184d0
    write (iwc, "(    10X,'TOTAL ENERGY            =',F17.5,      ' EV')") &
       & elect + enuclr + solv_energy
    write (iwc, "(    10X,'DIELECTRIC ENERGY       =',F17.5,      ' EV')") &
       & ediel
    write (iwc, "(10X,A,F14.2,A)") "COSMO AREA              =", area, " SQUARE ANGSTROMS"
    write (iwc, "(10X,A,F14.2,A)") "COSMO VOLUME            =", cosvol, " CUBIC ANGSTROMS"
    write (iwc, "(    10X,'TOTAL COSMO CHARGE      =',F17.5,a,/)") &
       & sum
    call web_message(iwc,"COSWRT.html")
    write (iwc, "(10X,A)") "ATOMIC DATA "
    write (iwc, "(2x,a10,4x,a24,10x,a13,2x,a12,A7,5x,a8)") " NR. ELEM.", "COORDINATES", &
       & "RADIUS  ", "COSMO-CHARGE", "AREA", "SIGMA"
    do i = 1, numat
     write (iwc, "(2i5,4x,7F12.6)") i, nat (i), &
       (coord(j, i), j=1, 3), srad(i), qscat(i), arat(i), qscat(i) / Max (arat(i), 1.d-10)
    end do
    write (iwc, "(/10X,A,i8)") " SEGMENT DATA: NPS=", nps
    write (iwc, "(a15,a36,5A12)") " NR. ATOM ELEM.", &
      & "COORDINATES (X, Y, Z)     ", "COSMO-CHARGE", " AREA ", &
      & "SIGMA", "POTENTIAL"
    do i = 1, nps
      k=0
      do j=1,3
        sum=cosurf(j,i)
        if(log10(abs(sum)).gt.4) k=max(k,int(log10(abs(sum))-3))
        if(sum.lt.0.0.and. log10(abs(sum)).gt.3) &
&                   k=max(k,int(log10(abs(sum))-2))
      end do
      if(k.eq.0) write (iwc, "(3i5,7F12.6)") i, iatsp (i), nat (iatsp(i)), &
&       (cosurf(j, i), j=1, 3), qscnet(i, 3), cosurf(4, i), &
&        qscnet(i, 3) / cosurf(4, i), phinet(i, 3)
      if(k.eq.1) write (iwc, "(3i5,7F12.5)") i, iatsp (i), nat (iatsp(i)), &
&       (cosurf(j, i), j=1, 3), qscnet(i, 3), cosurf(4, i), &
&        qscnet(i, 3) / cosurf(4, i), phinet(i, 3)
      if(k.eq.2) write (iwc, "(3i5,7F12.4)") i, iatsp (i), nat (iatsp(i)), &
&       (cosurf(j, i), j=1, 3), qscnet(i, 3), cosurf(4, i), &
&        qscnet(i, 3) / cosurf(4, i), phinet(i, 3)
    end do
  end subroutine coswrt
  subroutine write_paras_used(iwrt)
    use Common_arrays_C, only : nat
    USE parameters_C, only : guess1, guess2, guess3, gpp, gp2, hsp, gss, gsp, betas, betap, betad, &
      zs, zp, zd, uss, upp, udd, zsn, zpn, zdn, pocord, alpb, xfac, &
      f0sd, g2sd, main_group, alp
    use molkst_C, only : numat
    implicit none
    integer, intent(in) :: iwrt
!
!  Local
!
    integer :: i, j, k
    logical :: els(107) = .false.
    character :: elemnt(107)*2, num1*1
!-----------------------------------------------
      data (elemnt(i),i=1,107)/ 'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O '&
        , 'F ', 'NE', 'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', &
        'CA', 'SC', 'TI', 'V ', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA'&
        , 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y ', 'ZR', 'NB', 'MO', &
        'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I ', 'XE'&
        , 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', &
        'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA', 'W ', 'RE', 'OS', 'IR'&
        , 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', &
        'AC', 'TH', 'PA', 'U ', 'NP', 'PU', 'AM', 'CM', 'BK', 'MI', 'XX', 'FM'&
        , 'MD', 'CB', '++', '+', '--', '-', 'TV'/
      do i = 1, numat
        els(min(99,max(1,nat(i)))) = .true.
      end do
      write (iwrt, "(//,16X,A)") " Parameters used"
      write (iwrt, "(/5X, ' Parameter Type  Element    Parameter')")
      do i = 1, 98
        if (.not. els(i)) cycle
        if (uss(i) /= 0.D0)   write (iwrt, "(10x,A7,9x,A2,F16.6)") "USS  ", elemnt(i), uss(i)
        if (upp(i) /= 0.D0)   write (iwrt, "(10x,A7,9x,A2,F16.6)") "UPP  ", elemnt(i), upp(i)
        if (udd(i) /= 0.D0)   write (iwrt, "(10x,A7,9x,A2,F16.6)") "UDD  ", elemnt(i), udd(i)
        if (betas(i) /= 0.D0) write (iwrt, "(10x,A7,9x,A2,F16.6)") "BETAS", elemnt(i), betas(i)
        if (betap(i) /= 0.D0) write (iwrt, "(10x,A7,9x,A2,F16.6)") "BETAP", elemnt(i), betap(i)
        if (betad(i) /= 0.D0) write (iwrt, "(10x,A7,9x,A2,F16.6)") "BETAD", elemnt(i), betad(i)
        if (zs(i) /= 0.D0)    write (iwrt, "(10x,A7,9x,A2,F16.6)") "ZS   ", elemnt(i), zs(i)
        if (zp(i) /= 0.D0)    write (iwrt, "(10x,A7,9x,A2,F16.6)") "ZP   ", elemnt(i), zp(i)
        if (zd(i) /= 0.D0)    write (iwrt, "(10x,A7,9x,A2,F16.6)") "ZD   ", elemnt(i), zd(i)
        if (zsn(i) /= 0.D0)   write (iwrt, "(10x,A7,9x,A2,F16.6)") "ZSN  ", elemnt(i), zsn(i)
        if (zpn(i) /= 0.D0)   write (iwrt, "(10x,A7,9x,A2,F16.6)") "ZPN  ", elemnt(i), zpn(i)
        if (zdn(i) /= 0.D0)   write (iwrt, "(10x,A7,9x,A2,F16.6)") "ZDN  ", elemnt(i), zdn(i)
        if (alp(i) /= 0.D0)   write (iwrt, "(10x,A7,9x,A2,F16.6)") "ALP  ", elemnt(i), alp(i)
        if (gss(i) /= 0.D0)   write (iwrt, "(10x,A7,9x,A2,F16.6)") "GSS  ", elemnt(i), gss(i)
        if (gsp(i) /= 0.D0)   write (iwrt, "(10x,A7,9x,A2,F16.6)") "GSP  ", elemnt(i), gsp(i)
        if (gpp(i) /= 0.D0)   write (iwrt, "(10x,A7,9x,A2,F16.6)") "GPP  ", elemnt(i), gpp(i)
        if (gp2(i) /= 0.D0)   write (iwrt, "(10x,A7,9x,A2,F16.6)") "GP2  ", elemnt(i), gp2(i)
        if (hsp(i) /= 0.D0)   write (iwrt, "(10x,A7,9x,A2,F16.6)") "HSP  ", elemnt(i), hsp(i)
        if (pocord(i) /= 0.D0)write (iwrt, "(10x,A7,9x,A2,F16.6)") "POC  ", elemnt(i), pocord(i)
        if (.not. main_group(i) .and. f0sd(i) /= 0.D0) write (iwrt, "(10x,A7,9x,A2,F16.6)") "F0SD ", elemnt(i), f0sd(i)
        if (.not. main_group(i) .and. g2sd(i) /= 0.D0) write (iwrt, "(10x,A7,9x,A2,F16.6)") "G2SD ", elemnt(i), g2sd(i)
        if (i /= 99) then
          do j = 1, 4
            if (guess1(i,j) /= 0.D0) write (iwrt, '(12X,a3,I1,10x,a2,f16.6)')"FN1", j, elemnt(i), guess1(i,j)
            if (guess2(i,j) /= 0.D0) write (iwrt, '(12X,a3,I1,10x,a2,f16.6)')"FN2", j, elemnt(i), guess2(i,j)
            if (guess3(i,j) /= 0.D0) write (iwrt, '(12X,a3,I1,10x,a2,f16.6)')"FN3", j, elemnt(i), guess3(i,j)
          end do
        else
          do j = 1, 4
            if (guess1(i,j) /= 0.D0) then
              k = 1 + 3*(j - 1)
              if (k < 10) then
                num1 = "1"
              else
                num1 = "2"
              end if
              write (iwrt, &
                '(6X,"data gues1(",I3,",",I1,")/",          F17.6,"D0/ ! = PAR",i'//num1//')')&
                 i, j, guess1(i,j),k
            end if
            if (guess2(i,j) /= 0.D0) then
              k = 2 + 3*(j - 1)
              if (k < 10) then
                num1 = "1"
              else
                num1 = "2"
              end if
              write (iwrt, &
                '(6X,"data gues2(",I3,",",I1,")/",          F17.6,"D0/ ! = PAR",i'//num1//')')&
                 i, j, guess2(i,j),k
            end if
            if (guess3(i,j) /= 0.D0) then
              k = 3 + 3*(j - 1)
              if (k < 10) then
                num1 = "1"
              else
                num1 = "2"
              end if
              write (iwrt, &
                '(6X,"data gues3(",I3,",",I1,")/",          F17.6,"D0/ ! = PAR",i'//num1//')')&
                 i, j, guess3(i,j),k
            end if
          end do
        end if
        do j = 1, i
          if (.not.els(j)) cycle
          if (Abs(alpb(i,j)) > 0.01d0) write(iwrt,"(12x,a,a2,7x,a,f16.6)") "ALPB_",elemnt(j), elemnt(i), alpb(i,j)
          if (Abs(xfac(i,j)) > 0.01d0) write(iwrt,"(12x,a,a2,7x,a,f16.6)") "XFAC_",elemnt(j), elemnt(i), xfac(i,j)
        end do
      end do
    end subroutine write_paras_used


    subroutine PM7_TS()
    use common_arrays_C, only : grad, xparam
    use molkst_C, only : escf, moperr, mpack
    use chanel_C, only : iw
    implicit none
    integer :: i
    i = mpack
    call moldat(1)
    moperr = .false. ! An error in moldat is not important here.
    mpack = i
    call calpar
    call compfg (xparam, .TRUE., escf, .TRUE., grad, .FALSE.)
     write (iw, '(4/10X,''FINAL HEAT OF FORMATION ='',F17.5,'' KCAL/MOL''  ,'' ='',F14.5,'' KJ/MOL'')') escf, escf*4.184D0
    return
  end subroutine PM7_TS
