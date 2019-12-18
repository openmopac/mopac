      subroutine deriv(geo, gradnt)  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      USE molkst_C, ONLY: numat, norbs, nclose, nopen, fract, natoms, numcal, &
      & ndep, nvar, keywrd, cosine, moperr, mpack, isok, id, l123, line, &
      pressure, l1u, l2u, l3u, method_PM7, method_pm8
      use common_arrays_C, only : dxyz, loc, errfn, aicorr, tvec, nat
      USE symmetry_C, ONLY: locpar, idepfn 
      USE chanel_C, only : iw, ir, job_fn
      use elemts_C, only : elemnt
      use funcon_C, only : pi
      use derivs_C, only : aidref, work2
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  17:57:06  03/12/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use upcase_I 
      use mopend_I 
      use symtry_I 
      use gmetry_I 
      use dernvo_I 
      use dcart_I 
      use dfield_I 
      use volume_I
      use jcarin_I 
      use mxm_I
#if GPU
      use mxv_I
#endif
!              
      use dot_I 
      use geout_I 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double), intent(in)  :: geo(3,natoms) 
      real(double), intent(inout)  :: gradnt(nvar) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: idelta, icalcn, i, j, nstep, nw2, ncol
      real(double), dimension(3) :: change 
      real(double) :: coord(3*natoms), gold(3*natoms), xparam(3*natoms)


      real(double) :: grlim, sum, gnorm, step, press, summ, press1, &
      press2, press3
      logical :: scf1, halfe, slow, aifrst, debug, precis, intn, geochk, ci, &
        aic, noanci, field, saddle, DH_correction, large

      save change, scf1, halfe, idelta, slow, icalcn, aifrst, debug, &
        precis, intn, geochk, ci, aic, grlim, nw2, field, DH_correction, large
!-----------------------------------------------
!***********************************************************************
!
!    DERIV CALCULATES THE DERIVATIVES OF THE ENERGY WITH RESPECT TO THE
!          INTERNAL COORDINATES. THIS IS DONE BY FINITE DIFFERENCES.
!
!    THE MAIN ARRAYS IN DERIV ARE:
!        LOC    INTEGER ARRAY, LOC(1,I) CONTAINS THE ADDRESS OF THE ATOM
!               INTERNAL COORDINATE LOC(2,I) IS TO BE USED IN THE
!               DERIVATIVE CALCULATION.
!        GEO    ARRAY \GEO\ HOLDS THE INTERNAL COORDINATES.
!        gradnt   ON EXIT, CONTAINS THE DERIVATIVES
!
!***********************************************************************  
      data icalcn/ 0/  
      if (icalcn /= numcal) then          
        aifrst = index(keywrd,' RESTART') == 0 
        saddle = index(keywrd, " SADDLE") /= 0
        debug = index(keywrd,' DERIV') /= 0 
        field = index(keywrd,' FIELD') /= 0 
        large = index(keywrd,'LARGE') /= 0 
        precis = index(keywrd,' PREC') /= 0 
        DH_correction = (index(keywrd,' PM6-D') + index(keywrd,' PM6-H') /= 0 .or. &
          method_PM7 .or. method_pm8 )   
        intn = index(keywrd,'  XYZ') == 0 
        if (saddle) then
          nw2 = max(mpack, 9*natoms**2)
        else
          nw2 = max(mpack, 6*natoms*l123)
        end if
        errfn = 0.d0
        if (allocated(aidref)) deallocate(aidref)
        if (allocated(work2))  deallocate(work2)
        allocate(aidref(nvar), work2(nw2))
        aidref = 0.d0
!
!   GEOCHK is true if the system is a transition state in internal
!          coordinates with all coordinates marked for optimization.
!          It can be over-ridden by GEO-OK
!
        geochk = index(keywrd,' TS') + index(keywrd,' NLLSQ') + index(keywrd,&
          ' SIGMA') /= 0 
        geochk = geochk .and. intn .and. nvar>=numat*3-6 .and. id==0 .and. &
          index(keywrd,'GEO-OK')==0 
        geochk = geochk .and. index(keywrd,' XYZ')==0 
        ci = index(keywrd,' C.I.') /= 0 
        scf1 = index(keywrd,' 1SCF') /= 0 
        aic = index(keywrd,'AIDER') /= 0 
        if (aic .and. aifrst) then  
          open(unit=ir, file=job_fn, status='OLD', blank='ZERO',position='asis') 
          rewind ir 
!
!  ISOK IS SET FALSE: ONLY ONE SYSTEM ALLOWED
!
          isok = .FALSE. 
          do i = 1, 1000 
            read (ir, '(A)') line 
            call upcase (line, 80) 
            if (index(line,'AIDER') == 0) cycle  
            exit  
          end do 
          do j = 1, 1000 
            read (ir, '(A)', end=40, err=40) line 
            call upcase (line, 80) 
            if (index(line,'AIDER') /= 0) go to 60 
          end do 
   40     continue 
          write (iw, '(2/,A)') ' KEYWORD "AIDER" SPECIFIED, BUT NOT' 
          write (iw, '(A)') ' PRESENT AFTER Z-MATRIX.  JOB STOPPED' 
          call mopend (&
         'KEYWORD "AIDER" SPECIFIED, BUT NOT PRESENT AFTER Z-MATRIX.  JOB STOPPED')
          return  
   50     continue 
          write (iw, '(2/,A)') '  FAULT IN READ OF AB INITIO DERIVATIVES' 
          write (iw, '(A)') '  DERIVATIVES READ IN ARE AS FOLLOWS' 
          write (iw, '(6F12.6)') (aidref(j),j=1,i) 
          call mopend ('FAULT IN READ OF AB INITIO DERIVATIVES')  
          return  
   60     continue 
          if (natoms > 2) then 
            j = 3*natoms - 6 
          else 
            j = 1 
          endif 
          read (ir, *, end=50, err=50) (aidref(i),i=1,j) 
          write (iw, '(/,A,/)') &
            ' AB-INITIO DERIVATIVES IN KCAL/MOL/(ANGSTROM OR RADIAN)' 
          write (iw, '(5F12.6)') (aidref(i),i=1,j) 
          do i = 1, nvar 
            if (loc(1,i) > 3) then 
              j = 3*loc(1,i) + loc(2,i) - 9 
            else if (loc(1,i) == 3) then 
              j = loc(2,i) + 1 
            else 
              j = 1 
            endif 
            aidref(i) = aidref(j) 
          end do 
          write (iw, '(/,A,/)') ' AB-INITIO DERIVATIVES FOR VARIABLES' 
          write (iw, '(5F12.6)') (aidref(i),i=1,nvar) 
          if (ndep /= 0) then 
            do i = 1, nvar 
              sum = aidref(i) 
              aidref(i) = aidref(i) + count(loc(1,i)==locpar(:ndep) .and. &
              (loc(2,i)==idepfn(:ndep) .or. loc(2,i)==3 .and. &
              idepfn(:ndep)==14))*sum 
            end do 
            write (iw, '(/,A,/)') &
              ' AB-INITIO DERIVATIVES AFTER SYMMETRY WEIGHTING' 
            write (iw, '(5F12.6)') (aidref(j),j=1,nvar) 
          endif 
          close(ir, status='KEEP') 
        endif 
        grlim = 0.01D0 
        if (precis) grlim = 0.0001D0 
        halfe = nopen>nclose .and. Abs(fract - 2.d0) > 1.d-20 .and. Abs(fract) > 1.d-20 .or. ci 
        idelta = -7 
!
!   IDELTA IS A MACHINE-PRECISION DEPENDANT INTEGER
!
        change(1) = 10.D0**idelta 
        change(2) = 10.D0**idelta 
        change(3) = 10.D0**idelta 
!
!    CHANGE(I) IS THE STEP SIZE USED IN CALCULATING THE DERIVATIVES.
!    FOR "CARTESIAN" DERIVATIVES, CALCULATED USING DCART,AN
!    INFINITESIMAL STEP, HERE 0.000001, IS ACCEPTABLE. IN THE
!    HALF-ELECTRON METHOD A QUITE LARGE STEP IS NEEDED AS FULL SCF
!    CALCULATIONS ARE NEEDED, AND THE DIFFERENCE BETWEEN THE TOTAL
!    ENERGIES IS USED. THE STEP CANNOT BE VERY LARGE, AS THE SECOND
!    DERIVITIVE IN FLEPO IS CALCULATED FROM THE DIFFERENCES OF TWO
!    FIRST DERIVATIVES. CHANGE(1) IS FOR CHANGE IN BOND LENGTH,
!    (2) FOR ANGLE, AND (3) FOR DIHEDRAL.
!
      endif 
      if (nvar == 0) return  
      if (debug) then 
        write (iw, '('' GEO AT START OF DERIV'')') 
        call geout(-iw)
      endif 
      gnorm = 0.D0 
      do i = 1, nvar 
        gold(i) = gradnt(i) 
        xparam(i) = geo(loc(2,i),loc(1,i)) 
        gnorm = gnorm + gradnt(i)**2 
      end do 
      gnorm = sqrt(gnorm) 
      slow = .FALSE. 
      noanci = .FALSE. 
      if (halfe) then 
        noanci = (index(keywrd,'NOANCI') /= 0 .or. nopen == norbs)
        slow = (noanci .and. (gnorm<grlim .or. scf1))
      else
        slow = (index(keywrd,'NOANCI') /= 0)
      endif 
      if (ndep /= 0) call symtry 
      call gmetry (geo, coord) 
!
!  COORD NOW HOLDS THE CARTESIAN COORDINATES
!
      if (halfe .and. .not.noanci .and. numat > 1) then 
        if (debug) write (iw, *) 'DOING ANALYTICAL C.I. DERIVATIVES' 
        call dernvo () 
        if (moperr) return  
      else 
        if (debug) write (iw, *) 'DOING VARIATIONALLY OPTIMIZED DERIVATIVES' 
        call dcart (coord, dxyz)        
      endif  
      if (DH_correction) call post_scf_corrections(sum, .true.)
      if (DH_correction .and. debug) then  
        write (iw, '(2/10X,''CARTESIAN COORDINATE DERIVATIVES'')') 
        write (iw, '(11X,a)')"(Includes post-SCF corrections)" 
        write (iw, '(/3X,       ''NUMBER  ATOM '',5X,''X'',12X,''Y'',12X,''Z'',/)') 
        if (l123 == 1) then 
          write (iw, '(I6,4x,a2,F13.6,2F13.6)') (i,elemnt(nat(i)),(dxyz(j + (i - 1)*3),j=1,3),i=1,numat*l123) 
        else if (large) then 
          write (iw, '(I6,4x,a2,F13.6,2F13.6)') (i,elemnt(nat((i - 1)/l123+1)), &
            (dxyz(j + (i - 1)*3),j=1,3),i=1,numat*l123) 
        else 
          write (iw, '(I6,4x,a2,F13.6,2F13.6)') (i,elemnt(nat((i-1)/l123+1)),(dxyz(j + (i-1)*3) + &
            dxyz(j + i*3) + dxyz(j + (i + 1)*3),j=1,3),i=1, numat*l123 - 2,3) 
        endif
      end if 
!
!   THE CARTESIAN DERIVATIVES ARE IN DXYZ
!
      if (field) call dfield () 
      if (Abs (pressure) > 1.d-4) then 
        if (id == 1) then
      !
      !  Add in gradient tension contribution
      !
          i = 3 * l123 / 2 
        !
        !  For polymers, pressure is pull in Newtons for 1 mole
        !  1N = J/M = 10**(-3)/4.184 kcal/M = 4.184*10**(-3)*10**(-10) kcal/Angstrom
          press = pressure / Sqrt (dot(tvec(1, 1), tvec(1, 1), 3))
          do j = 1, 3
            dxyz(j + i) = dxyz(j + i) + tvec(j, 1) * press
            dxyz(j + i + 3) = dxyz(j + i + 3) - tvec(j, 1) * press
          end do
        else if (id == 3) then
          summ = volume (tvec, 3)
          press1 = summ / dot (tvec(1, 1), tvec(1, 1), 3) * pressure
          press2 = summ / dot (tvec(1, 2), tvec(1, 2), 3) * pressure
          press3 = summ / dot (tvec(1, 3), tvec(1, 3), 3) * pressure
  !
  ! Add in gradient pressure term
  !
  !  First, the central unit cell
  !
          i = 3*(l1u * (2*l2u+1) * (2*l3u+1) + l2u * (2*l3u+1) + l3u - 1)
          do j = 1, 3
            dxyz(j + i) = dxyz(j + i) + tvec(j, 1) * press1 
            dxyz(j + i) = dxyz(j + i) + tvec(j, 2) * press2
            dxyz(j + i) = dxyz(j + i) + tvec(j, 3) * press3
          end do
  !
  !  Cell in 0,0,1 position
  !
          i = 3*(l1u * (2*l2u+1) * (2*l3u+1) + l2u * (2*l3u+1) + l3u)
          do j = 1, 3
            dxyz(j + i) = dxyz(j + i) - tvec(j, 3) * press3
          end do
  !
  !  Cell in 0,1,0 position
  !
          i = 3*(l1u * (2*l2u+1) * (2*l3u+1) + (l2u+1) * (2*l3u+1) + l3u - 1)
          do j = 1, 3
            dxyz(j + i) = dxyz(j + i) - tvec(j, 2) * press2
          end do
  !
  !  Cell in 1,0,0 position
  !
          i = 3*((l1u+1) * (2*l2u+1) * (2*l3u+1) + l2u * (2*l3u+1) + l3u - 1)
          do j = 1, 3
            dxyz(j + i) = dxyz(j + i) - tvec(j, 1) * press1
          end do
        end if
      end if

      step = change(1) 
      nstep = nw2/(3*numat*l123)  
      do i = 1, nvar, nstep 
        j = min(i + nstep - 1,nvar) 
        call jcarin (xparam, step, precis, work2, ncol, i, j) 
#if GPU
        call mxv (work2, j - i + 1, dxyz, ncol, gradnt(i))         
#else       
        call mxm (work2, j - i + 1, dxyz, ncol, gradnt(i), 1)   
#endif       
               
!        
      end do 
      if (precis) then 
        step = 0.5D0/step 
      else 
        step = 1.0D0/step 
      endif 
      gradnt(:nvar) = gradnt(:nvar)*step 
!
!  NOW TO ENSURE THAT INTERNAL DERIVATIVES ACCURATELY REFLECT CARTESIAN
!  DERIVATIVES
!
      if (geochk) then 
        sum = dot(gradnt,gradnt,nvar) 
        if (sum<2.D0 .and. dot(dxyz,dxyz,3*numat)>max(4.D0,sum*4.D0)) then 
!
! OOPS, LOOKS LIKE AN ERROR.
!
          do i = 1, nvar 
            j = int(xparam(i)/3.141D0) 
            if (.not.(loc(2,i)==2 .and. loc(1,i)>3 .and. abs(xparam(i)-j*pi)<0.005D0)) cycle  
!
!  ERROR LOCATED, BUT CANNOT CORRECT IN THIS RUN
!
            write (iw, '(2/,3(A,/),I3,A)') &
              ' INTERNAL COORDINATE DERIVATIVES DO NOT REFLECT', &
              ' CARTESIAN COORDINATE DERIVATIVES', &
              ' TO CORRECT ERROR, INCREASE DIHEDRAL OF ATOM', loc(1,i), &
              ' BY 90 DEGREES' 
            write (iw, '(2/,A)') '     CURRENT GEOMETRY' 
            call geout (iw)
            call mopend (&
       ' INTERNAL COORDINATE DERIVATIVES DO NOT REFLECT CARTESIAN COORDINATE DERIVATIVES') 
            return  
          end do 
        endif 
      endif 
!
!  THIS CODE IS ONLY USED IF THE KEYWORD NOANCI IS SPECIFIED
      if (slow) then 
        if (debug) write (iw, *) 'DOING FULL SCF DERIVATIVES' 
        call deritr () 
        icalcn = numcal 
!
! THE ARRAY ERRFN HOLDS THE EXACT DERIVATIVES MINUS THE APPROXIMATE
! DERIVATIVES
        errfn(:nvar) = errfn(:nvar) - gradnt(:nvar) 
      endif 
!
!  AT THIS POINT, THE INTERNAL DERIVATIVES ARE KNOWN, AND ARE IN GRADNT.
!
      cosine = dot(gradnt,gold,nvar)/sqrt(dot(gradnt,gradnt,nvar)*dot(gold,gold,nvar)&
         + 1.D-20) 
      if (slow) gradnt(:nvar) = gradnt(:nvar) + errfn(:nvar) 
      if (aic) then 
        if (aifrst) then 
          aifrst = .FALSE. 
          aicorr(:nvar) = (-aidref(:nvar)) - gradnt(:nvar) 
        endif 
        gradnt(:nvar) = gradnt(:nvar) + aicorr(:nvar) 
      endif 
      if (debug) then 
        write (iw, '('' GRADIENTS'')') 
        write (iw, '(8F10.3)') (gradnt(i),i=1,nvar) 
        write (iw, '('' XPARAM'')') 
        write (iw, '(8F10.3)') (xparam(i),i=1,nvar) 
        if (slow) then 
          write (iw, '('' ERROR FUNCTION'')') 
          write (iw, '(8F10.3)') (errfn(i),i=1,nvar) 
        endif 
      endif 
      if (debug) write (iw, '('' COSINE OF SEARCH DIRECTION ='',F30.6)') cosine 
      return  
      end subroutine deriv 
