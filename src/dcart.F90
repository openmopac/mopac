      subroutine dcart(coord, dxyz) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!
      use common_arrays_C, only : nfirst, nlast, nat, p, pa, pb, tvec, &
      nbonds, ibonds, geoa, geo
!
      USE molkst_C, only : numat, numcal, keywrd, id, l1u, l2u, l3u, l123, use_ref_geo, &
      cutofp, method_pm6, method_PM7, mozyme, density, N_3_present, Si_O_H_present
!
      use MOZYME_C, only : iorbs, part_dxyz, mode, jopt
!
      use parameters_C, only : tore
!
      USE molmec_C, only : nnhco, nhco, htype
!
      USE funcon_C, only : fpc_9, a0, ev
!
      USE chanel_C, only : iw 
!
      USE elemts_C, only : elemnt 
!
      USE cosmo_C, only : useps 
!
!***********************************************************************
!DECK MOPAC
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dhc_I 
      use dihed_I 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision  :: coord(3,numat) 
      double precision  :: dxyz(3,numat*l123)  
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icalcn, icuc, numtot, i, j, ii, iii, im1, if, il, jj, &
        jjj, jf, jl, kkkk, ik, jk, kl, l, ij, k, loop, i2, j2, Si, O, H
      double precision, dimension(3,numat) :: work2 
      double precision, dimension(171) :: pdi, padi, pbdi 
      double precision, dimension(3,2) :: cdi 
      double precision, dimension(numat) :: q
      integer :: ndi(2), ione
      double precision :: chnge, chnge2, const, aa, ee, deriv, del, angle, refh, &
        heat, sum, sumx, sumy, sumz, half, rij, der, dstat(3) = 0.d0
      double precision, external :: derp, nsp2_atom_correction, Si_O_H_bond_correction
      logical :: debug, force, large, refeps, point
      integer, external :: ijbo
      save debug, force, large, chnge, chnge2, icalcn, ione, const
!***********************************************************************
!
!    DCART CALCULATES THE DERIVATIVES OF THE ENERGY WITH RESPECT TO THE
!          CARTESIAN COORDINATES. THIS IS DONE BY FINITE DIFFERENCES.
!
!    THE MAIN ARRAYS IN DCART ARE:
!        DXYZ   ON EXIT CONTAINS THE CARTESIAN DERIVATIVES.
!
!*********************************************************************** 
      data icalcn/ 0/  
      data chnge/ 1.D-4/  
      chnge2 = chnge*0.5D0 
      aa = 0.d0
!
! CHNGE IS A MACHINE-PRECISION DEPENDENT CONSTANT
! CHNGE2=CHNGE/2
!
      if (icalcn /= numcal) then 
        icalcn = numcal 
        const = fpc_9
        if (id == 0) then
          ione = 1
        else
          ione = 0
        end if
        large = index(keywrd,'LARGE') /= 0 
        debug = index(keywrd,'DERIV') + index(keywrd,'DCART') /= 0 
        force = index(keywrd,'PREC') + index(keywrd,'FORCE') /= 0 
      endif 
      icuc = (l123 + 1)/2 
      numtot = numat*l123       
      refeps = useps
      useps = .false.
      if (mozyme) then
   !
   !   MODE = 1:   ADD NEW DERIVATIVES ON TO OLD DERIVATIVES
   !          0:   CALCULATE ALL THE DERIVATIVES 'DE NOVO'
   !         -1:   CALCULATE 'OLD' DERIVATIVES, GIVEN ALL THE DERIVATIVES
        if (mode == 0) then
          dxyz(1:3, 1:numtot) = 0.d0
        else if (mode == 1) then
          dxyz(1:3, 1:numtot) = part_dxyz(1:3, 1:numtot)
        else if (mode ==-1) then
          dxyz(1:3, 1:numtot) = -dxyz(1:3, 1:numtot)
        end if
        call chrge_for_MOZYME(p, q)
      else
        dxyz(:,:numtot) = 0.D0 
        call chrge(p, q)
        mode = 0
      end if     
      q(:numat) = tore(nat(:numat)) - q(:numat)
      i2 = 1
      do ii = 1, numat 
        if (mozyme) then
          if (mode == 0 .or. jopt(i2) == ii) then
            i2 = i2 + 1
          else
            cycle
          end if
        end if
        iii = l123*(ii - 1) 
        im1 = ii - ione 
        if = nfirst(ii) 
        il = nlast(ii) 
        ndi(2) = nat(ii) 
        cdi(:,2) = coord(:,ii) 
        j2 = 1
        do jj = 1, im1 
          if (mozyme) then
            if (mode == 0 .or. jopt(j2) == jj) then
              j2 = j2 + 1
            else
              cycle
            end if
          end if
          if (ii == jj) then
            half = 0.5d0
          else
            half = 1.d0
          end if
          jjj = l123*(jj - 1) 
!  FORM DIATOMIC MATRICES
          jf = nfirst(jj) 
          jl = nlast(jj) 
!   GET FIRST ATOM
          ndi(1) = nat(jj) 
          if (jj == 1 .and. ii == 5) then
                   deriv = deriv
                   end if
          if (mozyme) then
            if (ijbo(ii, jj) >= 0) then
                  ! GET FIRST ATOM
              k = ijbo (jj, jj)
              ij = 0
              do i = 1, iorbs(jj)
                do j = 1, i
                  ij = ij + 1
                  k = k + 1
                  padi(ij) = p(k) * 0.5d0
                  pdi(ij) = p(k)
                end do
              end do
                  ! GET SECOND ATOM FIRST ATOM INTERSECTION
              if (ii == jj) then
                 ij = iorbs(jj)
                 do i = 1, iorbs(ii)
                  ij = ij + 1
                  l = (ij*(ij-1)) / 2
                  do j = 1, iorbs(jj)
                    l = l + 1
                    padi(l) = 0.d0
                    pdi(l) = 0.d0
                  end do
                end do
              else
                ij = iorbs(jj)
                k = ijbo (ii, jj)
                do i = 1, iorbs(ii)
                  ij = ij + 1
                  l = (ij*(ij-1)) / 2
                  do j = 1, iorbs(jj)
                    l = l + 1
                    k = k + 1
                    padi(l) = p(k) * 0.5d0
                    pdi(l) = p(k)
                  end do
                end do
              end if
                  ! GET SECOND ATOM
              k = ijbo (ii, ii)
              ij = iorbs(jj)
              do i = 1, iorbs(ii)
                ij = ij + 1
                l = (ij*(ij-1)) / 2 + iorbs(jj)
                do j = 1, i
                  k = k + 1
                  l = l + 1
                  padi(l) = p(k) * 0.5d0
                  pdi(l) = p(k)
                end do
              end do
              pbdi = padi
              point = .false.
            else
              point = .true.
            end if
          else
          point = .false.
          ij = 0 
            do i = jf, jl 
              k = (i*(i - 1))/2 + jf - 1 
              if (i - jf + 1 > 0) then 
                padi(ij+1:i-jf+1+ij) = pa(k+1:i-jf+1+k) 
                pbdi(ij+1:i-jf+1+ij) = pb(k+1:i-jf+1+k) 
                pdi(ij+1:i-jf+1+ij) = p(k+1:i-jf+1+k) 
                ij = i - jf + 1 + ij 
              endif 
            end do 
! GET SECOND ATOM FIRST ATOM INTERSECTION
            do i = if, il 
              l = (i*(i - 1))/2 
              k = l + jf - 1 
              if (jl - jf + 1 > 0) then 
                padi(ij+1:jl-jf+1+ij) = pa(k+1:jl-jf+1+k) 
                pbdi(ij+1:jl-jf+1+ij) = pb(k+1:jl-jf+1+k) 
                pdi(ij+1:jl-jf+1+ij) = p(k+1:jl-jf+1+k) 
                ij = jl - jf + 1 + ij 
              endif 
              k = l + if - 1 
              if (i - if + 1 > 0) then 
                padi(ij+1:i-if+1+ij) = pa(k+1:i-if+1+k) 
                pbdi(ij+1:i-if+1+ij) = pb(k+1:i-if+1+k) 
                pdi(ij+1:i-if+1+ij) = p(k+1:i-if+1+k) 
                ij = i - if + 1 + ij 
              endif 
            end do 
          end if
          kkkk = 0 
          do ik = -l1u, l1u 
            do jk = -l2u, l2u 
              do kl = -l3u, l3u 
                kkkk = kkkk + 1 
                cdi(:,1) = coord(:,jj) + tvec(:,1)*ik + tvec(:,2)*jk + tvec(:,3)*kl 
                if (id /= 0) then
                  rij = (cdi(1, 1)-cdi(1, 2)) ** 2 &
                       & + (cdi(2, 1)-cdi(2, 2)) ** 2 &
                       & + (cdi(3, 1)-cdi(3, 2)) ** 2
                  if (rij > (2.d0/3.d0*cutofp)**2) then
                        !
                        !   Use point-charge approximation
                        !
                    rij = Sqrt (rij)
                    der = derp (rij)
                    ee = q(ii) * q(jj) * fpc_9 * ev * a0 * der
                    do k = 1, 3
                      deriv = half * ee * (cdi(k, 1)-cdi(k, 2)) / rij
                      dxyz(k, iii+icuc) = dxyz(k, iii+icuc) - deriv
                      dxyz(k, jjj+kkkk) = dxyz(k, jjj+kkkk) + deriv
                    end do
                    cycle
                  end if
                end if
                if (point) then
                  do l = 1, 3
                    cdi(l, 1) = coord(l, jj) + tvec(l, 1) * ik &
                         & + tvec(l, 2) * jk + tvec(l, 3) * kl
                  end do
                  call delsta (nat, iorbs, p, cdi, dstat, ii, jj)  
                  dxyz(1:3, iii+icuc) = dxyz(1:3, iii+icuc) - dstat(1:3)
                  dxyz(1:3, jjj+kkkk) = dxyz(1:3, jjj+kkkk) + dstat(1:3) 
                else
                  if (.not.force) then 
                   cdi(1,1) = cdi(1,1) + chnge2 
                   cdi(2,1) = cdi(2,1) + chnge2 
                   cdi(3,1) = cdi(3,1) + chnge2 
                   call dhc (pdi, padi, pbdi, cdi, ndi, jf, jl, if, il, aa, 1) 
                  endif 
                  do k = 1, 3 
                    if (force) then 
                      cdi(k,2) = cdi(k,2) - chnge2 
                      call dhc (pdi, padi, pbdi, cdi, ndi, jf, jl, if, il, aa, 1) 
                    endif 
                    cdi(k,2) = cdi(k,2) + chnge 
                    call dhc (pdi, padi, pbdi, cdi, ndi, jf, jl, if, il, ee, 2) 
                    cdi(k,2) = cdi(k,2) - chnge2 
                    if (.not.force) cdi(k,2) = cdi(k,2) - chnge2 
                    deriv = half*(aa - ee)*const/chnge 
                    dxyz(k,iii+icuc) = dxyz(k,iii+icuc) - deriv 
                    dxyz(k,jjj+kkkk) = dxyz(k,jjj+kkkk) + deriv 
                  end do 
                end if
              end do 
            end do 
          end do 
        end do 
      end do 
      if (nnhco /= 0) then 
!
!   NOW ADD IN MOLECULAR-MECHANICS CORRECTION TO THE H-N-C=O TORSION
!
        del = 1.D-8 
        do i = 1, nnhco 
          do j = 1, 4 
            do k = 1, 3 
              coord(k,nhco(j,i)) = coord(k,nhco(j,i)) - del 
              call dihed (coord, nhco(1,i), nhco(2,i), nhco(3,i), nhco(4,i), &
                angle) 
              refh = htype*sin(angle)**2 
              coord(k,nhco(j,i)) = coord(k,nhco(j,i)) + del*2.D0 
              call dihed (coord, nhco(1,i), nhco(2,i), nhco(3,i), nhco(4,i), &
                angle) 
              coord(k,nhco(j,i)) = coord(k,nhco(j,i)) - del 
              heat = htype*sin(angle)**2 
              sum = (refh - heat)/(2.D0*del) 
              dxyz(k,nhco(j,i)) = dxyz(k,nhco(j,i)) - sum 
            end do 
          end do 
        end do 
      endif 
      if (method_pm6 .and. N_3_present) then
!
!   Add in the nitrogen sp2 correction
!
        del = 1.d-8 
        do i = 1, numat
          if (nat(i) == 7 .and. nbonds(i) == 4) then
            jj = 0
            if (nat(ibonds(2,i)) == 1) jj = 1
            if (nat(ibonds(3,i)) == 1) jj = jj + 1
            if (nat(ibonds(4,i)) == 1) jj = jj + 1
            if ( jj < 2) then
              do j = 1,4
                do k = 1, 3
                  coord(k,ibonds(j,i)) = coord(k,ibonds(j,i)) - del 
                  sum = nsp2_atom_correction(coord, i, ibonds(2,i), ibonds(3,i), ibonds(4,i))
                  coord(k,ibonds(j,i)) = coord(k,ibonds(j,i)) + 2.d0*del 
                  sum = (sum - nsp2_atom_correction(coord, i ,ibonds(2,i), ibonds(3,i), ibonds(4,i)))/(2.d0*del)
                  coord(k,ibonds(j,i)) = coord(k,ibonds(j,i)) - del 
                  dxyz(k,ibonds(j,i)) = dxyz(k,ibonds(j,i)) - sum 
                end do
              end do
            end if
          end if
        end do
      end if
      if (method_pm7 .and. Si_O_H_present) then
!
!   Add in the Si-O-H correction
!
      del = 1.d-8 
      do i = 1, numat
          if (nat(i) == 8) then
            O = i
            Si = 0
            H = 0
            do j = 1, nbonds(i)
              k = ibonds(j,i)
              if (nat(k) == 14) Si = k
              if (nat(k) == 1) H = k              
            end do
            if (Si /= 0 .and. H /= 0) then
              do k = 1, 3
!
! Si
!
                coord(k,Si) = coord(k,Si) - del 
                sum = Si_O_H_bond_correction(coord, Si, O, H)
                coord(k,Si) = coord(k,Si) + 2.d0*del 
                sum = (sum - Si_O_H_bond_correction(coord, Si, O, H))/(2.d0*del)
                coord(k,Si) = coord(k,Si) - del 
                dxyz(k,Si) = dxyz(k,Si) - sum 
!
! O
!
                coord(k,O) = coord(k,O) - del 
                sum = Si_O_H_bond_correction(coord, Si, O, H)
                coord(k,O) = coord(k,O) + 2.d0*del 
                sum = (sum - Si_O_H_bond_correction(coord, Si, O, H))/(2.d0*del)
                coord(k,O) = coord(k,O) - del 
                dxyz(k,O) = dxyz(k,O) - sum 
!
! H
!
                coord(k,H) = coord(k,H) - del 
                sum = Si_O_H_bond_correction(coord, Si, O, H)
                coord(k,H) = coord(k,H) + 2.d0*del 
                sum = (sum - Si_O_H_bond_correction(coord, Si, O, H))/(2.d0*del)
                coord(k,H) = coord(k,H) - del 
                dxyz(k,H) = dxyz(k,H) - sum 
              end do
            end if        
          end if
        end do
      end if
      if (mode ==-1) then
        part_dxyz(1:3, 1:numat*l123) = -dxyz(1:3, 1:numat*l123)
      end if
      useps = refeps
      if (useps) call diegrd (dxyz) 
      if (use_ref_geo) then
        do i = 1, numat
          do j = 1,3
            dxyz(j,i) = dxyz(j,i) + (geo(j,i) - geoa(j,i))*density*2.d0
          end do
        end do
      end if
      if (.not.debug) return 
      if (l123 == 1) then 
        write (iw, '(I6,4x,a2,F13.6,2F13.6)') (i,elemnt(nat(i)),(dxyz(j,i),j=1,3),i=1,numtot) 
      else if (large) then 
      call print_dxyz("(Does NOT include post-SCF corrections)")
      endif 
      if (id == 0) return  
      write (iw, &
      '(2/10X,"CARTESIAN COORDINATE DERIVATIVES FOR THE CENTRAL UNIT CELL",/,"  NO. AT.     X            Y            Z",/)') 
      if (l123 == 1) then 
        write (iw, '(I6,A2,3F13.6)') (i,elemnt(nat(i)),(dxyz(j,i),j=1,3),i=1,numtot) 
      else if (large) then 
        loop = 0 
        do i = 1, numat 
          sumx = 0.D0 
          sumy = 0.D0 
          sumz = 0.D0 
          do ik = -l1u, l1u 
            do jk = -l2u, l2u 
              do kl = -l3u, l3u 
                loop = loop + 1 
                sumx = sumx + dxyz(1,loop) 
                sumy = sumy + dxyz(2,loop) 
                sumz = sumz + dxyz(3,loop) 
              end do 
            end do 
          end do 
          work2(1,i) = sumx 
          work2(2,i) = sumy 
          work2(3,i) = sumz 
        end do 
        write (iw, '(I6,A2,F13.6,2F13.6)') (i,elemnt(nat(i)), &
        (work2(j,i),j=1,3), i = 1, numat) 
        if (id == 3) call xyzcry (tvec, numat, work2, iw) 
      else 
        write (iw, '(I6,A2,F13.6,2F13.6)') (i,elemnt(nat((i-1)/l123+1)), &
        (dxyz(j,i) + dxyz(j,i+1) + dxyz(j,i+2),j=1,3),i = 1, numtot - 2, 3) 
      endif
      return  
      end subroutine dcart 
!
      double precision function derp (r)
      use molkst_C, only: numcal, clower, cutofp, cupper
      implicit none
      double precision, intent (in) :: r
      integer, save :: icalcn = 0
      double precision, save :: bound1, bound2, c, cr, cr2, range
!
      if (icalcn /= numcal) then
      !
      ! Set constants for truncation function.
      !
      ! CLOWER = lower bound of truncation function, as a function of CUTOFP
        bound1 = clower / cutofp
      !
      ! CUPPER = upper bound of truncation function, as a function of CUTOFP
        bound2 = cupper / cutofp
        range = bound2 - bound1
      !
      ! Truncation function = C+CR*R+CR2*R**2
      !
        c = -0.5d0 * bound1 ** 2 * cutofp / range
        cr = 1.d0 + bound1 / range
        cr2 = -1 / (cutofp*2*range)
      !
      !   Above CUPPER function = constant
      !
        icalcn = numcal
      end if
   !
   !  Need a smooth function in the region of CUTOFP
   !
   !  Function has form:
   !    Up to CLOWER  R=R
   !    At CLOWER,  slope = 1.0
   !    Between CLOWER and CUPPER = Monotomic increase,
   !    with rate of increase dropping from 1.0 to 0.0
   !    At CUPPER, slope = 0
   !    Above CUTOFP   R=CUTOFP
      if (r > clower) then
        if (r > cupper) then
          derp = 0.d0
        else
          derp = -(cr+2*cr2*r) / (c+cr*r+cr2*r**2) ** 2
        end if
      else
        derp = -1 / r ** 2
      end if
  end function derp
  subroutine print_dxyz(header)
!
! Generic print of the array dxyz.  Use this wherever array dxyz needs to be printed
!
    use common_arrays_C, only : dxyz, nat
!
    USE molkst_C, only : numat, keywrd, l123, l11, l21, l31
    USE chanel_C, only : iw 
    USE elemts_C, only : elemnt 
!
    implicit none
    character :: header*(*)
    integer :: i1, j1, k1, i, k, j, l
    double precision :: sum
    logical :: large
    large = index(keywrd,'LARGE') /= 0 
    write (iw, '(2/10X,''CARTESIAN COORDINATE DERIVATIVES'')') 
    write (iw, '(7X,a)')trim(header)
    if (l123 == 1) then 
      write (iw, '(/1X, a, /)')" NUMBER ATOM           X                Y                Z              Total" 
      write (iw, '(I6,4x,a2,4F17.6)') (i, elemnt(nat(i)),(dxyz(i*3 - 3 + j), j = 1, 3), &
        sqrt(dxyz(i*3)**2 + dxyz(i*3 - 1)**2 + dxyz(i*3 - 2)**2), i = 1, numat)
    else if (large) then 
      write (iw, '(/1X, a, /)')"       CELL           ATOM            X                Y                Z            Total" 
      k = 0 
      l = 0
      do i1 = -l11, l11 
        do j1 = -l21, l21 
          do k1 = -l31, l31 
            l = l + 1
            do i = 1, numat              
              k = k + 1 
              sum = dxyz(k*3 - 2)**2 + dxyz(k*3 - 1)**2 + dxyz(k*3)**2
              if (sum > 0.1d0) &
              write (iw, '(I6, 2i4, i8, i4, 1x,a2,F13.6,3F17.6)')  i1, j1, k1, k, i, &
                elemnt(nat(i)), (dxyz(k*3 - 3 + j), j = 1, 3), sqrt(sum)
            end do
          end do
        end do
      end do
    endif 
    end subroutine print_dxyz

