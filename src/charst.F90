      real(kind(0.0d0)) function charst (vects, ntype, istate, ioper, r,nvecs, first) 
!-----------------------------------------------
!   M o d u l elem s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double
      USE meci_C, only : nalmat, microa, microb, lab, conf, nmos
      use chanel_C, only : iw
      use symmetry_C, only : elem, jelem
      use molkst_C, only : keywrd, numat, norbs
      use minv_I 
      use matout_I 
      implicit none
      integer , intent(in) :: istate 
      integer  :: ioper 
      integer , intent(in) :: nvecs 
      logical  :: first 
      integer , intent(in) :: ntype(norbs) 
      real(double) , intent(in) :: vects(nvecs,nmos) 
      real(double)  :: r(3,3) 
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
        ii, jj, k, l, ne, ia1, ia2, ia3, ib1, ib2, ib3, nai, ja1, ja2, ja3, jb1&
        , jb2, jb3, nbi 
      real(double), dimension(:,:), allocatable :: vect1, vect2 
      real(double) :: h(5), p(3), d(5) 
      real(double), dimension(:,:), allocatable :: t2, t4
      real(double), dimension(:), allocatable :: work
      real(double) :: sum, det, suma, sumb 
      logical :: posita, positb, debug 
      real(double), external :: ddot
!-----------------------------------------------
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
!    This routine is limited to systems having a maximum of the lesser
!    of three electrons of either spin, or three holes of either spin,
!    in other words all systems involving up to seven M.O.s
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
        endif 
        if (positb) then 
          microb(:nmos,:nstate) = 1 - microb(:nmos,:nstate)   
        endif 
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
      endif 
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
              endif 
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
              endif 
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
            endif 
            if (ibase /= 9) cycle  
!
!   'd' transform
!
            h = d
            call dtrans (d, ioper, first, r) 
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
        endif 
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
                endif 
                microa(i,j) = 1 - microa(i,j) 
              end do 
              iphase(j) = 1 - 2*mod(l,2) 
            end do 
          endif 
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
                endif 
                microb(i,j) = 1 - microb(i,j) 
              end do 
              iphase(j) = 1 - 2*mod(l,2) 
            end do 
          endif 
        endif 
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
        endif 
!
!   The big loop to fill T4
!
        do i = 1, nstate 
          ia1 = iperma(1,i) 
          ia2 = iperma(2,i) 
          ia3 = iperma(3,i) 
          ib1 = ipermb(1,i) 
          ib2 = ipermb(2,i) 
          ib3 = ipermb(3,i) 
          nai = nalmat(i) 
          do j = 1, nstate 
            ja1 = iperma(1,j) 
            ja2 = iperma(2,j) 
            ja3 = iperma(3,j) 
            jb1 = ipermb(1,j) 
            jb2 = ipermb(2,j) 
            jb3 = ipermb(3,j) 
            if (nalmat(j) /= nai) cycle  
!
!    NAI = Number of alpha electrons
!    NBI = Number of beta electrons
            nbi = ne - nai 
            select case (nai + 1)  
!
!  General case: for NAI greater than 3
!
            case default 
              k = 0 
              do ii = 1, nai 
                work(k+1:nai+k) = t2(iperma(ii,i),iperma(:nai,j)) 
                k = nai + k 
              end do 
              call minv (work, nai, suma) 
            case (1)  
              suma = 1.D0 
            case (2)  
              suma = t2(ja1,ia1) 
            case (3)  
              suma = t2(ja1,ia1)*t2(ja2,ia2) - t2(ja2,ia1)*t2(ja1,ia2) 
            case (4)  
              suma = t2(ja1,ia1)*t2(ja2,ia2)*t2(ja3,ia3) - t2(ja1,ia1)*t2(ja2,&
                ia3)*t2(ja3,ia2) - t2(ja1,ia2)*t2(ja2,ia1)*t2(ja3,ia3) + t2(ja1&
                ,ia2)*t2(ja2,ia3)*t2(ja3,ia1) + t2(ja1,ia3)*t2(ja2,ia1)*t2(ja3,&
                ia2) - t2(ja1,ia3)*t2(ja2,ia2)*t2(ja3,ia1) 
            end select 
            select case (nbi + 1)  
            case default 
              k = 0 
              do ii = 1, nbi 
                work(k+1:nbi+k) = t2(ipermb(ii,i),ipermb(:nbi,j)) 
                k = nbi + k 
              end do 
              call minv (work, nbi, sumb) 
            case (1)  
              sumb = 1.D0 
            case (2)  
              sumb = t2(jb1,ib1) 
            case (3)  
              sumb = t2(jb1,ib1)*t2(jb2,ib2) - t2(jb2,ib1)*t2(jb1,ib2) 
            case (4)  
              sumb = t2(jb1,ib1)*t2(jb2,ib2)*t2(jb3,ib3) - t2(jb1,ib1)*t2(jb2,&
                ib3)*t2(jb3,ib2) - t2(jb1,ib2)*t2(jb2,ib1)*t2(jb3,ib3) + t2(jb1&
                ,ib2)*t2(jb2,ib3)*t2(jb3,ib1) + t2(jb1,ib3)*t2(jb2,ib1)*t2(jb3,&
                ib2) - t2(jb1,ib3)*t2(jb2,ib2)*t2(jb3,ib1) 
            end select 
            t4(i,j) = suma*sumb*iphase(i)*iphase(j)*det 
          end do 
        end do 
        if (debug) then 
          write (iw, *) ' State Transform for State', istate, &
            ' under Operation', ioper 
          call matout (t4, t4, nstate, nstate, lab) 
        endif 
      endif 
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
