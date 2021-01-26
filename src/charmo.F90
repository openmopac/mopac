      double precision function charmo (vects, ntype, jorb, ioper, r, nvecs, first) 
      use symmetry_C, only : elem, jelem
      use molkst_C, only : norbs, numat
      implicit none
      integer , intent(in) :: jorb 
      integer  :: ioper 
      integer , intent(in) :: nvecs 
      logical, intent(inout)  :: first 
      integer , intent(in) :: ntype(*) 
      double precision , intent(in) :: vects(nvecs,nvecs) 
      double precision  :: r(3,3) 
      integer , dimension(2,3) :: ip 
      integer , dimension(2,5) :: id 
      integer , dimension(2,9) :: loc 
      integer :: i, iatom, jatom, ibase, kj, icheck, jcheck, ii, jj 
      double precision, dimension(norbs) :: vect1, vect2 
      double precision, dimension(5) :: h = 0.d0
      double precision, dimension(3) :: p 
      double precision, dimension(5) :: d 
      double precision :: sum 
      double precision, external :: ddot
!-----------------------------------------------
!
!  Trivial case:  Operation is 'E', the identity.
!
      charmo = 1.D0 
      if (ioper == 1) return  
!
!   Non-trivial case
!
      vect1(:nvecs) = 0.D0 
      vect2(:nvecs) = 0.D0 
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
        if (ibase > 0) then
!
!   's'-type basis function
!
          icheck = loc(1,1) 
          jcheck = loc(2,1) 
          vect1(icheck) = vects(icheck,jorb) 
          vect2(jcheck) = vects(icheck,jorb) 
        end if
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
            p(i-1) = vects(icheck,jorb) 
            ip(1,i-1) = loc(1,i) 
            ip(2,i-1) = loc(2,i) 
          else 
            d(i-4) = vects(icheck,jorb) 
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
            vect1(ii) = h(i) 
            vect2(jj) = p(i) 
          end do 
        endif 
        if (ibase /= 9) cycle  
!
!   'd' transform
!
        call dtrans (d, ioper, first, r) 
        do i = 1, 5 
          if (id(1,i) < 1) return  
          ii = id(1,i) 
          jj = id(2,i) 
          vect1(ii) = h(i) 
          vect2(jj) = d(i) 
        end do 
      end do 
!
!   Evaluate the term <psi(i)|R(theta)|psi(i)>
!
      sum = ddot(nvecs,vect1(:nvecs),1,vect2(:nvecs),1)
      charmo = sum 
      return  
      end function charmo 
