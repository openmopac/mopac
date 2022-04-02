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

  subroutine new_esp()
    use common_arrays_C, only : coord, c, p, nat, &
    nfirst, nlast
    use parameters_C, only : tore
    use molkst_C, only : numat, norbs, fract, nclose, nopen, &
      keywrd, jobnam
    use funcon_C, only : a0, pi
    use chanel_C, only : iesp
    use esp_C, only : ixn, iyn, izn, jxn, jyn, jzn
    USE overlaps_C, only : ccc, zzz
!

    implicit none
    logical dubl, sames, cube, espgrid
    integer :: i, j, ij, nx, ny, nz, iiz, iiy, iix, &
    mini, minj, maxi, maxj, lit, ljt, nshells, k, l, ig, jg, &
    isize, isizes, jsize, jsizes, jgmax, idx, nroots, n, mm, in, jn, &
    pij, kl, lmax, nxyz, nFit, m, el, nUse, nnx, nny, nnz

    real :: rtemp, dx, dy, dz, ep(81)

    double precision :: x, y, z, delx, dely, delz,  &
    iatomx, iatomy, iatomz, jatomx, jatomy, jatomz, px, py, pz, &
    icoeff, jcoeff, dtemp1, dtemp2, dtemp3, rr, iexp, jexp, &
    gamma, gammainverse, kfac, contractionDensity(100), tosp, yz, xyz, &
    rx, u(4), w(4), uu, ww, ttinverse, t, x0, y0, z0, &
    tt, xint, yint, zint, xin(64), yin(64), zin(64), &
     EUpperRange, &
    inorm, jnorm, &
    xmax = -1.d8, &
    xmin =  1.d8, &
    ymax = -1.d8, &
    ymin =  1.d8, &
    zmax = -1.d8, &
    zmin =  1.d8

    integer, dimension(:), allocatable :: jpvt, itypes, istart, istop, &
      & icoord
    real, dimension (:), allocatable :: aMatrix, qraux, bVector, aVector, &
                                      & rsd, qy, qty, xb, esp_array
    double precision, allocatable, dimension(:,:) :: norm


    integer ::  ijx(100), ijy(100), ijz(100), &
    itype, jtype
    real, dimension (:), allocatable :: planexy,  Scr2
    double precision, dimension (:,:), allocatable :: s, vecs
    double precision, external :: reada
    mini = 0
    minj = 0
    allocate (s(norbs, norbs), vecs(norbs, norbs), istart(norbs), istop(norbs))
    nFit = 0
    call get_minus_point_five_overlap(s)
    call mult (c, s, vecs, norbs)

! For GPU MOPAC
! GBR_new_addition

! ORG      call densit (vecs, norbs, norbs, nclose, 2.d0, nopen, fract, p, 2)
    ij = (norbs*(norbs+1)/2)
        call density_for_GPU (vecs,fract,nclose,nopen, &
                             & 2.d0,ij,norbs,2,p,5)
! this procedure was replaced with MKL/DSCAL calling
!    ij = 0
!    do i = 1, norbs
!      do j = 1, i - 1
!        ij = ij + 1
!        p(ij) = p(ij)*2.d0
!      end do
!      ij = ij + 1
!    end do

    call dscal(ij, 2.d0, p, 1)
    j = 0
    do i = 1,norbs
       j = i*(i+1)/2
       p(j) = p(j)/2.d0
    end do
!!!

    tosp = 2.d0/sqrt(pi)
!
!  Convert to AU and work out the upper and lower bounds of the box
!
    coord = coord/a0

    do i = 1, numat
      if (coord(1,i) > xmax) xmax = coord(1,i)
      if (coord(2,i) > ymax) ymax = coord(2,i)
      if (coord(3,i) > zmax) zmax = coord(3,i)
      if (coord(1,i) < xmin) xmin = coord(1,i)
      if (coord(2,i) < ymin) ymin = coord(2,i)
      if (coord(3,i) < zmin) zmin = coord(3,i)
    end do
    xmax = xmax + 4.d0
    ymax = ymax + 4.d0
    zmax = zmax + 4.d0
    xmin = xmin - 4.d0
    ymin = ymin - 4.d0
    zmin = zmin - 4.d0
  !  nnx = max(15, min(25, nint((xmax - xmin))))
  !  nnx = max(15, min(25, nint((ymax - ymin))), nnx)
  !  nnx = max(15, min(25, nint((zmax - zmin))), nnx)
  !  nny = nnx ! max(15, nint((ymax - ymin)*1.5d0))
  !  nnz = nnx ! max(15, nint((zmax - zmin)*1.5d0))


    nnx = max(15, min(25, nint((xmax - xmin))))
    nny = max(15, min(25, nint((ymax - ymin))))
    nnz = max(15, min(25, nint((zmax - zmin))))
    allocate (planexy(nnx*nny*nnz), scr2(nnx*nny*nnz))
    delx = (xmax - xmin)/(nnx - 1)
    dely = (ymax - ymin)/(nny - 1)
    delz = (zmax - zmin)/(nnz - 1)
!
!  Set up arrays for handling the STO6G orbitals
!
    call setupg
    allocate (itypes(norbs), icoord(norbs), norm(norbs,6))
    nshells = 0
    do i = 1, numat
      nshells = nshells + 1
      itypes(nshells) = 0   !  "s"-type
      icoord(nshells) = i
      istart(nshells) = nfirst(i)
      istop(nshells) = nfirst(i)
      norm(nshells,:) = (2.d0*zzz(nshells,:)/pi)**0.75
      if (nlast(i) == nfirst(i)) cycle
      nshells = nshells + 1
      itypes(nshells) = 1   !  "p"-type
      icoord(nshells) = i
      istart(nshells) = nfirst(i) + 1
      istop(nshells) = nfirst(i) + 3
      norm(nshells,:) = (128.d0/pi**3)**0.25*zzz(nshells,:)**1.25d0
    end do
!
!  Contribution to the ESP arising from the nuclii
!
    planexy = 0.e0
    pij = 0
    do iiz = 1, nnz
      z = zmin + delz*(iiz - 1)
      do iiy = 1, nny
        y = ymin + dely*(iiy - 1)
        do iix = 1, nnx
          x = xmin + delx*(iix - 1)
          pij = pij + 1
          do i = 1, numat
            dx = sngl(x - coord(1,i))
            dy = sngl(y - coord(2,i))
            dz = sngl(z - coord(3,i))
            planexy(pij) = planexy(pij) + sngl(tore(nat(i)))/sqrt(max(dx**2 + dy**2 + dz**2, 1.e-12))
          end do
        end do
      end do
    end do
!
!  Contribution to the ESP arising from the electrons
!
    do i = 1, nshells
      itype = itypes(i)
      iatomx = coord(1, icoord(i))
      iatomy = coord(2, icoord(i))
      iatomz = coord(3, icoord(i))
      select case (itype)
      case (0) !  "s"-type
        isize  = 1
        isizes = 1
        mini   = 1
        maxi   = 1
        lit    = 1
      case (1) !  "p"-type
        isize  = 3
        isizes = 3
        mini   = 2
        maxi   = 4
        lit    = 2
      end select

      do j = 1, i
        sames = (i == j)
        jtype = itypes(j)
        jatomx = coord(1, icoord(j))
        jatomy = coord(2, icoord(j))
        jatomz = coord(3, icoord(j))
        select case (jtype)
        case (0) !  "s"-type
          jsize  = 1
          jsizes = 1
          minj   = 1
          maxj   = 1
          ljt    = 1
        case (1) !  "p"-type
          jsize  = 3
          jsizes = 3
          minj   = 2
          maxj   = 4
          ljt    = 2
        end select
        dtemp1 = iatomx - jatomx
        dtemp2 = iatomy - jatomy
        dtemp3 = iatomz - jatomz
        rr = dtemp1**2 + dtemp2**2 + dtemp3**2
        if (sames) then
          ij = 0
          do k = mini, maxi
            do l = minj, k
              ij = ij + 1
              ijx(ij) = ixn(k) + jxn(l) + 1
              ijy(ij) = iyn(k) + jyn(l) + 1
              ijz(ij) = izn(k) + jzn(l) + 1
            end do
          end do
        else
          ij = 0
          do k = mini, maxi
            do l = minj, maxj
              ij = ij + 1
              ijx(ij) = ixn(k) + jxn(l) + 1
              ijy(ij) = iyn(k) + jyn(l) + 1
              ijz(ij) = izn(k) + jzn(l) + 1
            end do
          end do
        end if
        do ig = 1,6
          iexp = zzz(i,ig)
          icoeff = ccc(i,ig)
          inorm = norm(i,ig)
          if (sames) then
            jgmax = ig
          else
            jgmax = 6
          end if
          do jg = 1, jgmax
            dubl = (sames .and. ig /= jg)
            jexp = zzz(j,jg)
            jnorm = norm(j,jg)
            jcoeff = ccc(j,jg)
            gamma = iexp + jexp
            gammainverse = 1.d0/gamma
            dtemp1 = iexp*jexp*rr*gammainverse
            if (dtemp1 < 46.d0) then
              kfac = exp(-dtemp1)
              px = (iexp*iatomx + jexp*jatomx)*gammainverse
              py = (iexp*iatomy + jexp*jatomy)*gammainverse
              pz = (iexp*iatomz + jexp*jatomz)*gammainverse
              dtemp1 = icoeff*jcoeff*inorm*jnorm*kfac
              idx = 0
              if (sames) then
                if (dubl) then
                  do k = 1, isize
                    do l = 1, k
                      idx = idx + 1
                      contractionDensity(idx) = 2.d0*dtemp1
                    end do
                  end do
                else
                   do k = 1, isize
                    do l = 1, k
                      idx = idx + 1
                      contractionDensity(idx) = dtemp1
                    end do
                  end do
                end if
              else
                do k = 1, isize
                  do l = 1, jsize
                      idx = idx + 1
                      contractionDensity(idx) = dtemp1
                  end do
                end do
              end if
              contractionDensity(:ij) = contractionDensity(:ij)*tosp*gammainverse
              pij = 1
              do iiz = 1, nnz
                z = zmin + delz*(iiz - 1)
                dz = sngl(pz -  z)
                do iiy = 1, nny
                  y = ymin + dely*(iiy - 1)
                  dy = sngl(py - y)
                  yz = dy**2 + dz**2
                  do iix = 1, nnx
                    x = xmin + delx*(iix - 1)
                    isize = isizes
                    jsize = jsizes
                    dx = sngl(px - x)
                    xyz = max(dx**2 + yz, 1.d-2)
                    rx = gamma*xyz
                    nroots = (itype + jtype)/2 + 1
!
!  Solve the Rys polynomials.  See: J. Rys, M. Dupuis,
!  and H. F. King, J. Comput. Chem. 4,154 (1983).
!
                    call rys(rx, nroots, u, w)
                    mm = 0
                    do n = 1, nroots
                      uu = gamma*u(n)
                      ww = -w(n)
                      tt = gamma + uu
                      ttinverse = 1.d0/tt
                      t = sqrt(tt)
                      t = 1.d0/t
                      x0 = (px*gamma + uu*x)*ttinverse
                      y0 = (py*gamma + uu*y)*ttinverse
                      z0 = (pz*gamma + uu*z)*ttinverse
                      in = mm
                      do k = 1, lit
                        do l = 1, ljt
                          jn = in + l
!
!  Solve for the one-dimensional integrals using Gauss-Hermite Quadrature
!
                          call vint(xint, yint, zint, k,l, x0, y0, z0, &
                          iatomx, iatomy, iatomz, jatomx, jatomy, jatomz, t)
                          xin(jn) = xint
                          yin(jn) = yint
                          zin(jn) = zint * ww
                        end do
                        in = in + 4
                      end do
                      mm = mm + 16
                    end do
                    do k = 1, ij
                      nx = ijx(k)
                      ny = ijy(k)
                      nz = ijz(k)
                      dtemp1 = 0.d0
                      mm = 0
                      do l = 1,nroots
                        dtemp1 = dtemp1 + xin(nx + mm)*yin(ny + mm)*zin(nz + mm)
                        mm = mm + 16
                      end do
                      ep(k) = sngl(contractiondensity(k)*dtemp1)
                    end do
                    kl = 0
                    do k = istart(i), istop(i)
                      if (sames) then
                        lmax = k
                      else
                        lmax = istop(j)
                      end if
                      do l = istart(j), lmax
                        kl = kl + 1
                        planexy(pij) = planexy(pij) + ep(kl)*sngl(p((k*(k - 1))/2 + l))
                        if (pij == 231) then
                          pij = pij
                        end if
                      end do          !  "l" shell loop
                    end do            !  "k" shell loop
                  pij = pij + 1
                  end do              !  "x" loop
                end do                !  "y" loop
              end do                  !  "z" loop
            end if
          end do                      ! Outermost "jg" loop
        end do                        ! Outermost "ig" loop
      end do                          ! Outermost "j" loop
    end do
    nxyz = pij - 1                 ! Outermost "i" loop
    scr2 = planexy
    EUpperRange = 0.9d0
    do
      do i = 1, nxyz
        if (abs(scr2(i)) < EUpperrange .and. abs(scr2(i)) > 1.d-4) nFit = nFit + 1
      end do
      if (nFit > nxyz - 10) exit
      EUpperRange = EUpperRange*2.d0
    end do
    EUpperRange = EUpperRange*0.5d0
    pij = 1
    el = 7*numat
    allocate(amatrix(el*el), qraux(el), bVector(el), aVector(el), &
    rsd(el), qy(el), qty(el), xb(el), jpvt(el))
    jpvt = 0
    aMatrix = 0.e0
    bVector = 0.e0
    i = 0
    l = 0
    do iiz = 1, nnz
      z = zmin + delz*(iiz - 1)
      do iiy = 1, nny
        y = ymin + dely*(iiy - 1)
        do iix = 1, nnx
          x = xmin + delx*(iix - 1)
          l = l + 1
          if (abs(scr2(l)) <  EUpperrange .and. abs(scr2(l)) > 1.d-4) then
            call evec(aVector, x, y, z, coord, numat)
            do m = 1, el
              do n = 1, el
                aMatrix((m - 1)*el + n) =  aMatrix((m - 1)*el + n) + aVector(m)*aVector(n)
              end do
              bVector(m) =  bVector(m) + aVector(m)*Scr2(l)
            end do
            i = i + 1
          end if
          pij = pij + 1
        end do
      end do
    end do
!
!  Use the Householder transformation to compute the QR factorization of the aMatrix of size el.
!
    i = 1
    call sqrdc(aMatrix, el, el, el, qraux, jpvt, aVector, i)
    dtemp2 = Abs(aMatrix(1))*1.d-5
    do k = 1, el
      dtemp1 = Abs(aMatrix((k - 1)*el + k))
      if (dtemp1 < dtemp2) exit
    end do
    nUse = k - 1
!
!  SQRSL applies the output of SQRDC to compute coordinate
!     transformations, projections, and least squares solutions.
!
    i = 111
    call sqrsl(aMatrix, el, el, nUse, qraux, bVector, qy, qty, aVector, rsd, xb, i, j)
    do j = 1, el
      jpvt(j) = -jpvt(j)
      if (j > nUse) aVector(j) = 0.e0
    end do
!
! Untangle the vectors
!
    do j = 1, el
      if (jpvt(j) <= 0) then
        k = -jpvt(j)
        do while (k /= j)
          rtemp = aVector(j)
          aVector(j) = aVector(k)
          aVector(k) = rtemp
          jpvt(k) = -jpvt(k)
          k = jpvt(k)
        end do
      end if
    end do
    xmax = xmax + 4.d0
    xmin = xmin - 4.d0
    ymax = ymax + 4.d0
    ymin = ymin - 4.d0
    zmax = zmax + 4.d0
    zmin = zmin - 4.d0
    i = index(keywrd," ESPGR")
    if (i /= 0) then
      nnx = max( nint(reada(keywrd, i + 8)), 60)
      nny = nnx
      nnz = nnx
    else
      nnx = 60
      nny = 60
      nnz = 60
    end if
    allocate(esp_array(nny*nnz))
    delx = (xmax - xmin)/(nnx - 1)
    dely = (ymax - ymin)/(nny - 1)
    delz = (zmax - zmin)/(nnz - 1)

    cube = (index(keywrd, " CUBE") /= 0)
    espgrid = (index(keywrd, " ESPGRID") /= 0)
    if (cube) then
!
!  Write out ESP data in format for Jmol
!
      open (iesp, file=trim(jobnam)//".grd"  , status="UNKNOWN")
      write (iesp,"(a)") " 4 Density"
      write (iesp,"(a)") " Electron density from Total SCF Density"
      write (iesp, "(i5,3f12.6)")numat, xmin, ymin, zmin
      write (iesp, "(i5, 3f12.6)") nnx, delx, 0.d0, 0.d0
      write (iesp, "(i5, 3f12.6)") nny, 0.d0, dely, 0.d0
      write (iesp, "(i5, 3f12.6)") nnz, 0.d0, 0.d0, delz
      do i = 1, numat
        write(iesp,"(i5,4f12.6)") nat(i), nat(i)*1.d0, coord(1,i), coord(2,i), coord(3,i)
      end do
!
!   Write out electrostatic values
!
      do iix = 1, nnx
      x = xmin + delx*(iix - 1)
      pij = 0
      do iiy = 1, nny
        y = ymin + dely*(iiy - 1)
        do iiz = 1, nnz
          z = zmin + delz*(iiz - 1)
          pij = pij + 1
          call evec(bVector, x, y, z, coord, numat)
          rtemp = 0.e0
          do k = 1, el
            rtemp = rtemp + aVector(k)*bVector(k)
          end do
          esp_array(pij) = rtemp
        end do
      end do
      write(iesp, "(5e13.5)") esp_array(:pij)
    end do
    else if (espgrid) then
      open (iesp, file=trim(jobnam)//'.grd', status="UNKNOWN")
      write (iesp, "(3f12.6)") xmin*a0, ymin*a0, zmin*a0
      write (iesp, "(i5, 3f12.6)") nnx, delx*a0, delx*a0, delx*a0
      write (iesp, "(i5, 3f12.6)") nny, dely*a0, dely*a0, dely*a0
      write (iesp, "(i5, 3f12.6)") nnz, delz*a0, delz*a0, delz*a0
!
!   Write out electrostatic values
!
      do iiz = 1, nnz
      z = zmin + delz*(iiz - 1)
      pij = 0
      do iiy = 1, nny
        y = ymin + dely*(iiy - 1)
        do iix = 1, nnx
          x = xmin + delx*(iix - 1)
          pij = pij + 1
          call evec(bVector, x, y, z, coord, numat)
          rtemp = 0.e0
          do k = 1, el
            rtemp = rtemp + aVector(k)*bVector(k)
          end do
          esp_array(pij) = rtemp
        end do
      end do
      write(iesp, "(e13.5)") esp_array(:pij)
      end do
    end if
  end subroutine new_esp
