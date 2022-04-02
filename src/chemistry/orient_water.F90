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

module orient_water_C
  double precision :: near(3,30), O_coord(3)
  integer :: near_nat(30), n_atoms, O_copy, H1_copy, atom_nos(30)
end module orient_water_C
subroutine orient_water(O_of_H2O, H1, H2, fun)
use common_arrays_C,  only : coord, nat
use molkst_C,  only : numat
use orient_water_C, only: near,  O_coord, n_atoms, near_nat, O_copy, H1_copy, atom_nos
implicit none
integer, intent (in) :: O_of_H2O, H1, H2
double precision, intent (out) :: fun
double precision :: x(6)
integer   :: i, n
double precision, external :: distance
  n = 6
  x(1:3) = coord(:,H1)
  x(4:6) = coord(:,H2)
  O_coord = coord(:,O_of_H2O)
  O_copy = O_of_H2O
  H1_copy = H1
!
! Find nearby atoms
!
  n_atoms = 0
  do i = 1, numat
    if (distance(O_of_H2O, i) < 3.5d0) then
      if (i /= O_of_H2O .and. i /= H1 .and. i /= H2) then
        n_atoms = n_atoms + 1
        near(:, n_atoms) = coord(:, i)
        near_nat(n_atoms) = nat(i)
        atom_nos(n_atoms) = i
      end if
    end if
  end do
  if (n_atoms > 0) then
    call cobyla (n, x, fun)
    coord(:,H1) = x(1:3)
    coord(:,H2) = x(4:6)
  end if
  return
end subroutine orient_water

subroutine calcfc (x, fun)
!
! Evaluate the energy penalty due to steric interactions between water hydrogen atoms
! and nearby heavy atoms.
!
! Input:
! x:                Coordinates of the two hydrogen atoms
! O_coord:          Coordinates of the oxygen atom
! near(3, n_atoms): Coordinates of the n_atoms atoms near to oxygen atom
!
! On exit:
! fun:                Value of the energy-penalty, i.e., function to be minimized
!
use common_arrays_C,  only : coord
use orient_water_C, only : near, O_coord, n_atoms, near_nat, O_copy, H1_copy, atom_nos
implicit none
double precision, intent(in)   :: x(6)
double precision, intent(out)  :: fun
!
! Local
!
integer :: i, j
double precision :: h2(3,2), Rab, Rab2, min_O, c_O, angle
!
logical :: first = .true.
save
!
! This method of positioning the hydrogen atoms in a water molecule is "quick and dirty" -
! it works, but could do with being improved.
!
  if (first) then
!
!  Set the optimum H-X distances
!
    min_O = 2.d0
!
!  For a 12-6 potential, the precursor to the "6" term is 2/min_X**6
!
    c_O = 2.d0/min_O**6
    first = .false.
  end if
  h2(1:3,1) = x(1:3)
  h2(1:3,2) = x(4:6)
  fun = 0.d0
!
!  H - other atoms repulsion
!
  do j = 1, 2
    do i = 1, n_atoms
      Rab2 = (h2(1,j) - near(1,i))**2 + (h2(2,j) - near(2,i))**2 + (h2(3,j) - near(3,i))**2
      rab =  sqrt(rab2)
      select case (near_nat(i))
      case (1)  ! H-H repulsion
        fun = fun + 1.d0*Rab2**(-6)
      case (6, 7)  ! H-(C, and N) repulsion
        fun = fun + 100.d0*Rab2**(-6)
      case (8)  ! H-O interaction
        fun = fun + 100.d0*(Rab2**(-6) - c_O/Rab2**3)
        if (j == 1) then
          coord(:,H1_copy) = H2(:,1)
          call bangle (coord, O_copy, H1_copy, atom_nos(i), angle)
          fun = fun + 0.15d0*(3.14159d0 - angle)
        end if
      case (12)  ! H-Mg
        fun = fun + 200.d0*Rab2**(-6)
      case default  ! H-(everything else) repulsion
        fun = fun + 100.d0/Rab2**6
      end select
    end do
  end do
!
!  H - H distance
!
  Rab = sqrt((h2(1,1) - h2(1,2))**2 +(h2(2,1) - h2(2,2))**2 + (h2(3,1) - h2(3,2))**2)
  fun = fun + (Rab - 1.52d0)**2
!
!  O - H distance
!
  do i = 1, 2
    Rab = sqrt((h2(1,i) - O_coord(1))**2 +(h2(2,i) - O_coord(2))**2 + (h2(3,i) - O_coord(3))**2)
    fun = fun +    (Rab - 0.98d0)**2
  end do
return
end subroutine calcfc
subroutine cobyla (n,  x,  fun)
!
! Generic function minimized,  modified to work only for optimizing the positions of
! hydrogen atoms in water,  when ADD-H is used.
!
! n   = 6 = number of coordinates
! x   = Coordinates of the two hydrogen atoms,  updated on output
! fun = value of penalty (used in debugging only)
!
! This code is "black-box"
!
integer,  intent(in) :: n
double precision,  intent(in out)  :: x(6)
double precision :: con(2),  sim(n, n + 1),  simi(n, n),  datmat(2, n + 1),  a(n, 1),       &
             vsig(n),  veta(n),  sigbar(n),  dx(n),  w(n)
double precision :: alpha,  barmu,  beta,  cmin,  cmax,  cvmaxm,  cvmaxp,  delta,  denom,     &
             dxsign,  edgmax,  error,  fun,  gamma,  pareta,  parmu,  parsig,  phi,      &
             phimin,  prerec,  prerem,  ratio,  resmax,  resnew,  rho,  temp,  tempa,  &
             total,  trured,  vmnew,  vmold,  weta,  wsig,   rhobeg,  rhoend
integer   :: i,  ibrnch,  iflag,  ifull,  iptem,  iptemp,  j,  jdrop,  k,  l,  mp,   &
             nbest,  nfvals,  np,  iw = 6,  mpp,  m,  iprint,  maxfun
  rhobeg = 0.5d0
  rhoend = 1.d-3
  iprint = 0
  iflag = 0
  maxfun = 3500
  iptem = min(n, 5)
  iptemp = iptem + 1
  m = 0
  np = n + 1
  mp = m + 1
  mpp = m + 2
  alpha = 0.25d0
  beta = 2.1d0
  gamma = 0.5d0
  delta = 1.1d0
  rho = rhobeg
  parmu = 0.0d0
  parsig = 0.d0
  if (iprint >= 2) write(iw,  10) rho
  10 format (/'   The initial value of RHO is',  G13.6,    &
             '  and PARMU is set to zero.')
  nfvals = 0
  temp = 1.0d0/rho
  do i=1, n
    sim(i, np) = x(i)
    do j=1, n
      sim(i, j) = 0.0d0
      simi(i, j) = 0.0d0
    end do
    sim(i, i) = rho
    simi(i, i) = temp
  end do
  jdrop = np
  ibrnch = 0

  !  Make the next call of the user-supplied subroutine CALCFC. These
  !  instructions are also used for calling CALCFC during the iterations of
  !  the algorithm.

  40 if (nfvals >= maxfun .and. nfvals > 0) then
    if (iprint >= 1) write(iw,  50)
    50 format (/'   return from subroutine COBYLA because the ',   &
               'MAXFUN limit has been reached.')
    go to 600
  end if
  nfvals = nfvals + 1
  call calcfc (x,  fun)
  resmax = 0.0d0
  if (nfvals == iprint-1 .or. iprint == 3) then
    write(iw,  70) nfvals,  fun,  resmax,  x(1:iptem)
    70 format (/'   NFVALS = ',  i5,  '   fun = ',  G13.6,  '    MAXCV = ',   &
               G13.6/ ('   X = ',  5G14.6))
    if (iptem < n) write(iw,  80) x(iptemp:n)
    80 format (G19.6,  G15.6)
  end if
  con(mp) = fun
  con(mpp) = resmax
  if (ibrnch == 1) go to 440

  !  Set the recently calculated function values in a column of DATMAT. This
  !  array has a column for each vertex of the current simplex,  the entries of
  !  each column being the values of the constraint functions (if any)
  !  followed by the objective function and the greatest constraint violation
  !  at the vertex.

  do k=1, mpp
    datmat(k, jdrop) = con(k)
  end do
  if (nfvals > np) go to 130

  !  Exchange the new vertex of the initial simplex with the optimal vertex if
  !  necessary. then,  if the initial simplex is not complete,  pick its next
  !  vertex and calculate the function values there.

  if (jdrop <= n) then
    if (datmat(mp, np) <= fun) then
      x(jdrop) = sim(jdrop, np)
    else
      sim(jdrop, np) = x(jdrop)
      do k=1, mpp
        datmat(k, jdrop) = datmat(k, np)
        datmat(k, np) = con(k)
      end do
      do k=1, jdrop
        sim(jdrop, k) = -rho
        temp = -sum( simi(k:jdrop,  k) )
        simi(jdrop, k) = temp
      end do
    end if
  end if
  if (nfvals <= n) then
    jdrop = nfvals
    x(jdrop) = x(jdrop) + rho
    go to 40
  end if
  130 ibrnch = 1

  !  Identify the optimal vertex of the current simplex.

  140 phimin = datmat(mp, np) + parmu*datmat(mpp, np)
  nbest = np
  do j=1, n
    temp = datmat(mp, j) + parmu*datmat(mpp, j)
    if (temp < phimin) then
      nbest = j
      phimin = temp
    else if (abs(temp - phimin) < 1.d-6 .and. parmu == 0.0d0) then
      if (datmat(mpp, j) < datmat(mpp, nbest)) nbest = j
    end if
  end do

  !  Switch the best vertex into pole position if it is not there already,
  !  and also update SIM,  SIMI and DATMAT.

  if (nbest <= n) then
    do i=1, mpp
      temp = datmat(i, np)
      datmat(i, np) = datmat(i, nbest)
      datmat(i, nbest) = temp
    end do
    do i=1, n
      temp = sim(i, nbest)
      sim(i, nbest) = 0.0d0
      sim(i, np) = sim(i, np) + temp
      tempa = 0.0d0
      do k=1, n
        sim(i, k) = sim(i, k) - temp
        tempa = tempa - simi(k, i)
      end do
      simi(nbest, i) = tempa
    end do
  end if

  !  Make an error return if SIGI is a poor approximation to the inverse of
  !  the leading N by N submatrix of SIG.

  error = 0.0d0
  do i=1, n
    do j=1, n
      temp = 0.0d0
      if (i == j) temp = temp - 1.0d0
      temp = temp + dot_product( simi(i, 1:n),  sim(1:n, j) )
      error = max(error,  abs(temp))
    end do
  end do
  if (error > 0.1d0) then
    if (iprint >= 1) write(iw,  210)
    210 format (/'   return from subroutine COBYLA because ',   &
                'rounding errors are becoming damaging.')
    go to 600
  end if

  !  Calculate the coefficients of the linear approximations to the objective
  !  and constraint functions,  placing minus the objective function gradient
  !  after the constraint gradients in the array A. The vector W is used for
  !  working space.

  do k=1, mp
    con(k) = -datmat(k, np)
    do j=1, n
      w(j) = datmat(k, j) + con(k)
    end do
    do i=1, n
      temp = dot_product( w(1:n),  simi(1:n, i) )
      if (k == mp) temp = -temp
      a(i, k) = temp
    end do
  end do

  !  Calculate the values of sigma and eta,  and set IFLAG = 0 if the current
  !  simplex is not acceptable.

  iflag = 1
  parsig = alpha*rho
  pareta = beta*rho
  do j=1, n
    wsig = sum( simi(j, 1:n)**2 )
    weta = sum( sim(1:n, j)**2 )
    vsig(j) = 1.0d0/sqrt(wsig)
    veta(j) = sqrt(weta)
    if (vsig(j) < parsig .or. veta(j) > pareta) iflag = 0
  end do

  !  if a new vertex is needed to improve acceptability,  then decide which
  !  vertex to drop from the simplex.

  if (ibrnch == 1 .or. iflag == 1) go to 370
  jdrop = 0
  temp = pareta
  do j=1, n
    if (veta(j) > temp) then
      jdrop = j
      temp = veta(j)
    end if
  end do
  if (jdrop == 0) then
    do j=1, n
      if (vsig(j) < temp) then
        jdrop = j
        temp = vsig(j)
      end if
    end do
  end if

  !  Calculate the step to the new vertex and its sign.

  temp = gamma*rho*vsig(jdrop)
  dx(1:n) = temp*simi(jdrop, 1:n)
  cvmaxp = 0.0d0
  cvmaxm = 0.0d0
  do k=1, mp
    total = dot_product( a(1:n, k),  dx(1:n) )
    if (k < mp) then
      temp = datmat(k, np)
      cvmaxp = max(cvmaxp,  -total - temp)
      cvmaxm = max(cvmaxm,  total - temp)
    end if
  end do
  dxsign = 1.0d0
  if (parmu*(cvmaxp - cvmaxm) > total + total) dxsign = -1.0d0

  !  Update the elements of SIM and SIMI,  and set the next X.

  temp = 0.0d0
  do i=1, n
    dx(i) = dxsign*dx(i)
    sim(i, jdrop) = dx(i)
    temp = temp + simi(jdrop, i)*dx(i)
  end do
  simi(jdrop, 1:n) = simi(jdrop, 1:n) / temp
  do j=1, n
    if (j /= jdrop) then
      temp = dot_product( simi(j, 1:n),  dx(1:n) )
      simi(j, 1:n) = simi(j, 1:n) - temp*simi(jdrop, 1:n)
    end if
    x(j) = sim(j, np) + dx(j)
  end do
  go to 40

  !  Calculate DX = x(*)-x(0).
  !  Branch if the length of DX is less than 0.5*RHO.

  370 call trstlp (n,  m,  a,  con,  rho,  dx,  ifull)
  if (ifull == 0) then
    temp = sum( dx(1:n)**2 )
    if (temp < 0.25d0*rho*rho) then
      ibrnch = 1
      go to 550
    end if
  end if

  !  Predict the change to fun and the new maximum constraint violation if the
  !  variables are altered from x(0) to x(0) + DX.

  resnew = 0.0d0
  con(mp) = 0.0d0
  do k=1, mp
    total = con(k) - dot_product( a(1:n, k),  dx(1:n) )
    if (k < mp) resnew = max(resnew,  total)
  end do

  !  Increase PARMU if necessary and branch back if this change alters the
  !  optimal vertex. Otherwise PREREM and PREREC will be set to the predicted
  !  reductions in the merit function and the maximum constraint violation
  !  respectively.

  barmu = 0.0d0
  prerec = datmat(mpp, np) - resnew
  if (prerec > 0.0d0) barmu = total/prerec
  if (parmu < 1.5d0*barmu) then
    parmu = 2.0d0*barmu
    if (iprint >= 2) write(iw,  410) parmu
    410 format (/'   Increase in PARMU to',  G13.6)
    phi = datmat(mp, np) + parmu*datmat(mpp, np)
    do j=1, n
      temp = datmat(mp, j) + parmu*datmat(mpp, j)
      if (temp < phi) go to 140
      if (temp - phi < 1.d-6 .and. parmu == 0.0) then
        if (datmat(mpp, j) < datmat(mpp, np)) go to 140
      end if
    end do
  end if
  prerem = parmu*prerec - total

  !  Calculate the constraint and objective functions at x(*).
  !  then find the actual reduction in the merit function.

  x(1:n) = sim(1:n, np) + dx(1:n)
  ibrnch = 1
  go to 40

  440 vmold = datmat(mp, np) + parmu*datmat(mpp, np)
  vmnew = fun + parmu*resmax
  trured = vmold - vmnew
  if (abs(parmu) < 1.d-6 .and. abs(fun - datmat(mp, np)) < 1.d-6) then
    prerem = prerec
    trured = datmat(mpp, np) - resmax
  end if

  !  Begin the operations that decide whether x(*) should replace one of the
  !  vertices of the current simplex,  the change being mandatory if TRURED is
  !  positive. Firstly,  JDROP is set to the index of the vertex that is to be
  !  replaced.

  ratio = 0.0d0
  if (trured <= 0.0d0) ratio = 1.0d0
  jdrop = 0
  do j=1, n
    temp = dot_product( simi(j, 1:n),  dx(1:n) )
    temp = abs(temp)
    if (temp > ratio) then
      jdrop = j
      ratio = temp
    end if
    sigbar(j) = temp*vsig(j)
  end do

  !  Calculate the value of ell.

  edgmax = delta*rho
  l = 0
  do j=1, n
    if (sigbar(j) >= parsig .or. sigbar(j) >= vsig(j)) then
      temp = veta(j)
      if (trured > 0.0d0) then
        temp = sum( (dx(1:n) - sim(1:n, j))**2 )
        temp = sqrt(temp)
      end if
      if (temp > edgmax) then
        l = j
        edgmax = temp
      end if
    end if
  end do
  if (l > 0) jdrop = l
  if (jdrop == 0) go to 550

  !  Revise the simplex by updating the elements of SIM,  SIMI and DATMAT.

  temp = 0.0d0
  do i=1, n
    sim(i, jdrop) = dx(i)
    temp = temp + simi(jdrop, i)*dx(i)
  end do
  simi(jdrop, 1:n) = simi(jdrop, 1:n) / temp
  do j=1, n
    if (j /= jdrop) then
      temp = dot_product( simi(j, 1:n),  dx(1:n) )
      simi(j, 1:n) = simi(j, 1:n) - temp*simi(jdrop, 1:n)
    end if
  end do
  datmat(1:mpp, jdrop) = con(1:mpp)

  !  Branch back for further iterations with the current RHO.

  if (trured > 0.0d0 .and. trured >= 0.1d0*prerem) go to 140
  550 if (iflag == 0) then
    ibrnch = 0
    go to 140
  end if

  !  Otherwise reduce RHO if it is not at its least value and reset PARMU.

  if (rho > rhoend) then
    rho = 0.5d0*rho
    if (rho <= 1.5d0*rhoend) rho = rhoend
    if (parmu > 0.0d0) then
      denom = 0.0d0
      do k=1, mp
        cmin = datmat(k, np)
        cmax = cmin
        do i=1, n
          cmin = min(cmin,  datmat(k, i))
          cmax = max(cmax,  datmat(k, i))
        end do
        if (k <= m .and. cmin < 0.5d0*cmax) then
          temp = max(cmax, 0.0d0) - cmin
          if (denom <= 0.0d0) then
            denom = temp
          else
            denom = min(denom, temp)
          end if
        end if
      end do
      if (denom == 0.0d0) then
        parmu = 0.0d0
      else if (cmax - cmin < parmu*denom) then
        parmu = (cmax - cmin)/denom
      end if
    end if
    if (iprint >= 2) write(iw,  580) rho, parmu
    580 format ('   Reduction in RHO to ',  G13.6,  '  and PARMU = ',  G13.6)
    if (iprint == 2) then
      write(iw,  70) nfvals,  datmat(mp, np),  datmat(mpp, np),  sim(1:iptem, np)
      if (iptem < n) write(iw,  80) x(iptemp:n)
    end if
    go to 140
  end if

  !  return the best calculated values of the variables.

  if (iprint >= 1) write(iw,  590)
  590 format (/'   Normal return from subroutine COBYLA')
  if (ifull == 1) go to 620

  600 x(1:n) = sim(1:n, np)
  fun = datmat(mp, np)
  resmax = datmat(mpp, np)
  620 if (iprint >= 1) then
    write(iw,  70) nfvals,  fun,  resmax,  x(1:iptem)
    if (iptem < n) write(iw,  80) x(iptemp:n)
  end if
  maxfun = nfvals
  return
end subroutine cobyla
!------------------------------------------------------------------------------

subroutine trstlp (n,  m,  a,  b,  rho,  dx,  ifull)

! N.B. Arguments Z,  ZDOTA,  VMULTC,  SDIRN,  DXNEW,  VMULTD & IACT have been removed.

integer,  intent(in)     :: n
integer,  intent(in)     :: m
double precision,  intent(in)   :: a(n, 10)
double precision,  intent(in)   :: b(m)
double precision,  intent(in)   :: rho
double precision,  intent(out)  :: dx(10)
integer,  intent(out)    :: ifull

!

double precision :: z(n, n),  zdota(m + 1),  vmultc(m + 1),  sdirn(n),  dxnew(n),  vmultd(m + 1)
double precision :: acca,  accb,  alpha,  beta,  dd,  optnew,  optold,  ratio,  resmax,    &
             resold,  sd,  sp,  spabs,  ss,  step,  stpful,  sumabs,  temp,  tempa,  &
             tot,  total,  vsave,  zdotv,  zdotw,  zdvabs,  zdwabs
integer   :: i,  iact(m + 1),  icon,  icount,  isave,  k,  kk,  kl,  kp,  kw,  mcon,    &
             nact,  nactx

ifull = 1
icon = 0
mcon = m
nact = 0
nactx = 0
resold = 0.d0
resmax = 0.0d0
do i=1, n
  z(i, 1:n) = 0.0d0
  z(i, i) = 1.0d0
  dx(i) = 0.0d0
end do
if (m >= 1) then
  do k=1, m
    if (b(k) > resmax) then
      resmax = b(k)
      icon = k
    end if
  end do
  do k=1, m
    iact(k) = k
    vmultc(k) = resmax - b(k)
  end do
end if
if (resmax == 0.0d0) go to 480
sdirn(1:n) = 0.0d0

!  end the current stage of the calculation if 3 consecutive iterations
!  have either failed to reduce the best calculated value of the objective
!  function or to increase the number of active constraints since the best
!  value was calculated. This strategy prevents cycling,  but there is a
!  remote possibility that it will cause premature termination.

60 optold = 0.0d0
icount = 0
70 if (mcon == m) then
  optnew = resmax
else
  optnew = -dot_product( dx(1:n),  a(1:n, mcon) )
end if
if (icount == 0 .or. optnew < optold) then
  optold = optnew
  nactx = nact
  icount = 3
else if (nact > nactx) then
  nactx = nact
  icount = 3
else
  icount = icount - 1
  if (icount == 0) go to 490
end if

!  if ICON exceeds NACT,  then we add the constraint with index IACT(ICON) to
!  the active set. Apply Givens rotations so that the last N-NACT-1 columns
!  of Z are orthogonal to the gradient of the new constraint,  a scalar
!  product being set to zero if its nonzero value could be due to computer
!  rounding errors. The array DXNEW is used for working space.

if (icon <= nact) go to 260
kk = iact(icon)
dxnew(1:n) = a(1:n, kk)
tot = 0.0d0
k = n
100 if (k > nact) then
  sp = 0.0d0
  spabs = 0.0d0
  do i=1, n
    temp = z(i, k)*dxnew(i)
    sp = sp + temp
    spabs = spabs + abs(temp)
  end do
  acca = spabs + 0.1d0*abs(sp)
  accb = spabs + 0.2d0*abs(sp)
  if (spabs >= acca .or. acca >= accb) sp = 0.0d0
  if (tot == 0.0d0) then
    tot = sp
  else
    kp = k + 1
    temp = sqrt(sp*sp + tot*tot)
    alpha = sp/temp
    beta = tot/temp
    tot = temp
    do i=1, n
      temp = alpha*z(i, k) + beta*z(i, kp)
      z(i, kp) = alpha*z(i, kp) - beta*z(i, k)
      z(i, k) = temp
    end do
  end if
  k = k - 1
  go to 100
end if

!  Add the new constraint if this can be done without a deletion from the
!  active set.

if (tot /= 0.0d0) then
  nact = nact + 1
  zdota(nact) = tot
  vmultc(icon) = vmultc(nact)
  vmultc(nact) = 0.0d0
  go to 210
end if

!  The next instruction is reached if a deletion has to be made from the
!  active set in order to make room for the new active constraint,  because
!  the new constraint gradient is a linear combination of the gradients of
!  the old active constraints.  Set the elements of VMULTD to the multipliers
!  of the linear combination.  Further,  set IOUT to the index of the
!  constraint to be deleted,  but branch if no suitable index can be found.

ratio = -1.0d0
k = nact
130 zdotv = 0.0d0
zdvabs = 0.0d0
do i=1, n
  temp = z(i, k)*dxnew(i)
  zdotv = zdotv + temp
  zdvabs = zdvabs + abs(temp)
end do
acca = zdvabs + 0.1d0*abs(zdotv)
accb = zdvabs + 0.2d0*abs(zdotv)
if (zdvabs < acca .and. acca < accb) then
  temp = zdotv/zdota(k)
  if (temp > 0.0d0 .and. iact(k) <= m) then
    tempa = vmultc(k)/temp
    if (ratio < 0.0d0 .or. tempa < ratio) then
      ratio = tempa
    end if
  end if
  if (k >= 2) then
    kw = iact(k)
    dxnew(1:n) = dxnew(1:n) - temp*a(1:n, kw)
  end if
  vmultd(k) = temp
else
  vmultd(k) = 0.0d0
end if
k = k - 1
if (k > 0) go to 130
if (ratio < 0.0d0) go to 490

!  Revise the Lagrange multipliers and reorder the active constraints so
!  that the one to be replaced is at the end of the list. Also calculate the
!  new value of ZDOTA(NACT) and branch if it is not acceptable.

do k=1, nact
  vmultc(k) = max(0.0d0, vmultc(k) - ratio*vmultd(k))
end do
if (icon < nact) then
  isave = iact(icon)
  vsave = vmultc(icon)
  k = icon
  170 kp = k + 1
  kw = iact(kp)
  sp = dot_product( z(1:n, k),  a(1:n, kw) )
  temp = sqrt(sp*sp + zdota(kp)**2)
  alpha = zdota(kp)/temp
  beta = sp/temp
  zdota(kp) = alpha*zdota(k)
  zdota(k) = temp
  do i=1, n
    temp = alpha*z(i, kp) + beta*z(i, k)
    z(i, kp) = alpha*z(i, k) - beta*z(i, kp)
    z(i, k) = temp
  end do
  iact(k) = kw
  vmultc(k) = vmultc(kp)
  k = kp
  if (k < nact) go to 170
  iact(k) = isave
  vmultc(k) = vsave
end if
temp = dot_product( z(1:n, nact),  a(1:n, kk) )
if (temp == 0.0d0) go to 490
zdota(nact) = temp
vmultc(icon) = 0.0d0
vmultc(nact) = ratio

!  Update IACT and ensure that the objective function continues to be
!  treated as the last active constraint when MCON>M.

210 iact(icon) = iact(nact)
iact(nact) = kk
if (mcon > m .and. kk /= mcon) then
  k = nact - 1
  sp = dot_product( z(1:n, k),  a(1:n, kk) )
  temp = sqrt(sp*sp + zdota(nact)**2)
  alpha = zdota(nact)/temp
  beta = sp/temp
  zdota(nact) = alpha*zdota(k)
  zdota(k) = temp
  do i=1, n
    temp = alpha*z(i, nact) + beta*z(i, k)
    z(i, nact) = alpha*z(i, k) - beta*z(i, nact)
    z(i, k) = temp
  end do
  iact(nact) = iact(k)
  iact(k) = kk
  temp = vmultc(k)
  vmultc(k) = vmultc(nact)
  vmultc(nact) = temp
end if

!  if stage one is in progress,  then set SDIRN to the direction of the next
!  change to the current vector of variables.

if (mcon > m) go to 320
kk = iact(nact)
temp = dot_product( sdirn(1:n),  a(1:n, kk) )
temp = temp - 1.0d0
temp = temp/zdota(nact)
sdirn(1:n) = sdirn(1:n) - temp*z(1:n, nact)
go to 340

!  Delete the constraint that has the index IACT(ICON) from the active set.

260 if (icon < nact) then
  isave = iact(icon)
  vsave = vmultc(icon)
  k = icon
  do
    kp = k + 1
    kk = iact(kp)
    sp = dot_product( z(1:n, k),  a(1:n, kk) )
    temp = sqrt(sp*sp + zdota(kp)**2)
    alpha = zdota(kp)/temp
    beta = sp/temp
    zdota(kp) = alpha*zdota(k)
    zdota(k) = temp
    do i=1, n
      temp = alpha*z(i, kp) + beta*z(i, k)
      z(i, kp) = alpha*z(i, k) - beta*z(i, kp)
      z(i, k) = temp
    end do
    iact(k) = kk
    vmultc(k) = vmultc(kp)
    k = kp
    if (k >= nact) EXIT
  end do
  iact(k) = isave
  vmultc(k) = vsave
end if
nact = nact - 1

!  if stage one is in progress,  then set SDIRN to the direction of the next
!  change to the current vector of variables.

if (mcon > m) go to 320
temp = dot_product( sdirn(1:n),  z(1:n, nact + 1) )
sdirn(1:n) = sdirn(1:n) - temp*z(1:n, nact + 1)
go to 340

!  Pick the next search direction of stage two.

320 temp = 1.0d0/zdota(nact)
sdirn(1:n) = temp*z(1:n, nact)

!  Calculate the step to the boundary of the trust region or take the step
!  that reduces RESMAX to zero. The two statements below that include the
!  factor 1.0E-6 prevent some harmless underflows that occurred in a test
!  calculation. Further,  we skip the step if it could be zero within a
!  reasonable tolerance for computer rounding errors.

340 dd = rho*rho
sd = 0.0d0
ss = 0.0d0
do i=1, n
  if (abs(dx(i)) >= 1.0E-6*rho) dd = dd - dx(i)**2
  sd = sd + dx(i)*sdirn(i)
  ss = ss + sdirn(i)**2
end do
if (dd <= 0.0d0) go to 490
temp = sqrt(ss*dd)
if (abs(sd) >= 1.0E-6*temp) temp = sqrt(ss*dd + sd*sd)
stpful = dd/(temp + sd)
step = stpful
if (mcon == m) then
  acca = step + 0.1d0*resmax
  accb = step + 0.2d0*resmax
  if (step >= acca .or. acca >= accb) go to 480
  step = min(step, resmax)
end if

!  Set DXNEW to the new variables if STEP is the steplength,  and reduce
!  RESMAX to the corresponding maximum residual if stage one is being done.
!  Because DXNEW will be changed during the calculation of some Lagrange
!  multipliers,  it will be restored to the following value later.

dxnew(1:n) = dx(1:n) + step*sdirn(1:n)
if (mcon == m) then
  resold = resmax
  resmax = 0.0d0
  do k=1, nact
    kk = iact(k)
    temp = b(kk) - dot_product( a(1:n, kk),  dxnew(1:n) )
    resmax = max(resmax, temp)
  end do
end if

!  Set VMULTD to the VMULTC vector that would occur if DX became DXNEW. A
!  device is included to force VMULTD(K) = 0.0 if deviations from this value
!  can be attributed to computer rounding errors. First calculate the new
!  Lagrange multipliers.

k = nact
390 zdotw = 0.0d0
zdwabs = 0.0d0
do i=1, n
  temp = z(i, k)*dxnew(i)
  zdotw = zdotw + temp
  zdwabs = zdwabs + abs(temp)
end do
acca = zdwabs + 0.1d0*abs(zdotw)
accb = zdwabs + 0.2d0*abs(zdotw)
if (zdwabs >= acca .or. acca >= accb) zdotw = 0.0d0
vmultd(k) = zdotw / zdota(k)
if (k >= 2) then
  kk = iact(k)
  dxnew(1:n) = dxnew(1:n) - vmultd(k)*a(1:n, kk)
  k = k - 1
  go to 390
end if
if (mcon > m) vmultd(nact) = max(0.0d0, vmultd(nact))

!  Complete VMULTC by finding the new constraint residuals.

dxnew(1:n) = dx(1:n) + step*sdirn(1:n)
if (mcon > nact) then
  kl = nact + 1
  do k=kl, mcon
    kk = iact(k)
    total = resmax - b(kk)
    sumabs = resmax + abs(b(kk))
    do i=1, n
      temp = a(i, kk)*dxnew(i)
      total = total + temp
      sumabs = sumabs + abs(temp)
    end do
    acca = sumabs + 0.1*abs(total)
    accb = sumabs + 0.2*abs(total)
    if (sumabs >= acca .or. acca >= accb) total = 0.0
    vmultd(k) = total
  end do
end if

!  Calculate the fraction of the step from DX to DXNEW that will be taken.

ratio = 1.0d0
icon = 0
do k=1, mcon
  if (vmultd(k) < 0.0d0) then
    temp = vmultc(k)/(vmultc(k) - vmultd(k))
    if (temp < ratio) then
      ratio = temp
      icon = k
    end if
  end if
end do

!  Update DX,  VMULTC and RESMAX.

temp = 1.0d0 - ratio
dx(1:n) = temp*dx(1:n) + ratio*dxnew(1:n)
do k=1, mcon
  vmultc(k) = max(0.0d0, temp*vmultc(k) + ratio*vmultd(k))
end do
if (mcon == m) resmax = resold + ratio*(resmax - resold)

!  if the full step is not acceptable then begin another iteration.
!  Otherwise switch to stage two or end the calculation.

if (icon > 0) go to 70
if (abs(step - stpful) < 1.d-6) go to 500
480 mcon = m + 1
icon = mcon
iact(mcon) = mcon
vmultc(mcon) = 0.0d0
go to 60

!  We employ any freedom that may be available to reduce the objective
!  function before returning a DX whose length is less than RHO.

490 if (mcon == m) go to 480
ifull = 0

500 return
end subroutine trstlp
