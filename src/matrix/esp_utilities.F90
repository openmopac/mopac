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

   subroutine rys(x, nroots, u, w)
    double precision :: x, u(4), w(4)
    integer nroots
    double precision :: &
    pie4 = 7.85398163397448d-01, &
    r12  = 2.75255128608411d-01, &
    w22  = 9.17517095361369d-02, &
    r22  = 2.72474487139158d+00, &
    Ex, f1, y
    u = 0.d0
    w = 0.d0
    if (x > 33.d0) then
      w(1) = sqrt(pie4/x)
      if (nroots == 1) then
        u(1) = 0.5d0/(x - 0.5d0)
        w(1) = sqrt(pie4/x)
      else
        if (x > 40.d0) then
          u(1) = r12/(x-r12);
          u(2) = r22/(x-r22);
          w(2) = w22*w(1);
          w(1) = w(1)-w(2);
        else
          Ex =  exp(-x);
          u(1)=(-8.78947307498880d-01*x + 1.09243702330261d+01)*Ex + r12/(x-r12)
          u(2)=(-9.28903924275977d+00*x + 8.10642367843811d+01)*Ex + r22/(x-r22)
          w(2)=( 4.46857389308400d+00*x - 7.79250653461045d+01)*Ex + w22*w(1)
          w(1)=w(1)-w(2)
        end if
      end if
    else if (x > 15.d0) then
      ex = exp(-x)
      w(1) = ((( 1.9623264149430d-01/x -4.9695241464490d-01)/x - &
      6.0156581186481d-05)*ex + sqrt(pie4/x))
      f1 =(w(1) - ex)/(2.d0*x)
      if (nroots == 1) then
        u(1) = f1/(w(1) - f1)
      else
        u(1) = ((((-1.14906395546354d-06*x + 1.76003409708332d-04)*x - &
            1.71984023644904d-02)*x -1.37292644149838d-01)*x + &
            (-4.75742064274859d+01/x + 9.21005186542857d+00)/x - &
            2.31080873898939d-02)*Ex + r12/(x-r12)
        u(2) = ((( 3.64921633404158d-04*x - 9.71850973831558d-02)*x - &
            4.02886174850252d+00)*x + &
            (-1.35831002139173d+02/x -8.66891724287962d+01)/x + &
            2.98011277766958d+00)*Ex + r22/(x-r22)
         w(2) = ((f1 - w(1)) * u(1) + f1)*(1.0d0 + u(2))/(u(2) - u(1))
         w(1) = w(1)-w(2)
      end if
    else if (x > 10.d0) then
      ex = exp(-x)
      w(1) =((((-1.8784686463512d-01/x + 2.2991849164985d-01)/x - &
      4.9893752514047d-01)/x  -2.1916512131607d-05)*ex + &
      sqrt(pie4/x))
      f1 =(w(1) - ex)/(2.d0*x)
      if (nroots == 1) then
        u(1) = f1/(w(1) - f1)
      else
        u(1) = ((((-1.01041157064226d-05*x &
              +1.19483054115173d-03)*x &
               -6.73760231824074d-02)*x &
               +1.25705571069895d+00)*x &
               +(((-8.57609422987199d+03/x &
               +5.91005939591842d+03)/x &
            -1.70807677109425d+03)/x &
            +2.64536689959503d+02)/x &
            -2.38570496490846d+01)*Ex  &
            + r12/(x-r12)
         u(2) = ((( 3.39024225137123d-04*x &
             -9.34976436343509d-02)*x &
            -4.22216483306320d+00)*x  &
            +(((-2.08457050986847d+03/x &
            -1.04999071905664d+03)/x &
            +3.39891508992661d+02)/x &
            -1.56184800325063d+02)/x &
            +8.00839033297501d+00)*Ex  &
            + r22/(x-r22)
         w(2) = ((f1 - w(1))*u(1) + f1)*(1.d0 + u(2))/(u(2) - u(1))
         w(1) = w(1) - w(2)
      end if
    else if (x > 5.d0) then
      ex = exp(-x)
      w(1)= ((((((( 4.6897511375022d-01/x-6.9955602298985d-01)/x + &
      5.3689283271887d-01)/x-3.2883030418398d-01)/x - &
      2.4645596956002d-01)/x-4.9984072848436d-01)/x - 3.1501078774085d-06)*ex + &
      sqrt(pie4/x))
      f1 =(w(1) - ex)/(2.d0*x)
      if (nroots == 1) then
        u(1) = f1/(w(1) - f1)
      else
        y=x-7.5d0
        u(1)= (((((((((((((-1.43632730148572d-16*y+2.38198922570405d-16)*y + &
        1.358319618800d-14)*y-7.064522786879d-14)*y &
      -7.719300212748d-13 &
          )*y+7.802544789997d-12)*y+6.628721099436d-11)*y &
      -1.775564159743d-09 &
          )*y+1.713828823990d-08)*y-1.497500187053d-07)*y &
      +2.283485114279d-06 &
          )*y-3.76953869614706d-05)*y+4.74791204651451d-04 &
          )*y-4.60448960876139d-03)*y+3.72458587837249d-02
         u(2)= (((((((((((( 2.48791622798900d-14*y-1.36113510175724d-13 &
          )*y-2.224334349799d-12)*y+4.190559455515d-11)*y &
      -2.222722579924d-10 &
          )*y-2.624183464275d-09)*y+6.128153450169d-08)*y &
      -4.383376014528d-07 &
          )*y-2.49952200232910d-06)*y+1.03236647888320d-04 &
          )*y-1.44614664924989d-03)*y+1.35094294917224d-02 &
          )*y-9.53478510453887d-02)*y+5.44765245686790d-01
         w(2)=((f1-w(1))*u(1)+f1)*(1.0d0+u(2))/(u(2)-u(1))
         w(1)=w(1)-w(2)
      end if
    else if (x > 3.d0) then
      y = x - 4.d0;
      f1= ((((((((((-2.62453564772299d-11*y + 3.24031041623823d-10)*y - &
      3.614965656163d-09)*y + 3.760256799971d-08)*y - 3.553558319675d-07)*y + &
      3.022556449731d-06)*y - 2.290098979647d-05)*y + 1.526537461148d-04)*y - &
      8.81947375894379d-04)*y + 4.33207949514611d-03)*y - &
      1.75257821619926d-02)*y + 5.28406320615584d-02
      w(1) = 2.d0*x*f1 + exp(-x)
      if (nroots == 1) then
        w(1) = 2.d0*x*f1 + exp(-x)
        u(1) = f1/(w(1)-f1)
      else
        u(1)= ((((((((-4.11560117487296d-12*y+7.10910223886747d-11)*y &
        -1.73508862390291d-09)*y+5.93066856324744d-08)*y &
        -9.76085576741771d-07)*y+1.08484384385679d-05)*y &
        -1.12608004981982d-04)*y+1.16210907653515d-03)*y &
        -9.89572595720351d-03)*y+6.12589701086408d-02
           u(2)= (((((((((-1.80555625241001d-10*y+5.44072475994123d-10)*y &
        +1.603498045240d-08)*y-1.497986283037d-07)*y-7.017002532106d-07)*y &
        +1.85882653064034d-05)*y-2.04685420150802d-05)*y &
        -2.49327728643089d-03)*y+3.56550690684281d-02)*y &
        -2.60417417692375d-01)*y+1.12155283108289d+00
        w(2)=((f1-w(1))*u(1)+f1)*(1.0d0+u(2))/(u(2)-u(1))
        w(1)=w(1) - w(2)
      end if
    else if (x > 1.d0) then
      y=x-2.0d0
      f1= ((((((((((-1.61702782425558d-10*y + 1.96215250865776d-09)*y - &
      2.14234468198419d-08)*y + 2.17216556336318d-07)*y - &
      1.98850171329371d-06)*y + 1.62429321438911d-05)*y - &
      1.16740298039895d-04)*y + 7.24888732052332d-04)*y - &
      3.79490003707156d-03)*y + 1.61723488664661d-02)*y - &
      5.29428148329736d-02)*y + 1.15702180856167d-01
      w(1) = 2.d0*x*f1 + exp(-x)
      if (nroots == 1) then
        u(1) = f1/(w(1)-f1)
      else
        u(1) = (((((((((-6.36859636616415d-12*y+8.47417064776270d-11) * y &
        -5.152207846962d-10)*y-3.846389873308d-10)*y+8.472253388380d-08) * y &
        -1.85306035634293d-06)*y+2.47191693238413d-05) * y &
        -2.49018321709815d-04)*y+2.19173220020161d-03) * y &
        -1.63329339286794d-02)*y+8.68085688285261d-02
          u(2) = (((((((((1.45331350488343d-10*y+2.07111465297976d-09) * y &
        -1.878920917404d-08)*y - 1.725838516261d-07)*y+2.247389642339d-06) * y &
        +9.76783813082564d-06)*y - 1.93160765581969d-04) * y &
        -1.58064140671893d-03)*y + 4.85928174507904d-02) * y &
        -4.30761584997596d-01)*y + 1.80400974537950d+00
          w(2) = ((f1-w(1))*u(1)+f1)*(1.0d0+u(2))/(u(2)-u(1))
          w(1) = w(1)-w(2)
      end if
    else if (x > 3.d-7) then
      f1= ((((((((-8.36313918003957d-08*x + 1.21222603512827d-06)*x - &
      1.15662609053481d-05)*x + 9.25197374512647d-05)*x - &
      6.40994113129432d-04)*x + 3.78787044215009d-03)*x - &
      1.85185172458485d-02)*x + 7.14285713298222d-02)*x - &
      1.99999999997023d-01)*x + 3.33333333333318d-01
      w(1) = 2.d0*x*f1 + exp(-x)
      if (nroots == 1) then
        u(1) = f1/(w(1)-f1)
      else
        u(1) = (((((((-2.35234358048491d-09*x + 2.49173650389842d-08) &
             *x - 4.558315364581d-08)*x - 2.447252174587d-06)*x + 4.743292959463d-05) &
             *x - 5.33184749432408d-04)*x + 4.44654947116579d-03) &
             *x - 2.90430236084697d-02)*x + 1.30693606237085d-01
            u(2) = (((((((-2.47404902329170d-08*x + 2.36809910635906d-07) &
             *x + 1.835367736310d-06)*x - 2.066168802076d-05)*x - 1.345693393936d-04) &
             *x - 5.88154362858038d-05)*x + 5.32735082098139d-02) &
             *x - 6.37623643056745d-01)*x + 2.86930639376289d+00
            w(2)=((f1-w(1))*u(1)+f1)*(1.0d0+u(2))/(u(2)-u(1))
            w(1)=w(1) - w(2)
      end if
    else
!
!  NOT verified!
!
      if (nroots == 1) then
        u(1) = 0.5d0 - x/5.0d0
        w(1) = 1.0d0 - x/3.0d0
      else
        u(1) = 1.30693606237085d-01 - 2.90430236082028d-02*x
        u(2) = 2.86930639376291d+00 - 6.37623643058102d-01*x
        w(1) = 6.52145154862545d-01 - 1.22713621927067d-01*x
        w(2) = 3.47854845137453d-01 - 2.10619711404725d-01*x
      end if
    end if
    return
  end subroutine rys

  subroutine vint(xint, yint, zint, ni, nj, x0, y0, z0, xi, yi, zi, xj, yj, zj, t)
!
!  Gauss-Hermite Quadrature
!
    double precision, intent (in) :: x0, y0, z0, xi, yi, zi, xj, yj, zj, t
    double precision, intent (out) :: xint, yint, zint
    integer, intent (in) :: ni, nj
!
!  Local variables
!
    double precision :: dum, ptx, pty, ptz, ax, ay, az, bx, by, bz, px, py, pz
    integer ::  i, j, min, max, npt
    double precision :: Hermite_roots(21) = (/ &
                   0.0d0,                                                      & ! 1
                  -0.7071067811865d0,  0.7071067811865d0,                      & ! 2
                  -1.224744871391d0,   0.0d0,              1.2247448713915d0,  & ! 3
                  -1.6506801238857d0, -0.52464762327529d0, 0.52464762327529d0, & ! 4
                     1.6506801238857d0,                                        &
                  -2.020182870456d0,  -0.9585724646138d0,  0.0d0,              & ! 5
                        0.9585724646138d0,  2.020182870456d0,                  &
                  -2.350604973674d0,  -1.335849074014d0,  -0.436077411928d0,   & ! 6
                        0.436077411928d0,   1.335849074014d0,   2.350604973674d0 /)

    double precision ::   Hermite_weights(21)  = (/ &
            1.772453850905d0,                                                  & ! 1
            0.886226925452d0,     0.886226925452d0,                            & ! 2
            0.295408975150d0,     1.18163590060d0,      0.295408975150d0,      & ! 3
            8.131283544725d-02, 0.8049140900055d0,    8.049140900055d-01,      & ! 4
              8.131283544725d-02,                                              &
            1.995324205905d-02, 3.936193231522d-01, 9.453087204829d-01,        & ! 5
              3.936193231522d-01, 1.995324205905d-02,                          &
            4.530009905509d-03, 1.570673203229d-01, 7.246295952244d-01,        & ! 6
            7.246295952244d-01, 1.570673203229d-01, 4.530009905509d-03 /)
    integer :: start_index(6) = (/ 0, 1, 3, 6, 10, 15/), &
                 end_index(6) = (/ 0, 2, 5, 9, 14, 20/)


    xint = 0.0d0
    yint = 0.0d0
    zint = 0.0d0

    npt  = (ni + nj)/2 + 1      ! form the index into the start_index array
    min  = start_index(npt) + 1 ! get the starting index in the roots,weights arrays.
    max  = end_index(npt) + 1   ! get the ending index in the roots,weights arrays.

    do i = min, max
      px   = 1.0d0
      py   = 1.0d0
      pz   = 1.0d0
      dum  = Hermite_roots(i) * t
      ptx  = dum + x0
      pty  = dum + y0
      ptz  = dum + z0
      ax   = ptx - xi
      ay   = pty - yi
      az   = ptz - zi
      bx   = ptx - xj
      by   = pty - yj
      bz   = ptz - zj

      do j = 2, ni
        px = px*ax
        py = py*ay
        pz = pz*az
      end do

      do j = 2, nj
        px = px*bx
        py = py*by
        pz = pz*bz
      end do
      dum    = Hermite_weights(i)
      xint = xint + dum * px
      yint = yint + dum * py
      zint = zint + dum * pz
    end do
  end subroutine vint

  subroutine get_minus_point_five_overlap(s)
  use molkst_C, only : numat, norbs
  use common_arrays_C, only: nfirst, nlast, nat, h, &
                           & i1fact
  use parameters_C, only : betas, betap, betad
  implicit none
  double precision, dimension(norbs**2) :: s
  !
  !.. Local Scalars ..
  integer :: i, if, ii, ij, il, im1, j, jf, jj, jl, k, &
       & alloc_stat
  double precision :: bi, bj, sum
  !.. Local Arrays ..
  double precision, dimension (:), allocatable :: eigs, vecs
!

  !
  !
  !.. Intrinsic Functions ..
  intrinsic Abs, Sqrt
  !
  ! ... Executable Statements ...
  !
  !     PATAS
  !*********************************************************************
  !
  !  FIRST, RE-CALCULATE THE OVERLAP MATRIX
  !
  !*********************************************************************
  !
  !     PATAS

  allocate (eigs(norbs), vecs(norbs**2), stat=alloc_stat)
  if (alloc_stat /= 0) return  !WARNING
  if (.not. allocated(i1fact)) then

    allocate(i1fact(3 + norbs))
!
!   SET UP ARRAY OF LOWER HALF TRIANGLE INDICES (PASCAL'S TRIANGLE)
!
    do i = 1, norbs
      i1fact(i) = (i*(i + 1))/2
    end do
  end if
  do i = 1, numat
      if = nfirst(i)
      il = nlast(i)
      if (il >= if) then
        eigs(if) = betas(nat(i))
        if (il > if) then
          eigs(if+1) = betap(nat(i))
          eigs(if+2) = eigs(if+1)
          eigs(if+3) = eigs(if+1)
          if (il > if+3) then
            eigs(if+4) = betad(nat(i))
            eigs(if+5) = eigs(if+4)
            eigs(if+6) = eigs(if+4)
            eigs(if+7) = eigs(if+4)
            eigs(if+8) = eigs(if+4)
          end if
        end if
    else
      if = nfirst(i)
      il = nlast(i)
      if (il >= if) then
        eigs(if) = betas(nat(i))
        if (il > if) then
          eigs(if+1) = betap(nat(i))
          eigs(if+2) = eigs(if+1)
          eigs(if+3) = eigs(if+1)
          if (il > if+3) then
            eigs(if+4) = betad(nat(i))
            eigs(if+5) = eigs(if+4)
            eigs(if+6) = eigs(if+4)
            eigs(if+7) = eigs(if+4)
            eigs(if+8) = eigs(if+4)
          end if
        end if
      end if
    end if
    im1 = i - 1
    bi = betas(nat(i))
    do k = if, il
      bi = eigs(k)
      ii = (k*(k-1)) / 2
      do j = 1, im1
        jf = nfirst(j)
        jl = nlast(j)
        do jj = jf, jl
          bj = eigs(jj)
          ij = ii + jj
          h(ij) = 2.d0 * h(ij) / (bi+bj) + 1.d-14
          !  THE  +1.D-14 IS TO PREVENT POSSIBLE ERRORS IN THE DIAGONALIZATION.
        end do
      end do
      do jj = if, k
        ij = ii + jj
        h(ij) = 0.d0
      end do
    end do
  end do
  do i = 1, norbs
    h(i1fact(i)) = 1.d0
  end do
  call rsp (h, norbs, eigs, vecs)
  do i = 1, norbs
    eigs(i) = 1.d0 / Sqrt (Abs (eigs(i)))
  end do
  ij = 0
  do i = 1, norbs
    do j = 1, i
      ij = ij + 1
      sum = 0.d0
      do k = 1, norbs
        sum = sum + vecs(i+ (k-1)*norbs) * eigs(k) * vecs &
             & (j+ (k-1)*norbs)
      end do
      s(i+ (j-1)*norbs) = sum
      s(j+ (i-1)*norbs) = sum
    end do
  end do
  !
  end subroutine get_minus_point_five_overlap


  subroutine evec(aVector, x, y, z, coord, numat)
    integer :: numat
    double precision :: coord(3,numat), x, y, z
    real :: aVector(7*numat)
!
!
    integer :: i, j
    real :: shellRi, shellR2, shellR2i, shellR3i, &
    u, v, w, uu, vv, ww
    j = 1
    do i = 1, numat
      u = sngl(x - coord(1,i))
      v = sngl(y - coord(2,i))
      w = sngl(z - coord(3,i))
      uu = u * u
      vv = v * v
      ww = w * w
      shellR2 = max(uu + vv + ww, 1.e-2)
      shellR2i = 1.e0/(shellR2 + 1.e-7)
      shellRi = sqrt(shellR2i)

      aVector(j + 0) =  shellRi
      shellR3i = shellRi * shellR2i
      aVector(j + 1) = u * shellR3i
      aVector(j + 2) = v * shellR3i
      aVector(j + 3) = w * shellR3i
      aVector(j + 4) = shellR2i
      aVector(j + 5) = shellR3i
      aVector(j + 6) = shellR2i * shellR2i
      j = j + 7
    end do
  end subroutine evec

  subroutine sqrdc(x, ldx, n, p, qraux, jpvt, work, job)
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ldx
      integer  :: n
      integer , intent(in) :: p
      integer , intent(in) :: job
      integer , intent(inout) :: jpvt(ldx)
      real  :: x(ldx,ldx)
      real , intent(out) :: qraux(ldx)
      real , intent(inout) :: work(ldx)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, jp, l, lp1, lup, maxj, pl, pu, jj
      real :: maxnrm, tt, nrmxl, t
      logical :: negj, swapj
      real , external :: sdot, snrm2
!-----------------------------------------------
!
!     SQRDC USES HOUSEHOLDER TRANSFORMATIONS TO COMPUTE THE QR
!     FACTORIZATION OF AN N BY P MATRIX X.  COLUMN PIVOTING
!     BASED ON THE 2-NORMS OF THE REDUCED COLUMNS MAY BE
!     PERFORMED AT THE USERS OPTION.
!
!     ON ENTRY
!
!        X       REAL(LDX,P), WHERE LDX .GE. N.
!                X CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO BE
!                COMPUTED.
!
!        LDX     INTEGER.
!                LDX IS THE LEADING DIMENSION OF THE ARRAY X.
!
!        N       INTEGER.
!                N IS THE NUMBER OF ROWS OF THE MATRIX X.
!
!        P       INTEGER.
!                P IS THE NUMBER OF COLUMNS OF THE MATRIX X.
!
!        JPVT    INTEGER(P).
!                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION
!                OF THE PIVOT COLUMNS.  THE K-TH COLUMN X(K) OF X
!                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE
!                VALUE OF JPVT(K).
!
!                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL
!                                      COLUMN.
!
!                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE COLUMN.
!
!                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL COLUMN.
!
!                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL COLUMNS
!                ARE MOVED TO THE BEGINNING OF THE ARRAY X AND FINAL
!                COLUMNS TO THE END.  BOTH INITIAL AND FINAL COLUMNS
!                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY
!                FREE COLUMNS ARE MOVED.  AT THE K-TH STAGE OF THE
!                REDUCTION, IF X(K) IS OCCUPIED BY A FREE COLUMN
!                IT IS INTERCHANGED WITH THE FREE COLUMN OF LARGEST
!                REDUCED NORM.  JPVT IS NOT REFERENCED IF
!                JOB .EQ. 0.
!
!        WORK    REAL(P).
!                WORK IS A WORK ARRAY.  WORK IS NOT REFERENCED IF
!                JOB .EQ. 0.
!
!        JOB     INTEGER.
!                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING.
!                IF JOB .EQ. 0, NO PIVOTING IS DONE.
!                IF JOB .NE. 0, PIVOTING IS DONE.
!
!     ON RETURN
!
!        X       X CONTAINS IN ITS UPPER TRIANGLE THE UPPER
!                TRIANGULAR MATRIX R OF THE QR FACTORIZATION.
!                BELOW ITS DIAGONAL X CONTAINS INFORMATION FROM
!                WHICH THE ORTHOGONAL PART OF THE DECOMPOSITION
!                CAN BE RECOVERED.  NOTE THAT IF PIVOTING HAS
!                BEEN REQUESTED, THE DECOMPOSITION IS NOT THAT
!                OF THE ORIGINAL MATRIX X BUT THAT OF X
!                WITH ITS COLUMNS PERMUTED AS DESCRIBED BY JPVT.
!
!        QRAUX   REAL(P).
!                QRAUX CONTAINS FURTHER INFORMATION REQUIRED TO RECOVER
!                THE ORTHOGONAL PART OF THE DECOMPOSITION.
!
!        JPVT    JPVT(K) CONTAINS THE INDEX OF THE COLUMN OF THE
!                ORIGINAL MATRIX THAT HAS BEEN INTERCHANGED INTO
!                THE K-TH COLUMN, IF PIVOTING WAS REQUESTED.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
!
!     SQRDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
!
!     BLAS SAXPY,SDOT,SSCAL,SSWAP,SNRM2
!     FORTRAN ABS,AMAX1,MIN0,SQRT
!
!     INTERNAL VARIABLES
!
!
!
      pl = 1
      pu = 0
      if (job /= 0) then
!
!        PIVOTING HAS BEEN REQUESTED.  REARRANGE THE COLUMNS
!        ACCORDING TO JPVT.
!
        do j = 1, p
          swapj = jpvt(j) > 0
          negj = jpvt(j) < 0
          jpvt(j) = j
          if (negj) jpvt(j) = -j
          if (.not.swapj) cycle
          if (j /= pl) call sswap (n, x(1,pl), 1, x(1,j), 1)
          jpvt(j) = jpvt(pl)
          jpvt(pl) = j
          pl = pl + 1
        end do
        pu = p
        do jj = 1, p
          j = p - jj + 1
          if (jpvt(j) >= 0) cycle
          jpvt(j) = -jpvt(j)
          if (j /= pu) then
            call sswap (n, x(1,pu), 1, x(1,j), 1)
            jp = jpvt(pu)
            jpvt(pu) = jpvt(j)
            jpvt(j) = jp
          end if
          pu = pu - 1
        end do
      end if
      if (pu >= pl) then
        do j = pl, pu
          qraux(j) = snrm2(n,x(1,j),1)
          work(j) = qraux(j)
        end do
      end if
      lup = min0(n,p)
      do l = 1, lup
        if (l>=pl .and. l<pu) then
!
!           LOCATE THE COLUMN OF LARGEST NORM AND BRING IT
!           INTO THE PIVOT POSITION.
!
          maxnrm = 0.0E0
          maxj = l
          do j = l, pu
            if (qraux(j) <= maxnrm) cycle
            maxnrm = qraux(j)
            maxj = j
          end do
          if (maxj /= l) then
            call sswap (n, x(1,l), 1, x(1,maxj), 1)
            qraux(maxj) = qraux(l)
            work(maxj) = work(l)
            jp = jpvt(maxj)
            jpvt(maxj) = jpvt(l)
            jpvt(l) = jp
          end if
        end if
        qraux(l) = 0.0E0
        if (l == n) cycle
!
!           COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L.
!
        nrmxl = snrm2(n - l + 1,x(l,l),1)
        if (nrmxl == 0.0E0) cycle
        if (x(l,l) /= 0.0E0) nrmxl = sign(nrmxl,x(l,l))
        call sscal (n - l + 1, 1.0E0/nrmxl, x(l,l), 1)
        x(l,l) = 1.0E0 + x(l,l)
!
!              APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS,
!              UPDATING THE NORMS.
!
        lp1 = l + 1
        if (p >= lp1) then
          do j = lp1, p
            t = -sdot(n - l + 1,x(l,l),1,x(l,j),1)/x(l,l)
            call saxpy (n - l + 1, t, x(l,l), 1, x(l,j), 1)
            if (j<pl .or. j>pu) cycle
            if (qraux(j) == 0.0E0) cycle
            tt = 1.0E0 - (abs(x(l,j))/qraux(j))**2
            tt = amax1(tt,0.0E0)
            t = tt
            tt = 1.0E0 + 0.05E0*tt*(qraux(j)/work(j))**2
            if (Abs(tt - 1.0E0) > 1.e-10 ) then
              qraux(j) = qraux(j)*sqrt(t)
            else
              qraux(j) = snrm2(n - l,x(l+1,j),1)
              work(j) = qraux(j)
            end if
          end do
        end if
!
!              SAVE THE TRANSFORMATION.
!
        qraux(l) = x(l,l)
        x(l,l) = -nrmxl
      end do
      return
      end subroutine sqrdc


      subroutine saxpy(n, sa, sx, incx, sy, incy)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      integer , intent(in) :: incx
      integer , intent(in) :: incy
      real , intent(in) :: sa
      real , intent(in) :: sx(1)
      real , intent(inout) :: sy(1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ix, iy, m, mp1
!-----------------------------------------------
!
!     CONSTANT TIMES A VECTOR PLUS A VECTOR.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
!
      if (n <= 0) return
      if (sa == 0.0E0) return
      if (incx/=1 .or. incy/=1) then
!
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1
!
        ix = 1
        iy = 1
        if (incx < 0) ix = ((-n) + 1)*incx + 1
        if (incy < 0) iy = ((-n) + 1)*incy + 1
        sy(iy:(n-1)*incy+iy:incy) = sy(iy:(n-1)*incy+iy:incy) + sa*sx(ix:(n-1)*&
          incx+ix:incx)
        return
      end if
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
      m = mod(n,4)
      if (m /= 0) then
        sy(:m) = sy(:m) + sa*sx(:m)
        if (n < 4) return
      end if
      mp1 = m + 1
      sy(mp1:((n-mp1+4)/4)*4-1+mp1) = sy(mp1:((n-mp1+4)/4)*4-1+mp1) + sa*sx(mp1&
        :((n-mp1+4)/4)*4-1+mp1)
      return
      end subroutine saxpy


      real function sdot (n, sx, incx, sy, incy)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      integer , intent(in) :: incx
      integer , intent(in) :: incy
      real , intent(in) :: sx(1)
      real , intent(in) :: sy(1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ix, iy, m, mp1
      real :: temp
!-----------------------------------------------
!
!     FORMS THE DOT PRODUCT OF TWO VECTORS.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
!
      sdot = 0.0E0
      temp = 0.0E0
      if (n <= 0) return
      if (incx/=1 .or. incy/=1) then
!
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1
!
        ix = 1
        iy = 1
        if (incx < 0) ix = ((-n) + 1)*incx + 1
        if (incy < 0) iy = ((-n) + 1)*incy + 1
        temp = dot_product(sx(ix:(n-1)*incx+ix:incx),sy(iy:(n-1)*incy+iy:incy))
        sdot = temp
        return
      end if
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
      m = mod(n,5)
      if (m == 0) go to 40
      temp = dot_product(sx(:m),sy(:m))
      if (n < 5) go to 60
   40 continue
      mp1 = m + 1
      temp = temp + dot_product(sx(mp1:((n-mp1+5)/5)*5-1+mp1),sy(mp1:((n-mp1+5)&
        /5)*5-1+mp1))
   60 continue
      sdot = temp
      return
      end function sdot


      subroutine sscal(n, da, dx, incx)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      integer , intent(in) :: incx
      real , intent(in) :: da
      real , intent(inout) :: dx(1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: m, mp1, nincx
!-----------------------------------------------
!
!     SCALES A VECTOR BY A CONSTANT.
!     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
!
      if (n <= 0) return
      if (incx /= 1) then
!
!        CODE FOR INCREMENT NOT EQUAL TO 1
!
        nincx = n*incx
        dx(:nincx:incx) = da*dx(:nincx:incx)
        return
      end if
!
!        CODE FOR INCREMENT EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
      m = mod(n,5)
      if (m /= 0) then
        dx(:m) = da*dx(:m)
        if (n < 5) return
      end if
      mp1 = m + 1
      dx(mp1:((n-mp1+5)/5)*5-1+mp1) = da*dx(mp1:((n-mp1+5)/5)*5-1+mp1)
      return
      end subroutine sscal


      real function snrm2 (n, sx, incx)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      integer , intent(in) :: incx
      real , intent(in) :: sx(1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: next, nn, i, j
      real :: cutlo, cuthi, hitest, sum, xmax, zero, one
!-----------------------------------------------
      data zero, one/ 0.0E0, 1.0E0/
!
!     EUCLIDEAN NORM OF THE N-VECTOR STORED IN SX() WITH STORAGE
!     INCREMENT INCX .
!     IF    N .LE. 0 RETURN WITH RESULT = 0.
!     IF N .GE. 1 THEN INCX MUST BE .GE. 1
!
!           C.L.LAWSON, 1978 JAN 08
!
!     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
!     HOPEFULLY APPLICABLE TO ALL MACHINES.
!         CUTLO = MAXIMUM OF  SQRT(U/EPS)  OVER ALL KNOWN MACHINES.
!         CUTHI = MINIMUM OF  SQRT(V)      OVER ALL KNOWN MACHINES.
!     WHERE
!         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
!         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
!         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
!
!     BRIEF OUTLINE OF ALGORITHM..
!
!     PHASE 1    SCANS ZERO COMPONENTS.
!     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
!     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
!     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
!     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
!
!     VALUES FOR CUTLO AND CUTHI..
!     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
!     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
!     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
!                   UNIVAC AND DEC AT 2**(-103)
!                   THUS CUTLO = 2**(-51) = 4.44089E-16
!     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
!                   THUS CUTHI = 2**(63.5) = 1.30438E19
!     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
!                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
!     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
!     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
!     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      data cutlo, cuthi/ 4.441E-16, 1.304E19/
      xmax = 0.0
!
      if (n <= 0) then
        snrm2 = zero
      else
!
        next = 30
        sum = zero
        nn = n*incx
!                                                 BEGIN MAIN LOOP
        i = 1
   20   continue
        if (next == 30) go to 30
        if (next == 50) go to 50
        if (next == 70) go to 70
        if (next == 110) go to 110
   30   continue
        if (abs(sx(i)) > cutlo) go to 85
        next = 50
        xmax = zero
!
!                        PHASE 1.  SUM IS ZERO
!
   50   continue
        if (Abs(sx(i)) < 1.d-20) go to 200
        if (abs(sx(i)) > cutlo) go to 85
!
!                                PREPARE FOR PHASE 2.
        next = 70
        go to 105
!
!                                PREPARE FOR PHASE 4.
!
  100   continue
        i = j
        next = 110
        sum = (sum/sx(i))/sx(i)
  105   continue
        xmax = abs(sx(i))
        go to 115
!
!                   PHASE 2.  SUM IS SMALL.
!                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
!
   70   continue
        if (abs(sx(i)) > cutlo) go to 75
!
!                     COMMON CODE FOR PHASES 2 AND 4.
!                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
!
  110   continue
        if (abs(sx(i)) <= xmax) go to 115
        sum = one + sum*(xmax/sx(i))**2
        xmax = abs(sx(i))
        go to 200
!
  115   continue
        sum = sum + (sx(i)/xmax)**2
        go to 200
!
!
!                  PREPARE FOR PHASE 3.
!
   75   continue
        sum = (sum*xmax)*xmax
!
!
!     FOR REAL OR D.P. SET HITEST = CUTHI/N
!     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
!
   85   continue
        hitest = cuthi/float(n)
!
!                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
!
        do j = i, nn, incx
          if (abs(sx(j)) >= hitest) go to 100
          sum = sum + sx(j)**2
        end do
        snrm2 = sqrt(sum)
        go to 300
!
  200   continue
        i = i + incx
        if (i <= nn) go to 20
!
!              END OF MAIN LOOP.
!
!              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
!
        snrm2 = xmax*sqrt(sum)
      end if
  300 continue
      return
      end function snrm2


      subroutine sswap(n, sx, incx, sy, incy)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      integer , intent(in) :: incx
      integer , intent(in) :: incy
      real , intent(inout) :: sx(1)
      real , intent(inout) :: sy(1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ix, iy, m, mp1
      real :: stemp
!-----------------------------------------------
!
!     INTERCHANGES TWO VECTORS.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO 1.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
!
      if (n <= 0) return
      if (incx/=1 .or. incy/=1) then
!
!       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
!         TO 1
!
        ix = 1
        iy = 1
        if (incx < 0) ix = ((-n) + 1)*incx + 1
        if (incy < 0) iy = ((-n) + 1)*incy + 1
        do i = 1, n
          stemp = sx(ix)
          sx(ix) = sy(iy)
          sy(iy) = stemp
          ix = ix + incx
          iy = iy + incy
        end do
        return
      end if
!
!       CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!       CLEAN-UP LOOP
!
      m = mod(n,3)
      if (m /= 0) then
        do i = 1, m
          stemp = sx(i)
          sx(i) = sy(i)
          sy(i) = stemp
        end do
        if (n < 3) return
      end if
      mp1 = m + 1
      do i = mp1, n, 3
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
        stemp = sx(i+1)
        sx(i+1) = sy(i+1)
        sy(i+1) = stemp
        stemp = sx(i+2)
        sx(i+2) = sy(i+2)
        sy(i+2) = stemp
      end do
      return
      end subroutine sswap
            MODULE scopy_I
      INTERFACE
      SUBROUTINE scopy (N, SX, INCX, SY, INCY)
      integer, INTENT(IN) :: N
      real, DIMENSION(1), INTENT(IN) :: SX
      integer, INTENT(IN) :: INCX
      real, DIMENSION(1), INTENT(OUT) :: SY
      integer, INTENT(IN) :: INCY
      END SUBROUTINE
      END INTERFACE
      END MODULE
      MODULE saxpy_I
      INTERFACE
      SUBROUTINE saxpy (N, SA, SX, INCX, SY, INCY)
      integer, INTENT(IN) :: N
      real, INTENT(IN) :: SA
      real, DIMENSION(1), INTENT(IN) :: SX
      integer, INTENT(IN) :: INCX
      real, DIMENSION(1), INTENT(INOUT) :: SY
      integer, INTENT(IN) :: INCY
      END SUBROUTINE
      END INTERFACE
      END MODULE

      subroutine sqrsl(x, ldx, n, k, qraux, y, qy, qty, b, rsd, xb, job, info)
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ldx
      integer  :: n
      integer  :: k
      integer , intent(in) :: job
      integer , intent(out) :: info
      real  :: x(ldx,ldx)
      real , intent(in) :: qraux(ldx)
      real  :: y(ldx)
      real  :: qy(ldx)
      real  :: qty(ldx)
      real  :: b(ldx)
      real  :: rsd(ldx)
      real  :: xb(ldx)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, jj, ju, kp1
      real :: t, temp
      logical :: cb, cqy, cqty, cr, cxb
      real , external :: sdot
!-----------------------------------------------
!
!     SQRSL APPLIES THE OUTPUT OF SQRDC TO COMPUTE COORDINATE
!     TRANSFORMATIONS, PROJECTIONS, AND LEAST SQUARES SOLUTIONS.
!     FOR K .LE. MIN(N,P), LET XK BE THE MATRIX
!
!            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K)))
!
!     FORMED FROM COLUMNNS JPVT(1), ... ,JPVT(K) OF THE ORIGINAL
!     N X P MATRIX X THAT WAS INPUT TO SQRDC (IF NO PIVOTING WAS
!     DONE, XK CONSISTS OF THE FIRST K COLUMNS OF X IN THEIR
!     ORIGINAL ORDER).  SQRDC PRODUCES A FACTORED ORTHOGONAL MATRIX Q
!     AND AN UPPER TRIANGULAR MATRIX R SUCH THAT
!
!              XK = Q * (R)
!                       (0)
!
!     THIS INFORMATION IS CONTAINED IN CODED FORM IN THE ARRAYS
!     X AND QRAUX.
!
!     ON ENTRY
!
!        X      REAL(LDX,P).
!               X CONTAINS THE OUTPUT OF SQRDC.
!
!        LDX    INTEGER.
!               LDX IS THE LEADING DIMENSION OF THE ARRAY X.
!
!        N      INTEGER.
!               N IS THE NUMBER OF ROWS OF THE MATRIX XK.  IT MUST
!               HAVE THE SAME VALUE AS N IN SQRDC.
!
!        K      INTEGER.
!               K IS THE NUMBER OF COLUMNS OF THE MATRIX XK.  K
!               MUST NNOT BE GREATER THAN MIN(N,P), WHERE P IS THE
!               SAME AS IN THE CALLING SEQUENCE TO SQRDC.
!
!        QRAUX  REAL(P).
!               QRAUX CONTAINS THE AUXILIARY OUTPUT FROM SQRDC.
!
!        Y      REAL(N)
!               Y CONTAINS AN N-VECTOR THAT IS TO BE MANIPULATED
!               BY SQRSL.
!
!        JOB    INTEGER.
!               JOB SPECIFIES WHAT IS TO BE COMPUTED.  JOB HAS
!               THE DECIMAL EXPANSION ABCDE, WITH THE FOLLOWING
!               MEANING.
!
!                    IF A.NE.0, COMPUTE QY.
!                    IF B,C,D, OR E .NE. 0, COMPUTE QTY.
!                    IF C.NE.0, COMPUTE B.
!                    IF D.NE.0, COMPUTE RSD.
!                    IF E.NE.0, COMPUTE XB.
!
!               NOTE THAT A REQUEST TO COMPUTE B, RSD, OR XB
!               AUTOMATICALLY TRIGGERS THE COMPUTATION OF QTY, FOR
!               WHICH AN ARRAY MUST BE PROVIDED IN THE CALLING
!               SEQUENCE.
!
!     ON RETURN
!
!        QY     REAL(N).
!               QY CONNTAINS Q*Y, IF ITS COMPUTATION HAS BEEN
!               REQUESTED.
!
!        QTY    REAL(N).
!               QTY CONTAINS TRANS(Q)*Y, IF ITS COMPUTATION HAS
!               BEEN REQUESTED.  HERE TRANS(Q) IS THE
!               TRANSPOSE OF THE MATRIX Q.
!
!        B      REAL(K)
!               B CONTAINS THE SOLUTION OF THE LEAST SQUARES PROBLEM
!
!                    MINIMIZE NORM2(Y - XK*B),
!
!               IF ITS COMPUTATION HAS BEEN REQUESTED.  (NOTE THAT
!               IF PIVOTING WAS REQUESTED IN SQRDC, THE J-TH
!               COMPONENT OF B WILL BE ASSOCIATED WITH COLUMN JPVT(J)
!               OF THE ORIGINAL MATRIX X THAT WAS INPUT INTO SQRDC.)
!
!        RSD    REAL(N).
!               RSD CONTAINS THE LEAST SQUARES RESIDUAL Y - XK*B,
!               IF ITS COMPUTATION HAS BEEN REQUESTED.  RSD IS
!               ALSO THE ORTHOGONAL PROJECTION OF Y ONTO THE
!               ORTHOGONAL COMPLEMENT OF THE COLUMN SPACE OF XK.
!
!        XB     REAL(N).
!               XB CONTAINS THE LEAST SQUARES APPROXIMATION XK*B,
!               IF ITS COMPUTATION HAS BEEN REQUESTED.  XB IS ALSO
!               THE ORTHOGONAL PROJECTION OF Y ONTO THE COLUMN SPACE
!               OF X.
!
!        INFO   INTEGER.
!               INFO IS ZERO UNLESS THE COMPUTATION OF B HAS
!               BEEN REQUESTED AND R IS EXACTLY SINGULAR.  IN
!               THIS CASE, INFO IS THE INDEX OF THE FIRST ZERO
!               DIAGONAL ELEMENT OF R AND B IS LEFT UNALTERED.
!
!     THE PARAMETERS QY, QTY, B, RSD, AND XB ARE NOT REFERENCED
!     IF THEIR COMPUTATION IS NOT REQUESTED AND IN THIS CASE
!     CAN BE REPLACED BY DUMMY VARIABLES IN THE CALLING PROGRAM.
!     TO SAVE STORAGE, THE USER MAY IN SOME CASES USE THE SAME
!     ARRAY FOR DIFFERENT PARAMETERS IN THE CALLING SEQUENCE.  A
!     FREQUENTLY OCCURING EXAMPLE IS WHEN ONE WISHES TO COMPUTE
!     ANY OF B, RSD, OR XB AND DOES NOT NEED Y OR QTY.  IN THIS
!     CASE ONE MAY IDENTIFY Y, QTY, AND ONE OF B, RSD, OR XB, WHILE
!     PROVIDING SEPARATE ARRAYS FOR ANYTHING ELSE THAT IS TO BE
!     COMPUTED.  THUS THE CALLING SEQUENCE
!
!          CALL SQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO)
!
!     WILL RESULT IN THE COMPUTATION OF B AND RSD, WITH RSD
!     OVERWRITING Y.  MORE GENERALLY, EACH ITEM IN THE FOLLOWING
!     LIST CONTAINS GROUPS OF PERMISSIBLE IDENTIFICATIONS FOR
!     A SINGLE CALLINNG SEQUENCE.
!
!          1. (Y,QTY,B) (RSD) (XB) (QY)
!
!          2. (Y,QTY,RSD) (B) (XB) (QY)
!
!          3. (Y,QTY,XB) (B) (RSD) (QY)
!
!          4. (Y,QY) (QTY,B) (RSD) (XB)
!
!          5. (Y,QY) (QTY,RSD) (B) (XB)
!
!          6. (Y,QY) (QTY,XB) (B) (RSD)
!
!     IN ANY GROUP THE VALUE RETURNED IN THE ARRAY ALLOCATED TO
!     THE GROUP CORRESPONDS TO THE LAST MEMBER OF THE GROUP.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
!
!     SQRSL USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
!
!     BLAS SAXPY,SCOPY,SDOT
!     FORTRAN ABS,MIN0,MOD
!
!     INTERNAL VARIABLES
!
!
!
!     SET INFO FLAG.
!
      info = 0
!
!     DETERMINE WHAT IS TO BE COMPUTED.
!
      cqy = job/10000 /= 0
      cqty = mod(job,10000) /= 0
      cb = mod(job,1000)/100 /= 0
      cr = mod(job,100)/10 /= 0
      cxb = mod(job,10) /= 0
      ju = min0(k,n - 1)
!
!     SPECIAL ACTION WHEN N=1.
!
      if (ju == 0) then
        if (cqy) qy(1) = y(1)
        if (cqty) qty(1) = y(1)
        if (cxb) xb(1) = y(1)
        if (cb) then
          if (x(1,1) == 0.0E0) then
            info = 1
          else
            b(1) = y(1)/x(1,1)
          end if
        end if
        if (cr) rsd(1) = 0.0E0
      else
        if (cqy) call scopy (n, y, 1, qy, 1)
        if (cqty) call scopy (n, y, 1, qty, 1)
        if (cqy) then
!
!           COMPUTE QY.
!
          do jj = 1, ju
            j = ju - jj + 1
            if (qraux(j) == 0.0E0) cycle
            temp = x(j,j)
            x(j,j) = qraux(j)
            t = -sdot(n - j + 1,x(j,j),1,qy(j),1)/x(j,j)
            call saxpy (n - j + 1, t, x(j,j), 1, qy(j), 1)
            x(j,j) = temp
          end do
        end if
        if (cqty) then
!
!           COMPUTE TRANS(Q)*Y.
!
          do j = 1, ju
            if (qraux(j) == 0.0E0) cycle
            temp = x(j,j)
            x(j,j) = qraux(j)
            t = -sdot(n - j + 1,x(j,j),1,qty(j),1)/x(j,j)
            call saxpy (n - j + 1, t, x(j,j), 1, qty(j), 1)
            x(j,j) = temp
          end do
        end if
        if (cb) call scopy (k, qty, 1, b, 1)
        kp1 = k + 1
        if (cxb) call scopy (k, qty, 1, xb, 1)
        if (cr .and. k<n) call scopy (n - k, qty(kp1), 1, rsd(kp1), 1)
        if (.not.(.not.cxb .or. kp1>n)) then
          xb(kp1:n) = 0.0E0
        end if
        if (cr) then
          rsd(:k) = 0.0E0
        end if
        if (cb) then
!
!           COMPUTE B.
!
          do jj = 1, k
            j = k - jj + 1
            if (x(j,j) == 0.0E0) then
              info = j
!           ......EXIT
              exit
            end if
            b(j) = b(j)/x(j,j)
            if (j == 1) cycle
            t = -b(j)
            call saxpy (j - 1, t, x(1,j), 1, b, 1)
          end do
        end if
        if (.not.(.not.cr .and. .not.cxb)) then
!
!           COMPUTE RSD OR XB AS REQUIRED.
!
          do jj = 1, ju
            j = ju - jj + 1
            if (qraux(j) == 0.0E0) cycle
            temp = x(j,j)
            x(j,j) = qraux(j)
            if (cr) then
              t = -sdot(n - j + 1,x(j,j),1,rsd(j),1)/x(j,j)
              call saxpy (n - j + 1, t, x(j,j), 1, rsd(j), 1)
            end if
            if (cxb) then
              t = -sdot(n - j + 1,x(j,j),1,xb(j),1)/x(j,j)
              call saxpy (n - j + 1, t, x(j,j), 1, xb(j), 1)
            end if
            x(j,j) = temp
          end do
        end if
      end if
      return
      end subroutine sqrsl


      subroutine scopy(n, sx, incx, sy, incy)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      integer , intent(in) :: incx
      integer , intent(in) :: incy
      real , intent(in) :: sx(1)
      real , intent(out) :: sy(1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ix, iy, m, mp1
!-----------------------------------------------
!
!     COPIES A VECTOR, X, TO A VECTOR, Y.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO 1.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
!
      if (n <= 0) return
      if (incx/=1 .or. incy/=1) then
!
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1
!
        ix = 1
        iy = 1
        if (incx < 0) ix = ((-n) + 1)*incx + 1
        if (incy < 0) iy = ((-n) + 1)*incy + 1
        sy(iy:(n-1)*incy+iy:incy) = sx(ix:(n-1)*incx+ix:incx)
        return
      end if
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
      m = mod(n,7)
      if (m /= 0) then
        sy(:m) = sx(:m)
        if (n < 7) return
      end if
      mp1 = m + 1
      sy(mp1:((n-mp1+7)/7)*7-1+mp1) = sx(mp1:((n-mp1+7)/7)*7-1+mp1)
      return
      end subroutine scopy









