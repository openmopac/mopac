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

  subroutine prtpka(ipKa_sorted, pKa_sorted, ipKa_unsorted, pKa_unsorted,  no)
    use common_arrays_C, only : grad, xparam, p, nat, coord
    use molkst_C, only : numat, escf, mozyme, emin, mpack, moperr
    use linear_cosmo, only : ini_linear_cosmo, coscavz
    use parameters_C, only : tore
    use cosmo_C, only : iseps, useps, lpka
    use chanel_C, only : iw
    implicit none
    integer, intent (out) :: ipKa_sorted(numat), ipKa_unsorted(numat), no
    double precision, intent (out) :: pKa_sorted(numat), pKa_unsorted(numat)
!
! Calculates and prints the pKa_sorted for an organic molecule with an ionizable hydrogen
! attached to an oxygen atom.
!
   double precision :: sum, c1, c2, c3, dist(numat), sum_min, sum_H_min
   integer :: i, j, k, nh, loop, i_min
   double precision, dimension(numat) :: q
   logical :: store_iseps, store_useps
!
!  First, switch in parameters specific to the pKa_sorted calculation.
!
    call Parameters_for_PKA(c1, c2, c3)
!
!  Re-calculate the charges
!
    store_iseps = iseps
    store_useps = useps
    iseps = .true.
    useps = .true.
    emin = 0.d0
    call cosini(.true.)
    if (mozyme) call ini_linear_cosmo
    if (mozyme) then
      call coscavz(coord, nat)
      lpka = .true.
    else
       call coscav
       call mkbmat
    end if
    i = mpack
    call moldat(1)
    mpack = i
    moperr = .false. ! An error in moldat is not important here. only the electronics are affected.
    call calpar
    call compfg (xparam, .TRUE., escf, .TRUE., grad, .FALSE.)
    nh = 0
    if (moperr) return
    lpka = .false.
    call chrge (p, q)
    q(:numat) = tore(nat(:numat)) - q(:numat)
!
!  Find the ionizable hydrogens
!
    pKa_unsorted = 0.d0
    do i = 1, numat
      if (nat(i) /= 1) cycle
      sum_min = 2.d0
      i_min = 0
      do j = 1, numat
        if (nat(j) /= 8) cycle
        sum = (coord(1,i) - coord(1,j))**2 + &
            & (coord(2,i) - coord(2,j))**2 + &
            & (coord(3,i) - coord(3,j))**2
        if (sum < sum_min) then      ! Assume that a hydrogen within 1.4 Angstroms of an oxygen is
          sum_min = sum              ! attached to that atom
          i_min = j
        end if
      end do
      if (sum_min < 1.9999d0) then
!
! Exclude water
!
        sum_H_min = 2.d0
        do k = 1, numat
          if (nat(k) /= 1 .or. k == i) cycle
          sum = (coord(1,k) - coord(1,i_min))**2 + &
              & (coord(2,k) - coord(2,i_min))**2 + &
              & (coord(3,k) - coord(3,i_min))**2
! Assume that a hydrogen within 1.4 Angstroms of an oxygen is attached to that atom
          if (sum < sum_H_min) then
            sum_H_min = sum
          end if
        end do
        if (sum_H_min < 1.9999d0 .and. numat > 3) cycle
        nh = nh + 1
        ipKa_unsorted(nh) = i
        pKa_unsorted(nh)=q(i)
        dist(nh) = sqrt(sum_min)
      end if
    end do
!
! Now sort into decending order
!
    ipKa_sorted = 0
    pKa_sorted = 0.d0
    if (nh == 0) then
      write(iw,"(/,3(10x,a,/))")"A request was made to print the pKa_sorted values for this system,", &
      "but there are no hydrogen atoms attached to an oxygen atom,", &
      "so the pKa_sorted calculation cannot be completed."
      return
    end if
    do i = 1, nh !  Convert charges and bond-lengths into pKa_sorted values
     pKa_unsorted(i) = pKa_unsorted(i)*c1 + dist(i)*c2 + c3
    end do
    loop = 0
    do i = 1, nh  ! Sort pKa_sorted into order, select lowest four
      sum = 100.d0
      k = 0
      do j = 1, nh
        if (pKa_unsorted(j) < sum) then
          k = j
          sum = pKa_unsorted(j)
        end if
      end do
      if (k == 0) exit
      j = 1
      if (loop > 0) then
        do j = 1, loop
          if (ipKa_sorted(j) == ipKa_unsorted(k) .and. Abs(pKa_sorted(j) - pKa_unsorted(k)) < 0.2d0) exit
        end do
      end if
      if (j <= loop) cycle
      loop = loop + 1
      ipKa_sorted(loop) = ipKa_unsorted(k)
      pKa_sorted(loop) =  pKa_unsorted(k)
      pKa_unsorted(k) = 200.d0 + pKa_unsorted(k)
    end do
    pKa_unsorted(:loop) = pKa_unsorted(:loop) - 200.d0
    no = loop
    call switch
    iseps = store_iseps
    useps = store_useps
    i = mpack
    call moldat(1)
    mpack = i
    moperr = .false. ! An error in moldat is not important here. only the electronics are affected.
    call calpar
    emin = 0.d0
    call compfg (xparam, .TRUE., escf, .TRUE., grad, .FALSE.)
    return
  end subroutine prtpka
  subroutine Parameters_for_PKA(c1, c2, c3)
    use parameters_C, only : uss, upp, udd, zs, zp, zd, betas, &
    betap, betad, gss, gsp, gpp, gp2, hsp, zsn, zpn, zdn
    double precision, intent (out) :: c1, c2, c3
    c1   =    -288.0530769d0
    c2   =      28.6888717d0
    c3   =      89.1172382d0
!
!                    Data for pKa_sorted for Element  1         Hydrogen
!
    uss  ( 1)   =      -9.3555403d0
    betas( 1)   =      -2.6789813d0
    zs   ( 1)   =       1.2458568d0
    gss  ( 1)   =      14.5964443d0
!
!                    Data for pKa_sorted for Element  6           Carbon
!
    uss  ( 6)   =     -49.9390029d0
    upp  ( 6)   =     -43.8233347d0
    betas( 6)   =     -12.6446141d0
    betap( 6)   =      -9.4614222d0
    zs   ( 6)   =       1.6412075d0
    zp   ( 6)   =       1.5810984d0
    gss  ( 6)   =      16.7616538d0
    gsp  ( 6)   =      12.3925395d0
    gpp  ( 6)   =      11.2679258d0
    gp2  ( 6)   =      11.0158335d0
    hsp  ( 6)   =       0.4200022d0
!
!                    Data for pKa_sorted for Element  7         Nitrogen
!
    uss  ( 7)   =     -55.4307466d0
    upp  ( 7)   =     -49.8247298d0
    betas( 7)   =     -21.7612447d0
    betap( 7)   =     -16.9335869d0
    zs   ( 7)   =       1.4860676d0
    zp   ( 7)   =       2.0420069d0
    gss  ( 7)   =       7.1149144d0
    gsp  ( 7)   =       7.2242656d0
    gpp  ( 7)   =      14.8103679d0
    gp2  ( 7)   =      11.4500819d0
    hsp  ( 7)   =       4.5873409d0
!
!                    Data for pKa_sorted for Element  8           Oxygen
!
    uss  ( 8)   =     -89.9472647d0
    upp  ( 8)   =     -70.7249961d0
    betas( 8)   =     -66.7497177d0
    betap( 8)   =     -21.7795050d0
    zs   ( 8)   =       4.3991035d0
    zp   ( 8)   =       2.1614926d0
    gss  ( 8)   =      16.3875625d0
    gsp  ( 8)   =      16.0293626d0
    gpp  ( 8)   =      16.6098992d0
    gp2  ( 8)   =      11.0339705d0
    hsp  ( 8)   =       4.7920603d0
!
!                    Data for pKa_sorted for Element  9         Fluorine
!
    uss  ( 9)   =    -141.5218991d0
    upp  ( 9)   =     -98.4596662d0
    betas( 9)   =     -69.7628706d0
    betap( 9)   =     -30.1591819d0
    zs   ( 9)   =       5.9393528d0
    zp   ( 9)   =       4.7675390d0
    gss  ( 9)   =      12.5528331d0
    gsp  ( 9)   =      19.8940208d0
    gpp  ( 9)   =       8.8875462d0
    gp2  ( 9)   =      12.5543486d0
    hsp  ( 9)   =       3.5419742d0
!
!                    Data for pKa_sorted for Element 17         Chlorine
!
    uss  (17)   =     -61.5736619d0
    upp  (17)   =     -54.5254717d0
    udd  (17)   =     -38.2581550d0
    betas(17)   =      -0.5351767d0
    betap(17)   =     -11.7700638d0
    betad(17)   =      -4.0377510d0
    zs   (17)   =       4.2327248d0
    zp   (17)   =       1.2135427d0
    zd   (17)   =       1.3240330d0
    zsn  (17)   =       0.9562970d0
    zpn  (17)   =       2.4640670d0
    zdn  (17)   =       6.4103250d0
    gss  (17)   =       9.4258401d0
    gsp  (17)   =       5.3649222d0
    gpp  (17)   =      10.2118976d0
    gp2  (17)   =       9.1156399d0
    hsp  (17)   =       4.9183864d0
!
!                    Data for pKa_sorted for Element 35          Bromine
!
    uss  (35)   =     -46.3749516d0
    upp  (35)   =     -50.1929356d0
    udd  (35)   =       7.0867380d0
    betas(35)   =     -31.5606655d0
    betap(35)   =      -9.1283278d0
    betad(35)   =      -9.8391240d0
    zs   (35)   =       5.0204513d0
    zp   (35)   =       2.2579169d0
    zd   (35)   =       1.5210310d0
    zsn  (35)   =       3.0947770d0
    zpn  (35)   =       3.0657640d0
    zdn  (35)   =       2.8200030d0
    gss  (35)   =       8.1108378d0
    gsp  (35)   =       5.2179143d0
    gpp  (35)   =      10.3707055d0
    gp2  (35)   =       8.5089799d0
    hsp  (35)   =       4.8307027d0
!
!                    Data for pKa_sorted for Element 53           Iodine
!
    uss  (53)   =     -59.6492236d0
    upp  (53)   =     -56.2042321d0
    udd  (53)   =     -28.8226030d0
    betas(53)   =     -30.0143361d0
    betap(53)   =      -5.5364107d0
    betad(53)   =      -7.6761070d0
    zs   (53)   =       4.9453882d0
    zp   (53)   =       2.4202135d0
    zd   (53)   =       1.8751750d0
    zsn  (53)   =       9.1352440d0
    zpn  (53)   =       6.8881910d0
    zdn  (53)   =       3.7915230d0
    gss  (53)   =       7.7877000d0
    gsp  (53)   =       9.4756338d0
    gpp  (53)   =      10.5149673d0
    gp2  (53)   =       8.1832762d0
    hsp  (53)   =       4.7971922d0
    return
  end subroutine Parameters_for_PKA
