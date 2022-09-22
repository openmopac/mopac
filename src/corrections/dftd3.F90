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

  double precision function dftd3(l_grad, dxyz)
!
! dftd3 is the "D3" method of Grimme, et.al, as described in "A consistent and accurate ab initio
! parametrization of density functional dispersion correction (DFT-D) for the 94 elements H-Pu"
! S. Grimme, J. Antony, S. Ehrlich, H. Krieg, THE JOURNAL OF CHEMICAL PHYSICS 132, 154104 .2010.
!
!  This subroutine is based on materials provided by S. Grimme.
!
!
  use funcon_C, only : a0, fpc_9, fpc_2
  use molkst_C, only : numat, keywrd, E_hb, E_disp, numcal, method_pm8, method_pm6_org, l123
  use common_arrays_C, only: coord, nat, tvec
  use chanel_C, only : iw
  use elemts_C, only: elemnt
  use parameters_C, only: par7, par8, par9, par10
  implicit none
  double precision, intent (inout) ::  dxyz(3, numat)
  logical, intent (in):: l_grad
!
! local and dummy variables
!
  integer, parameter :: max_elem = 94, maxc = 5 ! maximum coordination number references per element
  integer :: i, mxc(max_elem), icalcn = -1
  logical :: D3H4, first = .true.
  double precision :: au_to_kcal, e6, e8, rs6, rs8, s6, s18, store_tvec(3,3), &
    alp6, alp8, rs18, alp, ehb, hbscale
  double precision :: &
    c6ab(max_elem, max_elem, maxc, maxc, 3), & ! C6 for all element pairs
    rcov(max_elem),                          & ! covalent radii
    r2r4(max_elem),                          & ! atomic <r^2>/<r^4> values
    r0ab(max_elem, max_elem)                   ! cut - off radii for all element pairs
    double precision, allocatable :: &
    dxyz_temp(:,:),                          &  ! Contribution to gradient
    store_coord(:,:)                                    ! Coordinates in au

  save
!
! PBE0/def2 - QZVP atomic values
!
      data r2r4 / &
       8.0589d0,  3.4698d0, 29.0974d0, 14.8517d0, 11.8799d0,  7.8715d0,  5.5588d0, &
       4.7566d0,  3.8025d0,  3.1036d0, 26.1552d0, 17.2304d0, 17.7210d0, 12.7442d0, &
       9.5361d0,  8.1652d0,  6.7463d0,  5.6004d0, 29.2012d0, 22.3934d0, 19.0598d0, &
      16.8590d0, 15.4023d0, 12.5589d0, 13.4788d0, 12.2309d0, 11.2809d0, 10.5569d0, &
      10.1428d0,  9.4907d0, 13.4606d0, 10.8544d0,  8.9386d0,  8.1350d0,  7.1251d0, &
       6.1971d0, 30.0162d0, 24.4103d0, 20.3537d0, 17.4780d0, 13.5528d0, 11.8451d0, &
      11.0355d0, 10.1997d0,  9.5414d0,  9.0061d0,  8.6417d0,  8.9975d0, 14.0834d0, &
      11.8333d0, 10.0179d0,  9.3844d0,  8.4110d0,  7.5152d0, 32.7622d0, 27.5708d0, &
      23.1671d0, 21.6003d0, 20.9615d0, 20.4562d0, 20.1010d0, 19.7475d0, 19.4828d0, &
      15.6013d0, 19.2362d0, 17.4717d0, 17.8321d0, 17.4237d0, 17.1954d0, 17.1631d0, &
      14.5716d0, 15.8758d0, 13.8989d0, 12.4834d0, 11.4421d0, 10.2671d0,  8.3549d0, &
       7.8496d0,  7.3278d0,  7.4820d0, 13.5124d0, 11.6554d0, 10.0959d0,  9.7340d0, &
       8.8584d0,  8.0125d0, 29.8135d0, 26.3157d0, 19.1885d0, 15.8542d0, 16.1305d0, &
      15.6161d0, 15.1226d0, 16.1576d0/
!
! covalent radii (taken from Pyykko and Atsumid0, Chem. Eur. J. 15d0, 2009d0, 188 - 197)
! values for metals decreased by 10 %
!
      data rcov/ &
      0.32d0, 0.46d0, 1.20d0, 0.94d0, 0.77d0, 0.75d0, 0.71d0, 0.63d0, 0.64d0, 0.67d0, &
      1.40d0, 1.25d0, 1.13d0, 1.04d0, 1.10d0, 1.02d0, 0.99d0, 0.96d0, 1.76d0, 1.54d0, &
      1.33d0, 1.22d0, 1.21d0, 1.10d0, 1.07d0, 1.04d0, 1.00d0, 0.99d0, 1.01d0, 1.09d0, &
      1.12d0, 1.09d0, 1.15d0, 1.10d0, 1.14d0, 1.17d0, 1.89d0, 1.67d0, 1.47d0, 1.39d0, &
      1.32d0, 1.24d0, 1.15d0, 1.13d0, 1.13d0, 1.08d0, 1.15d0, 1.23d0, 1.28d0, 1.26d0, &
      1.26d0, 1.23d0, 1.32d0, 1.31d0, 2.09d0, 1.76d0, 1.62d0, 1.47d0, 1.58d0, 1.57d0, &
      1.56d0, 1.55d0, 1.51d0, 1.52d0, 1.51d0, 1.50d0, 1.49d0, 1.49d0, 1.48d0, 1.53d0, &
      1.46d0, 1.37d0, 1.31d0, 1.23d0, 1.18d0, 1.16d0, 1.11d0, 1.12d0, 1.13d0, 1.32d0, &
      1.30d0, 1.30d0, 1.36d0, 1.31d0, 1.38d0, 1.42d0, 2.01d0, 1.81d0, 1.67d0, 1.58d0, &
      1.52d0, 1.53d0, 1.54d0, 1.55d0/
     if (icalcn /= numcal) then
        icalcn = numcal
        if (first) then
          first = .false.
          rcov = 4.d0/3.d0*rcov/a0
          au_to_kcal = fpc_9*fpc_2/a0
          do i = 1, max_elem
            r2r4(i) = sqrt(0.5d0*r2r4(i)*sqrt(float(i)))
          end do
          call setr0ab(max_elem, a0, r0ab)
          call copyc6(maxc, max_elem, c6ab, mxc)
        end if
      end if
      D3H4 = (method_pm6_org .or. method_PM8 .or. index(keywrd, "D3H4") + index(keywrd, "D3(H4)") /= 0)
      if (D3H4) then
! The D3H4 version of the dispersion
! Used in PM6-D3H4 and its variants PM6-D3H4X, PM6-D3(H4)
! I'VE CHECKED THAT THIS SETUP & THE PARAMETER VALUES READ FROM
! parameters_for_PM6 YIELD CORRECT PM6-D3H4 energies
        s6   = par7
        alp  = par8
        rs6 = par9
        s18 = par10
        rs18 = 1.d-10
      else ! hard-wired
!
! Grimme, S. (2012). "Supramolecular Binding Thermodynamics
! by Dispersion-Corrected Density Functional Theory." Chem. Eur. J.: 9955:9964.
!
        s6   = 1.0d0
        rs18 = 1.0d0
        alp = 14.0d0
        rs6 = 1.560d0   ! rs6 = s_(r,6) in Grimme's paper
        s18 = 1.009d0   ! s18 = s_8 in Grimme's paper
      end if
      hbscale = 1.301d0
      rs8   = rs18
      alp6  = alp
      alp8  = alp + 2.d0
!
!   Switch from MOPAC (convert coordinates from Angstroms to au)
!
      store_tvec = tvec
      allocate(dxyz_temp(3,numat*l123), store_coord(3,numat))
      store_coord = coord(:,:numat)
      coord(:,:numat) = coord(:,:numat)/a0
      tvec = tvec/a0
      call edisp(max_elem, maxc, numat, nat, c6ab, mxc, r2r4, r0ab, rcov, rs6, rs8, alp6, alp8, e6, e8)
      e6 = e6*s6
      e8 = e8*s18
      E_disp = (- e6 - e8)*au_to_kcal
      dxyz_temp = 0.0d0
!
! HBOND
!
      if (.not. D3H4) then
        call hbsimple(numat, nat, coord, hbscale, ehb, l_grad, dxyz_temp)
      else
        ehb = 0.d0
      end if
      E_hb   = ehb*au_to_kcal
      dftd3  = E_disp + E_hb
      if(l_grad)then
        call gdisp(r0ab, rs6, alp6, c6ab, s6, s18,mxc, r2r4, rcov, rs8, alp8, dxyz_temp)
        dxyz = dxyz + 2.d0*dxyz_temp*au_to_kcal
        if (index(keywrd, " DERIV") > 0) then
          write (iw, '(/16X,a)')"GRIMME'S D3 CORRECTIONS"
          write (iw, '(" NUMBER  ATOM  ",5X,"X",12X,"Y",12X,"Z",/)')
          write (iw, '(I6,4x,a2,F13.6,2F13.6)') (i, elemnt(nat(i)), dxyz_temp(:,i)*2.d0*au_to_kcal, i = 1,numat)
        end if
      end if
      coord(:,:numat) = store_coord
      tvec = store_tvec
      deallocate(dxyz_temp, store_coord)
  end function dftd3
