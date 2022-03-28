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

      module chanel_C
!dec$ attributes dllexport :: iend, end_fn
      integer :: &
     & iw0       = -1, & !    Abbrev. output, for GUI (By default, not used)
     & ifiles_1  =  1, &
     & input     =  2, & !    Input               <Filename>.dat         FORMATTED
     & iscr      =  3, & !    Scratch             <Filename>.scr         UNFORMATTED
     & irot      =  2, & !    Two-electron integrals used in analytical work
     & iend      =  4, & !    End                 <Filename>.end         FORMATTED
     & ires      =  9, & !    Restart             <Filename>.res         UNFORMATTED
     & iden      = 10, & !    Density             <Filename>.den         UNFORMATTED
     & ilog      = 31, & !    Log                 <Filename>.log         FORMATTED
     & iarc      = 12, & !    Archive             <Filename>.arc         FORMATTED
     & igpt      = 13, & !    Graphics            <Filename>.gpt         UNFORMATTED
     & iext      = 14, & !    EXTERNAL params      Defined by EXTERNAL    FORMATTED
     & iesr      = 15, & !    ESP restart         <Filename>.esr         UNFORMATTED
     & isyb      = 16, & !    SYBYL file          <Filename>.syb         FORMATTED
     & isol      = 17, & !    SOL map in MEP      <Filename>.sol         FORMATTED
     & ibrz      = 18, & !    Brillouin Zone      <Filename>.brz         FORMATTED
     & iump      = 20, & !    Grid map            <Filename>.ump         FORMATTED
     & iesp      = 21, & !    Electrostatic map   <Filename>.esp         FORMATTED
     & imep      = 22, & !    MEP map             <Filename>.mep         FORMATTED
     & ixyz      = 23, & !    XYZ output file     <Filename>.xyz         Formatted
     & param_out = 24, & !
     & ir        = 25, & !    Input (During run)  <Filename>.temp        FORMATTED
     & iw        = 26, & !    Output              <Filename>.out         FORMATTED
     & isetup    = 27, &
     & lbfgs_it  = 28, &
     & idaf      = 29, & !    DAF (polarization)  29   <Filename>.pol         FORMATTED
     & iwc       = 30    !    COSMO output (used by keyword COSWRT)

      logical :: log

      integer, parameter :: filename_len = 241

      character(len=filename_len) :: &
     & job_fn,       &
     & input_fn,     &
     & output_fn,    &
     & restart_fn,   &
     & density_fn,   &
     & log_fn,       &
     & end_fn,       &
     & ext_fn,       &
     & brillouin_fn, &
     & esp_fn,       &
     & esr_fn,       &
     & ump_fn,       &
     & mep_fn,       &
     & pol_fn,       &
     & gpt_fn,       &
     & xyz_fn,       &
     & syb_fn,       &
     & cosmo_fn,     &
     & archive_fn

      integer, dimension(145) :: ifilen, ioda
      integer :: irecln = 1023, irecst
      end module chanel_C
