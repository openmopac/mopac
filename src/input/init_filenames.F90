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

 subroutine init_filenames
    use molkst_C, only: jobnam, line, gui
    use chanel_C, only : output_fn, restart_fn, brillouin_fn, &
     & density_fn, log_fn, end_fn, archive_fn, esp_fn, ump_fn, &
     mep_fn, pol_fn, gpt_fn, esr_fn, input_fn, xyz_fn, syb_fn, &
     cosmo_fn
    implicit none
    integer :: text_length
    text_length = len_trim (jobnam)
    line = trim(jobnam)
    call upcase(line, len_trim(jobnam))
    if (text_length > 3) then
      if (jobnam(text_length - 3:text_length - 3) == ".") then
        if ((line(text_length - 3:text_length) == ".MOP") .or. &
            (line(text_length - 3:text_length) == ".DAT") .or. &
            (line(text_length - 3:text_length) == ".ARC") .or. &
            (line(text_length - 3:text_length) == ".PDB") .or. &
            (line(text_length - 3:text_length) == ".ENT") .or. &
            (line(text_length - 3:text_length) == ".NEW"))      &
        text_length = text_length - 4
      end if
    end if
!
!  Set up the name of the files that are related to the input file name
!
    input_fn      = jobnam(1:text_length) // ".temp"
    output_fn     = jobnam(1:text_length) // ".out"
    restart_fn    = jobnam(1:text_length) // ".res"
    density_fn    = jobnam(1:text_length) // ".den"
    log_fn        = jobnam(1:text_length) // ".log"
    end_fn        = jobnam(1:text_length) // ".end"
    archive_fn    = jobnam(1:text_length) // ".arc"
    brillouin_fn  = jobnam(1:text_length) // ".brz"
    esp_fn        = jobnam(1:text_length) // ".esp"
    ump_fn        = jobnam(1:text_length) // ".ump"
    mep_fn        = jobnam(1:text_length) // ".mep"
    pol_fn        = jobnam(1:text_length) // ".pol"
    gpt_fn        = jobnam(1:text_length) // ".gpt"
    esr_fn        = jobnam(1:text_length) // ".esr"
    xyz_fn        = jobnam(1:text_length) // ".xyz"
    syb_fn        = jobnam(1:text_length) // ".syb"
    cosmo_fn      = jobnam(1:text_length) // ".cos"
    if (gui) output_fn     = "OUTPUT file"
  end subroutine init_filenames
