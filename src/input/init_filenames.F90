! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

 subroutine init_filenames
    use molkst_C, only: jobnam, line
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
  end subroutine init_filenames
